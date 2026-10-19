#include "platform_observation.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <mach/processor_info.h>
#include <sys/sysctl.h>
#endif

namespace gadget::hardware_detail {
rfc::runtime::HardwareInventory darwin_inventory(
    const DarwinHardwareProbe &probe, const rfc::runtime::HardwareProbeLimits &limits) {
    using rfc::runtime::ObservationStatus;
    if (!limits.max_cpus)
        throw std::invalid_argument("hardware observation requires a positive CPU record limit");
    rfc::runtime::HardwareInventory result;
    result.platform = "Darwin";
    result.cpu_model = probe.cpu_model;
    result.ram_total_bytes = probe.ram_total_bytes;
    if (probe.scheduler_cpus) {
        if (probe.scheduler_cpus->size() > limits.max_cpus)
            throw std::length_error("observed processor count exceeds CPU record limit");
        auto cpus = *probe.scheduler_cpus;
        std::sort(cpus.begin(), cpus.end());
        if (std::adjacent_find(cpus.begin(), cpus.end()) != cpus.end())
            throw std::invalid_argument("duplicate observed processor slot");
        result.affinity_observed = !cpus.empty();
        for (const auto id : cpus)
            result.allowed_cpus.push_back({id, {}, {}, {}});
    }
    result.observations.push_back({
        "host_processor_info", result.affinity_observed ? ObservationStatus::Observed
                                                       : ObservationStatus::Unavailable,
        "online scheduler processor slots; no CPU placement or physical core topology"});
    if (probe.page_bytes && *probe.page_bytes && probe.free_pages && probe.inactive_pages) {
        const auto pages = rfc::machine::checked_size_add(*probe.free_pages, *probe.inactive_pages);
        auto bytes = rfc::machine::checked_size_mul(pages, *probe.page_bytes);
        if (probe.ram_total_bytes)
            bytes = std::min(bytes, *probe.ram_total_bytes);
        result.ram_available_bytes = bytes;
    }
    result.observations.push_back({
        "host_statistics64", result.ram_available_bytes ? ObservationStatus::Observed
                                                       : ObservationStatus::Unavailable,
        "advisory free plus inactive pages; speculative pages already belong to free pages"});
    return result;
}

#if defined(__APPLE__)
namespace {
struct HostPort {
    host_t value{mach_host_self()};
    HostPort() = default;
    HostPort(const HostPort &) = delete;
    HostPort &operator=(const HostPort &) = delete;
    ~HostPort() {
        if (MACH_PORT_VALID(value))
            (void)mach_port_deallocate(mach_task_self(), value);
    }
};
struct ProcessorInfo {
    processor_info_array_t data{};
    mach_msg_type_number_t count{};
    ProcessorInfo() = default;
    ProcessorInfo(const ProcessorInfo &) = delete;
    ProcessorInfo &operator=(const ProcessorInfo &) = delete;
    ~ProcessorInfo() {
        if (data)
            (void)vm_deallocate(mach_task_self(), reinterpret_cast<vm_address_t>(data),
                                static_cast<vm_size_t>(count) * sizeof(integer_t));
    }
};
DarwinHardwareProbe observe_darwin(bool cpu) {
    DarwinHardwareProbe probe;
    HostPort host;
    std::uint64_t total{};
    std::size_t length = sizeof(total);
    if (sysctlbyname("hw.memsize", &total, &length, nullptr, 0) == 0 &&
        length == sizeof(total) && total <= std::numeric_limits<std::size_t>::max())
        probe.ram_total_bytes = static_cast<std::size_t>(total);
    vm_size_t page_bytes{};
    vm_statistics64_data_t memory{};
    mach_msg_type_number_t memory_count = HOST_VM_INFO64_COUNT;
    if (host_page_size(host.value, &page_bytes) == KERN_SUCCESS &&
        host_statistics64(host.value, HOST_VM_INFO64,
                          reinterpret_cast<host_info64_t>(&memory), &memory_count) == KERN_SUCCESS &&
        memory_count >= (offsetof(vm_statistics64_data_t, inactive_count) +
                         sizeof(memory.inactive_count)) / sizeof(integer_t)) {
        probe.page_bytes = page_bytes;
        probe.free_pages = memory.free_count;
        probe.inactive_pages = memory.inactive_count;
    }
    if (!cpu)
        return probe;
    char model[256]{};
    length = sizeof(model);
    if (sysctlbyname("machdep.cpu.brand_string", model, &length, nullptr, 0) != 0) {
        length = sizeof(model);
        if (sysctlbyname("hw.model", model, &length, nullptr, 0) != 0)
            length = 0;
    }
    if (length && length <= sizeof(model))
        probe.cpu_model.assign(model, std::find(model, model + length, '\0'));
    ProcessorInfo processors;
    natural_t cpu_count{};
    if (host_processor_info(host.value, PROCESSOR_BASIC_INFO, &cpu_count,
                            &processors.data, &processors.count) == KERN_SUCCESS) {
        const auto expected = rfc::machine::checked_size_mul(cpu_count, PROCESSOR_BASIC_INFO_COUNT);
        const rfc::runtime::HardwareProbeLimits limits;
        if (processors.count != expected || cpu_count > limits.max_cpus)
            throw std::length_error("invalid observed processor record count");
        probe.scheduler_cpus.emplace();
        for (std::size_t i = 0; i < cpu_count; ++i) {
            const auto info = reinterpret_cast<const processor_basic_info *>(
                processors.data + i * PROCESSOR_BASIC_INFO_COUNT);
            if (info->running) {
                if (info->slot_num < 0)
                    throw std::invalid_argument("invalid observed processor slot");
                probe.scheduler_cpus->push_back(static_cast<std::size_t>(info->slot_num));
            }
        }
    }
    return probe;
}
} // namespace
#endif

rfc::runtime::HardwareInventory local_inventory() {
#if defined(__APPLE__)
    return darwin_inventory(observe_darwin(true));
#else
    return rfc::runtime::detect_hardware();
#endif
}
std::optional<std::size_t> darwin_memory_headroom() {
#if defined(__APPLE__)
    return darwin_inventory(observe_darwin(false)).advisory_memory_headroom();
#else
    return {};
#endif
}
} // namespace gadget::hardware_detail
