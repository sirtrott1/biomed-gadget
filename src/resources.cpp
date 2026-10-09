#include <gadget/resources.hpp>
#include <limits>
#include <mutex>

namespace gadget {
void CpuResources::process() const {
#if defined(__linux__)
    if (::getpid() != process_id_)
        throw std::logic_error("inherited CPU resources cannot be used after fork");
#endif
}
CpuResources::CpuResources(rfc::runtime::HardwareInventory hardware,
                           rfc::runtime::CpuResourcePolicy policy,
                           std::function<rfc::runtime::HardwareInventory()> refresh)
    : initial_(std::move(hardware)), refresh_(std::move(refresh)), coordinator_(initial_, policy) {
    if (!refresh_)
        throw std::invalid_argument("CPU resources require an observation refresh");
}
std::shared_ptr<CpuResources> CpuResources::local() {
    static std::mutex mutex;
    static std::weak_ptr<CpuResources> shared;
    std::lock_guard lock(mutex);
    if (auto current = shared.lock())
        return current;
    rfc::runtime::CpuResourcePolicy policy;
    auto hardware = rfc::runtime::detect_hardware();
    policy.worker_limit = hardware.allowed_cpus.size();
    policy.tile_limit = std::numeric_limits<std::size_t>::max();
    const auto headroom = hardware.advisory_memory_headroom();
    if (!headroom || !*headroom)
        throw std::invalid_argument("CPU resources require observed memory headroom");
    policy.declared_buffer_byte_limit = rfc::machine::checked_size(*headroom);
    auto current = std::make_shared<CpuResources>(std::move(hardware), policy,
                                                  [] { return rfc::runtime::detect_hardware(); });
    shared = current;
    return current;
}
CpuResources::CpuResources(const CpuResources &parent,
                           rfc::runtime::CpuResourceCoordinator assigned)
    : initial_(parent.initial_), refresh_(parent.refresh_), coordinator_(assigned.share()) {}
std::shared_ptr<CpuResources>
CpuResources::assigned(rfc::runtime::CpuResourceCoordinator &coordinator) const {
    process();
    return std::shared_ptr<CpuResources>(new CpuResources(*this, coordinator.share()));
}
rfc::runtime::HardwareInventory CpuResources::readmit() const {
    process();
    auto hardware = refresh_();
    coordinator_.readmit(hardware);
    return hardware;
}
rfc::runtime::CpuResourceCoordinator::Reservation CpuResources::reserve(std::size_t workers) {
    (void)readmit();
    return coordinator_.reserve(workers);
}
Json CpuResources::report() const {
    const auto hardware = readmit();
    auto result = Json::object();
    result.set("platform", Json::string(hardware.platform));
    result.set("cpu_model", Json::string(hardware.cpu_model));
    result.set("allowed_logical_threads",
               Json::integer(static_cast<std::int64_t>(hardware.allowed_cpus.size())));
    result.set("physical_cores",
               hardware.physical_cores()
                   ? Json::integer(static_cast<std::int64_t>(*hardware.physical_cores()))
                   : Json{});
    result.set("worker_capacity",
               Json::integer(static_cast<std::int64_t>(coordinator_.capacity())));
    result.set("worker_limit",
               Json::integer(static_cast<std::int64_t>(coordinator_.policy().worker_limit)));
    result.set("worker_basis",
               Json::string(coordinator_.policy().basis ==
                                    rfc::runtime::CpuWorkerBasis::AllowedLogicalThreads
                                ? "allowed_logical_threads"
                                : "observed_physical_core_count"));
    result.set("ram_total_bytes", hardware.ram_total_bytes
                                      ? Json::string(std::to_string(*hardware.ram_total_bytes))
                                      : Json{});
    result.set("advisory_memory_headroom_bytes",
               hardware.advisory_memory_headroom()
                   ? Json::string(std::to_string(*hardware.advisory_memory_headroom()))
                   : Json{});
    result.set("scope", Json::string("cooperating declared workers including each caller; no CPU "
                                     "pinning, external-process exclusion or heap quota"));
    result.set("gpu_compute_qualification", Json::string("not established by PCI observations"));
    return result;
}
} // namespace gadget
