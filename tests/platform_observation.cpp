#include "../src/platform_observation.hpp"
#include <gadget/resources.hpp>
#include <iostream>
#include <limits>
#if defined(__linux__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class E, class F> void rejects(F action, const char *message) {
    bool refused = false;
    try { action(); } catch (const E &) { refused = true; }
    check(refused, message);
}
void observations() {
    using namespace gadget::hardware_detail;
    DarwinHardwareProbe probe;
    probe.scheduler_cpus = std::vector<std::size_t>{11, 3, 7};
    probe.cpu_model = "authored scheduler observation";
    probe.ram_total_bytes = 128 * 4096;
    probe.page_bytes = 4096;
    probe.free_pages = 8;
    probe.inactive_pages = 16;
    const auto original = darwin_inventory(probe);
    check(original.platform == "Darwin" && original.affinity_observed &&
              original.allowed_cpus.size() == 3 && original.allowed_cpus[0].logical_id == 3 &&
              original.allowed_cpus[2].logical_id == 11,
          "scheduler observations lost original processor slots");
    check(!original.physical_cores() && original.caches.empty() && original.pci_devices.empty(),
          "scheduler observations invented physical topology or devices");
    check(original.advisory_memory_headroom() == 24 * 4096,
          "page observations changed the free plus inactive byte count");
    rfc::runtime::CpuResourcePolicy policy;
    policy.worker_limit = 2;
    policy.declared_buffer_byte_limit = 4096;
    gadget::CpuResources resources(original, policy, [&] { return darwin_inventory(probe); });
    check(resources.coordinator().capacity() == 2, "coordinator ignored worker policy");
    {
        auto reservation = resources.coordinator().reserve(1, 1024);
        check(resources.coordinator().available_workers() == 1 &&
                  resources.coordinator().available_declared_buffer_bytes() == 3072,
              "native observations bypassed canonical resource charges");
        rejects<rfc::runtime::CpuResourceBusy>([&] { (void)resources.reserve(2); },
                                              "occupied workers were admitted");
    }
    check(resources.coordinator().available_workers() == 2 &&
              resources.coordinator().available_declared_buffer_bytes() == 4096,
          "released resources retained their charges");
    probe.scheduler_cpus = std::vector<std::size_t>{3, 7};
    rejects<std::invalid_argument>([&] { (void)resources.readmit(); },
                                   "changed processor observations reused a coordinator");
    probe.scheduler_cpus.reset();
    const auto unknown = darwin_inventory(probe);
    check(!unknown.affinity_observed && unknown.allowed_cpus.empty(),
          "missing processor observation became a usable CPU set");
    rejects<std::invalid_argument>([&] { rfc::runtime::CpuResourceCoordinator cpu(unknown, policy); },
                                   "missing processors bypassed admission");
    probe.free_pages = 0;
    probe.inactive_pages = 0;
    check(darwin_inventory(probe).advisory_memory_headroom() == 0,
          "observed exhausted memory became unknown");
    probe.free_pages.reset();
    check(!darwin_inventory(probe).advisory_memory_headroom(),
          "missing page count became observed zero");
    probe.free_pages = 8;
    probe.page_bytes = 0;
    check(!darwin_inventory(probe).advisory_memory_headroom(), "zero page size was admitted");
    probe.page_bytes = 4096;
    probe.ram_total_bytes = 4096;
    check(darwin_inventory(probe).advisory_memory_headroom() == 4096,
          "advisory pages exceeded observed physical memory");
    probe.ram_total_bytes.reset();
    probe.free_pages = std::numeric_limits<std::size_t>::max();
    probe.inactive_pages = 1;
    rejects<std::overflow_error>([&] { (void)darwin_inventory(probe); },
                                "overflowed page sum was admitted");
    probe.inactive_pages = 0;
    rejects<std::overflow_error>([&] { (void)darwin_inventory(probe); },
                                "overflowed page bytes were admitted");
    probe.free_pages = 0;
    probe.scheduler_cpus = std::vector<std::size_t>{3, 3};
    rejects<std::invalid_argument>([&] { (void)darwin_inventory(probe); },
                                   "duplicate processor slots were admitted");
    probe.scheduler_cpus = std::vector<std::size_t>{3, 7};
    rfc::runtime::HardwareProbeLimits limits;
    limits.max_cpus = 1;
    rejects<std::length_error>([&] { (void)darwin_inventory(probe, limits); },
                              "processor record budget was ignored");
    limits.max_cpus = 0;
    rejects<std::invalid_argument>([&] { (void)darwin_inventory(probe, limits); },
                                   "zero processor record budget was admitted");
}
} // namespace

int main() try {
    observations();
    const auto cpu = gadget::CpuResources::local();
    const auto current = cpu->readmit();
    check(current.affinity_observed && !current.allowed_cpus.empty() &&
              current.advisory_memory_headroom().value_or(0) > 0,
          "native host observations are unavailable");
    {
        auto reservation = cpu->reserve(1);
        check(cpu->coordinator().available_workers() + 1 == cpu->coordinator().capacity(),
              "real host reservation bypassed admission");
    }
#if defined(__linux__) || defined(__APPLE__)
    const auto child = fork();
    if (child < 0)
        throw std::runtime_error("fork failed");
    if (!child) {
        try { (void)cpu->readmit(); }
        catch (const std::logic_error &) { _exit(0); }
        _exit(1);
    }
    int status{};
    check(waitpid(child, &status, 0) == child && WIFEXITED(status) && !WEXITSTATUS(status),
          "inherited resources were admitted after fork");
#endif
#if defined(__APPLE__)
    check(current.platform == "Darwin" &&
              gadget::hardware_detail::darwin_memory_headroom().has_value(),
          "macOS did not use native scheduler and memory observations");
#endif
    std::cout << "PASS platform CPU and memory observations, missing and exhausted support, "
                 "overflow refusals, canonical admission and live host execution\n";
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
