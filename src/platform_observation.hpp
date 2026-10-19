#pragma once
#include <rfc/runtime/hardware_inventory.hpp>

namespace gadget::hardware_detail {
struct DarwinHardwareProbe {
    std::optional<std::vector<std::size_t>> scheduler_cpus;
    std::string cpu_model;
    std::optional<std::size_t> ram_total_bytes, page_bytes, free_pages, inactive_pages;
};
rfc::runtime::HardwareInventory darwin_inventory(
    const DarwinHardwareProbe &probe, const rfc::runtime::HardwareProbeLimits &limits = {});
rfc::runtime::HardwareInventory local_inventory();
std::optional<std::size_t> darwin_memory_headroom();
} // namespace gadget::hardware_detail
