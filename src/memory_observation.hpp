#pragma once
#include <cstddef>
#include <optional>
#include <string>
namespace rfc::runtime { struct HardwareInventory; struct LinuxHardwareProbe; }
namespace gadget::hardware_detail {
bool memory_path(const std::string &path);
std::optional<std::string> read_observation(const std::string &path, std::size_t limit);
rfc::runtime::HardwareInventory memory_inventory(rfc::runtime::LinuxHardwareProbe probe);
std::optional<std::size_t> memory_headroom();
}
