#include "memory_observation.hpp"
#include "platform_observation.hpp"
#include <rfc/runtime/hardware_inventory.hpp>
#include <fstream>
namespace gadget::hardware_detail {
bool memory_path(const std::string &path) {
    return path == "/proc/meminfo" || path == "/proc/self/cgroup" ||
        path == "/proc/self/mountinfo" || path.ends_with("/memory.max") ||
        path.ends_with("/memory.current") || path.ends_with("/cpu.max");
}
std::optional<std::string> read_observation(const std::string &path, std::size_t limit) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return {};
    std::string text;
    char block[4096];
    while (file.read(block, sizeof(block)) || file.gcount()) {
        const auto count = static_cast<std::size_t>(file.gcount());
        if (text.size() > limit || count > limit - text.size())
            throw std::length_error("hardware file exceeds limit");
        text.append(block, count);
    }
    if (file.bad()) return {};
    return text;
}
rfc::runtime::HardwareInventory memory_inventory(rfc::runtime::LinuxHardwareProbe probe) {
    if (!probe.read_file) throw std::invalid_argument("memory observation requires a file reader");
    auto read = std::move(probe.read_file);
    probe.allowed_cpu_ids.reset();
    probe.visible_compute_nodes.clear();
    probe.read_file = [read = std::move(read)](const std::string &path) {
        return memory_path(path) ? read(path) : std::optional<std::string>{};
    };
    return rfc::runtime::probe_linux_hardware(probe);
}
std::optional<std::size_t> memory_headroom() {
#if defined(__linux__)
    rfc::runtime::LinuxHardwareProbe probe;
    probe.read_file = [limit = probe.limits.max_file_bytes](const std::string &path) {
        return read_observation(path, limit);
    };
    return memory_inventory(std::move(probe)).advisory_memory_headroom();
#elif defined(__APPLE__)
    return darwin_memory_headroom();
#else
    return rfc::runtime::detect_hardware().advisory_memory_headroom();
#endif
}
}
