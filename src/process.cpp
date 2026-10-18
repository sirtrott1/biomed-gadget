#include <gadget/process.hpp>
#include "memory_observation.hpp"
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#if defined(__linux__)
#include <sys/resource.h>
#endif
namespace gadget {
std::optional<std::size_t> limit_process_memory(std::optional<std::size_t> bytes) {
    if (bytes && !*bytes) throw std::invalid_argument("process memory budget must be positive");
#if defined(__linux__)
    const auto available = hardware_detail::memory_headroom();
    if (!bytes && available) bytes = *available / 2;
    if (!bytes) return {};
    if (available && *bytes > *available - *available / 4)
        throw std::invalid_argument("process memory budget exceeds available memory reserve");
    std::ifstream status("/proc/self/status");
    std::string key, unit;
    std::size_t kib{};
    bool observed = false;
    while (status >> key) {
        if (key == "VmData:") {
            observed = static_cast<bool>(status >> kib >> unit) && unit == "kB";
            break;
        }
        status.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    if (!observed || kib > std::numeric_limits<std::size_t>::max() / 1024)
        throw std::runtime_error("process data allocation observation failed");
    const auto current = kib * 1024;
    if (*bytes > std::numeric_limits<std::size_t>::max() - current)
        throw std::invalid_argument("process memory budget exceeds data allocation range");
    struct rlimit limit{};
    if (getrlimit(RLIMIT_DATA, &limit))
        throw std::runtime_error("process memory limit observation failed");
    const auto requested = static_cast<rlim_t>(current + *bytes);
    limit.rlim_cur = std::min(limit.rlim_cur, requested);
    if (setrlimit(RLIMIT_DATA, &limit))
        throw std::runtime_error(std::string("process memory limit failed: ") + std::strerror(errno));
    return limit.rlim_cur > current ? limit.rlim_cur - current : 0;
#else
    if (bytes) throw std::runtime_error("process memory limits require Linux");
    return {};
#endif
}
}
