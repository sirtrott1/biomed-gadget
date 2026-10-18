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
#include <unistd.h>
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
    std::ifstream status("/proc/self/statm");
    std::size_t pages{};
    const auto page_bytes = sysconf(_SC_PAGESIZE);
    if (!(status >> pages) || page_bytes <= 0 ||
        pages > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(page_bytes))
        throw std::runtime_error("process address space observation failed");
    const auto current = pages * static_cast<std::size_t>(page_bytes);
    if (*bytes > std::numeric_limits<std::size_t>::max() - current)
        throw std::invalid_argument("process memory budget exceeds address space");
    struct rlimit limit{};
    if (getrlimit(RLIMIT_AS, &limit))
        throw std::runtime_error("process memory limit observation failed");
    const auto requested = static_cast<rlim_t>(current + *bytes);
    limit.rlim_cur = std::min(limit.rlim_cur, requested);
    if (setrlimit(RLIMIT_AS, &limit))
        throw std::runtime_error(std::string("process memory limit failed: ") + std::strerror(errno));
    return limit.rlim_cur > current ? limit.rlim_cur - current : 0;
#else
    if (bytes) throw std::runtime_error("process memory limits require Linux");
    return {};
#endif
}
}
