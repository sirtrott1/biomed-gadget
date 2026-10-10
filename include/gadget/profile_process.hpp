#pragma once
#include <rfc/machine/checked_size.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string_view>

namespace gadget::profiling {
using Clock = std::chrono::steady_clock;
struct ProcessMemory {
    std::optional<std::size_t> resident, peak;
};
inline ProcessMemory process_memory() {
    ProcessMemory result;
    std::ifstream input("/proc/self/status");
    std::string line;
    while (std::getline(input, line)) {
        const auto read = [&]() {
            std::istringstream value(line.substr(line.find(':') + 1));
            std::size_t kib{};
            std::string unit;
            if (!(value >> kib >> unit) || unit != "kB")
                throw std::invalid_argument("invalid process memory observation");
            return rfc::machine::checked_size_mul(kib, 1024);
        };
        if (line.starts_with("VmRSS:")) result.resident = read();
        if (line.starts_with("VmHWM:")) result.peak = read();
    }
    return result;
}
inline void phase(std::string_view name, Clock::time_point start) {
    const auto seconds = std::chrono::duration<double>(Clock::now() - start).count();
    const auto memory = process_memory();
    std::cout << "phase\t" << name << '\t' << std::setprecision(9) << seconds << '\t';
    if (memory.resident) std::cout << *memory.resident; else std::cout << "unavailable";
    std::cout << '\t';
    if (memory.peak) std::cout << *memory.peak; else std::cout << "unavailable";
    std::cout << std::endl;
}
} // namespace gadget::profiling
