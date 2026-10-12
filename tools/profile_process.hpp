#pragma once
#include <gadget/profile_process.hpp>
#include <ctime>

namespace gadget::profiling {
struct ProcessClock {
    struct time_point {
        Clock::time_point wall;
        std::clock_t cpu;
    };
    static time_point now() { return {Clock::now(), std::clock()}; }
};
inline void cpu_clock() {
    std::cout << "value\tprocess.cpu_ticks_per_second\t" << CLOCKS_PER_SEC << std::endl;
}
inline void phase(std::string_view name, ProcessClock::time_point start) {
    const auto cpu = std::clock();
    phase(name, start.wall);
    if (start.cpu != std::clock_t{-1} && cpu != std::clock_t{-1} && cpu >= start.cpu) {
        static std::size_t index{};
        std::cout << "value\tphase.cpu." << index++ << '.' << name << ".ticks\t"
                  << cpu - start.cpu << std::endl;
    }
}
} // namespace gadget::profiling
