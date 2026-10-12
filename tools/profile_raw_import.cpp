#include "profile_process.hpp"
#include <gadget/assay.hpp>

int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("usage: gadget_raw_import_profile SPECIFICATION NEW_DIRECTORY");
        gadget::profiling::cpu_clock();
        const auto start = gadget::profiling::ProcessClock::now();
        auto stage = start;
        auto workspace = gadget::import_assay_workspace(argv[1], argv[2], [&](std::string_view name) {
            gadget::profiling::phase(name, stage);
            stage = gadget::profiling::ProcessClock::now();
        });
        gadget::profiling::phase("import", start);
        stage = gadget::profiling::ProcessClock::now();
        workspace.reset();
        gadget::profiling::phase("release", stage);
        gadget::profiling::phase("total", start);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
