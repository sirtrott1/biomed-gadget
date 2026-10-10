#include "profile_process.hpp"
#include <gadget/assay.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <algorithm>
#include <set>

int main(int argc, char **argv) {
    try {
        gadget::fs::path input, output;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            if ((option != "--input" && option != "--rewrite") ||
                !seen.insert(option).second || i + 1 == argc)
                throw std::invalid_argument("usage: gadget_workspace_load_profile --input FILE [--rewrite NEW_FILE]");
            (option == "--input" ? input : output) = argv[++i];
        }
        if (input.empty()) throw std::invalid_argument("workspace profile requires --input");
        {
            const rfc::storage::ReadOnlyMappedFile image(input);
            if (gadget::detect_native_data(image.bytes()) != gadget::NativeDataFormat::TemporalSeries)
                throw std::invalid_argument("workspace profile requires a canonical RFCXS input");
        }
        const auto start = gadget::profiling::Clock::now();
        auto stage = start;
        auto workspace = std::make_unique<gadget::AssayWorkspace>(input,
            [&](std::string_view name) {
                gadget::profiling::phase(name, stage);
                stage = gadget::profiling::Clock::now();
            });
        gadget::profiling::phase("admission", start);
        stage = gadget::profiling::Clock::now();
        {
            const auto assays = workspace->native_publication().at("manifest").at("assays");
            std::size_t partitions{};
            for (std::size_t i = 0; i < assays.size(); ++i)
                partitions = rfc::machine::checked_size_add(partitions, assays.at(i).at("partitions").size());
            std::cout << "value\tworkspace.assays\t" << assays.size() << '\n'
                      << "value\tworkspace.partitions\t" << partitions << '\n'
                      << "binding\tworkspace.digest\t" << workspace->digest() << std::endl;
        }
        gadget::profiling::phase("binding", stage);
        if (!output.empty()) {
            stage = gadget::profiling::Clock::now();
            workspace->write_native(output);
            gadget::profiling::phase("rewrite", stage);
            stage = gadget::profiling::Clock::now();
            const rfc::storage::ReadOnlyMappedFile before(input), after(output);
            if (before.size() != after.size() ||
                !std::equal(before.bytes().begin(), before.bytes().end(), after.bytes().begin()))
                throw std::invalid_argument("workspace rewrite differs from its original input");
            std::cout << "value\tworkspace.verified_bytes\t" << after.size() << std::endl;
            gadget::profiling::phase("verify_bytes", stage);
        }
        stage = gadget::profiling::Clock::now();
        workspace.reset();
        gadget::profiling::phase("release", stage);
        gadget::profiling::phase("total", start);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
