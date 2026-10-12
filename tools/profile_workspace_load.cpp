#include "profile_process.hpp"
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <algorithm>
#include <set>

namespace {
using Counts = std::map<std::string, std::size_t>;
void add(Counts &counts, const std::string &name, std::size_t value) {
    counts[name] = rfc::machine::checked_size_add(counts[name], value);
}
void byte_owners(Counts &counts, const std::string &prefix, const gadget::NativeRecord &record,
                 std::vector<gadget::RawByteSource> &owners) {
    for (const auto &source : record.raw_sources()) {
        add(counts, prefix + ".source_references", 1);
        if (std::none_of(owners.begin(), owners.end(),
                        [&](const auto &owner) { return source.shares_owner(owner); })) {
            owners.push_back(source);
            add(counts, prefix + ".unique_byte_owners", 1);
            add(counts, prefix + ".raw_owner_bytes", source.bytes().size());
        }
    }
}
void print(const Counts &counts) {
    for (const auto &[name, value] : counts)
        std::cout << "value\t" << name << '\t' << value << '\n';
}
void workspace_owners(const gadget::AssayWorkspace &workspace) {
    Counts counts;
    std::vector<gadget::RawByteSource> owners;
    byte_owners(counts, "workspace", workspace.native_publication(), owners);
    std::set<const rfc::storage::TemporalEndpoint *> endpoints;
    std::set<const std::vector<std::byte> *> components;
    const auto assays = workspace.native_publication().at("manifest").at("assays");
    for (std::size_t i = 0; i < assays.size(); ++i) {
        const auto assay = assays.at(i);
        const auto id = assay.at("declaration").at("id").text();
        for (std::size_t p = 0; p < assay.at("partitions").size(); ++p) {
            const auto study = workspace.study(id, p);
            byte_owners(counts, "workspace", study->native_publication(), owners);
            const auto endpoint = study->endpoint();
            if (endpoints.insert(endpoint.get()).second) {
                add(counts, "workspace.study_endpoint_owners", 1);
                for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k) {
                    const auto &grade = endpoint->source.grade({static_cast<std::uint16_t>(k)});
                    add(counts, "workspace.source.cells", grade.cell_count());
                    add(counts, "workspace.source.occurrences", grade.incidence_count());
                }
                for (const auto &field : endpoint->tensors) {
                    add(counts, "workspace.fields.tensors", 1);
                    add(counts, "workspace.fields.logical_coordinates", field.tensor.domain().logical_size());
                    add(counts, "workspace.fields.supported_coordinates", field.tensor.support().count());
                    add(counts, "workspace.fields.nonzero_entries", field.tensor.stored_nonzeros());
                }
            }
            study->native_history()->history().visit_components([&](const auto &part) {
                if (components.insert(part.bytes.get()).second) {
                    add(counts, "workspace.history.components", 1);
                    add(counts, "workspace.history.component_bytes", part.bytes->size());
                    add(counts, "workspace.history.component_capacity_bytes", part.bytes->capacity());
                }
            });
        }
    }
    print(counts);
}
} // namespace

int main(int argc, char **argv) {
    try {
        gadget::fs::path input, output, request_path, result_path;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            if ((option != "--input" && option != "--rewrite" &&
                 option != "--request" && option != "--result") ||
                !seen.insert(option).second || i + 1 == argc)
                throw std::invalid_argument("usage: gadget_workspace_load_profile --input FILE [--rewrite NEW_FILE] [--request FILE --result NEW_FILE]");
            const gadget::fs::path path(argv[++i]);
            if (option == "--input") input = path;
            else if (option == "--rewrite") output = path;
            else if (option == "--request") request_path = path;
            else result_path = path;
        }
        if (input.empty()) throw std::invalid_argument("workspace profile requires --input");
        if (request_path.empty() != result_path.empty())
            throw std::invalid_argument("workspace query profile requires both --request and --result");
        for (const auto &path : {output, result_path})
            if (!path.empty() && gadget::fs::exists(path))
                throw std::invalid_argument("profile publication destination already exists");
        if (!output.empty() && !result_path.empty() &&
            gadget::fs::absolute(output).lexically_normal() == gadget::fs::absolute(result_path).lexically_normal())
            throw std::invalid_argument("profile publications require distinct destinations");
        {
            const rfc::storage::ReadOnlyMappedFile image(input);
            if (gadget::detect_native_data(image.bytes()) != gadget::NativeDataFormat::TemporalSeries)
                throw std::invalid_argument("workspace profile requires a canonical RFCXS input");
            std::cout << "value\tworkspace.input_bytes\t" << image.size() << std::endl;
        }
        gadget::profiling::cpu_clock();
        using Clock = gadget::profiling::ProcessClock;
        const auto start = Clock::now();
        auto stage = start;
        auto workspace = std::make_unique<gadget::AssayWorkspace>(input,
            [&](std::string_view name) {
                gadget::profiling::phase(name, stage);
                stage = Clock::now();
            });
        gadget::profiling::phase("admission", start);
        stage = Clock::now();
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
        stage = Clock::now();
        workspace_owners(*workspace);
        gadget::profiling::phase("account_owners", stage);
        if (!request_path.empty()) {
            stage = Clock::now();
            auto request = gadget::read_record(request_path).native();
            gadget::profiling::phase("request_admission", stage);
            stage = Clock::now();
            {
                auto result = workspace->execute_record(request);
                gadget::profiling::phase("execute_record", stage);
                stage = Clock::now();
                {
                    Counts counts;
                    std::vector<gadget::RawByteSource> owners;
                    byte_owners(counts, "result", result, owners);
                    print(counts);
                }
                gadget::profiling::phase("account_result", stage);
                stage = Clock::now();
                result.write(result_path);
                gadget::profiling::phase("publication", stage);
                std::cout << "value\tresult.archive_bytes\t" << gadget::fs::file_size(result_path) << std::endl;
                stage = Clock::now();
            }
            gadget::profiling::phase("release_result", stage);
            stage = Clock::now();
        }
        if (!request_path.empty()) gadget::profiling::phase("release_request", stage);
        if (!output.empty()) {
            stage = Clock::now();
            workspace->write_native(output);
            gadget::profiling::phase("rewrite", stage);
            stage = Clock::now();
            const rfc::storage::ReadOnlyMappedFile before(input), after(output);
            if (before.size() != after.size() ||
                !std::equal(before.bytes().begin(), before.bytes().end(), after.bytes().begin()))
                throw std::invalid_argument("workspace rewrite differs from its original input");
            std::cout << "value\tworkspace.verified_bytes\t" << after.size() << std::endl;
            gadget::profiling::phase("verify_bytes", stage);
        }
        stage = Clock::now();
        workspace.reset();
        gadget::profiling::phase("release", stage);
        gadget::profiling::phase("total", start);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
