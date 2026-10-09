#include <gadget/native_record.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <algorithm>
#include <chrono>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <type_traits>

namespace {
using namespace rfc;
using Clock = std::chrono::steady_clock;
void require(bool value, const char *message) {
    if (!value) throw std::invalid_argument(message);
}
struct ProcessMemory {
    std::optional<std::size_t> resident, peak;
};
ProcessMemory process_memory() {
    ProcessMemory result;
    std::ifstream input("/proc/self/status");
    std::string line;
    while (std::getline(input, line)) {
        const auto read = [&]() {
            std::istringstream value(line.substr(line.find(':') + 1));
            std::size_t kib{};
            std::string unit;
            require(static_cast<bool>(value >> kib >> unit) && unit == "kB",
                    "invalid process memory observation");
            return machine::checked_size_mul(kib, 1024);
        };
        if (line.starts_with("VmRSS:")) result.resident = read();
        if (line.starts_with("VmHWM:")) result.peak = read();
    }
    return result;
}
void phase(std::string_view name, Clock::time_point start) {
    const auto seconds = std::chrono::duration<double>(Clock::now() - start).count();
    const auto memory = process_memory();
    std::cout << "phase\t" << name << '\t' << std::setprecision(9) << seconds << '\t';
    if (memory.resident) std::cout << *memory.resident; else std::cout << "unavailable";
    std::cout << '\t';
    if (memory.peak) std::cout << *memory.peak; else std::cout << "unavailable";
    std::cout << std::endl;
}
using Counts = std::map<std::string, std::size_t>;
void add(Counts &counts, const std::string &key, std::size_t value) {
    counts[key] = machine::checked_size_add(counts[key], value);
}
template<class T> void capacity(Counts &counts, const std::string &key, const std::vector<T> &values) {
    add(counts, key, machine::checked_size_mul(values.capacity(), sizeof(T)));
}
void tensor_counts(Counts &counts, const std::string &prefix,
                   const numeric::SparseRationalTensor &tensor) {
    add(counts, prefix + ".tensors", 1);
    add(counts, prefix + ".logical_coordinates", tensor.domain().logical_size());
    add(counts, prefix + ".supported_coordinates", tensor.support().count());
    add(counts, prefix + ".nonzero_entries", tensor.values().nonzeros());
    const auto* entries = tensor.values().materialized_entries();
    add(counts, prefix + ".retained_value_fields", tensor.values().retained_source() ? 1 : 0);
    add(counts, prefix + ".materialized_value_fields", entries ? 1 : 0);
    if (entries) capacity(counts, prefix + ".entry_capacity_bytes", *entries);
    else add(counts, prefix + ".entry_capacity_bytes", 0);
    capacity(counts, prefix + ".domain_capacity_bytes", tensor.domain().extents());
    add(counts, prefix + ".limb_capacity_bytes", 0);
    if (entries) for (const auto &entry : *entries)
        add(counts, prefix + ".limb_capacity_bytes", numeric::rational_limb_bytes(entry.value));
    if (const auto *indices = tensor.support().sparse_indices())
        capacity(counts, prefix + ".support_capacity_bytes", *indices);
    else if (const auto *bits = tensor.support().packed_bits())
        capacity(counts, prefix + ".support_capacity_bytes", bits->words());
    else add(counts, prefix + ".support_capacity_bytes", 0);
}
Counts endpoint_counts(const storage::TemporalEndpoint &endpoint) {
    Counts counts;
    add(counts, "source.grades", endpoint.source.grade_count());
    for (std::size_t k = 0; k < endpoint.source.grade_count(); ++k) {
        const auto &grade = endpoint.source.grade({static_cast<std::uint16_t>(k)});
        add(counts, "source.cells", grade.cell_count());
        add(counts, "source.occurrences", grade.incidence_count());
        capacity(counts, "source.identity_capacity_bytes", grade.identities);
        capacity(counts, "source.identity_capacity_bytes", grade.relation_identities);
        capacity(counts, "source.identity_capacity_bytes", grade.occurrence_identities);
        capacity(counts, "source.type_capacity_bytes", grade.relation_types);
        capacity(counts, "source.incidence_capacity_bytes", grade.offsets);
        capacity(counts, "source.incidence_capacity_bytes", grade.lower_index);
        capacity(counts, "source.incidence_capacity_bytes", grade.upper_index);
        capacity(counts, "source.reverse_index_capacity_bytes", grade.lower_offsets);
        capacity(counts, "source.reverse_index_capacity_bytes", grade.incidence_by_lower);
        capacity(counts, "source.scale_object_capacity_bytes", grade.cell_scales);
        capacity(counts, "source.lattice_object_capacity_bytes", grade.lattice_coefficient);
        for (const auto &value : grade.cell_scales)
            add(counts, "source.scale_limb_capacity_bytes", numeric::rational_limb_bytes(value));
        for (const auto &value : grade.lattice_coefficient)
            add(counts, "source.lattice_limb_capacity_bytes", numeric::integer_limb_bytes(value));
        for (std::size_t a = 0; a < grade.aspects.aspect_count(); ++a)
            capacity(counts, "source.aspect_capacity_bytes",
                     grade.aspects.plane(static_cast<composite::AspectId>(a)).words());
    }
    for (const auto &tensor : endpoint.tensors) tensor_counts(counts, "fields", tensor.tensor);
    if (endpoint.declarations) {
        capacity(counts, "declarations.owner_capacity_bytes", endpoint.declarations->values);
        std::vector<const composite::DeclaredValue *> pending;
        for (const auto &value : endpoint.declarations->values) pending.push_back(&value.value);
        while (!pending.empty()) {
            const auto *value = pending.back();
            pending.pop_back();
            add(counts, "declarations.nodes", 1);
            std::visit([&](const auto &item) {
                using T = std::decay_t<decltype(item)>;
                if constexpr (std::is_same_v<T, numeric::Integer>)
                    add(counts, "declarations.integer_limb_capacity_bytes", numeric::integer_limb_bytes(item));
                else if constexpr (std::is_same_v<T, numeric::Rational>)
                    add(counts, "declarations.rational_limb_capacity_bytes", numeric::rational_limb_bytes(item));
                else if constexpr (std::is_same_v<T, numeric::SparseRationalTensor>)
                    tensor_counts(counts, "declaration_tensors", item);
                else if constexpr (std::is_same_v<T, composite::DeclaredValue::Alternative> ||
                                   std::is_same_v<T, composite::DeclaredValue::Optional>) {
                    capacity(counts, "declarations.node_capacity_bytes", item.element);
                    for (const auto &child : item.element) pending.push_back(&child);
                } else if constexpr (std::is_same_v<T, composite::DeclaredValue::Product> ||
                                     std::is_same_v<T, composite::DeclaredValue::Repeated>) {
                    capacity(counts, "declarations.node_capacity_bytes", item.elements);
                    for (const auto &child : item.elements) pending.push_back(&child);
                }
            }, value->data);
        }
    }
    return counts;
}
std::uint64_t cutoff_value(const std::string &text) {
    std::uint64_t value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    require(!text.empty() && error == std::errc{} && end == text.data() + text.size(),
            "cutoff must be an unsigned publication sequence");
    return value;
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path input, output;
        std::optional<std::uint64_t> cutoff;
        bool record = false, accessions = false;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            require(seen.insert(option).second, "duplicate profile option");
            if (option == "--record") record = true;
            else if (option == "--accessions") accessions = true;
            else {
                require(option == "--input" || option == "--rewrite" || option == "--cutoff",
                        "usage: gadget_native_load_profile --input FILE [--cutoff N] [--record] [--accessions] [--rewrite FILE]");
                require(i + 1 < argc, "missing profile option value");
                const std::string value(argv[++i]);
                if (option == "--input") input = value;
                else if (option == "--rewrite") output = value;
                else cutoff = cutoff_value(value);
            }
        }
        require(!input.empty(), "native profile requires --input");
        auto start = Clock::now();
        auto history = [&] {
            const storage::ReadOnlyMappedFile image(input);
            std::cout << "value\tinput.bytes\t" << image.size() << std::endl;
            return storage::decode_native_series(image.bytes(), gadget::native_data_limits());
        }();
        phase("decode", start);
        start = Clock::now();
        const auto selected = storage::NativeSeriesPrefix::capture(
            history, cutoff.value_or(history.last_sequence())).history();
        const auto study = gadget::NativeStudy::from_history(selected);
        phase("select", start);
        start = Clock::now();
        Counts counts = endpoint_counts(*study.endpoint());
        add(counts, "history.archive_bytes", history.archive_bytes());
        std::set<const std::vector<std::byte> *> components;
        history.visit_components([&](const auto &part) {
            if (components.insert(part.bytes.get()).second) {
                add(counts, "history.components", 1);
                add(counts, "history.component_bytes", part.bytes->size());
                capacity(counts, "history.component_capacity_bytes", *part.bytes);
            }
        });
        phase("account", start);
        std::optional<gadget::NativeRecord> publication;
        if (record) {
            start = Clock::now();
            publication = gadget::NativeRecord::from_history(selected);
            phase("record", start);
            std::vector<gadget::RawByteSource> owners;
            for (const auto &source : publication->raw_sources()) {
                add(counts, "record.source_references", 1);
                const auto known = std::find_if(owners.begin(), owners.end(),
                    [&](const auto &owner) { return source.shares_owner(owner); });
                if (known == owners.end()) {
                    owners.push_back(source);
                    add(counts, "record.unique_byte_owners", 1);
                    add(counts, "record.raw_owner_bytes", source.bytes().size());
                }
            }
        }
        if (accessions) {
            start = Clock::now();
            const auto directory = publication ? publication->accession_directory()
                                               : study.accession_directory();
            phase("accessions", start);
            const auto stats = directory->stats();
            add(counts, "accessions.indexed_records", stats.indexed_records);
            add(counts, "accessions.binding_bytes", stats.binding_bytes);
            start = Clock::now();
            const auto reused = publication ? publication->accession_directory()
                                            : study.accession_directory();
            require(reused == directory && directory->pin_snapshot() == study.endpoint(),
                    "accession preparation changed the selected owner");
            phase("accessions_reused", start);
        }
        if (!output.empty()) {
            start = Clock::now();
            storage::publish_file_exact(output, storage::FilePublication::CreateNew,
                [&](std::ostream &stream) {
                    storage::stream_native_series(history, [&](auto bytes) {
                        storage::write_stream_exact(stream, bytes);
                    });
                });
            phase("rewrite", start);
            start = Clock::now();
            const storage::ReadOnlyMappedFile original(input), rewritten(output);
            require(original.size() == rewritten.size() &&
                    std::equal(original.bytes().begin(), original.bytes().end(), rewritten.bytes().begin()),
                    "rewritten native history differs from the complete original bytes");
            phase("verify_bytes", start);
            add(counts, "history.verified_bytes", original.size());
        }
        for (const auto &[key, value] : counts) std::cout << "value\t" << key << '\t' << value << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
