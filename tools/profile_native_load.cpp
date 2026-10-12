#include "profile_process.hpp"
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
using Clock = gadget::profiling::ProcessClock;
using gadget::profiling::process_memory;
using gadget::profiling::phase;
void require(bool value, const char *message) {
    if (!value) throw std::invalid_argument(message);
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
std::size_t verify_encoded(const storage::NativeEncodedTensor &encoded,
                           const storage::PinnedNativeAccession &pin) {
    require(encoded.publication().pin() == pin.endpoint() &&
                encoded.publication().source_namespace() == pin.source_namespace() &&
                encoded.publication().sequence() == pin.sequence(),
            "encoded tensor changed its selected publication");
    const auto tensor = pin.tensor();
    require(encoded.tensor().domain() == tensor->domain() &&
                storage::temporal_endpoint_detail::support_equal(encoded.tensor().support(),
                                                                tensor->support()),
            "encoded tensor changed its domain or support");
    const auto &values = tensor->values();
    auto cursor = values.retained_source() ? values.retained_source()->open() : nullptr;
    const auto *entries = values.materialized_entries();
    require(cursor || entries, "selected tensor has no retained coefficient reader");
    numeric::SparseRationalEntry entry;
    std::size_t count{};
    encoded.tensor().for_each_nonzero([&](std::size_t index, const numeric::Rational &value) {
        if (cursor) {
            require(cursor->next(entry) && entry.index == index && entry.value == value,
                    "encoded tensor coefficient differs from its selected publication");
        } else {
            require(count < entries->size() && entries->at(count).index == index &&
                        entries->at(count).value == value,
                    "encoded tensor coefficient differs from its selected publication");
        }
        count = machine::checked_size_add(count, 1);
    });
    require(count == values.nonzeros() && count == encoded.tensor().stored_nonzeros() &&
                (!cursor || !cursor->next(entry)),
            "encoded tensor omitted selected coefficients");
    return count;
}
} // namespace
int main(int argc, char **argv) {
    try {
        std::filesystem::path input, output;
        std::optional<std::uint64_t> cutoff, encoded_key;
        bool record = false, accessions = false, fields = false;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string option(argv[i]);
            require(seen.insert(option).second, "duplicate profile option");
            if (option == "--record") record = true;
            else if (option == "--accessions") accessions = true;
            else if (option == "--fields") fields = true;
            else {
                require(option == "--input" || option == "--rewrite" || option == "--cutoff" ||
                            option == "--encoded-tensor",
                        "usage: gadget_native_load_profile --input FILE [--cutoff N] [--record] [--accessions] [--fields] [--encoded-tensor KEY] [--rewrite FILE]");
                require(i + 1 < argc, "missing profile option value");
                const std::string value(argv[++i]);
                if (option == "--input") input = value;
                else if (option == "--rewrite") output = value;
                else if (option == "--cutoff") cutoff = cutoff_value(value);
                else encoded_key = cutoff_value(value);
            }
        }
        require(!input.empty(), "native profile requires --input");
        gadget::profiling::cpu_clock();
        const auto total = Clock::now();
        auto start = total;
        Counts counts;
        std::vector<storage::PinnedNativeAccession> retained_fields;
        std::optional<storage::PinnedNativeAccession> encoded_pin;
        std::optional<storage::NativeEncodedTensor> encoded;
        std::weak_ptr<const storage::TemporalEndpoint> selected_owner;
        {
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
        std::cout << "value\tpublication.source_namespace\t" << selected.source_namespace() << '\n';
        std::cout << "value\tpublication.cutoff\t" << study.endpoint()->sequence << '\n';
        start = Clock::now();
        counts = endpoint_counts(*study.endpoint());
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
        if (accessions || fields || encoded_key) {
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
            if (fields) {
                start = Clock::now();
                const auto endpoint = study.endpoint();
                selected_owner = endpoint;
                for (const auto *name : {"logical_coordinates", "supported_coordinates", "nonzero_entries"})
                    add(counts, std::string("field_access.") + name, 0);
                retained_fields.reserve(endpoint->tensors.size());
                for (const auto &field : endpoint->tensors) {
                    storage::NativeAccessionRequest request;
                    request.kind = storage::NativeAccessionKind::Tensor;
                    request.subject = field.key;
                    auto pin = directory->pin(request);
                    require(pin.endpoint() == endpoint && pin.tensor().get() == &field.tensor,
                            "field accession copied or changed the selected tensor");
                    add(counts, "field_access.logical_coordinates", field.tensor.domain().logical_size());
                    add(counts, "field_access.supported_coordinates", field.tensor.support().count());
                    add(counts, "field_access.nonzero_entries", field.tensor.stored_nonzeros());
                    retained_fields.push_back(std::move(pin));
                }
                add(counts, "field_access.pins", retained_fields.size());
                capacity(counts, "field_access.pin_capacity_bytes", retained_fields);
                phase("field_access", start);
                start = Clock::now();
                for (std::size_t i = 0; i < retained_fields.size(); ++i) {
                    const auto &pin = retained_fields[i];
                    require(pin.endpoint() == endpoint &&
                                pin.tensor().get() == &endpoint->tensors[i].tensor &&
                                pin.source_namespace() == selected.source_namespace() &&
                                pin.sequence() == endpoint->sequence,
                            "retained field pin changed source, cutoff or tensor ownership");
                }
                phase("fields_reused", start);
            }
            if (encoded_key) {
                start = Clock::now();
                encoded = storage::NativeEncodedTensor::from_series(selected, *encoded_key);
                phase("encoded_select", start);
                add(counts, "encoded.available", encoded ? 1 : 0);
                std::cout << "value\tencoded.key\t" << *encoded_key << '\n';
                if (encoded) {
                    storage::NativeAccessionRequest request;
                    request.kind = storage::NativeAccessionKind::Tensor;
                    request.subject = *encoded_key;
                    encoded_pin = directory->pin(request);
                    selected_owner = encoded_pin->endpoint();
                    start = Clock::now();
                    add(counts, "encoded.verified_nonzeros", verify_encoded(*encoded, *encoded_pin));
                    add(counts, "encoded.supported_coordinates", encoded->tensor().support().count());
                    std::cout << "value\tencoded.encoding_sequence\t"
                              << encoded->encoding_sequence() << '\n';
                    phase("encoded_verify", start);
                }
            }
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
        start = Clock::now();
        }
        phase("release", start);
        if (fields || encoded) {
            start = Clock::now();
            {
                const auto endpoint = selected_owner.lock();
                require(endpoint || (retained_fields.empty() && !encoded),
                        "retained field access lost its publication owner");
                add(counts, "field_access.materialized_value_fields_after", 0);
                for (std::size_t i = 0; i < retained_fields.size(); ++i) {
                    require(retained_fields[i].endpoint() == endpoint &&
                                retained_fields[i].tensor().get() == &endpoint->tensors[i].tensor,
                            "field pin did not survive intake owner destruction");
                    add(counts, "field_access.materialized_value_fields_after",
                        retained_fields[i].tensor()->values().materialized_entries() ? 1 : 0);
                }
                add(counts, "field_access.lifetime_verified_pins", retained_fields.size());
                if (encoded) {
                    add(counts, "encoded.lifetime_verified_nonzeros",
                        verify_encoded(*encoded, *encoded_pin));
                    add(counts, "encoded.materialized_field_after",
                        encoded_pin->tensor()->values().materialized_entries() ? 1 : 0);
                }
            }
            phase("fields_lifetime", start);
            start = Clock::now();
            std::vector<storage::PinnedNativeAccession>().swap(retained_fields);
            encoded.reset();
            encoded_pin.reset();
            require(selected_owner.expired(), "released field pins retained their publication");
            phase("release_fields", start);
        }
        for (const auto &[key, value] : counts) std::cout << "value\t" << key << '\t' << value << '\n';
        phase("total", total);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
