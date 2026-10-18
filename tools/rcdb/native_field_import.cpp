#include <gadget/data_field.hpp>
#include <gadget/record_view.hpp>
#include <gadget/process.hpp>
#include <rfc/storage/file_io.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/storage/rfcx.hpp>
#include <charconv>
#include <iostream>
#include <map>
#include <algorithm>

namespace {
gadget::RawByteSource assemble_source(const std::filesystem::path &manifest,
                                     gadget::Json &coordinates, std::size_t max_bytes) {
    using namespace rfc;
    const auto require = [](bool value, const char *message) {
        if (!value) throw std::invalid_argument(message);
    };
    const auto count = [&](const gadget::Json &value) {
        const auto n = value.integer();
        require(n >= 0, "native source component count must be nonnegative");
        return static_cast<std::size_t>(n);
    };
    const auto text = gadget::read_text(manifest, max_bytes);
    const auto binding = coordinates.at("native_source");
    require(gadget::sha256(text) == binding.at("components_sha256").text(), "native component manifest digest differs");
    const auto description = gadget::Json::parse(text);
    require(description.at("schema").text() == "gadget.native-data.components.v1", "native component manifest schema differs");
    const auto rows = count(coordinates.at("shape").at(0));
    const auto vertices = coordinates.at("vertex_map").size();
    const auto source_namespace = std::stoull(binding.at("namespace").text());
    const auto parts = description.at("parts");
    require(rows && vertices && parts.size(), "native source components must cover nonempty original axes");
    std::optional<core::RelationalComplexBuilder> builder;
    std::vector<std::byte> signature;
    std::optional<core::SemanticScale> scale;
    std::vector<bool> witnessed(vertices);
    std::size_t cursor = 0, next_occurrence = 1;
    auto retained = gadget::Json::array();
    for (std::size_t part_index = 0; part_index < parts.size(); ++part_index) {
        const auto part = parts.at(part_index);
        const auto offset = count(part.at("offset")), length = count(part.at("count"));
        require(offset == cursor && length > 0 && cursor <= rows && length <= rows - cursor,
                "native source components omit, overlap or reorder original rows");
        const std::filesystem::path relative(part.at("file").text());
        require(!relative.empty() && relative == relative.filename() && relative != "." && relative != "..",
                "native source component requires a local file name");
        const auto bytes = storage::read_file_exact(manifest.parent_path() / relative, max_bytes);
        const std::string_view raw(reinterpret_cast<const char *>(bytes.data()), bytes.size());
        require(gadget::sha256(raw) == part.at("sha256").text(), "native source component digest differs");
        const auto series = storage::decode_native_series(bytes, gadget::native_data_limits());
        const auto tail = series.tail();
        const auto &endpoint = tail.endpoint();
        const auto &source = endpoint.source;
        require(series.source_namespace() == source_namespace && source.grade_count() == 2 &&
                    endpoint.tensors.empty() && endpoint.carriers.empty() && endpoint.local_blocks.empty() &&
                    endpoint.finite_systems.empty() && !endpoint.codes && !endpoint.declarations,
                "native source component has a foreign namespace or additional mathematical owners");
        storage::rfcx_detail::Writer encoded_signature;
        storage::rfcx_detail::write_signature(encoded_signature, source.signature());
        if (!builder) {
            signature = encoded_signature.bytes();
            scale = source.semantic_scale();
            builder.emplace(source.signature(), *scale);
            builder->declare_grade(core::Grade{1});
            for (std::size_t i = 0; i < vertices; ++i) builder->declare_witness_identity(i + 1);
        } else {
            require(signature == encoded_signature.bytes() && source.semantic_scale().value == scale->value,
                    "native source component declaration differs");
        }
        const auto &base = source.grade(core::Grade{0});
        const auto &relations = source.grade(core::Grade{1});
        require(relations.cell_count() == length, "native component row carrier differs");
        for (const auto key : base.identities) {
            require(key > 0 && key <= vertices, "native component contains a foreign original vertex");
            witnessed[key - 1] = true;
        }
        const auto existence = source.signature().existence_aspect();
        const auto orientation = source.signature().orientation_aspect();
        for (std::size_t row = 0; row < length; ++row) {
            require(relations.relation_identities.at(row).value == offset + row + 1, "native component row identity differs");
            const auto [begin, end] = relations.incidence_range(row);
            std::vector<core::BoundaryOccurrenceInput> occurrences;
            occurrences.reserve(end - begin);
            for (std::size_t slot = begin; slot < end; ++slot) {
                require(relations.occurrence_identities.at(slot).value == next_occurrence++, "native component occurrence order differs");
                std::vector<core::AspectValue> aspects;
                for (std::size_t i = 0; i < source.signature().aspects().size(); ++i) {
                    const auto aspect = static_cast<composite::AspectId>(i);
                    if (aspect != existence && aspect != orientation)
                        aspects.push_back({aspect, relations.aspects.get(aspect, slot)});
                }
                auto magnitude = source.incidence_coefficient(core::Grade{1}, slot);
                if (magnitude < 0) magnitude = -magnitude;
                occurrences.emplace_back(base.identities.at(relations.lower_index.at(slot)),
                    relations.aspects.get(existence, slot), relations.aspects.get(orientation, slot),
                    std::move(magnitude), std::move(aspects), relations.occurrence_identities.at(slot));
            }
            builder->add_grade_one(relations.relation_identities.at(row), relations.relation_types.at(row), occurrences);
        }
        auto reference = gadget::Json::object();
        for (const auto &key : {"offset", "count", "sha256"}) reference.set(key, part.at(key));
        retained.push(reference);
        cursor += length;
    }
    require(cursor == rows && std::all_of(witnessed.begin(), witnessed.end(), [](bool value) { return value; }),
            "native source components leave original cells uncovered");
    storage::NativeTemporalSeries series(storage::TemporalEndpoint(builder->freeze()), source_namespace, gadget::native_data_limits());
    auto encoded = storage::encode_native_series(series);
    require(encoded.size() <= max_bytes, "assembled native source byte admission exceeded");
    coordinates.at("native_source").set("sha256", gadget::Json::string(gadget::sha256(
        std::string_view(reinterpret_cast<const char *>(encoded.data()), encoded.size()))));
    coordinates.at("native_source").set("components", retained);
    return gadget::RawByteSource::own(1, std::move(encoded));
}
}

int main(int argc, char **argv) {
    try {
        std::map<std::string, std::string> options;
        for (int i = 1; i < argc; i += 2) {
            const std::string key(argv[i]);
            if ((key != "--input" && key != "--coordinates" && key != "--namespace" && key != "--out" &&
                 key != "--max-bytes" && key != "--values" && key != "--support" && key != "--components") ||
                i + 1 == argc || !options.emplace(key, argv[i + 1]).second)
                throw std::invalid_argument("expected input, coordinates, namespace and fresh output");
        }
        for (const auto &key : {"--coordinates", "--namespace", "--out"})
            if (!options.contains(key)) throw std::invalid_argument("native field import options missing");
        if (options.contains("--input") == options.contains("--components"))
            throw std::invalid_argument("native field import requires one structural source input");
        std::uint64_t max_bytes = 64 * 1024 * 1024;
        if (options.contains("--max-bytes")) {
            const auto &text = options.at("--max-bytes");
            const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), max_bytes);
            if (error != std::errc{} || end != text.data() + text.size() || !max_bytes)
                throw std::invalid_argument("native field import byte admission must be positive");
        }
        (void)gadget::limit_process_memory();
        for (const auto &key : {"--input", "--components", "--coordinates"})
            if (options.contains(key) && std::filesystem::file_size(options.at(key)) > max_bytes)
                throw std::invalid_argument("native field import byte admission exceeded");
        std::optional<gadget::RawByteSource> values, support;
        for (const auto &key : {"--values", "--support"})
            if (options.contains(key)) {
                if (std::filesystem::file_size(options.at(key)) > max_bytes)
                    throw std::invalid_argument("native tensor input byte admission exceeded");
                auto raw = gadget::RawByteSource::own(1, rfc::storage::read_file_exact(options.at(key), max_bytes));
                if (std::string_view(key) == "--values") values.emplace(std::move(raw));
                else support.emplace(std::move(raw));
            }
        std::uint64_t source_namespace{};
        const auto &text = options.at("--namespace");
        const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), source_namespace);
        if (error != std::errc{} || end != text.data() + text.size() || !source_namespace)
            throw std::invalid_argument("native field import namespace must be positive");
        auto coordinates = gadget::Json::parse(gadget::read_text(options.at("--coordinates"), max_bytes));
        const auto source = options.contains("--components") ? assemble_source(options.at("--components"), coordinates, max_bytes) :
            gadget::RawByteSource::own(1, rfc::storage::read_file_exact(options.at("--input"), max_bytes));
        const auto control = gadget::admit_record_view(coordinates);
        const auto record = gadget::import_data_field(source, control, source_namespace,
            values ? &*values : nullptr, support ? &*support : nullptr);
        record.write(options.at("--out"));
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
