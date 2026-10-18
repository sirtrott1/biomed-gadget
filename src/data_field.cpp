#include <gadget/data_field.hpp>
#include <gadget/io.hpp>
#include <gadget/build.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <algorithm>
#include <limits>
#include <set>

namespace gadget {
NativeRecord import_data_field(const RawByteSource &source, const NativeRecord &coordinates,
                               std::uint64_t source_namespace, const RawByteSource *values,
                               const RawByteSource *support) {
    const auto require = [](bool value, const char *message) {
        if (!value) throw std::invalid_argument(message);
    };
    const auto count = [&](const NativeRecord &record) {
        require(record.kind() == NativeRecordKind::Integer, "data coordinate requires an integer");
        const auto value = record.integer();
        require(value >= 0 && value <= rfc::numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())),
                "data coordinate exceeds the host representation");
        return static_cast<std::size_t>(std::stoull(value.get_str()));
    };
    require(coordinates.at("schema").text() == "gadget.native-data.field.v1", "data field declaration schema");
    const std::set<std::string> allowed{"schema", "native_source", "dataset", "field", "field_id", "shape",
        "rows", "columns", "axes", "row_labels", "column_labels", "vertex_map", "row_map", "parent",
        "dependencies", "quantity", "units", "number_rule", "observed_count", "exporter", "identities_native"};
    for (const auto &name : coordinates.names())
        require(allowed.contains(name), "unknown native data coordinate declaration member");
    const auto binding = coordinates.at("native_source");
    const std::string_view raw(reinterpret_cast<const char *>(source.bytes().data()), source.bytes().size());
    require(binding.at("sha256").text() == sha256(raw), "data field native source digest differs");
    require(source_namespace != 0, "data field requires an explicit source namespace");
    const auto shape = coordinates.at("shape");
    require(shape.size() == 2, "data field requires two declared axes");
    const auto rows = count(shape.at(0)), columns = count(shape.at(1));
    require(rows > 0 && columns > 0, "data field axes must be nonempty");
    const rfc::numeric::IndexDomain domain({rows, columns});
    auto series = rfc::storage::decode_native_series(source.bytes(), native_data_limits());
    require(binding.at("namespace").text() == std::to_string(series.source_namespace()),
            "data field native source namespace differs");
    const auto tail = series.tail();
    auto endpoint = tail.endpoint();
    require(endpoint.source.grade(rfc::core::Grade{1}).cell_count() == rows, "data field native row carrier differs");
    require(coordinates.at("rows").size() == rows && coordinates.at("row_labels").size() == rows &&
                coordinates.at("columns").size() == columns && coordinates.at("column_labels").size() == columns,
            "data field original coordinate map differs");
    for (const auto &name : {"rows", "columns"}) {
        const auto values = coordinates.at(name);
        std::set<std::size_t> original;
        for (std::size_t i = 0; i < values.size(); ++i)
            require(original.insert(count(values.at(i))).second, "repeated original data coordinate");
    }
    if (values) {
        const auto decode = [&](const RawByteSource &input, const char *name) {
            const std::string_view raw(reinterpret_cast<const char *>(input.bytes().data()), input.bytes().size());
            require(binding.at(name).text() == sha256(raw), "native data tensor digest differs");
            const auto limits = native_data_limits();
            rfc::storage::temporal_codec_detail::Reader reader(input.bytes(), limits.temporal);
            auto tensor = rfc::storage::temporal_codec_detail::tensor(reader);
            require(reader.remaining() == 0 && tensor.domain().rank() == 1 &&
                        tensor.domain().logical_size() == domain.logical_size() &&
                        tensor.support().count() == domain.logical_size(), "native data tensor domain differs");
            return tensor;
        };
        auto tensor = decode(*values, "values_sha256");
        auto observed = rfc::numeric::StructuralSupport::full(domain.logical_size());
        if (support) {
            const auto indicator = decode(*support, "support_sha256");
            const auto kind = binding.at("support_kind").text();
            require(kind == "observed" || kind == "missing", "native data support indicator kind differs");
            std::vector<std::size_t> active;
            for (const auto &entry : indicator.values().entries())
                require(entry.value == 1, "native support indicator requires exact zero or one");
            if (kind == "observed") {
                active.reserve(indicator.values().entries().size());
                for (const auto &entry : indicator.values().entries()) active.push_back(entry.index);
            } else {
                active.reserve(domain.logical_size() - indicator.values().entries().size());
                auto missing = indicator.values().entries().begin();
                for (std::size_t i = 0; i < domain.logical_size(); ++i) {
                    if (missing != indicator.values().entries().end() && missing->index == i) ++missing;
                    else active.push_back(i);
                }
            }
            observed = rfc::numeric::StructuralSupport::packed(domain.logical_size(), active);
            for (const auto &entry : tensor.values().entries())
                require(observed.contains(entry.index), "native coefficient occurs outside observed support");
        }
        require(std::none_of(endpoint.carriers.begin(), endpoint.carriers.end(), [](const auto &value) { return value.key == 71; }) &&
                    std::none_of(endpoint.tensors.begin(), endpoint.tensors.end(), [](const auto &value) { return value.key == 72; }),
                "native data tensor source already owns the requested field");
        endpoint.carriers.push_back({71, rfc::storage::TemporalCarrierKind::Explicit, {}, {}, 71, tensor.domain()});
        rfc::storage::TemporalTensor field;
        field.key = 72; field.output = 71; field.role = rfc::storage::TemporalTensorRole::Field;
        field.tensor = {tensor.domain(), std::move(observed), tensor.values()};
        endpoint.tensors.push_back(std::move(field));
    } else {
        require(!support, "native support input requires native values");
    }
    auto field = std::find_if(endpoint.tensors.begin(), endpoint.tensors.end(), [](const auto &value) { return value.key == 72; });
    auto carrier = std::find_if(endpoint.carriers.begin(), endpoint.carriers.end(), [](const auto &value) { return value.key == 71; });
    require(field != endpoint.tensors.end() && carrier != endpoint.carriers.end() &&
                field->output == 71 && !field->input && field->tensor.domain().logical_size() == domain.logical_size(),
            "data field native tensor carrier differs");
    require(field->tensor.support().count() == count(coordinates.at("observed_count")),
            "data field observed support differs");
    carrier->domain = domain;
    field->tensor = {domain, field->tensor.support(), field->tensor.values()};
    field->interpretation = "Original selected row and measurement coordinates with declared observed support";
    const auto occupied = NativeRecord::source_keys(endpoint);
    std::uint64_t next_key = 1;
    for (const auto key : occupied) {
        require(key < std::numeric_limits<std::uint64_t>::max(), "data field byte owner identity space exhausted");
        next_key = std::max(next_key, key + 1);
    }
    std::map<std::uint64_t, std::uint64_t> source_map;
    std::vector<RawByteSource> owners;
    for (const auto &bytes : coordinates.raw_sources()) {
        require(next_key < std::numeric_limits<std::uint64_t>::max(), "data coordinate byte owner identity space exhausted");
        source_map.emplace(bytes.key(), next_key);
        owners.push_back(bytes.rekey(next_key++));
    }
    auto metadata = coordinates.unbound_value(source_map);
    auto &members = std::get<NativeRecordValue::Object>(metadata.data).members;
    members.emplace_back("gadget_sha256", NativeRecordValue(std::string(build::source_sha256())));
    members.emplace_back("rexgraph_sha256", NativeRecordValue(std::string(build::rexgraph_sha256())));
    return NativeRecord::create_from_endpoint(metadata, source_namespace,
                                             std::move(owners), std::move(endpoint));
}
} // namespace gadget
