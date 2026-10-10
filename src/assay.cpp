#include "assay_detail.hpp"
#include "assay_field_detail.hpp"
#include "assay_quality_detail.hpp"
#include "assay_coverage_detail.hpp"
#include "genomic_detail.hpp"
#include <algorithm>
#include <gadget/resources.hpp>
#include <limits>
#include <rfc/core/builder.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/runtime/work_graph.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <rfc/storage/mapped_file.hpp>

namespace gadget {
namespace {
using namespace assay_detail;
using namespace rfc::storage;
using rfc::numeric::SparseRationalEntry;
using rfc::numeric::SparseRationalTensor;
void admit_matrix_spec(const Json &spec) {
    members(spec,
            {"id", "kind", "format", "source", "features", "observations", "source_namespace",
             "carrier", "quantity", "units", "channel", "method", "arithmetic", "implicit_entries",
             "missing_tokens"},
            {"missing", "nan_policy", "relations", "reference", "observation_fields", "provenance",
             "partition_size", "metrics"});
    for (const auto *name : {"id", "kind", "quantity", "units", "method"})
        (void)text(spec.at(name), source_byte_limit());
    require(spec.at("carrier").text() == "features" || spec.at("carrier").text() == "observations",
            "assay carrier axis");
    const auto c = spec.at("channel").text();
    require(c == "G" || c == "E" || c == "T" || c == "observation", "assay channel");
    for (const auto *name : {"source", "features", "observations"})
        members(spec.at(name), {"path", "sha256"});
    if (spec.has("missing"))
        members(spec.at("missing"), {"path", "sha256"});
    if (spec.has("relations"))
        members(spec.at("relations"), {"path", "sha256"});
    if (spec.has("reference"))
        (void)reference_record(spec.at("reference"));
    (void)key(spec.at("source_namespace"));
}
void admit_resource(const Json &resource) {
    members(resource, {"id", "format", "source", "method"},
            {"reference", "weight_column", "reference_source", "provenance"});
    (void)text(resource.at("id"));
    (void)text(resource.at("method"));
    members(resource.at("source"), {"path", "sha256"});
    if (resource.has("provenance")) {
        require(resource.at("provenance").size() > 0, "resource provenance requires artifacts");
        for (std::size_t i = 0; i < resource.at("provenance").size(); ++i)
            members(resource.at("provenance").at(i), {"path", "sha256"});
    }
}
Json string_array(const std::vector<std::string> &values, std::size_t begin, std::size_t end) {
    auto out = Json::array();
    for (auto i = begin; i < end; ++i)
        out.push(Json::string(values[i]));
    return out;
}
struct RelationEvidence {
    std::string bytes;
    std::vector<std::pair<std::size_t, std::size_t>> relations, occurrences, fields;
};
// Render the already parsed declaration with the existing JSON codec. Retain
// byte spans as records are written, avoiding a copy of every relation for
// each of its occurrences. The original input remains in the workspace.
RelationEvidence relation_evidence(const Json &declaration) {
    RelationEvidence out;
    out.bytes = "{";
    bool first_root = true;
    json_object_object_foreach(declaration.raw(), root_name, root_value) {
        (void)root_value;
        if (!first_root)
            out.bytes += ',';
        first_root = false;
        out.bytes += Json::string(root_name).dump_compact() + ":";
        const std::string name(root_name);
        const auto array = declaration.at(name);
        if (name != "relations" && name != "fields") {
            out.bytes += array.dump_compact();
            continue;
        }
        out.bytes += '[';
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (i)
                out.bytes += ',';
            const auto start = out.bytes.size();
            const auto record = array.at(i);
            if (name == "fields") {
                out.bytes += record.dump_compact();
                out.fields.emplace_back(start, out.bytes.size());
                continue;
            }
            out.bytes += '{';
            bool first = true;
            json_object_object_foreach(record.raw(), member_name, member_value) {
                (void)member_value;
                if (!first)
                    out.bytes += ',';
                first = false;
                out.bytes += Json::string(member_name).dump_compact() + ":";
                const auto value = record.at(member_name);
                if (std::string(member_name) != "participants") {
                    out.bytes += value.dump_compact();
                    continue;
                }
                out.bytes += '[';
                for (std::size_t j = 0; j < value.size(); ++j) {
                    if (j)
                        out.bytes += ',';
                    const auto begin = out.bytes.size();
                    out.bytes += value.at(j).dump_compact();
                    out.occurrences.emplace_back(begin, out.bytes.size());
                }
                out.bytes += ']';
            }
            out.bytes += '}';
            out.relations.emplace_back(start, out.bytes.size());
        }
        out.bytes += ']';
    }
    out.bytes += '}';
    return out;
}
// The packet has one actual source carrier. Partitioning is a storage/execution
// policy, not an assertion that per-partition topology equals whole-carrier topology.
Json make_packet(const fs::path &destination, Matrix &matrix, const Json &spec,
                 const Json &design, std::size_t begin, std::size_t end,
                 std::uint64_t namespace_key, const Json &relation_data,
                 std::vector<NativeTemporalSeries> &produced,
                 const AssayImportObserver &observer) {
    const bool by_features = spec.at("carrier").text() == "features";
    const auto &axis = by_features ? matrix.features : matrix.observations;
    const auto &field_axis = by_features ? matrix.observations : matrix.features;
    const auto relations =
        relation_data.has("relations") ? relation_data.at("relations") : Json::array();
    rfc::core::RelationSignature signature;
    std::map<std::string, rfc::core::RelationTypeId> types;
    for (std::size_t i = 0; i < relations.size(); ++i) {
        const auto r = relations.at(i);
        members(r, {"id", "meaning", "participants"}, {"grade", "type"});
        const auto name = r.has("type") ? text(r.at("type"), 256) : "declared.assay.relation";
        if (!types.contains(name))
            types.emplace(name, signature.add_type(name));
    }
    rfc::core::RelationalComplexBuilder builder(signature);
    std::map<std::string, std::uint64_t> cells_by_id;
    for (auto i = begin; i < end; ++i) {
        builder.declare_witness_identity(i + 1);
        cells_by_id.emplace(axis.ids[i], i + 1);
    }
    std::map<std::string, rfc::core::RelationId> relation_ids;
    std::vector<std::size_t> relation_grades;
    std::vector<std::pair<std::size_t, std::uint64_t>> occurrence_grades;
    std::uint64_t occurrence = 1;
    std::size_t previous_grade = 1;
    for (std::size_t i = 0; i < relations.size(); ++i) {
        const auto r = relations.at(i);
        const auto id = text(r.at("id"));
        (void)text(r.at("meaning"));
        require(!relation_ids.contains(id), "duplicate authored relation identity");
        const auto grade = r.has("grade")
                               ? count(r.at("grade"), std::numeric_limits<std::uint16_t>::max(),
                                       "authored relation grade")
                               : 1;
        require(grade >= 1 && grade >= previous_grade,
                "authored relations must be ordered by increasing positive grade");
        previous_grade = grade;
        const auto name = r.has("type") ? r.at("type").text() : "declared.assay.relation";
        require(r.at("participants").size() > 0, "assay relation arity");
        std::vector<rfc::core::BoundaryOccurrenceInput> lower;
        std::vector<rfc::core::HigherIncidenceInput> higher;
        for (std::size_t j = 0; j < r.at("participants").size(); ++j) {
            const auto item = r.at("participants").at(j);
            members(item, {"cell", "exists", "negative", "magnitude"});
            const auto cell = text(item.at("cell"));
            const auto magnitude = exact_rational(text(item.at("magnitude")));
            require(magnitude >= 0,
                    "occurrence magnitude is nonnegative; its negative bit owns orientation");
            const rfc::core::OccurrenceKey key{occurrence};
            occurrence = rfc::machine::checked_size_add(occurrence, 1);
            occurrence_grades.emplace_back(grade, key.value);
            if (grade == 1) {
                require(cells_by_id.contains(cell), "relation cell outside its whole carrier");
                lower.push_back({cells_by_id.at(cell),
                                 item.at("exists").boolean(),
                                 item.at("negative").boolean(),
                                 magnitude,
                                 {},
                                 key});
            } else {
                require(relation_ids.contains(cell) &&
                            static_cast<std::size_t>(relation_ids.at(cell).grade.value) + 1 ==
                                grade,
                        "higher relation requires an authored adjacent lower cell");
                higher.push_back({relation_ids.at(cell).index,
                                  item.at("exists").boolean(),
                                  item.at("negative").boolean(),
                                  magnitude,
                                  {},
                                  key});
            }
        }
        const auto native =
            grade == 1
                ? builder.add_grade_one(rfc::core::RelationKey{1000000 + i}, types.at(name), lower)
                : builder.add_higher({static_cast<std::uint16_t>(grade)},
                                     rfc::core::RelationKey{1000000 + i}, types.at(name), higher);
        relation_ids.emplace(id, native);
        relation_grades.push_back(grade);
    }
    TemporalEndpoint endpoint(builder.freeze());
    const auto n = end - begin;
    endpoint.carriers.push_back(
        {1, TemporalCarrierKind::SourceGrade, {0}, {}, 0, rfc::numeric::IndexDomain::linear(n)});
    const auto axis_lines = lines(axis.bytes);
    std::string carrier_bytes =
        axis.bytes.substr(axis_lines[0].first, axis_lines[0].second - axis_lines[0].first) + "\n";
    for (auto i = begin; i < end; ++i)
        carrier_bytes +=
            axis.bytes.substr(axis.spans[i].first, axis.spans[i].second - axis.spans[i].first) +
            "\n";
    const auto carrier_axis = read_axis(carrier_bytes);
    auto cells = Json::array(), occurrences = Json::array(), fields = Json::array();
    const auto spec_bytes = spec.dump();
    const auto spec_evidence = span("assay_declaration", spec_bytes, 0, spec_bytes.size());
    const auto retained_relations = relation_evidence(relation_data);
    const auto &relations_bytes = retained_relations.bytes;
    require(retained_relations.occurrences.size() == occurrence_grades.size(),
            "complete relation occurrence evidence");
    for (std::size_t i = 0; i < n; ++i) {
        auto item = Json::object();
        item.set("grade", Json::integer(0));
        item.set("key", Json::string(std::to_string(begin + i + 1)));
        item.set("id", Json::string(axis.ids[begin + i]));
        item.set("label", Json::string(axis.ids[begin + i]));
        item.set("kind", Json::string(by_features ? "assay_feature" : "assay_observation"));
        item.set("evidence", span("carrier_axis", carrier_bytes, carrier_axis.spans[i].first,
                                  carrier_axis.spans[i].second));
        cells.push(item);
    }
    if (!relations.size())
        require(endpoint.source.grade_count() == 1, "unexpected assay source grade");
    else {
        for (std::size_t i = 0; i < relations.size(); ++i) {
            auto item = Json::object();
            item.set("grade", cardinality(relation_grades[i]));
            item.set("key", Json::string(std::to_string(1000000 + i)));
            item.set("id", Json::string("relation:" + relations.at(i).at("id").text()));
            item.set("label", relations.at(i).at("meaning"));
            item.set("kind", Json::string("declared_assay_relation"));
            const auto [start, stop] = retained_relations.relations.at(i);
            item.set("evidence", span("relations", relations_bytes, start, stop));
            cells.push(item);
        }
        for (std::size_t i = 0; i < occurrence_grades.size(); ++i) {
            const auto [grade, k] = occurrence_grades[i];
            auto item = Json::object();
            item.set("grade", cardinality(grade));
            item.set("key", Json::string(std::to_string(k)));
            const auto [start, stop] = retained_relations.occurrences.at(i);
            item.set("evidence", span("relations", relations_bytes, start, stop));
            occurrences.push(item);
        }
    }
    for (std::size_t f = 0; f < field_axis.ids.size(); ++f) {
        TemporalTensor tensor;
        tensor.key = 100 + f;
        tensor.role = TemporalTensorRole::Field;
        tensor.output = 1;
        tensor.interpretation = field_axis.ids[f];
        // These columns are private unpublished owners. A complete identity
        // carrier can transfer them; slicing and support normalization still
        // use the existing construction below.
        bool transfer = false;
        if (by_features && begin == 0 && end == axis.ids.size()) {
            const auto &support = matrix.columns[f].support();
            transfer = support.storage() == rfc::numeric::StructuralSupport::Storage::Full ||
                       (support.storage() == rfc::numeric::StructuralSupport::Storage::SparseIndices &&
                        support.count() < n);
        }
        if (transfer) {
            tensor.tensor = std::move(matrix.columns[f]);
        } else {
            std::vector<std::size_t> support;
            std::vector<SparseRationalEntry> entries;
            bool full = false;
            if (by_features) {
                const auto &column = matrix.columns[f];
                full = column.support().storage() == rfc::numeric::StructuralSupport::Storage::Full;
                if (!full) {
                    for (auto i = begin; i < end; ++i)
                        if (column.support().contains(i))
                            support.push_back(i - begin);
                    full = support.size() == n;
                }
                const auto &values = column.values().entries();
                auto it = std::lower_bound(
                    values.begin(), values.end(), begin,
                    [](const auto &entry, std::size_t index) { return entry.index < index; });
                for (; it != values.end() && it->index < end; ++it)
                    entries.push_back({it->index - begin, it->value});
            } else {
                for (auto i = begin; i < end; ++i) {
                    const auto state = matrix.columns[i].state(f);
                    if (state) {
                        support.push_back(i - begin);
                        if (*state != 0)
                            entries.push_back({i - begin, *state});
                    }
                }
                full = support.size() == n;
            }
            tensor.tensor = {rfc::numeric::IndexDomain::linear(n),
                             full ? rfc::numeric::StructuralSupport::full(n)
                                  : rfc::numeric::StructuralSupport::sparse(n, support),
                             rfc::numeric::SparseRationalVector(n, std::move(entries))};
        }
        endpoint.tensors.push_back(std::move(tensor));
        auto field = Json::object();
        field.set("id", Json::string("all." + field_axis.ids[f]));
        field.set("key", Json::string(std::to_string(100 + f)));
        field.set("context", Json::string("all"));
        for (const auto *name : {"channel", "quantity", "units", "method"})
            field.set(name, spec.at(name));
        field.set("value_type", Json::string("quantitative"));
        field.set("evidence", spec_evidence);
        fields.push(field);
    }
    if (spec.has("observation_fields")) {
        require(!by_features, "observation covariates require the observation carrier");
        const auto declared = spec.at("observation_fields");
        std::set<std::string> field_ids(field_axis.ids.begin(), field_axis.ids.end());
        for (std::size_t f = 0; f < declared.size(); ++f) {
            const auto cov = declared.at(f);
            members(cov, {"id", "column", "kind", "quantity", "units", "method", "missing_tokens"},
                    {"level"});
            const auto id = text(cov.at("id")), column = text(cov.at("column")),
                       kind = text(cov.at("kind"));
            require(field_ids.insert(id).second &&
                        (kind == "exact_value" || kind == "level_indicator" ||
                         kind == "categorical_code"),
                    "covariate identity/type");
            require((kind == "level_indicator") == cov.has("level"),
                    "indicator requires an explicit level");
            if (kind == "level_indicator")
                require(cov.at("units").text() == "dimensionless", "indicator covariate units");
            std::set<std::string> missing;
            for (std::size_t i = 0; i < cov.at("missing_tokens").size(); ++i)
                missing.insert(cov.at("missing_tokens").at(i).text());
            std::vector<std::size_t> support;
            std::vector<SparseRationalEntry> values;
            for (auto i = begin; i < end; ++i) {
                const auto record = axis.records.at(i);
                require(record.has(column), "covariate column is absent from observation metadata");
                const auto token = record.at(column).text();
                if (missing.contains(token))
                    continue;
                const auto value = kind == "level_indicator"
                                       ? Rational(token == cov.at("level").text() ? 1 : 0)
                                       : exact_rational(token);
                support.push_back(i - begin);
                if (value != 0)
                    values.push_back({i - begin, value});
            }
            TemporalTensor tensor;
            tensor.key = 100 + field_axis.ids.size() + f;
            tensor.role = TemporalTensorRole::Field;
            tensor.output = 1;
            tensor.interpretation = id;
            tensor.tensor = {rfc::numeric::IndexDomain::linear(n),
                             support.size() == n
                                 ? rfc::numeric::StructuralSupport::full(n)
                                 : rfc::numeric::StructuralSupport::sparse(n, support),
                             rfc::numeric::SparseRationalVector(n, std::move(values))};
            endpoint.tensors.push_back(std::move(tensor));
            auto field = Json::object();
            field.set("id", Json::string("all." + id));
            field.set("key", Json::string(std::to_string(100 + field_axis.ids.size() + f)));
            field.set("context", Json::string("all"));
            for (const auto *name : {"quantity", "units", "method"})
                field.set(name, cov.at(name));
            field.set("channel", Json::string("covariate"));
            field.set(
                "value_type",
                Json::string(kind == "level_indicator"
                                 ? "indicator"
                                 : (kind == "categorical_code" ? "categorical" : "quantitative")));
            field.set("evidence", spec_evidence);
            fields.push(field);
        }
    }
    if (relation_data.has("fields")) {
        const auto declared = relation_data.at("fields");
        std::map<std::uint16_t, std::uint64_t> carriers;
        for (std::size_t f = 0; f < declared.size(); ++f) {
            const auto field = declared.at(f);
            members(field, {"id", "grade", "channel", "quantity", "units", "method", "values"});
            const auto g = count(field.at("grade"), endpoint.source.grade_count() - 1,
                                 "relation measurement grade");
            require(g > 0, "relation fields require a positive authored grade");
            const rfc::core::Grade grade{static_cast<std::uint16_t>(g)};
            const auto size = endpoint.source.grade(grade).cell_count();
            if (!carriers.contains(grade.value)) {
                carriers.emplace(grade.value, g + 1);
                endpoint.carriers.push_back({g + 1,
                                             TemporalCarrierKind::SourceGrade,
                                             grade,
                                             {},
                                             0,
                                             rfc::numeric::IndexDomain::linear(size)});
            }
            std::set<std::size_t> reported;
            std::vector<std::size_t> support;
            std::vector<SparseRationalEntry> entries;
            require(field.at("values").size() <= size, "relation field coordinate budget");
            for (std::size_t j = 0; j < field.at("values").size(); ++j) {
                const auto v = field.at("values").at(j);
                members(v, {"relation", "value"});
                const auto id = text(v.at("relation"));
                require(relation_ids.contains(id) && relation_ids.at(id).grade == grade,
                        "relation measurement identity/grade");
                const auto index = relation_ids.at(id).index;
                require(reported.insert(index).second, "duplicate relation measurement");
                support.push_back(index);
                const auto value = exact_rational(text(v.at("value")));
                if (is_count_quantity(field.at("quantity").text()))
                    require(value >= 0 && value.get_den() == 1,
                            "relation counts must be nonnegative integers");
                if (value != 0)
                    entries.push_back({index, value});
            }
            TemporalTensor tensor;
            const auto ordinary_fields = rfc::machine::checked_size_add(
                field_axis.ids.size(),
                spec.has("observation_fields") ? spec.at("observation_fields").size() : 0);
            const auto relation_field_base =
                std::max<std::size_t>(100000, rfc::machine::checked_size_add(100, ordinary_fields));
            tensor.key = rfc::machine::checked_size_add(relation_field_base, f);
            tensor.role = TemporalTensorRole::Field;
            tensor.output = carriers.at(grade.value);
            tensor.interpretation = text(field.at("id"));
            tensor.tensor = {rfc::numeric::IndexDomain::linear(size),
                             support.size() == size
                                 ? rfc::numeric::StructuralSupport::full(size)
                                 : rfc::numeric::StructuralSupport::sparse(size, support),
                             rfc::numeric::SparseRationalVector(size, std::move(entries))};
            endpoint.tensors.push_back(std::move(tensor));
            auto meaning = Json::object();
            meaning.set("id", Json::string("all." + field.at("id").text()));
            meaning.set("key", Json::string(std::to_string(endpoint.tensors.back().key)));
            meaning.set("context", Json::string("all"));
            for (const auto *name : {"channel", "quantity", "units", "method"})
                meaning.set(name, field.at(name));
            meaning.set("value_type", Json::string("quantitative"));
            const auto [start, stop] = retained_relations.fields.at(f);
            meaning.set("evidence", span("relations", relations_bytes, start, stop));
            fields.push(meaning);
        }
    }
    if (observer)
        observer("source_construction");
    fs::create_directories(destination);
    {
        const NativeTemporalSeries series(endpoint, namespace_key,
                                          source_codec_limits(source_byte_limit()));
        // The immutable publication owns the complete admitted coefficient
        // lanes. Release the producer's private expanded field storage.
        std::vector<TemporalTensor>().swap(endpoint.tensors);
        publish_file_exact(destination / "source.rfcxs", FilePublication::CreateNew,
                           [&](std::ostream &output) {
                               stream_native_series(series, [&](auto chunk) {
                                   write_stream_exact(output, chunk);
                               });
                           });
        produced.push_back(series);
    }
    if (observer)
        observer("series_publication");
    auto declaration = Json::object();
    declaration.set("schema", Json::string("gadget.scientific.study.v1"));
    declaration.set("id", spec.at("id"));
    auto record = Json::object();
    record.set("path", Json::string("source.rfcxs"));
    record.set("sha256", Json::string("sha256:" + sha256_file(destination / "source.rfcxs")));
    record.set("source_namespace", Json::string(std::to_string(namespace_key)));
    declaration.set("native", record);
    declaration.set("design", design);
    declaration.set("cells", cells);
    declaration.set("occurrences", occurrences);
    declaration.set("fields", fields);
    declaration.set("field_sets", Json::array());
    auto contexts = Json::array(), context = Json::object(), metrics = Json::array();
    context.set("id", Json::string("all"));
    context.set("block", Json{});
    context.set("time", Json{});
    if (spec.has("metrics")) {
        require(spec.at("metrics").size() == endpoint.source.grade_count(),
                "one assay metric per native grade required");
        for (std::size_t g = 0; g < endpoint.source.grade_count(); ++g) {
            auto metric = Json::parse(spec.at("metrics").at(g).dump());
            if (metric.at("kind").text() == "diagonal" && g == 0) {
                members(metric, {"kind", "values"});
                const auto full = metric.at("values");
                require(full.size() == axis.ids.size(),
                        "declared assay metric must cover the entire source carrier");
                auto local = Json::array();
                for (auto i = begin; i < end; ++i)
                    local.push(full.at(i));
                metric.set("values", local);
            }
            metrics.push(metric);
        }
    } else {
        for (std::size_t g = 0; g < endpoint.source.grade_count(); ++g)
            metrics.push(Json::parse("{\"kind\":\"identity\"}"));
    }
    context.set("metrics", metrics);
    contexts.push(context);
    declaration.set("contexts", contexts);
    auto documents = Json::array();
    const auto write_doc = [&](const std::string &id, const std::string &path,
                               const std::string &bytes) {
        write_text(destination / path, bytes);
        auto doc = file_record(path, bytes);
        doc.set("id", Json::string(id));
        doc.set("format", Json::string("utf8"));
        documents.push(doc);
    };
    fs::create_directories(destination);
    write_doc("carrier_axis", "carrier.tsv", carrier_bytes);
    write_doc("assay_declaration", "assay.json", spec_bytes);
    write_doc("relations", "relations.json", relations_bytes);
    declaration.set("documents", documents);
    const auto declaration_bytes = declaration.dump() + "\n";
    write_text(destination / "study.json", declaration_bytes);
    auto packet = file_record("study.json", declaration_bytes);
    packet.set("source_namespace", Json::string(std::to_string(namespace_key)));
    packet.set("begin", cardinality(begin));
    packet.set("end", cardinality(end));
    packet.set("carrier_ids", string_array(axis.ids, begin, end));
    if (observer)
        observer("study_declaration");
    return packet;
}
Json unwrap(const Json &manifest) {
    auto out = Json::object();
    json_object_object_foreach(manifest.raw(), name, value) {
        (void)value;
        if (std::string(name) != "manifest_sha256")
            out.set(name, manifest.at(name));
    }
    return out;
}
void native_members(const NativeRecord &object, const std::set<std::string> &required,
                    const std::set<std::string> &optional = {}) {
    require(object.kind() == NativeRecordKind::Object, "expected object");
    for (const auto &name : required)
        require(object.has(name), "missing member: " + name);
    for (const auto &name : object.names())
        require(required.contains(name) || optional.contains(name), "unknown member: " + name);
}
NativeRecord admit_native_workspace(const fs::path &root,const Json &manifest,
                                   const std::map<std::string,std::vector<std::shared_ptr<const ScientificStudy>>> &supplied={}) {
    using N=NativeRecordValue;
    std::vector<RawByteSource> raw;N::Array originals,publications;
    std::set<std::string> paths;
    const auto capture=[&](const Json &record) {
        const auto path=record.at("path").text();if(!paths.insert(path).second)return;
        const auto data=local_artifact(fs::canonical(root),record);const auto key=raw.size()+1;
        raw.emplace_back(key,std::span<const std::byte>{reinterpret_cast<const std::byte *>(data.data()),data.size()});
        originals.elements.emplace_back(N::Object{{{"path",N(path)},{"bytes",N(N::ByteSpan{key,0,data.size(),false})}}});
    };
    for(const auto *name:{"import_spec","reimport","bindings"})if(manifest.has(name))capture(manifest.at(name));
    const auto assays=manifest.at("assays");
    std::map<fs::path,std::shared_ptr<const ScientificStudy>> owners;
    std::map<const ScientificStudy *,std::pair<std::uint64_t,std::size_t>> publication_bytes;
    for(std::size_t i=0;i<assays.size();++i) {
        const auto assay=assays.at(i);const auto id=assay.at("declaration").at("id").text();
        const auto inputs=assay.at("originals");
        json_object_object_foreach(inputs.raw(),name,value){(void)value;capture(inputs.at(name));}
        const auto partitions=assay.at("partitions");
        for(std::size_t p=0;p<partitions.size();++p) {
            const auto record=partitions.at(p);const auto path=record.at("path").text();
            std::shared_ptr<const ScientificStudy> owner;
            if(!supplied.empty())owner=supplied.at(id).at(p);
            else {
                const auto actual=fs::canonical(root/path);
                if(!owners.contains(actual))owners.emplace(actual,std::make_shared<const ScientificStudy>(actual));
                owner=owners.at(actual);
            }
            require(owner->digest()==record.at("source_digest").text(),"native workspace source binding");
            if(!publication_bytes.contains(owner.get())) {
                const auto key=raw.size()+1;auto bytes=owner->native_publication().archive_source(key);const auto size=bytes.bytes().size();
                raw.push_back(std::move(bytes));publication_bytes.emplace(owner.get(),std::pair{key,size});
            }
            const auto [key,size]=publication_bytes.at(owner.get());
            publications.elements.emplace_back(N::Object{{{"assay",N(id)},{"partition",N(rfc::numeric::Integer(std::to_string(p)))},
                {"path",N(path)},{"bytes",N(N::ByteSpan{key,0,size,false})}}});
        }
    }
    const auto resources=manifest.at("resources");
    for(std::size_t i=0;i<resources.size();++i) {
        const auto resource=resources.at(i);capture(resource.at("original"));
        if(resource.has("reference_original"))capture(resource.at("reference_original"));
        if(resource.has("provenance"))
            for(std::size_t p=0;p<resource.at("provenance").size();++p)
                capture(resource.at("provenance").at(p));
    }
    return NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.assay.workspace.native.v1"))},
        {"manifest",import_record_view(manifest)},{"originals",N(std::move(originals))},
        {"publications",N(std::move(publications))}}}),new_native_record_namespace(),std::move(raw));
}
NativeRecord publish_native_workspace(const fs::path &root,const Json &manifest,
                              const std::map<std::string,std::vector<std::shared_ptr<const ScientificStudy>>> &sources={}) {
    auto publication = admit_native_workspace(root,manifest,sources);
    publication.write(root/"workspace.rfcxs");
    return publication;
}
} // namespace
void import_assays(const fs::path &specification, const fs::path &destination) {
    import_assays(specification, destination, {});
}
void import_assays(const fs::path &specification, const fs::path &destination,
                   const AssayImportObserver &observer) {
    (void)import_assay_workspace(specification, destination, observer);
}
std::shared_ptr<const AssayWorkspace> import_assay_workspace(
    const fs::path &specification, const fs::path &destination,
    const AssayImportObserver &observer) {
    using namespace assay_detail;
    require(fs::is_regular_file(fs::symlink_status(specification)),
            "assay specification must be regular");
    require(!fs::exists(destination), "assay destination already exists");
    const auto path = fs::canonical(specification), root = path.parent_path();
    const auto bytes = read_text(path);
    const auto raw_spec = std::as_bytes(std::span(bytes.data(), bytes.size()));
    const bool native_spec = detect_native_data(raw_spec) == NativeDataFormat::TemporalSeries;
    const auto spec = static_cast<Json>(decode_record(raw_spec));
    const auto import_spec_path = native_spec ? "import-spec.rfcxs" : "import-spec.json";
    members(spec, {"schema", "id", "design", "assays", "resources", "links"});
    require(spec.at("schema").text() == "gadget.assays.import.v1", "assay import schema");
    (void)text(spec.at("id"));
    auto design = Json::parse(spec.at("design").dump());
    members(design, {"question", "population", "selection", "limitations"});
    design.set("covariates", Json::array());
    const auto assays = spec.at("assays"), resources = spec.at("resources"),
               links = spec.at("links");
    require(rfc::machine::checked_size_add(assays.size(), resources.size()) > 0,
            "workspace requires assays or resources");
    std::set<std::string> assay_ids, resource_ids;
    std::set<std::uint64_t> namespaces;
    std::map<std::string, std::set<std::string>> observations;
    auto admitted_assays = Json::array(), admitted_resources = Json::array();
    AssayWorkspace::AdmittedSources admitted_sources;
    std::map<std::string,std::vector<NativeTemporalSeries>> produced_histories;
    fs::create_directories(fs::absolute(destination).parent_path());
    require(fs::create_directory(destination), "workspace creation conflict");
    try {
        write_text(destination / import_spec_path, bytes);
        if (observer)
            observer("specification");
        for (std::size_t a = 0; a < assays.size(); ++a) {
            auto item = assays.at(a);
            admit_matrix_spec(item);
            const auto id = text(item.at("id"));
            require(assay_ids.insert(id).second, "duplicate assay identity");
            auto matrix = read_matrix(root, item);
            if (observer)
                observer("matrix_read");
            observations.emplace(id, std::set<std::string>(matrix.observations.ids.begin(),
                                                           matrix.observations.ids.end()));
            const auto dir = "assays/" + std::to_string(a);
            fs::create_directories(destination / dir);
            auto originals = Json::object();
            const auto retain = [&](const std::string &name, const std::string &content) {
                const auto rel = dir + "/source-" + name;
                write_text(destination / rel, content);
                originals.set(name, file_record(rel, content));
            };
            retain("matrix", matrix.bytes);
            // The complete original is retained on disk; native columns own
            // the parsed values independently of its text image.
            std::string().swap(matrix.bytes);
            retain("features", matrix.features.bytes);
            retain("observations", matrix.observations.bytes);
            if (item.has("provenance")) {
                require(item.at("provenance").size() > 0, "assay provenance artifact budget");
                for (std::size_t i = 0; i < item.at("provenance").size(); ++i)
                    retain("provenance-" + std::to_string(i),
                           local_artifact(root, item.at("provenance").at(i)));
            }
            if (item.has("missing"))
                retain("missing", local_artifact(root, item.at("missing")));
            Json relation_data = Json::object();
            if (item.has("relations")) {
                const auto raw = local_artifact(root, item.at("relations"));
                retain("relations", raw);
                relation_data = Json::parse(raw);
                members(relation_data, {"relations"}, {"fields"});
            }
            const bool by_features = item.at("carrier").text() == "features";
            const auto &carrier = by_features ? matrix.features : matrix.observations;
            const auto requested_chunk =
                item.has("partition_size")
                    ? count(item.at("partition_size"), std::numeric_limits<std::uint32_t>::max(),
                            "packet coordinate count")
                    : std::min<std::size_t>(std::numeric_limits<std::uint32_t>::max(),
                                            carrier.ids.size());
            require(requested_chunk > 0, "packet coordinate count must be positive");
            require(!relation_data.has("relations") || !item.has("partition_size") ||
                        requested_chunk >= carrier.ids.size(),
                    "declared relations require their complete carrier");
            const auto chunk = relation_data.has("relations")
                                   ? carrier.ids.size()
                                   : std::min(requested_chunk, carrier.ids.size());
            if (observer)
                observer("originals");
            auto packets = Json::array();
            const auto base = key(item.at("source_namespace"));
            for (std::size_t begin = 0, index = 0; begin < carrier.ids.size();
                 begin += chunk, ++index) {
                require(index <= std::numeric_limits<std::uint64_t>::max() - base &&
                            namespaces.insert(base + index).second,
                        "native assay namespace collision");
                const auto rel = dir + "/partitions/" + std::to_string(index);
                auto packet = make_packet(destination / rel, matrix, item, design, begin,
                                          std::min(begin + chunk, carrier.ids.size()), base + index,
                                          relation_data,produced_histories[item.at("id").text()],
                                          observer);
                packet.set("path", Json::string(rel + "/study.json"));
                packets.push(packet);
            }
            auto admitted = Json::object();
            admitted.set("declaration", item);
            admitted.set("originals", originals);
            admitted.set("features", matrix.features.records);
            admitted.set("observations", matrix.observations.records);
            admitted.set("partitions", packets);
            admitted_assays.push(admitted);
        }
        for (std::size_t r = 0; r < resources.size(); ++r) {
            const auto resource = resources.at(r);
            admit_resource(resource);
            require(resource_ids.insert(text(resource.at("id"))).second,
                    "duplicate physical resource");
            const auto raw = local_artifact(root, resource.at("source"));
            const auto reference = resource.has("reference_source")
                                       ? local_artifact(root, resource.at("reference_source"))
                                       : std::string{};
            if (!reference.empty())
                require(resource.has("reference") &&
                            "sha256:" + sha256(reference) ==
                                resource.at("reference").at("sequence_sha256").text(),
                        "reference FASTA hash binding");
            const auto inspected = inspect_resource(raw, resource, reference, true);
            const auto rel = "resources/" + std::to_string(r) + "/source";
            fs::create_directories((destination / rel).parent_path());
            write_text(destination / rel, raw);
            auto item = Json::object();
            item.set("declaration", resource);
            item.set("original", file_record(rel, raw));
            if (resource.has("reference_source")) {
                const auto ref_path = "resources/" + std::to_string(r) + "/reference.fa";
                write_text(destination / ref_path, reference);
                item.set("reference_original", file_record(ref_path, reference));
            }
            if (resource.has("provenance")) {
                auto originals = Json::array();
                for (std::size_t p = 0; p < resource.at("provenance").size(); ++p) {
                    const auto content = local_artifact(root, resource.at("provenance").at(p));
                    const auto provenance_path = "resources/" + std::to_string(r) +
                                                 "/provenance-" + std::to_string(p);
                    write_text(destination / provenance_path, content);
                    originals.push(file_record(provenance_path, content));
                }
                item.set("provenance", originals);
            }
            item.set("inspection", inspected);
            admitted_resources.push(item);
        }
        if (observer)
            observer("resources");
        std::set<std::string> link_ids;
        for (std::size_t i = 0; i < links.size(); ++i) {
            const auto link = links.at(i);
            members(link, {"id", "left_assay", "left_observation", "right_assay",
                           "right_observation", "relation", "evidence"});
            require(link_ids.insert(text(link.at("id"))).second, "duplicate assay link");
            for (const auto *side : {"left", "right"}) {
                const auto assay = text(link.at(std::string(side) + "_assay")),
                           obs = text(link.at(std::string(side) + "_observation"));
                require(observations.contains(assay) && observations.at(assay).contains(obs),
                        "link endpoint not admitted in its own assay");
            }
            (void)text(link.at("relation"));
            (void)text(link.at("evidence"), 16384);
        }
        // Producer matrices, endpoints and encoded bytes have left scope. Admit
        // each publication once, then hand its actual owner to final validation.
        for (std::size_t a = 0; a < admitted_assays.size(); ++a) {
            const auto item = admitted_assays.at(a);
            const auto id = item.at("declaration").at("id").text();
            const auto packets = item.at("partitions");
            for (std::size_t p = 0; p < packets.size(); ++p) {
                auto packet = packets.at(p);
                auto owner = ScientificStudy::admit_generated(destination / packet.at("path").text(),
                    produced_histories.at(id).at(p));
                owner->write_native((destination/packet.at("path").text()).parent_path()/"study.rfcxs");
                packet.set("source_digest", Json::string(owner->digest()));
                admitted_sources[id].push_back(std::move(owner));
            }
            produced_histories.erase(id);
        }
        if (observer)
            observer("study_admission");
        auto manifest = Json::object();
        manifest.set("schema", Json::string("gadget.assays.workspace.v1"));
        manifest.set("id", spec.at("id"));
        manifest.set("design", spec.at("design"));
        manifest.set("import_spec", file_record(import_spec_path, bytes));
        // Retain the original declaration and a distinct replay declaration
        // pointing only to artifacts inside this portable workspace.
        auto replay = Json::parse(spec.dump());
        auto replay_assays = Json::array(), replay_resources = Json::array();
        for (std::size_t i = 0; i < admitted_assays.size(); ++i) {
            auto item = Json::parse(assays.at(i).dump());
            const auto originals = admitted_assays.at(i).at("originals");
            item.set("source", originals.at("matrix"));
            for (const auto *name : {"features", "observations", "missing", "relations"})
                if (originals.has(name))
                    item.set(name, originals.at(name));
            if (item.has("provenance")) {
                auto records = Json::array();
                for (std::size_t p = 0; p < item.at("provenance").size(); ++p)
                    records.push(originals.at("provenance-" + std::to_string(p)));
                item.set("provenance", records);
            }
            replay_assays.push(item);
        }
        for (std::size_t i = 0; i < admitted_resources.size(); ++i) {
            auto item = Json::parse(resources.at(i).dump());
            item.set("source", admitted_resources.at(i).at("original"));
            if (admitted_resources.at(i).has("reference_original"))
                item.set("reference_source", admitted_resources.at(i).at("reference_original"));
            if (admitted_resources.at(i).has("provenance"))
                item.set("provenance", admitted_resources.at(i).at("provenance"));
            replay_resources.push(item);
        }
        replay.set("assays", replay_assays);
        replay.set("resources", replay_resources);
        const auto replay_bytes = replay.dump() + "\n";
        write_text(destination / "reimport.json", replay_bytes);
        manifest.set("reimport", file_record("reimport.json", replay_bytes));
        manifest.set("assays", admitted_assays);
        manifest.set("resources", admitted_resources);
        manifest.set("links", links);
        manifest.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
        manifest.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        manifest = seal(manifest, "manifest_sha256");
        write_text(destination / "workspace.json", manifest.dump() + "\n");
        auto publication = publish_native_workspace(destination,manifest,admitted_sources);
        if (observer)
            observer("workspace_publication");
        auto workspace = std::shared_ptr<const AssayWorkspace>(new AssayWorkspace(
            destination / "workspace.rfcxs", std::move(admitted_sources), std::move(publication)));
        if (observer)
            observer("workspace_validation");
        return workspace;
    } catch (...) {
        fs::remove_all(destination);
        throw;
    }
}
struct AssayWorkspace::State {
    fs::path root;
    RetainedRecord manifest;
    std::optional<NativeRecord> publication;
    std::map<std::string,RawByteSource> raw_artifacts;
    std::map<std::string, RetainedRecord> assays, physical;
    std::map<std::string, std::vector<std::shared_ptr<const ScientificStudy>>> sources;
    // Session reuse is scoped to the actual study and CPU coordinator owner.
    // No source tensors or result certificates are copied into this index.
    mutable std::mutex novelty_mutex;
    mutable std::shared_ptr<CpuResources> novelty_default_cpu;
    mutable std::map<std::tuple<std::string, std::size_t, CpuResources *>,
                     std::shared_ptr<ScientificSession>>
        novelty_sessions;
    explicit State(const fs::path &path, AssayWorkspace::AdmittedSources supplied = {},
                   std::optional<NativeRecord> admitted_publication = {}) {
        require(fs::is_regular_file(fs::symlink_status(path)),
                "workspace manifest must be regular");
        root = fs::canonical(path).parent_path();
        this->manifest = admitted_publication ? RetainedRecord(std::move(*admitted_publication))
                                             : read_record(path);
        if(this->manifest.has("schema") && (this->manifest.at("schema").text()=="gadget.assay.workspace.native.v1" || this->manifest.at("schema").text()=="gadget.assay.workspace.native.v2")) {
            publication=this->manifest.native();this->manifest=publication->at("manifest");
            const auto originals=publication->at("originals");
            for(std::size_t i=0;i<originals.size();++i)require(raw_artifacts.emplace(originals.at(i).at("path").text(),originals.at(i).at("bytes").byte_source()).second,"duplicate native input artifact");
        }
        using PublicationKey=std::pair<std::string,rfc::numeric::Integer>;
        std::map<PublicationKey,NativeRecord> publication_records;
        std::map<std::uint64_t,std::shared_ptr<const ScientificStudy>> publication_owners;
        if(publication) {
            const auto pubs=publication->at("publications");
            for(std::size_t i=0;i<pubs.size();++i) {
                const auto entry=pubs.at(i);
                require(publication_records.emplace(PublicationKey{entry.at("assay").text(),entry.at("partition").integer()},entry).second,"duplicate native biomedical publication binding");
            }
        }
        const auto manifest = this->manifest;
        const bool compact=manifest.at("schema").text()=="gadget.assays.workspace.v2";
        if(compact) {
            native_members(manifest.native(),{"schema","id","design","assays","resources","links","gadget_sha256","rexgraph_sha256"});
            require(publication && publication->at("schema").text()=="gadget.assay.workspace.native.v2","native collection publication required");
            require(manifest.at("resources").size()==0 && manifest.at("links").size()==0,"native study collection requires separately admitted resource/link bindings");
        } else {
        native_members(manifest.native(),
                {"schema", "id", "design", "import_spec", "reimport", "assays", "resources",
                 "links", "gadget_sha256", "rexgraph_sha256", "manifest_sha256"},
                {"bindings"});
        require(manifest.at("schema").text() == "gadget.assays.workspace.v1",
                "workspace schema");
        // The v1 checksum names the original compatibility representation.
        // Release that view before admitting any nested native publications.
        {
            const auto view = static_cast<Json>(manifest);
            require(manifest.at("manifest_sha256").text() ==
                        "sha256:" + sha256(unwrap(view).dump()),
                    "workspace manifest hash");
        }
        verify_artifact(manifest.at("import_spec").native());
        verify_artifact(manifest.at("reimport").native());
        if (manifest.has("bindings"))
            verify_artifact(manifest.at("bindings").native());
        }
        for (std::size_t i = 0; i < manifest.at("assays").size(); ++i) {
            const auto item = manifest.at("assays").at(i);
            native_members(item.native(), compact?std::set<std::string>{"declaration", "originals", "partitions"}:std::set<std::string>{"declaration", "originals", "features", "observations", "partitions"});
            const auto declaration = static_cast<Json>(item.at("declaration"));
            const auto id = text(declaration.at("id"), source_byte_limit());
            if(declaration.has("schema") && declaration.at("schema").text()=="gadget.assay.native.v1") {
                members(declaration,{"schema","id","kind","source_namespace","carrier"});
                require(declaration.at("carrier").text()=="source_grades","native assay requires its complete source grades");
                (void)key(declaration.at("source_namespace"));
            } else admit_matrix_spec(declaration);
            require(assays.emplace(id, this->manifest.member("assays").element(i)).second, "duplicate workspace assay");
            const auto originals = item.at("originals");
            for (const auto &name : originals.native().names())
                verify_artifact(originals.at(name).native());
            for (std::size_t p = 0; p < item.at("partitions").size(); ++p) {
                const auto packet = item.at("partitions").at(p);
                const auto packet_bytes = publication?std::string{}:local_artifact(root,packet);
                std::shared_ptr<const ScientificStudy> admitted;
                if (!supplied.empty()) {
                    require(supplied.contains(id) && supplied.at(id).size() ==
                                item.at("partitions").size(), "imported owner count");
                    admitted = supplied.at(id).at(p);
                    require(admitted && (publication || admitted->declaration().dump()==Json::parse(packet_bytes).dump()),"imported owner declaration");
                } else if(publication) {
                    const PublicationKey key{id,rfc::numeric::Integer(std::to_string(p))};
                    require(publication_records.contains(key),"native biomedical publication missing");
                    const auto entry=publication_records.at(key);
                    require(entry.at("path").text()==packet.at("path").text(),"native biomedical publication path binding");
                    const auto bytes=entry.at("bytes").byte_source();
                    if(!publication_owners.contains(bytes.key())) {
                        publication_owners.emplace(bytes.key(),std::make_shared<const ScientificStudy>(NativeRecord::decode(bytes.bytes())));
                    }
                    admitted=publication_owners.at(bytes.key());
                } else {
                    admitted = std::make_shared<const ScientificStudy>(root / packet.at("path").text());
                }
                require(admitted->digest() == packet.at("source_digest").text(),
                        "workspace native packet source digest");
                sources[id].push_back(std::move(admitted));
            }
        }
        for (std::size_t i = 0; i < manifest.at("resources").size(); ++i) {
            const auto item = manifest.at("resources").at(i);
            const auto id = text(item.at("declaration").at("id"), source_byte_limit());
            require(physical.emplace(id, this->manifest.member("resources").element(i)).second, "duplicate workspace resource");
            verify_artifact(item.at("original").native());
            if (item.has("reference_original"))
                verify_artifact(item.at("reference_original").native());
            if (item.has("provenance"))
                for (std::size_t p = 0; p < item.at("provenance").size(); ++p)
                    verify_artifact(item.at("provenance").at(p).native());
        }
        if(!publication) {
            publication=admit_native_workspace(root,static_cast<Json>(manifest),sources);this->manifest=publication->at("manifest");
            for(std::size_t i=0;i<manifest.at("assays").size();++i)assays.at(manifest.at("assays").at(i).at("declaration").at("id").text())=this->manifest.member("assays").element(i);
            for(std::size_t i=0;i<manifest.at("resources").size();++i)physical.at(manifest.at("resources").at(i).at("declaration").at("id").text())=this->manifest.member("resources").element(i);
        }
        if(raw_artifacts.empty()) {
            const auto originals=publication->at("originals");
            for(std::size_t i=0;i<originals.size();++i)raw_artifacts.emplace(originals.at(i).at("path").text(),originals.at(i).at("bytes").byte_source());
        }
    }
    template <class Record>
    static Json artifact_reference(const Record &record) {
        auto reference = Json::object();
        reference.set("path", Json::string(record.at("path").text()));
        reference.set("sha256", Json::string(record.at("sha256").text()));
        return reference;
    }
    template <class Record>
    const RawByteSource *retained_artifact(const Record &record) const {
        const auto path = record.at("path").text();
        const fs::path relative(path);
        require(!relative.empty() && !relative.is_absolute() &&
                    path.find('\0') == std::string::npos, "invalid native artifact path");
        for (const auto &part : relative)
            require(part != "..", "native artifact escapes workspace");
        const auto found = raw_artifacts.find(path);
        return found == raw_artifacts.end() ? nullptr : &found->second;
    }
    template <class Record>
    void verify_artifact(const Record &record) const {
        const auto *source = retained_artifact(record);
        if (!source) {
            (void)local_artifact(root, artifact_reference(record));
            return;
        }
        const auto data = source->bytes();
        require("sha256:" + sha256(std::string_view(
                    reinterpret_cast<const char *>(data.data()), data.size())) ==
                    record.at("sha256").text(), "native original artifact binding mismatch");
        if (fs::exists(root / record.at("path").text())) {
            const auto path = local_artifact_path(root, artifact_reference(record));
            const ReadOnlyMappedFile file(path);
            require(file.bytes().size() == data.size() &&
                        std::equal(data.begin(), data.end(), file.bytes().begin()),
                    "external artifact differs from native source");
        }
    }
    template <class Record>
    std::string artifact(const Record &record) const {
        const auto *source = retained_artifact(record);
        if (!source)
            return local_artifact(root, artifact_reference(record));
        verify_artifact(record);
        const auto data = source->bytes();
        return std::string(reinterpret_cast<const char *>(data.data()), data.size());
    }
};
AssayWorkspace::AssayWorkspace(const fs::path &path) : state_(std::make_unique<State>(path)) {}
AssayWorkspace::AssayWorkspace(const fs::path &path, AdmittedSources sources)
    : state_(std::make_unique<State>(path, std::move(sources))) {}
AssayWorkspace::AssayWorkspace(const fs::path &path, AdmittedSources sources,
                               NativeRecord publication)
    : state_(std::make_unique<State>(path, std::move(sources), std::move(publication))) {}
AssayWorkspace::~AssayWorkspace() = default;
const NativeRecord &AssayWorkspace::native_publication() const {return state_->publication.value();}
void AssayWorkspace::write_native(const fs::path &path) const {native_publication().write(path);}
Json AssayWorkspace::native_novelty_preparation_statistics() const {
    std::lock_guard lock(state_->novelty_mutex);
    auto out = Json::array();
    for (const auto &[key, session] : state_->novelty_sessions) {
        auto item = Json::object();
        item.set("assay", Json::string(std::get<0>(key)));
        item.set("partition", cardinality(std::get<1>(key)));
        item.set("preparation", session->native_novelty_preparation_statistics());
        out.push(item);
    }
    return out;
}
void AssayWorkspace::clear_native_novelty_preparations() {
    std::lock_guard lock(state_->novelty_mutex);
    state_->novelty_sessions.clear();
    state_->novelty_default_cpu.reset();
}
std::string AssayWorkspace::digest() const { return state_->manifest.has("manifest_sha256")?state_->manifest.at("manifest_sha256").text():native_declaration_digest(state_->manifest.native()); }
Json AssayWorkspace::catalog() const {
    auto result = Json::parse(state_->manifest.dump());
    result.set("capabilities", bioinformatics_capabilities());
    if(!result.has("manifest_sha256"))result.set("manifest_sha256",Json::string(digest()));
    return result;
}
Json AssayWorkspace::source_context(const Json &input) const {
    const auto request = Json::parse(input.dump());
    members(request, {"schema", "source_digest", "question", "targets"});
    require(request.at("schema").text() == "gadget.assays.source.context.v1" &&
                request.at("source_digest").text() == digest(),
            "source context schema/workspace mismatch");
    const auto question = text(request.at("question"), source_byte_limit());
    const auto entries = request.at("targets");
    source_array(entries, "source context targets");
    require(entries.size() > 0, "source context requires declared targets");
    std::set<std::string> ids;
    std::map<std::string, std::size_t> last_resource;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto target = entries.at(i).at("target");
        if (target.at("kind").text() == "resource")
            last_resource[target.at("resource").text()] = i;
    }
    std::map<std::string, std::map<std::string, Document>> resource_documents;
    std::map<std::string,std::map<std::string,std::string>> resource_bytes;
    std::map<std::string, Json> resource_records;
    auto targets = Json::array();
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto entry = entries.at(i), target = entry.at("target");
        members(entry, {"id", "target", "evidence"});
        require(ids.insert(text(entry.at("id"), source_byte_limit())).second,
                "duplicate source context target identity");
        const auto kind = target.at("kind").text();
        Json resolved;
        if (kind == "resource") {
            members(target, {"kind", "resource", "artifact", "sha256"}, {"record"});
            const auto id = text(target.at("resource"), source_byte_limit());
            require(state_->physical.contains(id), "source context resource is undeclared");
            const auto &resource = state_->physical.at(id);
            const auto artifact = target.at("artifact").text();
            require(artifact == "source" || artifact == "reference", "source resource artifact");
            const auto original =
                resource.at(artifact == "source" ? "original" : "reference_original");
            require(target.at("sha256").text() == original.at("sha256").text(),
                    "source resource owner/hash mismatch");
            auto &documents = resource_documents[id];
            auto load = [&](const std::string &name) -> const Document & {
                if (!documents.contains(name)) {
                    const auto record =
                        resource.at(name == "source" ? "original" : "reference_original");
                    auto bytes = state_->artifact( record);
                    require(bytes.find('\0') == std::string::npos,
                            "source evidence requires UTF-8 bytes");
                    (void)Json::parse(Json::string(bytes).dump());
                    auto metadata = Json::parse(record.dump());
                    metadata.set("format", Json::string("utf8"));
                    auto &retained=resource_bytes[id].emplace(name,std::move(bytes)).first->second;
                    documents.emplace(name, Document{metadata,retained});
                }
                return documents.at(name);
            };
            const auto &bytes = load(artifact).bytes;
            source_array(entry.at("evidence"), "resource evidence");
            for (std::size_t e = 0; e < entry.at("evidence").size(); ++e)
                require(entry.at("evidence").at(e).at("document").text() == artifact,
                        "resource evidence belongs to another artifact");
            const auto evidence = source_evidence(entry.at("evidence"), documents, true);
            auto descriptor = Json::parse(resource.at("declaration").dump());
            if (target.has("record")) {
                require(artifact == "source", "reference annotation requires a file target");
                // Reuse the format decoder; do not trust detached inspection metadata.
                if (!resource_records.contains(id)) {
                    const auto reference = resource.has("reference_original")
                                               ? load("reference").bytes
                                               : std::string_view{};
                    resource_records.emplace(
                        id, inspect_resource(std::string(bytes), resource.at("declaration"), std::string(reference))
                                .at("records"));
                }
                const auto records = resource_records.at(id);
                require(records.size() > 0, "source resource has no decoded records");
                const auto record = records.at(
                    count(target.at("record"), records.size() - 1, "source resource record"));
                require(record.has("source_start") && record.has("source_end"),
                        "decoded record has no original text span");
                for (std::size_t e = 0; e < evidence.size(); ++e)
                    require(evidence.at(e).at("start").integer() >=
                                    record.at("source_start").integer() &&
                                evidence.at(e).at("end").integer() <=
                                    record.at("source_end").integer(),
                            "annotation evidence is outside the selected resource record");
                descriptor.set("record", record);
            }
            resolved = Json::object();
            resolved.set("descriptor", descriptor);
            resolved.set("source_sha256", original.at("sha256"));
            resolved.set("evidence", evidence);
            if (last_resource.at(id) == i) {
                resource_documents.erase(id);
                resource_bytes.erase(id);
                resource_records.erase(id);
            }
        } else {
            members(target, {"kind", "assay", "partition", "source_digest"},
                    kind == "occurrence" ? std::set<std::string>{"grade", "key"}
                                         : std::set<std::string>{"id"});
            const auto assay = text(target.at("assay"), source_byte_limit());
            require(state_->sources.contains(assay) && !state_->sources.at(assay).empty(),
                    "source context assay is undeclared or empty");
            const auto partition =
                count(target.at("partition"), state_->sources.at(assay).size() - 1,
                      "source context partition");
            auto selector = Json::object();
            json_object_object_foreach(target.raw(), name, value) {
                (void)value;
                if (std::string(name) != "assay" && std::string(name) != "partition")
                    selector.set(name, target.at(name));
            }
            resolved = state_->sources.at(assay).at(partition)->source_target(selector,
                                                                              entry.at("evidence"));
        }
        auto bound = Json::object();
        bound.set("id", entry.at("id"));
        bound.set("target", target);
        bound.set("resolved", resolved);
        bound.set("evidence_verified", Json::boolean(true));
        targets.push(seal(bound, "binding_sha256"));
    }
    auto result = Json::object();
    result.set("schema", Json::string("gadget.assays.source.context.catalog.v1"));
    result.set("source_digest", Json::string(digest()));
    result.set("question", Json::string(question));
    result.set("question_sha256", Json::string("sha256:" + sha256(question)));
    result.set("request", request);
    result.set("targets", targets);
    return seal(result, "context_sha256");
}
Json AssayWorkspace::read_source(const Json &input) const {
    const auto request = Json::parse(input.dump());
    members(request, {"schema", "source_digest", "assay", "partition", "read"});
    require(request.at("schema").text() == "gadget.assays.source.read.v1" &&
                request.at("source_digest").text() == digest(),
            "assay source read schema/owner mismatch");
    const auto assay = text(request.at("assay"));
    require(state_->sources.contains(assay) && !state_->sources.at(assay).empty(),
            "source read assay undeclared or empty");
    const auto partition = count(request.at("partition"), state_->sources.at(assay).size() - 1,
                                 "source read partition");
    auto result = Json::object();
    result.set("schema", Json::string("gadget.assays.source.page.v1"));
    result.set("source_digest", Json::string(digest()));
    result.set("assay", request.at("assay"));
    result.set("partition", request.at("partition"));
    if (request.at("read").at("schema").text() == "gadget.scientific.source.index.v1") {
        members(request.at("read"), {"schema"});
        result.set("index", state_->sources.at(assay).at(partition)->source_index());
        result.set("next_request", Json{});
        return seal(result, "page_sha256");
    }
    const auto page = state_->sources.at(assay).at(partition)->read_source(request.at("read"));
    result.set("page", page);
    if (!page.at("complete").boolean()) {
        auto next = Json::parse(request.dump());
        next.set("read", page.at("next_request"));
        result.set("next_request", next);
    } else
        result.set("next_request", Json{});
    return seal(result, "page_sha256");
}
Json AssayWorkspace::bind_inputs(
    const std::map<std::string, std::shared_ptr<const AssayWorkspace>> &inputs, const Json &request,
    const fs::path &destination) {
    members(request,
            {"schema", "id", "design", "assays", "resources", "links", "method", "evidence"});
    require(request.at("schema").text() == "gadget.assays.bindings.v1", "assay bindings schema");
    (void)text(request.at("id"), source_byte_limit());
    (void)text(request.at("method"), source_byte_limit());
    (void)text(request.at("evidence"), source_byte_limit());
    members(request.at("design"), {"question", "population", "selection", "limitations"});
    for (const auto *name : {"question", "population", "selection"})
        (void)text(request.at("design").at(name), source_byte_limit());
    source_array(request.at("design").at("limitations"), "binding limitations");
    for (std::size_t i = 0; i < request.at("design").at("limitations").size(); ++i)
        (void)text(request.at("design").at("limitations").at(i), source_byte_limit());
    for (const auto *key : {"assays", "resources", "links"})
        source_array(request.at(key), key);
    require(request.at("assays").size() || request.at("resources").size(),
            "bindings require selected assays or resources");
    require(!fs::exists(destination), "bound workspace destination exists");
    auto receipt = Json::object(), sources = Json::array();
    std::set<std::string> used_inputs;
    const auto source_for = [&](const Json &selection) {
        const auto name = text(selection.at("input"), source_byte_limit());
        require(inputs.contains(name) && static_cast<bool>(inputs.at(name)),
                "unknown binding input: " + name);
        const auto &source = inputs.at(name);
        if (used_inputs.insert(name).second) {
            auto record = Json::object();
            record.set("input", selection.at("input"));
            record.set("workspace_sha256", Json::string(source->digest()));
            sources.push(record);
        }
        return source;
    };
    bool native_selection=false;
    for(std::size_t a=0;a<request.at("assays").size();++a) {
        const auto selection=request.at("assays").at(a);
        members(selection,{"input","assay","id"});
        const auto source=source_for(selection);
        require(source->state_->assays.contains(selection.at("assay").text()),"unknown bound assay");
        const auto declaration=source->state_->assays.at(selection.at("assay").text()).at("declaration");
        native_selection=native_selection || declaration.has("schema");
    }
    if(native_selection) {
        require(request.at("resources").size()==0 && request.at("links").size()==0,"native study bindings require separately admitted resource and link declarations");
        std::map<std::string,std::shared_ptr<const ScientificStudy>> selected;
        for(std::size_t a=0;a<request.at("assays").size();++a) {
            const auto selector=request.at("assays").at(a);const auto source=source_for(selector);
            require(source->state_->sources.at(selector.at("assay").text()).size()==1,"native collection binding requires an indivisible admitted study");
            require(selected.emplace(text(selector.at("id")),source->study(selector.at("assay").text())).second,"duplicate native assay alias");
        }
        const auto publication=publish_studies(admit_record_view(request),request.at("id").text(),selected);
        fs::create_directories(fs::absolute(destination).parent_path());
        const auto stage=private_directory(fs::absolute(destination).parent_path(),".gadget-native-bind-");
        struct Cleanup {fs::path path;~Cleanup(){std::error_code ec;fs::remove_all(path,ec);}} cleanup{stage};
        publication.write(stage/"workspace.rfcxs");
        auto result=render_record_view(publication.at("manifest"));
        result.set("manifest_sha256",Json::string(native_declaration_digest(publication.at("manifest"))));
        result.set("capabilities",bioinformatics_capabilities());fs::rename(stage,destination);return result;
    }
    // Admit all selectors before copying anything. Collection aliases are
    // request routing names; native primary identities and namespaces remain.
    std::set<std::string> assay_ids, resource_ids;
    std::set<std::uint64_t> namespaces;
    std::map<std::string, std::set<std::string>> observations;
    auto native_bindings = Json::array();
    for (std::size_t a = 0; a < request.at("assays").size(); ++a) {
        const auto selection = request.at("assays").at(a);
        members(selection, {"input", "assay", "id"});
        const auto source = source_for(selection);
        const auto alias = text(selection.at("id"), source_byte_limit());
        require(assay_ids.insert(alias).second &&
                    source->state_->assays.contains(selection.at("assay").text()),
                "duplicate alias or unknown selected assay");
        const auto item = source->state_->assays.at(selection.at("assay").text());
        for (std::size_t p = 0; p < item.at("partitions").size(); ++p)
            require(
                namespaces.insert(key(item.at("partitions").at(p).at("source_namespace"))).second,
                "bound collection requires distinct declared native namespaces for portable "
                "reimport");
        for (std::size_t i = 0; i < item.at("observations").size(); ++i)
            observations[alias].insert(item.at("observations").at(i).at("id").text());
        auto record = Json::parse(selection.dump());
        record.set("source_namespace", item.at("declaration").at("source_namespace"));
        record.set("partitions", item.at("partitions"));
        native_bindings.push(record);
    }
    for (std::size_t r = 0; r < request.at("resources").size(); ++r) {
        const auto selection = request.at("resources").at(r);
        members(selection, {"input", "resource", "id"});
        const auto source = source_for(selection);
        require(resource_ids.insert(text(selection.at("id"), source_byte_limit())).second &&
                    source->state_->physical.contains(selection.at("resource").text()),
                "duplicate alias or unknown selected resource");
    }
    std::set<std::string> link_ids;
    for (std::size_t i = 0; i < request.at("links").size(); ++i) {
        const auto link = request.at("links").at(i);
        members(link, {"id", "left_assay", "left_observation", "right_assay", "right_observation",
                       "relation", "evidence"});
        require(link_ids.insert(text(link.at("id"), source_byte_limit())).second,
                "duplicate binding link");
        for (const auto *side : {"left", "right"}) {
            const auto assay = link.at(std::string(side) + "_assay").text();
            require(observations.contains(assay) &&
                        observations.at(assay).contains(
                            link.at(std::string(side) + "_observation").text()),
                    "binding link endpoint absent from selected assay");
        }
        (void)text(link.at("relation"), source_byte_limit());
        (void)text(link.at("evidence"), source_byte_limit());
    }
    receipt.set("schema", Json::string("gadget.assays.bindings.receipt.v1"));
    receipt.set("request", request);
    receipt.set("sources", sources);
    receipt.set("native_assays", native_bindings);
    receipt.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
    receipt.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    receipt = seal(receipt, "certificate_sha256");
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage = private_directory(fs::absolute(destination).parent_path(), ".gadget-bind-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    const auto binding_bytes = receipt.dump() + "\n";
    write_text(stage / "bindings.json", binding_bytes);
    const auto binding_record = file_record("bindings.json", binding_bytes);
    auto admitted = Json::array(), physical = Json::array(), replay_assays = Json::array(),
         replay_resources = Json::array();
    AdmittedSources bound_sources;
    const auto copy_bytes = [&](const std::string &bytes, const Json &original, const std::string &path) {
        fs::create_directories((stage / path).parent_path());
        write_text(stage / path, bytes);
        auto copied = Json::parse(original.dump());
        copied.set("path", Json::string(path));
        return copied;
    };
    const auto copy_native = [&](const AssayWorkspace &source, const Json &original, const std::string &path) {
        return copy_bytes(source.state_->artifact(original), original, path);
    };
    for (std::size_t a = 0; a < request.at("assays").size(); ++a) {
        const auto selection = request.at("assays").at(a);
        const auto source = inputs.at(selection.at("input").text());
        auto item = Json::parse(source->state_->assays.at(selection.at("assay").text()).dump());
        const auto dir = "assays/" + std::to_string(a);
        auto originals = Json::object();
        json_object_object_foreach(item.at("originals").raw(), name, value) {
            (void)value;
            originals.set(name, copy_native(*source, item.at("originals").at(name),
                                     dir + "/source-" + name));
        }
        auto declaration = Json::parse(item.at("declaration").dump());
        declaration.set("id", selection.at("id"));
        declaration.set("source", originals.at("matrix"));
        for (const auto *name : {"features", "observations", "missing", "relations"})
            if (originals.has(name))
                declaration.set(name, originals.at(name));
        auto provenance = Json::array();
        if (declaration.has("provenance"))
            for (std::size_t i = 0; i < declaration.at("provenance").size(); ++i)
                provenance.push(originals.at("provenance-" + std::to_string(i)));
        originals.set("provenance-" + std::to_string(provenance.size()), binding_record);
        provenance.push(binding_record);
        declaration.set("provenance", provenance);
        auto packets = Json::array();
        for (std::size_t p = 0; p < item.at("partitions").size(); ++p) {
            const auto packet = item.at("partitions").at(p);
            const auto owner = source->study(selection.at("assay").text(), p);
            const auto partition_dir = dir + "/partitions/" + std::to_string(p) + "/";
            auto retained = Json::parse(packet.dump());
            fs::create_directories(stage / partition_dir);
            const auto path = partition_dir + "study.rfcxs";
            owner->native_publication().write(stage / path);
            const auto artifact = file_record(path, read_text(stage / path));
            retained.set("path", artifact.at("path"));
            retained.set("sha256", artifact.at("sha256"));
            bound_sources[selection.at("id").text()].push_back(owner);
            packets.push(retained);
        }
        item.set("declaration", declaration);
        item.set("originals", originals);
        item.set("partitions", packets);
        admitted.push(item);
        replay_assays.push(declaration);
    }
    for (std::size_t r = 0; r < request.at("resources").size(); ++r) {
        const auto selection = request.at("resources").at(r);
        const auto source = inputs.at(selection.at("input").text());
        auto item =
            Json::parse(source->state_->physical.at(selection.at("resource").text()).dump());
        auto declaration = Json::parse(item.at("declaration").dump());
        declaration.set("id", selection.at("id"));
        item.set("original", copy_native(*source, item.at("original"),
                                  "resources/" + std::to_string(r) + "/source"));
        declaration.set("source", item.at("original"));
        if (item.has("reference_original")) {
            item.set("reference_original", copy_native(*source, item.at("reference_original"),
                                                "resources/" + std::to_string(r) + "/reference"));
            declaration.set("reference_source", item.at("reference_original"));
        }
        if (item.has("provenance")) {
            auto originals = Json::array();
            for (std::size_t p = 0; p < item.at("provenance").size(); ++p)
                originals.push(copy_native(*source, item.at("provenance").at(p),
                                    "resources/" + std::to_string(r) +
                                    "/provenance-" + std::to_string(p)));
            item.set("provenance", originals);
            declaration.set("provenance", originals);
        }
        item.set("declaration", declaration);
        physical.push(item);
        replay_resources.push(declaration);
    }
    auto replay = Json::object();
    replay.set("schema", Json::string("gadget.assays.import.v1"));
    replay.set("id", request.at("id"));
    replay.set("design", request.at("design"));
    replay.set("assays", replay_assays);
    replay.set("resources", replay_resources);
    replay.set("links", request.at("links"));
    const auto replay_bytes = replay.dump() + "\n";
    write_text(stage / "reimport.json", replay_bytes);
    auto manifest = Json::object();
    manifest.set("schema", Json::string("gadget.assays.workspace.v1"));
    manifest.set("id", request.at("id"));
    manifest.set("design", request.at("design"));
    manifest.set("import_spec", binding_record);
    manifest.set("reimport", file_record("reimport.json", replay_bytes));
    manifest.set("bindings", binding_record);
    manifest.set("assays", admitted);
    manifest.set("resources", physical);
    manifest.set("links", request.at("links"));
    manifest.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
    manifest.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    manifest = seal(manifest, "manifest_sha256");
    write_text(stage / "workspace.json", manifest.dump() + "\n");
    publish_native_workspace(stage,manifest,bound_sources);
    const auto catalog = AssayWorkspace(stage / "workspace.rfcxs",bound_sources).catalog();
    fs::rename(stage, destination);
    return catalog;
}
std::shared_ptr<const ScientificStudy> AssayWorkspace::study(const std::string &assay,
                                                             std::size_t partition) const {
    require(state_->assays.contains(assay), "unknown workspace assay");
    return state_->sources.at(assay).at(partition);
}
Json AssayWorkspace::execute(const Json &request, std::shared_ptr<CpuResources> resources) const {
    members(request, {"schema", "operation"},
            {"assay",
             "query",
             "partitions",
             "resource",
             "thresholds",
             "strand_policy",
             "left",
             "right",
             "map",
             "metric",
             "motifs",
             "strands",
             "targets",
             "candidate_search",
             "kernel_certificate",
             "scale",
             "scale_provenance",
             "annotation",
             "feature_type",
             "group_attribute",
             "unit",
             "overlap",
             "filter",
             "samples",
             "contrasts",
             "absence_policy",
             "filter_policy",
             "quality_thresholds",
             "coordinates",
             "frame",
             "units",
             "metric_diagonal",
             "radius_squared",
             "relation_type",
             "evidence",
             "feature_groups",
             "observation_groups",
             "feature_design",
             "observation_design",
             "map_method",
             "map_evidence",
             "support_policy",
             "count_policy",
             "query_resource",
             "query_id",
             "reference_resource",
             "reference_id",
             "scoring",
             "ambiguity",
             "mode",
             "domain",
             "variables",
             "constraints",
             "objective",
             "witness_representation",
             "factor_certificate",
             "partition",
             "levels",
             "births",
             "filtration_meaning",
             "maximum_mismatches",
             "evidence_source",
             "coupling",
             "junction_policy",
             "junctions",
             "minimum_overhang",
             "k",
             "route",
             "lattice_snapshots",
             "maximum_edits",
             "minimum_exon_bases",
             "minimum_intron_bases",
             "maximum_intron_bases",
             "maximum_junctions",
             "motifs",
             "minimum_overlap",
             "queries",
             "concurrent_jobs",
             "workers",
             "controls",
             "intercept",
             "source_kernel",
             "execution_certificate",
             "views",
             "sections",
             "interfaces",
             "declaration",
             "sources",
             "embeddings",
             "boundary_policy"});
    require(request.at("schema").text() == "gadget.assays.query.v1", "assay query schema");
    const auto operation = request.at("operation").text();
    if (!resources) {
        if (operation == "native_query" &&
            request.at("query").at("operation").text() == "metric_relative_novelty") {
            std::lock_guard lock(state_->novelty_mutex);
            if (!state_->novelty_default_cpu)
                state_->novelty_default_cpu = CpuResources::local();
            resources = state_->novelty_default_cpu;
        } else
            resources = CpuResources::local();
    }
    auto result = Json::object();
    if (operation == "query_batch") {
        result = query_batch(*this, request, resources);
    } else if (operation == "inspect") {
        members(request, {"schema", "operation"});
        result = catalog();
    } else if (operation == "native_query") {
        members(request, {"schema", "operation", "assay", "query"}, {"partitions"});
        const auto id = text(request.at("assay"));
        require(state_->assays.contains(id), "unknown assay");
        const auto item = state_->assays.at(id);
        const auto packets = item.at("partitions");
        std::vector<std::size_t> selected;
        if (request.has("partitions")) {
            std::set<std::size_t> unique;
            for (std::size_t i = 0; i < request.at("partitions").size(); ++i) {
                const auto p =
                    count(request.at("partitions").at(i), packets.size() - 1, "assay partition");
                require(unique.insert(p).second, "duplicate assay partition");
                selected.push_back(p);
            }
            require(!selected.empty(), "empty assay partition selection");
        } else
            for (std::size_t i = 0; i < packets.size(); ++i)
                selected.push_back(i);
        // Native subspace/topology/Green queries require an indivisible carrier.
        const auto op = request.at("query").at("operation").text();
        require(packets.size() == 1 || op == "paired_field_contrasts" || op == "field_summary" ||
                    op == "coverage" || op == "tensor_window" || op == "inspect" ||
                    op == "provenance",
                "whole-carrier method cannot be applied independently to storage partitions");
        auto answers = Json::array();
        const bool reuse_novelty = op == "metric_relative_novelty" &&
                                   (!request.at("query").has("preparation") ||
                                    request.at("query").at("preparation").text() == "reuse");
        const auto session_for = [&](std::size_t p) {
            if (!reuse_novelty)
                return std::make_shared<ScientificSession>(study(id, p), resources);
            std::lock_guard lock(state_->novelty_mutex);
            const auto key = std::tuple{id, p, resources.get()};
            const auto found = state_->novelty_sessions.find(key);
            if (found != state_->novelty_sessions.end())
                return found->second;
            auto session = std::make_shared<ScientificSession>(study(id, p), resources);
            state_->novelty_sessions.emplace(key, session);
            return session;
        };
        // Preflight every leaf before publishing any answer.
        for (auto p : selected) {
            session_for(p)->validate(request.at("query"));
        }
        for (auto p : selected) {
            auto answer = Json::object();
            answer.set("partition", cardinality(p));
            if(packets.at(p).has("carrier_ids"))answer.set("carrier_ids",packets.at(p).at("carrier_ids"));
            else answer.set("source_index",study(id,p)->source_index());
            answer.set("certificate", session_for(p)->execute(request.at("query")));
            answers.push(answer);
        }
        result.set("assay", request.at("assay"));
        result.set("answers", answers);
        result.set("carrier_axis", item.at("declaration").at("carrier"));
        result.set(
            "partition_semantics",
            Json::string("Ordered disjoint carrier subsets with explicit identity witnesses; "
                         "whole-carrier topological operations require one partition"));
    } else if (operation == "constrained_allocations") {
        auto reservation = resources->reserve(1);
        result = constrained_allocations(request);
    } else if (operation == "relational_persistence") {
        const auto id = text(request.at("assay"));
        require(state_->assays.contains(id), "unknown persistence assay");
        require(state_->assays.at(id).at("partitions").size() == 1 &&
                    count(request.at("partition"), 0, "persistence partition") == 0,
                "persistence requires an indivisible source carrier");
        auto reservation = resources->reserve(1);
        result = relational_persistence(*study(id), request);
    } else if (operation == "aggregate_assay" || operation == "design_aggregation") {
        require(state_->assays.contains(text(request.at("assay"))), "unknown aggregate assay");
        const auto source = state_->assays.at(request.at("assay").text());
        auto declaration = Json::parse(source.at("declaration").dump());
        const auto originals = source.at("originals");
        declaration.set("source", originals.at("matrix"));
        for (const auto *name : {"features", "observations", "missing"})
            if (originals.has(name))
                declaration.set(name, originals.at(name));
        const auto workers = request.has("workers")
                                 ? count(request.at("workers"), resources->coordinator().capacity(),
                                         "aggregation workers")
                                 : 1;
        require(workers > 0, "aggregation workers must be positive");
        auto reservation = resources->reserve(workers);
        Matrix axes;
        axes.features = read_axis(state_->artifact( originals.at("features")));
        axes.observations = read_axis(state_->artifact( originals.at("observations")));
        const auto &axis = declaration.at("carrier").text() == "features" ? axes.features
                                                                         : axes.observations;
        std::vector<NativeAssayPartition> partitions;
        const auto packets = source.at("partitions");
        require(packets.size() > 0, "native aggregation requires stored packets");
        for (std::size_t p = 0; p < packets.size(); ++p) {
            const auto packet = packets.at(p);
            const auto begin = count(packet.at("begin"), axis.ids.size(), "native packet begin");
            const auto end = count(packet.at("end"), axis.ids.size(), "native packet end");
            require(end > begin, "native packet range");
            const auto ids = packet.at("carrier_ids");
            require(ids.size() == end - begin, "native packet identity count");
            for (std::size_t i = 0; i < ids.size(); ++i)
                require(ids.at(i).text() == axis.ids[begin + i], "native packet axis identity/order");
            const auto owner = study(request.at("assay").text(), p);
            require(owner->native_history()->source_namespace() == key(packet.at("source_namespace")) &&
                        owner->digest() == packet.at("source_digest").text(),
                    "native packet publication witness");
            partitions.push_back({begin, end, owner});
        }
        result = operation == "design_aggregation"
                     ? design_aggregation(axes, declaration, request, reservation, partitions)
                     : native_aggregate_analysis(axes.features, axes.observations, partitions,
                                                 declaration, request, reservation);
    } else if (operation == "spatial_relations") {
        require(state_->assays.contains(text(request.at("assay"))), "unknown spatial assay");
        auto reservation = resources->reserve(1);
        result = spatial_analysis(state_->assays.at(request.at("assay").text()), request);
    } else if (operation == "feature_counts" || operation == "expression_allocations" ||
               operation == "alignment_evidence" || operation == "variant_genotypes" ||
               operation == "annotation_features" || operation == "chromatin_contrast" ||
               operation == "controlled_coverage" || operation == "read_quality" ||
               operation == "exact_sequence_alignment" ||
               operation == "reference_sequence_matches" ||
               operation == "sequence_assembly_graph" ||
               operation == "indexed_sequence_alignment" ||
               operation == "sequence_splice_discovery" ||
               operation == "read_overlap_reconstruction") {
        ResourceImages images;
        const auto load = [&](const Json &id) {
            const auto name = text(id);
            require(state_->physical.contains(name), "unknown physical resource");
            if (images.contains(name))
                return;
            const auto item = state_->physical.at(name);
            images.emplace(name, ResourceImage{item.at("declaration"),
                                               state_->artifact( item.at("original")),
                                               item.has("reference_original")
                                                   ? state_->artifact(
                                                                    item.at("reference_original"))
                                                   : std::string{}});
        };
        if (request.has("resource"))
            load(request.at("resource"));
        if (request.has("annotation"))
            load(request.at("annotation"));
        if (request.has("query_resource"))
            load(request.at("query_resource"));
        if (request.has("reference_resource"))
            load(request.at("reference_resource"));
        if (request.has("samples"))
            for (std::size_t i = 0; i < request.at("samples").size(); ++i)
                load(request.at("samples").at(i).at("resource"));
        if (request.has("evidence_source"))
            for (const auto *name : {"annotation", "reference_resource"})
                if (request.at("evidence_source").has(name))
                    load(request.at("evidence_source").at(name));
        result = genomic_analysis(images, request, resources);
    } else if (operation == "assay_sections") {
        result = assay_sections(*this, request, resources);
    } else if (operation == "relational_gluing") {
        result = relational_gluing(*this, request, resources);
    } else if (operation == "field_transfer") {
        members(request, {"schema", "operation", "left", "right", "map", "metric"});
        require(request.at("metric").text() == "identity", "field transfer metric profile");
        const auto load = [&](const Json &record) {
            members(record, {"assay", "partition", "field"});
            const auto a = text(record.at("assay"));
            require(state_->assays.contains(a), "field transfer assay");
            const auto p =
                count(record.at("partition"), state_->assays.at(a).at("partitions").size() - 1,
                      "field transfer partition");
            auto input = load_assay_field(study(a, p), record.at("field").text());
            require(!input.values.empty(), "empty field transfer input");
            require(std::all_of(input.supported.begin(), input.supported.end(),
                                [](auto v) { return v; }),
                    "field transfer requires complete support; no implicit missing extension");
            return input;
        };
        const auto left = load(request.at("left")), right = load(request.at("right"));
        const auto map = request.at("map");
        members(map, {"domain_quantity", "domain_units", "codomain_quantity", "codomain_units",
                      "method", "evidence", "entries"});
        require(map.at("domain_quantity").text() == right.meaning.at("quantity").text() &&
                    map.at("domain_units").text() == right.meaning.at("units").text() &&
                    map.at("codomain_quantity").text() == left.meaning.at("quantity").text() &&
                    map.at("codomain_units").text() == left.meaning.at("units").text(),
                "field transfer domain/codomain meanings differ from measurements");
        (void)text(map.at("method"));
        (void)text(map.at("evidence"), 16384);
        std::map<std::string, std::size_t> li, ri;
        for (std::size_t i = 0; i < left.ids.size(); ++i)
            li.emplace(left.ids.at(i).text(), i);
        for (std::size_t i = 0; i < right.ids.size(); ++i)
            ri.emplace(right.ids.at(i).text(), i);
        rfc::numeric::ExactMatrix matrix(left.values.size(), right.values.size());
        std::set<std::pair<std::size_t, std::size_t>> seen;
        for (std::size_t i = 0; i < map.at("entries").size(); ++i) {
            const auto e = map.at("entries").at(i);
            members(e, {"left_id", "right_id", "coefficient"});
            const auto l = text(e.at("left_id")), r = text(e.at("right_id"));
            require(li.contains(l) && ri.contains(r) && seen.emplace(li.at(l), ri.at(r)).second,
                    "field transfer coordinate identity or duplicate entry");
            matrix.set(li.at(l), ri.at(r), exact_rational(text(e.at("coefficient"))));
        }
        auto reservation = resources->reserve(1);
        const auto transported = matrix.apply(right.values),
                   residual = rfc::numeric::sub(left.values, transported);
        const auto vector_json = [](const Vector &v) {
            auto out = Json::array();
            for (const auto &q : v)
                out.push(Json::string(q.get_str()));
            return out;
        };
        result.set("left", vector_json(left.values));
        result.set("right", vector_json(right.values));
        result.set("transported_right", vector_json(transported));
        result.set("residual", vector_json(residual));
        result.set("left_ids", left.ids);
        result.set("right_ids", right.ids);
        result.set("left_carrier_grade", Json::integer(left.grade));
        result.set("right_carrier_grade", Json::integer(right.grade));
        const auto rank = matrix.rank();
        result.set("map_rank", cardinality(rank));
        result.set("map_kernel_dimension", cardinality(matrix.cols() - rank));
        result.set("residual_quadrance",
                   Json::string(rfc::numeric::dot(residual, residual).get_str()));
        result.set("pairing", Json::string(rfc::numeric::dot(left.values, transported).get_str()));
        result.set("left_source_digest", Json::string(left.owner->digest()));
        result.set("right_source_digest", Json::string(right.owner->digest()));
        result.set("left_meaning", left.meaning);
        result.set("right_meaning", right.meaning);
        result.set("map", map);
        result.set(
            "scope",
            Json::string("Full field comparison under the declared measurement map and metric; the "
                         "map's biological validity is supplied by its evidence"));
    } else if (operation == "interval_coverage") {
        members(request, {"schema", "operation", "resource", "thresholds", "strand_policy"});
        const auto id = text(request.at("resource"));
        require(state_->physical.contains(id), "unknown interval resource");
        const auto item = state_->physical.at(id);
        const auto bytes = state_->artifact( item.at("original"));
        auto reservation = resources->reserve(1);
        result = interval_analysis(bytes, item.at("declaration"), request);
    } else if (operation == "sequence_motifs") {
        members(request, {"schema", "operation", "resource", "motifs", "strands"});
        const auto id = text(request.at("resource"));
        require(state_->physical.contains(id), "unknown sequence resource");
        const auto item = state_->physical.at(id);
        const auto bytes = state_->artifact( item.at("original"));
        auto reservation = resources->reserve(1);
        result = motif_analysis(bytes, item.at("declaration"), request);
    } else if (operation == "allocation_bounds") {
        members(request,
                {"schema", "operation", "resource", "targets", "candidate_search", "scale",
                 "scale_provenance"},
                {"kernel_certificate"});
        const auto id = text(request.at("resource"));
        require(state_->physical.contains(id), "unknown assignment resource");
        const auto item = state_->physical.at(id);
        require(item.at("declaration").at("format").text() == "salmon_eq",
                "assignment resource must be unweighted Salmon equivalence classes");
        const auto parsed = salmon_classes(state_->artifact( item.at("original")));
        auto allocation = Json::object();
        allocation.set("schema", Json::string("gadget.assay.allocations.v1"));
        allocation.set("operation", Json::string("allocation_bounds"));
        allocation.set("candidates", parsed.at("candidates"));
        allocation.set("classes", parsed.at("classes"));
        allocation.set("targets", request.at("targets"));
        if (request.has("kernel_certificate"))
            allocation.set("kernel_certificate", request.at("kernel_certificate"));
        allocation.set("contrasts", Json::array());
        allocation.set("class_partition", Json::string("disjoint_observations"));
        allocation.set("units", Json::string("assigned_fragment_count"));
        allocation.set("candidate_search", request.at("candidate_search"));
        auto evidence = Json::array();
        evidence.push(item.at("original").at("sha256"));
        allocation.set("evidence", evidence);
        auto samples = Json::array(), sample = Json::object();
        sample.set("id", request.at("resource"));
        sample.set("counts", parsed.at("counts"));
        sample.set("scale", request.at("scale"));
        sample.set("scale_provenance", request.at("scale_provenance"));
        samples.push(sample);
        allocation.set("samples", samples);
        result = analyze_allocations(allocation, resources).at("result");
    } else
        throw std::invalid_argument("unsupported assay workspace operation: " + operation);
    auto output = Json::object();
    output.set("schema", Json::string("gadget.assay.result.v1"));
    output.set("query", operation == "relational_gluing" ? Json::parse(request.dump()) : request);
    output.set("source_digest", Json::string(digest()));
    output.set("result", result);
    output.set("arithmetic", Json::string("rational"));
    output.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
    output.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    if (operation == "relational_gluing") {
        // The gluing adapter builds fresh native observations and uses detached
        // source pages/declarations. Capture the caller's query above, then seal
        // this already detached report without copying its whole result tree.
        output.set("certificate_sha256", Json::string("sha256:" + sha256(output.dump())));
        return output;
    }
    return seal(output, "certificate_sha256");
}
Json AssayWorkspace::materialize(const Json &request, const fs::path &destination,
                                 std::shared_ptr<CpuResources> resources) const {
    members(request, {"schema", "id", "source_namespace", "design", "recipe", "channel"},
            {"view", "component"});
    require(request.at("schema").text() == "gadget.assays.materialize.v1",
            "assay materialization schema");
    require(!fs::exists(destination), "derived assay destination exists");
    const auto recipe = request.at("recipe");
    const auto op = recipe.at("operation").text();
    require(op == "feature_counts" || op == "expression_allocations" || op == "variant_genotypes" ||
                op == "chromatin_contrast" || op == "controlled_coverage" ||
                op == "assay_sections" || op == "spatial_relations" || op == "aggregate_assay" ||
                op == "design_aggregation" || op == "exact_sequence_alignment" ||
                op == "sequence_assembly_graph" || op == "constrained_allocations",
            "unsupported materialization recipe");
    require(op == "assay_sections" ? request.has("view") && request.has("component")
                                   : !request.has("view") && !request.has("component"),
            "assay sections publication requires a view and component");
    if (!resources)
        resources = CpuResources::local();
    if(op=="controlled_coverage")return materialize_coverage(admit_record_view(request),destination,std::move(resources));
    const auto certificate = execute(recipe, resources), result = certificate.at("result");
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage =
        private_directory(fs::absolute(destination).parent_path(), ".gadget-derive-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    auto provenance = Json::array();
    std::size_t retained_inputs = 0;
    std::set<std::string> retained_hashes;
    const auto retain = [&](const std::string &bytes) {
        if (!retained_hashes.insert(sha256(bytes)).second)
            return;
        const auto path = "input-" + std::to_string(retained_inputs++);
        write_text(stage / path, bytes);
        provenance.push(file_record(path, bytes));
    };
    const auto lineage = certificate.dump() + "\n";
    retain(lineage);
    retain(request.dump() + "\n");
    auto descriptor = Json::object();
    descriptor.set("schema", Json::string("gadget.assay.lineage.v1"));
    descriptor.set("id", state_->manifest.at("id"));
    descriptor.set("design", state_->manifest.at("design"));
    descriptor.set("workspace_sha256", Json::string(digest()));
    auto source_records = Json::array();
    std::set<std::string> used_resources;
    const auto retain_resource = [&](const std::string &id) {
        if (!used_resources.insert(id).second)
            return;
        const auto source = state_->physical.at(id);
        retain(state_->artifact( source.at("original")));
        if (source.has("reference_original"))
            retain(state_->artifact( source.at("reference_original")));
        if (source.has("provenance"))
            for (std::size_t p = 0; p < source.at("provenance").size(); ++p)
                retain(state_->artifact(source.at("provenance").at(p)));
        source_records.push(source);
    };
    if (recipe.has("resource"))
        retain_resource(recipe.at("resource").text());
    if (recipe.has("annotation"))
        retain_resource(recipe.at("annotation").text());
    if (recipe.has("query_resource"))
        retain_resource(recipe.at("query_resource").text());
    if (recipe.has("reference_resource"))
        retain_resource(recipe.at("reference_resource").text());
    if (recipe.has("samples"))
        for (std::size_t i = 0; i < recipe.at("samples").size(); ++i)
            retain_resource(recipe.at("samples").at(i).at("resource").text());
    if (recipe.has("evidence_source"))
        for (const auto *name : {"annotation", "reference_resource"})
            if (recipe.at("evidence_source").has(name))
                retain_resource(recipe.at("evidence_source").at(name).text());
    descriptor.set("resources", source_records);
    auto assay = Json::object();
    std::string features_text = "id\n", observations_text = "id\n", matrix_text;
    std::vector<std::string> feature_ids, observation_ids;
    std::vector<std::vector<std::string>> rows;
    if (op == "spatial_relations") {
        const auto source = state_->assays.at(recipe.at("assay").text());
        require(source.at("declaration").at("carrier").text() == "observations",
                "spatial topology requires an observation-carrier assay");
        assay = Json::parse(source.at("declaration").dump());
        const auto originals = source.at("originals");
        for (const auto &name : originals.native().names()) {
            const auto bytes = state_->artifact( originals.at(name));
            retain(bytes);
            const auto target = std::string(name) == "matrix" ? "source" : std::string(name);
            if (target == "source" || target == "features" || target == "observations" ||
                target == "missing") {
                const auto path = "assay-" + target;
                write_text(stage / path, bytes);
                assay.set(target, file_record(path, bytes));
            }
        }
        descriptor.set("assay", source.at("declaration"));
        auto relations = Json::object();
        std::vector<Json> ordered;
        if (originals.has("relations")) {
            relations = Json::parse(state_->artifact( originals.at("relations")));
            for (std::size_t i = 0; i < relations.at("relations").size(); ++i)
                ordered.push_back(relations.at("relations").at(i));
        }
        for (std::size_t i = 0; i < result.at("relations").at("relations").size(); ++i)
            ordered.push_back(result.at("relations").at("relations").at(i));
        std::stable_sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) {
            return (a.has("grade") ? a.at("grade").integer() : 1) <
                   (b.has("grade") ? b.at("grade").integer() : 1);
        });
        auto all = Json::array();
        for (const auto &r : ordered)
            all.push(r);
        relations.set("relations", all);
        const auto bytes = relations.dump() + "\n";
        write_text(stage / "relations.json", bytes);
        assay.set("relations", file_record("relations.json", bytes));
    } else {
        std::string quantity, units;
        if (op == "assay_sections") {
            const auto view_id = request.at("view").text(),
                       component = request.at("component").text();
            require(component == "input" || component == "compatible" || component == "residual",
                    "assay section publication component");
            std::optional<std::size_t> selected;
            const auto views = recipe.at("views");
            auto originals = Json::array();
            std::set<std::string> seen;
            for (std::size_t v = 0; v < views.size(); ++v) {
                const auto view = views.at(v);
                if (view.at("id").text() == view_id)
                    selected = v;
                const auto source = state_->assays.at(view.at("assay").text());
                if (seen.insert(view.at("assay").text()).second) {
                    originals.push(source.at("declaration"));
                    const auto files = source.at("originals");
                    for (const auto &name : files.native().names()) {
                        retain(state_->artifact( files.at(name)));
                    }
                }
            }
            require(selected.has_value(), "published assay view not declared");
            descriptor.set("assays", originals);
            const auto sections = result.at("sections");
            const auto first = sections.at(0).at("fields").at(*selected);
            require(request.at("channel").text() == first.at("meaning").at("channel").text(),
                    "published assay view must retain its source channel");
            for (std::size_t f = 0; f < first.at("ids").size(); ++f) {
                feature_ids.push_back(first.at("ids").at(f).text());
                rows.emplace_back();
            }
            for (std::size_t s = 0; s < sections.size(); ++s) {
                observation_ids.push_back(sections.at(s).at("id").text());
                const auto values = sections.at(s).at("fields").at(*selected).at(component);
                for (std::size_t f = 0; f < values.size(); ++f)
                    rows[f].push_back(values.at(f).raw() ? values.at(f).text() : "NA");
            }
            quantity = component + "_section_of_" + first.at("meaning").at("quantity").text();
            units = first.at("meaning").at("units").text();
            auto metrics = Json::array(), metric = Json::object();
            metric.set("kind", Json::string("diagonal"));
            metric.set("values", views.at(*selected).at("metric").at("diagonal"));
            metrics.push(metric);
            assay.set("metrics", metrics);
        } else if (op == "expression_allocations") {
            require(result.at("expression").at("status").text() == "feasible",
                    "infeasible expression constraints have no derived bounds");
            const auto samples = result.at("expression").at("samples");
            observations_text = "id\tsample\tbound_endpoint\tscale\tscale_provenance_sha256\n";
            for (std::size_t s = 0; s < samples.size(); ++s)
                for (const auto *endpoint : {"lower", "upper"}) {
                    const auto id = "sample:" + std::to_string(s) + "." + endpoint;
                    observation_ids.push_back(id);
                    const auto declared = recipe.at("samples").at(s);
                    observations_text +=
                        id + "\t" + samples.at(s).at("id").text() + "\t" + endpoint + "\t" +
                        declared.at("scale").text() +
                        "\tsha256:" + sha256(declared.at("scale_provenance").text()) + "\n";
                }
            for (std::size_t t = 0; t < result.at("targets").size(); ++t) {
                feature_ids.push_back(result.at("targets").at(t).at("id").text());
                std::vector<std::string> row;
                for (std::size_t s = 0; s < samples.size(); ++s) {
                    const auto b = samples.at(s).at("bounds").at(t);
                    row.push_back(b.at("lower").text());
                    row.push_back(b.at("upper").text());
                }
                rows.push_back(std::move(row));
            }
            quantity = "allocation_expression_bound";
            units = result.at("units").text();
        } else if (op == "feature_counts") {
            observation_ids = {recipe.at("resource").text() + ".lower",
                               recipe.at("resource").text() + ".upper"};
            const auto bounds =
                result.at("allocation_certificate").at("result").at("samples").at(0).at("bounds");
            for (std::size_t i = 0; i < bounds.size(); ++i) {
                const auto b = bounds.at(i);
                feature_ids.push_back(b.at("id").text());
                rows.push_back({b.at("lower").text(), b.at("upper").text()});
            }
            observations_text = "id\tbound_endpoint\n" + observation_ids[0] + "\tlower\n" +
                                observation_ids[1] + "\tupper\n";
            quantity = "assigned_segment_count_bound";
            units = "segments";
        } else if (op == "variant_genotypes") {
            for (std::size_t i = 0; i < result.at("samples").size(); ++i)
                observation_ids.push_back(result.at("samples").at(i).text());
            for (std::size_t i = 0; i < result.at("allele_fields").size(); ++i) {
                const auto field = result.at("allele_fields").at(i);
                feature_ids.push_back(field.at("id").text());
                std::vector<std::string> row;
                for (std::size_t j = 0; j < field.at("values").size(); ++j) {
                    const auto value = field.at("values").at(j);
                    row.push_back(value.at("complete").boolean() ? value.at("copy_number").text()
                                                                 : "NA");
                }
                rows.push_back(row);
            }
            quantity = "called_allele_copy_number";
            units = "allele_copies";
        } else if (op == "controlled_coverage") {
            auto metrics = Json::array(), physical = Json::object();
            physical.set("kind", Json::string("diagonal"));
            physical.set("values", result.at("metric_diagonal"));
            metrics.push(physical);
            metrics.push(Json::parse("{\"kind\":\"identity\"}"));
            assay.set("metrics", metrics);
            features_text = "id\tcontig\tstart\tend\tstrand\n";
            observations_text = "id\tcontrast\trepresentation\n";
            for (std::size_t c = 0; c < result.at("contrasts").size(); ++c) {
                const auto contrast = result.at("contrasts").at(c).at("id").text();
                require(contrast.find_first_of("\t\r\n") == std::string::npos,
                        "contrast metadata contains separators");
                for (const auto *part : {"input", "represented", "residual"}) {
                    const auto id = "contrast:" + std::to_string(c) + ":" + part;
                    observation_ids.push_back(id);
                    observations_text += id + "\t" + contrast + "\t" + part + "\n";
                }
            }
            for (std::size_t i = 0; i < result.at("segments").size(); ++i) {
                const auto segment = result.at("segments").at(i);
                const auto id = segment.at("id").text();
                feature_ids.push_back(id);
                features_text += id + "\t" + segment.at("contig").text() + "\t" +
                                 std::to_string(segment.at("start").integer()) + "\t" +
                                 std::to_string(segment.at("end").integer()) + "\t" +
                                 segment.at("strand").text() + "\n";
                std::vector<std::string> row;
                for (std::size_t c = 0; c < result.at("contrasts").size(); ++c)
                    for (const auto *part : {"input", "represented", "residual"})
                        row.push_back(result.at("contrasts").at(c).at(part).at(i).text());
                rows.push_back(std::move(row));
            }
            quantity = result.at("quantity").text();
            units = result.at("units").text();
            const auto bytes = result.at("relation_declaration").dump() + "\n";
            write_text(stage / "relations.json", bytes);
            assay.set("relations", file_record("relations.json", bytes));
        } else if (op == "chromatin_contrast") {
            for (std::size_t i = 0; i < result.at("contrasts").size(); ++i)
                observation_ids.push_back(result.at("contrasts").at(i).at("id").text());
            features_text = "id\tcontig\tstart\tend\tstrand\n";
            for (std::size_t i = 0; i < result.at("segments").size(); ++i) {
                const auto s = result.at("segments").at(i);
                const auto id = "segment:" + std::to_string(i);
                feature_ids.push_back(id);
                features_text += id + "\t" + s.at("contig").text() + "\t" +
                                 std::to_string(s.at("start").integer()) + "\t" +
                                 std::to_string(s.at("end").integer()) + "\t" +
                                 s.at("strand").text() + "\n";
                std::vector<std::string> row;
                for (std::size_t j = 0; j < result.at("contrasts").size(); ++j)
                    row.push_back(result.at("contrasts").at(j).at("values").at(i).text());
                rows.push_back(row);
            }
            quantity = result.at("quantity").text();
            units = result.at("units").text();
        } else if (op == "exact_sequence_alignment") {
            features_text = "id\nminimum_prefix_cost\nminimum_suffix_cost\n";
            observations_text = "id\tquery_prefix\treference_prefix\n";
            feature_ids = {"minimum_prefix_cost", "minimum_suffix_cost"};
            rows = {{}, {}};
            for (std::size_t i = 0; i < result.at("prefix_suffix_witnesses").size(); ++i) {
                const auto p = result.at("prefix_suffix_witnesses").at(i);
                const auto q = std::to_string(p.at("query_prefix").integer()),
                           r = std::to_string(p.at("reference_prefix").integer());
                const auto id = "prefix:" + q + ":" + r;
                observation_ids.push_back(id);
                observations_text += id + "\t" + q + "\t" + r + "\n";
                rows[0].push_back(p.at("minimum_prefix_cost").text());
                rows[1].push_back(p.at("minimum_suffix_cost").text());
            }
            quantity = "declared_minimum_edit_cost";
            units = recipe.at("scoring").at("units").text();
            const auto bytes = result.at("relation_declaration").dump() + "\n";
            write_text(stage / "relations.json", bytes);
            assay.set("relations", file_record("relations.json", bytes));
        } else if (op == "constrained_allocations") {
            require(result.at("status").text() == "feasible",
                    "infeasible allocation has no derived measured bounds");
            observation_ids = {"feasible_lower", "feasible_upper", "optimal_lower",
                               "optimal_upper"};
            for (std::size_t i = 0; i < result.at("bounds").size(); ++i) {
                const auto b = result.at("bounds").at(i);
                feature_ids.push_back(b.at("id").text());
                rows.push_back({b.at("lower").text(), b.at("upper").text(),
                                b.at("optimal_lower").text(), b.at("optimal_upper").text()});
            }
            quantity = "constrained_allocation_bound";
            units = recipe.at("units").text();
        } else {
            for (std::size_t i = 0; i < result.at("feature_ids").size(); ++i)
                feature_ids.push_back(result.at("feature_ids").at(i).text());
            for (std::size_t i = 0; i < result.at("observation_ids").size(); ++i)
                observation_ids.push_back(result.at("observation_ids").at(i).text());
            for (std::size_t i = 0; i < result.at("rows").size(); ++i) {
                std::vector<std::string> row;
                for (std::size_t j = 0; j < result.at("rows").at(i).size(); ++j)
                    row.push_back(result.at("rows").at(i).at(j).text());
                rows.push_back(row);
            }
            quantity = result.at("quantity").text();
            units = result.at("units").text();
            if (op == "sequence_assembly_graph") {
                const auto bytes = result.at("relation_declaration").dump() + "\n";
                write_text(stage / "relations.json", bytes);
                assay.set("relations", file_record("relations.json", bytes));
            } else {
                const auto source = state_->assays.at(recipe.at("assay").text());
                descriptor.set("assay", source.at("declaration"));
                const auto originals = source.at("originals");
                for (const auto &name : originals.native().names()) {
                    retain(state_->artifact( originals.at(name)));
                }
            }
        }
        require(!feature_ids.empty() && !observation_ids.empty(),
                "derived assay has no measured axes");
        const auto token = [](const std::string &s) {
            require(!s.empty() && s.find_first_of("\t\r\n") == std::string::npos,
                    "derived TSV identity contains separators");
        };
        for (const auto &id : feature_ids)
            token(id);
        if (features_text == "id\n")
            for (const auto &id : feature_ids)
                features_text += id + "\n";
        if (observations_text == "id\n")
            for (const auto &id : observation_ids) {
                token(id);
                observations_text += id + "\n";
            }
        matrix_text = "id";
        for (const auto &id : observation_ids) {
            token(id);
            matrix_text += "\t" + id;
        }
        matrix_text += "\n";
        require(rows.size() == feature_ids.size(), "derived matrix feature count");
        for (std::size_t i = 0; i < rows.size(); ++i) {
            require(rows[i].size() == observation_ids.size(), "derived matrix observation count");
            matrix_text += feature_ids[i];
            for (const auto &v : rows[i])
                matrix_text += "\t" + v;
            matrix_text += "\n";
        }
        for (const auto &[name, bytes] :
             std::map<std::string, std::string>{{"source", matrix_text},
                                                {"features", features_text},
                                                {"observations", observations_text}}) {
            write_text(stage / name, bytes);
            assay.set(name, file_record(name, bytes));
        }
        assay.set("format", Json::string("dense_tsv"));
        assay.set("carrier",
                  Json::string(op == "sequence_assembly_graph" ? "observations" : "features"));
        if (op == "exact_sequence_alignment")
            assay.set("carrier", Json::string("observations"));
        assay.set("arithmetic", Json::string("exact_rational"));
        assay.set("implicit_entries", Json::string("unobserved"));
        assay.set("missing_tokens", Json::parse("[\"NA\"]"));
        assay.set("quantity", Json::string(quantity));
        assay.set("units", Json::string(units));
        assay.set("kind", Json::string(op));
        if (result.has("reference"))
            assay.set("reference", result.at("reference"));
    }
    retain(descriptor.dump() + "\n");
    assay.set("id", request.at("id"));
    assay.set("source_namespace", request.at("source_namespace"));
    assay.set("channel", request.at("channel"));
    assay.set(
        "method",
        Json::string("GADGET " + op +
                     "; full generating certificate and source files retained in provenance"));
    assay.set("provenance", provenance);
    auto import = Json::object(), assays = Json::array();
    assays.push(assay);
    import.set("schema", Json::string("gadget.assays.import.v1"));
    import.set("id", request.at("id"));
    import.set("design", request.at("design"));
    import.set("assays", assays);
    import.set("resources", Json::array());
    import.set("links", Json::array());
    write_text(stage / "import.json", import.dump() + "\n");
    import_assays(stage / "import.json", destination);
    return AssayWorkspace(destination / "workspace.rfcxs").catalog();
}
Json AssayWorkspace::replay_derivation(const std::string &assay_id, const fs::path &destination,
                                       std::shared_ptr<CpuResources> resources) const {
    require(!fs::exists(destination), "derivation replay destination exists");
    require(state_->assays.contains(assay_id), "unknown derived assay");
    const auto item = state_->assays.at(assay_id);
    if(item.at("declaration").has("schema") && item.at("declaration").at("schema").text()=="gadget.assay.native.v1")return replay_native_derivation(assay_id,destination,std::move(resources));
    require(item.at("declaration").has("provenance"), "assay has no retained generating lineage");
    std::map<std::string, Json> artifacts;
    Json request, expected, lineage;
    const auto originals = item.at("originals");
    for (std::size_t i = 0; i < item.at("declaration").at("provenance").size(); ++i) {
        const auto record = originals.at("provenance-" + std::to_string(i));
        artifacts.emplace(record.at("sha256").text(), record);
        const auto bytes = state_->artifact( record);
        if (i < 2) {
            const auto value = Json::parse(bytes);
            if (value.has("schema") && value.at("schema").text() == "gadget.assay.result.v1")
                expected = value;
            if (value.has("schema") && value.at("schema").text() == "gadget.assays.materialize.v1")
                request = value;
        }
        if (bytes.starts_with("{")) {
            try {
                const auto value = Json::parse(bytes);
                if (value.has("schema") && value.at("schema").text() == "gadget.assay.lineage.v1")
                    lineage = value;
            } catch (const std::invalid_argument &) {
            }
        }
    }
    require(request.raw() && expected.raw() && lineage.raw(),
            "incomplete or unsupported derivation lineage");
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage =
        private_directory(fs::absolute(destination).parent_path(), ".gadget-replay-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    std::map<std::string, Json> staged;
    const auto restore = [&](const Json &record) {
        const auto hash = record.at("sha256").text();
        require(artifacts.contains(hash), "generating input absent from retained lineage");
        if (!staged.contains(hash)) {
            const auto name = "source-" + std::to_string(staged.size());
            const auto bytes = state_->artifact( artifacts.at(hash));
            write_text(stage / name, bytes);
            staged.emplace(hash, file_record(name, bytes));
        }
        return staged.at(hash);
    };
    auto import = Json::object(), assays = Json::array(), images = Json::array();
    const auto restore_assay = [&](const Json &declaration) {
        auto source = Json::parse(declaration.dump());
        for (const auto *name : {"source", "features", "observations", "relations", "missing"})
            if (source.has(name))
                source.set(name, restore(source.at(name)));
        if (source.has("provenance")) {
            auto records = Json::array();
            for (std::size_t i = 0; i < source.at("provenance").size(); ++i)
                records.push(restore(source.at("provenance").at(i)));
            source.set("provenance", records);
        }
        assays.push(source);
    };
    if (lineage.has("assay"))
        restore_assay(lineage.at("assay"));
    if (lineage.has("assays"))
        for (std::size_t i = 0; i < lineage.at("assays").size(); ++i)
            restore_assay(lineage.at("assays").at(i));
    for (std::size_t i = 0; i < lineage.at("resources").size(); ++i) {
        auto source = Json::parse(lineage.at("resources").at(i).at("declaration").dump());
        source.set("source", restore(source.at("source")));
        if (source.has("reference_source"))
            source.set("reference_source", restore(source.at("reference_source")));
        if (source.has("provenance")) {
            auto originals = Json::array();
            for (std::size_t p = 0; p < source.at("provenance").size(); ++p)
                originals.push(restore(source.at("provenance").at(p)));
            source.set("provenance", originals);
        }
        images.push(source);
    }
    import.set("schema", Json::string("gadget.assays.import.v1"));
    import.set("id", lineage.at("id"));
    import.set("design", lineage.at("design"));
    import.set("assays", assays);
    import.set("resources", images);
    import.set("links", Json::array());
    Json catalog;
    if (assays.size() == 0 && images.size() == 0) {
        require(request.at("recipe").at("operation").text() == "constrained_allocations",
                "resource-free replay requires a fully declared constraint system");
        // Every generating input is in the retained request. A new empty source
        // workspace would violate admission and introduce a fictitious assay.
        catalog = materialize(request, destination, resources);
    } else {
        write_text(stage / "import.json", import.dump() + "\n");
        import_assays(stage / "import.json", stage / "workspace");
        catalog = AssayWorkspace(stage / "workspace/workspace.rfcxs")
                      .materialize(request, destination, resources);
    }
    try {
        const auto produced = catalog.at("assays").at(0).at("originals").at("provenance-0");
        const auto actual = Json::parse(local_artifact(fs::canonical(destination), produced));
        std::function<Json(const Json &)> payload = [&](const Json &value) -> Json {
            if (json_object_is_type(value.raw(), json_type_array)) {
                auto result = Json::array();
                for (std::size_t i = 0; i < value.size(); ++i)
                    result.push(payload(value.at(i)));
                return result;
            }
            if (json_object_is_type(value.raw(), json_type_object)) {
                auto result = Json::object();
                std::vector<std::string> keys;
                json_object_object_foreach(value.raw(), name, entry) {
                    (void)entry;
                    keys.emplace_back(name);
                }
                std::sort(keys.begin(), keys.end());
                for (const auto &key : keys) {
                    // Compiler/source bindings and native scheduling receipts
                    // remain in both full certificates. They are not numerical
                    // coordinates, units, support, relations or ordered history.
                    if (key == "certificate_sha256" || key == "gadget_sha256" ||
                        key == "gadget_sources_sha256" || key == "rexgraph_sha256" ||
                        key == "native_sha256" || key == "source_digest")
                        continue;
                    const auto child = value.at(key);
                    if (key == "execution" && json_object_is_type(child.raw(), json_type_object) &&
                        child.has("owner") &&
                        json_object_is_type(child.at("owner").raw(), json_type_string) &&
                        child.at("owner").text().starts_with("rfc::"))
                        continue;
                    result.set(key, payload(child));
                }
                return result;
            }
            return value;
        };
        require(payload(actual.at("result")).dump() == payload(expected.at("result")).dump(),
                "replayed generating result differs from retained mathematical payload");
        auto receipt = Json::object();
        receipt.set("schema", Json::string("gadget.assays.derivation.replay.v1"));
        receipt.set("status", Json::string("PASS"));
        receipt.set("source_digest", Json::string(digest()));
        receipt.set("expected_certificate", expected.at("certificate_sha256"));
        receipt.set("actual_certificate", actual.at("certificate_sha256"));
        receipt.set("comparison_profile", Json::string("mathematical_payload_v2"));
        receipt.set("result_payload_sha256",
                    Json::string("sha256:" + sha256(payload(actual.at("result")).dump())));
        receipt.set("workspace", catalog);
        return seal(receipt, "certificate_sha256");
    } catch (...) {
        fs::remove_all(destination);
        throw;
    }
}
NativeRecord AssayWorkspace::execute_coverage(const NativeRecord &request,
    std::shared_ptr<CpuResources> resources,std::optional<std::uint64_t> source_namespace) const {
    const auto query=render_record_view(request);
        using N=NativeRecordValue;
        require(query.at("schema").text()=="gadget.assays.query.v1","assay query schema");
        if(!resources)resources=CpuResources::local();
        auto reservation=resources->reserve(1);
        ResourceImages images;
        const auto samples=query.at("samples");
        for(std::size_t s=0;s<samples.size();++s) {
            const auto id=text(samples.at(s).at("resource"));
            require(state_->physical.contains(id),"unknown physical resource");
            if(images.contains(id))continue;
            const auto item=state_->physical.at(id);
            images.emplace(id,ResourceImage{item.at("declaration"),state_->artifact(item.at("original")),{}});
        }
        auto fields=controlled_coverage_source(images,query);
        std::uint64_t key=1;
        for(const auto &source:fields.sources) {
            require(source.key()!=std::numeric_limits<std::uint64_t>::max(),"raw source identity space exhausted");
            key=std::max(key,source.key()+1);
        }
        auto owner=native_publication().archive_source(key);
        const auto size=owner.bytes().size();
        auto certificate=N(N::Object{{
            {"schema",N(std::string("gadget.assay.result.native.v1"))},
            {"query",import_record_view(query)}, {"source_digest",N(digest())},
            {"result",std::move(fields.declaration)},
            {"workspace_publication",N(N::ByteSpan{key,0,size,false})},
            {"arithmetic",N(std::string("rational"))},
            {"gadget_sha256",N(std::string(GADGET_SOURCE_SHA256))},
            {"rexgraph_sha256",N(std::string(GADGET_REXGRAPH_SHA256))}
        }});
        fields.sources.push_back(std::move(owner));
        return retain_operation(std::move(fields.endpoint),source_namespace.value_or(new_native_record_namespace()),std::move(fields.sources),request,std::move(certificate));
}
NativeRecord AssayWorkspace::execute_record(const NativeRecord &request,
    std::shared_ptr<CpuResources> resources) const {
    const auto query = render_record_view(request);
    if (query.has("operation") && query.at("operation").text() == "controlled_coverage")
        return execute_coverage(request,std::move(resources));
    if (query.has("operation") && (query.at("operation").text() == "alignment_evidence" ||
                                   query.at("operation").text() == "feature_counts")) {
        using N = NativeRecordValue;
        require(query.at("schema").text() == "gadget.assays.query.v1", "assay query schema");
        if (!resources) resources = CpuResources::local();
        const bool counting = query.at("operation").text() == "feature_counts";
        ResourceImages images;
        const auto load = [&](const char *name) {
            const auto id = text(query.at(name));
            require(state_->physical.contains(id), "unknown physical resource");
            const auto item = state_->physical.at(id);
            images.emplace(id, ResourceImage{item.at("declaration"), state_->artifact(item.at("original")),
                item.has("reference_original") ? state_->artifact(item.at("reference_original")) : std::string{}});
        };
        load("resource");
        if (counting) load("annotation");
        auto fields = [&] {
            if (counting) return feature_count_fields_source(images, query, resources);
            auto reservation = resources->reserve(1);
            return alignment_fields_source(images.at(query.at("resource").text()), query);
        }();
        std::uint64_t key = 1;
        for (const auto &source : fields.sources) {
            require(source.key() != std::numeric_limits<std::uint64_t>::max(), "raw source identity space exhausted");
            key = std::max(key, source.key()+1);
        }
        auto owner = native_publication().archive_source(key);
        const auto bytes = owner.bytes().size();
        fields.sources.push_back(std::move(owner));
        N certificate(N::Object{{
            {"schema", N(std::string("gadget.assay.result.native.v1"))},
            {"query", import_record_view(query)},
            {"source_digest", N(digest())},
            {"result", std::move(fields.declaration)},
            {"workspace_publication", N(N::ByteSpan{key, 0, bytes, false})},
            {"arithmetic", N(std::string("rational"))},
            {"gadget_sha256", N(std::string(GADGET_SOURCE_SHA256))},
            {"rexgraph_sha256", N(std::string(GADGET_REXGRAPH_SHA256))}
        }});
        return retain_operation(std::move(fields.endpoint), new_native_record_namespace(),
                                std::move(fields.sources), request, std::move(certificate));
    }
    if (query.has("operation") && query.at("operation").text() == "read_quality") {
        using N = NativeRecordValue;
        require(query.at("schema").text() == "gadget.assays.query.v1", "assay query schema");
        const auto id = text(query.at("resource"));
        require(state_->physical.contains(id), "unknown physical resource");
        if (!resources) resources = CpuResources::local();
        auto reservation = resources->reserve(1);
        const auto resource = state_->physical.at(id);
        auto fields = quality_fields(ResourceImage{resource.at("declaration"),
                                     state_->artifact(resource.at("original")), {}}, query);
        const auto typed = NativeRecord::create(fields, new_native_record_namespace());
        auto certificate = Json::object();
        certificate.set("schema", Json::string("gadget.assay.result.v1"));
        certificate.set("query", query);
        certificate.set("source_digest", Json::string(digest()));
        certificate.set("result", quality_fields_view(typed));
        certificate.set("arithmetic", Json::string("rational"));
        certificate.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
        certificate.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        certificate.set("certificate_sha256", Json::string("sha256:" + sha256(certificate.dump())));
        // Keep the certificate's original member order without admitting its
        // numerical display arrays as a second representation.
        certificate.set("result", Json{});
        auto retained = import_record_view(certificate);
        for (auto &[name, value] : std::get<N::Object>(retained.data).members)
            if (name == "result") value = std::move(fields);
        return retain_operation(native_publication(), request, std::move(retained));
    }
    return retain_operation(native_publication(),request,
        import_record_view(execute(query,std::move(resources))));
}
Json AssayWorkspace::run(const Json &request,const fs::path &destination,
    std::shared_ptr<CpuResources> resources) const {
    return render_record_view(run_record(admit_record_view(request),destination,
        std::move(resources)).at("result"));
}
NativeRecord AssayWorkspace::run_record(const NativeRecord &retained_request,
    const fs::path &destination,std::shared_ptr<CpuResources> resources) const {
    const auto request=render_record_view(retained_request);
    members(request, {"schema", "operation", "assay", "query"}, {"partitions", "concurrent_jobs"});
    require(request.at("schema").text() == "gadget.assays.query.v1" &&
                request.at("operation").text() == "native_query",
            "partitioned run requires a native query");
    require(!fs::exists(destination), "assay result destination already exists");
    if (!resources)
        resources = CpuResources::local();
    const auto id = text(request.at("assay"));
    require(state_->assays.contains(id), "unknown assay");
    const auto packets = state_->assays.at(id).at("partitions");
    std::vector<std::size_t> selected;
    if (request.has("partitions")) {
        std::set<std::size_t> seen;
        for (std::size_t i = 0; i < request.at("partitions").size(); ++i) {
            const auto p =
                count(request.at("partitions").at(i), packets.size() - 1, "assay partition");
            require(seen.insert(p).second, "duplicate partition");
            selected.push_back(p);
        }
        require(!selected.empty(), "empty assay selection");
    } else
        for (std::size_t i = 0; i < packets.size(); ++i)
            selected.push_back(i);
    auto leaf = Json::parse(request.dump());
    if (leaf.has("concurrent_jobs")) {
        auto stripped = Json::object();
        json_object_object_foreach(leaf.raw(), name, value) {
            (void)value;
            if (std::string(name) != "concurrent_jobs")
                stripped.set(name, leaf.at(name));
        }
        leaf = stripped;
    }
    // Use the identical query admission/partition policy as execute().
    const auto op = request.at("query").at("operation").text();
    require(packets.size() == 1 || op == "paired_field_contrasts" || op == "field_summary" ||
                op == "coverage" || op == "tensor_window" || op == "inspect" || op == "provenance",
            "whole-carrier operation cannot be split across storage partitions");
    for (auto p : selected) {
        ScientificSession session(study(id, p), resources);
        session.validate(request.at("query"));
    }
    const auto q = request.at("query");
    const auto workers =
        q.has("workers")
            ? count(q.at("workers"), resources->coordinator().capacity(), "native query workers")
            : 1;
    require(workers > 0 && workers <= resources->coordinator().capacity(),
            "assay job worker demand");
    const auto requested_concurrency =
        request.has("concurrent_jobs")
            ? count(request.at("concurrent_jobs"), std::numeric_limits<std::size_t>::max(),
                    "concurrent assay jobs")
            : 1;
    require(requested_concurrency > 0, "assay concurrency must be positive");
    const auto concurrency = std::min(
        {requested_concurrency, selected.size(), resources->coordinator().capacity() / workers});
    fs::create_directories(fs::absolute(destination).parent_path());
    require(fs::create_directory(destination), "assay results creation conflict");
    auto entries = Json::array();
    try {
        // Joined jobs write their own certificates and retain only receipts.
        // The actual assigned concurrency determines each wave's size.
        for (std::size_t begin = 0; begin < selected.size();) {
            const auto end = begin + std::min(concurrency, selected.size() - begin);
            std::vector<rfc::runtime::WorkJob> jobs;
            std::vector<rfc::runtime::WorkKey> keys;
            const auto buffer = resources->coordinator().policy().declared_buffer_byte_limit /
                                std::min(concurrency, end - begin);
            for (auto i = begin; i < end; ++i) {
                const auto p = selected[i];
                const rfc::runtime::WorkKey key{i + 1, 1};
                keys.push_back(key);
                jobs.push_back(rfc::runtime::WorkJob::make<RetainedRecord>(
                    key, {workers, buffer, sizeof(RetainedRecord)}, {},
                    [this, p, id, q, resources, destination](const auto &, auto &assigned,
                                                             std::stop_token stop) {
                        if (stop.stop_requested())
                            throw rfc::runtime::ExecutionCancelled{};
                        ScientificSession session(study(id, p), resources->assigned(assigned));
                        const auto native_query=admit_record_view(q);
                        const auto native_answer=session.execute_record(native_query);
                        native_answer.write(destination/("partition-"+std::to_string(p)+".rfcxs"));
                        const auto answer=render_record_view(native_answer.at("result"));
                        const auto bytes = answer.dump() + "\n";
                        const auto path = "partition-" + std::to_string(p) + ".json";
                        write_text(destination / path, bytes);
                        const auto native_path="partition-"+std::to_string(p)+".rfcxs";
                        auto receipt = file_record(native_path, read_text(destination/native_path));
                        receipt.set("view_path",Json::string(path));
                        receipt.set("partition", cardinality(p));
                        receipt.set("source_digest", answer.at("source_digest"));
                        receipt.set("certificate_sha256", answer.at("certificate_sha256"));
                        receipt.set("native_path",Json::string(native_path));
                        return RetainedRecord(receipt);
                    }));
            }
            rfc::runtime::WorkGraphLimits limits;
            limits.jobs = jobs.size();
            limits.edges = 0;
            limits.depth = 1;
            limits.in_flight = concurrency;
            limits.retained_result_bytes = rfc::machine::checked_size_mul(jobs.size(),sizeof(RetainedRecord));
            const rfc::runtime::WorkGraph graph(std::move(jobs), limits);
            auto work = graph.submit(resources->coordinator());
            work.wait();
            for (std::size_t j = 0; j < keys.size(); ++j) {
                const auto completion = work.completion(keys[j]);
                const auto answer = completion.result<RetainedRecord>();
                const auto p = selected[begin + j];
                auto item = Json::parse(answer->dump());
                const auto receipt = completion.receipt();
                item.set("workers", cardinality(receipt.demand.workers));
                item.set("execution_ns", Json::string(std::to_string(receipt.execution_ns)));
                item.set("queue_ns", Json::string(std::to_string(receipt.queue_ns)));
                if (packets.at(p).has("carrier_ids"))
                    item.set("carrier_ids", packets.at(p).at("carrier_ids"));
                else
                    item.set("source_index", study(id, p)->source_index());
                entries.push(item);
            }
            begin = end;
        }
        auto result = Json::object();
        result.set("schema", Json::string("gadget.assays.run.v1"));
        result.set("source_digest", Json::string(digest()));
        result.set("query", request);
        result.set("results", entries);
        result.set("concurrent_jobs", cardinality(concurrency));
        result.set(
            "retention",
            Json::string(
                "certificate files written by joined native jobs; artifact receipts retained"));
        result.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
        result.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        result = seal(result, "certificate_sha256");
        auto publication=retain_operation(native_publication(),retained_request,import_record_view(result));
        publication.write(destination / "results.rfcxs");
        write_text(destination / "results.json", result.dump() + "\n");
        return publication;
    } catch (...) {
        fs::remove_all(destination);
        throw;
    }
}
Json bioinformatics_capabilities() {
    auto result = Json::parse(
        R"({"schema":"gadget.bioinformatics.capabilities.v1","matrix_formats":["dense_tsv","matrix_market"],"physical_formats":["bed","gtf","fasta","fastq"],"workspace_methods":["native_query","interval_coverage"],"allocation_method":"allocation_bounds","existing_native_methods":["paired_field_contrasts","controlled_pairings","trajectory_structure","structure_transition","local_sections","green_response","tensor_window","field_summary","boundary","adjoint","hodge","typed_moment"],"unavailable_services":["streaming_multigigabase_gapped_alignment","gapped_splice_alignment","error_correcting_consensus_assembly","unbounded_integer_optimization","population_significance_without_a_declared_null_model","automatic_cell_type_or_causal_inference","integral_or_zigzag_persistence"],"limits":{}})");
#if GADGET_HAS_HDF5
    auto formats = result.at("matrix_formats");
    formats.push(Json::string("h5ad"));
    result.set("h5ad_encodings", Json::parse("[\"dense\",\"csr_matrix\",\"csc_matrix\"]"));
#endif
    result.at("workspace_methods").push(Json::string("field_transfer"));
    result.at("workspace_methods").push(Json::string("assay_sections"));
    result.at("workspace_methods").push(Json::string("sequence_motifs"));
    result.at("workspace_methods").push(Json::string("allocation_bounds"));
    for (const auto *op : {"chain_certificate", "realized_path_comparison"})
        result.at("existing_native_methods").push(Json::string(op));
    result.at("physical_formats").push(Json::string("salmon_eq"));
    result.at("physical_formats").push(Json::string("gff3"));
    for (const auto *op :
         {"read_quality", "annotation_features", "chromatin_contrast", "controlled_coverage",
          "spatial_relations", "aggregate_assay", "design_aggregation", "exact_sequence_alignment",
          "reference_sequence_matches", "sequence_assembly_graph", "constrained_allocations",
          "relational_persistence", "indexed_sequence_alignment", "sequence_splice_discovery",
          "read_overlap_reconstruction", "expression_allocations", "relational_gluing", "query_batch"})
        result.at("workspace_methods").push(Json::string(op));
#if GADGET_HAS_HTSLIB
    for (const auto *format : {"sam", "bam", "cram", "vcf", "bcf"})
        result.at("physical_formats").push(Json::string(format));
    for (const auto *op : {"alignment_evidence", "feature_counts", "variant_genotypes"})
        result.at("workspace_methods").push(Json::string(op));
#endif
    result.set("htslib_enabled", Json::boolean(GADGET_HAS_HTSLIB));
    result.set(
        "source_admission",
        Json::parse(
            R"({"matrix_shape":"checked native IndexDomain; no application row, column or logical-entry ceiling","files":"complete retained input within string/stream representation","source_records":"complete decoded records and original bytes","relational_carriers":"whole authored carrier with native cell/index representation","structural_certificates":"complete native source; no application cell/path/level ceiling","query_jobs":"declared job count and observed native CPU capacity; owning file receipts","memory":"observed headroom supplies cooperative buffer declarations; not a GMP/process heap quota"})"));
    result.set("materialization", Json::string("derive-assays"));
    result.set("workspace_binding",
               Json::string("bind_assays pipeline recipe; explicit selected inputs and collection "
                            "aliases, unchanged native packets"));
    result.at("existing_native_methods").push(Json::string("green_action"));
    result.set("pipeline_profiles",
               Json::parse("[\"rna_seq\",\"single_cell\",\"atac_seq\",\"chip_seq\",\"spatial\","
                           "\"multi_omics\",\"variants\",\"sequence\",\"custom\"]"));
    result.set(
        "pair_alignment_execution",
        Json::parse(
            R"({"carrier":"Complete source-sized edit grid with native cell/index representation","witnesses":"Native minimum routes and reusable subpaths; complete prefix/suffix dual fields and all minimizing steps","history":"Every native witness route retained once with each coordinate's route and position","arithmetic":"Exact rational costs and native path/field actions"})"));
    result.set(
        "reference_matching_execution",
        Json::parse(
            R"({"literal_candidates":"Existing injective DNA seed index; disjoint-seed completeness and full native mismatch verification","IUPAC_and_junctions":"Complete original compatibility and declared-junction searches","output":"Every accepted placement, original query and source/junction identity","shape":"Declared source counts; exact counters with decimal strings beyond JSON int64"})"));
    result.set(
        "allocation_execution",
        Json::parse(
            R"({"bounds":"Sparse declared candidate-to-target memberships and native exact count actions","kernel":"Native class kernel frames joined through common target rows; complete rank and nullity","kernel_certificate":"Optional source-ordered kernel factors and target-image entries","shape":"Declared source counts and intrinsic native cell/index representation"})"));
    result.set(
        "indexed_sequence_execution",
        Json::parse(
            R"({"edit_budget":"Declared nonnegative budget strictly below each query length","index":"Complete requested literal seeds with injective uint64 two-bit keys","routes":"Complete band carriers and native dijkstra_exact endpoints","output":"All placements and one minimizing ordered realization per placement","shape":"Source-derived coordinates, checked counters and intrinsic native cell representation"})"));
    result.set(
        "sequence_profiles",
        Json::parse(
            R"({"indexed_sequence_alignment":"Complete unit-edit substring search for literal ACGT reads, with an injective DNA seed index and native exact routes","sequence_splice_discovery":"Complete exact exon chains with sequence-derived junctions, declared motif/length/junction constraints and complete source-sized output","read_overlap_reconstruction":"Every minimum suffix-prefix cover of all declared source reads, with native categorical section gluing and all alternative spellings"})"));
    result.set(
        "composite_sequence_execution",
        Json::parse(
            R"({"splice":"Complete exon extensions from intrinsic-sized literal seeds; exact gap interval joins and heap-based complete chain enumeration","read_cover":"Native impossible-start reachability witnesses, source-sized PackedBits subset keys, complete reachable prefix carrier and certified native suffix pruning","gluing":"Native categorical complete-overlap gluing and every original read restriction","comparison":"Native inherited Gram on every distinct-spelling representative; exact common endpoint certificates","assembly":"Complete original observed k-mer occurrence map, source walks and native boundary rank","shape":"Declared source/output counts and checked intrinsic native cell, identity and codec domains; no dataset/result policy ceilings"})"));
    result.set(
        "expression_allocation_execution",
        Json::parse(
            R"({"sources":["transcript_sequences","annotation_alignments","equivalence_classes"],"maps":"Declared shared candidate dictionary and complete transcript/feature memberships","outputs":"Sharp scaled sample and contrast bounds, complete evidence and native ambiguity kernels","coupling":"Authored constraints on unscaled target counts; existing native affine allocation solver with caller-selected complete witnesses","replay":"Portable generating source files and native endpoint bound fields"})"));
    result.set(
        "constrained_allocation_execution",
        Json::parse(
            R"({"coordinates":"Native affine equality kernel retaining every source variable","independence":"Native contrast-kernel indicators over reduced inequality and bound supports","search":"Complete local integer boxes or distinct rational active bases","witnesses":"Explicit complete Cartesian output by default; caller-requested factorized complete output","counts":"Exact nonnegative integers; decimal strings beyond signed 64-bit JSON representation","shape":"Declared source counts and intrinsic native cell/index representation"})"));
    result.set(
        "physical_geometry",
        Json::parse(
            R"({"intervals":"complete endpoint union and every source occurrence; sparse native sample actions", "controls":"declared multi-sample linear contrasts and optional intercept; native prepared Gram solve with reference-length metric", "spatial":"complete exact positive-diagonal neighborhoods using sorted-coordinate exclusion; original pair and occurrence order", "limits":"source counts and intrinsic native representation; no application sample, atom, membership, coordinate or relation ceiling"})"));
    result.set("assay_routes", Json::parse(R"({
        "bulk_expression":["feature_counts","allocation_bounds","paired_field_contrasts","controlled_pairings"],
        "splicing_evidence":["reference_sequence_matches","alignment_evidence","annotation_features"],
        "single_cell":["aggregate_assay","assay_sections","local_sections","controlled_pairings"],
        "chromatin":["interval_coverage","chromatin_contrast","controlled_coverage","sequence_motifs"],
        "spatial":["spatial_relations","green_response","hodge","local_sections"],
        "multi_omics":["field_transfer","assay_sections","aggregate_assay","typed_moment","local_sections"],
        "variant_calls":["variant_genotypes"],
        "sequence_validation":["read_quality","exact_sequence_alignment","reference_sequence_matches"],
        "assembly":["sequence_assembly_graph","hodge","green_response","relational_persistence"],
        "constrained_assignment":["constrained_allocations"],
        "persistent_structure":["relational_persistence","structure_transition"],
        "longitudinal_assays":["trajectory_structure","structure_transition"]
    })"));
    result.at("assay_routes")
        .at("sequence_validation")
        .push(Json::string("indexed_sequence_alignment"));
    result.at("assay_routes")
        .at("splicing_evidence")
        .push(Json::string("sequence_splice_discovery"));
    result.at("assay_routes").at("assembly").push(Json::string("read_overlap_reconstruction"));
#if !GADGET_HAS_HTSLIB
    result.at("assay_routes")
        .set("splicing_evidence", Json::parse("[\"annotation_features\",\"reference_sequence_"
                                              "matches\",\"sequence_splice_discovery\"]"));
    result.at("assay_routes").set("variant_calls", Json::array());
    result.at("assay_routes")
        .set("bulk_expression",
             Json::parse(
                 "[\"allocation_bounds\",\"paired_field_contrasts\",\"controlled_pairings\"]"));
#endif
    result.at("assay_routes").at("bulk_expression").push(Json::string("expression_allocations"));
    result.set("hdf5_enabled", Json::boolean(GADGET_HAS_HDF5));
    return result;
}
} // namespace gadget
