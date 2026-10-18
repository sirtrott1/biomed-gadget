#include <gadget/build.hpp>
#include "assay_detail.hpp"
#include <gadget/record_view.hpp>
#include <rfc/storage/file_io.hpp>
#include <rfc/storage/rfcx.hpp>

namespace gadget {
namespace {
using namespace assay_detail;
using namespace rfc::storage;
using N = NativeRecordValue;
N string(std::string_view value) { return N(std::string(value)); }
std::string field_key(const NativeRecord &key) { return key.integer().get_str(); }
N span_reference(const std::string &id, std::size_t bytes) {
    return N(N::Object{{{"document", N(id)},
                        {"start", N(rfc::numeric::Integer(0))},
                        {"end", N(rfc::numeric::Integer(std::to_string(bytes)))}}});
}
std::map<std::uint64_t, std::uint64_t> retained_request_sources(const NativeRecord &producer,
                                                             const NativeRecord &request) {
    std::map<std::uint64_t, std::uint64_t> bindings;
    const auto retained = producer.raw_sources();
    for (const auto &source : request.raw_sources()) {
        const auto found = std::find_if(retained.begin(), retained.end(),
                                       [&](const auto &item) { return item.shares_owner(source); });
        require(found != retained.end(), "generating request byte owner was lost");
        bindings.emplace(source.key(), found->key());
    }
    return bindings;
}
NativeRecord spatial_study(const NativeRecord &producer, const NativeRecord &request) {
    const auto spatial = producer.at("result").at("result");
    require(spatial.at("schema").text() == "gadget.spatial.fields.v1", "unsupported spatial producer");
    const auto source = spatial.at("source_study");
    auto docs = source.at("documents").unbound_value();
    auto raw = spatial.at("source_documents").unbound_value();
    const auto protocol = producer.at("request_publication");
    const auto bytes = protocol.bytes();
    const auto protocol_id = "spatial_protocol";
    for (std::size_t i = 0; i < source.at("documents").size(); ++i)
        require(source.at("documents").at(i).at("id").text() != protocol_id,
                "spatial protocol document identity already exists");
    std::get<N::Array>(docs.data).elements.emplace_back(N::Object{{{"id", string(protocol_id)},
        {"path", string(protocol_id)}, {"sha256", N("sha256:" + sha256(std::string_view(
            reinterpret_cast<const char *>(bytes.data()), bytes.size())))}, {"format", string("binary")}}});
    std::get<N::Array>(raw.data).elements.emplace_back(N::Object{{{"id", string(protocol_id)},
        {"bytes", protocol.materialize_value()}}});
    auto fields = source.at("fields").unbound_value();
    std::map<std::string, std::string> spatial_keys;
    const auto additions = spatial.at("fields");
    const N evidence(N::Array{{span_reference(protocol_id, bytes.size())}});
    for (std::size_t i = 0; i < additions.size(); ++i) {
        const auto f = additions.at(i);
        const auto name = f.at("id").text();
        spatial_keys.emplace(name, f.at("key").text());
        for (std::size_t j = 0; j < source.at("fields").size(); ++j)
            require(source.at("fields").at(j).at("id").text() != "all." + name,
                    "spatial field identity already exists");
        std::get<N::Array>(fields.data).elements.emplace_back(N::Object{{{"id", N("all." + name)},
            {"key", N(f.at("key").text())}, {"context", string("all")}, {"channel", string("covariate")},
            {"quantity", N(f.at("quantity").text())}, {"units", N(f.at("units").text())},
            {"method", string("Declared spatial coordinate map")}, {"value_type", string("quantitative")},
            {"evidence", evidence}}});
    }
    N::Array bindings, metrics, contexts;
    for (std::size_t k = 0; k < producer.endpoint()->source.grade_count(); ++k) {
        N binding = k == 0
            ? N(N::Object{{{"kind", string("spans")},
                {"document_index_field", N(spatial_keys.at("spatial.source_document"))},
                {"start_field", N(spatial_keys.at("spatial.source_start"))},
                {"end_field", N(spatial_keys.at("spatial.source_end"))},
                {"documents", N(N::Array{{string("carrier_axis")}})}}})
            : N(N::Object{{{"kind", string("documents")},
                {"documents", N(N::Array{{string(protocol_id)}})}}});
        bindings.elements.emplace_back(N::Object{{{"grade", N(rfc::numeric::Integer(std::to_string(k)))},
            {"id_prefix", N("grade" + std::to_string(k) + ":")},
            {"label_prefix", N(k == 0 ? "Observation " : "Relation ")},
            {"kind", N(k == 0 ? "source_observation" : "retained_relation")}, {"evidence", std::move(binding)}}});
        metrics.elements.emplace_back(N::Object{{{"kind", string("identity")}}});
    }
    // These source assays declare identity metrics. A different existing form
    // requires a supplied transport to the enlarged carrier.
    const auto original_contexts = source.at("contexts");
    require(original_contexts.size() == 1 && original_contexts.at(0).at("id").text() == "all",
            "spatial materialization requires the complete source context");
    const auto original_metrics = original_contexts.at(0).at("metrics");
    for (std::size_t i = 0; i < original_metrics.size(); ++i)
        require(original_metrics.at(i).at("kind").text() == "identity",
                "spatial materialization needs an explicit transport for a nonidentity source metric");
    contexts.elements.emplace_back(N::Object{{{"id", string("all")}, {"block", N{}},
        {"metrics", N(metrics)}, {"time", N{}}}});
    const auto charts = spatial.at("charts");
    for (std::size_t i = 0; i < charts.size(); ++i)
        contexts.elements.emplace_back(N::Object{{{"id", N(charts.at(i).at("id").text())},
            {"block", N(charts.at(i).at("block").text())}, {"metrics", N(metrics)}, {"time", N{}},
            {"witness_scope", string("block")}}});
    auto design = request.at("design").unbound_value(retained_request_sources(producer, request));
    require(!request.at("design").has("covariates"), "derived design covariates use field bindings");
    std::get<N::Object>(design.data).members.emplace_back("covariates", N(N::Array{}));
    require(request.at("channel").text() == source.at("fields").at(0).at("channel").text(),
            "spatial materialization must preserve the measured field channel");
    N study(N::Object{{{"schema", string("gadget.scientific.study.v3")}, {"id", N(request.at("id").text())},
        {"design", std::move(design)}, {"documents", std::move(docs)}, {"cell_bindings", N(std::move(bindings))},
        {"occurrence_evidence", string("parent_cell")}, {"contexts", N(std::move(contexts))},
        {"fields", std::move(fields)}, {"field_sets", source.at("field_sets").unbound_value()}}});
    return ScientificStudy::publish_annotations(producer, std::move(study), std::move(raw));
}
NativeRecord coverage_study(const NativeRecord &producer, const NativeRecord &request) {
    const auto coverage = producer.at("result").at("result");
    require(coverage.at("schema").text() == "gadget.coverage.fields.v1",
            "unsupported derived native source");
    const auto resources = coverage.at("resources"), samples = coverage.at("request").at("samples");
    N::Array docs, raw_docs, all_evidence, names;
    std::map<std::string, std::pair<std::string, std::size_t>> documents;
    for (std::size_t i = 0; i < resources.size(); ++i) {
        const auto resource = resources.at(i);
        const auto id = resource.at("id").text();
        const auto name = "resource:" + id;
        const auto bytes = resource.at("bytes").bytes();
        const auto hash =
            "sha256:" +
            sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
        docs.elements.emplace_back(N::Object{
            {{"id", N(name)}, {"path", N(name)}, {"sha256", N(hash)}, {"format", string("utf8")}}});
        raw_docs.elements.emplace_back(
            N::Object{{{"id", N(name)}, {"bytes", resource.at("bytes").materialize_value()}}});
        documents.emplace(id, std::pair{name, bytes.size()});
        if (!bytes.empty())
            all_evidence.elements.push_back(span_reference(name, bytes.size()));
    }
    for (std::size_t i = 0; i < samples.size(); ++i)
        names.elements.emplace_back(documents.at(samples.at(i).at("resource").text()).first);
    N::Array bindings;
    bindings.elements.emplace_back(
        N::Object{{{"grade", N(rfc::numeric::Integer(0))},
                   {"id_prefix", string("segment:")},
                   {"label_prefix", string("Segment ")},
                   {"kind", string("physical_interval")},
                   {"evidence", N(N::Object{{{"kind", string("cofaces")},
                                             {"upper_grade", N(rfc::numeric::Integer(1))}}})}}});
    const auto intervals = coverage.at("intervals");
    bindings.elements.emplace_back(N::Object{
        {{"grade", N(rfc::numeric::Integer(1))},
         {"id_prefix", string("interval:")},
         {"label_prefix", string("Original interval ")},
         {"kind", string("source_interval")},
         {"evidence", N(N::Object{{{"kind", string("spans")},
                                   {"document_index_field", N(field_key(intervals.at("sample")))},
                                   {"start_field", N(field_key(intervals.at("source_start")))},
                                   {"end_field", N(field_key(intervals.at("source_end")))},
                                   {"documents", N(std::move(names))}}})}}});
    N::Array metrics, contexts, fields;
    metrics.elements.emplace_back(
        N::Object{{{"kind", string("diagonal_field")},
                   {"key", N(field_key(coverage.at("metric_diagonal")))}}});
    metrics.elements.emplace_back(N::Object{{{"kind", string("identity")}}});
    contexts.elements.emplace_back(N::Object{{{"id", string("all")},
                                              {"block", N{}},
                                              {"metrics", N(std::move(metrics))},
                                              {"time", N{}}}});
    require(request.at("channel").text() == "G" || request.at("channel").text() == "E" ||
                request.at("channel").text() == "T" ||
                request.at("channel").text() == "observation",
            "derived field channel");
    for (const auto &tensor : producer.endpoint()->tensors) {
        if (tensor.role != TemporalTensorRole::Field)
            continue;
        const auto id = tensor.interpretation;
        const bool categorical = id == "interval.sample" || id == "interval.contig" ||
                                 id == "interval.strand" || id == "segment.contig" ||
                                 id == "segment.strand";
        const bool bytes = id == "interval.source_start" || id == "interval.source_end" ||
                           id == "interval.id_offset" || id == "interval.id_length";
        const bool coordinate = id == "interval.start" || id == "interval.end" ||
                                id == "segment.start" || id == "segment.end" ||
                                id == "reference_length";
        const bool structural =
            categorical || bytes || coordinate || id == "interval.source_record_index";
        const auto units = bytes        ? "byte"
                           : coordinate ? "reference_base"
                           : structural ? "dimensionless"
                                        : coverage.at("request").at("units").text();
        const auto quantity =
            id.starts_with("sample:") || id.starts_with("contrast:") ? "coverage"
            : id.starts_with("interval-sample:") || id.starts_with("interval-contrast:")
                ? "interval_weight"
                : id;
        fields.elements.emplace_back(
            N::Object{{{"id", N("all." + id)},
                       {"key", N(std::to_string(tensor.key))},
                       {"context", string("all")},
                       {"channel", N(structural ? "covariate" : request.at("channel").text())},
                       {"quantity", N(quantity)},
                       {"units", N(units)},
                       {"method", string("Retained controlled interval action")},
                       {"value_type", string(categorical ? "categorical" : "quantitative")},
                       {"evidence", N(all_evidence)}}});
    }
    auto design = request.at("design").unbound_value(retained_request_sources(producer, request));
    auto &design_members = std::get<N::Object>(design.data).members;
    require(!request.at("design").has("covariates"),
            "materialization design declares covariates through the study field bindings");
    design_members.emplace_back("covariates", N(N::Array{}));
    N study(N::Object{{{"schema", string("gadget.scientific.study.v3")},
                       {"id", N(request.at("id").text())},
                       {"design", std::move(design)},
                       {"documents", N(std::move(docs))},
                       {"cell_bindings", N(std::move(bindings))},
                       {"occurrence_evidence", string("parent_cell")},
                       {"contexts", N(std::move(contexts))},
                       {"fields", N(std::move(fields))},
                       {"field_sets", N(N::Array{})}}});
    return ScientificStudy::publish_annotations(producer, std::move(study), N(std::move(raw_docs)));
}
NativeRecord quantitative_study(const NativeRecord &producer, const NativeRecord &request) {
    const auto source = producer.at("result").at("result");
    const bool expression = source.at("schema").text() == "gadget.expression.fields.v1";
    const bool aggregation = source.at("schema").text() == "gadget.aggregation.fields.v1";
    require(expression || aggregation || source.at("schema").text() == "gadget.feature.count.fields.v1",
            "unsupported native quantitative producer");
    if (source.has("coupling"))
        require(source.at("coupling").at("status").text() == "feasible",
                "infeasible coupled expression has no derived measured bounds");
    N::Array docs, raw_docs, evidence, document_names;
    const auto document = [&](std::string id, const NativeRecord &original) {
        const auto bytes = original.bytes();
        docs.elements.emplace_back(N::Object{
            {{"id", N(id)},
             {"path", N(id)},
             {"sha256",
              N("sha256:" + sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),
                                                    bytes.size())))},
             {"format", string("binary")}}});
        raw_docs.elements.emplace_back(
            N::Object{{{"id", N(id)}, {"bytes", original.materialize_value()}}});
        if (!bytes.empty()) {
            evidence.elements.push_back(span_reference(id, bytes.size()));
            document_names.elements.emplace_back(id);
        }
    };
    if (aggregation) {
        const auto originals = source.at("original_documents");
        for (std::size_t i = 0; i < originals.size(); ++i)
            document(originals.at(i).at("id").text(), originals.at(i).at("bytes"));
        document("source_workspace", producer.at("result").at("workspace_publication"));
    } else if (expression) {
        const auto resources = source.at("resources");
        for (std::size_t i = 0; i < resources.size(); ++i) {
            const auto resource = resources.at(i);
            const auto prefix = resource.at("role").text() == "reference" ? "reference:" : "resource:";
            document(prefix + resource.at("id").text(), resource.at("bytes"));
        }
    } else {
        document("alignments", source.at("original_alignments"));
        document("annotation", source.at("original_annotation"));
        if (source.has("original_reference"))
            document("reference", source.at("original_reference"));
    }
    document(aggregation ? "aggregation_protocol" : "assignment_protocol",
             producer.at("request_publication"));
    N::Array bindings, metrics, contexts, fields;
    for (std::size_t k = 0; k < producer.endpoint()->source.grade_count(); ++k) {
        bindings.elements.emplace_back(N::Object{
            {{"grade", N(rfc::numeric::Integer(std::to_string(k)))},
             {"id_prefix", N(k == 0 ? (aggregation ? std::string("feature:")
                                                  : std::string("candidate_or_target:"))
                                    : std::string("assignment_relation:"))},
             {"label_prefix", N(k == 0 ? std::string("Candidate or declared target ")
                                       : std::string("Assignment relation "))},
             {"kind", N(k == 0 ? (aggregation ? std::string("aggregate_feature")
                                              : std::string("assignment_coordinate"))
                               : std::string("declared_assignment_relation"))},
             {"evidence",
              N(N::Object{{{"kind", string("documents")}, {"documents", N(document_names)}}})}}});
        if (aggregation) {
            auto &binding = std::get<N::Object>(bindings.elements.back().data);
            binding.members.emplace_back("ids", source.at("feature_ids").unbound_value());
            binding.members.emplace_back("labels", source.at("feature_ids").unbound_value());
        }
        metrics.elements.emplace_back(N::Object{{{"kind", string("identity")}}});
    }
    contexts.elements.emplace_back(N::Object{{{"id", string("all")},
                                              {"block", N{}},
                                              {"metrics", N(std::move(metrics))},
                                              {"time", N{}}}});
    const auto channel = request.at("channel").text();
    require(channel == "G" || channel == "E" || channel == "T" || channel == "observation",
            "derived field channel");
    for (const auto &tensor : producer.endpoint()->tensors) {
        if (tensor.role != TemporalTensorRole::Field)
            continue;
        const auto carrier = std::find_if(producer.endpoint()->carriers.begin(),
            producer.endpoint()->carriers.end(), [&](const auto &value) {return value.key == tensor.output;});
        require(carrier != producer.endpoint()->carriers.end(), "quantitative field carrier absent");
        if (carrier->kind != TemporalCarrierKind::SourceGrade) continue;
        auto id = tensor.interpretation;
        if (!expression && !aggregation && id.starts_with("sample:0."))
            id = source.at("request").at("resource").text() + id.substr(8);
        fields.elements.emplace_back(
            N::Object{{{"id", N("all." + id)},
                       {"key", N(std::to_string(tensor.key))},
                       {"context", string("all")},
                       {"channel", N(channel)},
                       {"quantity", N(aggregation ? source.at("summary").at("quantity").text()
                                         : tensor.interpretation.ends_with("class_count")
                                          ? std::string("observed_class_count")
                                          : std::string("assignment_count_bound"))},
                       {"units", N(aggregation ? source.at("summary").at("units").text()
                                               : source.at("allocation_fields").at("units").text())},
                       {"method", string(tensor.interpretation.starts_with("coupled.")
                                                     ? "Declared coupled constraints and complete exact local search"
                                                     : aggregation ? "Declared feature and observation maps"
                                                     : "Declared class incidence action")},
                       {"value_type", string("quantitative")},
                       {"evidence", N(evidence)}}});
    }
    auto design = request.at("design").unbound_value(retained_request_sources(producer, request));
    require(!request.at("design").has("covariates"), "derived design uses study field bindings");
    std::get<N::Object>(design.data).members.emplace_back("covariates", N(N::Array{}));
    N study(N::Object{{{"schema", string("gadget.scientific.study.v3")},
                       {"id", N(request.at("id").text())},
                       {"design", std::move(design)},
                       {"documents", N(std::move(docs))},
                       {"cell_bindings", N(std::move(bindings))},
                       {"occurrence_evidence", string("parent_cell")},
                       {"contexts", N(std::move(contexts))},
                       {"fields", N(std::move(fields))},
                       {"field_sets", N(N::Array{})}}});
    return ScientificStudy::publish_annotations(producer, std::move(study), N(std::move(raw_docs)));
}
NativeRecord structured_study(const NativeRecord &producer, const NativeRecord &request) {
    const auto source = producer.at("result").at("result");
    const auto schema = source.at("schema").text();
    require(schema == "gadget.sequence.alignment.fields.v1" || schema == "gadget.sequence.assembly.fields.v1" ||
                schema == "gadget.constraints.fields.v1" || schema == "gadget.genotype.fields.v1" || schema == "gadget.chromatin.fields.v1" ||
                schema == "gadget.sections.selected.v1",
            "unsupported structured native producer");
    if (schema == "gadget.constraints.fields.v1")
        require(source.at("summary").at("status").text() == "feasible",
                "infeasible allocation has no derived measured bounds");
    N::Array docs, raw_docs, evidence, document_names, bindings, metrics, fields;
    const auto document = [&](std::string id, const NativeRecord &original) {
        const auto bytes = original.bytes();
        docs.elements.emplace_back(N::Object{
            {{"id", N(id)},
             {"path", N(id)},
             {"sha256",
              N("sha256:" + sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),
                                                    bytes.size())))},
             {"format", string("binary")}}});
        raw_docs.elements.emplace_back(
            N::Object{{{"id", N(id)}, {"bytes", original.materialize_value()}}});
        if (!bytes.empty()) {
            evidence.elements.push_back(span_reference(id, bytes.size()));
            document_names.elements.emplace_back(id);
        }
    };
    document("source_workspace", producer.at("result").at("workspace_publication"));
    document("analysis_protocol", producer.at("request_publication"));
    if (source.has("sections_publication"))
        document("interface_sections", source.at("sections_publication"));
    for (std::size_t k = 0; k < producer.endpoint()->source.grade_count(); ++k) {
        N::Object binding{{{"grade", N(rfc::numeric::Integer(std::to_string(k)))},
                           {"id_prefix", N("grade" + std::to_string(k) + ":")},
                           {"label_prefix", N("Grade " + std::to_string(k) + " cell ")},
                           {"kind", string("retained_source_cell")},
                           {"evidence", N(N::Object{{{"kind", string("documents")},
                                                     {"documents", N(document_names)}}})}}};
        if (source.has("cell_ids") &&
            source.at("cell_ids").at(k).kind() != NativeRecordKind::Null) {
            binding.members.emplace_back("ids", source.at("cell_ids").at(k).unbound_value());
            binding.members.emplace_back("labels", source.at("cell_ids").at(k).unbound_value());
        }
        bindings.elements.emplace_back(std::move(binding));
        metrics.elements.push_back(source.has("metrics")
                                       ? source.at("metrics").at(k).unbound_value()
                                       : N(N::Object{{{"kind", string("identity")}}}));
    }
    const auto channel = request.at("channel").text();
    require(channel == "G" || channel == "E" || channel == "T" || channel == "observation",
            "derived field channel");
    std::set<std::string> ids;
    const auto declarations = source.at("fields");
    for (std::size_t f = 0; f < declarations.size(); ++f) {
        const auto declaration = declarations.at(f);
        const auto id = declaration.at("id").text();
        require(ids.insert(id).second, "derived field identity is ambiguous");
        fields.elements.emplace_back(
            N::Object{{{"id", N("all." + id)},
                       {"key", N(declaration.at("key").integer().get_str())},
                       {"context", string("all")},
                       {"channel", N(channel)},
                       {"quantity", N(declaration.at("quantity").text())},
                       {"units", N(declaration.at("units").text())},
                       {"method", N(declaration.at("method").text())},
                       {"value_type", string("quantitative")},
                       {"evidence", N(evidence)}}});
    }
    auto design = request.at("design").unbound_value(retained_request_sources(producer, request));
    require(!request.at("design").has("covariates"), "derived design uses study field bindings");
    std::get<N::Object>(design.data).members.emplace_back("covariates", N(N::Array{}));
    N study(N::Object{{{"schema", string("gadget.scientific.study.v3")},
                       {"id", N(request.at("id").text())},
                       {"design", std::move(design)},
                       {"documents", N(std::move(docs))},
                       {"cell_bindings", N(std::move(bindings))},
                       {"occurrence_evidence", string("parent_cell")},
                       {"contexts", N(N::Array{{N(N::Object{{{"id", string("all")},
                                                             {"block", N{}},
                                                             {"metrics", N(std::move(metrics))},
                                                             {"time", N{}}}})}})},
                       {"fields", N(std::move(fields))},
                       {"field_sets", N(N::Array{})}}});
    return ScientificStudy::publish_annotations(producer, std::move(study), N(std::move(raw_docs)));
}
void equal_source(const NativeRecord &expected, const NativeRecord &actual);
void equal_record_content(const NativeRecord &expected, const NativeRecord &actual,
                          bool software_receipt = false) {
    require(expected.kind() == actual.kind(), "replayed declaration constructor differs");
    switch (expected.kind()) {
    case NativeRecordKind::Null: return;
    case NativeRecordKind::Boolean:
        require(expected.boolean() == actual.boolean(), "replayed boolean differs"); return;
    case NativeRecordKind::Integer:
        require(expected.integer() == actual.integer(), "replayed integer differs"); return;
    case NativeRecordKind::Rational: case NativeRecordKind::ExactNumber:
        require(expected.rational() == actual.rational() &&
                    (expected.kind() != NativeRecordKind::ExactNumber ||
                     expected.approximate_view() == actual.approximate_view()),
                "replayed exact value differs"); return;
    case NativeRecordKind::Decimal:
        require(expected.rational() == actual.rational() && expected.text() == actual.text(),
                "replayed decimal differs"); return;
    case NativeRecordKind::Text: case NativeRecordKind::Bytes: {
        const auto a = expected.bytes(), b = actual.bytes();
        require(a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin()),
                "replayed original bytes differ"); return;
    }
    case NativeRecordKind::Object:
        require(expected.names() == actual.names(), "replayed declaration members differ");
        for (const auto &name : expected.names()) {
            if (software_receipt && (name == "gadget_sha256" || name == "rexgraph_sha256"))
                continue;
            const bool identity = name == "allocation_identity" && expected.has("schema") &&
                expected.at("schema").text() == "gadget.feature.count.fields.v1";
            if (name == "sections_publication" || name == "constraint_publication") {
                const auto old = NativeRecord::decode(expected.at(name).bytes());
                const auto now = NativeRecord::decode(actual.at(name).bytes());
                require(old.source_namespace() == now.source_namespace(),
                        "replayed retained publication namespace differs");
                equal_source(old, now);
            } else
                equal_record_content(expected.at(name), actual.at(name), identity);
        }
        return;
    case NativeRecordKind::Array:
        require(expected.size() == actual.size(), "replayed declaration cardinality differs");
        for (std::size_t i = 0; i < expected.size(); ++i) equal_record_content(expected.at(i), actual.at(i));
        return;
    case NativeRecordKind::Tensor:
        require(expected.tensor().domain() == actual.tensor().domain() &&
                    temporal_endpoint_detail::support_equal(expected.tensor().support(), actual.tensor().support()) &&
                    expected.tensor().values().entries() == actual.tensor().values().entries(),
                "replayed tensor value differs"); return;
    }
    throw std::invalid_argument("unsupported replay declaration constructor");
}
void equal_source(const NativeRecord &expected, const NativeRecord &actual) {
    const auto a = expected.endpoint(), b = actual.endpoint();
    require(encode_rfcx(a->source).bytes == encode_rfcx(b->source).bytes,
            "replayed complete relational source differs");
    require(a->local_blocks.size() == b->local_blocks.size(), "replayed local chart family differs");
    for (std::size_t i = 0; i < a->local_blocks.size(); ++i) {
        const auto &x = a->local_blocks[i], &y = b->local_blocks[i];
        require(x.block_key == y.block_key && x.kind == y.kind && x.cells == y.cells &&
                    x.occurrences == y.occurrences && x.payloads.empty() && y.payloads.empty(),
                "replayed local chart differs");
    }
    require(a->carriers.size() == b->carriers.size() && a->tensors.size() == b->tensors.size(),
            "replayed carrier/field family differs");
    for (std::size_t i = 0; i < a->carriers.size(); ++i) {
        const auto &x = a->carriers[i], &y = b->carriers[i];
        require(x.key == y.key && x.kind == y.kind && x.grade == y.grade && x.cell == y.cell &&
                    x.domain == y.domain && x.coordinate_system == y.coordinate_system,
                "replayed carrier differs");
    }
    for (std::size_t i = 0; i < a->tensors.size(); ++i) {
        const auto &x = a->tensors[i], &y = b->tensors[i];
        require(x.key == y.key && x.role == y.role && x.input == y.input && x.output == y.output &&
                    x.interpretation == y.interpretation &&
                    temporal_endpoint_detail::tensor_equal(x, y),
                "replayed field/support differs for tensor " + std::to_string(x.key) + " (" +
                    x.interpretation + ")");
    }
    const auto old = expected.has("result") ? expected.at("result").at("result") : expected,
               now = actual.has("result") ? actual.at("result").at("result") : actual;
    require(old.at("schema").text() == now.at("schema").text(),
            "replayed declaration schema differs");
    if (old.at("schema").text() != "gadget.coverage.fields.v1") {
        equal_record_content(old, now);
        return;
    }
    const auto originals = old.at("resources"), replayed = now.at("resources");
    require(originals.size() == replayed.size(), "replayed original resource count differs");
    for (std::size_t i = 0; i < originals.size(); ++i) {
        require(originals.at(i).at("id").text() == replayed.at(i).at("id").text(),
                "replayed resource identity differs");
        const auto x = originals.at(i).at("bytes").bytes(), y = replayed.at(i).at("bytes").bytes();
        require(x.size() == y.size() && std::equal(x.begin(), x.end(), y.begin()),
                "replayed original bytes differ");
    }
}
} // namespace
NativeRecord AssayWorkspace::publish_studies(
    const NativeRecord &request, const std::string &id,
    const std::map<std::string, std::shared_ptr<const ScientificStudy>> &studies) {
    require(!id.empty() && !studies.empty(), "native study collection needs an id and studies");
    N::Array assays, publications;
    std::vector<RawByteSource> raw;
    std::map<const ScientificStudy *, std::pair<std::uint64_t, std::size_t>> retained;
    const auto key = std::uint64_t{1};
    auto request_bytes = request.owning_publication().archive_source(key);
    const auto request_size = request_bytes.bytes().size();
    raw.push_back(std::move(request_bytes));
    std::map<std::uint64_t, std::uint64_t> request_sources;
    for (const auto &source : request.raw_sources()) {
        const auto retained_key = raw.size() + 1;
        raw.push_back(source.rekey(retained_key));
        request_sources.emplace(source.key(), retained_key);
    }
    auto design = request.at("design").unbound_value(request_sources);
    for (const auto &[name, study] : studies) {
        require(study && !name.empty(), "native collection study is absent");
        if (!retained.contains(study.get())) {
            const auto source = raw.size() + 1;
            auto bytes = study->native_publication().archive_source(source);
            const auto size = bytes.bytes().size();
            raw.push_back(std::move(bytes));
            retained.emplace(study.get(), std::pair{source, size});
        }
        const auto [source, size] = retained.at(study.get());
        const auto path = "studies/" + name + ".rfcxs";
        N::Array packets;
        packets.elements.emplace_back(N::Object{
            {{"path", N(path)},
             {"source_namespace", N(std::to_string(study->native_history()->source_namespace()))},
             {"source_digest", N(study->digest())}}});
        assays.elements.emplace_back(N::Object{
            {{"declaration",
              N(N::Object{{{"schema", string("gadget.assay.native.v1")},
                           {"id", N(name)},
                           {"kind", string("native_study")},
                           {"source_namespace",
                            N(std::to_string(study->native_history()->source_namespace()))},
                           {"carrier", string("source_grades")}}})},
             {"originals", N(N::Object{})},
             {"partitions", N(std::move(packets))}}});
        publications.elements.emplace_back(
            N::Object{{{"assay", N(name)},
                       {"partition", N(rfc::numeric::Integer(0))},
                       {"path", N(path)},
                       {"bytes", N(N::ByteSpan{source, 0, size, false})}}});
    }
    N manifest(N::Object{{{"schema", string("gadget.assays.workspace.v2")},
                          {"id", N(id)},
                          {"design", std::move(design)},
                          {"assays", N(std::move(assays))},
                          {"resources", N(N::Array{})},
                          {"links", N(N::Array{})},
                          {"gadget_sha256", string(gadget::build::source_sha256())},
                          {"rexgraph_sha256", string(gadget::build::rexgraph_sha256())}}});
    return NativeRecord::create(
        N(N::Object{{{"schema", string("gadget.assay.workspace.native.v2")},
                     {"manifest", std::move(manifest)},
                     {"originals", N(N::Array{})},
                     {"publications", N(std::move(publications))},
                     {"request_publication", N(N::ByteSpan{key, 0, request_size, false})}}}),
        new_native_record_namespace(), std::move(raw));
}
Json AssayWorkspace::materialize_record(const NativeRecord &request, const fs::path &destination,
                                        std::shared_ptr<CpuResources> resources) const {
    const auto workspace = materialize_workspace(request, destination, std::move(resources));
    auto result = render_record_view(workspace->native_publication().at("manifest"));
    result.set("manifest_sha256", Json::string(workspace->digest()));
    result.set("capabilities", bioinformatics_capabilities());
    return result;
}
bool AssayWorkspace::supports_native_materialization(std::string_view operation) {
    return operation == "controlled_coverage" || operation == "spatial_relations" ||
           operation == "feature_counts" || operation == "expression_allocations" ||
           operation == "aggregate_assay" || operation == "design_aggregation" ||
           operation == "variant_genotypes" || operation == "chromatin_contrast" ||
           operation == "assay_sections" || operation == "constrained_allocations" ||
           operation == "exact_sequence_alignment" || operation == "sequence_assembly_graph";
}
bool AssayWorkspace::supports_native_materialization(const NativeRecord &recipe) {
    const auto operation = recipe.at("operation").text();
    return supports_native_materialization(operation);
}
std::shared_ptr<const AssayWorkspace>
AssayWorkspace::materialize_workspace(const NativeRecord &request, const fs::path &destination,
                                      std::shared_ptr<CpuResources> resources) const {
    require(supports_native_materialization(request.at("recipe")),
            "operation has no native materialization producer");
    return materialize_native_fields(request, destination, std::move(resources));
}
std::shared_ptr<const AssayWorkspace>
AssayWorkspace::materialize_native_fields(const NativeRecord &request, const fs::path &destination,
                                          std::shared_ptr<CpuResources> resources,
                                          std::optional<std::uint64_t> section_namespace,
                                          const NativeRecord *previous_coupling) const {
    const auto view = render_record_view(request);
    members(view, {"schema", "id", "source_namespace", "design", "recipe", "channel"},
            {"view", "component"});
    const bool section = request.at("recipe").at("operation").text() == "assay_sections";
    require(section ? request.has("view") && request.has("component")
                    : !request.has("view") && !request.has("component"),
            "section selection scope");
    require(request.at("schema").text() == "gadget.assays.materialize.v1",
            "materialization schema");
    require(!fs::exists(destination), "derived assay destination exists");
    const auto ns = key(view.at("source_namespace"));
    require(ns != native_publication().source_namespace() && ns != request.source_namespace(),
            "derived source must have its own namespace");
    const bool spatial = request.at("recipe").at("operation").text() == "spatial_relations";
    const bool coverage = request.at("recipe").at("operation").text() == "controlled_coverage";
    const auto name = request.at("recipe").at("operation").text();
    const bool aggregation = name == "aggregate_assay" || name == "design_aggregation";
    const bool structured = section || name == "variant_genotypes" || name == "chromatin_contrast" || name == "constrained_allocations" || name == "exact_sequence_alignment" || name == "sequence_assembly_graph";
    const auto operation =
        structured    ? execute_structured(request.at("recipe"), std::move(resources), ns, true,
                                           section ? &request : nullptr, section_namespace)
        : spatial     ? execute_spatial(request.at("recipe"), std::move(resources), ns)
        : coverage    ? execute_coverage(request.at("recipe"), std::move(resources), ns)
        : aggregation ? execute_aggregation(request.at("recipe"), std::move(resources), ns)
                      : execute_assignment(request.at("recipe"), std::move(resources), ns, previous_coupling);
    auto admitted = std::make_shared<const ScientificStudy>(
        structured ? structured_study(operation, request)
        : spatial  ? spatial_study(operation, request)
        : coverage ? coverage_study(operation, request)
                   : quantitative_study(operation, request));
    const auto collection =
        publish_studies(request, request.at("id").text(), {{request.at("id").text(), admitted}});
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage =
        private_directory(fs::absolute(destination).parent_path(), ".gadget-native-derive-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    auto workspace = std::shared_ptr<const AssayWorkspace>(
        new AssayWorkspace(destination / "workspace.rfcxs",
                           {{request.at("id").text(), {std::move(admitted)}}}, collection));
    collection.write(stage / "workspace.rfcxs");
    fs::rename(stage, destination);
    return workspace;
}
Json AssayWorkspace::replay_native_derivation(const std::string &id, const fs::path &destination,
                                              std::shared_ptr<CpuResources> resources) const {
    const auto original = study(id);
    const auto producer = original->producer_publication();
    require(producer.at("schema").text() == "gadget.operation.native.v1",
            "missing generating operation");
    const auto workspace_bytes = producer.at("result").at("workspace_publication").bytes();
    const auto request_bytes = producer.at("request_publication").bytes();
    const auto request = NativeRecord::decode(request_bytes);
    require(request.at("schema").text() == "gadget.assays.materialize.v1" &&
                supports_native_materialization(request.at("recipe")),
            "generating request is not a native field materialization");
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage =
        private_directory(fs::absolute(destination).parent_path(), ".gadget-native-replay-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    write_file_exact(stage / "original-workspace.rfcxs", workspace_bytes);
    AssayWorkspace generating(stage / "original-workspace.rfcxs");
    require(generating.digest() == producer.at("result").at("source_digest").text(),
            "generating workspace binding differs");
    std::optional<std::uint64_t> section_namespace;
    const auto generated_fields = producer.at("result").at("result");
    if (generated_fields.has("sections_publication"))
        section_namespace =
            NativeRecord::decode(generated_fields.at("sections_publication").bytes())
                .source_namespace();
    const auto previous_coupling = generated_fields.has("coupling")
        ? std::optional<NativeRecord>(generated_fields.at("coupling")) : std::nullopt;
    const auto regenerated_workspace = generating.materialize_native_fields(
        request, stage / "replayed", std::move(resources), section_namespace,
        previous_coupling ? &*previous_coupling : nullptr);
    const auto regenerated = regenerated_workspace->study(request.at("id").text());
    const auto actual = regenerated->producer_publication();
    equal_source(producer, actual);
    auto receipt = Json::object();
    receipt.set("schema", Json::string("gadget.assay.derivation.replay.v2"));
    receipt.set("status", Json::string("PASS"));
    receipt.set("source_digest", Json::string(original->digest()));
    receipt.set("replayed_source_digest",
                Json::string(regenerated->digest()));
    receipt.set("source_cells", original->source_index().at("catalog_sections").at("cells"));
    receipt.set("source_occurrences",
                original->source_index().at("catalog_sections").at("occurrences"));
    const auto publication =
        retain_operation(native_publication(), request, import_record_view(receipt));
    publication.write(stage / "replayed/replay.rfcxs");
    fs::rename(stage / "replayed", destination);
    return receipt;
}
Json AssayWorkspace::derivation_certificate(const std::string &id) const {
    const auto source = study(id);
    const auto producer = source->producer_publication();
    const auto generated = producer.at("result"), coverage = generated.at("result");
    const bool spatial = coverage.at("schema").text() == "gadget.spatial.fields.v1";
    const bool assignment = coverage.at("schema").text() == "gadget.feature.count.fields.v1" ||
                            coverage.at("schema").text() == "gadget.expression.fields.v1" ||
                            coverage.at("schema").text() == "gadget.aggregation.fields.v1";
    const auto schema = coverage.at("schema").text();
    const bool structured = schema == "gadget.sequence.alignment.fields.v1" || schema == "gadget.sequence.assembly.fields.v1" ||
                            schema == "gadget.constraints.fields.v1" || schema == "gadget.genotype.fields.v1" ||
                            schema == "gadget.chromatin.fields.v1" ||
                            schema == "gadget.sections.selected.v1";
    require(spatial || assignment || structured || schema == "gadget.coverage.fields.v1",
            "unsupported native generating result");
    auto result = (assignment || structured) ? render_record_view(coverage) : Json::object();
    if (!assignment && !structured) {
        result.set("schema", Json::string(spatial ? "gadget.spatial.fields.reference.v1"
                                                  : "gadget.coverage.fields.reference.v1"));
        result.set("native", source->declaration().at("native"));
        result.set("source_index", source->source_index());
        for (const auto *name : {"intervals", "segments", "sample_fields", "metric_diagonal",
                                 "contrasts", "units", "thresholds", "fields", "frame", "charts",
                                 "neighbor_begin", "neighbor_count"})
            if (coverage.has(name))
                result.set(name, render_record_view(coverage.at(name)));
        if (!spatial)
            result.set("control_rank", render_record_view(coverage.at("controls").at("rank")));
        result.set("native_declaration_sha256", Json::string(native_declaration_digest(generated)));
    }
    auto certificate = Json::object();
    certificate.set("schema", Json::string("gadget.assay.result.v1"));
    certificate.set("query", render_record_view(generated.at("query")));
    certificate.set("source_digest", Json::string(generated.at("source_digest").text()));
    certificate.set("result", result);
    return seal(certificate, "certificate_sha256");
}
} // namespace gadget
