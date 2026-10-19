#include "result_answer_witness.hpp"
#include "scientific_result_fields.hpp"
#include "pipeline_answers.hpp"
#include <gadget/scientific.hpp>
#include <gadget/biomedical.hpp>
#include <algorithm>
#include <set>

namespace gadget::pipeline_detail {
namespace {
void require(bool ok, const char *why) { if (!ok) throw std::invalid_argument(why); }
void keys(const Json &value, const std::set<std::string> &names) {
    require(json_object_is_type(value.raw(), json_type_object), "native result witness object required");
    require(static_cast<std::size_t>(json_object_object_length(value.raw())) == names.size(), "native result witness members differ");
    for (const auto &name : names) require(value.has(name), "native result witness member absent");
}
std::string digest(std::span<const std::byte> bytes) {
    return "sha256:" + sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
}
std::string archive_digest(const NativeRecord &record) {
    require(record.is_root(), "native result archive requires actual owning root");
    const auto history = record.publication_history();
    require(history.has_value(), "native result owning history absent");
    return "sha256:" + sha256_chunks([&](const auto &sink) {
        rfc::storage::stream_native_series(*history, sink);
    });
}
NativeRecord checked(const fs::path &directory, const Json &reference) {
    const auto path = checked_answer_path(directory, reference);
    const auto record = NativeRecord::read(path);
    require(archive_digest(record) == reference.at("sha256").text(), "native result answer archive hash mismatch");
    return record;
}
void same(const NativeRecord &a, const NativeRecord &b) {
    require(a.kind() == b.kind(), "native result wrapper changed producer result kind");
    switch (a.kind()) {
    case NativeRecordKind::Null: return;
    case NativeRecordKind::Boolean: require(a.boolean() == b.boolean(), "native result result boolean differs"); return;
    case NativeRecordKind::Integer: require(a.integer() == b.integer(), "native result result integer differs"); return;
    case NativeRecordKind::Rational: require(a.rational() == b.rational(), "native result result rational differs"); return;
    case NativeRecordKind::Decimal: require(a.rational() == b.rational() && a.text() == b.text(), "native result result decimal differs"); return;
    case NativeRecordKind::ExactNumber: require(a.rational() == b.rational() && a.approximate_view() == b.approximate_view(), "native result exact result differs"); return;
    case NativeRecordKind::Text: case NativeRecordKind::Bytes:
        require(std::ranges::equal(a.bytes(), b.bytes()), "native result result bytes differ"); return;
    case NativeRecordKind::Object:
        require(a.names() == b.names(), "native result result members differ");
        for (const auto &name : a.names()) same(a.at(name), b.at(name));
        return;
    case NativeRecordKind::Array:
        require(a.size() == b.size(), "native result result cardinality differs");
        for (std::size_t i = 0; i < a.size(); ++i) same(a.at(i), b.at(i));
        return;
    case NativeRecordKind::Tensor:
        throw std::invalid_argument("native result results require retained native tensor references");
    }
}
NativeRecord previous(const NativeRecord &record) {
    const auto history = record.publication_history();
    require(record.is_root() && history && history->segment_count() && history->last_sequence() > history->first_sequence(), "native result answer previous history absent");
    const auto parent = record.parent();
    require(parent.has_value(), "native result answer parent declaration absent");
    const auto owner = record.historical_publication(history->last_sequence() - 1);
    same(*parent, owner);
    return owner;
}
void control(const NativeRecord &actual, const Json &expected) {
    const auto rendered = render_record_view(actual);
    require(json_object_equal(rendered.raw(), expected.raw()), "native result actual request differs from effective request");
}
Json seal(Json value) { value.set("certificate_sha256", Json::string("sha256:" + sha256(value.dump()))); return value; }
} // namespace
Json result_answer_view(const fs::path &directory, const Json &reference) {
    const auto selected = reference.at("native_selection");
    keys(selected, {"kind", "method", "request", "analysis_sha256"});
    require(selected.at("kind").text() == "result_operation", "unsupported native result selection");
    keys(selected.at("request"), {"path", "sha256"});
    const auto method = selected.at("method").text();
    require(method == "pathway_overrepresentation" || method == "feature_evidence_context", "unreviewed native result answer method");
    const auto operation = checked(directory, reference), effective = checked(directory, selected.at("request"));
    require(operation.at("schema").text() == "gadget.operation.native.v1", "native result answer requires actual native operation");
    require(digest(operation.at("request_publication").bytes()) == selected.at("request").at("sha256").text(), "native result answer substituted request owner");
    const auto native_result = operation.at("result");
    const auto result = native_result.at("schema").text() == "gadget.scientific.fields.v1"
        ? native_result.at("view") : native_result;
    auto pins = Json::object();
    pins.set("archive_sha256", reference.at("sha256"));
    pins.set("request_archive_sha256", selected.at("request").at("sha256"));
    pins.set("request_declaration_sha256", Json::string(native_declaration_digest(operation.at("request"))));
    pins.set("result_declaration_sha256", Json::string(native_declaration_digest(result)));
    pins.set("source_namespace", Json::string(std::to_string(operation.source_namespace())));
    require(operation.publication_history().has_value(), "native result answer history absent");
    pins.set("cutoff", Json::string(std::to_string(operation.publication_history()->last_sequence())));
    std::string source_digest;
    {
        control(operation.at("request"), render_record_view(effective));
        require(effective.at("operation").text() == method, "native result dispatch method differs");
        const auto producer = previous(operation);
        require(producer.at("schema").text() == "gadget.operation.native.v1", "native result answer lost producer history");
        const auto analysis = result.has("analysis_publication")
            ? NativeRecord::decode(result.at("analysis_publication").bytes()) : previous(producer);
        source_digest = archive_digest(analysis);
        require(source_digest == selected.at("analysis_sha256").text() && result.at("analysis_sha256").text() == source_digest, "native result answer analysis owner differs");
        auto query = render_record_view(effective.at("query"));
        const auto mode = effective.at("analysis_binding").text();
        require(mode == "exact" || mode == "actual_result", "native result analysis binding mode differs");
        if (mode == "actual_result") {
            require(!query.has("analysis_sha256"), "native result template already pins analysis");
            query.set("analysis_sha256", Json::string(source_digest));
        } else require(query.at("analysis_sha256").text() == source_digest, "native result exact analysis pin differs");
        control(producer.at("request"), query);
        const auto biological_request = NativeRecord::decode(producer.at("request_publication").bytes());
        if (mode == "exact") {
            require(archive_digest(biological_request) == selected.at("request").at("sha256").text(),
                    "native result exact request owner differs");
            same(producer.at("request"), biological_request.at("query"));
        } else {
            control(biological_request, query);
            const auto history = biological_request.publication_history();
            require(biological_request.is_root() && history &&
                        biological_request.source_namespace() == effective.source_namespace(),
                    "native result bound query owner differs");
            // Child derive retains the containing publication as its parent.
            // A fresh wrapper has no temporal segment and derive retains its
            // original declaration without inventing a historical prefix.
            const auto parent_wrapper = biological_request.parent();
            require(parent_wrapper.has_value(), "native result bound query parent absent");
            same(*parent_wrapper, effective);
            require(parent_wrapper->type() == effective.type() &&
                        rfc::composite::declared_value_equal(parent_wrapper->value(), effective.value()),
                    "native result bound query original declaration differs");
            const auto &retained_registry = biological_request.endpoint()->declarations->registry;
            for (const auto &definition : effective.endpoint()->declarations->registry.definitions())
                require(retained_registry.at(definition.id) == definition,
                        "native result bound query original schema differs");
            if (history->segment_count()) {
                require(history->last_sequence() > history->first_sequence(),
                        "native result bound query history cutoff differs");
                const auto original_request = biological_request.historical_publication(history->last_sequence() - 1);
                require(archive_digest(original_request) == selected.at("request").at("sha256").text(),
                        "native result bound query original wrapper differs");
            }
        }
        pins.set("biological_request_archive_sha256", Json::string(digest(producer.at("request_publication").bytes())));
        same(native_result, producer.at("result"));
        const auto x = operation.endpoint(), y = producer.endpoint();
        require(x->carriers == y->carriers && x->tensors.size() == y->tensors.size(), "native result wrapper changed native fields");
        for (std::size_t i = 0; i < x->tensors.size(); ++i)
            require(rfc::storage::temporal_endpoint_detail::tensor_equal(x->tensors[i], y->tensors[i]), "native result wrapper changed native field");
        const auto source_name = method == "pathway_overrepresentation" ? "gene_set_source_archive" : "mapping_publication";
        const auto hash_name = method == "pathway_overrepresentation" ? "gene_sets_sha256" : "mapping_sha256";
        const auto source_hash = digest(result.at(source_name).bytes());
        require(source_hash == effective.at("source_record").at("sha256").text() && result.at(hash_name).text() == source_hash && query.at(hash_name).text() == source_hash, "native result released source owner differs");
        require(result.at("schema").text() == (method == "pathway_overrepresentation" ? "gadget.pathways.overrepresentation.result.v1" : "gadget.features.evidence-context.v1"), "native result result schema differs");
        pins.set("producer_archive_sha256", Json::string(archive_digest(producer)));
        pins.set("source_archive_sha256", Json::string(source_hash));
        if (method == "feature_evidence_context") {
            require(result.at("evidence_digest").text() == query.at("evidence_digest").text(), "native result context owner differs");
            const auto evidence = NativeRecord::decode(result.at("evidence_publication").bytes());
            require(evidence.at("receipt").at("source_digest").text() == query.at("evidence_digest").text(), "native result original evidence digest differs");
            pins.set("evidence_archive_sha256", Json::string(digest(result.at("evidence_publication").bytes())));
        }
    }
    auto certificate = Json::object();
    certificate.set("schema", Json::string("gadget.pipeline.result.certificate.v1"));
    certificate.set("source_digest", Json::string(source_digest));
    certificate.set("query", render_record_view(operation.at("request")));
    certificate.set("result", scientific_detail::result_fields_view(native_result));
    certificate.set("native_result", pins);
    certificate.set("scope", Json::string("Displays retained native results with their original request, source and history owners."));
    return seal(certificate);
}
Json result_answer_reference(const fs::path &directory, const Json &result,
    const Json &request, const std::string &method, const std::string &analysis_sha256) {
    auto reference = Json::object(), selected = Json::object();
    reference.set("path", result.at("path")); reference.set("sha256", result.at("sha256"));
    selected.set("kind", Json::string("result_operation")); selected.set("method", Json::string(method));
    selected.set("request", request); selected.set("analysis_sha256", Json::string(analysis_sha256));
    reference.set("native_selection", selected);
    const auto certificate = result_answer_view(directory, reference);
    reference.set("pointer", Json::string("/result"));
    for (const auto *name : {"source_digest", "query", "certificate_sha256"}) reference.set(name, certificate.at(name));
    return reference;
}
} // namespace gadget::pipeline_detail
