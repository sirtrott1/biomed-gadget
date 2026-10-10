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
    std::map<std::uint64_t, std::uint64_t> request_sources;
    const auto producer_sources = producer.raw_sources();
    for (const auto &source : request.raw_sources()) {
        const auto retained =
            std::find_if(producer_sources.begin(), producer_sources.end(),
                         [&](const auto &item) { return item.shares_owner(source); });
        require(retained != producer_sources.end(), "generating request byte owner was lost");
        request_sources.emplace(source.key(), retained->key());
    }
    auto design = request.at("design").unbound_value(request_sources);
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
void equal_source(const NativeRecord &expected, const NativeRecord &actual) {
    const auto a = expected.endpoint(), b = actual.endpoint();
    require(encode_rfcx(a->source).bytes == encode_rfcx(b->source).bytes,
            "replayed complete relational source differs");
    require(a->local_blocks.empty() && b->local_blocks.empty(),
            "coverage replay must retain its ambient source without reconstructed blocks");
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
                "replayed field/support differs");
    }
    const auto old = expected.at("result").at("result"), now = actual.at("result").at("result");
    require(old.at("schema").text() == now.at("schema").text(),
            "replayed declaration schema differs");
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
                          {"gadget_sha256", string(GADGET_SOURCE_SHA256)},
                          {"rexgraph_sha256", string(GADGET_REXGRAPH_SHA256)}}});
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
    if (request.at("recipe").at("operation").text() == "controlled_coverage")
        return materialize_coverage(request, destination, std::move(resources));
    return materialize(render_record_view(request), destination, std::move(resources));
}
Json AssayWorkspace::materialize_coverage(const NativeRecord &request, const fs::path &destination,
                                          std::shared_ptr<CpuResources> resources,
                                          std::shared_ptr<const ScientificStudy> *retained) const {
    const auto view = render_record_view(request);
    members(view, {"schema", "id", "source_namespace", "design", "recipe", "channel"});
    require(request.at("schema").text() == "gadget.assays.materialize.v1",
            "materialization schema");
    require(!fs::exists(destination), "derived assay destination exists");
    const auto ns = key(view.at("source_namespace"));
    require(ns != native_publication().source_namespace() && ns != request.source_namespace(),
            "derived source must have its own namespace");
    const auto operation = execute_coverage(request.at("recipe"), std::move(resources), ns);
    auto admitted = std::make_shared<const ScientificStudy>(coverage_study(operation, request));
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
    collection.write(stage / "workspace.rfcxs");
    auto result = render_record_view(collection.at("manifest"));
    result.set("manifest_sha256",
               Json::string(native_declaration_digest(collection.at("manifest"))));
    result.set("capabilities", bioinformatics_capabilities());
    fs::rename(stage, destination);
    if (retained)
        *retained = std::move(admitted);
    return result;
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
                request.at("recipe").at("operation").text() == "controlled_coverage",
            "generating request is not a native coverage materialization");
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
    std::shared_ptr<const ScientificStudy> regenerated;
    generating.materialize_coverage(request, stage / "replayed", std::move(resources),
                                   &regenerated);
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
    require(coverage.at("schema").text() == "gadget.coverage.fields.v1",
            "unsupported native generating result");
    auto result = Json::object();
    result.set("schema", Json::string("gadget.coverage.fields.reference.v1"));
    result.set("native", source->declaration().at("native"));
    result.set("source_index", source->source_index());
    for (const auto *name : {"intervals", "segments", "sample_fields", "metric_diagonal",
                             "contrasts", "units", "thresholds"})
        if (coverage.has(name))
            result.set(name, render_record_view(coverage.at(name)));
    result.set("control_rank", render_record_view(coverage.at("controls").at("rank")));
    result.set("native_declaration_sha256", Json::string(native_declaration_digest(generated)));
    auto certificate = Json::object();
    certificate.set("schema", Json::string("gadget.assay.result.v1"));
    certificate.set("query", render_record_view(generated.at("query")));
    certificate.set("source_digest", Json::string(generated.at("source_digest").text()));
    certificate.set("result", result);
    return seal(certificate, "certificate_sha256");
}
} // namespace gadget
