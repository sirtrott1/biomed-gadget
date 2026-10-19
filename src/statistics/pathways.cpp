#include "pathway_ora.hpp"
#include "../archive_identity.hpp"
#include <algorithm>
#include <gadget/statistics/pathways.hpp>
#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <limits>
#include <set>

namespace gadget::statistics {
namespace {
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("native pathway operation cancelled");
}
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::invalid_argument("native pathway admission: " + message);
}
void members(const NativeRecord &record, const std::set<std::string> &names) {
    require(record.kind() == NativeRecordKind::Object, "expected native object");
    for (const auto &name : names)
        require(record.has(name), "missing source member: " + name);
    for (const auto &name : record.names())
        require(names.contains(name), "unknown source member: " + name);
}
std::string text(const NativeRecord &record) {
    require(record.kind() == NativeRecordKind::Text && record.bytes().size() <= 4096 &&
                !record.bytes().empty(),
            "native source text representation");
    auto out = record.text();
    require(out.find('\0') == std::string::npos, "NUL in source declaration");
    return out;
}
std::string_view original(const NativeRecord &record, std::size_t limit) {
    require(record.kind() == NativeRecordKind::Text || record.kind() == NativeRecordKind::Bytes,
            "actual immutable source bytes required");
    const auto bytes = record.bytes();
    require(bytes.size() <= limit, "native source byte budget");
    return {reinterpret_cast<const char *>(bytes.data()), bytes.size()};
}
} // namespace
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources) {
    return statistics::pathway_overrepresentation(analysis, gene_sets, request, std::move(resources), {});
}
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources,
                                        std::stop_token stop) {
    cancelled(stop);
    require(analysis.is_root() && gene_sets.is_root(),
            "input archive pins require owning root records");
    const auto query = render_record_view(request); // Request only; never render molecular results.
    const auto p = pathway_detail::admit_request(query);
    if (!resources)
        resources = CpuResources::local();
    (void)resources->readmit();
    auto reservation = resources->coordinator().reserve(1, rfc::machine::checked_size_add(
        p.declared_buffer_bytes,rfc::machine::checked_size_mul(p.archive_bytes,2)));
    require(archive_identity(analysis, p.archive_bytes, stop) == p.analysis_sha256,
            "molecular archive SHA mismatch");
    cancelled(stop);
    require(archive_identity(gene_sets, p.archive_bytes, stop) == p.gene_sets_sha256,
            "released gene set archive SHA mismatch");
    cancelled(stop);
    members(gene_sets, {"schema", "id", "release", "identifier_namespace", "format", "sha256",
                        "data", "mapping"});
    require(text(gene_sets.at("schema")) == "gadget.pathways.source.v1" &&
                text(gene_sets.at("format")) == "gmt",
            "released source schema/format");
    const auto mapping = gene_sets.at("mapping");
    members(mapping, {"id", "release", "origin", "feature_namespace", "identifier_namespace",
                      "format", "sha256", "data"});
    require(text(mapping.at("format")) == "tsv", "correspondence source format");
    pathway_detail::Source source;
    source.id = text(gene_sets.at("id"));
    source.release = text(gene_sets.at("release"));
    source.identifier_namespace = text(gene_sets.at("identifier_namespace"));
    source.sha256 = text(gene_sets.at("sha256"));
    require(text(mapping.at("identifier_namespace")) == source.identifier_namespace,
            "mapping target namespace differs from released gene sets");
    source.mapping_id = text(mapping.at("id"));
    source.mapping_release = text(mapping.at("release"));
    source.mapping_origin = text(mapping.at("origin"));
    source.feature_namespace = text(mapping.at("feature_namespace"));
    source.mapping_sha256 = text(mapping.at("sha256"));
    source.gmt = original(gene_sets.at("data"), p.source_bytes);
    source.mapping = original(mapping.at("data"), p.source_bytes - source.gmt.size());
    const auto features = select_scientific_result(analysis, request.at("features_path"));
    require(features.kind() == NativeRecordKind::Array && features.size() <= p.features,
            "actual native feature array budget");
    require(features.size() <= p.declared_buffer_bytes / 512, "native feature working byte admission");
    std::vector<std::string> identities;
    identities.reserve(features.size());
    std::size_t identity_bytes{};
    for (std::size_t i = 0; i < features.size(); ++i) {
        cancelled(stop);
        const auto identity = p.feature_id_path.size() ? select_scientific_result(features.at(i), request.at("feature_id_path")) : features.at(i);
        require(identity.kind() == NativeRecordKind::Text &&
                    identity.bytes().size() <= p.field_bytes,
                "native feature identity representation");
        require(identity.bytes().size() <= p.declared_buffer_bytes / 8 - identity_bytes,
                "native feature identity working byte admission");
        identity_bytes += identity.bytes().size();
        identities.push_back(identity.text());
    }
    scientific_detail::ResultFields fields(analysis);
    auto result = pathway_detail::overrepresentation(identities, source, query, fields, stop);
    cancelled(stop);
    const auto analysis_key=fields.next_source_key();
    auto analysis_archive=analysis.archive_source(analysis_key);
    require(analysis_archive.bytes().size()<=p.archive_bytes,"retained analysis archive byte admission");
    result.set("analysis_publication",scientific_detail::ResultValue(NativeRecordValue(
        NativeRecordValue::ByteSpan{analysis_key,0,analysis_archive.bytes().size(),false})));
    fields.retain_sources({std::move(analysis_archive)});
    const auto source_key=fields.next_source_key();
    auto archive=gene_sets.archive_source(source_key);
    require(archive.bytes().size()<=p.archive_bytes,"retained released source archive byte admission");
    result.set("gene_set_source_archive",scientific_detail::ResultValue(NativeRecordValue(
        NativeRecordValue::ByteSpan{source_key,0,archive.bytes().size(),false})));
    fields.retain_sources({std::move(archive)});
    cancelled(stop);
    auto publication=fields.finish_operation(std::move(result),{},"",{},request);
    cancelled(stop);
    return publication;
}
} // namespace gadget::statistics
