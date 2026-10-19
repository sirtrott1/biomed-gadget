#include <gadget/feature_context.hpp>
#include <gadget/record_view.hpp>
#include "feature_context_model.hpp"
#include "archive_identity.hpp"
#include <rfc/database/interpretation_codec.hpp>
#include <algorithm>
#include <set>

namespace gadget {
namespace {
using N=NativeRecordValue;
using namespace feature_context_detail;
void require(bool ok,const char *message) {if(!ok)throw std::invalid_argument(message);}
void bound(std::size_t n,std::size_t cap,const char *message) {if(n>cap)throw std::length_error(message);}
N number(std::size_t n) {return N(rfc::numeric::Integer(std::to_string(n)));}
std::string digest(std::span<const std::byte> bytes) {
    return "sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size()));
}
std::vector<std::string> feature_ids(const NativeRecord &analysis,const Request &q,Work &work) {
    const auto features=select_scientific_result(analysis,admit_record_view(q.feature_path));
    require(features.kind()==NativeRecordKind::Array,"analysis feature path is not an actual native array");
    bound(features.size(),q.limits.features,"analysis feature occurrence budget");
    std::vector<std::string> ids;ids.reserve(features.size());
    std::size_t identity_bytes{};
    const auto identity_path=admit_record_view(q.feature_id_path);
    for(std::size_t i=0;i<features.size();++i) {
        work.take();const auto identity=q.feature_id_path.size()
            ? select_scientific_result(features.at(i),identity_path) : features.at(i);
        require(identity.kind()==NativeRecordKind::Text,"feature identity must be native text");
        bound(identity.bytes().size(),q.limits.output_bytes/8-identity_bytes,
              "feature identity working byte admission");
        identity_bytes+=identity.bytes().size();const auto id=identity.text();
        require(!id.empty()&&id.find('\0')==std::string::npos,"analysis feature ID is empty or NUL");
        ids.push_back(id);
    }
    return ids;
}
Json mapping_declaration(const NativeRecord &source,std::size_t working_bytes) {
    const std::set<std::string> names{"schema","id","origin","release","source_url","license","format","sha256","bytes"};
    require(source.names().size()==names.size(),"undeclared mapping source member");
    auto declaration=Json::object();
    std::size_t text_bytes{};
    for(const auto &name:source.names()) {
        require(names.contains(name),"undeclared mapping source member");
        if(name!="bytes") {
            const auto value=source.at(name);
            require(value.kind()==NativeRecordKind::Text,"mapping declaration must be native text");
            bound(value.bytes().size(),working_bytes/8-text_bytes,
                  "mapping declaration working byte admission");
            text_bytes+=value.bytes().size();declaration.set(name,Json::string(value.text()));
        }
    }
    return declaration;
}
} // namespace
NativeRecord contextualize_feature_evidence(const NativeRecord &analysis,const EvidenceStore &evidence,
    const NativeRecord &mapping_source,const NativeRecord &request,std::shared_ptr<CpuResources> resources) {
    return contextualize_feature_evidence(analysis,evidence,mapping_source,request,std::move(resources),{});
}
NativeRecord contextualize_feature_evidence(const NativeRecord &analysis,const EvidenceStore &evidence,
    const NativeRecord &mapping_source,const NativeRecord &request,std::shared_ptr<CpuResources> resources,
    std::stop_token stop) {
    using namespace feature_context_detail;
    if(stop.stop_requested())throw std::runtime_error("feature evidence context cancelled");
    // Rendering is restricted to explicit control declarations. Analysis values
    // and source assertions are selected through their actual native owners.
    const auto q=admit_request(render_record_view(request));Work work(q.limits.work,stop);
    if(!resources)resources=CpuResources::local();
    (void)resources->readmit();auto reservation=resources->coordinator().reserve(1,q.limits.output_bytes);
    work.check();
    const auto owner=analysis.owning_publication(),mapping_owner=mapping_source.owning_publication();
    require(mapping_source.is_root(),"mapping must be the complete original source publication");
    require(q.evidence_digest==evidence.digest(),"feature context evidence digest mismatch");
    const auto &evidence_owner=evidence.native_publication();
    require(archive_identity(owner,q.limits.output_bytes,stop)==q.analysis_sha256,
            "feature context actual analysis archive SHA-256 mismatch");
    const auto source_roles=roles(evidence.interpreted_source().schema());
    const auto declaration=mapping_declaration(mapping_source,q.limits.output_bytes);
    const auto raw=mapping_source.at("bytes");
    require(raw.kind()==NativeRecordKind::Bytes,"mapping source must retain original opaque bytes");
    const auto original=raw.byte_source(); // Refuses slices of a different/partial byte owner.
    bound(original.bytes().size(),q.limits.source_bytes,"mapping original byte budget");
    const std::string source(reinterpret_cast<const char *>(original.bytes().data()),original.bytes().size());
    const auto mapping=admit_mapping(declaration,source,source_roles,q.limits,work);
    const auto ids=feature_ids(owner,q,work);
    // Bound immutable source payload before invoking bounded native
    // archive encoders. This is cooperative admission, not a process heap quota.
    for(const auto *publication:{&owner,&mapping_owner,&evidence_owner}) {
        std::size_t raw_bytes{};
        for(const auto &bytes:publication->raw_sources()) {
            work.take();bound(bytes.bytes().size(),q.limits.output_bytes-raw_bytes,"context owner raw-byte output budget");raw_bytes+=bytes.bytes().size();
        }
    }
    std::uint64_t next=1;
    for(const auto key:NativeRecord::source_keys(*owner.endpoint())) {
        require(key<std::numeric_limits<std::uint64_t>::max(),"context byte source key exhausted");next=std::max(next,key+1);
    }
    std::vector<RawByteSource> added;std::size_t archived_bytes{};
    const auto fresh=[&]() {require(next<std::numeric_limits<std::uint64_t>::max(),"context byte source key exhausted");return next++;};
    const auto archive=[&](const NativeRecord &publication) {
        work.check();
        (void)archive_identity(publication,q.limits.output_bytes-archived_bytes,stop);
        const auto key=fresh();auto bytes=publication.archive_source(key);work.check();
        bound(bytes.bytes().size(),q.limits.output_bytes-archived_bytes,"combined native context archive budget");
        archived_bytes+=bytes.bytes().size();const auto size=bytes.bytes().size();added.push_back(std::move(bytes));
        return N(N::ByteSpan{key,0,size,false});
    };
    work.check();
    auto mapping_archive=archive(mapping_owner);require(digest(added.back().bytes())==q.mapping_sha256,"feature context actual mapping archive SHA-256 mismatch");
    auto evidence_archive=archive(evidence_owner);
    const auto assertions=project_assertions(evidence.interpreted_source(),q,source_roles,work);
    const auto joined=join_features(ids,mapping,assertions,q.limits,work);
    const auto mapping_key=fresh();added.push_back(original.rekey(mapping_key));
    N::Array mapping_rows,features;
    for(std::size_t r=0;r<mapping.size();++r) {
        work.take();const auto &row=mapping[r];
        mapping_rows.elements.emplace_back(N::Object{{{"mapping_row",number(r)},{"source_record",number(row.source_record)},
            {"feature_id",N(row.feature)},{"evidence_id",N(row.identity)},{"role",N(row.role)},{"status",N(joined.row_status[r])},
            {"source_span",N(N::ByteSpan{mapping_key,row.begin,row.length,false})}}});
    }
    for(std::size_t f=0;f<joined.features.size();++f) {
        work.take();const auto &feature=joined.features[f];N::Array rows,matches;
        for(const auto row:feature.rows){work.take();rows.elements.push_back(number(row));}
        for(const auto &match:feature.matches) {
            work.take();const auto &assertion=assertions.at(match.assertion);const auto &participant=assertion.participants.at(match.participant);
            matches.elements.emplace_back(N::Object{{{"mapping_row",number(match.mapping_row)},
                {"assertion_key",N(std::to_string(assertion.key))},{"assertion_kind",N(std::to_string(assertion.kind))},
                {"assertion_selection_index",number(match.assertion)},{"native_occurrence_key",N(std::to_string(participant.occurrence))},
                {"role_position",number(participant.position)},{"active",N(participant.active)},
                {"role",N(participant.role)},{"evidence_identity",N(participant.identity)}}});
        }
        features.elements.emplace_back(N::Object{{{"feature_index",number(f)},{"feature_id",N(feature.id)},{"status",N(feature.status)},
            {"candidate_identities",number(feature.candidates)},{"mapping_rows",N(std::move(rows))},{"correspondences",N(std::move(matches))}}});
    }
    N result(N::Object{{{"schema",N(std::string("gadget.features.evidence-context.v1"))},
        {"analysis_sha256",N(q.analysis_sha256)},{"mapping_sha256",N(q.mapping_sha256)},{"evidence_digest",N(q.evidence_digest)},
        {"question",N(q.question)},{"mapping_schema",N(declaration.at("schema").text())},
        {"mapping_id",N(declaration.at("id").text())},{"mapping_origin",N(declaration.at("origin").text())},
        {"mapping_release",N(declaration.at("release").text())},{"mapping_source_sha256",N(declaration.at("sha256").text())},
        {"features",N(std::move(features))},{"mapping_rows",N(std::move(mapping_rows))},
        {"mapping_bytes",N(N::ByteSpan{mapping_key,0,source.size(),false})},
        {"analysis_identity",N(N::Object{{{"namespace",N(std::to_string(owner.source_namespace()))},
            {"cutoff",number(owner.endpoint()->sequence)}}})},{"mapping_publication",std::move(mapping_archive)},
        {"evidence_publication",std::move(evidence_archive)},
        {"scope",N(std::string("Explicit source correspondences only. Original analysis and evidence remain unchanged; no inferred biological equivalence, efficacy or causal relationship."))}}});
    work.check();auto operation=annotate_operation(owner,std::move(added),request,std::move(result));
    work.check();(void)archive_identity(operation,q.limits.output_bytes,stop);work.check();return operation;
}
} // namespace gadget
