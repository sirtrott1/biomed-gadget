#pragma once
#include <gadget/feature_context.hpp>
#include <gadget/record_view.hpp>

namespace gadget::test {
struct FeatureContextFixture {
    NativeRecord analysis,mapping_source;
    std::shared_ptr<EvidenceStore> evidence;
    NativeRecord request;
};
inline std::string context_archive_sha(const NativeRecord &record) {
    const auto bytes=record.owning_publication().encode();
    return "sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size()));
}
inline FeatureContextFixture make_feature_context_fixture(const std::filesystem::path &root) {
    namespace fs=std::filesystem;
    fs::create_directories(root/"snapshot");
    auto body=Json::parse(R"({"activities":[{"activity_id":1,"molecule_chembl_id":"CHEMBL_C1","target_chembl_id":"CHEMBL_T1","assay_chembl_id":"CHEMBL_A1","document_chembl_id":"CHEMBL_D1","standard_type":"Ki","standard_relation":"<","standard_units":"nM","standard_value":"0"},{"activity_id":2,"molecule_chembl_id":"CHEMBL_C2","target_chembl_id":"CHEMBL_T2","assay_chembl_id":"CHEMBL_A2","document_chembl_id":"CHEMBL_D2","standard_type":"Kd","standard_relation":"=","standard_units":"uM","standard_value":null}]})");
    const auto bytes=body.dump(),headers=std::string("authored offline mapping fixture\n");
    write_text(root/"snapshot/activity.json",bytes);write_text(root/"snapshot/activity.headers.txt",headers);
    auto manifest=Json::parse(R"({"schema":"gadget.database.snapshot.v1","id":"authored-context-evidence","records":[{"request":{"id":"activity","provider":"chembl","resource":"activity","identifier":"P00519"},"url":"offline:authored-context","path":"activity.json","sha256":"","headers_path":"activity.headers.txt","headers_sha256":""}]})");
    manifest.at("records").at(0).set("sha256",Json::string("sha256:"+sha256(bytes)));
    manifest.at("records").at(0).set("headers_sha256",Json::string("sha256:"+sha256(headers)));
    manifest.set("snapshot_sha256",Json::string("sha256:"+sha256(manifest.dump())));
    write_text(root/"snapshot/snapshot.json",manifest.dump());
    import_evidence(root/"snapshot",root/"store",3);
    auto evidence=std::make_shared<EvidenceStore>(root/"store");
    const auto analysis=admit_record_view(Json::parse(R"({"schema":"authored.molecular.result.v1","result":{"view":{"features":[{"id":"ENSG000001.1","coefficient":{"rational":"-3/7"}},{"id":"ENSG000001.1","coefficient":{"rational":"0"}},{"id":"ENSG000001.2"},{"id":"ENSG000001.1_PAR_Y"},{"id":"ambiguous"},{"id":"unresolved"}]}}})"));
    const std::string mapping="feature_id\tevidence_id\trole\nENSG000001.1\tchembl:CHEMBL_T1\ttarget_or_right\nENSG000001.1\tchembl:CHEMBL_T1\ttarget_or_right\nENSG000001.1_PAR_Y\tchembl:CHEMBL_T2\ttarget_or_right\nambiguous\tchembl:CHEMBL_T1\ttarget_or_right\nambiguous\tchembl:CHEMBL_T2\ttarget_or_right\nunresolved\tchembl:missing\ttarget_or_right\nunused\tchembl:CHEMBL_C1\tcompound_or_left\n";
    using N=NativeRecordValue;const RawByteSource raw(1,std::as_bytes(std::span(mapping)));
    const auto mapping_source=NativeRecord::create(N(N::Object{{
        {"schema",N(std::string("gadget.features.mapping.native.v1"))},{"id",N(std::string("authored.feature-map"))},
        {"origin",N(std::string("explicit_declared"))},{"release",N(std::string("authored.v1"))},
        {"source_url",N(std::string("offline:authored-mapping"))},{"license",N(std::string("authored test data"))},
        {"format",N(std::string("tsv"))},{"sha256",N("sha256:"+sha256(mapping))},
        {"bytes",N(N::ByteSpan{1,0,mapping.size(),false})}}}),new_native_record_namespace(),{raw});
    auto query=Json::parse(R"({"schema":"gadget.features.evidence-context.request.v1","analysis_sha256":"","mapping_sha256":"","evidence_digest":"","feature_path":["result","view","features"],"feature_id_path":["id"],"assertion_keys":[],"question":"Which explicitly mapped targets occur in these original source assertions?","limits":{"source_bytes":65536,"features":100,"mapping_rows":100,"assertions":10,"participants":100,"matches":100,"work":10000000,"output_bytes":67108864}})");
    query.set("analysis_sha256",Json::string(context_archive_sha(analysis)));
    query.set("mapping_sha256",Json::string(context_archive_sha(mapping_source)));
    query.set("evidence_digest",Json::string(evidence->digest()));
    const auto records=evidence->native_publication().at("records");
    for(std::size_t i=0;i<records.size();++i) {
        const auto value=records.at(i).materialize_value();const auto &owned=std::get<N::Owned>(value.data);
        query.at("assertion_keys").push(Json::string(std::to_string(owned.cell->key)));
    }
    return {analysis,mapping_source,evidence,admit_record_view(query)};
}
} // namespace gadget::test
