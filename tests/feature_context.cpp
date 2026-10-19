#include "feature_context_fixture.hpp"
#include "feature_identity_fixture.hpp"
#include <algorithm>
#include <iostream>
#include "pubchem_evidence.hpp"
using namespace gadget::examples;
#define GADGET_CONTEXT_PUBCHEM_AVAILABLE 1

using namespace gadget;
namespace {
std::size_t checks{};
void check(bool ok,const char *message){++checks;if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F call){bool refused{};try{call();}catch(const std::exception &){refused=true;}check(refused,"invalid native feature context accepted");}
void contracts(const std::filesystem::path &root) {
    using N=NativeRecordValue;
    std::filesystem::remove_all(root);const auto f=gadget::test::make_feature_context_fixture(root);
    const auto before=f.analysis.encode(),mapping_before=f.mapping_source.encode(),evidence_before=f.evidence->native_publication().encode(),request_before=f.request.encode();
    const auto resources=CpuResources::local();
    const auto operation=contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,f.request,resources);
    const auto result=operation.at("result");
    check(result.at("schema").text()=="gadget.features.evidence-context.v1","context result schema");
    check(f.analysis.encode()==before&&f.mapping_source.encode()==mapping_before&&f.evidence->native_publication().encode()==evidence_before&&f.request.encode()==request_before,"context mutated an input publication");
    const auto parent=operation.parent();
    check(parent&&parent->type()==f.analysis.type()&&
          rfc::composite::declared_value_equal(parent->value(),f.analysis.value()),
          "actual unchanged analysis parent declaration lost");
    const auto retained_sources=operation.raw_sources();
    for(const auto &original:f.analysis.raw_sources())
        check(std::any_of(retained_sources.begin(),retained_sources.end(),[&](const auto &retained){
            return retained.key()==original.key()&&retained.shares_owner(original);
        }),"analysis original byte owner was replaced");
    check(rfc::storage::encode_rfcx(operation.endpoint()->source).bytes==rfc::storage::encode_rfcx(f.analysis.endpoint()->source).bytes,"evidence complex bindings transferred to analysis complex");
    check(result.at("analysis_identity").at("namespace").text()==std::to_string(f.analysis.source_namespace()),"analysis identity changed");
    check(NativeRecord::decode(result.at("mapping_publication").bytes()).encode()==mapping_before,"mapping archive changed");
    check(NativeRecord::decode(result.at("evidence_publication").bytes()).encode()==evidence_before,"evidence source roles/values/units/censoring changed");
    check(result.at("mapping_bytes").byte_source().shares_owner(f.mapping_source.at("bytes").byte_source()),"original mapping immutable byte owner replaced");
    check(render_record_view(operation.at("request")).dump()==render_record_view(f.request).dump(),"actual admitted request not retained");
    const auto features=result.at("features");
    check(features.size()==6&&features.at(0).at("correspondences").size()==2&&features.at(1).at("correspondences").size()==2,"repeated feature/map occurrences dropped");
    check(features.at(2).at("status").text()=="unmapped"&&features.at(3).at("status").text()=="mapped","feature versions or PAR_Y normalized");
    check(features.at(4).at("status").text()=="ambiguous"&&features.at(4).at("correspondences").size()==2&&features.at(5).at("status").text()=="unresolved","ambiguity or unresolved mapping hidden");
    const auto rows=result.at("mapping_rows");check(rows.size()==7&&rows.at(6).at("status").text()=="unused","unused original mapping occurrence hidden");
    const auto raw=f.mapping_source.at("bytes").bytes();
    for(std::size_t r=0;r<rows.size();++r) {
        const auto bytes=rows.at(r).at("source_span").bytes();
        const auto original=std::string_view(reinterpret_cast<const char *>(raw.data()),raw.size());
        const auto span=std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size());
        check(original.find(span)!=std::string_view::npos&&span.ends_with('\n'),"mapping span lost original row bytes");
    }
    const auto retained=NativeRecord::decode(result.at("evidence_publication").bytes());
    for(std::size_t i=0;i<features.size();++i)for(std::size_t m=0;m<features.at(i).at("correspondences").size();++m) {
        const auto match=features.at(i).at("correspondences").at(m);
        const auto cell=rfc::storage::StableCellKey{{1},std::stoull(match.at("assertion_key").text())};
        check(rfc::storage::StableSourceIndex(retained.endpoint()->source).cell(cell).has_value(),
              "correspondence lost its actual assertion owner");
    }
    operation.write(root/"context.rfcxs");const auto reopened=NativeRecord::read(root/"context.rfcxs");
    const auto reopened_parent=reopened.parent();
    check(reopened.encode()==operation.encode()&&reopened_parent&&
          reopened_parent->type()==f.analysis.type()&&
          rfc::composite::declared_value_equal(reopened_parent->value(),f.analysis.value()),
          "native context save/reopen or retained parent declaration changed");
    const auto replay=contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,f.request,resources);
    // New operations receive fresh declaration identities. Replay compares the
    // source bound result and retained owners, not a newly authored layout hash.
    const auto replay_result=replay.at("result");
    for(const auto *member:{"schema","analysis_sha256","mapping_sha256","evidence_digest",
                           "question","mapping_schema","mapping_id","mapping_origin",
                           "mapping_release","mapping_source_sha256","scope"})
        check(replay_result.at(member).text()==result.at(member).text(),"replay changed source pins or provenance");
    for(const auto *member:{"features","mapping_rows"})
        check(render_record_view(replay_result.at(member)).dump_compact()==
              render_record_view(result.at(member)).dump_compact(),"replay changed feature/mapping occurrences");
    check(std::ranges::equal(replay_result.at("mapping_bytes").bytes(),raw)&&
          replay_result.at("mapping_bytes").byte_source().shares_owner(f.mapping_source.at("bytes").byte_source()),
          "replay changed original mapping bytes/owner");
    const auto replay_rows=replay_result.at("mapping_rows");
    for(std::size_t r=0;r<rows.size();++r) {
        const auto original_span=std::get<N::ByteSpan>(rows.at(r).at("source_span").materialize_value().data);
        const auto replay_span=std::get<N::ByteSpan>(replay_rows.at(r).at("source_span").materialize_value().data);
        check(original_span.offset==replay_span.offset&&original_span.length==replay_span.length&&
              std::ranges::equal(rows.at(r).at("source_span").bytes(),replay_rows.at(r).at("source_span").bytes()),
              "replay changed an original mapping source interval");
    }
    for(const auto *member:{"mapping_publication","evidence_publication"})
        check(std::ranges::equal(replay_result.at(member).bytes(),result.at(member).bytes()),
              "replay changed an exact original publication archive");
    check(std::ranges::equal(operation.at("request_publication").bytes(),request_before)&&
          std::ranges::equal(replay.at("request_publication").bytes(),request_before),
          "replay replaced the actual original request owner");
    // Exact ownership pins refuse another native publication even with equal views.
    const auto same_view=admit_record_view(render_record_view(f.analysis));
    rejects([&]{(void)contextualize_feature_evidence(same_view,*f.evidence,f.mapping_source,f.request,resources);});
    for(const auto *name:{"analysis_sha256","mapping_sha256","evidence_digest"}) {
        auto q=render_record_view(f.request);q.set(name,Json::string("sha256:"+sha256("wrong actual owner")));
        rejects([&]{(void)contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,admit_record_view(q),resources);});
    }
    for(const auto *name:{"source_bytes","features","mapping_rows","assertions","participants","matches","work","output_bytes"}) {
        auto q=render_record_view(f.request);q.at("limits").set(name,Json::integer(0));
        rejects([&]{(void)contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,admit_record_view(q),resources);});
    }
    auto q=render_record_view(f.request);auto paths=Json::array();paths.push(Json::string("not_a_feature_array"));q.set("feature_path",paths);
    rejects([&]{(void)contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,admit_record_view(q),resources);});
    auto changed=f.mapping_source.materialize_value();for(auto &[name,value]:std::get<N::Object>(changed.data).members)if(name=="sha256")value=N("sha256:"+sha256("wrong original mapping"));
    const auto wrong_source=f.mapping_source.derive(changed);q=render_record_view(f.request);q.set("mapping_sha256",Json::string(test::context_archive_sha(wrong_source)));
    rejects([&]{(void)contextualize_feature_evidence(f.analysis,*f.evidence,wrong_source,admit_record_view(q),resources);});
    auto metadata=f.mapping_source.materialize_value();
    auto sources=f.mapping_source.raw_sources();std::uint64_t text_key=1;
    for(const auto &source:sources)text_key=std::max(text_key,source.key()+1);
    const std::string shared(4096,'X');
    sources.emplace_back(text_key,std::as_bytes(std::span(shared)));
    for(auto &[name,value]:std::get<N::Object>(metadata.data).members)
        if(name=="id"||name=="release"||name=="source_url"||name=="license")
            value=N(N::ByteSpan{text_key,0,shared.size(),true});
    const auto aliased_mapping=NativeRecord::create(metadata,new_native_record_namespace(),std::move(sources));
    auto declaration_limit=render_record_view(f.request);
    declaration_limit.set("mapping_sha256",Json::string(test::context_archive_sha(aliased_mapping)));
    declaration_limit.at("limits").set("output_bytes",Json::integer(65536));
    bool metadata_refusal{};
    try {(void)contextualize_feature_evidence(f.analysis,*f.evidence,aliased_mapping,admit_record_view(declaration_limit),resources);}
    catch(const std::exception &error) {
        metadata_refusal=std::string_view(error.what()).find("mapping declaration working byte admission")!=std::string_view::npos;
    }
    check(metadata_refusal,"shared mapping declaration spans expanded beyond working admission");
    const auto repeated=test::repeated_feature_identity_source();
    auto bounded=render_record_view(f.request);
    bounded.set("analysis_sha256",Json::string(test::context_archive_sha(repeated)));
    bounded.set("feature_path",Json::parse(R"(["features"])"));
    bounded.at("limits").set("features",Json::integer(512));
    bounded.at("limits").set("output_bytes",Json::integer(1048576));
    bool early_refusal{};
    try {(void)contextualize_feature_evidence(repeated,*f.evidence,f.mapping_source,admit_record_view(bounded),resources);}
    catch(const std::exception &error) {
        early_refusal=std::string_view(error.what()).find("identity working byte admission")!=std::string_view::npos;
    }
    check(early_refusal,"shared native feature labels expanded beyond working admission");
    std::stop_source stopped;stopped.request_stop();rejects([&]{(void)contextualize_feature_evidence(f.analysis,*f.evidence,f.mapping_source,f.request,resources,stopped.get_token());});
}
#ifdef GADGET_CONTEXT_PUBCHEM_AVAILABLE
void typed_pubchem_contracts(const fs::path &root) {
    fs::create_directories(root);
    using N=NativeRecordValue;
    const std::string assay=R"({"PC_AssaySubmit":{"assay":{"descr":{"aid":{"id":100},"name":"authored source","results":[{"tid":1,"name":"modifier","type":4},{"tid":2,"name":"reported number","type":1},{"tid":3,"name":"Boolean","type":3}]}},"data":[{"sid":100,"outcome":1,"data":[{"tid":1,"value":{"sval":"<"}},{"tid":2,"value":{"fval":0.000}},{"tid":3,"value":{"bval":false}}]},{"sid":101,"data":[]}]}})";
    const std::string concise="AID,SID,CID,Target Accession\n100,100,100,NP_1.1\n";
    const auto input=NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.pubchem.sources.native.v1"))},
        {"assay",N(N::ByteSpan{1,0,assay.size(),false})},{"concise",N(N::ByteSpan{2,0,concise.size(),false})}}}),
        new_native_record_namespace(),{RawByteSource(1,std::as_bytes(std::span(assay))),RawByteSource(2,std::as_bytes(std::span(concise)))});
    auto admission=Json::parse(R"({"schema":"gadget.pubchem.admission.request.v1","source_id":"context.authored.pubchem","aid":"100","assay_sha256":"","concise_sha256":"","provenance":{"release":"authored.v1","selection":"complete authored source","license":"authored test data","assay_url":"offline:assay","concise_url":"offline:concise"},"limits":{"source_bytes":1048576,"definitions":100,"substances":100,"readout_slots":10000,"concise_rows":100,"links":1000,"assertions":10000,"numeric_digits":100,"numeric_scale":1000,"native_records":100000,"native_bytes":67108864,"native_work":10000000}})");
    admission.set("assay_sha256",Json::string("sha256:"+sha256(assay)));admission.set("concise_sha256",Json::string("sha256:"+sha256(concise)));
    const auto publication=admit_pubchem_evidence(input,admit_record_view(admission));const EvidenceStore evidence(publication,pubchem_evidence_adapter());
    const auto analysis=admit_record_view(Json::parse(R"({"features":[{"id":"A"},{"id":"B"},{"id":"C"}]})"));
    const std::string mapping="feature_id\tevidence_id\trole\nA\tpubchem:sid:100\tleft\nA\tpubchem:cid:100\tright\nB\tpubchem:sid:404\tleft\nC\tpubchem:target_accession:NP_1.1\tright\nall.negative\tpubchem:sid:100\tleft\nall.positive\tpubchem:target_accession:NP_1.1\tright\n";
    const auto source=NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.features.mapping.native.v1"))},
        {"id",N(std::string("explicit.pubchem-map"))},{"origin",N(std::string("explicit_declared"))},
        {"release",N(std::string("authored.v1"))},{"source_url",N(std::string("offline:map"))},{"license",N(std::string("authored test data"))},
        {"format",N(std::string("tsv"))},{"sha256",N("sha256:"+sha256(mapping))},{"bytes",N(N::ByteSpan{1,0,mapping.size(),false})}}}),
        new_native_record_namespace(),{RawByteSource(1,std::as_bytes(std::span(mapping)))});
    auto query=Json::parse(R"({"schema":"gadget.features.evidence-context.request.v1","analysis_sha256":"","mapping_sha256":"","evidence_digest":"","feature_path":["features"],"feature_id_path":["id"],"assertion_keys":[],"question":"Which exact typed source participants were explicitly mapped?","limits":{"source_bytes":65536,"features":100,"mapping_rows":100,"assertions":100,"participants":1000,"matches":1000,"work":10000000,"output_bytes":67108864}})");
    query.set("analysis_sha256",Json::string(test::context_archive_sha(analysis)));query.set("mapping_sha256",Json::string(test::context_archive_sha(source)));query.set("evidence_digest",Json::string(evidence.digest()));
    for(std::size_t i=0;i<publication.at("records").size();++i){const auto value=publication.at("records").at(i).materialize_value();query.at("assertion_keys").push(Json::string(std::to_string(std::get<N::Owned>(value.data).cell->key)));}
    const auto operation=contextualize_feature_evidence(analysis,evidence,source,admit_record_view(query));const auto result=operation.at("result");
    check(result.at("features").at(0).at("status").text()=="ambiguous"&&result.at("features").at(1).at("status").text()=="unresolved"&&result.at("features").at(2).at("status").text()=="mapped","typed source exact context statuses changed");
    const auto matches=result.at("features").at(0).at("correspondences");bool sid{},cid{};
    for(std::size_t i=0;i<matches.size();++i){const auto row=matches.at(i);const auto role=row.at("role").text(),identity=row.at("evidence_identity").text();sid|=role=="left"&&identity=="pubchem:sid:100";cid|=role=="right"&&identity=="pubchem:cid:100";check((role=="left"&&identity=="pubchem:sid:100")||(role=="right"&&identity=="pubchem:cid:100"),"typed context silently identified SID with CID or ignored actual role");}
    check(sid&&cid,"explicit SID and CID source links missing");
    check(NativeRecord::decode(result.at("evidence_publication").bytes()).encode()==publication.encode(),"typed readouts/modifiers/outcomes changed during context join");
    const auto reopened=NativeRecord::decode(operation.encode());check(reopened.encode()==operation.encode(),"typed context native save/reopen changed");
    publication.write(root/"evidence.rfcxs");analysis.write(root/"analysis.rfcxs");source.write(root/"mapping.rfcxs");
    json_object_object_del(query.raw(),"analysis_sha256");
    auto dispatch=Json::object();dispatch.set("schema",Json::string("gadget.result-operation.v1"));
    dispatch.set("operation",Json::string("feature_evidence_context"));dispatch.set("analysis_binding",Json::string("actual_result"));
    dispatch.set("analysis_bytes",Json::integer(67108864));dispatch.set("query",query);
    auto mapping_record=Json::object();mapping_record.set("path",Json::string(fs::absolute(root/"mapping.rfcxs").string()));
    mapping_record.set("sha256",Json::string(test::context_archive_sha(source)));mapping_record.set("max_bytes",Json::integer(67108864));
    dispatch.set("source_record",mapping_record);write_text(root/"dispatch.json",dispatch.dump());
    auto retrieval=Json::parse(R"({"schema":"gadget.pubchem.query.v1","source_digest":"","operation":"readouts","max_records":100})");
    retrieval.set("source_digest",Json::string(evidence.digest()));write_text(root/"query.json",retrieval.dump());
    write_text(root/"expected.json",render_record_view(result.at("features")).dump());

}
#endif

}
int main(int argc,char **argv){try{if(argc!=2)throw std::invalid_argument("native context tests require a test-data path");contracts(argv[1]);
#ifdef GADGET_CONTEXT_PUBCHEM_AVAILABLE
    typed_pubchem_contracts(fs::path(argv[1])/"typed");
#endif
    std::cout<<checks<<" native feature context owner/reopen/replay checks passed\n";return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
