#include "pubchem_evidence.hpp"
using namespace gadget::examples;
#include <gadget/record_view.hpp>
#include <iostream>

using namespace gadget;
using N=NativeRecordValue;
namespace {
std::size_t checks{};
void check(bool value,const char *message) {++checks;if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F callback) {bool refused{};try{callback();}catch(const std::exception &){refused=true;}check(refused,"invalid native PubChem publication/request was admitted");}
const std::string assay=R"({"PC_AssaySubmit":{"assay":{"descr":{"aid":{"id":100,"version":1},"name":"all typed values","description":["unchanged conditions"],"literal_metadata":{"rational":"literal provider label"},"results":[{"tid":1,"name":"float","type":1,"unit":5},{"tid":2,"name":"integer","type":2},{"tid":3,"name":"Boolean","type":3},{"tid":4,"name":"text","type":4}]}},"data":[{"sid":123,"outcome":1,"data":[{"tid":1,"value":{"fval":0.00}},{"tid":2,"value":{"ival":0}},{"tid":3,"value":{"bval":false}},{"tid":4,"value":{"sval":""}}]},{"sid":124}],"revoke":[124]}})";
const std::string csv="AID,SID,CID,Target Accession\n100,123,500,NP_1.1\n100,123,501,\n";
NativeRecord source() {
    auto a=std::as_bytes(std::span(assay.data(),assay.size())),c=std::as_bytes(std::span(csv.data(),csv.size()));
    return NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.pubchem.sources.native.v1"))},
        {"assay",N(N::ByteSpan{1,0,a.size(),false})},{"concise",N(N::ByteSpan{2,0,c.size(),false})}}}),new_native_record_namespace(),
        {RawByteSource(1,a),RawByteSource(2,c)});
}
Json request() {
    auto control=Json::parse(R"({"schema":"gadget.pubchem.admission.request.v1","source_id":"publication.case","aid":"100","assay_sha256":"","concise_sha256":"","provenance":{"release":"authored","selection":"two complete source rows","license":"authored test fixture","assay_url":"https://example.invalid/assay","concise_url":"https://example.invalid/concise"},"limits":{"source_bytes":1048576,"definitions":100,"substances":100,"readout_slots":10000,"concise_rows":100,"links":1000,"assertions":10000,"numeric_digits":100,"numeric_scale":1000,"native_records":100000,"native_bytes":67108864,"native_work":10000000}})");
    control.set("assay_sha256",Json::string("sha256:"+sha256(assay)));control.set("concise_sha256",Json::string("sha256:"+sha256(csv)));return control;
}
void contracts(const fs::path &parent) {
    fs::create_directories(parent);const auto directory=private_directory(parent,"pubchem-native-");
    const auto input=source(),admission=admit_record_view(request());const auto publication=admit_pubchem_evidence(input,admission);
    check(publication.at("schema").text()=="gadget.evidence.pubchem.native.v1","typed native publication missing");
    check(publication.at("records").at(0).at("original").at("literal_metadata").at("rational").text()=="literal provider label","provider metadata was reinterpreted as a numerical shorthand");
    check(publication.at("assay_bytes").byte_source().shares_owner(input.at("assay").byte_source()),"original immutable assay owner copied or replaced");
    check(NativeRecord::decode(publication.at("source_input").bytes()).encode()==input.encode(),"original input/history archive changed");
    check(NativeRecord::decode(publication.at("admission_request").bytes()).encode()==admission.encode(),"original admission request archive changed");
    publication.write(directory/"store.rfcxs");const EvidenceStore store(NativeRecord::read(directory/"store.rfcxs"),pubchem_evidence_adapter()),in_memory(publication,pubchem_evidence_adapter());
    check(store.digest()==in_memory.digest()&&store.native_publication().encode()==publication.encode(),"native store save/reopen changed source");
    const EvidenceStore repeated(admit_pubchem_evidence(input,admission),pubchem_evidence_adapter());
    check(repeated.digest()!=store.digest() && repeated.native_publication().source_namespace()!=publication.source_namespace(),"new publications reused a source identity");
    const auto independent_source=source();const EvidenceStore independent(admit_pubchem_evidence(independent_source,admission),pubchem_evidence_adapter());check(independent.native_publication().source_namespace()!=publication.source_namespace(),"independent sources share a mathematical namespace");check(independent.digest()!=store.digest(),"equal byte contents identified different actual native source owners");
    auto query=Json::parse(R"({"schema":"gadget.pubchem.query.v1","source_digest":"","operation":"readouts","max_records":8})");query.set("source_digest",Json::string(store.digest()));
    const auto operation=store.execute_record(admit_record_view(query));const auto rows=operation.at("result").at("records");
    check(rows.size()==8,"missing coordinates excluded from native query");std::size_t absent{},zero{},false_values{},empty{};
    for(std::size_t i=0;i<rows.size();++i) {
        const auto row=rows.at(i);check(row.at("source_bytes").byte_source().shares_owner(operation.parent()->at("assay_bytes").byte_source()),"native query lost raw source binding");
        if(row.at("support").text()=="absent")++absent;
        else {const auto value=row.at("original").at("value");if(value.has("fval")){++zero;check(row.at("number_token").text()=="0.00"&&value.at("fval").rational()==0,"zero exact value/literal changed");}
            if(value.has("ival"))check(value.at("ival").integer()==0,"integer zero lost");
            if(value.has("bval"))false_values+=!value.at("bval").boolean();
            if(value.has("sval"))empty+=value.at("sval").text().empty();}
    }
    check(absent==4&&zero==1&&false_values==1&&empty==1,"native typed value support changed");
    check(render_record_view(operation.at("request")).dump()==query.dump(),"native query request/history changed");
    const auto replayed=NativeRecord::decode(operation.encode());check(render_record_view(replayed).dump()==render_record_view(operation).dump(),"native query save/replay changed");
    check(operation.parent()&&operation.parent()->at("schema").text()=="gadget.evidence.pubchem.native.v1","query detached actual source history");
    auto context=Json::parse(R"({"schema":"gadget.pubchem.context.request.v1","source_digest":"","question":"Which values were actually reported?","assertion_keys":[],"max_records":1})");context.set("source_digest",Json::string(store.digest()));context.at("assertion_keys").push(Json::string(rows.at(0).at("native_assertion_key").text()));
    const auto selected=store.source_context(admit_record_view(context));check(selected.at("result").at("records").size()==1&&!selected.at("result").at("availability").at("historical_availability_established").boolean(),"typed source context changed selection or invented historical time");
    query.set("tid",Json::string("1"));query.set("sid",Json::string("123"));check(store.execute_record(admit_record_view(query)).at("result").at("records").size()==1,"native SID/TID selection failed");
    query.set("max_records",Json::integer(0));rejects([&]{(void)store.execute_record(admit_record_view(query));});
    auto control=request();control.set("assay_sha256",Json::string("wrong"));rejects([&]{(void)admit_pubchem_evidence(input,admit_record_view(control));});
    auto small=request();small.at("limits").set("native_bytes",Json::integer(1024));rejects([&]{(void)admit_pubchem_evidence(input,admit_record_view(small));});
    std::stop_source stopped;stopped.request_stop();rejects([&]{(void)admit_pubchem_evidence(input,admission,stopped.get_token());});
    auto changed=publication.materialize_value();auto &members=std::get<N::Object>(changed.data).members;
    for(auto &[name,value]:members)if(name=="receipt")for(auto &[field,item]:std::get<N::Object>(value.data).members)if(field=="source_digest")item=N(std::string("wrong source"));
    const auto forged=publication.derive(std::move(changed));rejects([&]{const EvidenceStore refused(forged,pubchem_evidence_adapter());});
}
} // namespace
int main(int argc,char **argv) {try{if(argc!=2)throw std::invalid_argument("test requires generated-data parent");contracts(argv[1]);std::cout<<checks<<" native PubChem publication/replay checks passed\n";return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
