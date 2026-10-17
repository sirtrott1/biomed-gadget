#include <gadget/homology_diagram.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/inference/native_homology_diagram_codec.hpp>
#include <fstream>
#include <iostream>

using namespace gadget;
namespace {
using N=NativeRecordValue;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
NativeRecord source(std::uint64_t space,long order){
    rfc::core::RelationSignature signature;
    const auto loop=signature.add_type("loop",rfc::core::BoundaryPolicy::DeclaredExact,2,2);
    const auto face=signature.add_type("face",rfc::core::BoundaryPolicy::DeclaredExact,1,3);
    rfc::core::RelationalComplexBuilder b(signature);
    b.add_grade_one(loop,{{42,true,true,false,1},{42,true,false,false,1}});
    b.add_higher(rfc::core::Grade{2},face,{{0,true,false,false,rfc::numeric::Rational{order}},
        {0,false,false,false,1},{0,true,false,false,0}});
    rfc::storage::TemporalEndpoint endpoint(b.freeze());
    for(std::size_t k=0;k<endpoint.source.grade_count();++k){rfc::storage::TemporalCarrier c;c.key=10+k;
        c.kind=rfc::storage::TemporalCarrierKind::SourceGrade;c.grade=rfc::core::Grade{static_cast<std::uint16_t>(k)};
        c.domain=rfc::numeric::IndexDomain::linear(endpoint.source.grade(c.grade).cell_count());endpoint.carriers.push_back(c);}
    const std::vector<std::byte> original{std::byte{0},std::byte{255},std::byte{42},std::byte{10}};
    return NativeRecord::create_from_endpoint(N(N::Object{{{"schema",N(std::string("gadget.homology.fixture.v1"))},
        {"original",N(N::ByteSpan{1,0,original.size(),false})}}}),space,{RawByteSource(1,original)},std::move(endpoint));
}
std::string specification(const std::string& product){return
    "{\"schema\":\"gadget.homology.diagram.v1\",\"product\":\""+product+"\",\"grade\":1,"
    "\"vertices\":[{\"source\":\"four.rfcxs\",\"carriers\":[10,11,12]},{\"source\":\"two.rfcxs\",\"carriers\":[10,11,12]}],"
    "\"arrows\":[{\"forward\":true,\"maps\":[{\"rows\":1,\"columns\":1,\"entries\":[[0,0,\"1\"]]},"
    "{\"rows\":1,\"columns\":1,\"entries\":[[0,0,\"1\"]]},{\"rows\":1,\"columns\":1,\"entries\":[[0,0,\"2\"]]}]}],"
    "\"evidence\":\"Authored exact torsion chain map, including retained inactive and zero occurrences\"}";}
}
int main(int argc,char** argv){try{
    check(argc==2,"destination required");const fs::path destination=fs::absolute(argv[1]);fs::create_directories(destination);
    fs::path root;
    do {root=destination/std::to_string(new_native_record_namespace());} while(!fs::create_directory(root));
    const std::vector<NativeRecord> sources{source(900,4),source(901,2)};
    sources[0].write(root/"four.rfcxs");sources[1].write(root/"two.rfcxs");
    auto resources=CpuResources::local();
    for(const std::string product:{"rational_zigzag","integral_modules"}){
        const auto text=specification(product);{std::ofstream output(root/(product+".json"));output<<text;}
        const auto request=admit_record_view(Json::parse(text));
        const auto record=homology_diagram_record(request,sources,resources);record.write(root/(product+".rfcxs"));
        const auto packet=record.at("view").at("provider").bytes();
        const auto provider=rfc::inference::decode_native_homology_diagram(packet,native_data_limits().temporal,rfc::runtime::ExecutionLease{},native_accession_limits());
        const auto encoded=rfc::inference::encode_native_homology_diagram(provider,source_byte_limit());
        check(std::equal(packet.begin(),packet.end(),encoded.begin(),encoded.end()),"generic SDK cold replay changed complete application source channels");
        const auto view=homology_diagram_view(record);
        check(view.at("product").text()==product,"product selector differs");
        if(product=="integral_modules"){
            check(view.at("vertices").at(0).at("torsion").at(0).integer()==4 && view.at("vertices").at(1).at("torsion").at(0).integer()==2,"integral torsion lost");
            check(view.at("arrows").at(0).at("smith_map").at("entries").size()==1,"integral map lost");
            auto endpoint=*record.endpoint();const auto key=record.at("view").at("arrows").at(0).at("smith_map").at("tensor").integer().get_ui();bool changed=false;
            for(auto& tensor:endpoint.tensors)if(tensor.key==key){
                tensor.tensor=rfc::numeric::SparseRationalTensor(tensor.tensor.domain(),tensor.tensor.support(),
                    rfc::numeric::SparseRationalVector(1,{{0,rfc::numeric::Rational{2}}}));changed=true;}
            check(changed,"published induced map field absent");
            endpoint.declarations.reset();endpoint.codes.reset();
            const auto false_result=NativeRecord::create_from_endpoint(record.unbound_value(),new_native_record_namespace(),record.raw_sources(),std::move(endpoint));
            bool refused=false;try{(void)replay_homology_diagram(false_result,resources);}catch(const std::invalid_argument& error){
                if(std::string_view(error.what()).find("native result field differs")==std::string_view::npos)throw;
                refused=true;}
            check(refused,"false numerical result accepted against its original sources");
        }else check(view.at("intervals").size()==0,"rational torsion quotient differs");
        const auto cold=replay_homology_diagram(NativeRecord::read(root/(product+".rfcxs")),resources);
        check(homology_diagram_view(cold).dump_compact()==view.dump_compact(),"complete replay result differs");
        auto changed=record.unbound_value();
        for(auto& [key,value]:std::get<N::Object>(changed.data).members)if(key=="view")
            for(auto& [name,member]:std::get<N::Object>(value.data).members)if(name=="grade")member=N(rfc::numeric::Integer{0});
        bool contradicted=false;try{(void)replay_homology_diagram(record.derive(std::move(changed)),resources);}catch(const std::invalid_argument& error){
            if(std::string_view(error.what()).find("typed result declaration differs")==std::string_view::npos)throw;
            contradicted=true;}
        check(contradicted,"changed result declaration accepted");
        const auto routed=run_homology_diagram(root/(product+".json"),resources);
        check(homology_diagram_view(routed).dump_compact()==view.dump_compact(),"relative input paths differ");
        for(std::size_t i=0;i<sources.size();++i){const auto archive=record.at("view").at("vertices").at(i).at("source").at("publication").bytes();
            const auto original=sources[i].encode();check(std::equal(original.begin(),original.end(),archive.begin(),archive.end()),"full original native source changed");}
        auto stop=std::make_shared<std::atomic_bool>(true);bool refused=false;
        try{(void)homology_diagram_record(request,sources,resources,rfc::runtime::ExecutionLease({},stop));}catch(const rfc::runtime::ExecutionCancelled&){refused=true;}check(refused,"cancellation ignored");
        rfc::runtime::ExecutionHostPolicy tiny;tiny.limits.rational_slots=1;refused=false;
        try{(void)homology_diagram_record(request,sources,resources,rfc::runtime::ExecutionLease(tiny));}catch(const std::length_error&){refused=true;}check(refused,"resource refusal ignored");
    }
    auto positive=Json::parse(specification("rational_zigzag"));positive.set("grade",Json::integer(0));
    {std::ofstream output(root/"rational-positive.json");output<<positive.dump_compact();}
    const auto persistent=homology_diagram_record(admit_record_view(positive),sources,resources);
    const auto intervals=homology_diagram_view(persistent).at("intervals");
    check(intervals.size()==1 && intervals.at(0).at("begin").integer()==0 && intervals.at(0).at("end").integer()==1 && intervals.at(0).at("multiplicity").integer()==1 && intervals.at(0).at("generators").size()==2,"native interval generator publication differs");
    (void)replay_homology_diagram(persistent,resources);
    const auto hardware=resources->readmit();rfc::runtime::CpuResourcePolicy policy;policy.worker_limit=1;
    policy.declared_buffer_byte_limit=resources->coordinator().available_declared_buffer_bytes();
    const auto one=std::make_shared<CpuResources>(hardware,policy,[hardware]{return hardware;});
    const auto request=admit_record_view(Json::parse(specification("integral_modules")));
    {auto held=one->reserve(1);bool refused=false;
        try{(void)homology_diagram_record(request,sources,one);}catch(const rfc::runtime::CpuResourceBusy&){refused=true;}check(refused,"busy caller reservation ignored");}
    (void)homology_diagram_record(request,sources,one);check(one->coordinator().available_workers()==1,"caller worker not released");
    auto receipt=Json::object();receipt.set("status",Json::string("PASS"));receipt.set("fixtures",Json::string(root.string()));
    std::cout<<receipt.dump_compact()<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
