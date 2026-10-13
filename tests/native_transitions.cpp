#include "../src/research_factored.hpp"
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <gadget/resources.hpp>
#include <rfc/storage/rfcx.hpp>
#include <atomic>
#include <iostream>

using namespace gadget;
using Rational=rfc::numeric::Rational;
using Matrix=rfc::numeric::ExactMatrix;
namespace {
std::size_t checks{};
void check(bool value,const char *message) {++checks;if(!value)throw std::runtime_error(message);}
const rfc::storage::TemporalTensor &tensor(const NativeRecord &owner,const NativeRecord &reference) {
    const auto key=reference.at("tensor").integer();
    for(const auto &field:owner.endpoint()->tensors)if(field.key==key)return field;
    throw std::runtime_error("native transition tensor absent");
}
Vector values(const NativeRecord &owner,const NativeRecord &reference) {
    const auto &field=tensor(owner,reference).tensor;
    check(field.support().count()==field.domain().logical_size(),"transition complete support lost");
    Vector result;
    for(std::size_t i=0;i<field.domain().logical_size();++i)result.push_back(*field.state(i));
    return result;
}
Matrix matrix(const NativeRecord &owner,const NativeRecord &reference) {
    const auto &field=tensor(owner,reference);
    check(field.tensor.support().count()==field.tensor.domain().logical_size(),"complete map support lost");
    return rfc::storage::temporal_tensor_matrix(field);
}
void audit(const NativeRecord &owner,const NativeRecord &value) {
    if(value.kind()==NativeRecordKind::Object) {
        if(value.has("schema")&&value.at("schema").kind()==NativeRecordKind::Text&&
           value.at("schema").text()=="gadget.scientific.tensor.view.v1") {
            const auto &field=tensor(owner,value);
            check(field.output.has_value(),"native transition output carrier absent");
            check(field.tensor.support().count()==field.tensor.domain().logical_size(),"transition field omitted supported zero");
            return;
        }
        for(const auto &name:value.names())audit(owner,value.at(name));
    } else if(value.kind()==NativeRecordKind::Array)
        for(std::size_t i=0;i<value.size();++i)audit(owner,value.at(i));
}
void verify(const NativeRecord &publication,const ScientificStudy &study) {
    const auto owner=publication.at("result").at("fields"),view=owner.at("view"),binding=owner.at("bindings").at(0);
    audit(owner,view);audit(owner,owner.at("bindings"));
    const auto original=NativeRecord::decode(publication.at("result").at("source_study_publication").bytes());
    check(original.encode()==study.native_publication().encode(),"transition original source history changed");
    for(const auto &name:{"old_input","new_input"}) {
        const auto input=binding.at(name);
        const auto pin=study.field_accession(input.at("id").text());
        check(input.at("source_namespace").text()==std::to_string(pin.source_namespace())&&
              input.at("source_cutoff").integer()==pin.sequence(),"transition field producer cutoff changed");
        const auto lift=matrix(owner,input.at("realization"));
        Vector full;
        for(std::size_t i=0;i<pin.tensor()->domain().logical_size();++i)
            full.push_back(pin.tensor()->state(i).value_or(Rational{0}));
        check(lift.apply(full)==values(owner,input.at("values")),"transition source lift differs");
    }
    const auto old_snapshot=NativeRecord::decode(binding.at("old_source").at("publication").bytes());
    const auto new_snapshot=NativeRecord::decode(binding.at("new_source").at("publication").bytes());
    check(rfc::storage::encode_rfcx(owner.endpoint()->source).bytes==rfc::storage::encode_rfcx(new_snapshot.endpoint()->source).bytes,
          "new actual context source differs");
    check(old_snapshot.source_namespace()!=new_snapshot.source_namespace(),"independent context snapshots merged");
    const auto old=view.at("old"),next=view.at("new");
    check(tensor(owner,old.at("input")).output!=tensor(owner,next.at("input")).output,"old and new axes identified by dimension");
    for(const auto &endpoint:{old,next}) {
        auto sum=rfc::numeric::add(values(owner,endpoint.at("gradient")),values(owner,endpoint.at("curl")));
        sum=rfc::numeric::add(sum,values(owner,endpoint.at("harmonic")));
        check(sum==values(owner,endpoint.at("input")),"native sectors do not reconstruct");
    }
    auto total=rfc::numeric::add(values(owner,view.at("realization_term")),values(owner,view.at("innovation_term")));
    total=rfc::numeric::add(total,values(owner,view.at("operator_term")));
    check(total==values(owner,view.at("total")),"native Green terms do not reconstruct");
    Vector transport(values(owner,next.at("input")).size());
    for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b)
        transport=rfc::numeric::add(transport,values(owner,view.at("sector_transport").at(a).at(b)));
    check(transport==values(owner,view.at("old_transported")),"complete nine-sector field transport differs");
    if(!view.has("factored_maps")) {
        const auto j=matrix(owner,view.at("correspondence"));Matrix sum(j.rows(),j.cols());
        for(std::size_t a=0;a<3;++a)for(std::size_t b=0;b<3;++b) {
            const auto map=view.at("sector_maps").at(a).at(b);
            check(tensor(owner,map).input==tensor(owner,old.at("input")).output&&
                  tensor(owner,map).output==tensor(owner,next.at("input")).output,"sector map endpoint axes differ");
            sum=sum.add(matrix(owner,map));
        }
        check(sum==j,"complete sector maps do not reconstruct correspondence");
    } else {
        const auto nodes=view.at("factored_maps").at("nodes");
        for(std::size_t i=0;i<nodes.size();++i) {
            const auto raw=nodes.at(i);
            const auto node=raw.has("additions")?raw.at("additions"):raw;
            check(node.has("input_carrier")&&node.has("output_carrier"),"native factored axes absent");
            if(node.at("kind").text()=="matrix") {
                const auto base=raw.at("base");
                const auto &field=tensor(owner,base);
                check(field.input==node.at("input_carrier").integer().get_ui()&&
                      field.output==node.at("output_carrier").integer().get_ui(),"factored matrix actual carrier differs");
            }
            for(const auto &name:{"input","system"})if(node.has(name))check(node.at(name).integer()<i,"factored declaration order differs");
            for(const auto &name:{"inputs","application_order"})if(node.has(name))
                for(std::size_t j=0;j<node.at(name).size();++j)check(node.at(name).at(j).integer()<i,"factored declaration cycle");
        }
    }
}
}
int main(int argc,char **argv) try {
    check(argc==4,"transition study, requests and output required");
    auto initial=std::make_shared<const ScientificStudy>(argv[1]);
    auto study=std::make_shared<const ScientificStudy>(NativeRecord::annotate(initial->native_publication().unbound_value(),initial->native_publication()));
    initial.reset();const auto requests=admit_record_view(Json::parse(read_text(argv[2])));fs::create_directories(argv[3]);
    auto session=std::make_unique<ScientificSession>(study);std::vector<NativeRecord> retained;
    {
        scientific_detail::ResultFields fields(study->native_publication());
        research_detail::NativeFactoredTemporalMaps maps({},&fields,{});
        const auto a=fields.axis(2),b=fields.axis(2),c=fields.axis(2);
        const auto first=maps.matrix(Matrix::identity(2),a,b),second=maps.matrix(Matrix::identity(2),c,a);
        for(std::size_t scenario=0;scenario<2;++scenario) {
            bool refused=false;
            try {
                if(scenario==0)(void)maps.compose({first,second});
                if(scenario==1)(void)maps.sum({first,maps.transpose(first)});
            } catch(const std::invalid_argument &){refused=true;}
            check(refused,"equal dimensions identified independent factored carriers");
        }
        const auto inverse=maps.solve(first);
        check(inverse.at("input_carrier").count_value()==b&&inverse.at("output_carrier").count_value()==a,
              "declared inverse did not reverse its independent carriers");
    }
    for(const auto &name:requests.names()) {
        std::cerr<<"native transition "<<name<<" begin\n"<<std::flush;
        const auto result=session->execute_record(requests.at(name));verify(result,*study);
        const fs::path path=fs::path(argv[3])/(name+".rfcxs");result.write(path);
        const auto reopened=NativeRecord::read(path);verify(reopened,*study);
        const auto replay=replay_scientific_record(reopened);verify(replay,*study);
        check(reopened.endpoint()->tensors.size()==replay.endpoint()->tensors.size(),"replay field count changed");
        for(std::size_t i=0;i<reopened.endpoint()->tensors.size();++i)
            check(rfc::storage::temporal_endpoint_detail::tensor_equal(reopened.endpoint()->tensors[i],replay.endpoint()->tensors[i]),"replay complete tensor changed");
        retained.push_back(reopened);std::cerr<<"native transition "<<name<<" replay PASS\n"<<std::flush;
    }
    const auto request=render_record_view(requests.at(requests.names().front()));
    const auto maximum=std::numeric_limits<std::size_t>::max();
    for(std::size_t scenario=0;scenario<7;++scenario) {
        rfc::runtime::ExecutionHostPolicy policy;policy.name="native-transition-contract";policy.limits={maximum,maximum,maximum,maximum,maximum};
        auto stop=std::make_shared<std::atomic_bool>(scenario==0);
        if(scenario==1)policy.limits.records=0;
        if(scenario==2)policy.limits.rational_slots=0;
        if(scenario==3)policy.limits.sparse_entries=0;
        if(scenario==4)policy.limits.scalar_work=0;
        if(scenario==5)policy.limits.integer_bits=0;
        if(scenario==6)policy.allow_allocation=false;
        bool refused=false;try{session->execute_structure_transition(request,rfc::runtime::ExecutionLease(policy,stop));}
        catch(const std::exception &){refused=true;}
        check(refused,"native transition resource/cancellation ignored");
    }
    rfc::runtime::ExecutionHostPolicy admitted;admitted.name="native-transition-admitted";
    admitted.limits={maximum,maximum,maximum,maximum,maximum};
    const auto allowed=session->execute_structure_transition(request,rfc::runtime::ExecutionLease(admitted));
    check(allowed.at("result").has("green_response_defect"),"admitted structure lease bypassed native publication");
    session.reset();study.reset();for(const auto &result:retained)audit(result.at("result").at("fields"),result.at("result").at("fields").at("view"));
    std::cout<<"{\"status\":\"PASS\",\"requests\":"<<requests.names().size()<<",\"checks\":"<<checks<<"}\n";return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
