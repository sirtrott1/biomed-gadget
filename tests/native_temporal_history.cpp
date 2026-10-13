#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/rfcx.hpp>
#include <atomic>
#include <iostream>
#include <set>
using namespace gadget;
using Matrix=rfc::numeric::ExactMatrix;
namespace {
std::size_t checks{};
void check(bool value,const char *message) {++checks;if(!value)throw std::runtime_error(message);}
const rfc::storage::TemporalTensor &tensor(const NativeRecord &owner,const NativeRecord &reference) {
    const auto key=reference.at("tensor").integer();
    for(const auto &field:owner.endpoint()->tensors)if(field.key==key)return field;
    throw std::runtime_error("history field absent");
}
Vector values(const NativeRecord &owner,const NativeRecord &reference) {
    const auto &field=tensor(owner,reference).tensor;
    check(field.support().count()==field.domain().logical_size(),"history complete support lost");
    Vector result;for(std::size_t i=0;i<field.domain().logical_size();++i)result.push_back(*field.state(i));return result;
}
void audit(const NativeRecord &owner,const NativeRecord &value) {
    if(value.kind()==NativeRecordKind::Object) {
        if(value.has("schema")&&value.at("schema").kind()==NativeRecordKind::Text&&value.at("schema").text()=="gadget.scientific.tensor.view.v1") {
            const auto &field=tensor(owner,value);check(field.output.has_value(),"history axis absent");
            if(!value.has("coordinates"))check(field.tensor.support().count()==field.tensor.domain().logical_size(),"history supported zero omitted");
            return;
        }
        for(const auto &name:value.names())audit(owner,value.at(name));
    } else if(value.kind()==NativeRecordKind::Array)for(std::size_t i=0;i<value.size();++i)audit(owner,value.at(i));
}
void verify(const NativeRecord &operation,const ScientificStudy &study) {
    const auto fields=operation.at("result").at("fields"),view=fields.at("view"),binding=fields.at("bindings").at(0);
    audit(fields,view);audit(fields,fields.at("bindings"));
    check(NativeRecord::decode(operation.at("result").at("source_study_publication").bytes()).encode()==study.native_publication().encode(),"original coordinator changed");
    const auto endpoints=binding.at("endpoints");check(endpoints.size()==view.at("endpoints").size(),"history endpoint family lost");
    std::set<std::uint64_t> axes;
    for(std::size_t t=0;t<endpoints.size();++t) {
        const auto endpoint=endpoints.at(t),decl=view.at("endpoints").at(t).at("declaration");
        const auto original=NativeRecord::decode(endpoint.at("original_study_publication").at("publication").bytes());
        const auto context=NativeRecord::decode(endpoint.at("context_source").at("publication").bytes());
        const ScientificStudy original_study(original);
        for(std::size_t k=0;k<endpoint.at("grades").size();++k) {
            const auto grade=endpoint.at("grades").at(k);const auto axis=grade.at("coordinate_carrier").integer().get_ui();
            check(axes.insert(axis).second,"equal dimension or grade identified independent endpoint axes");
            const auto &arena=context.endpoint()->source.grade(rfc::core::Grade{static_cast<std::uint16_t>(k)});
            check(tensor(fields,grade.at("original_coordinates")).tensor.domain().logical_size()==arena.cell_count(),"history source lift extent changed");
        }
        auto inputs=endpoint.at("diagonal_inputs");
        std::vector<NativeRecord> retained{endpoint.at("input")};for(std::size_t i=0;i<inputs.size();++i)retained.push_back(inputs.at(i));
        for(const auto &input:retained) {
            const auto pin=original_study.field_accession(input.at("id").text());
            check(input.at("source_namespace").text()==std::to_string(pin.source_namespace())&&input.at("source_cutoff").integer()==pin.sequence(),"history producer cutoff changed");
            const auto descriptor=pin.descriptor();
            check(input.at("source_tensor").integer()==descriptor.request().subject,"history original tensor changed");
            const auto lift=rfc::storage::temporal_tensor_matrix(tensor(fields,input.at("realization")));
            check(lift.cols()==pin.tensor()->domain().logical_size()&&lift.rows()==values(fields,input.at("execution_values")).size(),"history realization carrier changed");
            const auto &reported=tensor(fields,input.at("reported_values")).tensor;
            for(std::size_t i=0;i<lift.rows();++i) {
                std::optional<rfc::numeric::Rational> original;
                for(const auto &[position,coefficient]:lift.entries())if(position.first==i&&coefficient==1)original=pin.tensor()->state(position.second);
                check(reported.state(i)==original,"reported support or original value changed under completion");
            }
        }
        const auto observed=decl.at("observations");const auto x=values(fields,view.at("endpoints").at(t).at("field"));
        check(observed.size()==x.size(),"history original observation population changed");
        for(std::size_t i=0;i<x.size();++i) {
            check(values(fields,observed.at(i).at("execution_value"))[0]==x[i],"history execution observation changed");
            if(observed.at(i).at("reported").boolean())check(values(fields,observed.at(i).at("observed_value"))[0]==x[i],"reported observation changed");
            else check(observed.at(i).at("observed_value").kind()==NativeRecordKind::Null,"missing observation became supported zero");
        }
    }
    for(std::size_t t=0;t<view.at("steps").size();++t) {
        const auto step=view.at("steps").at(t);Vector sum(values(fields,step.at("field_delta")).size());
        for(std::size_t i=0;i<step.at("field_delta_channels").size();++i)sum=rfc::numeric::add(sum,values(fields,step.at("field_delta_channels").at(i).at("field")));
        check(sum==values(fields,step.at("field_delta")),"ordered history channels fail reconstruction");
    }
    Vector terminal(values(fields,view.at("endpoints").at(endpoints.size()-1).at("field")).size());
    for(std::size_t i=0;i<view.at("terminal_channels").size();++i)terminal=rfc::numeric::add(terminal,values(fields,view.at("terminal_channels").at(i).at("field")));
    check(terminal==values(fields,view.at("endpoints").at(endpoints.size()-1).at("field")),"terminal history reconstruction lost");
}
}
int main(int argc,char **argv) try {
    check(argc==4,"history study, requests and owned output required");
    auto initial=std::make_shared<const ScientificStudy>(argv[1]);
    auto study=std::make_shared<const ScientificStudy>(NativeRecord::annotate(initial->native_publication().unbound_value(),initial->native_publication()));initial.reset();
    const auto requests=admit_record_view(Json::parse(read_text(argv[2])));fs::create_directories(argv[3]);
    auto session=std::make_unique<ScientificSession>(study);std::vector<NativeRecord> retained;
    for(const auto &name:requests.names()) {
        std::cerr<<"native history "<<name<<" begin\n"<<std::flush;
        const auto result=session->execute_record(requests.at(name));verify(result,*study);const auto path=fs::path(argv[3])/(name+".rfcxs");result.write(path);
        const auto reopened=NativeRecord::read(path);verify(reopened,*study);const auto replay=replay_scientific_record(reopened);verify(replay,*study);
        const auto a=reopened.at("result").at("fields"),b=replay.at("result").at("fields");check(a.endpoint()->tensors.size()==b.endpoint()->tensors.size(),"replay history field count differs");
        for(std::size_t i=0;i<a.endpoint()->tensors.size();++i)check(rfc::storage::temporal_endpoint_detail::tensor_equal(a.endpoint()->tensors[i],b.endpoint()->tensors[i]),"replay history complete tensor differs");
        retained.push_back(reopened);std::cerr<<"native history "<<name<<" replay PASS\n"<<std::flush;
    }
    const auto request=render_record_view(requests.at(requests.names().front()));const auto maximum=std::numeric_limits<std::size_t>::max();
    for(std::size_t scenario=0;scenario<7;++scenario) {
        rfc::runtime::ExecutionHostPolicy policy;policy.limits={maximum,maximum,maximum,maximum,maximum};auto stop=std::make_shared<std::atomic_bool>(scenario==0);
        if(scenario==1)policy.limits.records=0;
        if(scenario==2)policy.limits.rational_slots=0;
        if(scenario==3)policy.limits.sparse_entries=0;
        if(scenario==4)policy.limits.scalar_work=0;
        if(scenario==5)policy.limits.integer_bits=0;
        if(scenario==6)policy.allow_allocation=false;
        bool refused=false;try{session->execute_temporal_history(request,rfc::runtime::ExecutionLease(policy,stop));}catch(const std::exception &){refused=true;}check(refused,"history resource/cancellation ignored");
    }
    rfc::runtime::ExecutionHostPolicy policy;policy.limits={maximum,maximum,maximum,maximum,maximum};
    check(session->execute_temporal_history(request,rfc::runtime::ExecutionLease(policy)).at("result").at("identity_verified").boolean(),"admitted history lease failed");
    session.reset();study.reset();for(const auto &result:retained) {const auto fields=result.at("result").at("fields");audit(fields,fields.at("view"));}
    std::cout<<"{\"status\":\"PASS\",\"requests\":"<<requests.names().size()<<",\"checks\":"<<checks<<"}\n";return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
