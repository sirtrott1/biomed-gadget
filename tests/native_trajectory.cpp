#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <rfc/storage/rfcx.hpp>
#include <rfc/algebra/materialize.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <atomic>
#include <iostream>
#include <limits>
using namespace gadget;
using namespace rfc::storage;
namespace {
std::size_t checks{};
void check(bool value,const char *message){++checks;if(!value)throw std::runtime_error(message);}
const TemporalTensor &tensor(const NativeRecord &owner,const NativeRecord &reference) {
    const auto id=reference.at("tensor").integer();
    for(const auto &t:owner.endpoint()->tensors)if(t.key==id)return t;
    throw std::runtime_error("trajectory tensor absent");
}
Vector vector(const NativeRecord &owner,const NativeRecord &reference) {
    const auto &t=tensor(owner,reference).tensor;
    check(t.support().count()==t.domain().logical_size(),"trajectory complete support lost");
    return t.values().dense();
}
rfc::numeric::ExactMatrix matrix(const NativeRecord &owner,const NativeRecord &reference) {
    const auto &t=tensor(owner,reference);check(t.tensor.support().count()==t.tensor.domain().logical_size(),"trajectory map supported zero lost");return temporal_tensor_matrix(t);
}
void audit(const NativeRecord &owner,const NativeRecord &value) {
    if(value.kind()==NativeRecordKind::Object) {
        if(value.has("schema")&&value.at("schema").kind()==NativeRecordKind::Text&&value.at("schema").text()=="gadget.scientific.tensor.view.v1") {
            const auto &t=tensor(owner,value);check(t.input||t.output,"trajectory unbound numerical field");return;
        }
        for(const auto &name:value.names())audit(owner,value.at(name));
    } else if(value.kind()==NativeRecordKind::Array)for(std::size_t i=0;i<value.size();++i)audit(owner,value.at(i));
}
void verify(const NativeRecord &operation,const ScientificStudy &study) {
    const auto result=operation.at("result"),owner=result.at("fields"),view=owner.at("view"),bindings=owner.at("bindings");
    check(NativeRecord::decode(result.at("source_study_publication").bytes()).encode()==study.native_publication().encode(),"trajectory original study bytes changed");
    check(encode_rfcx(owner.endpoint()->source).bytes==encode_rfcx(study.endpoint()->source).bytes,"trajectory original complex changed");
    audit(owner,view);audit(owner,bindings);
    std::vector<NativeRecord> intervals;std::size_t derived{},partial{};
    for(std::size_t i=0;i<bindings.size();++i) {
        const auto b=bindings.at(i);if(!b.has("kind"))continue;const auto kind=b.at("kind").text();
        if(kind=="declared_comparison_source") {
            const auto pins=b.at("fields");for(std::size_t j=0;j<pins.size();++j) {
                const auto pin=pins.at(j);const auto actual=study.field_accession(pin.at("field").text());
                check(pin.at("source_namespace").text()==std::to_string(actual.source_namespace())&&pin.at("source_cutoff").integer()==actual.sequence(),"trajectory original field producer pin changed");
                check(temporal_endpoint_detail::tensor_equal(tensor(owner,pin.at("source_field")).tensor,*actual.tensor()),"trajectory original field coefficients/support changed");
            }
        } else if(b.has("publication")) {
            const auto source=NativeRecord::decode(b.at("publication").bytes());++derived;
            check(std::to_string(source.source_namespace())==b.at("execution").at("source_namespace").text(),"trajectory derived source identity changed");
            const auto grades=b.at("grades");
            for(std::size_t k=0;k<source.endpoint()->source.grade_count();++k) {
                const auto n=source.endpoint()->source.grade({static_cast<std::uint16_t>(k)}).cell_count();
                check(grades.at(k).at("source_coordinates").size()==n,"trajectory original grade coordinates omitted");
                check(vector(owner,grades.at(k).at("identity_metric_diagonal"))==Vector(n,Rational{1}),"trajectory declared identity metric changed");
                if(grades.at(k).has("ambient_lift")) {
                    const auto lift=matrix(owner,grades.at(k).at("ambient_lift"));check(lift.cols()==n,"trajectory subcomplex map shape changed");
                    for(std::size_t j=0;j<n;++j)check(lift.get(grades.at(k).at("ambient_coordinates").at(j).integer().get_ui(),j)==1,"trajectory original subcomplex coordinate changed");
                }
            }
            if(kind=="trajectory_interval")intervals.push_back(source);
        } else if(kind=="trajectory_endpoint_reported_support") {
            const auto reference=b.at("reported_values");const auto &t=tensor(owner,reference);if(t.tensor.support().count()!=t.tensor.domain().logical_size())++partial;
        }
    }
    check(derived>=1+2*view.at("intervals").size(),"trajectory derived source family incomplete");
    check(intervals.size()==view.at("intervals").size(),"trajectory interval source count changed");
    for(std::size_t i=0;i<intervals.size();++i) {
        const auto &source=intervals[i].endpoint()->source;const auto interval=view.at("intervals").at(i);
        check(matrix(owner,interval.at("B1"))==rfc::field::materialize_boundary_matrix(source,{1}),"trajectory derived boundary differs from retained occurrences");
        const auto fields=interval.at("fields");for(std::size_t j=0;j<fields.size();++j) {
            const auto f=fields.at(j);if(f.has("complete_support")&&!f.at("complete_support").boolean()){check(!f.has("harmonic"),"trajectory missing endpoint invented harmonic values");continue;}
            const auto x=vector(owner,f.at("input")),g=vector(owner,f.at("gradient")),c=vector(owner,f.at("curl")),h=vector(owner,f.at("harmonic"));
            check(rfc::numeric::add(rfc::numeric::add(g,c),h)==x,"trajectory native full sector reconstruction failed");
            check(rfc::numeric::dot(g,c)==0&&rfc::numeric::dot(g,h)==0&&rfc::numeric::dot(c,h)==0,"trajectory native sectors not orthogonal");
        }
    }
    (void)partial;
}
}
int main(int argc,char **argv)try {
    check(argc==4,"trajectory study, request list and fresh output required");auto study=std::make_shared<const ScientificStudy>(argv[1]);
    const auto requests=admit_record_view(Json::parse(read_text(argv[2])));auto session=std::make_unique<ScientificSession>(study);fs::create_directories(argv[3]);std::vector<NativeRecord> retained;
    for(const auto &name:requests.names()) {
        const auto result=session->execute_record(requests.at(name));verify(result,*study);const auto path=fs::path(argv[3])/(name+".rfcxs");result.write(path);
        const auto reopened=NativeRecord::read(path),replay=replay_scientific_record(reopened);verify(reopened,*study);verify(replay,*study);
        const auto a=reopened.at("result").at("fields"),b=replay.at("result").at("fields");check(a.endpoint()->tensors.size()==b.endpoint()->tensors.size(),"trajectory replay family changed");
        for(std::size_t i=0;i<a.endpoint()->tensors.size();++i)check(temporal_endpoint_detail::tensor_equal(a.endpoint()->tensors[i],b.endpoint()->tensors[i]),"trajectory replay tensor/support/carrier changed");
        retained.push_back(reopened);
    }
    const auto change_member=[](NativeRecordValue &v,const std::string &name)->NativeRecordValue & {
        for(auto &[key,item]:std::get<NativeRecordValue::Object>(v.data).members)if(key==name)return item;
        throw std::runtime_error("trajectory negative control member missing");
    };
    for(const auto &change:{"source_namespace","source_cutoff"}) {
        const auto base=retained.front(),fields=base.at("result").at("fields");auto changed=fields.unbound_value();
        auto &bindings=std::get<NativeRecordValue::Array>(change_member(changed,"bindings").data).elements;
        bool altered=false;
        for(std::size_t index=0;index<bindings.size();++index) {
            const auto original_binding=fields.at("bindings").at(index);
            if(!original_binding.has("schema") || original_binding.at("schema").text()!="gadget.trajectory.source.binding.v1")continue;
            auto &binding=bindings[index];
            if(std::string(change)=="source_namespace")change_member(change_member(binding,"execution"),change)=NativeRecordValue(std::string("foreign-derived-owner"));
            else change_member(binding,change)=NativeRecordValue(rfc::numeric::Integer(999));
            altered=true;break;
        }
        check(altered,"trajectory source negative control absent");
        const auto bad_fields=fields.derive(std::move(changed));std::uint64_t key=1;for(const auto &source:base.raw_sources())key=std::max(key,rfc::machine::checked_size_add(source.key(),1));
        auto archive=bad_fields.archive_source(key);auto value=base.unbound_value();change_member(change_member(value,"result"),"fields")=NativeRecordValue(NativeRecordValue::Object{
            {{"schema",NativeRecordValue(std::string("gadget.scientific.publication.view.v1"))},
             {"publication",NativeRecordValue(NativeRecordValue::ByteSpan{key,0,archive.bytes().size(),false})},
             {"path",NativeRecordValue(NativeRecordValue::Array{})}}});
        const auto bad=base.derive(std::move(value),{std::move(archive)});bool refused=false;try{(void)replay_scientific_record(bad);}catch(const std::exception &){refused=true;}check(refused,"trajectory changed source provenance replayed");
    }
    const auto request=render_record_view(requests.at(requests.names().front()));const auto maximum=std::numeric_limits<std::size_t>::max();
    for(std::size_t scenario=0;scenario<7;++scenario) {
        rfc::runtime::ExecutionHostPolicy p;p.limits={maximum,maximum,maximum,maximum,maximum};auto stop=std::make_shared<std::atomic_bool>(scenario==0);
        if(scenario==1)p.limits.records=0;
        if(scenario==2)p.limits.rational_slots=0;
        if(scenario==3)p.limits.sparse_entries=0;
        if(scenario==4)p.limits.scalar_work=0;
        if(scenario==5)p.limits.integer_bits=0;
        if(scenario==6)p.allow_allocation=false;
        bool refused=false;try{session->execute_trajectory(request,rfc::runtime::ExecutionLease(p,stop));}catch(const std::exception &){refused=true;}check(refused,"trajectory current resource/cancellation ignored");
    }
    rfc::runtime::ExecutionHostPolicy p;p.limits={maximum,maximum,maximum,maximum,maximum};check(session->execute_trajectory(request,rfc::runtime::ExecutionLease(p)).at("result").has("intervals"),"trajectory explicit admitted lease failed");
    const auto hardware=rfc::runtime::detect_hardware();
    rfc::runtime::CpuResourcePolicy one_worker;one_worker.worker_limit=1;
    one_worker.declared_buffer_byte_limit=maximum;
    auto resources=std::make_shared<CpuResources>(hardware,one_worker,[hardware]{return hardware;});
    ScientificSession single(study,resources);
    check(single.execute_trajectory(request,rfc::runtime::ExecutionLease(p)).at("result").has("intervals"),"trajectory caller reservation charged twice");
    session.reset();study.reset();for(const auto &r:retained){const auto owner=r.at("result").at("fields");audit(owner,owner.at("view"));}
    std::cout<<"{\"status\":\"PASS\",\"checks\":"<<checks<<",\"requests\":"<<requests.size()<<"}\n";
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
