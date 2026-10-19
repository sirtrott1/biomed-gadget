#include <gadget/data_study.hpp>
#include <gadget/record_view.hpp>
#include "archive_identity.hpp"
#include "scientific_detail.hpp"
#include "source_grade_fields.hpp"
#include <rfc/core/builder.hpp>

namespace gadget {
NativeRecord bind_measurement_study(const NativeRecord &data, const NativeRecord &request,
                                   std::stop_token stop) {
    using namespace scientific_detail;
    using namespace rfc::storage;
    using N = NativeRecordValue;
    const auto cancelled = [&] {if (stop.stop_requested()) throw std::runtime_error("measurement study binding cancelled");};
    cancelled();
    const auto admission=request.at("max_bytes").integer();
    require(admission>0 && admission<=rfc::numeric::Integer(std::to_string(source_byte_limit())),
            "measurement study byte admission");
    const auto budget=static_cast<std::size_t>(admission.get_ui());
    (void)archive_identity(request.owning_publication(),budget,stop);
    const auto q = render_record_view(request);
    members(q,{"schema","data_sha256","id","design","responses","covariates","max_bytes"});
    require(q.at("schema").text()=="gadget.data.study.v1" && data.is_root() &&
            data.at("schema").text()=="gadget.native-data.field.v1", "native measurement study source/schema");
    require(budget>0 && archive_identity(data,budget,stop)==q.at("data_sha256").text(),"measurement study source archive differs");
    const auto endpoint=data.endpoint();
    const auto found=std::find_if(endpoint->tensors.begin(),endpoint->tensors.end(),[](const auto &field){return field.key==72;});
    require(found!=endpoint->tensors.end() && found->output==71 && !found->input,
            "native measurement tensor missing");
    const auto rows=data.at("shape").at(0).integer(), columns=data.at("shape").at(1).integer();
    require(rows>0 && columns>0 && found->tensor.domain().rank()==2 &&
            rows==rfc::numeric::Integer(std::to_string(found->tensor.domain().extents()[0])) &&
            columns==rfc::numeric::Integer(std::to_string(found->tensor.domain().extents()[1])),"native measurement shape differs");
    const auto n=found->tensor.domain().extents()[1], features=found->tensor.domain().extents()[0];
    const auto responses=q.at("responses"),covariates=q.at("covariates");
    source_array(responses,"response coordinate declarations");source_array(covariates,"covariate declarations");
    require(responses.size()>0,"measurement study requires explicit response coordinates");
    const auto slots=rfc::machine::checked_size_mul(n,rfc::machine::checked_size_add(responses.size(),covariates.size()));
    require(rfc::machine::checked_size_mul(slots,512)<=budget,"measurement chart working byte admission");
    auto original=data.archive_source(1);
    const auto original_size=original.bytes().size();require(original_size<=budget,"original measurement archive byte admission");
    std::size_t coefficient_bytes{};
    std::map<std::size_t,std::size_t> selected_rows;
    for(std::size_t r=0;r<responses.size();++r)
        ++selected_rows[count(responses.at(r).at("row"),features-1,"response row")];
    for(const auto &entry:found->tensor.values().entries()) {
        cancelled();const auto row=selected_rows.find(entry.index/n);
        if(row==selected_rows.end())continue;
        const auto bits=rfc::machine::checked_size_add(mpz_sizeinbase(entry.value.get_num_mpz_t(),2),
            mpz_sizeinbase(entry.value.get_den_mpz_t(),2));
        coefficient_bytes=rfc::machine::checked_size_add(coefficient_bytes,
            rfc::machine::checked_size_mul(bits/8+2,row->second));
    }
    std::size_t covariate_bytes{};
    for(std::size_t c=0;c<covariates.size();++c) {
        const auto values=request.at("covariates").at(c).at("values");
        require(values.kind()==NativeRecordKind::Array && values.size()==n,
                "covariate must declare every original observation column");
        for(std::size_t i=0;i<n;++i) {
            cancelled();const auto value=values.at(i);
            if(value.kind()==NativeRecordKind::Null)continue;
            std::size_t size{};
            if(value.kind()==NativeRecordKind::Text)size=value.bytes().size();
            else {
                const auto coefficient=value.rational();
                size=rfc::machine::checked_size_add(mpz_sizeinbase(coefficient.get_num_mpz_t(),2),
                    mpz_sizeinbase(coefficient.get_den_mpz_t(),2))/8+2;
            }
            covariate_bytes=rfc::machine::checked_size_add(covariate_bytes,size);
        }
    }
    coefficient_bytes=rfc::machine::checked_size_add(coefficient_bytes,covariate_bytes);
    const auto working_bytes=rfc::machine::checked_size_add(rfc::machine::checked_size_mul(slots,512),
        rfc::machine::checked_size_add(rfc::machine::checked_size_mul(coefficient_bytes,16),
                                      rfc::machine::checked_size_mul(original_size,2)));
    require(working_bytes<=budget,"exact measurement coefficient working admission");
    rfc::core::RelationalComplexBuilder builder(rfc::core::RelationSignature{});
    for(std::size_t i=0;i<n;++i) builder.declare_witness_identity(i+1);
    TemporalEndpoint chart(builder.freeze());
    chart.carriers.push_back({1,TemporalCarrierKind::SourceGrade,{0},{},0,rfc::numeric::IndexDomain::linear(n)});
    assay_detail::SourceGradeFields fields(chart);
    N::Array meanings,coordinates;
    const auto source_view = [&](N::Array path) {
        return N(N::Object{{{"schema",N(std::string("gadget.scientific.publication.view.v1"))},
            {"publication",N(N::ByteSpan{1,0,original_size,false})},{"path",N(std::move(path))}}});
    };
    const auto evidence=Json::parse("[{\"document\":\"data\",\"start\":0,\"end\":"+std::to_string(original_size)+"}]");
    std::set<std::string> ids;
    const auto declaration=[&](const std::string &id, const N &key,const std::string &quantity,
                               const std::string &units,const std::string &value_type) {
        require(!id.empty() && ids.insert(id).second,"duplicate or empty measurement chart field");
        meanings.elements.emplace_back(N::Object{{{"id",N(id)},{"key",N(std::get<rfc::numeric::Integer>(key.data).get_str())},{"context",N(std::string("all"))},
            {"channel",N(std::string("observation"))},{"quantity",N(quantity)},{"units",N(units)},
            {"method",N(std::string("Explicit native measurement chart"))},{"value_type",N(value_type)},
            {"evidence",import_record_view(evidence)}}});
    };
    for(std::size_t r=0;r<responses.size();++r) {
        cancelled();const auto response=responses.at(r);members(response,{"id","row"});
        const auto row=count(response.at("row"),features-1,"response row");const auto id=text(response.at("id"));
        Vector values(n);std::vector<std::size_t> support;
        const auto base=row*n;
        for(std::size_t i=0;i<n;++i)
            if(found->tensor.support().contains(base+i))support.push_back(i);
        const auto &source_entries=found->tensor.values().entries();
        auto first=std::lower_bound(source_entries.begin(),source_entries.end(),base,
            [](const auto &entry,std::size_t index){return entry.index<index;});
        for(;first!=source_entries.end() && first->index-base<n;++first)
            values[first->index-base]=first->value;
        TemporalTensor field;field.role=TemporalTensorRole::Field;field.output=1;
        field.interpretation="Declared source row on the original ordered observation columns";
        std::vector<rfc::numeric::SparseRationalEntry> nonzeros;
        for(std::size_t i=0;i<n;++i)if(values[i]!=0)nonzeros.push_back({i,values[i]});
        field.tensor={rfc::numeric::IndexDomain::linear(n),rfc::numeric::StructuralSupport::sparse(n,support),
            rfc::numeric::SparseRationalVector(n,std::move(nonzeros))};
        declaration(id,fields.insert(std::move(field)),data.at("quantity").text(),data.at("units").text(),"quantitative");
        coordinates.elements.emplace_back(N::Object{{{"field",N(id)},{"source_row",N(data.at("rows").at(row).integer())},
            {"source_label",source_view(N::Array{{N(std::string("row_labels")),N(rfc::numeric::Integer(std::to_string(row)))}})},{"selected_row",N(rfc::numeric::Integer(std::to_string(row)))}}});
    }
    const auto response_fields=meanings.elements.size();
    for(std::size_t c=0;c<covariates.size();++c) {
        cancelled();const auto column=covariates.at(c);members(column,{"id","values","quantity","units","value_type"});
        const auto values=column.at("values");require(values.size()==n,"covariate must declare every original observation column");
        const auto type=column.at("value_type").text();require(type=="quantitative"||type=="indicator","categorical values require an explicit indicator declaration");
        Vector amplitudes(n);std::vector<std::size_t> support;
        for(std::size_t i=0;i<n;++i) {
            const auto original=request.at("covariates").at(c).at("values").at(i);
            if(original.kind()==NativeRecordKind::Null)continue;
            const auto value=original.kind()==NativeRecordKind::Text?exact_rational(original.text()):original.rational();
            require(type!="indicator"||value==0||value==1,"indicator values must be zero or one");amplitudes[i]=value;support.push_back(i);
        }
        TemporalTensor field;field.role=TemporalTensorRole::Field;field.output=1;
        field.interpretation="Authored covariate values on the explicit original observation order";
        std::vector<rfc::numeric::SparseRationalEntry> nonzeros;
        for(std::size_t i=0;i<n;++i)if(amplitudes[i]!=0)nonzeros.push_back({i,amplitudes[i]});
        field.tensor={rfc::numeric::IndexDomain::linear(n),rfc::numeric::StructuralSupport::sparse(n,support),
            rfc::numeric::SparseRationalVector(n,std::move(nonzeros))};
        declaration(text(column.at("id")),fields.insert(std::move(field)),text(column.at("quantity")),text(column.at("units")),type);
    }
    cancelled();
    const std::string binding="Native measurement coordinates and authored observation chart. ";
    const auto note=RawByteSource(2,std::as_bytes(std::span(binding.data(),binding.size())));
    N result(N::Object{{{"schema",N(std::string("gadget.data.study.binding.v1"))},
        {"data_sha256",N(q.at("data_sha256").text())},{"source",N(N::ByteSpan{1,0,original_size,false})},
        {"response_coordinates",N(std::move(coordinates))},{"observation_coordinates",source_view(N::Array{{N(std::string("columns"))}})},
        {"observation_labels",source_view(N::Array{{N(std::string("column_labels"))}})}}});
    const auto producer=retain_operation(std::move(chart),new_native_record_namespace(),{std::move(original),note},request,std::move(result));
    N::Array docs{{N(N::Object{{{"id",N(std::string("binding"))},{"path",N(std::string("binding"))},
        {"sha256",N("sha256:"+sha256(binding))},{"format",N(std::string("utf8"))}}})}};
    docs.elements.emplace_back(N::Object{{{"id",N(std::string("data"))},{"path",N(std::string("data"))},
        {"sha256",N(q.at("data_sha256").text())},{"format",N(std::string("binary"))}}});
    N::Array raw_docs{{N(N::Object{{{"id",N(std::string("binding"))},{"bytes",N(N::ByteSpan{2,0,binding.size(),false})}}}),
        N(N::Object{{{"id",N(std::string("data"))},{"bytes",N(N::ByteSpan{1,0,original_size,false})}}})}};
    const auto protocol=producer.at("request_publication");const auto protocol_bytes=protocol.bytes();
    docs.elements.emplace_back(N::Object{{{"id",N(std::string("protocol"))},{"path",N(std::string("protocol"))},
        {"sha256",N("sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(protocol_bytes.data()),protocol_bytes.size())))},
        {"format",N(std::string("binary"))}}});
    raw_docs.elements.emplace_back(N::Object{{{"id",N(std::string("protocol"))},{"bytes",protocol.materialize_value()}}});
    auto covariate_evidence=Json::parse("[{\"document\":\"protocol\",\"start\":0,\"end\":"+std::to_string(protocol_bytes.size())+"}]");
    for(std::size_t c=response_fields;c<meanings.elements.size();++c)
        for(auto &[name,value]:std::get<N::Object>(meanings.elements[c].data).members)
            if(name=="evidence")value=import_record_view(covariate_evidence);
    auto design=import_record_view(q.at("design"));
    const auto control=q.at("design");members(control,{"question","population","selection","limitations"});
    std::get<N::Object>(design.data).members.emplace_back("covariates",N(N::Array{}));
    N study(N::Object{{{"schema",N(std::string("gadget.scientific.study.v3"))},{"id",N(text(q.at("id")))},
        {"design",std::move(design)},{"documents",N(std::move(docs))},
        {"cell_bindings",import_record_view(Json::parse(R"([{"grade":0,"id_prefix":"observation:","label_prefix":"Observation ","kind":"source_observation","evidence":{"kind":"documents","documents":["binding"]}}])"))},
        {"occurrence_evidence",N(std::string("parent_cell"))},
        {"contexts",import_record_view(Json::parse(R"([{"id":"all","block":null,"metrics":[{"kind":"identity"}],"time":null}])"))},
        {"fields",N(std::move(meanings))},{"field_sets",N(N::Array{})}}});
    cancelled();return ScientificStudy::publish_annotations(producer,std::move(study),N(std::move(raw_docs)));
}
}
