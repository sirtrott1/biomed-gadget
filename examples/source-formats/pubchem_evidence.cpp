#include "pubchem_evidence_detail.hpp"
#include "pubchem_evidence.hpp"
#include <gadget/record_view.hpp>
#include <rfc/database/interpretation_codec.hpp>
#include <rfc/database/native_query.hpp>
#include <set>
#include <algorithm>

namespace gadget::examples {
namespace {
using N = NativeRecordValue;
using namespace rfc::database;
void require(bool value,const char *message) {if(!value)throw std::invalid_argument(message);}
void cancelled(std::stop_token stop) {if(stop.stop_requested())throw std::runtime_error("PubChem native admission cancelled");}
std::string bytes_text(const NativeRecord &record) {
    require(record.kind()==NativeRecordKind::Bytes,"PubChem sources must be opaque original bytes");
    const auto b=record.bytes();return {reinterpret_cast<const char *>(b.data()),b.size()};
}
struct Input {
    std::string assay;
    std::optional<std::string> concise;
};
std::size_t source_cap(const Json &request) {
    const auto limit=request.at("limits").at("source_bytes").integer();
    require(limit>=0 && static_cast<std::uint64_t>(limit)<=std::numeric_limits<std::size_t>::max(),"PubChem source byte budget is not addressable");
    return static_cast<std::size_t>(limit);
}
std::size_t working_cap(const Json &request) {
    const auto limit=request.at("limits").at("native_bytes").integer();
    require(limit>0 && static_cast<std::uint64_t>(limit)<=source_byte_limit(),
            "PubChem native working byte admission");
    return static_cast<std::size_t>(limit);
}
Input source_bytes(const NativeRecord &source,std::size_t limit,std::size_t working_bytes) {
    require(source.is_root() && source.at("schema").text()=="gadget.pubchem.sources.native.v1","unsupported PubChem original source owner");
    for(const auto &name:source.names())require(name=="schema"||name=="assay"||name=="concise","undeclared PubChem byte source");
    (void)source.at("assay").byte_source();
    const auto assay_bytes=source.at("assay").bytes().size();
    const auto concise_bytes=source.has("concise")?source.at("concise").bytes().size():0;
    require(assay_bytes<=working_bytes/32 && concise_bytes<=working_bytes/32-assay_bytes,
            "PubChem decoded source working admission");
    if(source.at("assay").bytes().size()>limit)throw std::length_error("PubChem original assay byte budget");
    Input input{bytes_text(source.at("assay")),{}};
    if(source.has("concise")) {(void)source.at("concise").byte_source();if(source.at("concise").bytes().size()>limit)throw std::length_error("PubChem original concise byte budget");input.concise=bytes_text(source.at("concise"));}
    return input;
}
std::vector<std::byte> archive(const InterpretedSnapshot &native) {return encode_interpreted_snapshot(native);}
std::string byte_digest(std::span<const std::byte> bytes) {
    return "sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size()));
}
Json publication_receipt(const PubChemEvidencePacket &packet,std::span<const std::byte> source,
                         std::span<const std::byte> request,std::span<const std::byte> encoded) {
    auto receipt=Json::parse(packet.receipt.dump_compact());const auto source_hash=byte_digest(source),request_hash=byte_digest(request);
    receipt.set("source_format_digest",receipt.at("source_digest"));
    receipt.set("source_publication_sha256",Json::string(source_hash));
    receipt.set("admission_request_sha256",Json::string(request_hash));
    receipt.set("source_digest",Json::string("sha256:"+sha256(packet.logical.revision+"\n"+source_hash+"\n"+request_hash+"\n"+byte_digest(encoded))));
    receipt.set("native_sha256",Json::string(byte_digest(encoded)));return receipt;
}
// Compile literal source metadata with the existing authored native value
// constructors. Generic application presentation shorthands (such as an object
// named "rational") do not reinterpret a provider's original JSON members.
N literal_source_value(const Json &value,std::size_t depth=1) {
    if(depth>native_data_limits().temporal.declarations.max_depth)throw std::length_error("PubChem native metadata nesting profile");
    switch(json_object_get_type(value.raw())) {
    case json_type_null:return {};
    case json_type_boolean:return N(value.boolean());
    case json_type_int:return N(rfc::numeric::Integer(value.dump_compact()));
    case json_type_double:{const auto token=value.dump_compact();return N(N::Decimal{exact_decimal(token),token});}
    case json_type_string:return N(value.text());
    case json_type_array:{N::Array out;out.elements.reserve(value.size());for(std::size_t i=0;i<value.size();++i)out.elements.push_back(literal_source_value(value.at(i),depth+1));return N(std::move(out));}
    case json_type_object:{N::Object out;auto current=json_object_iter_begin(value.raw()),end=json_object_iter_end(value.raw());while(!json_object_iter_equal(&current,&end)) {const std::string key=json_object_iter_peek_name(&current);out.members.emplace_back(key,literal_source_value(value.at(key),depth+1));json_object_iter_next(&current);}return N(std::move(out));}
    }
    throw std::invalid_argument("unsupported PubChem source JSON value");
}
N number(std::size_t value) {return N(rfc::numeric::Integer(std::to_string(value)));}
std::map<std::uint64_t,NativeRecord> records_by_key(const NativeRecord &publication,const InterpretedSnapshot &native) {
    std::map<std::string,NativeRecord> ids;const auto rows=publication.at("records");
    for(std::size_t i=0;i<rows.size();++i)require(ids.emplace(rows.at(i).at("id").text(),rows.at(i)).second,"duplicate typed PubChem observation identity");
    std::map<std::uint64_t,NativeRecord> result;
    for(const auto &binding:native.lineage().assertions)result.emplace(binding.native_key,ids.at(binding.origin.token));
    return result;
}
N selected(const NativeRecord &record,std::uint64_t key) {
    auto value=record.materialize_value();auto *owner=std::get_if<N::Owned>(&value.data);
    require(owner && owner->cell && owner->cell->grade==rfc::core::Grade{1} && owner->cell->key==key,
            "PubChem selected record lost its actual assertion-cell owner");
    auto payload=*owner->payload;auto &members=std::get<N::Object>(payload.data).members;
    members.emplace_back("native_assertion_key",N(std::to_string(key)));
    owner->payload=std::make_shared<const N>(std::move(payload));return value;
}
std::size_t cap(const Json &request) {return static_cast<std::size_t>(request.at("max_records").integer());}
} // namespace
bool is_pubchem_evidence(const NativeRecord &publication) {
    return publication.has("schema") && publication.at("schema").kind()==NativeRecordKind::Text &&
        publication.at("schema").text()=="gadget.evidence.pubchem.native.v1";
}
void verify_pubchem_evidence(const NativeRecord &publication,const InterpretedSnapshot &native) {
    require(is_pubchem_evidence(publication) && publication.is_root(),
            "unsupported native PubChem publication");
    const auto source=NativeRecord::decode(publication.at("source_input").bytes());
    const auto request=NativeRecord::decode(publication.at("admission_request").bytes());
    const auto control=render_record_view(request);const auto input=source_bytes(source,source_cap(control),working_cap(control));
    const auto packet=prepare_pubchem_evidence(control,input.assay,input.concise,publication.source_namespace());
    const auto expected=interpret_pubchem_evidence(packet);const auto expected_archive=archive(expected);
    const auto receipt=publication_receipt(packet,publication.at("source_input").bytes(),publication.at("admission_request").bytes(),expected_archive);
    require(render_record_view(publication.at("receipt")).dump()==receipt.dump() &&
            render_record_view(publication.at("snapshot")).dump()==packet.manifest.dump(),"PubChem receipt/source identity differs from retained admission");
    require(std::ranges::equal(publication.at("interpretation").bytes(),expected_archive) && archive(native)==expected_archive,
            "typed PubChem native interpretation differs from retained source");
    require(bytes_text(publication.at("assay_bytes"))==input.assay && publication.has("concise_bytes")==input.concise.has_value(),
            "PubChem retained assay/concise byte owners differ");
    if(input.concise)require(bytes_text(publication.at("concise_bytes"))==*input.concise,"PubChem concise bytes differ from retained source");
    const auto rows=publication.at("records");require(rows.size()==packet.records.size(),"PubChem native record count differs");
    auto normalized=Json::array();std::map<std::string,std::uint64_t> assertions,subjects;
    for(const auto &entry:native.lineage().assertions)assertions.emplace(entry.origin.token,entry.native_key);
    for(const auto &entry:native.lineage().subjects)subjects.emplace(entry.origin.token,entry.native_key);
    for(std::size_t i=0;i<rows.size();++i) {
        const auto row=rows.at(i);auto view=render_record_view(row);json_object_object_del(view.raw(),"source_bytes");normalized.push(view);
        const auto key=assertions.at(row.at("id").text());
        (void)selected(row,key);
        const auto owner=row.at("source_role").text()=="assay"?publication.at("assay_bytes"):publication.at("concise_bytes");
        require(row.at("source_bytes").byte_source().shares_owner(owner.byte_source()),"PubChem record byte binding belongs to another source owner");
    }
    require(normalized.dump_compact()==packet.records.dump_compact(),"PubChem native declarations differ from original readout occurrences");
    auto entities=Json::object();for(const auto &[id,value]:packet.entities)entities.set(id,value);
    require(render_record_view(publication.at("entities")).dump()==entities.dump(),"PubChem native subject declarations differ");
    for(const auto &[id,key]:subjects) {
        const auto value=publication.at("entities").at(id).materialize_value();const auto *owner=std::get_if<N::Owned>(&value.data);
        require(owner && owner->cell && owner->cell->grade==rfc::core::Grade{0} && owner->cell->key==key,
                "PubChem entity lost its actual subject-cell owner");
    }
    require(rfc::storage::encode_rfcx(publication.endpoint()->source).bytes==rfc::storage::encode_rfcx(native.endpoint()->source).bytes,
            "PubChem native complex differs from its typed source interpretation");
    const auto &declarations=*publication.endpoint()->declarations;
    for(const auto &definition:native.endpoint()->declarations->registry.definitions())
        require(declarations.registry.at(definition.id)==definition,"PubChem original payload type was altered");
    std::map<std::uint64_t,const rfc::storage::TemporalDeclaredValue *> values;
    for(const auto &value:declarations.values)values.emplace(value.key,&value);
    for(const auto &original:native.endpoint()->declarations->values) {
        require(values.contains(original.key),"PubChem original native payload is missing");const auto &value=*values.at(original.key);
        require(value.type==original.type && rfc::storage::temporal_declaration_detail::owner_equal(value,original) &&
                rfc::composite::declared_value_equal(value.value,original.value),"PubChem original typed payload was altered");
    }
}
NativeRecord query_pubchem_evidence(const NativeRecord &publication,const InterpretedSnapshot &native,const NativeRecord &request) {
    const auto query=render_record_view(request);const auto digest=publication.at("receipt").at("source_digest").text();
    validate_pubchem_query(query,digest);const auto op=query.at("operation").text();
    const std::map<std::string,std::uint64_t> kinds{{"descriptions",40},{"definitions",41},{"outcomes",42},{"readouts",43},{"revocations",44},{"links",45}};
    auto selected_keys=std::vector<std::uint64_t>{};
    if(op=="inspect")for(const auto &binding:native.lineage().assertions)selected_keys.push_back(binding.native_key);
    else {
        const auto kind=kinds.at(op);
        if(query.has("tid"))selected_keys=select_exact(native,kind,{{1},ExactPredicateKind::Equal,
            rfc::composite::DeclaredValue::integer(rfc::numeric::Integer(query.at("tid").text()))},native.assertions(kind).size());
        else selected_keys=native.assertions(kind);
    }
    if(query.has("sid")) {
        std::optional<std::uint64_t> subject;for(const auto &binding:native.lineage().subjects)
            if(binding.origin.token=="pubchem:sid:"+query.at("sid").text())subject=binding.native_key;
        std::erase_if(selected_keys,[&](auto assertion) {
            if(!subject)return true;
            for(const auto &occurrence:native.participants(assertion,11)) {
                if(occurrence.participant().key==*subject)return false;
            }
            return true;
        });
    }
    if(selected_keys.size()>cap(query))throw std::length_error("PubChem query result record limit");
    const auto rows=records_by_key(publication,native);N::Array records;
    for(const auto key:selected_keys)records.elements.push_back(selected(rows.at(key),key));
    return retain_operation(publication,request,N(N::Object{{
        {"schema",N(std::string("gadget.pubchem.query.result.v1"))},{"source_digest",N(digest)},
        {"operation",N(op)},{"selected",number(selected_keys.size())},{"records",N(std::move(records))},
        {"scope",N(std::string("Reported typed source assertions only. Missing values, outcomes, revocations and modifier readouts remain separate; no activity threshold, identifier equivalence or assay comparability is inferred."))}}}));
}
NativeRecord context_pubchem_evidence(const NativeRecord &publication,const InterpretedSnapshot &native,const NativeRecord &request) {
    const auto control=render_record_view(request);const auto digest=publication.at("receipt").at("source_digest").text();
    validate_pubchem_context(control,digest);const auto rows=records_by_key(publication,native);const auto keys=control.at("assertion_keys");N::Array records;
    for(std::size_t i=0;i<keys.size();++i) {
        const auto key=static_cast<std::uint64_t>(std::stoull(keys.at(i).text()));
        require(rows.contains(key),"PubChem context references an unknown native assertion");records.elements.push_back(selected(rows.at(key),key));
    }
    const auto question=control.at("question").text();
    return retain_operation(publication,request,N(N::Object{{
        {"schema",N(std::string("gadget.pubchem.context.v1"))},{"source_digest",N(digest)},
        {"question",N(question)},{"question_sha256",N("sha256:"+sha256(question))},
        {"records",N(std::move(records))},{"selected",number(keys.size())},
        {"availability",N(N::Object{{{"kind",N(std::string("frozen_source_capture"))},{"historical_availability_established",N(false)}}})},
        {"scope",N(std::string("Explicitly selected original native source assertions; interpretations remain unverified and source capture does not establish biological or historical observation time."))}}}));
}
} // namespace gadget::examples
namespace gadget::examples {
NativeRecord admit_pubchem_evidence(const NativeRecord &source,const NativeRecord &request,std::stop_token stop) {
    using namespace gadget::examples;using N=NativeRecordValue;
    require(request.is_root(),"PubChem admission requires an owning native request");
    const auto control=render_record_view(request);const auto input=source_bytes(source,source_cap(control),working_cap(control));cancelled(stop);
    const auto packet=prepare_pubchem_evidence(control,input.assay,input.concise,new_native_record_namespace(),stop);
    const auto native=interpret_pubchem_evidence(packet,stop);auto encoded=archive(native);
    std::vector<RawByteSource> raw;
    raw.push_back(source.archive_source(1));raw.push_back(request.archive_source(2));raw.push_back(RawByteSource::own(3,std::move(encoded)));
    const auto receipt=publication_receipt(packet,raw.at(0).bytes(),raw.at(1).bytes(),raw.at(2).bytes());
    raw.push_back(source.at("assay").byte_source().rekey(4));
    if(input.concise)raw.push_back(source.at("concise").byte_source().rekey(5));
    std::map<std::string,std::uint64_t> assertions,subjects;
    for(const auto &entry:native.lineage().assertions)assertions.emplace(entry.origin.token,entry.native_key);
    for(const auto &entry:native.lineage().subjects)subjects.emplace(entry.origin.token,entry.native_key);
    N::Array records;N::Object entities;
    for(std::size_t i=0;i<packet.records.size();++i) {
        cancelled(stop);const auto row=packet.records.at(i);auto value=literal_source_value(row);auto &fields=std::get<N::Object>(value.data).members;
        const bool assay=row.at("source_role").text()=="assay";
        fields.emplace_back("source_bytes",N(N::ByteSpan{assay?4u:5u,0,assay?input.assay.size():input.concise->size(),false}));
        records.elements.emplace_back(N::Owned{0,rfc::storage::StableCellKey{rfc::core::Grade{1},assertions.at(row.at("id").text())},{},std::make_shared<const N>(std::move(value))});
    }
    for(const auto &[id,value]:packet.entities)entities.members.emplace_back(id,N(N::Owned{0,rfc::storage::StableCellKey{rfc::core::Grade{0},subjects.at(id)},{},std::make_shared<const N>(literal_source_value(value))}));
    N::Object root{{{"schema",N(std::string("gadget.evidence.pubchem.native.v1"))},
        {"receipt",literal_source_value(receipt)},{"snapshot",literal_source_value(packet.manifest)},
        {"records",N(std::move(records))},{"entities",N(std::move(entities))},
        {"assay_bytes",N(N::ByteSpan{4,0,input.assay.size(),false})},
        {"source_input",N(N::ByteSpan{1,0,raw.at(0).bytes().size(),false})},
        {"admission_request",N(N::ByteSpan{2,0,raw.at(1).bytes().size(),false})},
        {"interpretation",N(N::ByteSpan{3,0,raw.at(2).bytes().size(),false})}}};
    if(input.concise)root.members.emplace_back("concise_bytes",N(N::ByteSpan{5,0,input.concise->size(),false}));
    cancelled(stop);auto publication=NativeRecord::create(N(std::move(root)),packet.schema.source_namespace,std::move(raw),native.endpoint().get());
    auto represented=Json::array();const auto admitted_rows=publication.at("records");
    for(std::size_t i=0;i<admitted_rows.size();++i) {
        cancelled(stop);auto view=render_record_view(admitted_rows.at(i));json_object_object_del(view.raw(),"source_bytes");represented.push(view);
    }
    require(represented.dump_compact()==packet.records.dump_compact(),"PubChem source metadata exceeds the literal native representation");
    if(publication.encode().size()>packet.limits.max_archive_bytes)throw std::length_error("PubChem complete publication archive byte budget");
    cancelled(stop);return publication;
}
} // namespace gadget

namespace gadget::examples {
namespace {
class PubChemEvidenceAdapter final : public EvidenceAdapter {
  public:
    void verify(const NativeRecord &source, const rfc::database::InterpretedSnapshot &native) const override {
        verify_pubchem_evidence(source, native);
    }
    void validate(const Json &query, const std::string &digest) const override {
        validate_pubchem_query(query, digest);
    }
    NativeRecord execute(const NativeRecord &source, const rfc::database::InterpretedSnapshot &native,
                         const NativeRecord &query) const override {
        return query_pubchem_evidence(source, native, query);
    }
    NativeRecord context(const NativeRecord &source, const rfc::database::InterpretedSnapshot &native,
                         const NativeRecord &query) const override {
        return context_pubchem_evidence(source, native, query);
    }
};
}
std::shared_ptr<const EvidenceAdapter> pubchem_evidence_adapter() {
    static const auto adapter = std::make_shared<const PubChemEvidenceAdapter>();
    return adapter;
}
}
