#include "pubchem_model.hpp"
#include "pubchem.hpp"
#include <charconv>
#include <set>
#include <rfc/database/interpretation_codec.hpp>

namespace gadget::examples {
namespace {
using namespace rfc::database;
using D = rfc::composite::Declaration;
using V = rfc::composite::DeclaredValue;
using I = rfc::numeric::Integer;
void require(bool condition, const char *message) {
    if (!condition) throw std::invalid_argument(message);
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested()) throw std::runtime_error("PubChem native admission cancelled");
}
void members(const Json &value, const std::set<std::string> &required,
             const std::set<std::string> &optional = {}) {
    require(json_object_is_type(value.raw(), json_type_object), "PubChem declaration must be an object");
    for (const auto &name : required) require(value.has(name), "missing PubChem declaration member");
    auto current = json_object_iter_begin(value.raw()), end = json_object_iter_end(value.raw());
    while (!json_object_iter_equal(&current, &end)) {
        const std::string name = json_object_iter_peek_name(&current);
        require(required.contains(name) || optional.contains(name), "undeclared PubChem request member");
        json_object_iter_next(&current);
    }
}
std::string identity(const Json &value) {
    auto token = value.text();
    require(!token.empty() && token.find('\0') == std::string::npos, "empty or NUL PubChem source identity");
    return token;
}
std::string public_id(const Json &value) {
    const auto token = identity(value);
    std::uint64_t number{};
    const auto parsed = std::from_chars(token.data(), token.data() + token.size(), number);
    require(parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size() && number > 0 &&
                number <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) &&
                token == std::to_string(number), "PubChem identifier is not a canonical positive integer");
    return token;
}
std::size_t budget(const Json &limits, const char *name) {
    const auto n = limits.at(name).integer();
    require(n >= 0 && static_cast<std::uint64_t>(n) <= std::numeric_limits<std::size_t>::max(),
            "PubChem resource limit is not addressable");
    return static_cast<std::size_t>(n);
}
void bound(std::size_t n, std::size_t maximum, const char *message) {
    if (n > maximum) throw std::length_error(message);
}
void numeric_literal(const std::string &token, std::size_t digits_limit, std::size_t scale_limit) {
    std::size_t digits{}, fractional{}, exponent{};
    bool point{}, negative_exponent{};
    std::size_t i{};
    if (i < token.size() && token[i] == '-') ++i;
    for (; i < token.size() && token[i] != 'e' && token[i] != 'E'; ++i) {
        if (token[i] == '.') point = true;
        else { ++digits; if (point) ++fractional; }
    }
    bound(digits, digits_limit, "PubChem numeric coefficient digit budget");
    if (i < token.size()) {
        ++i;
        if (i < token.size() && (token[i] == '+' || token[i] == '-')) negative_exponent = token[i++] == '-';
        const auto parsed = std::from_chars(token.data() + i, token.data() + token.size(), exponent);
        require(parsed.ec == std::errc{} && parsed.ptr == token.data() + token.size(), "PubChem numeric exponent is not addressable");
    }
    if (negative_exponent) {
        bound(fractional, scale_limit, "PubChem fractional scale budget");
        bound(exponent, scale_limit - fractional, "PubChem numeric scale budget");
    } else bound(exponent >= fractional ? exponent - fractional : fractional - exponent,
                 scale_limit, "PubChem numeric scale budget");
}
void numeric_metadata(const Json &root, std::size_t digits, std::size_t scale, std::stop_token stop) {
    std::vector<Json> pending{root};
    while (!pending.empty()) {
        cancelled(stop);
        auto value = pending.back(); pending.pop_back();
        if (json_object_is_type(value.raw(), json_type_int) || json_object_is_type(value.raw(), json_type_double))
            numeric_literal(value.dump_compact(), digits, scale);
        else if (json_object_is_type(value.raw(), json_type_array)) {
            for (std::size_t i = 0; i < value.size(); ++i) pending.push_back(value.at(i));
        } else if (json_object_is_type(value.raw(), json_type_object)) {
            auto current = json_object_iter_begin(value.raw()), end = json_object_iter_end(value.raw());
            while (!json_object_iter_equal(&current, &end)) {
                pending.emplace_back(json_object_get(json_object_iter_peek_value(&current)));
                json_object_iter_next(&current);
            }
        }
    }
}
void metadata_budget(const Json &root,std::size_t byte_limit,std::size_t work_limit,std::stop_token stop) {
    std::size_t bytes{},work{};std::vector<Json> pending{root};
    const auto add=[&](std::size_t n) {if(n>byte_limit-bytes)throw std::length_error("PubChem source metadata byte budget");bytes+=n;};
    const auto text=[&](const std::string &value) {add(2);for(const unsigned char c:value)add(c<32?6:(c=='"'||c=='\\'||c=='/'?2:1));};
    while(!pending.empty()) {
        cancelled(stop);if(work==work_limit)throw std::length_error("PubChem source metadata work budget");++work;
        auto value=pending.back();pending.pop_back();const auto type=json_object_get_type(value.raw());
        if(type==json_type_object) {
            add(2);auto current=json_object_iter_begin(value.raw()),end=json_object_iter_end(value.raw());
            while(!json_object_iter_equal(&current,&end)) {text(json_object_iter_peek_name(&current));add(2);pending.emplace_back(json_object_get(json_object_iter_peek_value(&current)));json_object_iter_next(&current);}
        } else if(type==json_type_array) {add(2);for(std::size_t i=0;i<value.size();++i) {add(1);pending.push_back(value.at(i));}}
        else if(type==json_type_string)text(value.text());else add(value.dump_compact().size());
    }
}
V text_value(const std::string &text) {
    std::vector<V> bytes; bytes.reserve(text.size());
    for (const unsigned char byte : text) bytes.push_back(V::symbol(byte));
    return V::repeated(std::move(bytes));
}
InterpretationSchema schema(std::size_t source_bytes,std::uint64_t source_namespace) {
    const auto id = [](std::uint64_t key) { return rfc::composite::DeclarationId{pubchem_declaration_namespace,key,1}; };
    std::vector<std::uint64_t> bytes; for (std::uint64_t i=0;i<256;++i) bytes.push_back(i);
    InterpretationSchema out;
    out.source_namespace = source_namespace; out.binding_namespace = pubchem_binding_namespace;
    out.values = rfc::composite::DeclarationRegistry({
        D::boolean(id(1)), D::integer(id(2)), D::rational(id(3)),
        D::finite(id(4), {1,2,3,4}), D::finite(id(5), std::move(bytes)),
        D::repeated(id(6), id(5), 0, source_bytes),
        D::alternative(id(7), {{1,id(3),0},{2,id(2),0},{3,id(1),0},{4,id(6),0}}),
        D::optional(id(8),id(7)), D::optional(id(9),id(2)),
        D::product(id(10),{{1,id(2),0},{2,id(4),0},{3,id(9),0}}),
        D::product(id(11),{{1,id(2),0},{2,id(4),0},{3,id(8),0}}),
        D::finite(id(12),{1,2,3,4,5}), D::optional(id(13),id(12)),
        D::finite(id(14),{1,2,3})});
    out.subjects={{1,"separate reported PubChem source identities","namespaced_identity","UTF-8 token"}};
    const std::vector<ExternalRole> roles{{11,1,0,1,true,false},{12,1,0,1,true,false},
        {13,1,1,1,false,false},{14,1,0,1,true,false},{15,1,1,1,false,false}};
    const std::vector<ExternalFieldLabel> references{{11,"left"},{12,"right"},{13,"assay"},{14,"definition"},{15,"capture"}};
    const auto kind = [&](std::uint64_t key,const char *name,std::uint64_t payload,std::vector<ExternalFieldLabel> fields={}) {
        return ExternalAssertionKind{key,name,id(payload),roles,std::move(fields),"capture_occurrence","exact source capture and pointer",references};
    };
    out.assertions={kind(40,"complete assay description",1),
        kind(41,"reported readout definition",10,{{1,"tid"},{2,"type"},{3,"unit"}}),
        kind(42,"reported substance outcome",13),
        kind(43,"typed reported or absent readout",11,{{1,"tid"},{2,"type"},{3,"reported_value"}}),
        kind(44,"reported substance revocation",2), kind(45,"explicit reported identity link",14)};
    return out;
}
} // namespace
PubChemEvidencePacket prepare_pubchem_evidence(const Json &input,
    const std::string &assay_bytes, const std::optional<std::string> &concise_bytes, std::uint64_t source_namespace, std::stop_token stop) {
    cancelled(stop);
    PubChemEvidencePacket p; p.request = Json::parse(input.dump_compact());
    members(p.request,{"schema","source_id","aid","assay_sha256","provenance","limits"},{"concise_sha256"});
    require(p.request.at("schema").text()=="gadget.pubchem.admission.request.v1", "unsupported PubChem admission request");
    const auto source_id=identity(p.request.at("source_id")), aid=public_id(p.request.at("aid"));
    const auto provenance=p.request.at("provenance"), limits=p.request.at("limits");
    members(provenance,{"release","selection","license","assay_url"},{"concise_url"});
    for(const auto *name:{"release","selection","license","assay_url"}) (void)identity(provenance.at(name));
    require(provenance.has("concise_url")==concise_bytes.has_value(),"PubChem concise provenance differs from source selection");
    if(concise_bytes) (void)identity(provenance.at("concise_url"));
    members(limits,{"source_bytes","definitions","substances","readout_slots","concise_rows","links","assertions","numeric_digits","numeric_scale","native_records","native_bytes","native_work"});
    PubChemSourceLimits reader{budget(limits,"source_bytes"),budget(limits,"definitions"),budget(limits,"substances"),budget(limits,"readout_slots"),budget(limits,"concise_rows"),budget(limits,"links")};
    const auto numeric_digits=budget(limits,"numeric_digits"),numeric_scale=budget(limits,"numeric_scale"),native_bytes=budget(limits,"native_bytes"),assertions=budget(limits,"assertions");
    // Decoder admissions account for absent coordinates before projection.
    reader.readout_slots=std::min(reader.readout_slots,native_bytes/1024);
    reader.definitions=std::min(reader.definitions,native_bytes/512);
    reader.substances=std::min(reader.substances,native_bytes/1024);
    bound(numeric_digits,native_bytes,"PubChem digit budget exceeds native byte budget");
    bound(numeric_scale,native_bytes,"PubChem numeric scale budget exceeds native byte budget");
    bound(assay_bytes.size(),reader.source_bytes,"PubChem assay source byte budget");
    if(concise_bytes)bound(concise_bytes->size(),reader.source_bytes,"PubChem concise source byte budget");
    require(p.request.at("assay_sha256").text()=="sha256:"+sha256(assay_bytes),"PubChem assay source SHA-256 mismatch");
    require(p.request.has("concise_sha256")==concise_bytes.has_value(),"PubChem concise source/hash presence differs");
    if(concise_bytes) require(p.request.at("concise_sha256").text()=="sha256:"+sha256(*concise_bytes),"PubChem concise source SHA-256 mismatch");
    const auto assay=inspect_pubchem_assay(assay_bytes,aid,reader,stop);
    numeric_metadata(assay.original,numeric_digits,numeric_scale,stop);
    const auto concise=concise_bytes?std::optional{inspect_pubchem_concise(*concise_bytes,aid,reader,stop)}:std::nullopt;
    p.records=Json::array();
    p.manifest=Json::object();p.manifest.set("schema",Json::string("gadget.pubchem.source.v1"));
    for(const auto *name:{"source_id","aid","assay_sha256","provenance"})p.manifest.set(name,p.request.at(name));
    if(concise_bytes)p.manifest.set("concise_sha256",p.request.at("concise_sha256"));
    p.logical.revision="sha256:"+sha256(p.request.dump_compact());
    p.manifest.set("source_digest",Json::string(p.logical.revision));
    p.schema=schema(reader.source_bytes,source_namespace);
    auto &profile=p.limits;profile.max_records=budget(limits,"native_records");
    profile.max_origin_bytes=profile.max_metadata_bytes=profile.max_archive_bytes=native_bytes;
    profile.max_work=budget(limits,"native_work");
    profile.declarations.max_components=profile.max_records;profile.declarations.max_value_nodes=profile.max_work;
    profile.declarations.max_encoded_bytes=profile.declarations.max_materialized_bytes=native_bytes;
    profile.declarations.max_integer_bits=rfc::machine::checked_size_mul(native_bytes,8);
    profile.atlas.max_records=profile.max_records;profile.atlas.max_work=profile.max_work;
    profile.atlas.max_archive_bytes=profile.atlas.max_integer_bytes=native_bytes;
    auto &access=profile.atlas.accession;
    access.max_handles=access.max_selection_cells=profile.max_records;
    access.max_bucket_entries=access.max_index_entries=access.max_action_work=profile.max_work;
    access.max_descriptor_bytes=access.max_total_descriptor_bytes=access.max_binding_bytes=access.max_integer_bytes=native_bytes;
    const auto entity=[&](const std::string &token,const char *kind) {
        if(token.empty()||p.entities.contains(token))return;
        auto value=Json::object();value.set("id",Json::string(token));value.set("kind",Json::string(kind));
        p.entities.emplace(token,value);p.logical.subjects.push_back({1,token});
    };
    const auto capture="pubchem:source:"+source_id, assay_id="pubchem:aid:"+aid;
    entity(capture,"declared_source_capture");entity(assay_id,"assay_record");
    const auto add=[&](std::uint64_t kind,const char *record_type,const char *source_role,
                       const std::string &pointer,const Json &original,V payload,
                       const std::string &left,const std::string &right,const std::string &definition) {
        cancelled(stop);bound(p.records.size()+1,assertions,"PubChem native assertion budget");
        const auto token=source_id+"/"+source_role+pointer;
        auto row=Json::object();row.set("id",Json::string(token));row.set("kind",Json::integer(static_cast<std::int64_t>(kind)));
        row.set("record_type",Json::string(record_type));row.set("provider",Json::string("pubchem"));
        row.set("source_role",Json::string(source_role));row.set("source_pointer",Json::string(pointer));row.set("original",original);
        row.set("left",Json::string(left));row.set("right",Json::string(right));row.set("assay",Json::string(assay_id));
        row.set("definition",Json::string(definition));row.set("capture",Json::string(capture));
        entity(left,"reported_identity");entity(right,"reported_identity");entity(definition,"readout_definition");
        ExternalAssertion native{{kind,token},std::move(payload),{}};
        std::uint64_t role=11;
        for(const auto &participant:{left,right,assay_id,definition,capture}) {
            ExternalRoleValue value;value.role=role++;value.state=participant.empty()?ReferenceState::Missing:ReferenceState::Present;
            if(!participant.empty())value.occurrences.push_back({std::to_string(value.role),{1,participant},true});
            native.roles.push_back(std::move(value));
        }
        p.records.push(row);p.logical.assertions.push_back(std::move(native));return row;
    };
    const auto context = [](const Json &original, const std::string &separate) {
        auto out=Json::object();
        json_object_object_foreach(original.raw(), name, value) {
            (void)value;
            if(name!=separate) out.set(name,original.at(name));
        }
        return out;
    };
    const auto view=assay.projection, description=view.at("description");
    add(40,"description","assay",description.at("source_pointer").text(),context(description.at("original"),"results"),V::boolean(true),assay_id,"","");
    const auto definitions=view.at("definitions");std::map<std::string,Json> descriptors;
    for(std::size_t i=0;i<definitions.size();++i) {
        const auto d=definitions.at(i),original=d.at("original");const auto tid=d.at("tid").text();descriptors.emplace(tid,original);
        auto definition_row=add(41,"definition","assay",d.at("source_pointer").text(),original,
            V::product({V::integer(I(tid)),V::symbol(static_cast<std::uint64_t>(original.at("type").integer())),
                original.has("unit")?V::present(V::integer(I(std::to_string(original.at("unit").integer())))):V::absent()}),
            assay_id+":tid:"+tid,"",assay_id+":tid:"+tid);
        definition_row.set("tid",Json::string(tid));
    }
    const auto substances=view.at("substances");
    for(std::size_t i=0;i<substances.size();++i) {
        const auto occurrence=substances.at(i);const auto sid=occurrence.at("sid").text(),substance="pubchem:sid:"+sid,path=occurrence.at("source_pointer").text();
        entity(substance,"substance_record");
        auto row=add(42,"outcome","assay",path,context(occurrence.at("original"),"data"),occurrence.at("outcome_reported").boolean()
            ?V::present(V::symbol(static_cast<std::uint64_t>(occurrence.at("outcome").integer()))):V::absent(),substance,"","");
        row.set("sid",Json::string(sid));row.set("outcome_reported",occurrence.at("outcome_reported"));
        if(occurrence.has("outcome"))row.set("outcome",occurrence.at("outcome"));
        const auto readouts=occurrence.at("readouts");
        for(std::size_t j=0;j<readouts.size();++j) {
            const auto readout=readouts.at(j);const auto tid=readout.at("tid").text();const auto d=descriptors.at(tid);
            const auto type=d.at("type").integer();auto value=V::absent();
            if(readout.at("support").text()=="reported") {
                const auto scalar=readout.at("original").at("value").at(readout.at("value_member").text());
                auto typed= type==1?V::rational(exact_decimal(readout.at("number_token").text()))
                    :type==2?V::integer(I(readout.at("number_token").text()))
                    :type==3?V::boolean(scalar.boolean()):text_value(scalar.text());
                value=V::present(V::alternative(static_cast<std::uint64_t>(type),std::move(typed)));
            }
            row=add(43,"readout","assay",path+"/readout/"+tid,
                readout.has("original")?readout.at("original"):Json{},
                V::product({V::integer(I(tid)),V::symbol(static_cast<std::uint64_t>(type)),std::move(value)}),substance,"",assay_id+":tid:"+tid);
            row.set("source_pointer",occurrence.at("source_pointer"));
            if(readout.has("source_pointer"))row.set("source_item_pointer",readout.at("source_pointer"));
            row.set("sid",Json::string(sid));row.set("tid",Json::string(tid));
            row.set("definition_id",Json::string(assay_id+":tid:"+tid));row.set("support",readout.at("support"));
            if(readout.has("number_token"))row.set("number_token",readout.at("number_token"));
        }
    }
    const auto revoked=view.at("revocations");
    for(std::size_t i=0;i<revoked.size();++i) {
        const auto r=revoked.at(i);const auto sid=r.at("sid").text();
        auto row=add(44,"revocation","assay",r.at("source_pointer").text(),r.at("original"),V::integer(I(sid)),"pubchem:sid:"+sid,"","");row.set("sid",Json::string(sid));
    }
    if(concise) {
        const auto links=concise->projection.at("links");
        for(std::size_t i=0;i<links.size();++i) {
            const auto link=links.at(i);const auto meaning=link.at("relationship").text();
            const std::uint64_t code=meaning=="reported_substance_compound_link"?1:meaning=="reported_assay_target"?2:3;
            auto row=add(45,"link","concise","/records/"+std::to_string(link.at("source_record").integer())+"/links/"+std::to_string(i),link,V::symbol(code),link.at("left").text(),link.at("right").text(),"");
            row.set("source_record",link.at("source_record"));row.set("relationship",link.at("relationship"));
        }
    }
    metadata_budget(p.records,native_bytes,profile.max_work,stop);
    p.receipt=Json::object();p.receipt.set("schema",Json::string("gadget.pubchem.evidence.store.v1"));
    p.receipt.set("source_digest",Json::string(p.logical.revision));p.receipt.set("source_id",Json::string(source_id));
    p.receipt.set("records_sha256",Json::string("sha256:"+sha256(p.records.dump_compact())));
    p.receipt.set("assertions",Json::integer(static_cast<std::int64_t>(p.logical.assertions.size())));
    cancelled(stop);return p;
}
InterpretedSnapshot interpret_pubchem_evidence(const PubChemEvidencePacket &p,std::stop_token stop) {
    cancelled(stop);auto native=InterpretedSnapshot::interpret(p.schema,p.logical,p.limits);cancelled(stop);return native;
}
void validate_pubchem_query(const Json &query,const std::string &digest) {
    members(query,{"schema","source_digest","operation","max_records"},{"sid","tid"});
    require(query.at("schema").text()=="gadget.pubchem.query.v1" && query.at("source_digest").text()==digest,"PubChem query schema/source mismatch");
    const auto op=query.at("operation").text();
    require(std::set<std::string>{"inspect","descriptions","definitions","outcomes","readouts","revocations","links"}.contains(op),"unsupported typed PubChem query");
    (void)budget(query,"max_records");
    if(query.has("sid")) {require(op=="readouts"||op=="outcomes"||op=="revocations","SID filter not declared for this query");(void)public_id(query.at("sid"));}
    if(query.has("tid")) {require(op=="readouts"||op=="definitions","TID filter not declared for this query");(void)public_id(query.at("tid"));}
}
void validate_pubchem_context(const Json &request,const std::string &digest) {
    members(request,{"schema","source_digest","question","assertion_keys","max_records"});
    require(request.at("schema").text()=="gadget.pubchem.context.request.v1" && request.at("source_digest").text()==digest,"PubChem context schema/source mismatch");
    (void)identity(request.at("question"));(void)budget(request,"max_records");
    const auto keys=request.at("assertion_keys");require(json_object_is_type(keys.raw(),json_type_array)&&keys.size()>0,"PubChem context requires assertion keys");
    bound(keys.size(),budget(request,"max_records"),"PubChem context result limit");std::set<std::string> distinct;
    for(std::size_t i=0;i<keys.size();++i)require(distinct.insert(public_id(keys.at(i))).second,"duplicate PubChem context key");
}
} // namespace gadget::examples
