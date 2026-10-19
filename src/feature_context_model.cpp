#include "feature_context_model.hpp"
#include <algorithm>
#include <charconv>
#include <set>

namespace gadget::feature_context_detail {
namespace {
void require(bool condition,const char *message) {if(!condition)throw std::invalid_argument(message);}
void bound(std::size_t n,std::size_t limit,const char *message) {if(n>limit)throw std::length_error(message);}
void members(const Json &value,const std::set<std::string> &required) {
    require(json_object_is_type(value.raw(),json_type_object),"context declaration must be an object");
    require(static_cast<std::size_t>(json_object_object_length(value.raw()))==required.size(),"undeclared context member");
    for(const auto &name:required)require(value.has(name),"missing context member");
}
std::string text(const Json &value) {
    auto out=value.text();require(!out.empty()&&out.find('\0')==std::string::npos,"empty or NUL context identity");
    (void)Json::parse(Json::string(out).dump_compact());return out;
}
std::size_t count(const Json &value) {
    const auto n=value.integer();require(n>=0&&static_cast<std::uint64_t>(n)<=std::numeric_limits<std::size_t>::max(),"context count is not addressable");return static_cast<std::size_t>(n);
}
std::string digest(const Json &value) {
    const auto out=text(value);require(out.starts_with("sha256:")&&out.size()==71&&out.substr(7).find_first_not_of("0123456789abcdef")==std::string::npos,"invalid context SHA-256");return out;
}
Json path(const Json &value, bool empty = false) {
    require(json_object_is_type(value.raw(),json_type_array)&&(empty||value.size()>0),"context feature path must be an array");
    for(std::size_t i=0;i<value.size();++i) {
        const auto item=value.at(i);
        require(json_object_is_type(item.raw(),json_type_string) ||
            (json_object_is_type(item.raw(),json_type_int)&&item.integer()>=0),"native feature path component");
    }
    return value;
}
std::uint64_t key(const Json &value) {
    const auto token=text(value);std::uint64_t n{};const auto parsed=std::from_chars(token.data(),token.data()+token.size(),n);
    require(parsed.ec==std::errc{}&&parsed.ptr==token.data()+token.size()&&n>0&&std::to_string(n)==token,"context assertion key must be positive canonical text");return n;
}
} // namespace
void Work::check() const {if(stop_.stop_requested())throw std::runtime_error("feature evidence context cancelled");}
void Work::take(std::size_t n) {check();if(n>limit_-used_)throw std::length_error("feature evidence context work budget");used_+=n;}
Request admit_request(const Json &q) {
    members(q,{"schema","analysis_sha256","mapping_sha256","evidence_digest","feature_path","feature_id_path","assertion_keys","question","limits"});
    require(q.at("schema").text()=="gadget.features.evidence-context.request.v1","unsupported feature context request");
    Request out;out.analysis_sha256=digest(q.at("analysis_sha256"));out.mapping_sha256=digest(q.at("mapping_sha256"));out.evidence_digest=digest(q.at("evidence_digest"));out.question=text(q.at("question"));
    out.feature_path=path(q.at("feature_path"));out.feature_id_path=path(q.at("feature_id_path"),true);
    const auto limits=q.at("limits");members(limits,{"source_bytes","features","mapping_rows","assertions","participants","matches","work","output_bytes"});
    out.limits={count(limits.at("source_bytes")),count(limits.at("features")),count(limits.at("mapping_rows")),count(limits.at("assertions")),count(limits.at("participants")),count(limits.at("matches")),count(limits.at("work")),count(limits.at("output_bytes"))};
    const auto keys=q.at("assertion_keys");require(json_object_is_type(keys.raw(),json_type_array)&&keys.size()>0,"context requires explicit selected assertions");bound(keys.size(),out.limits.assertions,"context selected assertion budget");std::set<std::uint64_t> distinct;
    for(std::size_t i=0;i<keys.size();++i){const auto value=key(keys.at(i));require(distinct.insert(value).second,"duplicate context assertion selection");out.assertions.push_back(value);}return out;
}
Roles roles(const rfc::database::InterpretationSchema &schema) {
    Roles out;
    for (const auto &kind : schema.assertions)
        for (const auto &role : kind.reference_labels) {
            const auto previous = out.find(role.label);
            require(previous == out.end() || previous->second == role.factor,
                    "one role label names different native factors; explicit source profile required");
            out.emplace(role.label, role.factor);
        }
    require(!out.empty(), "native evidence has no declared role labels");
    return out;
}
std::vector<MappingRow> admit_mapping(const Json &d,const std::string &source,const Roles &declared_roles,const Limits &limits,Work &work) {
    members(d,{"schema","id","origin","release","source_url","license","format","sha256"});
    require(d.at("schema").text()=="gadget.features.mapping.native.v1"&&d.at("format").text()=="tsv","unsupported retained feature mapping source");
    for(const auto *name:{"id","release","source_url","license"})(void)text(d.at(name));
    const auto origin=d.at("origin").text();require(origin=="released"||origin=="explicit_declared","mapping requires released or explicit declared origin");
    bound(source.size(),limits.source_bytes,"mapping source byte budget");work.take(source.size());
    require(digest(d.at("sha256"))=="sha256:"+sha256(source),"mapping original source SHA-256 differs");
    const auto table=read_csv(source,'\t');
    for(const auto *column:{"feature_id","evidence_id","role"}) {
        require(table.columns.contains(column)&&table.header.at(table.columns.at(column))==column,"mapping source columns missing or renamed");
    }
    bound(table.rows.size(),limits.mapping_rows,"mapping source occurrence budget");
    // read_csv owns grammar/decoding. Audit original record intervals only,
    // respecting quoted fields, doubled quotes and embedded line endings.
    std::vector<std::pair<std::size_t,std::size_t>> spans;bool quoted{};std::size_t begin{};
    for(std::size_t i=0;i<source.size();++i) {
        work.take();const auto c=source[i];
        if(c=='"') {if(quoted&&i+1<source.size()&&source[i+1]=='"'){++i;work.take();}else quoted=!quoted;}
        else if(!quoted&&(c=='\n'||c=='\r')) {if(c=='\r'&&i+1<source.size()&&source[i+1]=='\n'){++i;work.take();}spans.emplace_back(begin,i+1-begin);begin=i+1;}
    }
    if(begin<source.size())spans.emplace_back(begin,source.size()-begin);
    require(!quoted&&spans.size()==table.rows.size()+1,"mapping source contains omitted/blank or unmatched records");
    const auto &declared=declared_roles;std::vector<MappingRow> result;
    for(std::size_t i=0;i<table.rows.size();++i) {
        work.take();const auto &row=table.rows[i];const auto feature=text(Json::string(table.get(row,"feature_id"))),identity=text(Json::string(table.get(row,"evidence_id"))),role=text(Json::string(table.get(row,"role")));
        const auto colon=identity.find(':');require(colon!=std::string::npos&&colon>0&&colon+1<identity.size(),"mapping evidence identity must retain an explicit namespace");
        require(declared.contains(role),"mapping role is unsupported by the actual evidence profile");
        result.push_back({row.record,spans.at(i+1).first,spans.at(i+1).second,feature,identity,role});
    }
    work.check();return result;
}
std::vector<Assertion> project_assertions(const rfc::database::InterpretedSnapshot &source,const Request &request,const Roles &declared_roles,Work &work) {
    const auto &declared=declared_roles;
    const std::set<std::uint64_t> selected_keys(request.assertions.begin(),request.assertions.end());
    std::map<std::uint64_t,std::uint64_t> assertions;
    for(const auto &binding:source.lineage().assertions) {
        work.take();
        if(selected_keys.contains(binding.native_key)) assertions.emplace(binding.native_key,binding.origin.kind);
    }
    using Position=std::pair<std::size_t,std::size_t>;
    std::map<std::uint64_t,std::vector<Position>> subjects;
    std::vector<Assertion> result;std::size_t participants{};
    for(const auto selected:request.assertions) {
        work.take();require(assertions.contains(selected),"context selects an unknown actual native assertion");
        result.push_back({selected,assertions.at(selected),{}});auto &row=result.back();
        const auto &kinds=source.schema().assertions;
        const auto found=std::find_if(kinds.begin(),kinds.end(),[&](const auto &kind){return kind.key==row.kind;});
        require(found!=kinds.end(),"native assertion kind is undeclared");
        for(const auto &[name,role]:declared) {
            work.take();if(std::none_of(found->roles.begin(),found->roles.end(),[&](const auto &r){return r.key==role;}))continue;
            std::size_t position{};
            for(const auto &actual:source.participants(selected,role)) {
                work.take();bound(++participants,request.limits.participants,"context actual participant budget");
                subjects[actual.participant().key].emplace_back(result.size()-1,row.participants.size());
                row.participants.push_back({{},name,actual.key().key.value,position++,actual.aspect("exists")});
            }
        }
    }
    // Retain an index only for the selected participant occurrences.
    for(const auto &binding:source.lineage().subjects) {
        work.take();const auto found=subjects.find(binding.native_key);
        if(found==subjects.end())continue;
        for(const auto &[a,p]:found->second) {
            work.take();result[a].participants[p].identity=binding.origin.token;
        }
        subjects.erase(found);
    }
    require(subjects.empty(),"native role references an unknown subject");
    return result;
}
Join join_features(const std::vector<std::string> &features,const std::vector<MappingRow> &mapping,const std::vector<Assertion> &assertions,const Limits &limits,Work &work) {
    bound(features.size(),limits.features,"context analysis feature occurrence budget");bound(mapping.size(),limits.mapping_rows,"context mapping occurrence budget");bound(assertions.size(),limits.assertions,"context assertion budget");
    using Candidate=std::pair<std::string,std::string>;using Reference=std::pair<std::size_t,std::size_t>;
    std::map<Candidate,std::vector<Reference>> index;std::size_t participant_count{};
    for(std::size_t a=0;a<assertions.size();++a)for(std::size_t p=0;p<assertions[a].participants.size();++p) {work.take();bound(++participant_count,limits.participants,"context actual participant budget");const auto &value=assertions[a].participants[p];index[{value.identity,value.role}].push_back({a,p});}
    std::map<std::string,std::vector<std::size_t>> rows;for(std::size_t r=0;r<mapping.size();++r){work.take();rows[mapping[r].feature].push_back(r);}
    Join out;out.row_status.assign(mapping.size(),"unused");std::size_t matches{};
    for(const auto &id:features) {
        work.take();require(!id.empty()&&id.find('\0')==std::string::npos,"analysis feature identity is empty or NUL");Feature feature;feature.id=id;const auto found=rows.find(id);
        if(found!=rows.end())feature.rows=found->second;
        std::set<Candidate> candidates;
        for(const auto r:feature.rows) {
            work.take();const auto &row=mapping[r];const Candidate candidate{row.identity,row.role};candidates.insert(candidate);const auto matched=index.find(candidate);
            if(matched==index.end())out.row_status[r]="unresolved";
            else {out.row_status[r]="matched";for(const auto &[a,p]:matched->second){work.take();bound(++matches,limits.matches,"context expanded correspondence budget");feature.matches.push_back({r,a,p});}}
        }
        feature.candidates=candidates.size();feature.status=feature.rows.empty()?"unmapped":candidates.size()>1?"ambiguous":feature.matches.empty()?"unresolved":"mapped";out.features.push_back(std::move(feature));
    }
    work.check();return out;
}
} // namespace gadget::feature_context_detail
