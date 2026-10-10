#include <gadget/record_view.hpp>
#include "assay_quality_detail.hpp"
#include "assay_coverage_detail.hpp"
#include "genomic_detail.hpp"
#include <rfc/storage/mapped_file.hpp>
#include <rfc/storage/declaration_codec.hpp>
#include <cmath>

namespace gadget {
namespace {
using N=NativeRecordValue;
N import(const Json &value,const RecordByteBindings &bindings,std::size_t depth=1) {
    if(depth>native_data_limits().temporal.declarations.max_depth)
        throw std::length_error("input exceeds the native composite nesting profile");
    switch(json_object_get_type(value.raw())) {
    case json_type_null:return {};
    case json_type_boolean:return N(value.boolean());
    case json_type_int:return N(rfc::numeric::Integer(value.dump_compact()));
    case json_type_double:{const auto literal=value.dump_compact();return N(N::Decimal{exact_decimal(literal),literal});}
    case json_type_string:return N(value.text());
    case json_type_array:{
        N::Array out;out.elements.reserve(value.size());
        for(std::size_t i=0;i<value.size();++i)out.elements.push_back(import(value.at(i),bindings,depth+1));
        return N(std::move(out));
    }
    case json_type_object:{
        // This is the existing exact_number external representation. Only the
        // exact rational is mathematical data; a float is an optional view.
        const auto members=json_object_object_length(value.raw());
        if(value.has("rational") && (members==1 || (members==2 && value.has("approximation"))) &&
           json_object_is_type(value.at("rational").raw(),json_type_string)) {
            const auto exact=exact_rational(value.at("rational").text());
            if(value.has("approximation") && value.at("approximation").dump()!=Json::number(exact.get_d()).dump())
                throw std::invalid_argument("approximation differs from its declared exact rational");
            return N(N::ExactNumber{exact,value.has("approximation")});
        }
        N::Object out;
        json_object_object_foreach(value.raw(),key,entry) {
            (void)entry;const auto item=value.at(key);
            if(std::string_view(key)=="text" && value.has("document") && value.has("start") && value.has("end") &&
               bindings.documents.contains(value.at("document").text())) {
                const auto source=bindings.documents.at(value.at("document").text());
                const auto start=value.at("start").integer(),end=value.at("end").integer();
                if(start<0 || end<start)throw std::invalid_argument("invalid original evidence byte span");
                const auto owner=std::find_if(bindings.sources.begin(),bindings.sources.end(),[&](const auto &s){return s.key()==source;});
                if(owner==bindings.sources.end() || static_cast<std::uint64_t>(end)>owner->bytes().size())
                    throw std::invalid_argument("evidence byte source unavailable or span out of range");
                const auto span=owner->bytes().subspan(static_cast<std::size_t>(start),static_cast<std::size_t>(end-start));
                if(item.text()!=std::string(reinterpret_cast<const char *>(span.data()),span.size()))
                    throw std::invalid_argument("evidence quotation differs from original raw bytes");
                out.members.emplace_back(key,N(N::ByteSpan{source,static_cast<std::size_t>(start),static_cast<std::size_t>(end-start)}));
            } else out.members.emplace_back(key,import(item,bindings,depth+1));
        }
        return N(std::move(out));
    }
    }
    throw std::invalid_argument("unsupported external record constructor");
}
} // namespace
NativeRecordValue import_record_view(const Json &input,const RecordByteBindings &bytes) {return import(input,bytes);}
NativeRecord admit_record_view(const Json &input,RecordByteBindings bytes,const rfc::storage::TemporalEndpoint *source,
                              std::optional<std::uint64_t> source_namespace) {
    const auto data=import(input,bytes);
    return NativeRecord::create(data,source_namespace.value_or(new_native_record_namespace()),std::move(bytes.sources),source);
}
std::string native_declaration_digest(const NativeRecord &record) {
    const auto &declaration=record.endpoint()->declarations.value();
    const rfc::composite::CompiledLayout layout(declaration.registry,record.type(),
        declaration.profile,declaration.resource_limits);
    const auto schema=rfc::storage::encode_declarations(layout);
    const auto bits=rfc::storage::encode_declared_bits(layout,record.value());
    return "sha256:"+sha256_chunks([&](const auto &sink) {
        const auto count=std::to_string(schema.size())+"\n"+std::to_string(bits.bit_count)+"\n";
        sink(std::as_bytes(std::span(count.data(),count.size())));sink(schema);sink(bits.bytes);
        std::vector<NativeRecord> pending{record};
        while(!pending.empty()) {
            const auto item=pending.back();pending.pop_back();
            if(item.kind()==NativeRecordKind::Text || item.kind()==NativeRecordKind::Bytes)sink(item.bytes());
            else if(item.kind()==NativeRecordKind::Object) {
                const auto names=item.names();for(auto name=names.rbegin();name!=names.rend();++name)pending.push_back(item.at(*name));
            } else if(item.kind()==NativeRecordKind::Array) {
                for(auto i=item.size();i>0;--i)pending.push_back(item.at(i-1));
            }
        }
    });
}
Json render_record_view(const NativeRecord &record) {
    switch(record.kind()) {
    case NativeRecordKind::Null:return Json{};
    case NativeRecordKind::Boolean:return Json::boolean(record.boolean());
    case NativeRecordKind::Integer:{
        const auto n=record.integer();
        const rfc::numeric::Integer low(std::to_string(std::numeric_limits<std::int64_t>::min())),
                                   high(std::to_string(std::numeric_limits<std::int64_t>::max()));
        if(n<low || n>high)return Json::string(n.get_str());
        return Json::integer(std::stoll(n.get_str()));
    }
    case NativeRecordKind::Rational:return Json::string(record.rational().get_str());
    case NativeRecordKind::ExactNumber:{
        auto out=Json::object();out.set("rational",Json::string(record.rational().get_str()));
        if(record.approximate_view())out.set("approximation",Json::number(record.rational().get_d()));
        return out;
    }
    case NativeRecordKind::Decimal:{
        const auto literal=record.text();
        if(exact_decimal(literal)!=record.rational())throw std::invalid_argument("native decimal literal/exact value mismatch");
        return Json::parse(literal+" ");
    }
    case NativeRecordKind::Text:return Json::string(record.text());
    case NativeRecordKind::Bytes:{auto out=Json::object();out.set("byte_length",Json::string(std::to_string(record.bytes().size())));return out;}
    case NativeRecordKind::Object:{
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.alignment.fields.v1")
            return assay_detail::alignment_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.feature.count.fields.v1")
            return assay_detail::feature_count_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()==assay_detail::coverage_fields_schema)
            return assay_detail::coverage_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.assay.result.native.v1") {
            auto out=Json::object();out.set("schema",Json::string("gadget.assay.result.v1"));
            for(const auto *name:{"query","source_digest","result","arithmetic","gadget_sha256","rexgraph_sha256"})
                out.set(name,render_record_view(record.at(name)));
            out.set("certificate_sha256",Json::string("sha256:"+sha256(out.dump())));
            return out;
        }
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()==assay_detail::quality_fields_schema)
            return assay_detail::quality_fields_view(record);
        auto out=Json::object();for(const auto &name:record.names())out.set(name,render_record_view(record.at(name)));return out;
    }
    case NativeRecordKind::Array:{auto out=Json::array();for(std::size_t i=0;i<record.size();++i)out.push(render_record_view(record.at(i)));return out;}
    case NativeRecordKind::Tensor:{
        const auto &tensor=record.tensor();auto out=Json::object(),shape=Json::array(),support=Json::array(),values=Json::array();
        for(const auto n:tensor.domain().extents())shape.push(Json::string(std::to_string(n)));
        tensor.support().for_each_active([&](std::size_t index){support.push(Json::string(std::to_string(index)));});
        for(const auto &e:tensor.values().entries()){auto row=Json::object();row.set("index",Json::string(std::to_string(e.index)));row.set("value",Json::string(e.value.get_str()));values.push(row);}
        out.set("shape",shape);out.set("support",support);out.set("values",values);return out;
    }
    }
    throw std::invalid_argument("unknown native record view");
}
RetainedRecord read_record(const fs::path &path) {
    const rfc::storage::ReadOnlyMappedFile file(path);
    return decode_record(file.bytes());
}
namespace {
template<class Mutation>
NativeRecordValue update_native_container(NativeRecordValue value,const Mutation &mutation) {
    if(auto *owned=std::get_if<NativeRecordValue::Owned>(&value.data)) {
        owned->payload=std::make_shared<const NativeRecordValue>(
            update_native_container(*owned->payload,mutation));
    } else {
        mutation(value);
    }
    return value;
}
}
void RetainedRecord::set(const std::string &name,const Json &value) {
    auto updated=update_native_container(native().materialize_value(),[&](auto &container) {
        auto *object=std::get_if<NativeRecordValue::Object>(&container.data);
        if(!object)throw std::invalid_argument("native set requires a product");
        auto replacement=import_record_view(value);
        for(auto &member:object->members)if(member.first==name) {
            member.second=std::move(replacement);return;
        }
        object->members.emplace_back(name,std::move(replacement));
    });
    record_=native().derive(updated);
}
void RetainedRecord::push(const Json &value) {
    auto updated=update_native_container(native().materialize_value(),[&](auto &container) {
        auto *array=std::get_if<NativeRecordValue::Array>(&container.data);
        if(!array)throw std::invalid_argument("native push requires a repeated value");
        array->elements.push_back(import_record_view(value));
    });
    record_=native().derive(updated);
}
RetainedRecord decode_record(std::span<const std::byte> raw) {
    if(detect_native_data(raw)==NativeDataFormat::TemporalSeries)return NativeRecord::decode(raw);
    const std::string bytes(reinterpret_cast<const char *>(raw.data()),raw.size());
    return admit_record_view(Json::parse(bytes+" "),{{RawByteSource(1,raw)},{} });
}
void write_record(const fs::path &path,const Json &view,const rfc::storage::TemporalEndpoint *source) {
    admit_record_view(view,{},source).write(path);
}
} // namespace gadget
