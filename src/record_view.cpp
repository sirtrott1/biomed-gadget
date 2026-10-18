#include "assay_coverage_detail.hpp"
#include "assay_field_detail.hpp"
#include "assay_result_fields.hpp"
#include "assay_quality_detail.hpp"
#include "genomic_detail.hpp"
#include <cmath>
#include <gadget/record_view.hpp>
#include "scientific_result_fields.hpp"
#include <rfc/storage/declaration_codec.hpp>
#include <rfc/storage/mapped_file.hpp>

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
               json_object_is_type(value.at("document").raw(),json_type_string) &&
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
    const rfc::composite::PreparedLayoutRegistry registry(declaration.registry, declaration.resource_limits);
    std::map<rfc::composite::DeclarationId, rfc::composite::CompiledLayout> layouts;
    const auto layout_for = [&](rfc::composite::DeclarationId type) -> const rfc::composite::CompiledLayout & {
        if (!layouts.contains(type))
            layouts.emplace(type, rfc::composite::CompiledLayout(registry, type, declaration.profile));
        return layouts.at(type);
    };
    const auto &layout = layout_for(record.type());
    const auto schema=rfc::storage::encode_declarations(layout);
    const auto bits=rfc::storage::encode_declared_bits(layout,record.value());
    return "sha256:"+sha256_chunks([&](const auto &sink) {
        const auto count=std::to_string(schema.size())+"\n"+std::to_string(bits.bit_count)+"\n";
        sink(std::as_bytes(std::span(count.data(),count.size())));sink(schema);sink(bits.bytes);
        struct Frame {
            NativeRecord record;
            std::vector<std::string> names;
            std::size_t next{};
            bool entered{};
            bool packed_scope{};
        };
        std::vector<Frame> pending{{record, {}, 0, false, false}};
        std::map<rfc::composite::DeclarationId, std::vector<std::byte>> packed_schemas;
        while(!pending.empty()) {
            auto &frame = pending.back();
            const auto kind = frame.record.kind();
            if (!frame.entered) {
                frame.entered = true;
                if (kind == NativeRecordKind::Text || kind == NativeRecordKind::Bytes)
                    sink(frame.record.bytes());
                if (kind == NativeRecordKind::Object) frame.names = frame.record.names();
                else if (kind != NativeRecordKind::Array) {
                    pending.pop_back();
                    continue;
                }
            }
            const auto size = kind == NativeRecordKind::Object ? frame.names.size() : frame.record.size();
            if (frame.next == size) { pending.pop_back(); continue; }
            const auto index = frame.next++;
            const auto &container = declaration.registry.at(frame.record.type());
            auto stored_type = kind == NativeRecordKind::Object ? container.factors.at(index).type
                                                               : container.element;
            auto packed_scope = frame.packed_scope;
            if (kind == NativeRecordKind::Array) {
                if (const auto encoding = frame.record.packed_encoding(index)) {
                    stored_type = encoding->type;
                    packed_scope = true;
                    if (!packed_schemas.contains(encoding->type)) {
                        packed_schemas.emplace(encoding->type,
                            rfc::storage::encode_declarations(layout_for(encoding->type)));
                    }
                    const auto &schema = packed_schemas.at(encoding->type);
                    const auto framing = "gadget.packed.value\n" + std::to_string(schema.size()) +
                        "\n" + std::to_string(encoding->bit_count) + "\n";
                    sink(std::as_bytes(std::span(framing.data(), framing.size())));
                    sink(schema);
                    sink(encoding->bytes);
                } else if (declaration.registry.at(stored_type).kind == rfc::composite::DeclarationKind::Alternative) {
                    const auto &elements = std::get<rfc::composite::DeclaredValue::Repeated>(frame.record.value().data).elements;
                    const auto &alternative = std::get<rfc::composite::DeclaredValue::Alternative>(elements.at(index).data);
                    const auto ordinal = declaration.registry.alternative_ordinal(stored_type, alternative.key);
                    stored_type = declaration.registry.at(stored_type).factors.at(ordinal).type;
                }
            }
            auto child = kind == NativeRecordKind::Object ? frame.record.at(frame.names.at(index))
                                                         : frame.record.at(index);
            // A packet may name an owning declaration outside its encoded
            // tree. Bind that resolved value as well as its reference address.
            if (packed_scope && stored_type != child.type()) {
                if (!packed_schemas.contains(child.type()))
                    packed_schemas.emplace(child.type(), rfc::storage::encode_declarations(layout_for(child.type())));
                const auto &resolved_schema = packed_schemas.at(child.type());
                const auto resolved_bits = rfc::storage::encode_declared_bits(layout_for(child.type()), child.value());
                const auto framing = "gadget.resolved.value\n" + std::to_string(resolved_schema.size()) +
                    "\n" + std::to_string(resolved_bits.bit_count) + "\n";
                sink(std::as_bytes(std::span(framing.data(), framing.size())));
                sink(resolved_schema);
                sink(resolved_bits.bytes);
            }
            pending.push_back({std::move(child), {}, 0, false, packed_scope});
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
        if (record.has("schema") && record.at("schema").kind() == NativeRecordKind::Text) {
            const auto schema = record.at("schema").text();
            if (schema == "gadget.scientific.result.fields.v1") {
                auto out = render_record_view(record.at("certificate"));
                const auto fields = record.at("fields");
                const auto view = fields.at("view");
                const auto rendered = render_record_view(fields);
                if (view.kind() == NativeRecordKind::Object && view.has("schema") &&
                    view.at("schema").text() == "gadget.scientific.direct.result.v1") {
                    out.set("result", rendered.at("result"));
                    out.set("output_grade", rendered.at("output_grade"));
                    for (const auto &name : {"green_iterations", "green_execution", "residual_verified"})
                        if (rendered.has(name)) out.set(name, rendered.at(name));
                } else if (view.kind() == NativeRecordKind::Object && view.has("schema") &&
                           view.at("schema").text() == "gadget.scientific.native.boundary.result.v1") {
                    for (const auto &name : {"result", "field_set", "retained_occurrences", "mathematical_definition"})
                        out.set(name, rendered.at(name));
                } else out.set("result", rendered);
                if (fields.has("execution"))
                    out.set("execution", render_record_view(fields.at("execution")));
                return scientific_detail::seal(out, "certificate_sha256");
            }
            if (schema == "gadget.scientific.plan.fields.v1") {
                auto out = Json::object(), steps = Json::array();
                for (std::size_t i = 0; i < record.at("steps").size(); ++i) {
                    const auto retained = record.at("steps").at(i);
                    auto step = Json::object();
                    step.set("id", render_record_view(retained.at("id")));
                    const auto publication = NativeRecord::decode(retained.at("publication").bytes());
                    step.set("certificate", render_record_view(publication.at("result")));
                    steps.push(step);
                }
                const auto plan = render_record_view(record.at("plan"));
                out.set("schema", Json::string("gadget.scientific.plan.certificate.v1"));
                out.set("status", Json::string("complete"));
                out.set("plan", plan);
                out.set("plan_sha256", Json::string("sha256:" + sha256(plan.dump())));
                out.set("source_digest", render_record_view(record.at("source_digest")));
                out.set("steps", steps);
                return scientific_detail::seal(out, "certificate_sha256");
            }
            if (schema == "gadget.scientific.fields.v1")
                return scientific_detail::result_fields_view(record);
            if (schema == "gadget.scientific.tensor.view.v1")
                return scientific_detail::result_tensor_view(record);
            if (schema == "gadget.sequence.alignment.fields.v1")
                return assay_detail::sequence_alignment_fields_view(record);
            if (schema == "gadget.sequence.assembly.fields.v1")
                return assay_detail::sequence_assembly_fields_view(record);
            if (schema == "gadget.constraints.fields.v1")
                return assay_detail::constraint_fields_view(record);
            if (schema == "gadget.genotype.fields.v1")
                return assay_detail::genotype_fields_view(record);
            if (schema == "gadget.chromatin.fields.v1")
                return assay_detail::chromatin_fields_view(record);
            if (schema == "gadget.intervals.coverage.fields.v1")
                return assay_detail::interval_fields_view(record);
            if (schema == "gadget.annotation.features.fields.v1")
                return assay_detail::annotation_fields_view(record);
            if (schema == "gadget.sequence.motifs.fields.v1")
                return assay_detail::motif_fields_view(record);
            if (schema == "gadget.assay.batch.fields.v1")
                return assay_detail::query_batch_fields_view(record);
            if (schema == "gadget.persistence.fields.v1")
                return assay_detail::persistence_fields_view(record);
            if (schema == "gadget.sequence.indexed.fields.v1")
                return assay_detail::indexed_alignment_fields_view(record);
            if (schema == "gadget.sequence.reference.fields.v1")
                return assay_detail::reference_match_fields_view(record);
            if (schema == "gadget.sequence.splice.fields.v1")
                return assay_detail::splice_search_fields_view(record);
            if (schema == "gadget.sequence.read_cover.fields.v1")
                return assay_detail::read_reconstruction_fields_view(record);
            if (schema == "gadget.allocation.fields.v1")
                return assay_detail::allocation_bounds_fields_view(record);
            if (schema == "gadget.assays.gluing.fields.v1")
                return assay_detail::gluing_fields_view(record);
            if (schema == "gadget.assays.transfer.fields.v1")
                return assay_detail::transfer_fields_view(record);
            if (schema == "gadget.sequence.consensus.fields.v1")
                return assay_detail::read_consensus_fields_view(record);
            if (schema == "gadget.sections.fields.v1" || schema == "gadget.sections.selected.v1")
                return assay_detail::section_fields_view(record);
        }
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.alignment.fields.v1")
            return assay_detail::alignment_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.feature.count.fields.v1")
            return assay_detail::feature_count_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.expression.fields.v1")
            return assay_detail::expression_fields_view(record);
        if(record.has("schema") && record.at("schema").kind()==NativeRecordKind::Text &&
           record.at("schema").text()=="gadget.aggregation.fields.v1")
            return assay_detail::aggregation_fields_view(record);
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
    index_.reset();
}
void RetainedRecord::push(const Json &value) {
    auto updated=update_native_container(native().materialize_value(),[&](auto &container) {
        auto *array=std::get_if<NativeRecordValue::Array>(&container.data);
        if(!array)throw std::invalid_argument("native push requires a repeated value");
        array->elements.push_back(import_record_view(value));
    });
    record_=native().derive(updated);
    index_.reset();
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
