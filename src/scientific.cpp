#include "scientific_detail.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <limits>
#include <rfc/storage/mapped_file.hpp>
#include <rfc/storage/native_atlas_types.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::scientific_detail {
void require(bool valid, const std::string &message) {
    if (!valid)
        throw std::invalid_argument("scientific study: " + message);
}
void members(const Json &object, const std::set<std::string> &required,
             const std::set<std::string> &optional) {
    require(json_object_is_type(object.raw(), json_type_object), "expected object");
    for (const auto &name : required)
        require(object.has(name), "missing member: " + name);
    json_object_object_foreach(object.raw(), name, value) {
        (void)value;
        require(required.contains(name) || optional.contains(name),
                "unknown member: " + std::string(name));
    }
}
std::size_t count(const Json &value, std::size_t limit, const std::string &name) {
    const auto n = value.integer();
    require(n >= 0 && static_cast<std::uint64_t>(n) <= limit, name + " exceeds profile budget");
    return static_cast<std::size_t>(n);
}
Json cardinality(std::size_t value) {
    require(value <= static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()),
            "count cannot be represented in JSON");
    return Json::integer(static_cast<std::int64_t>(value));
}
Json exact_count(const rfc::numeric::Integer &value) {
    require(value >= 0, "negative exact count");
    static const rfc::numeric::Integer largest{
        std::to_string(std::numeric_limits<std::int64_t>::max())};
    return value <= largest ? Json::integer(std::stoll(value.get_str()))
                            : Json::string(value.get_str());
}
rfc::numeric::SparseRationalTensor declared_comparison(const Json &map, std::size_t rows,
                                                       std::size_t columns) {
    members(map, {"rows", "columns", "entries"});
    require(count(map.at("rows"), rows, "comparison rows") == rows &&
                count(map.at("columns"), columns, "comparison columns") == columns,
            "comparison dimensions differ from actual carriers");
    source_array(map.at("entries"), "comparison entries");
    const rfc::numeric::IndexDomain domain({rows, columns});
    std::vector<rfc::numeric::SparseRationalEntry> values;
    std::set<std::size_t> seen;
    values.reserve(map.at("entries").size());
    for (std::size_t i = 0; i < map.at("entries").size(); ++i) {
        const auto e = map.at("entries").at(i);
        members(e, {"row", "column", "value"});
        const auto row = count(e.at("row"), rows, "comparison row"),
                   column = count(e.at("column"), columns, "comparison column");
        require(row < rows && column < columns, "comparison entry outside actual carrier");
        const std::array<std::size_t, 2> coordinates{row, column};
        const auto index = domain.flatten(coordinates);
        require(seen.insert(index).second, "duplicate comparison entry");
        auto value = exact_rational(e.at("value").text());
        if (value != 0)
            values.push_back({index, std::move(value)});
    }
    return {domain, rfc::numeric::StructuralSupport::full(domain.logical_size()),
            rfc::numeric::SparseRationalVector(domain.logical_size(), std::move(values))};
}
rfc::storage::NativeSeriesLimits source_codec_limits(std::size_t bytes) {
    auto limits = native_data_limits();
    limits.max_archive_bytes = bytes;
    limits.max_segments = bytes;
    limits.max_checkpoints = bytes;
    limits.max_analysis_steps = bytes;
    limits.max_analysis_dimension = std::numeric_limits<std::size_t>::max();
    limits.max_analysis_scalars = std::numeric_limits<std::size_t>::max();
    limits.max_analysis_work = std::numeric_limits<std::size_t>::max();
    limits.temporal.max_file_bytes = bytes;
    // Encoded records, rank words and integer bytes cannot outnumber input
    // bytes. Sparse logical dimensions instead use the native index domain.
    limits.temporal.max_records = bytes;
    limits.temporal.max_rank = bytes;
    limits.temporal.max_integer_bytes = bytes;
    limits.temporal.max_dimension = std::numeric_limits<std::size_t>::max();
    limits.temporal.source_declarations.max_occurrences = bytes;
    limits.temporal.source_declarations.max_aspect_bytes = bytes;
    return limits;
}
std::uint64_t key(const Json &value) {
    const auto s = value.text();
    std::uint64_t n{};
    const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), n);
    require(!s.empty() && s.front() != '0' && error == std::errc{} && end == s.data() + s.size(),
            "native key must be a positive canonical decimal string");
    return n;
}
std::string text(const Json &value, std::size_t limit) {
    auto s = value.text();
    require(!s.empty() && s.size() <= limit && s.find('\0') == std::string::npos,
            "invalid text or text budget");
    return s;
}
Json seal(Json result, const std::string &name) {
    const auto serialized = result.dump();
    result = Json::parse(serialized);
    result.set(name, Json::string("sha256:" + sha256(result.dump())));
    return result;
}
fs::path local_artifact_path(const fs::path &root, const Json &record) {
    const fs::path relative(text(record.at("path")));
    require(!relative.is_absolute(), "artifact path must be relative");
    for (const auto &part : relative)
        require(part != "..", "artifact path escapes study root");
    auto path = root / relative;
    require(fs::is_regular_file(fs::symlink_status(path)), "artifact must be a regular file");
    auto canonical = fs::canonical(path).lexically_relative(fs::canonical(root));
    require(!canonical.empty() && *canonical.begin() != "..", "artifact escapes study root");
    return path;
}
std::string local_artifact(const fs::path &root, const Json &record, std::size_t limit) {
    const auto path = local_artifact_path(root, record);
    auto bytes = read_text(path, limit);
    require(record.at("sha256").text() == "sha256:" + sha256(bytes),
            "artifact SHA-256 mismatch: " + record.at("path").text());
    return bytes;
}
void source_array(const Json &value, const std::string &name) {
    require(json_object_is_type(value.raw(), json_type_array), name + " must be an array");
}
} // namespace gadget::scientific_detail

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc::storage;

bool same_metadata(const NativeRecord &record, const Json &value) {
    if (record.kind() == NativeRecordKind::Array) {
        if (!json_object_is_type(value.raw(), json_type_array) || record.size() != value.size())
            return false;
        for (std::size_t i = 0; i < value.size(); ++i)
            if (!same_metadata(record.at(i), value.at(i))) return false;
        return true;
    }
    if (record.kind() == NativeRecordKind::Object) {
        if (!json_object_is_type(value.raw(), json_type_object)) return false;
        const auto names = record.names();
        if (names.size() != static_cast<std::size_t>(json_object_object_length(value.raw()))) return false;
        std::size_t i{};
        json_object_object_foreach(value.raw(), name, member) {
            (void)member;
            if (names.at(i++) != name || !same_metadata(record.at(name), value.at(name)))
                return false;
        }
        return true;
    }
    return render_record_view(record).dump() == value.dump();
}
bool boundary(std::string_view bytes, std::size_t position) {
    return position == bytes.size() ||
           (static_cast<unsigned char>(bytes.at(position)) & 0xc0) != 0x80;
}
} // namespace
Json scientific_detail::source_evidence(const Json &value,
                                        const std::map<std::string, Document> &documents,
                                        bool require_utf8, bool allow_references) {
    source_array(value, "evidence");
    require(value.size() > 0, "every semantic/source record requires evidence");
    auto result = Json::array();
    for (std::size_t i = 0; i < value.size(); ++i) {
        auto span = value.at(i);
        if (allow_references)
            members(span, {"document", "start", "end"}, {"text"});
        else
            members(span, {"document", "start", "end", "text"});
        const auto document = text(span.at("document"));
        require(documents.contains(document), "unknown evidence document");
        const auto &source = documents.at(document);
        require(!require_utf8 || source.metadata.at("format").text() == "utf8",
                "source annotations require a UTF-8 document");
        const auto start = count(span.at("start"), source.bytes.size(), "span start");
        const auto end = count(span.at("end"), source.bytes.size(), "span end");
        require(start < end, "evidence span must be nonempty");
        if (span.has("text"))
            require(span.at("text").text() ==
                        std::string_view(source.bytes).substr(start, end - start),
                    "evidence quote/span mismatch");
        if (source.metadata.at("format").text() == "utf8")
            require(boundary(source.bytes, start) && boundary(source.bytes, end),
                    "span splits a UTF-8 code point");
        auto item = Json::parse(span.dump());
        item.set("document_sha256", source.metadata.at("sha256"));
        result.push(item);
    }
    return result;
}
Json ScientificStudy::State::evidence(const Json &value) const {
    return scientific_detail::source_evidence(
        value, documents, false, declaration.at("schema").text() != "gadget.scientific.study.v1");
}
ScientificStudy::State::State(const fs::path &path,std::optional<NativeRecord> publication,
                             std::optional<NativeTemporalSeries> producer,
                             const ScientificAdmissionObserver &observer) {
    const bool generated=producer.has_value();
    if(!publication) {
        require(fs::is_regular_file(fs::symlink_status(path)), "declaration must be a regular file");
        const ReadOnlyMappedFile input(path);
        if(detect_native_data(input.bytes())==NativeDataFormat::TemporalSeries)
            publication=NativeRecord::read(path);
        else {
            const auto sibling=path.parent_path()/"study.rfcxs";
            if(sibling!=path && fs::is_regular_file(fs::symlink_status(sibling))) {
                auto candidate=NativeRecord::read_if_record(sibling);
                if(candidate && candidate->has("schema") &&
                    (candidate->at("schema").text()=="gadget.biomedical.native.v1" ||
                     candidate->at("schema").text()=="gadget.biomedical.native.v2")) {
                    const auto original=candidate->at("declaration").bytes();
                    if(original.size()==input.bytes().size() &&
                        std::equal(original.begin(),original.end(),input.bytes().begin())) {
                        const auto input_root=fs::canonical(path).parent_path();
                        const auto view=candidate->at("study");
                        const auto check_original=[&](const NativeRecord &reference,std::span<const std::byte> expected) {
                            if(!fs::exists(input_root/reference.at("path").text()))return;
                            const auto source_path=local_artifact_path(input_root,render_record_view(reference));
                            const ReadOnlyMappedFile source_file(source_path);
                            require(source_file.bytes().size()==expected.size() &&
                                std::equal(expected.begin(),expected.end(),source_file.bytes().begin()),
                                "external source bytes differ from the admitted native publication");
                        };
                        if(candidate->at("schema").text()=="gadget.biomedical.native.v1")
                            check_original(view.at("native"),candidate->at("history").bytes());
                        else if(fs::exists(input_root/view.at("native").at("path").text())) {
                            const auto external=local_artifact_path(input_root,render_record_view(view.at("native")));
                            require("sha256:"+sha256_file(external)==view.at("native").at("sha256").text(),
                                    "external source bytes differ from the admitted native publication");
                        }
                        const auto documents=view.at("documents");
                        for(std::size_t i=0;i<documents.size();++i)
                            check_original(documents.at(i),candidate->at("documents").at(i).at("bytes").bytes());
                        publication=std::move(candidate);
                    }
                }
            }
        }
    }
    if(publication && publication->at("schema").text()=="gadget.biomedical.native.v3") {
        root=fs::absolute(path).parent_path();admit_compact(std::move(*publication), observer);return;
    }
    if(publication)require(publication->at("schema").text()=="gadget.biomedical.native.v1" ||
        publication->at("schema").text()=="gadget.biomedical.native.v2","unsupported native biomedical declaration");
    const bool embedded=publication && publication->at("schema").text()=="gadget.biomedical.native.v2";
    const auto absolute = publication?fs::absolute(path):fs::canonical(path);
    // Cycles are rejected by source identity, without a numerical nesting cap.
    struct AdmissionPath {
        std::set<fs::path> &paths;
        fs::path path;
        ~AdmissionPath() { paths.erase(path); }
    };
    static thread_local std::set<fs::path> admitting;
    require(admitting.insert(absolute).second, "cyclic study publication references");
    AdmissionPath active{admitting, absolute};
    root = absolute.parent_path();
    const auto bytes = publication?[&] {
        const auto b=publication->at("declaration").bytes();
        return std::string(reinterpret_cast<const char *>(b.data()),b.size());
    }():read_text(absolute);
    declaration_digest = "sha256:" + sha256(bytes);
    const auto declaration = Json::parse(bytes);
    if(publication)require(same_metadata(publication->at("study"), declaration),"native biomedical declaration differs from its retained original source");
    members(declaration,
            {"schema", "id", "native", "design", "documents", "cells", "occurrences", "contexts",
             "fields", "field_sets"},
            {"publications", "operators"});
    const auto schema = declaration.at("schema").text();
    require(schema == "gadget.scientific.study.v1" || schema == "gadget.scientific.study.v2",
            "unsupported schema");
    (void)text(declaration.at("id"));
    const auto native = declaration.at("native");
    members(native, {"path", "sha256", "source_namespace"}, {"cutoff"});
    native_digest = native.at("sha256").text();
    std::vector<std::byte> original_history;
    std::optional<NativeTemporalSeries> producer_history;
    {
        auto series = [&] {
            if(producer) {
                require(!publication,"generated source cannot replace an existing publication");
                const auto hash=sha256_chunks([&](const auto &sink){stream_native_series(*producer,sink);});
                require(native_digest=="sha256:"+hash,"generated history differs from its emitted bytes");
                return *producer;
            }
            if(publication) {
                if(embedded) {
                    const auto retained=publication->publication_history();
                    require(retained.has_value(),"embedded study has no owning publication history");
                    const auto binding=publication->at("history");
                    require(binding.at("kind").text()=="publication_prefix","unsupported native history binding");
                    const auto cutoff=binding.at("cutoff").integer();
                    require(cutoff>=0 && cutoff<=rfc::numeric::Integer(std::to_string(std::numeric_limits<std::uint64_t>::max())),
                            "embedded history cutoff outside native sequence domain");
                    const auto prefix=NativeSeriesPrefix::capture(*retained,std::stoull(cutoff.get_str()));
                    const auto prefix_tail=prefix.tail();
                    require(prefix.cutoff()<retained->last_sequence(),"annotation must follow its bound source prefix");
                    for(auto sequence=prefix.cutoff();sequence<retained->last_sequence();++sequence) {
                        const auto index=static_cast<std::size_t>(sequence-retained->first_sequence());
                        const auto update=decode_temporal_segment(retained->segment_bytes(index),
                            retained->limits().temporal);
                        require(update.source.grades.empty() &&
                            update.source.target_grade_count==prefix_tail.endpoint().source.grade_count() &&
                            update.removed_blocks.empty() && update.blocks.empty() &&
                            update.removed_carriers.empty() && update.carriers.empty() &&
                            update.removed_tensors.empty() && update.tensors.empty() &&
                            update.correspondences.empty() && update.removed_finite_systems.empty() &&
                            update.finite_systems.empty() && !update.codes && !update.source_migration,
                            "biomedical annotation changed its bound mathematical source");
                    }
                    const auto series=prefix.history();
                    const auto hash=sha256_chunks([&](const auto &sink){stream_native_series(series,sink);});
                    require(native_digest=="sha256:"+hash,"embedded native input differs from its integrity binding");
                    return series;
                }
                const auto b=publication->at("history").bytes();
                require(native_digest=="sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(b.data()),b.size())),"original native input differs from its compatibility integrity binding");
                return decode_native_series(b,native_data_limits());
            }
            const auto native_path = local_artifact_path(root, native);
            require(fs::file_size(native_path) <= source_byte_limit(), "input size budget exceeded");
            const rfc::storage::ReadOnlyMappedFile archive(native_path);
            require(archive.size() <= source_byte_limit(), "input size budget exceeded");
            const std::string_view bytes(reinterpret_cast<const char *>(archive.data()), archive.size());
            require(native_digest == "sha256:" + sha256(bytes),
                    "artifact SHA-256 mismatch: " + native.at("path").text());
            if(!generated)original_history.assign(archive.bytes().begin(),archive.bytes().end());
            return decode_native_series(archive.bytes(), source_codec_limits(archive.size()));
        }();
        require(series.source_namespace() == key(native.at("source_namespace")),
                "native namespace mismatch");
        require(native.has("cutoff") || series.segment_count() == 0,
                "a native history requires an explicit publication cutoff");
        auto cutoff = series.first_sequence();
        if (native.has("cutoff")) {
            const auto s = native.at("cutoff").text();
            const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), cutoff);
            require(!s.empty() && (s.size() == 1 || s.front() != '0') && error == std::errc{} &&
                        end == s.data() + s.size(),
                    "native cutoff must be a canonical unsigned decimal string");
        }
        history =
            std::make_shared<const NativeSeriesPrefix>(NativeSeriesPrefix::capture(series, cutoff));
        if(generated)producer_history=history->history();
    }
    const auto tail = history->tail();
    // The history already owns this immutable endpoint. Application admission
    // supplies the execution directory below; making another directory here
    // would copy every field only to replace it with the publication's owner.
    const auto endpoint = tail.pin();
    const auto docs = declaration.at("documents");
    source_array(docs, "documents");
    std::size_t document_bytes = 0;
    RecordByteBindings raw;
    const auto retained_source_keys=NativeRecord::source_keys(*endpoint);
    std::set<std::uint64_t> used_source_keys(retained_source_keys.begin(),retained_source_keys.end());
    std::uint64_t next_source_key=1;
    const auto allocate_source_key=[&] {
        while(used_source_keys.contains(next_source_key)) {
            require(next_source_key!=std::numeric_limits<std::uint64_t>::max(),"raw source identity space exhausted");
            ++next_source_key;
        }
        const auto key=next_source_key;used_source_keys.insert(key);return key;
    };
    raw.sources.emplace_back(allocate_source_key(),std::as_bytes(std::span(bytes.data(),bytes.size())));
    for (std::size_t i = 0; i < docs.size(); ++i) {
        const auto doc = docs.at(i);
        members(doc, {"id", "path", "sha256", "format"});
        const auto id = text(doc.at("id"));
        const auto format = doc.at("format").text();
        require(format == "utf8" || format == "opaque", "unsupported document format");
        auto content=publication?[&] {
            const auto item=publication->at("documents").at(i);
            require(item.at("id").text()==id,"native document identity/order mismatch");
            const auto b=item.at("bytes").bytes();
            return std::string(reinterpret_cast<const char *>(b.data()),b.size());
        }():local_artifact(root, doc);
        require(doc.at("sha256").text()=="sha256:"+sha256(content),"native evidence bytes differ from source binding");
        const auto source_key=allocate_source_key();
        raw.documents.emplace(id,source_key);
        raw.sources.emplace_back(source_key,std::as_bytes(std::span(content.data(),content.size())));
        document_bytes = rfc::machine::checked_size_add(document_bytes, content.size());
        if (format == "utf8") {
            require(content.find('\0') == std::string::npos, "NUL in UTF-8 evidence");
            (void)Json::parse(Json::string(content).dump());
        }
        const auto owned=raw.sources.back().bytes();
        require(documents.emplace(id, Document{doc, std::string_view(reinterpret_cast<const char *>(owned.data()),owned.size())}).second,
                "duplicate document identity");
    }
    // Biomedical meanings and original evidence spans now inhabit native
    // declaration channels on the actual retained relational complex.
    const auto make_publication=[&] {
        using N=NativeRecordValue;
        RecordByteBindings bound=raw;
        const auto history_key=generated?0:allocate_source_key();
        if(!generated)bound.sources.emplace_back(history_key,original_history);
        N::Array document_values,publication_values;
        for(std::size_t i=0;i<docs.size();++i)
            document_values.elements.emplace_back(N::Object{{{"id",N(docs.at(i).at("id").text())},
                {"bytes",N(N::ByteSpan{raw.sources.at(i+1).key(),0,raw.sources.at(i+1).bytes().size(),false})}}});
        std::map<const ScientificStudy *,std::uint64_t> owners;
        for(const auto &[id,owner]:publications) {
            auto found=owners.find(owner.get());
            if(found==owners.end()) {
                const auto key=allocate_source_key();
                bound.sources.push_back(owner->native_publication().archive_source(key));
                found=owners.emplace(owner.get(),key).first;
            }
            const auto source=std::find_if(bound.sources.begin(),bound.sources.end(),[&](const auto &s){return s.key()==found->second;});
            publication_values.elements.emplace_back(N::Object{{{"id",N(id)},
                {"source_key",N(rfc::numeric::Integer(std::to_string(found->second)))},
                {"bytes",N(N::ByteSpan{found->second,0,source->bytes().size(),false})}}});
        }
        N meanings(N::Object{});
        auto &meaning_members=std::get<N::Object>(meanings.data).members;
        json_object_object_foreach(declaration.raw(),name,entry) {
            (void)entry;
            const auto array=declaration.at(name);
            if(std::string_view(name)!="cells" && std::string_view(name)!="occurrences") {
                meaning_members.emplace_back(name,import_record_view(array,bound));
                continue;
            }
            const bool cells=std::string_view(name)=="cells";
            meaning_members.emplace_back(name,N(N::PackedArray{array.size(),
                [bindings=bound,array,cells,grade_count=endpoint->source.grade_count()](const auto &emit) {
                for(std::size_t i=0;i<array.size();++i) {
                    const auto row=array.at(i);
                    const auto grade=rfc::core::Grade{static_cast<std::uint16_t>(count(row.at("grade"),grade_count-1,"source grade"))};
                    N::Owned owned;
                    owned.payload=std::make_shared<const N>(import_record_view(row,bindings));
                    owned.packed=true;
                    if(cells)owned.cell=StableCellKey{grade,key(row.at("key"))};
                    else owned.occurrence=StableOccurrenceKey{grade,rfc::core::OccurrenceKey{key(row.at("key"))}};
                    emit(N(std::move(owned)));
                }
            }}));
        }
        N history_value=generated
            ?N(N::Object{{{"kind",N(std::string("publication_prefix"))},
                {"cutoff",N(rfc::numeric::Integer(std::to_string(history->cutoff())))}}})
            :N(N::ByteSpan{history_key,0,original_history.size(),false});
        N body(N::Object{{{"schema",N(std::string(generated?"gadget.biomedical.native.v2":"gadget.biomedical.native.v1"))},
            {"study",std::move(meanings)},
            {"declaration",N(N::ByteSpan{raw.sources.front().key(),0,bytes.size(),false})},
            {"history",std::move(history_value)},
            {"documents",N(std::move(document_values))},
            {"publications",N(std::move(publication_values))}}});
        if(generated)return NativeRecord::annotate(body,*producer_history,std::move(bound.sources));
        return NativeRecord::create(body,history->source_namespace(),std::move(bound.sources),endpoint.get());
    };
    native_publication=publication?publication:std::optional<NativeRecord>(make_publication());
    this->declaration=native_publication->at("study");
    // Numerical execution belongs to the explicitly bound source cutoff.
    // Annotation publication order is a separate, later native endpoint.
    directory=(generated || embedded)
        ?std::make_shared<const NativeAccessionDirectory>(tail,
            std::vector<RelationalCodeObservationQuery>{},
            [](auto value){return rfcx_detail::fnv1a64(value);},native_accession_limits())
        :native_publication->accession_directory();
    for(std::size_t i=0;i<docs.size();++i) {
        auto &document=documents.at(docs.at(i).at("id").text());
        document.metadata=this->declaration.member("documents").element(i);
        const auto b=native_publication->at("documents").at(i).at("bytes").bytes();
        document.bytes=std::string_view(reinterpret_cast<const char *>(b.data()),b.size());
    }
    for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k)
        grade_pins.push_back(directory->pin(NativeAccessionRequest::sparse_chain(
            rfc::core::Grade{static_cast<std::uint16_t>(k)}, {}, {})));
    const auto design = declaration.at("design");
    members(design, {"question", "population", "selection", "covariates", "limitations"});
    for (const auto &name : {"question", "population", "selection"})
        (void)text(design.at(name));
    source_array(design.at("limitations"), "limitations");
    for (std::size_t i = 0; i < design.at("limitations").size(); ++i)
        (void)text(design.at("limitations").at(i));
    source_array(design.at("covariates"), "covariates");
    std::set<std::string> covariates;
    for (std::size_t i = 0; i < design.at("covariates").size(); ++i) {
        auto item = design.at("covariates").at(i);
        members(item, {"name", "description", "field", "evidence"});
        require(covariates.insert(text(item.at("name"))).second, "duplicate covariate");
        (void)text(item.at("description"));
        (void)evidence(item.at("evidence"));
    }
    const auto declared_cells = declaration.at("cells");
    source_array(declared_cells, "cells");
    const StableSourceIndex index(endpoint->source);
    std::size_t native_cells = 0;
    for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k) {
        const auto size =
            endpoint->source.grade(rfc::core::Grade{static_cast<std::uint16_t>(k)}).cell_count();
        native_cells = rfc::machine::checked_size_add(native_cells, size);
        require(native_cells <= declared_cells.size(), "incomplete native cell evidence map");
        coordinate_cells.emplace_back(size);
    }
    for (std::size_t i = 0; i < declared_cells.size(); ++i) {
        auto item = declared_cells.at(i);
        members(item, {"grade", "key", "id", "kind", "label", "evidence"});
        const auto cell =
            StableCellKey{rfc::core::Grade{static_cast<std::uint16_t>(count(
                              item.at("grade"), endpoint->source.grade_count() - 1, "grade"))},
                          key(item.at("key"))};
        const auto coordinate = index.cell(cell);
        require(coordinate.has_value(), "semantic cell absent from native source");
        require(cells.emplace(cell, this->declaration.member("cells").element(i)).second &&
                    cell_ids.emplace(text(item.at("id")), cell).second,
                "duplicate cell identity");
        coordinate_cells.at(coordinate->grade.value).at(coordinate->index) = cell;
        (void)text(item.at("kind"));
        (void)text(item.at("label"));
        (void)evidence(item.at("evidence"));
    }
    const auto declared_occurrences = declaration.at("occurrences");
    source_array(declared_occurrences, "occurrences");
    for (std::size_t i = 0; i < declared_occurrences.size(); ++i) {
        auto item = declared_occurrences.at(i);
        members(item, {"grade", "key", "evidence"});
        const auto occurrence = StableOccurrenceKey{
            rfc::core::Grade{static_cast<std::uint16_t>(
                count(item.at("grade"), endpoint->source.grade_count() - 1, "grade"))},
            rfc::core::OccurrenceKey{key(item.at("key"))}};
        require(index.occurrence(occurrence).has_value(),
                "semantic occurrence absent from native source");
        require(occurrences.emplace(occurrence, this->declaration.member("occurrences").element(i)).second, "duplicate occurrence identity");
        (void)evidence(item.at("evidence"));
    }
    std::size_t cell_count = 0, occurrence_count = 0;
    for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k) {
        const auto &arena = endpoint->source.grade(rfc::core::Grade{static_cast<std::uint16_t>(k)});
        cell_count = rfc::machine::checked_size_add(cell_count, arena.cell_count());
        occurrence_count =
            rfc::machine::checked_size_add(occurrence_count, arena.incidence_count());
    }
    require(cells.size() == cell_count && occurrences.size() == occurrence_count,
            "semantic map must cover every retained cell and occurrence");
    const auto declared_contexts = declaration.at("contexts");
    source_array(declared_contexts, "contexts");
    require(declared_contexts.size() > 0, "at least one declared context required");
    for (std::size_t i = 0; i < declared_contexts.size(); ++i) {
        auto item = declared_contexts.at(i);
        members(item, {"id", "block", "metrics", "time"}, {"witness_scope"});
        if (item.has("witness_scope"))
            require(item.at("witness_scope").text() == "block" && item.at("block").raw(),
                    "block witness scope requires a declared strict subcomplex");
        Context context{this->declaration.member("contexts").element(i), {}};
        if (!json_object_is_type(item.at("block").raw(), json_type_null)) {
            context.block = key(item.at("block"));
            require(std::any_of(endpoint->local_blocks.begin(), endpoint->local_blocks.end(),
                                [&](const auto &b) {
                                    return b.block_key == *context.block &&
                                           b.kind == LocalBlockKind::StrictSubcomplex;
                                }),
                    "context requires a declared strict subcomplex block");
        }
        const auto metrics = item.at("metrics");
        require(metrics.size() == endpoint->source.grade_count(),
                "one explicit metric per native grade required");
        for (std::size_t k = 0; k < metrics.size(); ++k) {
            auto metric = metrics.at(k);
            const auto kind = metric.at("kind").text();
            if (kind == "identity")
                members(metric, {"kind"});
            else {
                require(kind == "diagonal", "only declared identity/diagonal metrics admitted");
                members(metric, {"kind", "values"});
                auto values = metric.at("values");
                require(values.size() ==
                            endpoint->source.grade(rfc::core::Grade{static_cast<std::uint16_t>(k)})
                                .cell_count(),
                        "ambient metric dimension mismatch");
                for (std::size_t j = 0; j < values.size(); ++j)
                    require(exact_decimal(values.at(j).text()) > 0, "nonpositive metric");
            }
        }
        auto time = item.at("time");
        if (!json_object_is_type(time.raw(), json_type_null)) {
            members(time, {"axis", "position", "unit", "population"});
            for (const auto &name : {"axis", "unit", "population"})
                (void)text(time.at(name));
            (void)exact_decimal(time.at("position").text());
        }
        require(contexts.emplace(text(item.at("id")), std::move(context)).second,
                "duplicate context");
    }
    const auto tensors = temporal_endpoint_detail::tensors(*endpoint);
    const auto carriers = temporal_endpoint_detail::carriers(*endpoint);
    const auto declared_fields = declaration.at("fields");
    source_array(declared_fields, "fields");
    for (std::size_t i = 0; i < declared_fields.size(); ++i) {
        auto item = declared_fields.at(i);
        members(item, {"id", "key", "context", "channel", "quantity", "units", "method",
                       "value_type", "evidence"});
        const auto tensor = key(item.at("key"));
        require(tensors.contains(tensor) && tensors.at(tensor)->role == TemporalTensorRole::Field &&
                    tensors.at(tensor)->output,
                "field must name a native field tensor");
        require(contexts.contains(text(item.at("context"))), "field context unknown");
        const auto channel = item.at("channel").text();
        require(channel == "G" || channel == "E" || channel == "T" || channel == "covariate" ||
                    channel == "observation",
                "unknown independent evidence channel");
        for (const auto &name : {"quantity", "units", "method"})
            (void)text(item.at(name));
        const auto value_type = item.at("value_type").text();
        require(value_type == "quantitative" || value_type == "indicator" ||
                    value_type == "categorical",
                "unknown field value type");
        if (value_type == "indicator")
            for (const auto &value : tensors.at(tensor)->tensor.values().entries())
                require(value.value == 0 || value.value == 1,
                        "indicator field must contain exact zero/one values");
        (void)evidence(item.at("evidence"));
        require(carriers.contains(*tensors.at(tensor)->output), "field carrier missing");
        if (!field_pins.contains(tensor)) {
            NativeAccessionRequest accession;
            accession.kind = NativeAccessionKind::Field;
            accession.subject = tensor;
            field_pins.emplace(tensor, directory->pin(accession));
        }
        require(fields.emplace(text(item.at("id")), Field{this->declaration.member("fields").element(i), tensor}).second,
                "duplicate field identity");
    }
    const auto sets=declaration.at("field_sets");
    admit_scoped_actions(declaration.has("operators")?declaration.at("operators"):Json{},sets);
    for (std::size_t i = 0; i < design.at("covariates").size(); ++i) {
        auto field = design.at("covariates").at(i).at("field");
        if (!json_object_is_type(field.raw(), json_type_null))
            require(fields.contains(field.text()) &&
                        fields.at(field.text()).declaration.at("channel").text() == "covariate",
                    "covariate must name a declared covariate field");
    }
    digest = "sha256:" + sha256(schema + "\n" + bytes + "\n" + native_digest);
    if (declaration.has("publications")) {
        const auto records = declaration.at("publications");
        source_array(records, "study publications");
        std::map<fs::path, std::shared_ptr<const ScientificStudy>> loaded;
        std::map<std::uint64_t,std::shared_ptr<const ScientificStudy>> native_loaded;
        for (std::size_t i = 0; i < records.size(); ++i) {
            const auto record = records.at(i);
            members(record, {"id", "study", "source_digest", "meaning", "evidence"});
            const auto id = text(record.at("id"), source_byte_limit());
            require(id != "self" && !publications.contains(id),
                    "duplicate/reserved publication id");
            (void)text(record.at("meaning"), source_byte_limit());
            (void)evidence(record.at("evidence"));
            const auto artifact = record.at("study");
            members(artifact, {"path", "sha256"});
            std::shared_ptr<const ScientificStudy> owner;
            if(publication) {
                const auto entries=publication->at("publications");
                std::optional<NativeRecord> selected;
                std::uint64_t owner_key{};
                for(std::size_t e=0;e<entries.size();++e)if(entries.at(e).at("id").text()==id) {
                    owner_key=std::stoull(entries.at(e).at("source_key").integer().get_str());
                    if(native_loaded.contains(owner_key)) {owner=native_loaded.at(owner_key);break;}
                    const auto b=entries.at(e).at("bytes").bytes();
                    const auto series=decode_native_series(b,native_data_limits());
                    selected=NativeRecord::from_history(series);
                    break;
                }
                require(owner || selected.has_value(),"referenced native publication bytes are missing");
                if(!owner) {owner=std::make_shared<const ScientificStudy>(std::move(*selected));native_loaded.emplace(owner_key,owner);}
            } else {
                const auto retained = local_artifact(root, artifact);
                const auto file = fs::canonical(root / artifact.at("path").text());
                if (!loaded.contains(file))loaded.emplace(file,std::make_shared<const ScientificStudy>(file));
                owner=loaded.at(file);
                require(owner->state_->declaration_digest == "sha256:" + sha256(retained),"publication declaration changed during admission");
            }
            require(owner->digest() == record.at("source_digest").text(),
                    "publication admitted source digest mismatch");
            publications.emplace(id, owner);
            publication_declarations.emplace(id, this->declaration.member("publications").element(i));
        }
    }
    if(!publication && !publications.empty()) {
        native_publication=make_publication();
        this->declaration=native_publication->at("study");
        directory=generated
            ?std::make_shared<const NativeAccessionDirectory>(tail,
                std::vector<RelationalCodeObservationQuery>{},
                [](auto value){return rfcx_detail::fnv1a64(value);},native_accession_limits())
            :native_publication->accession_directory();
        grade_pins.clear();field_pins.clear();operator_pins.clear();
        for(std::size_t k=0;k<endpoint->source.grade_count();++k)
            grade_pins.push_back(directory->pin(NativeAccessionRequest::sparse_chain(rfc::core::Grade{static_cast<std::uint16_t>(k)}, {}, {})));
        for(std::size_t i=0;i<docs.size();++i) {
            auto &document=documents.at(docs.at(i).at("id").text());
            document.metadata=this->declaration.member("documents").element(i);
            const auto b=native_publication->at("documents").at(i).at("bytes").bytes();
            document.bytes=std::string_view(reinterpret_cast<const char *>(b.data()),b.size());
        }
        for(std::size_t i=0;i<declared_cells.size();++i) {
            const auto item=declared_cells.at(i);
            const auto cell=cell_ids.at(item.at("id").text());cells.at(cell)=this->declaration.member("cells").element(i);
        }
        for(std::size_t i=0;i<declared_occurrences.size();++i) {
            const auto item=declared_occurrences.at(i);
            const StableOccurrenceKey key{rfc::core::Grade{static_cast<std::uint16_t>(item.at("grade").integer())},rfc::core::OccurrenceKey{scientific_detail::key(item.at("key"))}};
            occurrences.at(key)=this->declaration.member("occurrences").element(i);
        }
        for(std::size_t i=0;i<declared_contexts.size();++i)
            contexts.at(declared_contexts.at(i).at("id").text()).declaration=this->declaration.member("contexts").element(i);
        for(std::size_t i=0;i<declared_fields.size();++i) {
            auto &field=fields.at(declared_fields.at(i).at("id").text());field.declaration=this->declaration.member("fields").element(i);
            NativeAccessionRequest request;request.kind=NativeAccessionKind::Field;request.subject=field.tensor;
            field_pins.emplace(field.tensor,directory->pin(request));
        }
        for(std::size_t i=0;i<sets.size();++i)field_sets.at(sets.at(i).at("id").text()).declaration=this->declaration.member("field_sets").element(i);
        if(declaration.has("operators"))for(std::size_t i=0;i<declaration.at("operators").size();++i) {
            const auto item=declaration.at("operators").at(i);const auto id=item.at("id").text();operators.at(id)=this->declaration.member("operators").element(i);
            NativeAccessionRequest request;request.subject=scientific_detail::key(item.at("key"));
            const auto role=tensors.at(request.subject)->role;
            request.kind=role==TemporalTensorRole::Metric || role==TemporalTensorRole::CrossForm?NativeAccessionKind::Tensor:NativeAccessionKind::Operator;
            operator_pins.emplace(id,directory->pin(request));
        }
        for(std::size_t i=0;i<declaration.at("publications").size();++i)
            publication_declarations.at(declaration.at("publications").at(i).at("id").text())=this->declaration.member("publications").element(i);
    }
}
void ScientificStudy::State::admit_scoped_actions(const Json &declared_operators,
                                                  const Json &sets) {
    const auto endpoint = directory->pin_snapshot();
    const auto tensors = temporal_endpoint_detail::tensors(*endpoint);
    if (declared_operators.raw()) {
        source_array(declared_operators, "operators");
        for (std::size_t i = 0; i < declared_operators.size(); ++i) {
            const auto item = declared_operators.at(i);
            members(item, {"id", "key", "quantity", "units", "method", "evidence"});
            const auto id = text(item.at("id"));
            for (const auto *name : {"quantity", "units", "method"})
                (void)text(item.at(name));
            (void)evidence(item.at("evidence"));
            NativeAccessionRequest request;
            request.subject = key(item.at("key"));
            require(tensors.contains(request.subject), "operator tensor key unknown");
            const auto role = tensors.at(request.subject)->role;
            // A stored metric is a bilinear form. Its tensor accession does
            // not implicitly turn it into an endomorphism.
            request.kind =
                role == TemporalTensorRole::Metric || role == TemporalTensorRole::CrossForm
                    ? NativeAccessionKind::Tensor
                    : NativeAccessionKind::Operator;
            const auto pin = directory->pin(request);
            const auto tensor = pin.tensor_declaration();
            require(tensor->input && tensor->output && tensor->tensor.domain().rank() == 2,
                    "operator declaration requires an actual map with both carrier bindings");
            require(operators.emplace(id, this->declaration.member("operators").element(i)).second,
                    "duplicate operator identity");
            operator_pins.emplace(id, pin);
        }
    }
    source_array(sets, "field sets");
    for (std::size_t i = 0; i < sets.size(); ++i) {
        auto item = sets.at(i);
        members(item, {"id", "grade", "fields", "carriers", "maps"}, {"coordinate_system"});
        NativeFieldRequest request;
        request.upper_grade = rfc::core::Grade{static_cast<std::uint16_t>(
            count(item.at("grade"), endpoint->source.grade_count() - 1, "field set grade"))};
        source_array(item.at("fields"), "field set fields");
        for (std::size_t j = 0; j < item.at("fields").size(); ++j) {
            const auto name = item.at("fields").at(j).text();
            require(fields.contains(name), "unknown field set member");
            require(fields.at(name).declaration.at("value_type").text() != "categorical",
                    "categorical fields cannot be linear native amplitudes");
            require(!contexts.at(fields.at(name).declaration.at("context").text()).block,
                    "native field action requires the whole declared source context");
            NativeAccessionRequest accession;
            accession.kind = NativeAccessionKind::Field;
            accession.subject = fields.at(name).tensor;
            request.fields.push_back(directory->pin(accession));
        }
        source_array(item.at("carriers"), "field set carriers");
        for (std::size_t j = 0; j < item.at("carriers").size(); ++j)
            request.stalk_carriers.push_back(key(item.at("carriers").at(j)));
        NativeFieldRealization realization;
        const auto maps = item.at("maps").text();
        if (maps == "constant") {
            require(item.has("coordinate_system"),
                    "constant fiber requires a declared coordinate system");
            realization = NativeFieldRealization::constant_fiber(key(item.at("coordinate_system")));
        } else
            require(maps == "explicit" && !item.has("coordinate_system"),
                    "explicit maps cannot declare implicit identity transport");
        // Canonical admission owns all map roles, complete coordinates, support,
        // repeated/inactive occurrences, matching fibers and source identities.
        (void)NativeFieldPlan::inspect(
            directory->pin(NativeAccessionRequest::sparse_chain(request.upper_grade, {}, {})),
            request, realization);
        require(
            field_sets
                .emplace(text(item.at("id")),
                         FieldSet{this->declaration.member("field_sets").element(i), realization})
                .second,
            "duplicate field set");
    }
}
Json ScientificStudy::State::carrier_metadata(const TemporalCarrier &carrier) const {
    auto item = Json::object(), extents = Json::array();
    item.set("key", Json::string(std::to_string(carrier.key)));
    item.set("kind", Json::string(carrier.kind == TemporalCarrierKind::SourceGrade ? "source_grade"
                                  : carrier.kind == TemporalCarrierKind::CellStalk ? "cell_stalk"
                                                                                   : "explicit"));
    item.set("grade", Json::integer(carrier.grade.value));
    item.set("coordinate_system", Json::string(std::to_string(carrier.coordinate_system)));
    for (const auto extent : carrier.domain.extents())
        extents.push(cardinality(extent));
    item.set("extents", extents);
    item.set("coordinates", cardinality(carrier.domain.logical_size()));
    if (carrier.kind == TemporalCarrierKind::CellStalk) {
        auto cell = Json::object();
        cell.set("grade", Json::integer(carrier.cell->grade.value));
        cell.set("key", Json::string(std::to_string(carrier.cell->key)));
        item.set("cell", cell);
    }
    return item;
}
Json ScientificStudy::State::field_metadata(const Field &field) const {
    auto item = Json::object();
    for (const auto &member :
         {"id", "key", "context", "channel", "quantity", "units", "method", "value_type"})
        item.set(member, field.declaration.at(member));
    const auto pin = field_pins.at(field.tensor);
    const auto &tensor = *pin.tensor_declaration();
    const auto &carrier = *pin.carrier();
    item.set("grade", Json::integer(carrier.grade.value));
    item.set("carrier", Json::string(std::to_string(carrier.key)));
    item.set("coordinate_system", Json::string(std::to_string(carrier.coordinate_system)));
    auto axes = Json::array();
    for (const auto extent : tensor.tensor.domain().extents())
        axes.push(cardinality(extent));
    item.set("extents", axes);
    item.set("coordinates", cardinality(tensor.tensor.domain().logical_size()));
    item.set("reported", cardinality(tensor.tensor.support().count()));
    item.set("profile",
             Json::string(carrier.kind == TemporalCarrierKind::SourceGrade ? "source_grade"
                          : carrier.kind == TemporalCarrierKind::Explicit  ? "explicit_coordinates"
                                                                           : "native_stalk"));
    return item;
}
Json ScientificStudy::State::operator_metadata(const std::string &id) const {
    auto item = Json::parse(operators.at(id).dump());
    const auto t = operator_pins.at(id).tensor_declaration();
    item.set("native_role", Json::integer(static_cast<std::int64_t>(t->role)));
    item.set("input_carrier", Json::string(std::to_string(*t->input)));
    item.set("output_carrier", Json::string(std::to_string(*t->output)));
    auto extents = Json::array();
    for (auto extent : t->tensor.domain().extents())
        extents.push(cardinality(extent));
    item.set("extents", extents);
    item.set("full_support",
             Json::boolean(t->tensor.support().count() == t->tensor.domain().logical_size()));
    return item;
}
Json ScientificStudy::State::history_metadata() const {
    auto binding = Json::object();
    binding.set("archive_sha256", Json::string(native_digest));
    binding.set("first_sequence", Json::string(std::to_string(history->first_sequence())));
    binding.set("cutoff", Json::string(std::to_string(history->cutoff())));
    binding.set("retained_segments", cardinality(history->segment_count()));
    return binding;
}
NativeRecordValue ScientificStudy::State::catalog_value(std::vector<RawByteSource> &sources) const {
    using N = NativeRecordValue;
    const auto count_value = [](std::size_t value) { return N(rfc::numeric::Integer(std::to_string(value))); };
    const auto text_value = [](std::string value) { return N(std::move(value)); };
    const auto put = [](N &value, std::string name, N member) {
        auto &members = std::get<N::Object>(value.data).members;
        for (auto &[key, current] : members)
            if (key == name) { current = std::move(member); return; }
        members.emplace_back(std::move(name), std::move(member));
    };
    std::map<std::uint64_t, std::uint64_t> byte_map;
    const auto retain_value = [&](const NativeRecord &value) {
        for (const auto &raw : value.raw_sources()) {
            if (byte_map.contains(raw.key())) continue;
            const auto found = std::find_if(sources.begin(), sources.end(),
                [&](const auto &candidate) { return candidate.shares_owner(raw); });
            if (found != sources.end()) byte_map.emplace(raw.key(), found->key());
            else {
                std::uint64_t key = 1;
                for (const auto &candidate : sources)
                    key = std::max<std::uint64_t>(key, rfc::machine::checked_size_add(candidate.key(), 1));
                sources.push_back(raw.rekey(key)); byte_map.emplace(raw.key(), key);
            }
        }
        return value.unbound_value(byte_map);
    };
    const auto declared = [&](std::string_view name) { return retain_value(declaration.native().at(name)); };
    const auto endpoint = directory->pin_snapshot();
    N result(N::Object{});
    put(result, "schema", text_value("gadget.scientific.catalog.v1"));
    put(result, "study_id", declared("id"));
    put(result, "source_digest", text_value(digest));
    put(result, "source_sequence", text_value(std::to_string(endpoint->sequence)));
    N carriers(N::PackedArray{endpoint->carriers.size(), [endpoint, put, count_value, text_value](const auto &emit) {
      for (const auto &carrier : endpoint->carriers) {
        N item(N::Object{});
        put(item, "key", text_value(std::to_string(carrier.key)));
        put(item, "kind", text_value(carrier.kind == TemporalCarrierKind::SourceGrade ? "source_grade"
            : carrier.kind == TemporalCarrierKind::CellStalk ? "cell_stalk" : "explicit"));
        put(item, "grade", count_value(carrier.grade.value));
        put(item, "coordinate_system", text_value(std::to_string(carrier.coordinate_system)));
        N::Array extents;
        for (const auto extent : carrier.domain.extents()) extents.elements.push_back(count_value(extent));
        put(item, "extents", N(std::move(extents)));
        put(item, "coordinates", count_value(carrier.domain.logical_size()));
        if (carrier.kind == TemporalCarrierKind::CellStalk)
            put(item, "cell", N(N::Object{{{"grade", count_value(carrier.cell->grade.value)},
                {"key", text_value(std::to_string(carrier.cell->key))}}}));
        emit(std::move(item));
      }
    }});
    put(result, "carriers", std::move(carriers));
    if (declaration.at("native").has("cutoff"))
        put(result, "native_history", N(N::Object{{{"archive_sha256", text_value(native_digest)},
            {"first_sequence", text_value(std::to_string(history->first_sequence()))},
            {"cutoff", text_value(std::to_string(history->cutoff()))},
            {"retained_segments", count_value(history->segment_count())}}}));
    put(result, "design", declared("design"));
    N::Array bindings;
    for (const auto &[id, owner] : publications) {
        auto binding = retain_value(publication_declarations.at(id).native());
        put(binding, "catalog", owner->state_->catalog_value(sources));
        bindings.elements.push_back(std::move(binding));
    }
    put(result, "publications", N(std::move(bindings)));
    put(result, "contexts", declared("contexts"));
    N::Array list;
    for (const auto &[name, field] : fields) {
        (void)name;
        N item(N::Object{});
        for (const auto *member : {"id", "key", "context", "channel", "quantity", "units", "method", "value_type"})
            put(item, member, retain_value(field.declaration.native().at(member)));
        const auto pin = field_pins.at(field.tensor);
        const auto &tensor = *pin.tensor_declaration(); const auto &carrier = *pin.carrier();
        put(item, "grade", count_value(carrier.grade.value));
        put(item, "carrier", text_value(std::to_string(carrier.key)));
        put(item, "coordinate_system", text_value(std::to_string(carrier.coordinate_system)));
        N::Array extents;
        for (const auto extent : tensor.tensor.domain().extents()) extents.elements.push_back(count_value(extent));
        put(item, "extents", N(std::move(extents)));
        put(item, "coordinates", count_value(tensor.tensor.domain().logical_size()));
        put(item, "reported", count_value(tensor.tensor.support().count()));
        put(item, "profile", text_value(carrier.kind == TemporalCarrierKind::SourceGrade ? "source_grade"
            : carrier.kind == TemporalCarrierKind::Explicit ? "explicit_coordinates" : "native_stalk"));
        list.elements.push_back(std::move(item));
    }
    put(result, "fields", N(std::move(list)));
    if (declaration.has("operators")) {
        N::Array list_operators;
        for (const auto &[id, declaration_] : operators) {
            auto item = retain_value(declaration_.native());
            const auto tensor = operator_pins.at(id).tensor_declaration();
            put(item, "native_role", count_value(static_cast<std::size_t>(tensor->role)));
            put(item, "input_carrier", text_value(std::to_string(*tensor->input)));
            put(item, "output_carrier", text_value(std::to_string(*tensor->output)));
            N::Array extents;
            for (const auto extent : tensor->tensor.domain().extents()) extents.elements.push_back(count_value(extent));
            put(item, "extents", N(std::move(extents)));
            put(item, "full_support", N(tensor->tensor.support().count() == tensor->tensor.domain().logical_size()));
            list_operators.elements.push_back(std::move(item));
        }
        put(result, "operators", N(std::move(list_operators)));
    }
    put(result, "field_sets", declared("field_sets"));
    N::Array operations;
    for (const auto *operation : {"inspect", "coverage", "provenance", "boundary", "adjoint", "hodge",
        "typed_moment", "temporal_moment_change", "native_boundary", "field_summary", "support_mask",
        "controlled_pairings", "trajectory_structure", "structure_transition", "paired_field_contrasts",
        "local_sections", "tensor_window", "green_response", "chain_certificate", "realized_path_comparison",
        "green_action", "linear_model_comparison", "finite_reference_calibration", "relational_field_ensemble",
        "structural_reference_calibration", "source_structure_ablation", "observed_response_model",
        "constrained_response_model", "temporal_field_history", "censored_endpoint_bounds",
        "declared_response_transition", "hodge_frame_decomposition", "linear_response_observation", "metric_relative_novelty"})
        operations.elements.push_back(text_value(operation));
    put(result, "operations", N(std::move(operations)));
    N::Array docs;
    for (const auto &[name, document] : documents) {
        (void)name;
        auto item = retain_value(document.metadata.native());
        if (document.metadata.at("format").text() == "utf8") {
            auto end = std::min<std::size_t>(document.bytes.size(), 2048);
            while (end && !boundary(document.bytes, end)) --end;
            put(item, "text_prefix", text_value(std::string(document.bytes.substr(0, end))));
            put(item, "prefix_bytes", count_value(end));
        }
        docs.elements.push_back(std::move(item));
    }
    put(result, "documents", N(std::move(docs)));
    auto identity_records = std::make_shared<const std::map<StableCellKey, RetainedRecord>>(cells);
    N identities(N::PackedArray{identity_records->size(), [identity_records, count_value, text_value, put](const auto &emit) {
        for (const auto &[cell, item] : *identity_records) {
            const auto original = item.native();
            N identity(N::Object{});
            put(identity, "grade", count_value(cell.grade.value));
            for (const auto *member : {"key", "id", "kind", "label"})
                put(identity, member, text_value(original.at(member).text()));
            emit(std::move(identity));
        }
    }});
    put(result, "identities", std::move(identities));
    put(result, "identity_catalog_complete", N(!compact));
    return result;
}
Json ScientificStudy::State::catalog() const {
    auto sources = native_publication->raw_sources();
    auto value = catalog_value(sources);
    return render_record_view(NativeRecord::create(std::move(value), new_native_record_namespace(),
                                                   std::move(sources)));
}
ScientificStudy::ScientificStudy(const fs::path &path) : ScientificStudy(std::make_unique<State>(path)) {}
ScientificStudy::ScientificStudy(const fs::path &path, const ScientificAdmissionObserver &observer)
    : ScientificStudy(std::make_unique<State>(path, std::nullopt, std::nullopt, observer)) {
    if (observer) observer("study_complete");
}
std::shared_ptr<const ScientificStudy> ScientificStudy::admit_generated(const fs::path &path,NativeTemporalSeries history) {
    return std::shared_ptr<const ScientificStudy>(new ScientificStudy(std::make_unique<State>(path,std::nullopt,std::move(history))));
}
ScientificStudy::ScientificStudy(NativeRecord publication)
    :ScientificStudy(std::make_unique<State>(fs::path("native-publication-"+std::to_string(new_native_record_namespace())),std::move(publication))) {}
ScientificStudy::ScientificStudy(NativeRecord publication, const ScientificAdmissionObserver &observer)
    : ScientificStudy(std::make_unique<State>(fs::path("native-publication-" +
          std::to_string(new_native_record_namespace())), std::move(publication), std::nullopt, observer)) {
    if (observer) observer("study_complete");
}
ScientificStudy::ScientificStudy(std::unique_ptr<State> state) :state_(std::move(state)) {
    auto index = [&](const std::string &section, const auto &items) {
        auto &keys = state_->catalog_keys[section];
        keys.reserve(items.size());
        for (const auto &[id, item] : items) {
            (void)item;
            keys.push_back(id);
        }
    };
    index("fields", state_->fields);
    index("contexts", state_->contexts);
    index("operators", state_->operators);
    index("field_sets", state_->field_sets);
    index("documents", state_->documents);
    index("publications", state_->publications);
    for (const auto &[cell, item] : state_->cells) {
        (void)item;
        state_->cell_keys.push_back(cell);
    }
    for (const auto &[occurrence, item] : state_->occurrences) {
        (void)item;
        state_->occurrence_keys.push_back(occurrence);
    }
}
ScientificStudy::~ScientificStudy() = default;
const NativeRecord &ScientificStudy::native_publication() const {return state_->native_publication.value();}
void ScientificStudy::write_native(const fs::path &path) const {native_publication().write(path);}
rfc::storage::PinnedNativeAccession ScientificStudy::field_accession(const std::string &id) const {
    return state_->field_pins.at(state_->fields.at(id).tensor);
}
rfc::storage::NativeAtlasSource ScientificStudy::native_atlas_source(std::uint64_t key) const {
    return rfc::storage::NativeAtlasSource(key, *state_->directory);
}
Json ScientificStudy::catalog() const { return Json::parse(state_->catalog().dump()); }
Json ScientificStudy::declaration() const { return Json::parse(state_->declaration.dump()); }
rfc::storage::StableCellKey ScientificStudy::cell_identity(rfc::core::Grade grade,
                                                         const std::string &id) const {
    const auto cell=state_->resolve_cell(id);
    require(cell.grade==grade,"source cell identity has a different grade");
    return cell;
}
Json ScientificStudy::field_declaration(const std::string &id) const {
    require(state_->fields.contains(id), "source field identity is absent");
    const auto &source = state_->fields.at(id).declaration;
    auto result = Json::object();
    for (const auto *name : {"id", "key", "value_type", "quantity", "units", "channel"})
        result.set(name, source.at(name));
    return Json::parse(result.dump());
}
Json ScientificStudy::source_binding() const {
    auto binding = Json::object();
    binding.set("native", Json::parse(state_->declaration.at("native").dump()));
    binding.set("declaration_sha256", Json::string(state_->declaration_digest));
    return binding;
}
std::string ScientificStudy::digest() const { return state_->digest; }
std::shared_ptr<const rfc::storage::NativeSeriesPrefix> ScientificStudy::native_history() const {
    return state_->history;
}
std::shared_ptr<const rfc::storage::TemporalEndpoint> ScientificStudy::endpoint() const {
    return state_->directory->pin_snapshot();
}
} // namespace gadget
