#include <algorithm>
#include <gadget/snapshot.hpp>
#include <rfc/core/axioms.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <set>
#include <stdexcept>

namespace gadget {
namespace {
using namespace rfc::storage;
constexpr std::uint64_t study_namespace = 0x4741444745540001ULL;
constexpr std::size_t artifact_limit = std::numeric_limits<std::size_t>::max();
void require(bool valid, const std::string &message) {
    if (!valid)
        throw std::invalid_argument("snapshot: " + message);
}
void members(const Json &j, const std::set<std::string> &keys) {
    require(json_object_is_type(j.raw(), json_type_object), "expected object");
    require(static_cast<std::size_t>(json_object_object_length(j.raw())) == keys.size(),
            "missing or extra object members");
    json_object_object_foreach(j.raw(), key, value) {
        (void)value;
        require(keys.contains(key), std::string("unknown member: ") + key);
    }
}
template <class F> void entries(const Json &j, F &&f) {
    require(json_object_is_type(j.raw(), json_type_object), "expected object map");
    json_object_object_foreach(j.raw(), key, value) {
        (void)value;
        f(std::string(key), j.at(key));
    }
}
std::size_t count(const Json &j, std::size_t limit) {
    auto n = j.integer();
    require(n >= 0 && static_cast<std::uint64_t>(n) <= limit, "integer out of range");
    return static_cast<std::size_t>(n);
}
bool digest_valid(const std::string &s) {
    return s.starts_with("sha256:") && s.size() == 71 &&
           std::all_of(s.begin() + 7, s.end(),
                       [](char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); });
}
std::string compact(const Json &j) {
    return json_object_to_json_string_ext(j.raw(), JSON_C_TO_STRING_NOSLASHESCAPE);
}
std::span<const std::byte> bytes(const std::string &s) {
    return {reinterpret_cast<const std::byte *>(s.data()), s.size()};
}
std::string digest(std::span<const std::byte> b) {
    return "sha256:" + sha256(std::string(reinterpret_cast<const char *>(b.data()), b.size()));
}
std::string study_digest(const std::string &native, const std::string &metadata) {
    return "sha256:" + sha256("gadget.snapshot.v1\n" + native + "\n" + metadata + "\n");
}
Json completion_view(const Json &metadata,std::span<const std::byte> native,const Json &creator) {
    auto completion=Json::object(),artifacts=Json::object();const auto view=compact(metadata)+"\n";
    for(const auto &[name,data]:std::map<std::string,std::span<const std::byte>>{{"study.rfcxs",native},{"metadata.json",bytes(view)}}) {
        auto record=Json::object();record.set("bytes",Json::integer(data.size()));record.set("sha256",Json::string(digest(data)));artifacts.set(name,record);
    }
    completion.set("schema",Json::string("gadget.snapshot.v1"));completion.set("profile",Json::string("pager_get.v1"));
    completion.set("codec",Json::string("rfc.native_series.v1"));completion.set("source_digest",metadata.at("source_manifest").at("source_digest"));
    completion.set("study_digest",Json::string(study_digest(artifacts.at("study.rfcxs").at("sha256").text(),artifacts.at("metadata.json").at("sha256").text())));
    completion.set("artifacts",artifacts);completion.set("rexgraph_headers_sha256",creator.at("rexgraph_headers_sha256"));completion.set("gadget_sources_sha256",creator.at("gadget_sources_sha256"));return completion;
}
std::vector<std::byte> artifact(const fs::path &root, const std::string &name, std::size_t limit) {
    auto path = root / name;
    require(fs::is_regular_file(fs::symlink_status(path)), "artifact must be a regular file");
    return read_file_exact(path, limit);
}
Json json_bytes(const std::vector<std::byte> &b) {
    return Json::parse(std::string(reinterpret_cast<const char *>(b.data()), b.size()));
}
template<class Record> Json json_map(const std::map<std::string, Record> &values) {
    auto j = Json::object();
    for (const auto &[name, value] : values)
        j.set(name, value);
    return j;
}
rfc::numeric::SparseRationalVector stored_coordinates(const Vector &values) {
    std::vector<rfc::numeric::SparseRationalEntry> nonzero;
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] != 0)
            nonzero.push_back({i, values[i]});
    return rfc::numeric::SparseRationalVector(values.size(), std::move(nonzero));
}
// Admission for this biomedical profile. General source serialization remains
// the native library's responsibility; unsupported study profiles are explicit.
void validate_source(const Workspace &w) {
    require(w.complex && w.complex->grade_count() == 2, "PAGER profile requires grades 0 and 1");
    const auto &c = *w.complex;
    const auto &zero = c.grade(rfc::core::Grade{0});
    const auto &one = c.grade(rfc::core::Grade{1});
    require(zero.cell_count() == w.genes.size() && one.cell_count() == w.assertions.size(),
            "source coordinate map disagrees with native carrier");
    const auto &sig = c.signature();
    require(sig.type_count() == 1 && sig.type(0).name == "PAGER.symmetric_participant_set" &&
                sig.aspects().size() == 5,
            "unsupported source signature");
    const auto exists = sig.existence_aspect(), orientation = sig.orientation_aspect();
    std::vector<rfc::composite::AspectId> aspects;
    for (const auto &name : {"undirected", "participant", "rp_observed"})
        aspects.push_back(sig.aspects().id(name));
    std::set<std::string> gene_names, assertion_names;
    std::map<std::string, std::size_t> external_assertions;
    for (std::size_t i = 0; i < w.genes.size(); ++i)
        require(!w.genes[i].empty() && gene_names.insert(w.genes[i]).second &&
                    zero.identities[i] == i + 1 && w.gene_index.at(w.genes[i]) == i,
                "gene identity/coordinate mismatch");
    for (const auto &f : w.files) {
        fs::path path(f.path);
        require(!path.empty() && !path.is_absolute() && f.path.find('\0') == std::string::npos &&
                    digest_valid(f.hash) && !f.parser.empty(),
                "invalid original source-file record");
        for (const auto &part : path)
            require(part != "..", "source-file path escapes original root");
    }
    require(w.contexts.size() == 2 && w.contexts.contains("AD") && w.contexts.contains("HD"),
            "unsupported context profile");
    for (const auto &[name, mask] : w.contexts) {
        (void)name;
        require(mask && mask->boundary_closed(c), "context is not a native subcomplex");
    }
    std::size_t occurrences = 0;
    for (std::size_t i = 0; i < w.assertions.size(); ++i) {
        const auto &a = w.assertions[i];
        ++external_assertions[a.external];
        const auto [begin, end] = one.incidence_range(i);
        require(one.relation_identities[i].value == i + 1 && one.relation_types[i] == 0 &&
                    end - begin == a.participants.size() && !a.participants.empty() &&
                    !a.external.empty() && !a.contexts.empty() && !a.files.empty(),
                "assertion disagrees with native source");
        std::string canonical = std::to_string(a.external.size()) + ":" + a.external;
        std::string previous;
        for (std::size_t j = 0; j < a.participants.size(); ++j) {
            const auto &p = a.participants[j];
            const auto slot = begin + j;
            require(p.gene < w.genes.size() && one.lower_index[slot] == p.gene &&
                        one.occurrence_identities[slot].value == p.occurrence && p.rp >= 0 &&
                        exact_decimal(p.literal) == p.rp && !p.origins.empty(),
                    "participant/observation disagrees with native occurrence");
            require(c.incidence_coefficient(rfc::core::Grade{1}, slot) == 1 &&
                        one.aspects.get(exists, slot) && !one.aspects.get(orientation, slot),
                    "primary assertion must retain positive unit support");
            for (auto aspect : aspects)
                require(one.aspects.get(aspect, slot), "missing participant aspect");
            const auto &gene = w.genes[p.gene];
            require(previous.empty() || previous < gene, "participant identity order");
            previous = gene;
            const auto number = p.rp.get_str();
            canonical += std::to_string(gene.size()) + ":" + gene + std::to_string(number.size()) +
                         ":" + number;
            for (auto origin : p.origins)
                require(origin.file < w.files.size() && origin.record > 0 &&
                            std::find(a.files.begin(), a.files.end(), origin.file) != a.files.end(),
                        "invalid participant source-row origin");
        }
        require(a.id == a.external + "@sha256:" + sha256(canonical) &&
                    assertion_names.insert(a.id).second,
                "assertion content identity mismatch");
        for (auto file : a.files)
            require(file < w.files.size(), "invalid assertion source file");
        for (const auto &name : a.contexts)
            require(w.contexts.contains(name), "unknown assertion context");
        for (const auto &[name, mask] : w.contexts)
            require(mask->contains({rfc::core::Grade{1}, static_cast<std::uint32_t>(i)}) ==
                        a.contexts.contains(name),
                    "context membership disagrees with native block");
        for (const auto &[name, value] : a.observations) {
            (void)value;
            require(a.contexts.contains(name), "statistics belong to absent assertion context");
        }
        occurrences += a.participants.size();
    }
    require(occurrences == one.incidence_count(), "retained occurrence count mismatch");
    auto audit_count = [&](const std::string &name, std::size_t actual) {
        require(count(w.audit.at(name), actual) == actual, "source audit mismatch: " + name);
    };
    audit_count("ambient_grade0", w.genes.size());
    audit_count("ambient_grade1", w.assertions.size());
    audit_count("occurrences", occurrences);
    audit_count("chain_issues", 0);
    for (const auto &[name, mask] : w.contexts)
        for (std::uint16_t grade = 0; grade < 2; ++grade)
            audit_count(name + "_grade" + std::to_string(grade),
                        mask->grade(rfc::core::Grade{grade}).count());
    const auto common = meet(*w.contexts.at("AD"), *w.contexts.at("HD"));
    audit_count("shared_grade0", common.grade(rfc::core::Grade{0}).count());
    audit_count("strict_shared_grade1", common.grade(rfc::core::Grade{1}).count());
    audit_count("same_external_id_different_assertion",
                w.assertions.size() - external_assertions.size());
    require(w.fields.size() == w.field_semantics.size() &&
                w.fields.size() == w.field_presence.size(),
            "field declarations are incomplete");
    for (const auto &context : {"AD", "HD"})
        for (char channel : std::string("GET"))
            require(w.fields.contains(std::string(context) + "." + channel),
                    "missing independent GET channel");
    for (const auto &[name, values] : w.fields) {
        require((name.starts_with("AD.") || name.starts_with("HD.")) && name.size() > 3 &&
                    values.size() == w.genes.size() &&
                    w.field_presence.at(name).size() == w.genes.size(),
                "field carrier/namespace mismatch");
        const auto &meaning = w.field_semantics.at(name);
        require(meaning.at("context").text() == name.substr(0, 2) &&
                    meaning.at("grade").integer() == 0,
                "field meaning/carrier mismatch");
        const auto &present = w.field_presence.at(name);
        for (std::size_t i = 0; i < values.size(); ++i)
            require(present[i] || values[i] == 0, "nonzero value outside reported support");
    }
    for (const auto &[name, rows] : w.measurements) {
        (void)name;
        for (std::size_t i = 0; i < rows.size(); ++i) {
            auto row = rows.at(i);
            (void)count(row.at("source_file"), w.files.size() - 1);
            require(row.at("source_record").integer() > 0, "invalid measurement origin");
        }
    }
    require(rfc::core::bianchi_through(c, c.max_grade()).satisfied(), "native Bianchi failure");
}
Json capture_metadata(const Workspace &w, TemporalEndpoint &endpoint) {
    auto metadata = Json::object();
    metadata.set("schema", Json::string("gadget.snapshot.metadata.v1"));
    metadata.set("source_manifest", w.manifest());
    auto genes = Json::array(), assertions = Json::array();
    for (const auto &name : w.genes)
        genes.push(Json::string(name));
    Vector rp;
    for (const auto &a : w.assertions) {
        auto item = Json::object(), participants = Json::array(), files = Json::array(),
             contexts = Json::array();
        item.set("external", Json::string(a.external));
        item.set("id", Json::string(a.id));
        for (auto file : a.files)
            files.push(Json::integer(file));
        for (const auto &name : a.contexts)
            contexts.push(Json::string(name));
        for (const auto &p : a.participants) {
            auto row = Json::object(), origins = Json::array();
            row.set("gene", Json::integer(p.gene));
            row.set("literal", Json::string(p.literal));
            row.set("occurrence", Json::integer(p.occurrence));
            for (auto origin : p.origins) {
                auto pair = Json::array();
                pair.push(Json::integer(origin.file));
                pair.push(Json::integer(origin.record));
                origins.push(pair);
            }
            row.set("origins", origins);
            participants.push(row);
            rp.push_back(p.rp);
        }
        item.set("participants", participants);
        item.set("files", files);
        item.set("contexts", contexts);
        item.set("observations", json_map(a.observations));
        assertions.push(item);
    }
    metadata.set("genes", genes);
    metadata.set("assertions", assertions);
    metadata.set("measurements", json_map(w.measurements));
    metadata.set("field_semantics", json_map(w.field_semantics));
    auto fields = Json::object(), contexts = Json::object();
    endpoint.carriers.push_back({1,
                                 TemporalCarrierKind::SourceGrade,
                                 rfc::core::Grade{0},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(w.genes.size())});
    std::uint64_t key = 1;
    for (const auto &[name, values] : w.fields) {
        std::vector<std::size_t> active;
        for (std::size_t i = 0; i < values.size(); ++i)
            if (w.field_presence.at(name)[i])
                active.push_back(i);
        TemporalTensor tensor;
        tensor.key = key;
        tensor.role = TemporalTensorRole::Field;
        tensor.output = 1;
        tensor.interpretation = name;
        tensor.tensor = rfc::numeric::SparseRationalTensor(
            rfc::numeric::IndexDomain::linear(values.size()),
            rfc::numeric::StructuralSupport::sparse(values.size(), std::move(active)),
            stored_coordinates(values));
        endpoint.tensors.push_back(std::move(tensor));
        fields.set(name, Json::integer(key++));
    }
    TemporalTensor observation;
    observation.key = key;
    observation.role = TemporalTensorRole::OpaqueExact;
    observation.interpretation = "PAGER.RP in retained grade-one occurrence order";
    observation.tensor = rfc::numeric::SparseRationalTensor(
        rfc::numeric::IndexDomain::linear(rp.size()),
        rfc::numeric::StructuralSupport::full(rp.size()), stored_coordinates(rp));
    endpoint.tensors.push_back(std::move(observation));
    metadata.set("rp_tensor", Json::integer(key));
    metadata.set("fields", fields);
    key = 1;
    for (const auto &[name, mask] : w.contexts) {
        LocalBlockDefinition block;
        block.block_key = key;
        block.kind = LocalBlockKind::StrictSubcomplex;
        for (std::uint16_t grade = 0; grade < w.complex->grade_count(); ++grade)
            mask->grade(rfc::core::Grade{grade}).for_each_set([&](std::size_t i) {
                block.cells.push_back(stable_cell_key(
                    *w.complex, {rfc::core::Grade{grade}, static_cast<std::uint32_t>(i)}));
            });
        endpoint.local_blocks.push_back(std::move(block));
        contexts.set(name, Json::integer(key++));
    }
    metadata.set("contexts", contexts);
    return metadata;
}
} // namespace

NativeSnapshot::NativeSnapshot(const rfc::storage::TemporalEndpoint &endpoint,
                               std::uint64_t source_namespace,
                               std::map<std::string, std::uint64_t> fields)
    : directory_(
          std::make_shared<const rfc::storage::NativeAccessionDirectory>(endpoint, source_namespace,
              std::vector<rfc::storage::RelationalCodeObservationQuery>{},
              [](auto bytes){return rfc::storage::rfcx_detail::fnv1a64(bytes);},native_accession_limits())),
      fields_(std::move(fields)) {}
NativeSnapshot::NativeSnapshot(const NativeRecord &publication,std::map<std::string,std::uint64_t> fields)
    :directory_(publication.accession_directory()),fields_(std::move(fields)) {}
std::shared_ptr<const rfc::storage::TemporalEndpoint> NativeSnapshot::endpoint() const {
    return directory_->pin_snapshot();
}
rfc::storage::PinnedNativeAccession NativeSnapshot::field(const std::string &name) const {
    rfc::storage::NativeAccessionRequest request;
    request.kind = rfc::storage::NativeAccessionKind::Field;
    request.subject = fields_.at(name);
    return directory_->pin(request);
}
rfc::storage::PinnedNativeAccession NativeSnapshot::source_grade(rfc::core::Grade grade) const {
    auto request = rfc::storage::NativeAccessionRequest::chain(0, grade);
    request.inline_value = true;
    return directory_->pin(request);
}

NativeRecord capture_workspace(const Workspace &w) {
    using namespace rfc::storage;
    validate_source(w);
    TemporalEndpoint endpoint(*w.complex);
    auto metadata = capture_metadata(w, endpoint);
    using N=NativeRecordValue;
    auto creator=Json::object();creator.set("rexgraph_headers_sha256",Json::string(GADGET_REXGRAPH_SHA256));creator.set("gadget_sources_sha256",Json::string(GADGET_SOURCE_SHA256));
    N::Array original_files;std::vector<RawByteSource> raw_sources;
    require(w.raw_sources.size()==w.files.size(),"original source bytes were not retained at admission");
    for(std::size_t i=0;i<w.files.size();++i) {
        const auto &file=w.files[i];const auto &raw=w.raw_sources.at(i);raw_sources.push_back(raw);
        original_files.elements.emplace_back(N::Object{{{"path",N(file.path)},
            {"bytes",N(N::ByteSpan{raw.key(),0,raw.bytes().size(),false})}}});
    }
    std::optional<NativeRecord> publication;
    if(w.native_publication) {
        auto old=*w.native_publication;
        auto old_numerical=*old.endpoint();old_numerical.declarations.reset();
        if(render_record_view(old.at("study")).dump()==metadata.dump() &&
           encode_temporal_endpoint(old_numerical)==encode_temporal_endpoint(endpoint))publication=std::move(old);
    }
    if(publication)return *publication;
    auto meanings=import_record_view(metadata);
    auto &members=std::get<N::Object>(meanings.data).members;
    for(auto &[name,value]:members)if(name=="genes" || name=="assertions") {
        auto &rows=std::get<N::Array>(value.data).elements;
        for(std::size_t i=0;i<rows.size();++i) {
            if(name=="assertions") {
                auto &attributes=std::get<N::Object>(rows[i].data).members;
                for(auto &[attribute,array]:attributes)if(attribute=="participants") {
                    auto &participants=std::get<N::Array>(array.data).elements;
                    for(std::size_t j=0;j<participants.size();++j)participants[j]=N(N::Owned{0,{},StableOccurrenceKey{rfc::core::Grade{1},rfc::core::OccurrenceKey{w.assertions[i].participants[j].occurrence}},std::make_shared<const N>(std::move(participants[j]))});
                }
            }
            const auto owner=stable_cell_key(*w.complex,{rfc::core::Grade{static_cast<std::uint16_t>(name=="genes"?0:1)},static_cast<std::uint32_t>(i)});
            rows[i]=N(N::Owned{0,owner,{},std::make_shared<const N>(std::move(rows[i]))});
        }
    }
    if(!publication)publication=NativeRecord::create(N(N::Object{{
        {"schema",N(std::string("gadget.workspace.native.v1"))},
        {"study",std::move(meanings)},
        {"creator",import_record_view(creator)},
        {"original_files",N(std::move(original_files))}}}),study_namespace,std::move(raw_sources),&endpoint);
    return *publication;
}
void admit_workspace(Workspace &w) {
    w.native_publication=capture_workspace(w);
    const auto &publication=*w.native_publication;
    std::map<std::string,std::uint64_t> fields;
    for(const auto &tensor:publication.endpoint()->tensors)
        if(tensor.role==rfc::storage::TemporalTensorRole::Field)fields.emplace(tensor.interpretation,tensor.key);
    w.native_snapshot=std::make_shared<const NativeSnapshot>(publication,std::move(fields));
    const auto study=publication.at("study");
    w.audit=study.at("source_manifest").at("audit");
    for(auto &[name,value]:w.measurements)value=study.at("measurements").at(name);
    for(auto &[name,value]:w.field_semantics)value=study.at("field_semantics").at(name);
    for(std::size_t i=0;i<w.assertions.size();++i)
        for(auto &[name,value]:w.assertions[i].observations)value=study.at("assertions").at(i).at("observations").at(name);
}
Json save_snapshot(const Workspace &w,const fs::path &requested,bool export_json_views) {
    using namespace rfc::storage;
    validate_source(w);
    validate_file_path(requested);
    require(!requested.filename().empty() && requested.filename() != "." &&
                requested.filename() != "..",
            "destination must name a new directory");
    const auto destination = fs::absolute(requested);
    require(!fs::exists(fs::symlink_status(destination)), "destination already exists");
    const auto stage = destination.parent_path() / (destination.filename().string() + ".partial");
    fs::create_directories(destination.parent_path());
    require(fs::create_directory(stage), "staging directory already exists");
    fs::permissions(stage, fs::perms::owner_all);
    const auto publication=capture_workspace(w);
    const auto metadata=render_record_view(publication.at("study"));
    const auto text=compact(metadata)+"\n";
    const auto creator=render_record_view(publication.at("creator"));
    const auto native=publication.encode();
    write_file_exact(stage / "study.rfcxs", native, FilePublication::CreateNew);
    const auto completion=completion_view(metadata,native,creator);
    if(export_json_views) {
        write_file_exact(stage/"metadata.json",bytes(text),FilePublication::CreateNew);
        const auto manifest=completion.dump()+"\n";
        write_file_exact(stage/"snapshot.json",bytes(manifest),FilePublication::CreateNew);
    }
    require(!fs::exists(fs::symlink_status(destination)), "destination appeared during capture");
    fs::rename(stage, destination);
    return completion;
}

Workspace load_snapshot(const fs::path &requested) {
    using namespace rfc::storage;
    validate_file_path(requested);
    require(fs::is_directory(fs::symlink_status(requested)), "expected regular snapshot directory");
    const auto root = fs::canonical(requested);
    const auto publication=NativeRecord::read(root/"study.rfcxs");
    require(publication.at("schema").text()=="gadget.workspace.native.v1","unsupported native workspace declaration");
    auto endpoint=publication.endpoint();
    require(publication.source_namespace()==study_namespace && endpoint->sequence==0,"snapshot is not a single source publication");
    auto metadata=render_record_view(publication.at("study"));
    const auto source_bytes=artifact(root,"study.rfcxs",artifact_limit);
    const auto creator=render_record_view(publication.at("creator"));
    const auto completion=completion_view(metadata,source_bytes,creator);
    const auto completion_text=completion.dump()+"\n";
    const auto completion_bytes=std::vector<std::byte>(bytes(completion_text).begin(),bytes(completion_text).end());
    if(fs::exists(root/"snapshot.json"))require(json_bytes(artifact(root,"snapshot.json",artifact_limit)).dump()==completion.dump(),"optional completion view differs from native publication");
    if(fs::exists(root/"metadata.json"))require(json_bytes(artifact(root,"metadata.json",artifact_limit)).dump()==metadata.dump(),"optional biomedical view differs from canonical declaration");
    require(endpoint->finite_systems.empty() && endpoint->carriers.size()==1 && endpoint->carriers[0].key==1 && endpoint->carriers[0].kind==TemporalCarrierKind::SourceGrade && endpoint->carriers[0].grade==rfc::core::Grade{0},"unsupported native workspace channels");
    members(metadata, {"schema", "source_manifest", "genes", "assertions", "measurements",
                       "field_semantics", "fields", "rp_tensor", "contexts"});
    require(metadata.at("schema").text() == "gadget.snapshot.metadata.v1",
            "unsupported metadata schema");
    Workspace w;w.native_publication=publication;
    const auto manifest = metadata.at("source_manifest");
    members(manifest, {"schema", "source_digest", "root", "audit", "files", "field_semantics",
                       "boundary_contract"});
    require(manifest.at("schema").text() == "gadget.source.v1",
            "unsupported original source manifest");
    w.root = manifest.at("root").text();
    validate_file_path(w.root);
    require(w.root.is_absolute(), "original source root must be absolute");
    w.audit = publication.at("study").at("source_manifest").at("audit");
    auto files = manifest.at("files"), genes = metadata.at("genes");
    require(files.size() > 0 && genes.size() > 0,
            "source file/gene budget");
    for (std::size_t i = 0; i < files.size(); ++i) {
        auto f = files.at(i);
        members(f, {"index", "path", "sha256", "bytes", "parser"});
        require(count(f.at("index"), files.size()) == i, "original file index mismatch");
        w.files.push_back({f.at("path").text(), f.at("sha256").text(), f.at("parser").text(),
                           count(f.at("bytes"), artifact_limit)});
    }
    const auto original_files=publication.at("original_files");
    require(original_files.size()==w.files.size(),"incomplete original source byte map");
    for(std::size_t i=0;i<original_files.size();++i) {
        const auto row=original_files.at(i);const auto path=row.at("path").text();
        require(path==w.files[i].path,"original file identity differs from biomedical map");
        const auto raw=row.at("bytes").byte_source();const auto data=raw.bytes();
        require(data.size()==w.files[i].bytes && digest(data)==w.files[i].hash,"original raw byte binding mismatch");
        w.raw_sources.push_back(raw);
    }
    refresh_source_digest(w);
    require(w.source_digest == manifest.at("source_digest").text() &&
                w.source_digest == completion.at("source_digest").text(),
            "original source digest mismatch");
    for (std::size_t i = 0; i < genes.size(); ++i) {
        auto name = genes.at(i).text();
        require(w.gene_index.emplace(name, i).second, "duplicate gene label");
        w.genes.push_back(std::move(name));
    }
    std::map<std::uint64_t, const TemporalTensor *> tensors;
    for (const auto &tensor : endpoint->tensors)
        require(tensors.emplace(tensor.key, &tensor).second, "duplicate native tensor key");
    const auto rp_key = count(metadata.at("rp_tensor"), std::numeric_limits<std::uint64_t>::max());
    require(tensors.contains(rp_key), "missing RP observation tensor");
    const auto &rp = *tensors.at(rp_key);
    require(endpoint->source.grade_count() == 2, "unsupported source grade profile");
    const auto occurrences = endpoint->source.grade(rfc::core::Grade{1}).incidence_count();
    require(rp.role == TemporalTensorRole::OpaqueExact && !rp.input && !rp.output &&
                rp.owner.kind == TemporalOwnerKind::Global &&
                rp.interpretation == "PAGER.RP in retained grade-one occurrence order" &&
                rp.tensor.domain() == rfc::numeric::IndexDomain::linear(occurrences) &&
                rp.tensor.support().count() == occurrences,
            "invalid RP observation declaration");
    auto rp_values = rp.tensor.values().dense();
    auto assertions = metadata.at("assertions");
    require(assertions.size() == endpoint->source.grade(rfc::core::Grade{1}).cell_count(),
            "assertion map dimension mismatch");
    std::size_t slot = 0;
    for (std::size_t i = 0; i < assertions.size(); ++i) {
        auto item = assertions.at(i);
        members(item, {"external", "id", "participants", "files", "contexts", "observations"});
        Assertion a;
        a.external = item.at("external").text();
        a.id = item.at("id").text();
        auto participants = item.at("participants");
        const auto [begin, end] = endpoint->source.grade(rfc::core::Grade{1}).incidence_range(i);
        require(participants.size() == end - begin && slot == begin,
                "participant map dimension mismatch");
        for (std::size_t j = 0; j < participants.size(); ++j) {
            const auto row = participants.at(j);
            members(row, {"gene", "literal", "occurrence", "origins"});
            Participant p;
            p.gene = static_cast<std::uint32_t>(count(row.at("gene"), w.genes.size() - 1));
            p.literal = row.at("literal").text();
            p.occurrence = count(row.at("occurrence"), INT64_MAX);
            p.rp = rp_values.at(slot++);
            const auto origins = row.at("origins");
            for (std::size_t n = 0; n < origins.size(); ++n) {
                const auto origin = origins.at(n);
                require(origin.size() == 2, "invalid source-row origin pair");
                p.origins.push_back(
                    {count(origin.at(0), w.files.size() - 1), count(origin.at(1), INT64_MAX)});
            }
            a.participants.push_back(std::move(p));
        }
        auto input_files = item.at("files"), contexts = item.at("contexts");
        for (std::size_t j = 0; j < input_files.size(); ++j)
            a.files.push_back(count(input_files.at(j), w.files.size() - 1));
        for (std::size_t j = 0; j < contexts.size(); ++j)
            require(a.contexts.insert(contexts.at(j).text()).second, "duplicate assertion context");
        entries(item.at("observations"),
                [&](const auto &name, const auto &) { a.observations.emplace(name, publication.at("study").at("assertions").at(i).at("observations").at(name)); });
        w.assertions.push_back(std::move(a));
    }
    assertions = Json{};
    Vector().swap(rp_values);
    entries(metadata.at("measurements"),
            [&](const auto &name, const auto &) { w.measurements.emplace(name, publication.at("study").at("measurements").at(name)); });
    entries(metadata.at("field_semantics"), [&](const auto &name, const auto &) {
        w.field_semantics.emplace(name, publication.at("study").at("field_semantics").at(name));
    });
    std::map<std::string, std::uint64_t> field_keys;
    std::set<std::uint64_t> used{rp_key};
    const auto declared_fields = metadata.at("fields");
    require(json_object_is_type(declared_fields.raw(), json_type_object), "expected field map");
    const auto field_count =
        static_cast<std::size_t>(json_object_object_length(declared_fields.raw()));
    (void)rfc::machine::checked_size_mul(field_count,w.genes.size());
    entries(metadata.at("fields"), [&](const auto &name, const auto &value) {
        auto key = count(value, std::numeric_limits<std::uint64_t>::max());
        require(tensors.contains(key) && used.insert(key).second, "missing/aliased field tensor");
        const auto &tensor = *tensors.at(key);
        require(tensor.role == TemporalTensorRole::Field && tensor.output == 1 && !tensor.input &&
                    tensor.owner.kind == TemporalOwnerKind::Global && tensor.interpretation == name,
                "field declaration mismatch");
        auto values = tensor.tensor.values().dense();
        require(values.size() == w.genes.size(), "field dimension mismatch");
        std::vector<bool> present(w.genes.size(), false);
        tensor.tensor.support().for_each_active([&](std::size_t i) { present.at(i) = true; });
        w.fields.emplace(name, std::move(values));
        w.field_presence.emplace(name, std::move(present));
        field_keys.emplace(name, key);
    });
    require(used.size() == tensors.size(), "unbound native tensor");
    w.complex = std::shared_ptr<const rfc::core::RelationalComplex>(endpoint, &endpoint->source);
    const StableSourceIndex index(*w.complex);
    std::set<std::uint64_t> used_blocks;
    entries(metadata.at("contexts"), [&](const auto &name, const auto &value) {
        auto key = count(value, 2);
        require(used_blocks.insert(key).second, "aliased context block");
        auto block = std::find_if(endpoint->local_blocks.begin(), endpoint->local_blocks.end(),
                                  [&](const auto &b) { return b.block_key == key; });
        require(block != endpoint->local_blocks.end() &&
                    block->kind == LocalBlockKind::StrictSubcomplex && block->payloads.empty(),
                "missing or unsupported context block");
        auto mask = std::make_unique<rfc::algebra::SubcomplexMask>(*w.complex);
        for (auto cell : block->cells) {
            auto coordinate = index.cell(cell);
            require(coordinate.has_value(), "context references missing source cell");
            mask->select(*coordinate);
        }
        w.contexts.emplace(name, std::move(mask));
    });
    require(used_blocks.size() == endpoint->local_blocks.size(), "unbound native context block");
    metadata = Json{};
    genes = Json{};
    validate_source(w);
    require(compact(w.manifest()) == compact(manifest), "source manifest did not round-trip");
    w.native_snapshot =
        std::make_shared<NativeSnapshot>(publication,std::move(field_keys));
    const auto pinned = w.native_snapshot->endpoint();
    w.complex = std::shared_ptr<const rfc::core::RelationalComplex>(pinned, &pinned->source);
    w.snapshot_origin = Json::parse(completion.dump());
    w.snapshot_origin.set("directory", Json::string(root.string()));
    w.snapshot_origin.set("completion_sha256", Json::string(digest(completion_bytes)));
    return w;
}
} // namespace gadget
