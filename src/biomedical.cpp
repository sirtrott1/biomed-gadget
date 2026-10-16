#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <cmath>
#include <gadget/biomedical.hpp>
#include <gadget/table.hpp>
#include <rfc/database/interpretation_codec.hpp>
#include <rfc/database/native_query.hpp>
#include <rfc/database/native_role_join.hpp>

namespace gadget {
NativeRecord retain_study_sources(const ScientificStudy &study, const NativeRecord &request) {
    using namespace scientific_detail;
    using N = NativeRecordValue;
    const auto declaration = render_record_view(request);
    members(declaration, {"schema", "source_digest", "method", "sources"});
    require(declaration.at("schema").text() == "gadget.biomedical.sources.request.v1" &&
                declaration.at("source_digest").text() == study.digest(),
            "original source request differs from admitted study");
    const auto method = text(declaration.at("method"));
    const auto &original = study.native_publication();
    require(!original.has("source_provenance"), "study already has declared source provenance");
    auto value = original.materialize_value();
    auto &root = std::get<N::Object>(value.data).members;
    std::uint64_t next = 1;
    for (const auto key : NativeRecord::source_keys(*original.endpoint()))
        next = rfc::machine::checked_size_add(std::max(next, key), 1);
    std::vector<RawByteSource> added;
    N::Array sources;
    std::set<std::string> identities;
    const auto inputs = declaration.at("sources");
    source_array(inputs, "original source declarations");
    require(inputs.size() > 0, "original source declarations are empty");
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        const auto item = inputs.at(i);
        members(item, {"id", "path", "sha256", "source_url", "license"});
        const auto id = text(item.at("id"));
        require(identities.insert(id).second, "duplicate original source identity");
        const fs::path path(text(item.at("path")));
        require(fs::is_regular_file(fs::symlink_status(path)), "original source must be regular");
        auto bytes = RawByteSource::read(next, path);
        const auto span = bytes.bytes();
        const auto digest = "sha256:" + sha256(std::string_view(
            reinterpret_cast<const char *>(span.data()), span.size()));
        require(item.at("sha256").text() == digest, "original source SHA-256 mismatch");
        sources.elements.emplace_back(N::Object{{
            {"id", N(id)}, {"sha256", N(digest)},
            {"source_url", N(text(item.at("source_url")))},
            {"license", N(text(item.at("license")))},
            {"bytes", N(N::ByteSpan{next, 0, span.size(), false})}}});
        added.push_back(std::move(bytes));
        next = rfc::machine::checked_size_add(next, 1);
    }
    auto archive = request.owning_publication().archive_source(next);
    const auto archive_size = archive.bytes().size();
    added.push_back(std::move(archive));
    root.emplace_back("source_provenance", N(N::Object{{
        {"schema", N(std::string("gadget.biomedical.sources.v1"))},
        {"source_digest", N(study.digest())}, {"method", N(method)},
        {"sources", N(std::move(sources))},
        {"request_publication", N(N::ByteSpan{next, 0, archive_size, false})}}}));
    return original.derive(std::move(value), std::move(added));
}
namespace {
using namespace scientific_detail;
using namespace rfc::database;
using V = rfc::composite::DeclaredValue;
constexpr std::uint64_t source_namespace = 72101, binding_namespace = 72102;
const std::vector<std::string> relations{"", "=", "<", "<=", ">", ">=", "~", "other"};
const std::vector<std::string> endpoints{"", "Kd", "Ki", "IC50", "EC50", "other"};
const std::vector<std::string> units{"", "nM", "uM", "M", "other"};
std::uint64_t code(const std::vector<std::string> &values, const std::string &v) {
    const auto i = std::find(values.begin(), values.end(), v);
    return i == values.end() ? values.size() - 1
                             : static_cast<std::uint64_t>(std::distance(values.begin(), i));
}
NativeRoleJoinRequest role_join_request(const Json &query, const InterpretedSnapshot &source) {
    const std::map<std::string, std::uint64_t> kinds{{"binding_observation",10},
        {"identifier_cross_reference",20},{"source_relationship",30}};
    const std::map<std::string, std::uint64_t> roles{{"compound_or_left",11},
        {"target_or_right",12},{"assay",13},{"document",14},{"response",15}};
    const auto side = [&](const char *kind_name, const char *role_name) {
        const auto kind = query.at(kind_name).text(), role = query.at(role_name).text();
        require(kinds.contains(kind) && roles.contains(role), "unknown declared biomedical join kind or role");
        const auto &types = source.schema().assertions;
        const auto type = std::find_if(types.begin(),types.end(),[&](const auto &entry) { return entry.key==kinds.at(kind); });
        require(type!=types.end(), "biomedical join kind is absent from this source adapter");
        require(std::any_of(type->roles.begin(),type->roles.end(),[&](const auto &entry) { return entry.key==roles.at(role); }),
                "biomedical join role is absent from its original declaration");
        return std::pair{kinds.at(kind),roles.at(role)};
    };
    const auto left = side("left_kind","left_role"), right = side("right_kind","right_role");
    const auto activity = query.at("activity").text();
    require(activity=="retained" || activity=="active_only", "unknown biomedical join activity policy");
    return {left.first,left.second,right.first,right.second,
        activity=="retained" ? NativeRoleActivity::Retained : NativeRoleActivity::ActiveOnly};
}
std::string string_or_empty(const Json &j, const std::string &name) {
    if (!j.has(name) || !j.at(name).raw())
        return {};
    if (json_object_is_type(j.at(name).raw(), json_type_string))
        return text(j.at(name));
    if (json_object_is_type(j.at(name).raw(), json_type_int))
        return std::to_string(j.at(name).integer());
    throw std::invalid_argument("source identity/value must be a string or integer");
}
void verify_seal(const Json &record, const std::string &member) {
    auto copy = Json::parse(record.dump());
    const auto expected = copy.at(member).text();
    json_object_object_del(copy.raw(), member.c_str());
    require(expected == "sha256:" + sha256(copy.dump()), "evidence manifest digest mismatch");
}
InterpretationSchema schema(std::uint32_t version) {
    using D = rfc::composite::Declaration;
    const auto id = [](std::uint64_t n) { return rfc::composite::DeclarationId{72103, n, 1}; };
    InterpretationSchema out;
    out.source_namespace = source_namespace;
    out.binding_namespace = binding_namespace;
    out.values = rfc::composite::DeclarationRegistry(
        {D::rational(id(1)), D::optional(id(2), id(1)), D::finite(id(3), {0, 1, 2, 3, 4, 5, 6, 7}),
         D::finite(id(4), {0, 1, 2, 3, 4, 5}), D::finite(id(5), {0, 1, 2, 3, 4}),
         D::finite(id(6), {1, 2, 3, 4, 5}),
         D::product(id(7),
                    {{1, id(2), 0}, {2, id(3), 0}, {3, id(4), 0}, {4, id(5), 0}, {5, id(6), 0}})});
    out.subjects = {{1, "namespaced biomedical source identities", "database_identifier",
                     "UTF-8 namespaced token"}};
    const std::vector<ExternalRole> roles{{11, 1, 0, 1, true, false},
                                          {12, 1, 0, 1, true, false},
                                          {13, 1, 0, 1, true, false},
                                          {14, 1, 0, 1, true, false},
                                          {15, 1, 1, 1, false, false}};
    const std::vector<ExternalFieldLabel> fields{
        {1, "reported_number"}, {2, "relation"}, {3, "endpoint"}, {4, "units"}, {5, "provider"}};
    const std::vector<ExternalFieldLabel> labels{{11, "compound_or_left"},
                                                 {12, "target_or_right"},
                                                 {13, "assay"},
                                                 {14, "document"},
                                                 {15, "response"}};
    out.assertions = {{10, "reported_binding_observation", id(7), roles, fields, "source_record",
                       "provider identifier or response-local record pointer", labels},
                      {20, "reported_identifier_cross_reference", id(7), roles, fields,
                       "source_record", "response-local mapping pointer", labels}};
    if (version >= 3)
        out.assertions.push_back({30, "reported_source_relationship", id(7), roles, fields,
            "source_record", "provider record or response-local statement pointer", labels});
    return out;
}
struct Packet {
    std::uint32_t adapter_version{2};
    Json manifest, records{Json::array()};
    std::map<std::string, Json> entities;
};
Packet normalize(const NativeRecord &database, std::uint32_t version, const fs::path &directory={}) {
    require(version == 1 || version == 2 || version == 3, "unsupported biomedical adapter version");
    Packet p;
    p.adapter_version = version;
    p.manifest=render_record_view(database.at("snapshot"));
    const auto retained_artifact=[&](const Json &record) {
        const auto items=database.at("artifacts");const auto path=record.at("path").text();
        for(std::size_t i=0;i<items.size();++i)if(items.at(i).at("path").text()==path) {
            const auto b=items.at(i).at("bytes").bytes();const std::string bytes(reinterpret_cast<const char *>(b.data()),b.size());
            require("sha256:"+sha256(bytes)==record.at("sha256").text(),"native database byte binding mismatch");
            if(!directory.empty() && fs::exists(directory/path))require(local_artifact(directory,record)==bytes,"external database artifact differs from native source");
            return bytes;
        }
        throw std::invalid_argument("database source bytes missing from native publication");
    };
    require(p.manifest.at("schema").text() == "gadget.database.snapshot.v1",
            "unsupported frozen database snapshot");
    verify_seal(p.manifest, "snapshot_sha256");
    const auto responses = p.manifest.at("records");
    require(responses.size(), "evidence snapshot requires responses");
    std::set<std::string> ids;
    const auto entity = [&](const std::string &id, const std::string &kind,
                            const std::string &label) {
        if (id.empty())
            return;
        if (!p.entities.contains(id)) {
            auto j = Json::object();
            j.set("id", Json::string(id));
            j.set("kind", Json::string(kind));
            j.set("label", Json::string(label.empty() ? id : label));
            p.entities.emplace(id, j);
        }
    };
    for (std::size_t i = 0; i < responses.size(); ++i) {
        const auto receipt = responses.at(i), request = receipt.at("request");
        const auto rid = request.at("id").text(), provider = request.at("provider").text(),
                   resource = request.at("resource").text();
        require(ids.insert(rid).second, "duplicate response identity");
        const auto bytes = retained_artifact(receipt);
        auto headers = Json::object();
        headers.set("path", receipt.at("headers_path"));
        headers.set("sha256", receipt.at("headers_sha256"));
        (void)retained_artifact(headers);
        const auto raw = Json::parse(bytes);
        const auto response = "response:" + rid + ":" + receipt.at("sha256").text();
        entity(response, "frozen_api_response", rid);
        const auto add =
            [&](const Json &original, const std::string &pointer, const std::string &identity,
                const std::string &compound, const std::string &target, const std::string &assay,
                const std::string &document, const std::string &value, const std::string &relation,
                const std::string &endpoint, const std::string &unit, std::uint64_t kind) {
                auto item = Json::object();
                item.set("id", Json::string(rid + ":" + identity));
                item.set("kind", Json::integer(static_cast<std::int64_t>(kind)));
                item.set("provider", Json::string(provider));
                item.set("response", Json::string(response));
                item.set("compound", Json::string(compound));
                item.set("target", Json::string(target));
                item.set("assay", Json::string(assay));
                item.set("document", Json::string(document));
                item.set("value_token", Json::string(value));
                item.set("relation", Json::string(relation));
                item.set("endpoint", Json::string(endpoint));
                item.set("units", Json::string(unit));
                item.set("source_pointer", Json::string(pointer));
                item.set("source_receipt", receipt);
                item.set("original", original);
                if (!value.empty())
                    item.set("exact_reported_number", exact_number(exact_decimal(value)));
                p.records.push(item);
            };
        const auto mapping = [&](const std::string &left, const std::string &right,
                                 const std::string &pointer, const Json &original) {
            entity(left, "source_identity", left);
            entity(right, "source_identity", right);
            add(original, pointer, pointer, left, right, "", "", "", "", "", "", 20);
        };
        const auto relationship = [&](const Json &original, const std::string &pointer,
                                      const std::string &identity, const std::string &left,
                                      const std::string &right, const std::string &meaning) {
            require(version >= 3, "source relationships require biomedical adapter version 3");
            add(original, pointer, identity, left, right, "", "", "", "", "", "", 30);
            p.records.at(p.records.size()-1).set("relationship", Json::string(meaning));
        };
        if (provider == "chembl" && resource == "activity") {
            const auto rows = raw.at("activities");
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const auto row = rows.at(r);
                const auto compound = "chembl:" + string_or_empty(row, "molecule_chembl_id"),
                           target = "chembl:" + string_or_empty(row, "target_chembl_id"),
                           assay = "chembl:" + string_or_empty(row, "assay_chembl_id"),
                           document = "chembl:" + string_or_empty(row, "document_chembl_id");
                require(compound != "chembl:" && target != "chembl:" && assay != "chembl:" &&
                            document != "chembl:",
                        "activity lacks declared producer identities");
                entity(compound, "compound_record", string_or_empty(row, "molecule_pref_name"));
                entity(target, "target_record", string_or_empty(row, "target_pref_name"));
                entity(assay, "assay_record", string_or_empty(row, "assay_description"));
                entity(document, "publication_record", document);
                const auto activity = string_or_empty(row, "activity_id");
                require(!activity.empty(), "ChEMBL activity lacks stable record ID");
                add(row, "/activities/" + std::to_string(r), "activity:" + activity, compound,
                    target, assay, document, string_or_empty(row, "standard_value"),
                    string_or_empty(row, "standard_relation"),
                    string_or_empty(row, "standard_type"), string_or_empty(row, "standard_units"),
                    10);
            }
        } else if (provider == "bindingdb" && resource == "uniprot") {
            const auto rows = raw.at("getLindsByUniprotsResponse").at("affinities");
            const auto target = "uniprot:" + request.at("identifier").text();
            entity(target, "protein_record", target);
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const auto row = rows.at(r);
                auto value = trim(string_or_empty(row, "affinity")), relation = std::string("=");
                require(!value.empty(), "BindingDB compact response lacks an affinity token");
                if (value.starts_with("<") || value.starts_with(">")) {
                    relation = value.substr(0, value.size() > 1 && value[1] == '=' ? 2 : 1);
                    value = trim(value.substr(relation.size()));
                }
                const auto compound = "bindingdb:monomer:" + string_or_empty(row, "monomerid");
                require(compound != "bindingdb:monomer:", "BindingDB record lacks monomer ID");
                entity(compound, "compound_record", compound);
                const auto doi = string_or_empty(row, "doi");
                const auto document = doi.empty() ? std::string{} : "doi:" + doi;
                entity(document, "publication_record", doi);
                const auto pointer = "/getLindsByUniprotsResponse/affinities/" + std::to_string(r);
                // This response has no measurement/assay ID. Its pointer is local
                // to these exact bytes, never a stable database measurement key.
                add(row, pointer, response + pointer, compound, target, "", document, value,
                    relation, string_or_empty(row, "affinity_type"), "nM", 10);
            }
        } else if (provider == "chembl" && resource == "target") {
            const auto left = "chembl:" + raw.at("target_chembl_id").text();
            entity(left, "target_record", string_or_empty(raw, "pref_name"));
            const auto components = raw.at("target_components");
            for (std::size_t r = 0; r < components.size(); ++r) {
                const auto accession = string_or_empty(components.at(r), "accession");
                if (!accession.empty())
                    mapping(left, "uniprot:" + accession,
                            "/target_components/" + std::to_string(r) + "/accession",
                            components.at(r));
            }
        } else if (provider == "chembl" && resource == "molecule") {
            const auto left = "chembl:" + raw.at("molecule_chembl_id").text();
            entity(left, "compound_record", string_or_empty(raw, "pref_name"));
            const auto refs = raw.at("cross_references");
            for (std::size_t r = 0; r < refs.size(); ++r)
                if (string_or_empty(refs.at(r), "xref_src") == "PubChem")
                    mapping(left, "pubchem:sid:" + string_or_empty(refs.at(r), "xref_id"),
                            "/cross_references/" + std::to_string(r), refs.at(r));
            // ChEMBL PubChem xrefs are Substance IDs, not Compound IDs.
            if (version >= 2 && raw.has("molecule_structures") &&
                raw.at("molecule_structures").raw()) {
                const auto key =
                    string_or_empty(raw.at("molecule_structures"), "standard_inchi_key");
                if (!key.empty())
                    mapping(left, "inchikey:" + key, "/molecule_structures/standard_inchi_key",
                            raw.at("molecule_structures"));
            }
        } else if (provider == "uniprot" && resource == "protein") {
            const auto left = "uniprot:" + raw.at("primaryAccession").text();
            entity(left, "protein_record", string_or_empty(raw, "uniProtkbId"));
            const auto refs = raw.at("uniProtKBCrossReferences");
            for (std::size_t r = 0; r < refs.size(); ++r) {
                const auto db = refs.at(r).at("database").text();
                if (db == "PDB")
                    mapping(left, "pdb:" + refs.at(r).at("id").text(),
                            "/uniProtKBCrossReferences/" + std::to_string(r), refs.at(r));
            }
            if (version >= 3 && raw.has("comments")) {
                const auto comments = raw.at("comments");
                for (std::size_t r = 0; r < comments.size(); ++r) {
                    const auto comment = comments.at(r);
                    if (string_or_empty(comment, "commentType") != "DISEASE") continue;
                    const auto pointer = "/comments/" + std::to_string(r);
                    std::string right;
                    if (comment.has("disease") && comment.at("disease").raw() &&
                        !string_or_empty(comment.at("disease"), "diseaseAccession").empty()) {
                        const auto disease = comment.at("disease");
                        right = "uniprot:disease:" + disease.at("diseaseAccession").text();
                        entity(right, "disease_record", string_or_empty(disease, "diseaseId"));
                        if (disease.has("diseaseCrossReference") && disease.at("diseaseCrossReference").raw()) {
                            const auto reference = disease.at("diseaseCrossReference");
                            if (string_or_empty(reference, "database") == "MIM" &&
                                !string_or_empty(reference, "id").empty())
                                mapping(right, "mim:" + reference.at("id").text(),
                                        pointer + "/disease/diseaseCrossReference", reference);
                        }
                    } else {
                        right = response + pointer;
                        entity(right, "source_statement", "Original disease statement without one declared disease accession");
                    }
                    relationship(comment, pointer, pointer, left, right, "reported_disease_annotation");
                }
            }
        } else if (provider == "chembl" && resource == "drug_indication") {
            require(version >= 3, "drug indications require biomedical adapter version 3");
            const auto rows = raw.at("drug_indications");
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const auto row = rows.at(r);
                require(row.at("molecule_chembl_id").text()==request.at("identifier").text(),
                        "drug indication molecule differs from the explicitly requested source");
                const auto compound = "chembl:" + row.at("molecule_chembl_id").text();
                const auto efo = string_or_empty(row, "efo_id"), mesh = string_or_empty(row, "mesh_id");
                require(!efo.empty() || !mesh.empty(), "drug indication lacks an original condition identifier");
                entity(compound, "compound_record", compound);
                const auto efo_id = efo.empty() ? std::string{} : "efo:" + efo;
                const auto mesh_id = mesh.empty() ? std::string{} : "mesh:" + mesh;
                if (!efo.empty()) entity(efo_id, "condition_record", string_or_empty(row, "efo_term"));
                if (!mesh.empty()) entity(mesh_id, "condition_record", string_or_empty(row, "mesh_heading"));
                const auto pointer = "/drug_indications/" + std::to_string(r);
                relationship(row, pointer, pointer, compound, efo.empty() ? mesh_id : efo_id,
                             "reported_drug_indication");
                if (!efo.empty() && !mesh.empty())
                    mapping(efo_id, mesh_id, pointer + "/condition_cross_reference", row);
            }
        } else if (provider == "pubchem" && resource == "compound") {
            const auto rows = raw.at("PropertyTable").at("Properties");
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const auto id = "pubchem:cid:" + string_or_empty(rows.at(r), "CID");
                entity(id, "compound_record", id);
                mapping(id, "inchikey:" + rows.at(r).at("InChIKey").text(),
                        "/PropertyTable/Properties/" + std::to_string(r), rows.at(r));
            }
        } else if (provider == "pdb" && resource == "polymer") {
            const auto left = "pdb:polymer:" + raw.at("rcsb_id").text();
            entity(left, "structure_polymer_entity", left);
            const auto container = raw.at("rcsb_polymer_entity_container_identifiers");
            if (container.has("reference_sequence_identifiers") &&
                container.at("reference_sequence_identifiers").raw()) {
                const auto refs = container.at("reference_sequence_identifiers");
                for (std::size_t r = 0; r < refs.size(); ++r)
                    if (string_or_empty(refs.at(r), "database_name") == "UniProt")
                        mapping(left, "uniprot:" + refs.at(r).at("database_accession").text(),
                                "/rcsb_polymer_entity_container_identifiers/"
                                "reference_sequence_identifiers/" +
                                    std::to_string(r),
                                refs.at(r));
            }
        } else if (provider == "pdb" && resource == "entry") {
            const auto left = "pdb:" + raw.at("rcsb_id").text();
            entity(left, "experimental_structure_entry", left);
            const auto ids_ = raw.at("rcsb_entry_container_identifiers").at("polymer_entity_ids");
            for (std::size_t r = 0; r < ids_.size(); ++r)
                mapping(left, "pdb:polymer:" + raw.at("rcsb_id").text() + "_" + ids_.at(r).text(),
                        "/rcsb_entry_container_identifiers/polymer_entity_ids/" + std::to_string(r),
                        raw);
        } else
            require(provider == "chembl" && (resource == "assay" || resource == "status"),
                    "no admitted evidence adapter for response");
    }
    require(p.records.size() > 0, "no admitted biomedical assertions");
    return p;
}
Packet normalize(const fs::path &input_directory,std::uint32_t version=2) {
    const auto directory=fs::canonical(input_directory);
    return normalize(capture_database_snapshot(directory),version,directory);
}
InterpretedSnapshot interpret(const Packet &p) {
    LogicalSnapshot logical;
    logical.revision = p.manifest.at("snapshot_sha256").text();
    if (p.adapter_version >= 2)
        logical.revision += "/adapter/" + std::to_string(p.adapter_version);
    for (const auto &[id, j] : p.entities) {
        (void)j;
        logical.subjects.push_back({1, id});
    }
    for (std::size_t i = 0; i < p.records.size(); ++i) {
        const auto row = p.records.at(i);
        const auto provider = row.at("provider").text();
        const std::uint64_t pc = provider == "chembl"      ? 1
                                 : provider == "bindingdb" ? 2
                                 : provider == "uniprot"   ? 3
                                 : provider == "pubchem"   ? 4
                                                           : 5;
        auto value = row.at("value_token").text().empty()
                         ? V::absent()
                         : V::present(V::rational(exact_decimal(row.at("value_token").text())));
        ExternalAssertion a;
        a.origin = {static_cast<std::uint64_t>(row.at("kind").integer()), row.at("id").text()};
        a.payload = V::product({value, V::symbol(code(relations, row.at("relation").text())),
                                V::symbol(code(endpoints, row.at("endpoint").text())),
                                V::symbol(code(units, row.at("units").text())), V::symbol(pc)});
        std::uint64_t role = 11;
        for (const auto *name : {"compound", "target", "assay", "document", "response"}) {
            const auto id = row.at(name).text();
            ExternalRoleValue r;
            r.role = role++;
            r.state = id.empty() ? ReferenceState::Missing : ReferenceState::Present;
            if (!id.empty())
                r.occurrences.push_back({name, {1, id}, true});
            a.roles.push_back(std::move(r));
        }
        logical.assertions.push_back(std::move(a));
    }
    // These are complete in-memory producer declarations. Admission uses the
    // addressable native containers and byte representation, not the defaults
    // intended for bounded external decoders. Intrinsic declaration depth and
    // tensor-rank representations remain the canonical SDK's responsibility.
    InterpretationLimits limits;
    const auto addressable = std::numeric_limits<std::size_t>::max();
    limits.max_records = addressable;
    limits.max_origin_bytes = source_byte_limit();
    limits.max_metadata_bytes = source_byte_limit();
    limits.max_archive_bytes = source_byte_limit();
    limits.max_work = addressable;
    limits.declarations.max_components = addressable;
    limits.declarations.max_value_nodes = addressable;
    limits.declarations.max_encoded_bytes = source_byte_limit();
    limits.declarations.max_integer_bits = addressable;
    limits.declarations.max_materialized_bytes = source_byte_limit();
    limits.atlas.accession = {addressable, addressable, addressable, addressable, addressable,
                              addressable, addressable, addressable, addressable, addressable};
    return InterpretedSnapshot::interpret(schema(p.adapter_version), logical, limits);
}
std::string archive(const InterpretedSnapshot &native) {
    const auto bytes = encode_interpreted_snapshot(native);
    return {reinterpret_cast<const char *>(bytes.data()), bytes.size()};
}
std::string reported_point_quantity(const std::string &endpoint) {
    return "reported standard " + endpoint +
           "; exact point (=) only; upstream normalization retained; "
           "missing/censored/other units excluded from this field";
}
void scientific_view(const Packet &p, const InterpretedSnapshot &native, const fs::path &root) {
    std::map<std::string, std::map<std::string, std::size_t>> categories;
    const std::vector<std::string> cats{"compound", "target",   "assay", "endpoint",
                                        "relation", "provider", "units", "document"};
    std::vector<Json> rows;
    for (std::size_t i = 0; i < p.records.size(); ++i)
        if (p.records.at(i).at("kind").integer() == 10 &&
            p.records.at(i).at("provider").text() == "chembl")
            rows.push_back(p.records.at(i));
    if (rows.empty())
        return;
    for (const auto &name : cats) {
        std::set<std::string> values;
        for (const auto &row : rows)
            if (!row.at(name).text().empty())
                values.insert(row.at(name).text());
        for (const auto &value : values)
            categories[name].emplace(value, categories[name].size() + 1);
    }
    std::string table, dictionary;
    auto columns = Json::array();
    const auto column = [&](const std::string &name, const std::string &channel,
                            const std::string &quantity, const std::string &unit,
                            const std::string &type) {
        const auto start = dictionary.size();
        dictionary += name + ": " + quantity + "; " + unit + ".\n";
        if (categories.contains(name))
            for (const auto &[id, n] : categories[name])
                dictionary += std::to_string(n) + " = " + id + "\n";
        auto j = Json::object();
        j.set("id", Json::string(name));
        j.set("channel", Json::string(channel));
        j.set("quantity", Json::string(quantity));
        j.set("units", Json::string(unit));
        j.set("value_type", Json::string(type));
        j.set("dictionary_start", cardinality(start));
        j.set("dictionary_end", cardinality(dictionary.size()));
        columns.push(j);
    };
    for (const auto &name : cats)
        column(name, "covariate",
               "reported source " + name + " identity; category code has no amplitude",
               "category_code", "categorical");
    for (const auto &e : std::vector<std::string>{"Kd", "Ki", "IC50", "EC50"})
        column(e + "_nM", "T",
               reported_point_quantity(e),
               "nM", "quantitative");
    auto mapping = Json::array();
    std::map<std::string, std::uint64_t> native_keys;
    for (const auto &entry : native.lineage().assertions)
        native_keys.emplace(entry.origin.token, entry.native_key);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto row = rows[i];
        for (const auto &name : cats) {
            if (!table.empty() && table.back() != '\n')
                table += '\t';
            table += row.at(name).text().empty()
                         ? "?"
                         : std::to_string(categories.at(name).at(row.at(name).text()));
        }
        for (const auto &e : std::vector<std::string>{"Kd", "Ki", "IC50", "EC50"}) {
            table += '\t';
            const auto point = row.at("relation").text() == "=" && row.at("endpoint").text() == e &&
                               row.at("units").text() == "nM" &&
                               !row.at("value_token").text().empty();
            table += point ? row.at("value_token").text() : "?";
        }
        table += '\n';
        auto item = Json::object();
        item.set("scientific_record_key", Json::string(std::to_string(i + 1)));
        item.set("assertion_key",
                 Json::string(std::to_string(native_keys.at(row.at("id").text()))));
        item.set("record_id", row.at("id"));
        mapping.push(item);
    }
    auto spec = Json::object();
    spec.set("schema", Json::string("gadget.table.import.v1"));
    spec.set("id", Json::string(p.manifest.at("id").text() + "-binding-observations"));
    spec.set("source_namespace", Json::string("72201"));
    auto file = Json::object();
    file.set("path", Json::string("observations.tsv"));
    file.set("sha256", Json::string("sha256:" + sha256(table)));
    file.set("delimiter", Json::string("\t"));
    file.set("header", Json::boolean(false));
    spec.set("source", file);
    file = Json::object();
    file.set("path", Json::string("dictionary.txt"));
    file.set("sha256", Json::string("sha256:" + sha256(dictionary)));
    spec.set("dictionary", file);
    spec.set("columns", columns);
    spec.set("missing_tokens", Json::parse("[\"?\"]"));
    auto rel = Json::array();
    for (const auto &name : cats)
        rel.push(Json::string(name));
    spec.set("relation_columns", rel);
    spec.set("question", Json::string("Describe retained reported binding measurements and their "
                                      "assay/target/source strata."));
    spec.set("population",
             Json::string("ChEMBL activity records in explicit frozen API pages; each row is a "
                          "source record, not an independent experiment or patient."));
    spec.set(
        "selection",
        Json::string("Exact saved query/page selection; BindingDB compact results remain in native "
                     "evidence store and do not enter this assay-complete numeric view."));
    spec.set("limitations",
             Json::parse(
                 "[\"No pooled affinity, efficacy, mechanism, causal, significance or population "
                 "claim.\",\"Kd, Ki, IC50 and EC50 remain distinct; numeric fields admit only "
                 "reported equality in nM, not censored thresholds.\",\"Upstream "
                 "duplicate/validity flags and assay conditions stay in original records; "
                 "different assays, constructs and species are not automatically comparable.\"]"));
    spec.set("source_url", Json::string(p.manifest.at("snapshot_sha256").text()));
    spec.set("attribution", Json::string("ChEMBL; original source documents and API metadata "
                                         "retained in adjacent evidence snapshot."));
    spec.set("license",
             Json::string("Source data terms retained per response; ChEMBL CC BY-SA 3.0."));
    write_text(root / "observations.tsv", table);
    write_text(root / "dictionary.txt", dictionary);
    write_text(root / "table-import.json", spec.dump() + "\n");
    write_text(root / "record-map.json", mapping.dump() + "\n");
    import_table(root / "table-import.json", root / "study");
    // Bind the mathematical view to the complete immutable source packet and
    // explicit native-assertion correspondence, without substituting its carrier.
    auto study = Json::parse(read_text(root / "study/study.json"));
    auto docs = study.at("documents");
    const auto attach = [&](const std::string &id, const std::string &name,
                            const std::string &bytes) {
        write_text(root / "study" / name, bytes);
        auto d = Json::object();
        d.set("id", Json::string(id));
        d.set("path", Json::string(name));
        d.set("format", Json::string("utf8"));
        d.set("sha256", Json::string("sha256:" + sha256(bytes)));
        docs.push(d);
    };
    attach("database_snapshot", "database-snapshot.json", p.manifest.dump() + "\n");
    attach("native_record_correspondence", "record-map.json", mapping.dump() + "\n");
    study.set("documents", docs);
    write_text(root / "study/study.json", study.dump() + "\n");
    ScientificStudy admitted(root / "study/study.json");
    write_text(root / "study/admission.json", admitted.catalog().dump() + "\n");
    rfc::storage::write_file_exact(root/"study/study.rfcxs",admitted.native_publication().encode(),rfc::storage::FilePublication::ReplaceRegular);
}
NativeRecord evidence_publication(const fs::path &snapshot,const Packet &packet,
                                  const InterpretedSnapshot &native,const Json &receipt) {
    using N=NativeRecordValue;auto database=capture_database_snapshot(snapshot);
    auto bytes=database.encode();const auto count=bytes.size();
    auto archive_bytes=rfc::database::encode_interpreted_snapshot(native);const auto archive_count=archive_bytes.size();
    std::vector<RawByteSource> raw;raw.push_back(RawByteSource::own(1,std::move(bytes)));
    raw.push_back(RawByteSource::own(2,std::move(archive_bytes)));
    N::Array observations;N::Object entities;
    std::map<std::string,std::uint64_t> assertion_keys,subject_keys;
    for(const auto &entry:native.lineage().assertions)assertion_keys.emplace(entry.origin.token,entry.native_key);
    for(const auto &entry:native.lineage().subjects)subject_keys.emplace(entry.origin.token,entry.native_key);
    for(std::size_t i=0;i<packet.records.size();++i) {
        const auto row=packet.records.at(i);
        observations.elements.emplace_back(N::Owned{0,rfc::storage::StableCellKey{rfc::core::Grade{1},assertion_keys.at(row.at("id").text())},{},std::make_shared<const N>(import_record_view(row))});
    }
    for(const auto &[id,value]:packet.entities)entities.members.emplace_back(id,N(N::Owned{0,rfc::storage::StableCellKey{rfc::core::Grade{0},subject_keys.at(id)},{},std::make_shared<const N>(import_record_view(value))}));
    const auto endpoint=native.endpoint();
    return NativeRecord::create(N(N::Object{{{"schema",N(std::string("gadget.evidence.native.v1"))},
        {"receipt",import_record_view(receipt)},{"snapshot",import_record_view(packet.manifest)},
        {"records",N(std::move(observations))},{"entities",N(std::move(entities))},
        {"source_database",N(N::ByteSpan{1,0,count,false})},
        {"interpretation",N(N::ByteSpan{2,0,archive_count,false})}}}),source_namespace,std::move(raw),endpoint.get());
}
} // namespace
void import_evidence(const fs::path &input_snapshot, const fs::path &destination,
                     std::uint32_t version) {
    const auto snapshot = fs::canonical(input_snapshot);
    require(!fs::exists(destination), "evidence destination already exists");
    const auto original = normalize(snapshot, version);
    const auto native = interpret(original);
    const auto bytes = archive(native);
    fs::create_directories(destination.parent_path());
    require(fs::create_directory(destination), "evidence creation conflict");
    try {
        fs::create_directory(destination / "snapshot");
        const auto database=capture_database_snapshot(snapshot);database.write(destination/"snapshot/snapshot.rfcxs");
        write_text(destination / "snapshot/snapshot.json", render_record_view(database.at("snapshot")).dump()+"\n");
        write_text(destination / "evidence.rfcxi", bytes);
        auto receipt = Json::object();
        receipt.set("schema", Json::string("gadget.biomedical.store.v1"));
        receipt.set("adapter_version", Json::integer(version));
        receipt.set("snapshot_sha256", original.manifest.at("snapshot_sha256"));
        receipt.set("native_sha256", Json::string("sha256:" + sha256(bytes)));
        receipt.set("records_sha256", Json::string("sha256:" + sha256(original.records.dump())));
        receipt.set("assertions", cardinality(original.records.size()));
        receipt.set("entities", cardinality(original.entities.size()));
        receipt = seal(receipt, "source_digest");
        write_text(destination / "store.json", receipt.dump() + "\n");
        evidence_publication(snapshot,original,native,receipt).write(destination/"store.rfcxs");
        scientific_view(original, native, destination);
    } catch (...) {
        fs::remove_all(destination);
        throw;
    }
}
struct EvidenceStore::State {
    NativeRecord publication;
    RetainedRecord receipt,manifest;
    InterpretedSnapshot native;
    std::map<std::uint64_t,RetainedRecord> records,entities;
    static NativeRecord open(const fs::path &directory) {
        if(fs::is_regular_file(directory/"store.rfcxs"))return NativeRecord::read(directory/"store.rfcxs");
        const auto receipt=Json::parse(read_text(directory/"store.json"));
        const auto packet=normalize(directory/"snapshot",receipt.has("adapter_version")?static_cast<std::uint32_t>(count(receipt.at("adapter_version"),3,"adapter version")):1);
        const auto native=interpret(packet);return evidence_publication(directory/"snapshot",packet,native,receipt);
    }
    static InterpretedSnapshot interpret_native(const NativeRecord &publication) {
        auto limits=InterpretationLimits{};const auto addressable=std::numeric_limits<std::size_t>::max();
        limits.max_records=limits.max_origin_bytes=limits.max_metadata_bytes=limits.max_archive_bytes=limits.max_work=addressable;
        limits.declarations=native_data_limits().temporal.declarations;
        limits.atlas.accession=native_accession_limits();
        return decode_interpreted_snapshot(publication.at("interpretation").bytes(),limits);
    }
    explicit State(const fs::path &directory)
        :publication(open(directory)),receipt(publication.at("receipt")),manifest(publication.at("snapshot")),native(interpret_native(publication)) {
        require(publication.at("schema").text()=="gadget.evidence.native.v1","unsupported native biomedical evidence");
        verify_seal(receipt,"source_digest");
        const auto normalized=publication.at("records");
        require(receipt.at("snapshot_sha256").text()==manifest.at("snapshot_sha256").text() && receipt.at("records_sha256").text()=="sha256:"+sha256(render_record_view(normalized).dump()),"evidence source/records mismatch");
        const auto bytes=publication.at("interpretation").bytes();
        require(receipt.at("native_sha256").text()=="sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size())),"native evidence/source mismatch");
        if(fs::exists(directory/"evidence.rfcxi"))require(read_text(directory/"evidence.rfcxi")==archive(native),"external native evidence differs from publication");
        const auto source_series=rfc::storage::decode_native_series(publication.at("source_database").bytes(),native_data_limits());
        const auto source=NativeRecord::from_history(source_series);
        const auto expected=normalize(source,static_cast<std::uint32_t>(count(receipt.at("adapter_version"),3,"adapter version")));
        auto expected_entities=Json::object();for(const auto &[id,value]:expected.entities)expected_entities.set(id,value);
        require(expected.manifest.dump()==manifest.dump() && expected.records.dump()==render_record_view(normalized).dump() && expected_entities.dump()==render_record_view(publication.at("entities")).dump(),"biomedical declarations differ from retained original responses");
        require(archive(interpret(expected))==archive(native),"native interpretation differs from retained original responses");
        require(rfc::storage::encode_rfcx(publication.endpoint()->source).bytes==rfc::storage::encode_rfcx(native.endpoint()->source).bytes,"biomedical complex differs from retained native interpretation");
        const auto &declared=*publication.endpoint()->declarations;
        for(const auto &definition:native.endpoint()->declarations->registry.definitions())
            require(declared.registry.at(definition.id)==definition,"original biomedical type was altered");
        std::map<std::uint64_t,const rfc::storage::TemporalDeclaredValue *> original_rows;
        for(const auto &row:declared.values)original_rows.emplace(row.key,&row);
        for(const auto &original:native.endpoint()->declarations->values) {
            require(original_rows.contains(original.key),"original biomedical payload is missing");
            const auto &row=*original_rows.at(original.key);
            require(row.type==original.type && rfc::storage::temporal_declaration_detail::owner_equal(row,original) && rfc::composite::declared_value_equal(row.value,original.value),"original biomedical payload was altered");
        }
        std::map<std::string,NativeRecord> by_id;
        for(std::size_t i=0;i<normalized.size();++i)require(by_id.emplace(normalized.at(i).at("id").text(),normalized.at(i)).second,"duplicate observation identity");
        for(const auto &entry:native.lineage().assertions)records.emplace(entry.native_key,by_id.at(entry.origin.token));
        for(const auto &entry:native.lineage().subjects)entities.emplace(entry.native_key,publication.at("entities").at(entry.origin.token));
    }
    static ResultValue selected_value(const NativeRecord &record, const std::string &key_name,
                                      std::uint64_t k,
                                      const std::optional<std::string> &criterion = {}) {
        using N = NativeRecordValue;
        auto value = record.materialize_value();
        auto *owner = std::get_if<N::Owned>(&value.data);
        auto payload = owner ? ResultValue(*owner->payload) : ResultValue(value);
        payload.set(key_name, ResultValue::string(std::to_string(k)));
        if (criterion) payload.set("criterion", ResultValue::string(*criterion));
        if (!owner) return payload;
        owner->payload = std::make_shared<const N>(std::move(payload).release());
        return ResultValue(std::move(value));
    }
    ResultValue record(std::uint64_t k, const std::optional<std::string> &criterion = {}) const {
        return selected_value(records.at(k).native(), "native_assertion_key", k, criterion);
    }
};
EvidenceStore::EvidenceStore(const fs::path &directory)
    : state_(std::make_unique<State>(directory)) {}
EvidenceStore::~EvidenceStore() = default;
const NativeRecord &EvidenceStore::native_publication() const {return state_->publication;}
std::string EvidenceStore::digest() const { return state_->receipt.at("source_digest").text(); }
NativeRecord EvidenceStore::bind_reported_field(const ScientificStudy &study,
                                               const NativeRecord &request) const {
    using N = NativeRecordValue;
    const auto control = render_record_view(request);
    members(control, {"schema", "source_digest", "study_digest", "field", "endpoint", "units"});
    require(control.at("schema").text() == "gadget.biomedical.field.request.v1" &&
                control.at("source_digest").text() == digest() &&
                control.at("study_digest").text() == study.digest(),
            "biomedical field source/study/schema mismatch");
    const auto endpoint = text(control.at("endpoint")), unit = text(control.at("units"));
    require(endpoint == "Kd" || endpoint == "Ki" || endpoint == "IC50" || endpoint == "EC50",
            "reported field requires an explicit binding endpoint");
    require(unit == "nM", "reported field projection requires its declared nM units");
    const auto field = text(control.at("field"));
    const auto declaration = study.field_declaration(field);
    require(declaration.at("value_type").text() == "quantitative" &&
                declaration.at("units").text() == unit &&
                declaration.at("quantity").text() == reported_point_quantity(endpoint),
            "reported field endpoint/units/quantity declaration mismatch");
    const auto &publication = study.native_publication();
    const auto document = [&](const std::string &id) {
        const auto documents = publication.at("documents");
        for (std::size_t i = 0; i < documents.size(); ++i) {
            const auto item = documents.at(i);
            if (item.at("id").text() == id) {
                const auto bytes = item.at("bytes").bytes();
                return Json::parse(std::string(reinterpret_cast<const char *>(bytes.data()),
                                               bytes.size()));
            }
        }
        throw std::invalid_argument("reported field source document is absent");
    };
    require(document("database_snapshot").dump() == state_->manifest.dump(),
            "reported field snapshot differs from evidence owner");
    const auto mapping = document("native_record_correspondence");
    source_array(mapping, "reported field correspondence");
    std::vector<std::uint64_t> original_order;
    const auto original = state_->publication.at("records");
    std::map<std::string, std::uint64_t> assertion_ids;
    for (const auto &[assertion, row] : state_->records)
        assertion_ids.emplace(row.at("id").text(), assertion);
    for (std::size_t i = 0; i < original.size(); ++i) {
        const auto row = original.at(i);
        if (row.at("kind").integer() == 10 && row.at("provider").text() == "chembl")
            original_order.push_back(assertion_ids.at(row.at("id").text()));
    }
    const auto pin = study.field_accession(field);
    const auto tensor = pin.tensor_declaration();
    const auto carrier = pin.carrier();
    require(tensor->role == rfc::storage::TemporalTensorRole::Field && !tensor->input &&
                tensor->output && *tensor->output == carrier->key &&
                carrier->kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                carrier->grade == rfc::core::Grade{0} &&
                tensor->tensor.domain() == carrier->domain &&
                mapping.size() == original_order.size() &&
                mapping.size() == tensor->tensor.domain().logical_size(),
            "reported field correspondence requires the complete original record carrier");
    N::Array correspondence;
    std::size_t supported = 0;
    for (std::size_t i = 0; i < mapping.size(); ++i) {
        const auto entry = mapping.at(i);
        members(entry, {"scientific_record_key", "assertion_key", "record_id"});
        const auto assertion = key(entry.at("assertion_key"));
        require(assertion == original_order.at(i) &&
                    key(entry.at("scientific_record_key")) ==
                        rfc::storage::stable_cell_key(pin.endpoint()->source,
                            {{0}, static_cast<std::uint32_t>(i)}).key,
                "reported field correspondence changed original order or identity");
        const auto row = state_->records.at(assertion);
        require(row.at("id").text() == entry.at("record_id").text(),
                "reported field record identity mismatch");
        const auto payload = state_->native.payload(assertion, {1}).declared_value();
        const auto &value = std::get<V::Optional>(payload->data).element;
        const bool point = row.at("relation").text() == "=" &&
                           row.at("endpoint").text() == endpoint &&
                           row.at("units").text() == unit && !value.empty();
        const auto actual = tensor->tensor.state(i);
        require(actual.has_value() == point &&
                    (!point || *actual == std::get<Rational>(value.at(0).data)),
                "reported field native support/value differs from original assertion");
        supported += point;
        correspondence.elements.emplace_back(N::Object{{
            {"coordinate", N(rfc::numeric::Integer(std::to_string(i)))},
            {"scientific_record_key", N(entry.at("scientific_record_key").text())},
            {"assertion_key", N(std::to_string(assertion))},
            {"record", state_->record(assertion).release()}, {"supported", N(point)}}});
    }
    const auto number = [](std::uint64_t value) {
        return N(rfc::numeric::Integer(std::to_string(value)));
    };
    const auto nullable = [&](const std::optional<std::uint64_t> &value) {
        return value ? number(*value) : N{};
    };
    // Separate archives preserve the study declarations and actual field cutoff.
    // Their raw source keys do not establish a map to the evidence complex.
    std::vector<RawByteSource> sources;
    std::uint64_t next = 1;
    for (const auto source_key : NativeRecord::source_keys(*state_->publication.endpoint()))
        next = rfc::machine::checked_size_add(std::max(next, source_key), 1);
    const auto study_key = next, field_key = rfc::machine::checked_size_add(next, 1);
    sources.push_back(publication.archive_source(study_key));
    sources.push_back(publication.historical_publication(pin.sequence()).archive_source(field_key));
    const auto study_bytes = sources.at(0).bytes().size(), field_bytes = sources.at(1).bytes().size();
    N selector(N::Object{{{"schema", N(std::string("gadget.rcdb.field.query.v1"))},
        {"source_namespace", N(std::to_string(pin.source_namespace()))},
        {"cutoff", number(pin.sequence())}, {"tensor", number(tensor->key)},
        {"input_carrier", nullable(tensor->input)}, {"output_carrier", nullable(tensor->output)},
        {"interpretation", N(tensor->interpretation)}}});
    auto source = state_->publication.derive(state_->publication.materialize_value(),
                                             std::move(sources));
    return retain_operation(source, request, N(N::Object{{
        {"schema", N(std::string("gadget.biomedical.field.binding.v1"))},
        {"source_digest", N(digest())}, {"study_digest", N(study.digest())},
        {"field", N(field)}, {"endpoint", N(endpoint)}, {"units", N(unit)},
        {"relation", N(std::string("="))}, {"supported", number(supported)},
        {"selector", std::move(selector)}, {"correspondence", N(std::move(correspondence))},
        {"study_publication", N(N::ByteSpan{study_key, 0, study_bytes, false})},
        {"field_publication", N(N::ByteSpan{field_key, 0, field_bytes, false})},
        {"scope", N(std::string("Complete retained assay-identified ChEMBL record projection; "
            "reported equality points in declared nM units. Separate sources, assays and "
            "endpoints do not imply biological comparability."))}}}));
}
NativeRecord EvidenceStore::source_context(const NativeRecord &request) const {
    using N = NativeRecordValue;
    const auto control = render_record_view(request);
    members(control, {"schema", "source_digest", "snapshot_sha256", "question", "assertion_keys"});
    require(control.at("schema").text() == "gadget.biomedical.context.v1" &&
                text(control.at("source_digest")) == digest() &&
                text(control.at("snapshot_sha256")) == state_->manifest.at("snapshot_sha256").text(),
            "biomedical context source/snapshot/schema mismatch");
    const auto question = text(control.at("question"));
    const auto keys = control.at("assertion_keys");
    source_array(keys, "biomedical assertion selection");
    require(keys.size() > 0, "biomedical context requires selected assertions");
    N::Array records;
    std::set<std::uint64_t> selected;
    for (std::size_t i = 0; i < keys.size(); ++i) {
        const auto k = key(keys.at(i));
        require(selected.insert(k).second && state_->records.contains(k),
                "biomedical context assertion is unknown or repeated");
        records.elements.emplace_back(N::Object{
            {{"assertion_key", N(std::to_string(k))},
             {"record", state_->records.at(k).native().materialize_value()}}});
    }
    N::Object context{{{"schema", N(std::string("gadget.biomedical.context.catalog.v1"))},
        {"source_digest", N(digest())},
        {"snapshot_sha256", state_->manifest.native().at("snapshot_sha256").materialize_value()},
        {"question", N(question)},
        {"question_sha256", N("sha256:" + sha256(question))},
        {"context_sha256", N("sha256:" + sha256(digest() + "\n" + native_declaration_digest(request)))},
        {"availability", N(N::Object{{{"kind", N(std::string("frozen_snapshot"))},
            {"historical_availability_established", N(false)}}})},
        {"records", N(std::move(records))}}};
    return retain_operation(state_->publication, request, N(std::move(context)));
}
NativeRecord EvidenceStore::admit_proposal(const NativeRecord &context_request,
                                         const NativeRecord &input,
                                         const std::string &question) const {
    using N = NativeRecordValue;
    const auto context = source_context(context_request);
    const auto catalog = context.at("result");
    require(question == catalog.at("question").text(), "biomedical proposal question mismatch");
    const auto proposal = render_record_view(input);
    members(proposal, {"schema", "source_digest", "context_sha256", "question_sha256", "action"},
            {"query", "mentions", "reason"});
    require(proposal.at("schema").text() == "gadget.biomedical.proposal.v1" &&
                text(proposal.at("source_digest")) == digest() &&
                text(proposal.at("context_sha256")) == catalog.at("context_sha256").text() &&
                text(proposal.at("question_sha256")) == catalog.at("question_sha256").text(),
            "biomedical proposal source/context/question/schema mismatch");
    std::set<std::uint64_t> selected;
    for (std::size_t i = 0; i < catalog.at("records").size(); ++i)
        selected.insert(key(render_record_view(catalog.at("records").at(i).at("assertion_key"))));
    const auto action = text(proposal.at("action"));
    auto mentions = Json::array();
    if (action == "clarify") {
        require(proposal.has("reason") && !proposal.has("query") && !proposal.has("mentions"),
                "biomedical clarification requires only a reason");
        (void)text(proposal.at("reason"));
    } else if (action == "query") {
        require(proposal.has("query") && !proposal.has("reason") && !proposal.has("mentions"),
                "biomedical query proposal requires only a query");
        const auto query = proposal.at("query");
        validate(query);
        const auto operation = query.at("operation").text();
        if (operation == "compare")
            require(selected.contains(key(query.at("left"))) &&
                        selected.contains(key(query.at("right"))),
                    "biomedical comparison lies outside the selected context");
        if (operation == "role_join") {
            const auto join = role_join_request(query,state_->native);
            for (const auto kind : {join.left_kind,join.right_kind})
                for (const auto assertion : state_->native.assertions(kind))
                    require(selected.contains(assertion),
                            "biomedical role join requires complete selected source kinds");
        }
        if (operation == "neighbors") {
            const auto subject = key(query.at("subject_key"));
            const auto id = state_->entities.at(subject).at("id").text();
            bool found = false;
            for (const auto k : selected)
                for (const auto *role : {"compound", "target", "assay", "document", "response"})
                    found = found || state_->records.at(k).at(role).text() == id;
            require(found, "biomedical subject lies outside the selected context");
        }
        if (operation == "affinity_evidence" || operation == "select") {
            bool found = false;
            for (const auto k : selected) {
                const auto &record = state_->records.at(k);
                if (operation == "affinity_evidence")
                    found = found || (record.at("endpoint").text() == query.at("endpoint").text() &&
                                      record.at("units").text() == query.at("units").text());
                else
                    found = found || record.at(query.at("field").text()).text() ==
                                         query.at("value").text();
            }
            require(found, "biomedical query predicates lack selected source evidence");
        }
    } else {
        require(action == "annotate" && proposal.has("mentions") && !proposal.has("query") &&
                    !proposal.has("reason"), "unknown biomedical proposal action");
        source_array(proposal.at("mentions"), "biomedical annotation mentions");
        require(proposal.at("mentions").size() > 0, "biomedical annotation requires mentions");
        const auto database = NativeRecord::decode(state_->publication.at("source_database").bytes());
        std::map<std::string, Document> documents;
        for (std::size_t i = 0; i < database.at("snapshot").at("records").size(); ++i) {
            const auto record = database.at("snapshot").at("records").at(i);
            const auto artifacts = database.at("artifacts");
            for (std::size_t j = 0; j < artifacts.size(); ++j)
                if (artifacts.at(j).at("path").text() == record.at("path").text()) {
                    const auto bytes = artifacts.at(j).at("bytes").bytes();
                    auto metadata = Json::object();
                    metadata.set("format", Json::string("utf8"));
                    metadata.set("sha256", render_record_view(record.at("sha256")));
                    documents.emplace(record.at("request").at("id").text(), Document{
                        metadata, {reinterpret_cast<const char *>(bytes.data()), bytes.size()}});
                }
        }
        std::set<std::string> ids;
        for (std::size_t i = 0; i < proposal.at("mentions").size(); ++i) {
            auto mention = Json::parse(proposal.at("mentions").at(i).dump());
            members(mention, {"id", "kind", "assertion_key", "polarity", "interpretation", "evidence"},
                    {"role", "entity", "endpoint", "units", "relation", "value_token"});
            require(ids.insert(text(mention.at("id"))).second, "duplicate biomedical mention");
            const auto k = key(mention.at("assertion_key"));
            require(selected.contains(k), "annotation assertion lies outside the selected context");
            const auto &record = state_->records.at(k);
            const auto kind = text(mention.at("kind"));
            if (kind == "entity") {
                require(mention.has("role") && mention.has("entity") && !mention.has("endpoint") &&
                            !mention.has("units") && !mention.has("relation") && !mention.has("value_token"),
                        "entity annotation requires only a role and identifier");
                const auto role = text(mention.at("role"));
                require(role == "compound" || role == "target" || role == "assay" ||
                            role == "document" || role == "response", "unknown biomedical entity role");
                require(text(mention.at("entity")) == record.at(role).text(),
                        "biomedical annotation identifier differs from the original assertion role");
            } else {
                require(kind == "measurement" && record.at("kind").integer() == 10 &&
                            !mention.has("role") && !mention.has("entity"),
                        "measurement annotation requires a binding observation");
                for (const auto *name : {"endpoint", "units", "relation", "value_token"})
                    require(mention.has(name) && mention.at(name).text() == record.at(name).text(),
                            "measurement annotation changes original endpoint/units/relation/value");
            }
            const auto polarity = text(mention.at("polarity"));
            require(polarity == "asserted" || polarity == "negated" || polarity == "uncertain",
                    "biomedical annotation polarity must be explicit");
            (void)text(mention.at("interpretation"));
            const auto evidence = source_evidence(mention.at("evidence"), documents, true);
            for (std::size_t j = 0; j < evidence.size(); ++j)
                require(evidence.at(j).at("document").text() ==
                            record.at("source_receipt").at("request").at("id").text(),
                        "biomedical annotation quote belongs to another response");
            mention.set("evidence", evidence);
            mention.set("interpretation_verified", Json::boolean(false));
            mentions.push(mention);
        }
    }
    N::Object receipt{{{"schema", N(std::string("gadget.biomedical.proposal.receipt.v1"))},
        {"status", N(std::string(action == "clarify" ? "clarification" : "prospective"))},
        {"context", catalog.materialize_value()},
        {"mentions", import_record_view(mentions)},
        {"scope", N(std::string("Source, question, identifiers, literal measurements and quotes admitted. "
            "Interpretations remain unverified; queries require separate execution. "
            "This snapshot establishes no historical availability or biological observation time."))}}};
    return retain_operation(context, input, N(std::move(receipt)));
}
void EvidenceStore::validate(const Json &query) const {
    const auto op = query.at("operation").text();
    if (op == "inspect" || op == "coverage")
        members(query, {"operation"});
    else if (op == "select") {
        members(query, {"operation", "field", "comparison", "value"});
        require(query.at("comparison").text() == "equal", "categorical predicate has no order");
        const auto field = query.at("field").text(), value = query.at("value").text();
        require(field == "endpoint" || field == "relation" || field == "units" ||
                    field == "provider",
                "numeric predicates need a separate assay/unit admission profile");
        if (field == "provider")
            require(value == "chembl" || value == "bindingdb", "unknown binding provider");
        else {
            const auto &dictionary = field == "endpoint"   ? endpoints
                                     : field == "relation" ? relations
                                                           : units;
            require(std::find(dictionary.begin(), dictionary.end(), value) != dictionary.end(),
                    "unknown finite predicate value");
        }
    } else if (op == "neighbors") {
        members(query, {"operation", "subject_key"});
        require(state_->entities.contains(key(query.at("subject_key"))), "unknown native subject");
    } else if (op == "role_join") {
        members(query,{"operation","source_digest","left_kind","left_role","right_kind","right_role","activity","max_matches"});
        require(query.at("source_digest").text()==digest(), "biomedical join source digest mismatch");
        (void)role_join_request(query,state_->native);
        (void)count(query.at("max_matches"),std::numeric_limits<std::size_t>::max(),"role join result matches");
    } else if (op == "affinity_evidence") {
        members(query, {"operation", "endpoint", "units", "maximum"});
        const auto endpoint = query.at("endpoint").text(), unit = query.at("units").text();
        require(endpoint == "Kd" || endpoint == "Ki" || endpoint == "IC50" || endpoint == "EC50",
                "threshold requires an explicit binding endpoint");
        require(unit == "nM" || unit == "uM" || unit == "M",
                "threshold requires an explicit molar unit");
        const auto maximum = exact_decimal(text(query.at("maximum")));
        require(maximum >= 0, "affinity threshold must be nonnegative");
    } else if (op == "compare") {
        members(query, {"operation", "left", "right"});
        const auto l = key(query.at("left")), r = key(query.at("right"));
        require(l != r && state_->records.contains(l) && state_->records.contains(r),
                "comparison requires two admitted assertion keys");
        const auto a = state_->records.at(l), b = state_->records.at(r);
        for (const auto *name : {"provider", "compound", "target", "assay", "endpoint", "units"})
            require(a.at(name).text() == b.at(name).text(),
                    "comparison requires identical compound/target/assay/endpoint/units");
        require(a.at("kind").integer() == 10 && b.at("kind").integer() == 10 &&
                    !a.at("assay").text().empty() && a.at("relation").text() == "=" &&
                    b.at("relation").text() == "=" && !a.at("value_token").text().empty() &&
                    !b.at("value_token").text().empty(),
                "comparison requires assay-identified reported points");
    } else
        throw std::invalid_argument("unknown biomedical evidence query");
}
NativeRecord EvidenceStore::execute_record(const NativeRecord &input) const {
    using R = ResultValue;
    using N = NativeRecordValue;
    const auto number = [](Rational value) {
        const bool finite = std::isfinite(value.get_d());
        return R(N(N::ExactNumber{std::move(value), finite}));
    };
    const auto query = render_record_view(input);
    validate(query);
    const auto op = query.at("operation").text();
    auto result = R::object();
    if (op == "inspect") {
        members(query, {"operation"});
        result.set("assertions", R::count(state_->records.size()));
        result.set("selected", R::count(state_->records.size()));
        auto entities = R::array();
        for (const auto &[k, j] : state_->entities)
            entities.push(State::selected_value(j.native(), "native_subject_key", k));
        result.set("entities", entities);
        result.set("database_snapshot", R(state_->manifest.native().materialize_value()));
    } else if (op == "select") {
        members(query, {"operation", "field", "comparison", "value"});
        const auto field = query.at("field").text();
        ExactScalarPredicate p;
        if (field == "endpoint" || field == "relation" || field == "units" || field == "provider") {
            const auto &dictionary = field == "endpoint"   ? endpoints
                                     : field == "relation" ? relations
                                                           : units;
            std::uint64_t symbol = 0;
            if (field == "provider") {
                const auto v = query.at("value").text();
                require(v == "chembl" || v == "bindingdb", "binding provider must be explicit");
                symbol = v == "chembl" ? 1 : 2;
            } else {
                const auto v = query.at("value").text();
                require(std::find(dictionary.begin(), dictionary.end(), v) != dictionary.end(),
                        "unknown finite predicate value");
                symbol = code(dictionary, v);
            }
            require(query.at("comparison").text() == "equal", "categorical predicate has no order");
            p.factors = {field == "endpoint"   ? 3u
                         : field == "relation" ? 2u
                         : field == "units"    ? 4u
                                               : 5u};
            p.literal = V::symbol(symbol);
        } else
            throw std::invalid_argument(
                "numeric threshold predicates require an endpoint/unit/assay admission profile; "
                "this query selects declared categories only");
        const auto selected =
            select_exact(state_->native, 10, p, state_->native.assertions(10).size());
        auto rows = R::array();
        for (const auto k : selected)
            rows.push(state_->record(k));
        result.set("selected", R::count(selected.size()));
        result.set("records", rows);
    } else if (op == "neighbors") {
        members(query, {"operation", "subject_key"});
        const auto subject = key(query.at("subject_key"));
        require(state_->entities.contains(subject), "unknown native biomedical subject");
        auto rows = R::array();
        std::size_t selected = 0;
        for (const auto &kind : state_->native.schema().assertions)
            for (const auto assertion : state_->native.assertions(kind.key)) {
                bool found = false;
                for (const auto role : {11u, 12u, 13u, 14u, 15u})
                    for (const auto &occurrence : state_->native.participants(assertion, role))
                        found = found || occurrence.participant().key == subject;
                if (found) {
                    rows.push(state_->record(assertion));
                    ++selected;
                }
            }
        result.set("selected", R::count(selected));
        result.set("records", rows);
        result.set("subject", R(state_->entities.at(subject).native().materialize_value()));
    } else if (op == "role_join") {
        const auto request = role_join_request(query,state_->native);
        const auto lease = source_execution_lease("gadget-biomedical-native-role-join");
        const rfc::runtime::WorkAvailability availability{state_->native.schema().source_namespace,state_->native.sequence()};
        const auto prepared = PreparedNativeRoleJoin::prepare(state_->native,request,availability,lease);
        const auto maximum = count(query.at("max_matches"),std::numeric_limits<std::size_t>::max(),"role join result matches");
        if (prepared.match_count()>maximum)
            throw rfc::runtime::ExecutionLimitExceeded("role join result matches",prepared.match_count(),maximum);
        const auto joined = prepared.execute(availability,lease);
        auto rows = R::array();
        for (const auto &match : joined.matches()) {
            lease.check_cancelled();
            auto row = R::object();
            row.set("participant_key",R::string(std::to_string(match.participant.key)));
            row.set("participant",R::string(state_->entities.at(match.participant.key).at("id").text()));
            row.set("left_assertion",R::string(std::to_string(match.left_assertion)));
            row.set("right_assertion",R::string(std::to_string(match.right_assertion)));
            row.set("left_record_id",R::string(state_->records.at(match.left_assertion).at("id").text()));
            row.set("right_record_id",R::string(state_->records.at(match.right_assertion).at("id").text()));
            row.set("left_occurrence",R::string(std::to_string(match.left_occurrence.key.value)));
            row.set("right_occurrence",R::string(std::to_string(match.right_occurrence.key.value)));
            row.set("left_position",R::count(match.left_position));
            row.set("right_position",R::count(match.right_position));
            rows.push(std::move(row));
        }
        result.set("selected",R::count(joined.matches().size()));
        result.set("matches",std::move(rows));
        result.set("missing_role_fibers",R::count(prepared.missing_fibers()));
        result.set("null_role_fibers",R::count(prepared.null_fibers()));
        result.set("source_namespace",R::string(std::to_string(availability.source_namespace)));
        result.set("cutoff",R::count(availability.max_sequence));
        result.set("scope",R::string("Complete native role occurrence pairs on this exact retained source. Cross references do not assert identity equivalence, assay comparability, disease causality or efficacy."));
    } else if (op == "coverage") {
        members(query, {"operation"});
        const auto missing = select_exact(state_->native, 10, {{1}, ExactPredicateKind::IsNull, {}},
                                          state_->native.assertions(10).size());
        const auto points =
            select_exact(state_->native, 10, {{2}, ExactPredicateKind::Equal, V::symbol(1)},
                         state_->native.assertions(10).size());
        const auto all = state_->native.assertions(10).size();
        result.set("reported", R::count(all - missing.size()));
        result.set("missing", R::count(missing.size()));
        result.set("equality_relation_records", R::count(points.size()));
        result.set("binding_records", R::count(all));
        result.set(
            "scope",
            R::string(
                "Reported numerical tokens include censored bounds; equality count is a relation "
                "count and does not imply comparable assays or a present numeric value."));
    } else if (op == "affinity_evidence") {
        const auto maximum = exact_decimal(query.at("maximum").text());
        const auto as_set = [](std::vector<std::uint64_t> v) {
            return std::set<std::uint64_t>(v.begin(), v.end());
        };
        const auto endpoint =
            as_set(select_exact(state_->native, 10,
                                {{3},
                                 ExactPredicateKind::Equal,
                                 V::symbol(code(endpoints, query.at("endpoint").text()))},
                                state_->native.assertions(10).size()));
        const auto unit = as_set(select_exact(
            state_->native, 10,
            {{4}, ExactPredicateKind::Equal, V::symbol(code(units, query.at("units").text()))},
            state_->native.assertions(10).size()));
        const auto le = as_set(select_exact(
            state_->native, 10, {{1}, ExactPredicateKind::LessEqual, V::rational(maximum)},
            state_->native.assertions(10).size()));
        const auto ge = as_set(select_exact(
            state_->native, 10, {{1}, ExactPredicateKind::GreaterEqual, V::rational(maximum)},
            state_->native.assertions(10).size()));
        const auto gt = as_set(select_exact(
            state_->native, 10, {{1}, ExactPredicateKind::Greater, V::rational(maximum)},
            state_->native.assertions(10).size()));
        const auto missing =
            as_set(select_exact(state_->native, 10, {{1}, ExactPredicateKind::IsNull, {}},
                                state_->native.assertions(10).size()));
        std::size_t yes = 0, no = 0, unknown = 0;
        auto rows = R::array();
        for (const auto k : state_->native.assertions(10))
            if (endpoint.contains(k) && unit.contains(k)) {
                const auto relation = state_->records.at(k).at("relation").text();
                auto classification = std::string("indeterminate");
                if (!missing.contains(k)) {
                    if (relation == "=")
                        classification = le.contains(k) ? "entailed" : "contradicted";
                    else if ((relation == "<" || relation == "<=") && le.contains(k))
                        classification = "entailed";
                    else if ((relation == ">" && ge.contains(k)) ||
                             (relation == ">=" && gt.contains(k)))
                        classification = "contradicted";
                }
                if (classification == "entailed")
                    ++yes;
                else if (classification == "contradicted")
                    ++no;
                else
                    ++unknown;
                rows.push(state_->record(k, classification));
            }
        result.set("selected", R::count(yes));
        result.set("contradicted", R::count(no));
        result.set("indeterminate", R::count(unknown));
        result.set("eligible", R::count(yes + no + unknown));
        result.set("other_endpoint_or_units",
                   R::count(state_->native.assertions(10).size() - (yes + no + unknown)));
        result.set("maximum", number(maximum));
        result.set("records", rows);
        result.set(
            "scope",
            R::string("Per-record reported interval entails, contradicts or leaves unknown x <= "
                         "maximum under the explicit endpoint/unit. No pooled selectivity, "
                         "biological affinity certainty, efficacy or assay-comparability claim."));
    } else if (op == "compare") {
        members(query, {"operation", "left", "right"});
        const auto l = key(query.at("left")), r = key(query.at("right"));
        require(l != r, "comparison requires two distinct observations");
        const auto a = state_->records.at(l), b = state_->records.at(r);
        for (const auto *name : {"provider", "compound", "target", "assay", "endpoint", "units"})
            require(a.at(name).text() == b.at(name).text(),
                    "comparison requires identical declared compound/target/assay/endpoint/units");
        require(
            !a.at("assay").text().empty() && a.at("relation").text() == "=" &&
                b.at("relation").text() == "=" && !a.at("value_token").text().empty() &&
                !b.at("value_token").text().empty(),
            "comparison requires assay-identified reported points; censored/missing values reject");
        const auto av = state_->native.payload(l, {1}).declared_value(),
                   bv = state_->native.payload(r, {1}).declared_value();
        const auto x = std::get<Rational>(std::get<V::Optional>(av->data).element.at(0).data),
                   y = std::get<Rational>(std::get<V::Optional>(bv->data).element.at(0).data);
        result.set("left", state_->record(l));
        result.set("right", state_->record(r));
        result.set("reported_difference", number(Rational{x - y}));
        result.set("reported_ratio", y == 0 ? R{} : number(Rational{x / y}));
        result.set("scope",
                   R::string("Exact arithmetic on two reported standard point tokens; shared "
                                "assay ID does not establish independent experiments, full "
                                "construct/condition comparability or biological disagreement."));
    } else
        throw std::invalid_argument("unknown biomedical evidence query");
    auto cert = R::object();
    cert.set("schema", R::string("gadget.biomedical.certificate.v1"));
    cert.set("source_digest", R::string(digest()));
    cert.set("snapshot_sha256", R(state_->receipt.native().at("snapshot_sha256").materialize_value()));
    cert.set("native_sha256", R(state_->receipt.native().at("native_sha256").materialize_value()));
    auto execution = R::object();
    execution.set("gadget_sources_sha256", R::string(GADGET_SOURCE_SHA256));
    execution.set("rexgraph_headers_sha256", R::string(GADGET_REXGRAPH_SHA256));
    cert.set("compiled_sources", execution);
    cert.set("query", R(query));
    cert.set("result", result);
    cert.set("arithmetic",
             R::string("RexGraph native role traversal, three-valued exact predicates and GMP "
                          "rationals; no inference or identity coercion"));
    cert.set("interpretation_limits",
             R::string("Source assertions and explicit cross-references only. Mapping is not "
                          "identity equivalence, binding is not efficacy, omitted edges are not "
                          "negative evidence. API arrival is not biological time."));
    return retain_operation(state_->publication, input, std::move(cert).release());
}
Json EvidenceStore::execute(const Json &query) const {
    const auto operation = execute_record(admit_record_view(query));
    return seal(render_record_view(operation.at("result")), "certificate_sha256");
}
} // namespace gadget
