#include <algorithm>
#include <gadget/source.hpp>
#include <gadget/snapshot.hpp>
#include <limits>
#include <rfc/core/axioms.hpp>
#include <rfc/core/builder.hpp>
#include <stdexcept>

namespace gadget {
namespace {
Json row_json(const CsvTable &table, const CsvRow &row) {
    auto j = Json::object();
    for (std::size_t i = 0; i < table.header.size(); ++i)
        j.set(table.header[i], Json::string(row.cells[i]));
    j.set("source_record", Json::integer(row.record));
    return j;
}
std::uint32_t intern(Workspace &w, const std::string &gene) {
    if (gene.empty())
        throw std::invalid_argument("empty gene identity");
    auto [it, added] = w.gene_index.emplace(gene, static_cast<std::uint32_t>(w.genes.size()));
    if (added) {
        if (w.genes.size() >= std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("gene identity capacity");
        w.genes.push_back(gene);
    }
    return it->second;
}
std::pair<std::size_t, std::string> capture(Workspace &w, const fs::path &relative,
                                            const std::string &parser) {
    auto path = w.root / relative;
    auto canonical = fs::canonical(path);
    auto local = canonical.lexically_relative(w.root);
    if (local.empty() || *local.begin() == "..")
        throw std::invalid_argument("source escapes dataset root");
    auto bytes = read_text(path);
    auto index = w.files.size();
    w.files.push_back({relative.generic_string(), "sha256:" + sha256(bytes), parser, bytes.size()});
    w.raw_sources.emplace_back(index+1,std::span<const std::byte>{reinterpret_cast<const std::byte *>(bytes.data()),bytes.size()});
    return {index, std::move(bytes)};
}
} // namespace
Workspace load_pager(const fs::path &root) {
    Workspace w;
    w.root = fs::canonical(root);
    std::map<std::string, std::vector<std::pair<std::uint32_t, unsigned>>> annotations;
    std::map<std::string, std::map<std::string, Json>> metadata;
    const std::map<std::string, std::string> get_paths = {
        {"AD", "gene lists/Alzheimer_Disease_gene_analysis 1.csv"},
        {"HD", "gene lists/Huntington_Disease_gene_analysis_new 1.csv"}};
    const std::map<std::string, std::string> meta_paths = {
        {"AD", "pag lists/alzheimer_pag_list.csv"}, {"HD", "pag lists/huntington_pag_list.csv"}};
    for (const auto &[context, path] : get_paths) {
        auto [file, bytes] = capture(w, path, "gadget.legacy_get.v1");
        auto table = read_csv(bytes);
        auto records = Json::array();
        for (const auto &row : table.rows) {
            auto gene = trim(table.get(row, "Gene"));
            auto id = intern(w, gene);
            unsigned bits = 0;
            for (char c : trim(table.get(row, "GET"))) {
                if (c == 'G')
                    bits |= 1;
                else if (c == 'E')
                    bits |= 2;
                else if (c == 'T')
                    bits |= 4;
                else if (c != '+' && c != ' ' && c != ',')
                    throw std::invalid_argument("invalid GET channel");
            }
            annotations[context].push_back({id, bits});
            auto observation = row_json(table, row);
            observation.set("source_file", Json::integer(file));
            records.push(observation);
        }
        w.measurements.emplace(context, std::move(records));
    }
    for (const auto &[context, path] : meta_paths) {
        auto [file, bytes] = capture(w, path, "gadget.legacy_pager_metadata.v1");
        auto table = read_csv(bytes);
        for (const auto &row : table.rows) {
            auto id = trim(table.get(row, "PAG ID"));
            auto observation = row_json(table, row);
            observation.set("source_file", Json::integer(file));
            auto numbers = Json::object();
            for (const auto &column :
                 {"Size", "Overlap", "Similarity score", "CoCo score", "p-value", "FDR"}) {
                auto literal = trim(table.get(row, column));
                if (!literal.empty() && literal != "N/A")
                    numbers.set(column, exact_number(exact_decimal(literal)));
            }
            observation.set("exact_statistics", numbers);
            if (!metadata[context].emplace(id, observation).second)
                throw std::invalid_argument("duplicate PAGER metadata ID: " + id);
        }
    }
    const std::map<std::string, std::string> directories = {
        {"AD", "pag data/Alzheimers"}, {"HD", "pag data/Huntingtons/huntington_all_pager_files"}};
    // These four zero-byte exports are part of the supplied legacy source.
    const std::set<std::string> expected_empty = {
        "pag data/Alzheimers/MOX010352_genes.csv", "pag data/Alzheimers/TAX019115_genes.csv",
        "pag data/Alzheimers/TAX026391_genes.csv",
        "pag data/Huntingtons/huntington_all_pager_files/MOX007392_genes.csv"};
    std::map<std::string, std::size_t> assertion_by_content;
    std::map<std::string, std::set<std::string>> external_content;
    std::uint64_t occurrence = 1;
    std::size_t total_occurrences = 0;
    auto empty = Json::array();
    for (const auto &[context, directory] : directories) {
        std::vector<fs::path> paths;
        for (const auto &entry : fs::directory_iterator(w.root / directory))
            if (entry.is_regular_file() && entry.path().extension() == ".csv")
                paths.push_back(entry.path().lexically_relative(w.root));
        std::sort(paths.begin(), paths.end());
        if (paths.empty())
            throw std::invalid_argument("no PAGER assertions: " + directory);
        std::size_t duplicate_count = 0, raw_count = 0, empty_count = 0;
        for (const auto &path : paths) {
            auto [file, bytes] = capture(w, path, "gadget.legacy_pager_assertion.v1");
            if (bytes.empty()) {
                if (!expected_empty.contains(path.generic_string()))
                    throw std::invalid_argument("unexpected empty assertion: " + path.string());
                empty.push(Json::string(path.generic_string()));
                ++empty_count;
                continue;
            }
            auto table = read_csv(bytes);
            std::map<std::string, std::pair<Participant, std::vector<Origin>>> parsed;
            for (const auto &row : table.rows) {
                ++raw_count;
                auto gene = trim(table.get(row, "Gene symbol"));
                auto literal = trim(table.get(row, "RP_score"));
                auto value = exact_decimal(literal);
                if (value < 0)
                    throw std::invalid_argument("negative RP magnitude: " + path.string());
                auto [it, added] =
                    parsed.emplace(gene, std::make_pair(Participant{0, value, literal, 0, {}},
                                                        std::vector<Origin>{}));
                if (!added) {
                    ++duplicate_count;
                    if (it->second.first.rp != value)
                        throw std::invalid_argument("conflicting duplicate RP: " + path.string() +
                                                    ":" + std::to_string(row.record));
                }
                it->second.second.push_back({file, row.record});
            }
            if (parsed.empty())
                throw std::invalid_argument("empty retained participant set: " + path.string());
            auto external = path.stem().string();
            if (external.ends_with("_genes"))
                external.resize(external.size() - 6);
            // Exact retained content is the equality criterion. Digest collisions
            // cannot merge assertions: this map is keyed by the complete string.
            std::string canonical = std::to_string(external.size()) + ":" + external;
            for (const auto &[gene, observation] : parsed) {
                auto number = observation.first.rp.get_str();
                canonical += std::to_string(gene.size()) + ":" + gene +
                             std::to_string(number.size()) + ":" + number;
            }
            auto [entry, added] = assertion_by_content.emplace(canonical, w.assertions.size());
            auto rid = entry->second;
            if (added) {
                Assertion assertion;
                assertion.external = external;
                assertion.id = external + "@sha256:" + sha256(canonical);
                for (auto &[gene, observation] : parsed) {
                    auto p = std::move(observation.first);
                    p.gene = intern(w, gene);
                    p.occurrence = occurrence++;
                    p.origins = std::move(observation.second);
                    assertion.participants.push_back(std::move(p));
                }
                total_occurrences += assertion.participants.size();
                w.assertions.push_back(std::move(assertion));
            } else {
                auto &assertion = w.assertions.at(rid);
                std::size_t pos = 0;
                for (const auto &[gene, observation] : parsed) {
                    (void)gene;
                    auto &origins = assertion.participants.at(pos++).origins;
                    origins.insert(origins.end(), observation.second.begin(),
                                   observation.second.end());
                }
            }
            auto &assertion = w.assertions.at(rid);
            assertion.files.push_back(file);
            assertion.contexts.insert(context);
            external_content[external].insert(assertion.id);
            if (auto it = metadata[context].find(external); it != metadata[context].end())
                assertion.observations.emplace(context, it->second);
        }
        w.audit.set(context + "_files", Json::integer(paths.size()));
        w.audit.set(context + "_empty_files", Json::integer(empty_count));
        w.audit.set(context + "_raw_rows", Json::integer(raw_count));
        w.audit.set(context + "_duplicate_rows", Json::integer(duplicate_count));
    }
    rfc::core::RelationSignature signature;
    auto undirected = signature.aspects().add("undirected");
    auto participant = signature.aspects().add("participant");
    auto rp_observed = signature.aspects().add("rp_observed");
    auto type = signature.add_type("PAGER.symmetric_participant_set");
    rfc::core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < w.genes.size(); ++i)
        builder.declare_witness_identity(rfc::core::IdentityKey{i + 1});
    for (std::size_t i = 0; i < w.assertions.size(); ++i) {
        std::vector<rfc::core::BoundaryOccurrenceInput> incidences;
        for (const auto &p : w.assertions[i].participants)
            incidences.emplace_back(rfc::core::IdentityKey{p.gene + 1u}, true, false, Rational{1},
                                    std::vector<rfc::core::AspectValue>{{undirected, true},
                                                                        {participant, true},
                                                                        {rp_observed, true}},
                                    rfc::core::OccurrenceKey{p.occurrence});
        builder.add_grade_one(rfc::core::RelationKey{i + 1}, type, incidences);
    }
    w.complex = std::make_unique<rfc::core::RelationalComplex>(builder.freeze());
    for (const auto &context : {"AD", "HD"}) {
        auto mask = std::make_unique<rfc::algebra::SubcomplexMask>(*w.complex);
        for (std::size_t i = 0; i < w.assertions.size(); ++i)
            if (w.assertions[i].contexts.contains(context))
                mask->select({rfc::core::Grade{1}, static_cast<std::uint32_t>(i)});
        mask->close_boundaries(*w.complex);
        for (auto [gene, bits] : annotations[context]) {
            (void)bits;
            mask->select({rfc::core::Grade{0}, gene});
        }
        w.contexts.emplace(context, std::move(mask));
        for (std::size_t channel = 0; channel < 3; ++channel) {
            Vector values(w.genes.size(), Rational{0});
            std::vector<bool> present(w.genes.size(), false);
            for (auto [gene, bits] : annotations[context]) {
                present[gene] = true;
                if (bits & (1u << channel))
                    values[gene] = 1;
            }
            auto name = std::string(context) + "." + std::string(1, "GET"[channel]);
            auto semantics = Json::object();
            semantics.set("kind", Json::string("binary_annotation"));
            semantics.set("units", Json::string("binary"));
            semantics.set("method", Json::string("retained source GET evidence channel"));
            semantics.set("null_model", Json::string("not supplied"));
            semantics.set("missing_policy",
                          Json::string("declared zero extension to common ambient genes"));
            semantics.set("grade", Json::integer(0));
            semantics.set("context", Json::string(context));
            w.field_semantics.emplace(name, std::move(semantics));
            w.field_presence.emplace(name, std::move(present));
            w.fields.emplace(name, std::move(values));
        }
    }
    auto common = meet(*w.contexts.at("AD"), *w.contexts.at("HD"));
    std::size_t split = 0;
    for (const auto &[external, ids] : external_content) {
        (void)external;
        split += ids.size() - 1;
    }
    auto &audit=w.audit;
    audit.set("ambient_grade0", Json::integer(w.genes.size()));
    audit.set("ambient_grade1", Json::integer(w.assertions.size()));
    audit.set("occurrences", Json::integer(total_occurrences));
    for (const auto &context : {"AD", "HD"})
        for (std::uint16_t k = 0; k < 2; ++k)
            audit.set(std::string(context) + "_grade" + std::to_string(k),
                        Json::integer(w.contexts.at(context)->grade(rfc::core::Grade{k}).count()));
    audit.set("shared_grade0", Json::integer(common.grade(rfc::core::Grade{0}).count()));
    audit.set("strict_shared_grade1", Json::integer(common.grade(rfc::core::Grade{1}).count()));
    audit.set("same_external_id_different_assertion", Json::integer(split));
    auto report = rfc::core::bianchi_through(*w.complex, w.complex->max_grade());
    if (!report.satisfied())
        throw std::logic_error("canonical source Bianchi failure");
    audit.set("chain_issues", Json::integer(0));
    audit.set("empty_exports", empty);
    refresh_source_digest(w);
    admit_workspace(w);
    return w;
}
Json Workspace::manifest() const {
    auto j = Json::object();
    j.set("schema", Json::string("gadget.source.v1"));
    j.set("source_digest", Json::string(source_digest));
    j.set("root", Json::string(root.string()));
    j.set("audit", audit);
    auto records = Json::array();
    for (std::size_t i = 0; i < files.size(); ++i) {
        const auto &f = files[i];
        auto r = Json::object();
        r.set("index", Json::integer(i));
        r.set("path", Json::string(f.path));
        r.set("sha256", Json::string(f.hash));
        r.set("bytes", Json::integer(f.bytes));
        r.set("parser", Json::string(f.parser));
        records.push(r);
    }
    j.set("files", records);
    auto semantics = Json::object();
    for (const auto &[name, value] : field_semantics) {
        auto item = Json::parse(value.dump());
        auto count =
            std::count(field_presence.at(name).begin(), field_presence.at(name).end(), true);
        item.set("reported_coordinates", Json::integer(count));
        item.set("ambient_coordinates", Json::integer(genes.size()));
        semantics.set(name, item);
    }
    j.set("field_semantics", semantics);
    j.set("boundary_contract", Json::string("symmetric positive unit participant occurrences; RP "
                                            "and PAGER statistics retained separately"));
    return j;
}
Json Workspace::provenance(std::uint32_t gene) const {
    auto j = Json::object();
    j.set("gene", Json::string(genes.at(gene)));
    j.set("identity", Json::integer(gene + 1u));
    auto contributions = Json::array();
    for (std::size_t rid = 0; rid < assertions.size(); ++rid) {
        const auto &a = assertions[rid];
        for (const auto &p : a.participants)
            if (p.gene == gene) {
                auto c = Json::object();
                c.set("relation", Json::string(a.id));
                c.set("relation_key", Json::integer(rid + 1));
                c.set("external_id", Json::string(a.external));
                c.set("occurrence_key", Json::integer(p.occurrence));
                c.set("rp_literal", Json::string(p.literal));
                c.set("rp", exact_number(p.rp));
                auto origins = Json::array();
                for (auto o : p.origins) {
                    auto r = Json::object();
                    r.set("path", Json::string(files.at(o.file).path));
                    r.set("sha256", Json::string(files.at(o.file).hash));
                    r.set("record", Json::integer(o.record));
                    origins.push(r);
                }
                c.set("origins", origins);
                auto retained_contexts = Json::array();
                for (const auto &context : a.contexts)
                    retained_contexts.push(Json::string(context));
                c.set("contexts", retained_contexts);
                auto metadata = Json::object();
                for (const auto &[context, observation] : a.observations)
                    metadata.set(context, observation);
                c.set("source_statistics", metadata);
                contributions.push(c);
            }
    }
    j.set("contributions", contributions);
    auto observed = Json::array();
    for (const auto &[dataset, records] : measurements) {
        for (std::size_t i = 0; i < records.size(); ++i) {
            auto record = records.at(i);
            std::string identity;
            if (record.has("Gene"))
                identity = record.at("Gene").text();
            else if (record.has("gene"))
                identity = record.at("gene").text();
            if (identity == genes.at(gene)) {
                auto row = Json::parse(record.dump());
                row.set("dataset", Json::string(dataset));
                auto file = static_cast<std::size_t>(row.at("source_file").integer());
                row.set("path", Json::string(files.at(file).path));
                row.set("sha256", Json::string(files.at(file).hash));
                observed.push(row);
            }
        }
    }
    j.set("measurements", observed);
    return j;
}
} // namespace gadget
