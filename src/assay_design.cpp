#include "genomic_detail.hpp"
#include <algorithm>
#include <gadget/record_view.hpp>
#include <iterator>

namespace gadget::assay_detail {
namespace {
struct CompiledDesign {
    Json groups = Json::array(), certificate = Json::object();
};
CompiledDesign compile_design(const Axis &axis, const Json &design, const std::string &name) {
    members(design, {"mode"}, {"groups", "missing_tokens", "unassigned_policy", "overlap_policy"});
    CompiledDesign out;
    auto witnesses = Json::array();
    out.certificate.set("axis", Json::string(name));
    out.certificate.set("dictionary_sha256", Json::string("sha256:" + sha256(axis.bytes)));
    out.certificate.set("declaration", design);
    const bool identity = design.at("mode").text() == "identity";
    if (identity) {
        require(!design.has("groups") && !design.has("missing_tokens") &&
                    !design.has("unassigned_policy") && !design.has("overlap_policy"),
                "identity design has no grouping policies");
        for (std::size_t i = 0; i < axis.ids.size(); ++i) {
            auto group = Json::object(), member = Json::object(), entries = Json::array();
            group.set("id", Json::string(axis.ids[i]));
            member.set("id", Json::string(axis.ids[i]));
            member.set("coefficient", Json::string("1"));
            entries.push(member);
            group.set("members", entries);
            out.groups.push(group);
        }
    } else {
        require(design.at("mode").text() == "groups" && design.has("groups") &&
                    design.has("missing_tokens") && design.has("unassigned_policy") &&
                    design.has("overlap_policy"),
                "group design requires explicit missingness, exclusion and overlap policies");
        require(design.at("unassigned_policy").text() == "error" ||
                    design.at("unassigned_policy").text() == "exclude",
                "unknown design unassigned policy");
        require(design.at("overlap_policy").text() == "error" ||
                    design.at("overlap_policy").text() == "allow",
                "unknown design overlap policy");
        source_array(design.at("groups"), "design groups");
        require(design.at("groups").size() > 0, "design requires named groups");
    }
    std::set<std::string> missing;
    if (!identity) {
        source_array(design.at("missing_tokens"), "dictionary missing tokens");
        for (std::size_t i = 0; i < design.at("missing_tokens").size(); ++i) {
            const auto value = design.at("missing_tokens").at(i).text();
            require(missing.insert(value).second, "duplicate dictionary missing token");
        }
    }
    struct Predicate {
        std::string column;
        std::set<std::string> values;
    };
    std::vector<std::vector<Predicate>> predicates;
    std::set<std::string> targets;
    if (!identity)
        for (std::size_t g = 0; g < design.at("groups").size(); ++g) {
            const auto declared = design.at("groups").at(g);
            members(declared, {"id", "where", "coefficient"});
            require(targets.insert(text(declared.at("id"), source_byte_limit())).second,
                    "duplicate design group");
            (void)exact_rational(text(declared.at("coefficient"), source_byte_limit()));
            source_array(declared.at("where"), "design predicates");
            require(declared.at("where").size() > 0, "design group requires predicates");
            std::vector<Predicate> selection;
            std::set<std::string> columns;
            for (std::size_t p = 0; p < declared.at("where").size(); ++p) {
                const auto entry = declared.at("where").at(p);
                members(entry, {"column", "values"});
                Predicate predicate{text(entry.at("column"), source_byte_limit()), {}};
                require(columns.insert(predicate.column).second &&
                            axis.records.at(0).has(predicate.column),
                        "unknown or repeated dictionary predicate column");
                source_array(entry.at("values"), "predicate values");
                require(entry.at("values").size() > 0, "predicate requires values");
                for (std::size_t v = 0; v < entry.at("values").size(); ++v) {
                    const auto value = entry.at("values").at(v).text();
                    require(!missing.contains(value) && predicate.values.insert(value).second,
                            "predicate value is missing or repeated");
                }
                selection.push_back(std::move(predicate));
            }
            predicates.push_back(std::move(selection));
            auto group = Json::object();
            group.set("id", declared.at("id"));
            group.set("members", Json::array());
            out.groups.push(group);
        }
    // Build categorical postings once for the queried columns. Intersect
    // source-ordered row lists for each authored conjunction; do not rescan
    // every row for every specimen or feature group.
    std::map<std::string, std::map<std::string, std::vector<std::size_t>>> postings;
    for (const auto &group : predicates)
        for (const auto &p : group)
            for (const auto &value : p.values)
                postings[p.column].try_emplace(value);
    std::vector<std::vector<std::size_t>> row_groups(axis.ids.size());
    std::vector<std::vector<std::string>> row_missing(axis.ids.size());
    for (auto &[column, values] : postings)
        for (std::size_t i = 0; i < axis.ids.size(); ++i) {
            const auto value = axis.records.at(i).at(column).text();
            if (missing.contains(value))
                row_missing[i].push_back(column);
            else if (values.contains(value))
                values.at(value).push_back(i);
        }
    for (std::size_t g = 0; g < predicates.size(); ++g) {
        std::vector<std::size_t> selected;
        bool first = true;
        for (const auto &p : predicates[g]) {
            std::vector<std::size_t> candidates;
            for (const auto &value : p.values) {
                const auto &rows = postings.at(p.column).at(value);
                candidates.insert(candidates.end(), rows.begin(), rows.end());
            }
            std::sort(candidates.begin(), candidates.end());
            if (first) {
                selected = std::move(candidates);
                first = false;
            } else {
                std::vector<std::size_t> intersection;
                std::set_intersection(selected.begin(), selected.end(), candidates.begin(),
                                      candidates.end(), std::back_inserter(intersection));
                selected = std::move(intersection);
            }
        }
        for (const auto i : selected) {
            auto member = Json::object();
            member.set("id", Json::string(axis.ids[i]));
            member.set("coefficient", design.at("groups").at(g).at("coefficient"));
            out.groups.at(g).at("members").push(member);
            row_groups[i].push_back(g);
        }
    }
    for (std::size_t i = 0; i < axis.ids.size(); ++i) {
        auto witness = Json::object(), matched = Json::array(), absent = Json::array();
        if (identity)
            matched.push(Json::string(axis.ids[i]));
        else
            for (const auto g : row_groups[i])
                matched.push(out.groups.at(g).at("id"));
        if (!identity) {
            require(matched.size() > 0 || design.at("unassigned_policy").text() == "exclude",
                    "dictionary row has no declared group: " + axis.ids[i]);
            require(matched.size() <= 1 || design.at("overlap_policy").text() == "allow",
                    "dictionary row belongs to overlapping groups: " + axis.ids[i]);
        }
        for (const auto &column : row_missing[i])
            absent.push(Json::string(column));
        witness.set("id", Json::string(axis.ids[i]));
        witness.set("record", axis.records.at(i));
        witness.set("groups", matched);
        witness.set("excluded", Json::boolean(matched.size() == 0));
        witness.set("missing_columns", absent);
        witness.set("evidence", span(name, axis.bytes, axis.spans[i].first, axis.spans[i].second));
        witnesses.push(witness);
    }
    for (std::size_t g = 0; g < out.groups.size(); ++g)
        require(out.groups.at(g).at("members").size() > 0,
                "declared design group has no dictionary members");
    out.certificate.set("rows", witnesses);
    out.certificate.set("compiled_groups", out.groups);
    out.certificate.set("indexed_columns", cardinality(postings.size()));
    return out;
}
} // namespace

Json design_aggregation(const Matrix &matrix, const Json &declaration, const Json &request,
                        const rfc::runtime::CpuResourceCoordinator::Reservation &reservation,
                        const std::vector<NativeAssayPartition> &partitions,
                        std::optional<AssayFieldsSource> *retained) {
    members(request,
            {"schema", "operation", "assay", "feature_design", "observation_design", "map_method",
             "map_evidence", "support_policy"},
            {"count_policy", "workers"});
    const auto features = compile_design(matrix.features, request.at("feature_design"), "features");
    const auto observations =
        compile_design(matrix.observations, request.at("observation_design"), "observations");
    auto compiled = Json::object();
    for (const auto *key : {"schema", "assay", "map_method", "map_evidence", "support_policy",
                            "count_policy", "workers"})
        if (request.has(key))
            compiled.set(key, request.at(key));
    compiled.set("operation", Json::string("aggregate_assay"));
    compiled.set("feature_groups", features.groups);
    compiled.set("observation_groups", observations.groups);
    auto out = partitions.empty()
                   ? aggregate_analysis(matrix, declaration, compiled, reservation)
                   : native_aggregate_analysis(matrix.features, matrix.observations, partitions,
                                               declaration, compiled, reservation, retained);
    auto design = Json::object();
    design.set("features", features.certificate);
    design.set("observations", observations.certificate);
    design.set("compiled_request", compiled);
    out.set("design", design);
    if (retained) {
        require(retained->has_value(), "design aggregation fields missing");
        auto &value = std::get<NativeRecordValue::Object>(retained->value().declaration.data);
        for (auto &[name, summary] : value.members)
            if (name == "summary")
                summary = import_record_view(out);
    }
    return out;
}
} // namespace gadget::assay_detail
