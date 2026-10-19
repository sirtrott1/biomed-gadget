#include "pathway_ora.hpp"
#include <algorithm>
#include "scientific_reference.hpp"
#include <numeric>
#include <rfc/machine/checked_size.hpp>
#include <limits>
#include <map>

#include <set>

namespace gadget::pathway_detail {
namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::invalid_argument("pathway admission: " + message);
}
void members(const Json &object, const std::set<std::string> &names) {
    require(json_object_is_type(object.raw(), json_type_object), "expected object");
    for (const auto &name : names)
        require(object.has(name), "missing member: " + name);
    json_object_object_foreach(object.raw(), name, value) {
        (void)value;
        require(names.contains(name), "unknown member: " + std::string(name));
    }
}
std::size_t count(const Json &value, std::size_t maximum, bool positive = false) {
    const auto n = value.integer();
    require(n >= (positive ? 1 : 0) && static_cast<std::uint64_t>(n) <= maximum,
            "count outside declared budget");
    return static_cast<std::size_t>(n);
}
Json cardinality(std::size_t n) {
    require(n <= static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()),
            "JSON count overflow");
    return Json::integer(static_cast<std::int64_t>(n));
}
std::string text(const Json &value) {
    const auto out = value.text();
    require(!out.empty() && out.size() <= source_byte_limit() && out.find('\0') == std::string::npos,
            "empty/oversized identity or declaration");
    return out;
}
bool digest(const std::string &value) {
    if (value.size() != 71 || !value.starts_with("sha256:"))
        return false;
    return std::all_of(value.begin() + 7, value.end(),
                       [](char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); });
}
std::vector<std::string> strings(const Json &array, std::size_t maximum, bool unique) {
    require(json_object_is_type(array.raw(), json_type_array) && array.size() <= maximum,
            "string array budget");
    std::vector<std::string> out;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < array.size(); ++i) {
        auto value = text(array.at(i));
        require(!unique || seen.insert(value).second, "duplicate declared identity");
        out.push_back(std::move(value));
    }
    return out;
}
std::vector<std::size_t> indices(const Json &array, std::size_t maximum) {
    require(json_object_is_type(array.raw(), json_type_array) && array.size() <= maximum,
            "feature selector budget");
    std::vector<std::size_t> out;
    std::set<std::size_t> seen;
    for (std::size_t i = 0; i < array.size(); ++i) {
        const auto n = count(array.at(i), maximum);
        require(seen.insert(n).second, "duplicate feature selector");
        out.push_back(n);
    }
    return out;
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("pathway operation cancelled");
}
struct Row {
    std::size_t ordinal{}, start{}, end{};
    std::vector<std::string_view> cells;
};
std::vector<Row> rows(std::string_view bytes, std::size_t max_rows, std::size_t max_cells,
                      std::size_t field_bytes, std::stop_token stop) {
    require(bytes.find('\0') == std::string_view::npos, "NUL in source table");
    std::vector<Row> out;
    std::size_t start{}, cells{};
    while (start < bytes.size()) {
        cancelled(stop);
        require(out.size() < max_rows, "source row budget");
        auto end = bytes.find('\n', start);
        if (end == std::string_view::npos)
            end = bytes.size();
        const auto next = end + static_cast<std::size_t>(end < bytes.size());
        if (end > start && bytes[end - 1] == '\r')
            --end;
        require(end > start, "empty source record");
        Row row{out.size(), start, next, {}};
        for (std::size_t at = start;;) {
            require(cells++ < max_cells, "source membership/cell budget");
            auto tab = bytes.find('\t', at);
            if (tab == std::string_view::npos || tab > end)
                tab = end;
            require(tab - at <= field_bytes, "source field byte budget");
            const auto cell = bytes.substr(at, tab - at);
            require(cell.find('\r') == std::string_view::npos, "embedded CR in table");
            row.cells.push_back(cell);
            if (tab == end)
                break;
            at = tab + 1;
        }
        out.push_back(std::move(row));
        start = next;
    }
    require(!out.empty(), "empty source table");
    return out;
}
Json span(const Row &row) {
    auto out = Json::object();
    out.set("start", cardinality(row.start));
    out.set("end", cardinality(row.end));
    out.set("ordinal", cardinality(row.ordinal));
    return out;
}
Json names(const std::set<std::string> &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(Json::string(value));
    return out;
}
} // namespace
Parameters admit_request(const Json &request) {
    members(request,
            {"schema", "operation", "analysis_sha256", "gene_sets_sha256", "feature_namespace",
             "features_path", "feature_id_path", "background_indices", "selected_indices",
             "mapping_policy", "set_ids", "size_filter", "membership_policy", "correction",
             "selection_method", "null_assumptions", "limits"});
    require(request.at("schema").text() == "gadget.pathways.overrepresentation.v1" &&
                request.at("operation").text() == "pathway_overrepresentation",
            "request schema/operation");
    Parameters out;
    out.analysis_sha256 = text(request.at("analysis_sha256"));
    out.gene_sets_sha256 = text(request.at("gene_sets_sha256"));
    require(digest(out.analysis_sha256) && digest(out.gene_sets_sha256),
            "exact archive SHA pins required");
    out.feature_namespace = text(request.at("feature_namespace"));
    out.features_path = request.at("features_path");
    out.feature_id_path = request.at("feature_id_path");
    for (const auto &path : {out.features_path, out.feature_id_path}) {
        require(json_object_is_type(path.raw(), json_type_array), "native feature path requires an array");
        for (std::size_t i=0;i<path.size();++i) {
            const auto item=path.at(i);
            require(json_object_is_type(item.raw(),json_type_string) ||
                    (json_object_is_type(item.raw(),json_type_int) && item.integer()>=0), "native feature path component");
        }
    }
    require(out.features_path.size()>0, "native feature path must be nonempty");
    const auto limits = request.at("limits");
    members(limits, {"source_bytes", "mapping_rows", "membership_cells", "sets", "features",
                     "field_bytes", "kernel_work", "archive_bytes", "declared_buffer_bytes"});
    const auto largest = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
    out.source_bytes = count(limits.at("source_bytes"), source_byte_limit(), true);
    out.mapping_rows = count(limits.at("mapping_rows"), largest, true);
    out.membership_cells = count(limits.at("membership_cells"), largest, true);
    out.sets = count(limits.at("sets"), largest, true);
    out.features = count(limits.at("features"), largest, true);
    out.field_bytes = count(limits.at("field_bytes"), out.source_bytes, true);
    out.kernel_work = count(limits.at("kernel_work"), largest, true);
    out.archive_bytes = count(limits.at("archive_bytes"), source_byte_limit(), true);
    out.declared_buffer_bytes = count(limits.at("declared_buffer_bytes"), largest, true);
    out.background = indices(request.at("background_indices"), out.features);
    out.selected = indices(request.at("selected_indices"), out.features);
    require(!out.background.empty(), "explicit background feature selectors required");
    const std::set<std::size_t> background(out.background.begin(), out.background.end());
    for (const auto index : out.selected)
        require(background.contains(index), "selected feature outside declared background");
    out.set_ids = strings(request.at("set_ids"), out.sets, true);
    require(!out.set_ids.empty(), "explicit tested set family required");
    const auto filter = request.at("size_filter");
    members(filter, {"minimum", "maximum", "basis"});
    require(filter.at("basis").text() == "members_in_background",
            "set size basis must be explicit background intersection");
    out.min_size = count(filter.at("minimum"), largest);
    out.max_size = count(filter.at("maximum"), largest);
    require(out.min_size <= out.max_size, "reversed set size filter");
    const auto policy = request.at("mapping_policy");
    members(policy, {"ambiguous", "unmapped", "many_to_one", "duplicate_occurrences"});
    out.ambiguous = text(policy.at("ambiguous"));
    out.unmapped = text(policy.at("unmapped"));
    out.collisions = text(policy.at("many_to_one"));
    require((out.ambiguous == "exclude" || out.ambiguous == "refuse") &&
                (out.unmapped == "exclude" || out.unmapped == "refuse"),
            "explicit ambiguous/unmapped mapping policies");
    require((out.collisions == "refuse" || out.collisions == "count_member_once") &&
                policy.at("duplicate_occurrences").text() == "retain_count_once",
            "explicit mapping collision/occurrence policy");
    const auto membership = request.at("membership_policy");
    members(membership, {"duplicates", "empty", "outside_background"});
    require(membership.at("duplicates").text() == "retain_count_once" &&
                membership.at("empty").text() == "retain_exclude" &&
                membership.at("outside_background").text() == "retain_exclude",
            "explicit gene set occurrence policies");
    out.correction = text(request.at("correction"));
    require(out.correction == "none" || out.correction == "benjamini_hochberg" || out.correction == "benjamini_yekutieli" ||
                out.correction == "bonferroni",
            "unsupported multiplicity correction");
    (void)text(request.at("selection_method"));
    require(!strings(request.at("null_assumptions"), 64, false).empty(),
            "declared null assumptions required");
    return out;
}
scientific_detail::ResultValue overrepresentation(const std::vector<std::string> &feature_ids, const Source &source,
                        const Json &request, scientific_detail::ResultFields &fields, std::stop_token stop) {
    cancelled(stop);
    const auto p = admit_request(request);
    require(feature_ids.size() <= p.features, "actual feature occurrence budget");
    for (const auto &id : feature_ids)
        require(!id.empty() && id.size() <= p.field_bytes && id.find('\0') == std::string::npos,
                "feature identity representation");
    require(!source.id.empty() && !source.release.empty() && !source.identifier_namespace.empty() &&
                !source.mapping_id.empty() && !source.mapping_release.empty(),
            "released gene set and mapping source declarations required");
    require(source.mapping_origin == "released" || source.mapping_origin == "explicit_declared",
            "mapping origin must be explicitly declared");
    require(source.feature_namespace == p.feature_namespace, "mapping feature namespace mismatch");
    require(source.gmt.size() <= p.source_bytes &&
                source.mapping.size() <= p.source_bytes - source.gmt.size(),
            "aggregate original source byte budget");
    require(source.sha256 == "sha256:" + sha256(source.gmt) &&
                source.mapping_sha256 == "sha256:" + sha256(source.mapping),
            "original source SHA mismatch");
    for (const auto i : p.background)
        require(i < feature_ids.size(), "background selector outside actual native features");
    std::size_t retained_bytes = rfc::machine::checked_size_add(
        rfc::machine::checked_size_mul(source.gmt.size() + source.mapping.size(), 32),
        rfc::machine::checked_size_mul(feature_ids.size(), 512));
    for(const auto &id:feature_ids)
        retained_bytes=rfc::machine::checked_size_add(retained_bytes,
            rfc::machine::checked_size_mul(id.size(),8));
    require(retained_bytes <= p.declared_buffer_bytes, "mapping and feature working byte admission");

    const auto mapping =
        rows(source.mapping, p.mapping_rows, p.membership_cells, p.field_bytes, stop);
    require(mapping[0].cells == std::vector<std::string_view>{"feature_id", "member", "evidence"},
            "explicit correspondence TSV header");
    std::map<std::string, std::set<std::string>> choices;
    auto mapping_occurrences = Json::array();
    std::map<std::string, std::vector<std::size_t>> occurrences;
    for (std::size_t i = 0; i < feature_ids.size(); ++i)
        occurrences[feature_ids[i]].push_back(i);
    for (std::size_t i = 1; i < mapping.size(); ++i) {
        cancelled(stop);
        const auto &row = mapping[i];
        require(row.cells.size() == 3 && !row.cells[0].empty() && !row.cells[2].empty(),
                "mapping occurrence width/identity/evidence");
        const std::string id(row.cells[0]), member(row.cells[1]);
        if (!member.empty())
            choices[id].insert(member);
        auto item = span(row);
        item.set("feature_id", Json::string(id));
        item.set("member", Json::string(member));
        item.set("evidence", Json::string(std::string(row.cells[2])));
        auto matched = Json::array();
        for (const auto index : occurrences[id])
            matched.push(cardinality(index));
        item.set("actual_feature_indices", matched);
        mapping_occurrences.push(item);
    }
    const std::set<std::size_t> selected(p.selected.begin(), p.selected.end());
    std::set<std::string> background, foreground;
    auto admission = Json::array();
    std::map<std::string, std::vector<std::size_t>> member_features;
    for (const auto index : p.background) {
        cancelled(stop);
        const auto &id = feature_ids[index];
        const auto &mapped = choices[id];
        const auto reason = mapped.empty()      ? "unmapped"
                            : mapped.size() > 1 ? "ambiguous"
                                                : "mapped";
        require(mapped.size() != 0 || p.unmapped == "exclude",
                "unmapped background feature refused");
        require(mapped.size() <= 1 || p.ambiguous == "exclude",
                "ambiguous background feature refused");
        auto item = Json::object();
        item.set("feature_index", cardinality(index));
        item.set("feature_id", Json::string(id));
        item.set("selected", Json::boolean(selected.contains(index)));
        item.set("status", Json::string(reason));
        item.set("candidate_members", names(mapped));
        if (mapped.size() == 1) {
            const auto &member = *mapped.begin();
            require(p.collisions == "count_member_once" || !background.contains(member),
                    "many-to-one mapping refused");
            background.insert(member);
            if (selected.contains(index))
                foreground.insert(member);
            member_features[member].push_back(index);
        }
        admission.push(item);
    }
    require(!background.empty(), "mapped background must be nonempty");
    const auto sets = rows(source.gmt, p.sets, p.membership_cells, p.field_bytes, stop);
    std::map<std::string, const Row *> by_id;
    for (const auto &row : sets) {
        require(row.cells.size() >= 2 && !row.cells[0].empty(),
                "GMT requires a name and description");
        require(by_id.emplace(std::string(row.cells[0]), &row).second,
                "duplicate GMT set identity");
    }
    using Value = scientific_detail::ResultValue;
    using Rational = rfc::numeric::Rational;
    using Integer = rfc::numeric::Integer;
    auto tests = Value::array();
    auto excluded = Json::array();
    std::vector<Rational> probabilities;
    std::vector<std::size_t> tested_indices;
    std::size_t work{};
    for (const auto &id : p.set_ids) {
        cancelled(stop);
        require(by_id.contains(id), "declared set absent from released source");
        const auto &row = *by_id.at(id);
        std::set<std::string> members_in_background, all_members, hits;
        auto membership = Json::array();
        for (std::size_t c = 2; c < row.cells.size(); ++c) {
            const std::string member(row.cells[c]);
            const bool repeated = !member.empty() && !all_members.insert(member).second;
            auto occurrence = Json::object();
            occurrence.set("column", cardinality(c));
            occurrence.set("member", Json::string(member));
            occurrence.set("status",
                           Json::string(member.empty()                 ? "empty"
                                        : repeated                     ? "repeated"
                                        : !background.contains(member) ? "outside_background"
                                                                       : "in_background"));
            membership.push(occurrence);
            if (!member.empty() && background.contains(member)) {
                members_in_background.insert(member);
                if (foreground.contains(member))
                    hits.insert(member);
            }
        }
        auto test = Json::object();
        test.set("id", Json::string(id));
        test.set("description", Json::string(std::string(row.cells[1])));
        test.set("source_span", span(row));
        test.set("membership_occurrences", membership);
        test.set("source_unique_members", cardinality(all_members.size()));
        test.set("background_members", names(members_in_background));
        test.set("selected_members", names(hits));
        test.set("N", cardinality(background.size()));
        test.set("K", cardinality(members_in_background.size()));
        test.set("draws", cardinality(foreground.size()));
        test.set("hits", cardinality(hits.size()));
        if (members_in_background.size() < p.min_size ||
            members_in_background.size() > p.max_size) {
            test.set("status", Json::string("outside_size_filter"));
            excluded.push(test);
            continue;
        }
        const auto N = background.size(), K = members_in_background.size();
        const auto draws = foreground.size(), observed = hits.size();
        const auto lower = draws > N - K ? draws - (N - K) : 0;
        const auto upper = std::min(K, draws);
        require(lower <= observed && observed <= upper, "observed intersection outside finite law");
        const auto states = upper - lower + 1;
        using rfc::machine::checked_size_add;
        using rfc::machine::checked_size_mul;
        // A binomial coefficient has at most N + 1 bits. Charge retained law,
        // finite calibration maps, exact result fields and their construction.
        const auto limb_bytes = checked_size_add(N / 8, 1);
        const auto demand_bytes = checked_size_mul(states, checked_size_add(4096,
            checked_size_mul(limb_bytes, 256)));
        require(demand_bytes <= p.declared_buffer_bytes - retained_bytes, "exact finite law working byte admission");
        retained_bytes += demand_bytes;
        const auto demand = checked_size_add(checked_size_mul(N, 3), checked_size_mul(states, 32));
        require(demand <= p.kernel_work - work, "declared exact finite law work budget");
        work += demand;
        const auto choose = [&](std::size_t n, std::size_t k) {
            cancelled(stop);
            require(k <= n, "binomial domain");
            k = std::min(k, n - k);
            require(k <= std::numeric_limits<unsigned long>::max(), "binomial index addressability");
            Integer value;
            const Integer owner(std::to_string(n));
            mpz_bin_ui(value.get_mpz_t(), owner.get_mpz_t(), static_cast<unsigned long>(k));
            return value;
        };
        const auto denominator = choose(N, draws);
        Rational mass(choose(K, lower) * choose(N - K, draws - lower), denominator);
        mass.canonicalize();
        std::vector<Vector> scores;
        std::vector<std::vector<bool>> defined;
        Vector weights;
        scores.reserve(states); defined.reserve(states); weights.reserve(states);
        for (std::size_t k = lower;; ++k) {
            cancelled(stop);
            scores.push_back({Rational(std::to_string(k))});
            defined.push_back({true}); weights.push_back(mass);
            if (k == upper) break;
            // Adjacent intersection counts compress the uniform subset law
            // without enumerating subsets or identifying independent carriers.
            mass *= Rational(Integer(std::to_string(K - k)) * Integer(std::to_string(draws - k)),
                Integer(std::to_string(k + 1)) * Integer(std::to_string((N - K) - (draws - k - 1))));
            mass.canonicalize();
        }
        auto family = Json::array(), hypothesis = Json::object();
        hypothesis.set("id", Json::string(id)); family.push(hypothesis);
        scientific_detail::FiniteScoreValues values;
        const auto state_axis = fields.axis(states), hypothesis_axis = fields.axis(1);
        auto native_test = Value(test);
        native_test.set("status", Value::string("tested"));
        native_test.set("intersection_counts", fields.vector([&] {
            Vector counts;
            for (std::size_t k = lower;; ++k) {
                counts.emplace_back(std::to_string(k)); if (k == upper) break;
            }
            return counts;
        }(), state_axis));
        native_test.set("reference_masses", fields.vector(weights, state_axis));
        native_test.set("observed_state", Value::count(observed - lower));
        native_test.set("calibration", scientific_detail::calibrate_finite_scores_fields(
            scores, defined, weights, observed - lower, family, {}, fields,
            state_axis, hypothesis_axis, &values));
        require(values.reference_tail.size() == 1 && values.reference_tail[0].has_value(),
                "finite law tail unavailable");
        native_test.set("p_value", fields.scalar(*values.reference_tail[0]));
        tested_indices.push_back(tests.size());
        probabilities.push_back(*values.reference_tail[0]);
        tests.push(std::move(native_test));
    }
    // These are declared testing policies over exact marginal tails. They do
    // not manufacture a joint reference law from independent marginal laws.
    const auto m = probabilities.size();
    std::vector<Rational> adjusted = probabilities;
    Rational multiplier(std::to_string(m));
    if (p.correction == "benjamini_yekutieli") {
        Rational harmonic{0};
        for (std::size_t i = 1; i <= m; ++i) {
            cancelled(stop); harmonic += Rational(Integer{1}, Integer(std::to_string(i)));
        }
        const auto bits=rfc::machine::checked_size_add(mpz_sizeinbase(harmonic.get_num_mpz_t(),2),
            mpz_sizeinbase(harmonic.get_den_mpz_t(),2));
        const auto additional=rfc::machine::checked_size_mul(m,rfc::machine::checked_size_mul(bits/8+2,16));
        require(additional<=p.declared_buffer_bytes-retained_bytes,"exact multiplicity working byte admission");
        retained_bytes+=additional;
        multiplier *= harmonic;
    }
    if (p.correction == "bonferroni") {
        for (auto &value : adjusted) value = std::min(Rational{1}, Rational{value * multiplier});
    } else if (p.correction != "none") {
        std::vector<std::size_t> order(m); std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) {return probabilities[a] < probabilities[b];});
        Rational following{1};
        for (std::size_t j = m; j > 0; --j) {
            cancelled(stop);
            following = std::min(following, Rational{probabilities[order[j - 1]] * multiplier / Rational(std::to_string(j))});
            adjusted[order[j - 1]] = following;
        }
    }
    auto complete_tests = Value::array();
    for (std::size_t i = 0; i < m; ++i) {
        auto test = tests.at(tested_indices[i]);
        test.set("adjusted_p_value", fields.scalar(adjusted[i]));
        complete_tests.push(std::move(test));
    }
    auto collisions = Json::array();
    for (const auto &[member, indices] : member_features)
        if (indices.size() > 1) {
            auto item = Json::object(), ids = Json::array();
            item.set("member", Json::string(member));
            for (const auto index : indices)
                ids.push(cardinality(index));
            item.set("feature_indices", ids);
            collisions.push(item);
        }
    auto out = Json::object();
    out.set("schema", Json::string("gadget.pathways.overrepresentation.result.v1"));
    out.set("method", Json::string("exact_finite_intersection_law"));
    out.set("analysis_sha256", Json::string(p.analysis_sha256));
    out.set("gene_sets_sha256", Json::string(p.gene_sets_sha256));
    out.set("source_id", Json::string(source.id));
    out.set("release", Json::string(source.release));
    out.set("identifier_namespace", Json::string(source.identifier_namespace));
    out.set("gmt_sha256", Json::string(source.sha256));
    out.set("mapping_id", Json::string(source.mapping_id));
    out.set("mapping_release", Json::string(source.mapping_release));
    out.set("mapping_origin", Json::string(source.mapping_origin));
    out.set("mapping_sha256", Json::string(source.mapping_sha256));
    out.set("feature_namespace", Json::string(p.feature_namespace));
    out.set("feature_admission", admission);
    out.set("mapping_occurrences", mapping_occurrences);
    out.set("many_to_one_occurrences", collisions);
    out.set("background_members", names(background));
    out.set("selected_members", names(foreground));
    out.set("declared_set_family", request.at("set_ids"));
    out.set("tested_family_size", cardinality(tests.size()));
    out.set("excluded_sets", excluded);

    out.set("correction", Json::string(p.correction));
    out.set("kernel_work_bound", cardinality(work));
    out.set("numerical_kind", Json::string("exact rational native fields"));
    out.set("multiplicity_scope", Json::string("Marginal uniform subset laws are explicit model assumptions. Bonferroni uses the union bound; BH requires valid null tails and independence or suitable positive dependence; BY permits arbitrary dependence with valid null tails. No joint pathway law is inferred."));
    out.set("null_hypothesis", Json::string("Uniform subset of the mapped declared background with "
                                            "the observed selected-member count"));
    out.set("null_assumptions", request.at("null_assumptions"));
    out.set("selection_method", request.at("selection_method"));
    out.set("interpretation",
            Json::string("Competitive membership overrepresentation conditional on the explicit "
                         "mapping, background, selection and tested family; not pathway activity, "
                         "effect direction or causation"));
    cancelled(stop);
    auto native = Value(out);
    native.set("tests", std::move(complete_tests));
    return native;
}
} // namespace gadget::pathway_detail
