#include "genomic_detail.hpp"
#include "allocation_detail.hpp"
#include "source_grade_fields.hpp"
#include <algorithm>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail {
namespace {
using rfc::numeric::ExactMatrix;
std::string label(const Json &value) { return text(value, source_byte_limit()); }
Json values(const Vector &v) {
    auto out = Json::array();
    for (const auto &q : v)
        out.push(Json::string(q.get_str()));
    return out;
}
Json bound(const Json &id, const Rational &lo, const Rational &hi, bool contrast = false) {
    auto b = Json::object();
    b.set("id", id);
    b.set("lower", Json::string(lo.get_str()));
    b.set("upper", Json::string(hi.get_str()));
    b.set("identifiable", Json::boolean(lo == hi));
    if (contrast)
        b.set("direction", Json::string(lo > 0               ? "positive_for_every_allocation"
                                        : hi < 0             ? "negative_for_every_allocation"
                                        : lo == 0 && hi == 0 ? "zero_for_every_allocation"
                                                             : "unresolved"));
    return b;
}
Json sequence_evidence(const ResourceImages &images, const Json &source, const Json &sample,
                       Vector *native_counts = nullptr) {
    members(source, {"mode", "reference_resource", "strands", "ambiguity", "maximum_mismatches",
                     "evidence"});
    auto query = Json::object();
    query.set("schema", Json::string("gadget.assays.query.v1"));
    query.set("operation", Json::string("reference_sequence_matches"));
    query.set("query_resource", sample.at("resource"));
    for (const auto *name :
         {"reference_resource", "strands", "ambiguity", "maximum_mismatches", "evidence"})
        query.set(name, source.at(name));
    query.set("minimum_overhang",
              Json::integer(1)); // No genomic junctions in transcript coordinates.
    query.set("junctions", Json::array());
    const auto search = reference_matches(images, query);
    const auto reads = sequence_records(images.at(label(sample.at("resource"))));
    const auto references = sequence_records(images.at(label(source.at("reference_resource"))));
    auto candidates = Json::array(), classes = Json::array(), counts = Json::array(),
         observations = Json::array();
    std::map<std::string, std::size_t> ids;
    for (const auto &r : references) {
        require(r.id != "__unassigned__" && ids.emplace(r.id, ids.size()).second,
                "unique transcript identity; __unassigned__ is reserved");
        candidates.push(Json::string(r.id));
    }
    const auto unassigned = ids.size();
    candidates.push(Json::string("__unassigned__"));
    std::vector<std::set<std::size_t>> memberships(reads.size());
    std::vector<Json> witnesses;
    for (std::size_t i = 0; i < reads.size(); ++i)
        witnesses.push_back(Json::array());
    for (std::size_t i = 0; i < search.at("matches").size(); ++i) {
        const auto m = search.at("matches").at(i);
        const auto row =
            count(m.at("query_evidence").at("record_index"), reads.size() - 1, "read source index");
        memberships[row].insert(ids.at(label(m.at("reference_id"))));
        witnesses[row].push(cardinality(i));
    }
    std::map<std::vector<std::size_t>, std::size_t> groups;
    std::vector<std::size_t> class_counts;
    for (std::size_t i = 0; i < reads.size(); ++i) {
        if (memberships[i].empty())
            memberships[i].insert(unassigned);
        const std::vector<std::size_t> p(memberships[i].begin(), memberships[i].end());
        if (!groups.contains(p)) {
            const auto index = groups.size();
            groups.emplace(p, index);
            class_counts.push_back(0);
            auto c = Json::object(), members = Json::array();
            c.set("id", Json::string("class:" + std::to_string(index)));
            for (auto candidate : p)
                members.push(candidates.at(candidate));
            c.set("candidates", members);
            classes.push(c);
        }
        const auto index = groups.at(p);
        class_counts[index] = rfc::machine::checked_size_add(class_counts[index], 1);
        auto o = Json::object();
        o.set("query", Json::string(reads[i].id));
        o.set("source_evidence", reads[i].evidence);
        o.set("class_index", cardinality(index));
        o.set("match_indices", witnesses[i]);
        observations.push(o);
    }
    for (auto n : class_counts) {
        counts.push(Json::string(std::to_string(n)));
        if (native_counts)
            native_counts->emplace_back(std::to_string(n));
    }
    auto out = Json::object();
    out.set("candidates", candidates);
    out.set("classes", classes);
    out.set("counts", counts);
    out.set("observations", observations);
    out.set("search_request", query);
    out.set("search_result", search);
    out.set(
        "reference_sha256",
        Json::string("sha256:" + sha256(images.at(label(source.at("reference_resource"))).bytes)));
    return out;
}
struct Measure {
    Json id;
    std::size_t output, target;
    bool contrast;
    Vector coefficients;
};


NativeRecordValue &member(NativeRecordValue &value, std::string_view name) {
    for (auto &entry : std::get<NativeRecordValue::Object>(value.data).members)
        if (entry.first == name) return entry.second;
    throw std::invalid_argument("native expression declaration member absent");
}
NativeRecordValue coupled_expression(const Json &allocation, const Json &declaration,
    const std::vector<Vector> &counts, AssayFieldsSource &native,
    const std::shared_ptr<CpuResources> &resources, const NativeRecord *previous) {
    using N = NativeRecordValue;
    using allocation_detail::Problem;
    members(declaration, {"constraints", "witness_representation", "evidence"});
    (void)label(declaration.at("evidence"));
    const auto samples = allocation.at("samples"), classes = allocation.at("classes"),
               targets = allocation.at("targets"), contrasts = allocation.at("contrasts");
    std::map<std::string, std::size_t> sample_ids, target_ids;
    std::vector<std::set<std::string>> target_members;
    Vector scales;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        sample_ids.emplace(label(samples.at(i).at("id")), i);
        scales.push_back(exact_rational(label(samples.at(i).at("scale"))));
    }
    for (std::size_t i = 0; i < targets.size(); ++i) {
        target_ids.emplace(label(targets.at(i).at("id")), i);
        std::set<std::string> members;
        for (std::size_t j = 0; j < targets.at(i).at("candidates").size(); ++j)
            members.insert(label(targets.at(i).at("candidates").at(j)));
        target_members.push_back(std::move(members));
    }
    Problem problem;
    problem.domain = "bounded_integer";
    problem.names = Json::array();
    auto variables = Json::array(), origins = Json::array(), rows = Json::array();
    std::vector<std::vector<std::size_t>> conservation;
    Vector totals;
    std::vector<std::pair<std::size_t, std::string>> coordinates;
    std::vector<std::size_t> class_indices;
    for (std::size_t i = 0; i < samples.size(); ++i)
        for (std::size_t c = 0; c < classes.size(); ++c) {
            std::vector<std::size_t> columns;
            const auto candidates = classes.at(c).at("candidates");
            for (std::size_t j = 0; j < candidates.size(); ++j) {
                const auto id = "allocation:" + std::to_string(i) + ":" +
                                std::to_string(c) + ":" + std::to_string(j);
                auto v = Json::object(), origin = Json::object();
                v.set("id", Json::string(id)); variables.push(v); problem.names.push(v.at("id"));
                columns.push_back(problem.lower.size());
                problem.lower.emplace_back(0); problem.upper.push_back(counts.at(i).at(c));
                const Rational width = counts.at(i).at(c) + 1;
                problem.declared_box *= width.get_num();
                coordinates.emplace_back(i, label(candidates.at(j)));
                class_indices.push_back(c);
                origin.set("id", v.at("id")); origin.set("sample", samples.at(i).at("id"));
                origin.set("class", classes.at(c).at("id")); origin.set("candidate", candidates.at(j));
                origins.push(origin);
            }
            conservation.push_back(std::move(columns)); totals.push_back(counts.at(i).at(c));
        }
    const auto n = problem.lower.size();
    ExactMatrix target_map(rfc::machine::checked_size_mul(samples.size(), targets.size()), n);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t t = 0; t < targets.size(); ++t)
            if (target_members[t].contains(coordinates[j].second))
                target_map.set(coordinates[j].first * targets.size() + t, j, 1);
    const auto transpose = target_map.transpose();
    const auto target_row = [&](std::size_t i) {
        Vector selector(target_map.rows(), 0); selector.at(i) = 1;
        return transpose.apply(selector);
    };
    const auto declared = declaration.at("constraints");
    problem.constraints = ExactMatrix(rfc::machine::checked_size_add(conservation.size(), declared.size()), n);
    for (std::size_t i = 0; i < conservation.size(); ++i) {
        for (auto j : conservation[i]) problem.constraints.set(i, j, 1);
        problem.senses.emplace_back("eq"); problem.rhs.push_back(totals[i]);
        auto row = Json::object(); row.set("id", Json::string("conservation:" + std::to_string(i))); rows.push(row);
    }
    std::set<std::string> constraint_ids;
    for (std::size_t i = 0; i < declared.size(); ++i) {
        const auto row = declared.at(i);
        members(row, {"id", "terms", "sense", "rhs", "evidence"});
        require(constraint_ids.insert(label(row.at("id"))).second, "duplicate expression constraint identity");
        (void)label(row.at("evidence"));
        Vector coefficients(n, 0);
        for (std::size_t j = 0; j < row.at("terms").size(); ++j) {
            const auto term = row.at("terms").at(j);
            members(term, {"sample", "target", "coefficient"});
            require(sample_ids.contains(label(term.at("sample"))) && target_ids.contains(label(term.at("target"))),
                    "coupling term requires declared sample and target");
            coefficients = rfc::numeric::add(coefficients,
                rfc::numeric::scale(exact_rational(label(term.at("coefficient"))),
                    target_row(sample_ids.at(label(term.at("sample"))) * targets.size() +
                               target_ids.at(label(term.at("target"))))));
        }
        for (std::size_t j = 0; j < n; ++j) problem.constraints.set(conservation.size() + i, j, coefficients[j]);
        problem.senses.push_back(label(row.at("sense"))); problem.rhs.push_back(exact_rational(label(row.at("rhs"))));
        auto translated = Json::object(); translated.set("id", Json::string("declared:" + label(row.at("id")))); rows.push(translated);
    }
    problem.declaration = Json::object();
    for (const auto &[name, value] : std::vector<std::pair<std::string, std::string>>{
            {"schema", "gadget.assays.query.v1"}, {"operation", "constrained_allocations"}, {"domain", problem.domain}})
        problem.declaration.set(name, Json::string(value));
    problem.declaration.set("variables", variables); problem.declaration.set("constraints", rows);
    problem.declaration.set("units", allocation.at("units")); problem.declaration.set("evidence", declaration.at("evidence"));
    problem.declaration.set("witness_representation", declaration.at("witness_representation"));
    problem.declaration.set("typed_source", Json::boolean(true));
    std::vector<Measure> measures;
    for (std::size_t i = 0; i < samples.size(); ++i)
        for (std::size_t t = 0; t < targets.size(); ++t)
            measures.push_back({targets.at(t).at("id"), i, t, false,
                rfc::numeric::scale(scales[i], target_row(i * targets.size() + t))});
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto contrast = contrasts.at(c);
        const auto left = sample_ids.at(label(contrast.at("left"))), right = sample_ids.at(label(contrast.at("right")));
        for (std::size_t t = 0; t < targets.size(); ++t)
            measures.push_back({targets.at(t).at("id"), c, t, true,
                rfc::numeric::sub(rfc::numeric::scale(scales[left], target_row(left * targets.size() + t)),
                                  rfc::numeric::scale(scales[right], target_row(right * targets.size() + t)))});
    }
    SourceGradeFields store(native.endpoint);
    const auto assignment_axis = store.explicit_carrier(n), output_axis = store.explicit_carrier(target_map.rows());
    const auto target_map_key = store.map("coupled.target_map", assignment_axis, output_axis, target_map);
    std::map<std::string, std::size_t> candidate_indices;
    const auto candidate_ids = allocation.at("candidates");
    for (std::size_t i = 0; i < candidate_ids.size(); ++i) candidate_indices.emplace(label(candidate_ids.at(i)), i);
    ExactMatrix candidate_map(native.endpoint.source.grade({0}).cell_count(), n),
                class_map(native.endpoint.source.grade({1}).cell_count(), n);
    Vector sample_coordinates, variable_coordinates;
    for (std::size_t j = 0; j < n; ++j) {
        candidate_map.set(candidate_indices.at(coordinates[j].second), j, 1);
        class_map.set(class_indices[j], j, 1);
        sample_coordinates.emplace_back(std::to_string(coordinates[j].first));
        variable_coordinates.emplace_back(std::to_string(j));
    }
    const auto candidate_map_key = store.map("coupled.assignment_candidates", assignment_axis, 1, candidate_map),
        class_map_key = store.map("coupled.assignment_classes", assignment_axis, 2, class_map),
        sample_coordinates_key = store.vector("coupled.assignment_samples", assignment_axis, sample_coordinates),
        variable_coordinates_key = store.vector("coupled.solver_variables", assignment_axis, variable_coordinates);

    struct Solved { std::size_t index; std::optional<Rational> optimum; };
    std::map<std::pair<bool, Vector>, Solved> solved;
    N::Array solvers, certificates, sample_bounds, contrast_bounds;
    Vector sample_lower(target_map.rows()), sample_upper(target_map.rows()),
           contrast_lower(rfc::machine::checked_size_mul(contrasts.size(), targets.size())),
           contrast_upper(contrast_lower.size());
    bool feasible = true;
    for (const auto &measure : measures) {
        for (bool maximize : {false, true}) {
            const auto identity = std::make_pair(maximize, measure.coefficients);
            if (!solved.contains(identity)) {
                problem.direction = maximize ? "maximize" : "minimize";
                auto objective = Json::object(); objective.set("sense", Json::string(problem.direction));
                problem.declaration.set("objective", objective);
                problem.objective = ExactMatrix(1, n);
                for (std::size_t j = 0; j < n; ++j) problem.objective.set(0, j, measure.coefficients[j]);
                auto reservation = resources->reserve(1);
                auto result = allocation_detail::solve_fields(problem);
                const auto index = solvers.elements.size();
                const auto ns = previous ? NativeRecord::decode(previous->at("solvers").at(index).at("constraint_publication").bytes()).source_namespace()
                                         : new_native_record_namespace();
                const auto publication = NativeRecord::create_from_endpoint(std::move(result.fields.declaration), ns,
                    std::move(result.fields.sources), std::move(result.fields.endpoint));
                const auto key = static_cast<std::uint64_t>(native.sources.size() + 1);
                auto owner = publication.archive_source(key); const auto size = owner.bytes().size();
                native.sources.push_back(std::move(owner));
                solvers.elements.emplace_back(N::Object{{{"constraint_publication", N(N::ByteSpan{key, 0, size, false})},
                    {"variable_binding", N(N::Object{{{"source_namespace", N(rfc::numeric::Integer(std::to_string(ns)))},
                        {"grade", N(rfc::numeric::Integer(0))}, {"source_coordinates", variable_coordinates_key},
                        {"meaning", N(std::string("Original ordered assignment variables of this solver source"))}}})}}});
                solved.emplace(identity, Solved{index, result.optimum});
            }
            const auto &result = solved.at(identity);
            certificates.elements.emplace_back(N::Object{{{"solver", N(rfc::numeric::Integer(std::to_string(result.index)))},
                {"kind", N(std::string(measure.contrast ? "contrast" : "sample"))},
                {"output_index", N(rfc::numeric::Integer(std::to_string(measure.output)))},
                {"target_index", N(rfc::numeric::Integer(std::to_string(measure.target)))}}});
            if (!result.optimum) {feasible = false; break;}
            auto &values = measure.contrast ? (maximize ? contrast_upper : contrast_lower)
                                            : (maximize ? sample_upper : sample_lower);
            values.at(measure.output * targets.size() + measure.target) = *result.optimum;
        }
        if (!feasible) break;
    }
    if (previous) require(previous->at("solvers").size() == solvers.elements.size(), "replayed coupled solver family differs");
    const auto &coordinates_tensor = std::get<rfc::numeric::SparseRationalTensor>(member(native.declaration, "target_coordinates").data);
    std::vector<std::size_t> target_coordinates;
    for (std::size_t t = 0; t < targets.size(); ++t)
        target_coordinates.push_back(coordinates_tensor.state(t)->get_num().get_ui());
    const auto source_size = native.endpoint.source.grade({0}).cell_count();
    const auto field = [&](std::string name, const Vector &values, std::size_t begin) {
        Vector mapped(source_size);
        if (feasible)
            for (std::size_t t = 0; t < target_coordinates.size(); ++t) mapped.at(target_coordinates[t]) = values.at(begin + t);
        auto active = target_coordinates; std::sort(active.begin(), active.end());
        active.erase(std::unique(active.begin(), active.end()), active.end());
        if (!feasible) active.clear();
        return store.column(std::move(name), 0, source_size, [&](auto i) {return mapped[i];},
            rfc::numeric::StructuralSupport::sparse(source_size, std::move(active)));
    };
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto prefix = "coupled.sample:" + std::to_string(i);
        sample_bounds.elements.emplace_back(N::Object{{{"id", N(label(samples.at(i).at("id")))},
            {"lower", field(prefix + ".lower", sample_lower, i * targets.size())},
            {"upper", field(prefix + ".upper", sample_upper, i * targets.size())}}});
    }
    for (std::size_t i = 0; i < contrasts.size(); ++i) {
        const auto prefix = "coupled.contrast:" + std::to_string(i);
        N::Object item;
        for (const auto *name : {"id", "left", "right"}) item.members.emplace_back(name, import_record_view(contrasts.at(i).at(name)));
        item.members.emplace_back("lower", field(prefix + ".lower", contrast_lower, i * targets.size()));
        item.members.emplace_back("upper", field(prefix + ".upper", contrast_upper, i * targets.size()));
        contrast_bounds.elements.emplace_back(std::move(item));
    }
    return N(N::Object{{{"status", N(std::string(feasible ? "feasible" : "infeasible"))},
        {"variable_origins", import_record_view(origins)}, {"declaration", import_record_view(declaration)},
        {"target_map", target_map_key}, {"assignment_candidates", candidate_map_key},
        {"assignment_classes", class_map_key}, {"assignment_samples", sample_coordinates_key},
        {"solvers", N(std::move(solvers))},
        {"certificates", N(std::move(certificates))}, {"samples", N(std::move(sample_bounds))},
        {"contrasts", N(std::move(contrast_bounds))}}});
}

} // namespace

Json expression_allocations(const ResourceImages &images, const Json &request,
                            std::shared_ptr<CpuResources> resources,
                            std::optional<AssayFieldsSource> *retained, const NativeRecord *previous_coupling) {
    using N = NativeRecordValue;
    if (!retained) {
        std::optional<AssayFieldsSource> fields;
        expression_allocations(images, request, resources, &fields);
        const auto publication = NativeRecord::create_from_endpoint(std::move(fields->declaration),
            new_native_record_namespace(), std::move(fields->sources), std::move(fields->endpoint));
        return expression_fields_view(publication);
    }
    members(request,
            {"schema", "operation", "samples", "evidence_source", "targets", "contrasts",
             "map_evidence"},
            {"kernel_certificate", "coupling"});
    (void)label(request.at("map_evidence"));
    const auto source = request.at("evidence_source"), samples = request.at("samples");
    const auto mode = label(source.at("mode"));
    require(mode == "transcript_sequences" || mode == "annotation_alignments" ||
                mode == "equivalence_classes",
            "expression evidence mode");
    require(samples.size() > 0, "expression requires samples");
    std::set<std::string> sample_ids, resources_seen, contrast_ids;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto sample = samples.at(i);
        members(sample, {"id", "resource", "scale", "scale_provenance"});
        require(label(sample.at("id")).find_first_of("\t\r\n") == std::string::npos,
                "sample identity must be one table token");
        require(sample_ids.insert(label(sample.at("id"))).second &&
                    resources_seen.insert(label(sample.at("resource"))).second,
                "distinct sample identities and resources; repeated evidence cannot be an "
                "independent sample");
        require(exact_rational(label(sample.at("scale"))) > 0, "positive expression sample scale");
        (void)label(sample.at("scale_provenance"));
    }
    for (std::size_t i = 0; i < request.at("contrasts").size(); ++i) {
        const auto c = request.at("contrasts").at(i);
        members(c, {"id", "left", "right"});
        require(contrast_ids.insert(label(c.at("id"))).second,
                "distinct expression contrast identities");
    }
    auto evidence_records = Json::array(), evidence_hashes = Json::array(),
         canonical_candidates = Json::array();
    std::vector<Json> local;
    std::vector<Vector> local_counts(samples.size()), shared_counts;
    N::Array native_maps;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        Json adapted;
        {
            auto reservation = resources->reserve(1);
            if (mode == "transcript_sequences")
                adapted = sequence_evidence(images, source, samples.at(s), &local_counts[s]);
            else if (mode == "annotation_alignments") {
                members(source, {"mode", "annotation", "feature_type", "group_attribute", "overlap",
                                 "strand_policy", "filter", "candidate_search", "junction_policy"});
                auto q = Json::parse(source.dump());
                json_object_object_del(q.raw(), "mode");
                q.set("schema", request.at("schema"));
                q.set("operation", Json::string("feature_counts"));
                q.set("resource", samples.at(s).at("resource"));
                q.set("unit", Json::string("sequenced_segment"));
                adapted = feature_allocation_evidence(images, q, &local_counts[s]);
                const auto a = adapted.at("allocation_request");
                adapted.set("candidates", a.at("candidates"));
                adapted.set("classes", a.at("classes"));
                adapted.set("counts", a.at("samples").at(0).at("counts"));
                adapted.set("reference_sha256",
                            source.at("candidate_search").at("reference_sha256"));
            } else {
                members(source, {"mode", "candidate_search"});
                const auto &image = images.at(label(samples.at(s).at("resource")));
                require(image.declaration.at("format").text() == "salmon_eq",
                        "unweighted equivalence-class resource required");
                adapted = salmon_classes(image.bytes, &local_counts[s]);
                adapted.set("reference_sha256",
                            source.at("candidate_search").at("reference_sha256"));
            }
        }
        if (s == 0)
            canonical_candidates = adapted.at("candidates");
        // Matching names prove this explicitly declared shared dictionary; order
        // is preserved in every local receipt and the first source defines output order.
        std::set<std::string> expected, actual;
        for (std::size_t i = 0; i < canonical_candidates.size(); ++i)
            expected.insert(label(canonical_candidates.at(i)));
        for (std::size_t i = 0; i < adapted.at("candidates").size(); ++i)
            actual.insert(label(adapted.at("candidates").at(i)));
        require(expected == actual && actual.size() == adapted.at("candidates").size(),
                "identical declared expression candidate dictionaries");
        auto record = Json::object();
        record.set("sample", samples.at(s));
        record.set("evidence", adapted);
        evidence_records.push(record);
        local.push_back(adapted);
        evidence_hashes.push(
            Json::string("sha256:" + sha256(images.at(label(samples.at(s).at("resource"))).bytes)));
    }
    if (source.has("annotation"))
        evidence_hashes.push(
            Json::string("sha256:" + sha256(images.at(label(source.at("annotation"))).bytes)));
    if (source.has("reference_resource"))
        evidence_hashes.push(Json::string(
            "sha256:" + sha256(images.at(label(source.at("reference_resource"))).bytes)));
    std::map<std::string, std::size_t> candidate_ids;
    for (std::size_t i = 0; i < canonical_candidates.size(); ++i)
        candidate_ids.emplace(label(canonical_candidates.at(i)), i);
    std::map<std::vector<std::size_t>, std::size_t> groups;
    auto shared_classes = Json::array(), class_maps = Json::array();
    std::vector<std::vector<std::size_t>> maps;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        std::vector<std::size_t> map;
        auto receipt = Json::array();
        const auto classes = local[s].at("classes");
        for (std::size_t c = 0; c < classes.size(); ++c) {
            std::vector<std::size_t> p;
            const auto candidates = classes.at(c).at("candidates");
            for (std::size_t i = 0; i < candidates.size(); ++i)
                p.push_back(candidate_ids.at(label(candidates.at(i))));
            std::sort(p.begin(), p.end());
            if (!groups.contains(p)) {
                const auto index = groups.size();
                groups.emplace(p, index);
                auto shared = Json::object(), members = Json::array();
                shared.set("id", Json::string("class:" + std::to_string(index)));
                for (auto id : p)
                    members.push(canonical_candidates.at(id));
                shared.set("candidates", members);
                shared_classes.push(shared);
            }
            map.push_back(groups.at(p));
            receipt.push(cardinality(groups.at(p)));
        }
        maps.push_back(std::move(map));
        class_maps.push(receipt);
    }
    auto shared_samples = Json::array();
    for (std::size_t s = 0; s < samples.size(); ++s) {
        ExactMatrix map(shared_classes.size(), maps[s].size());
        require(local_counts[s].size() == maps[s].size(), "local expression count shape");
        for (std::size_t c = 0; c < maps[s].size(); ++c)
            map.set(maps[s][c], c, 1);
        shared_counts.push_back(map.apply(local_counts[s]));
        if (retained) {
            native_maps.elements.emplace_back(N::Object{
                {{"action", SourceGradeFields::matrix(map)},
                 {"local_counts",
                  N(SourceGradeFields::tensor(
                      local_counts[s].size(), [&](auto i) { return local_counts[s][i]; },
                      rfc::numeric::StructuralSupport::full(local_counts[s].size())))},
                 {"class_indices",
                  N(SourceGradeFields::tensor(
                      maps[s].size(), [&](auto i) { return Rational(std::to_string(maps[s][i])); },
                      rfc::numeric::StructuralSupport::full(maps[s].size())))}}});
        }
        auto sample = Json::parse(samples.at(s).dump());
        json_object_object_del(sample.raw(), "resource");
        sample.set("counts", values(shared_counts.back()));
        shared_samples.push(sample);
    }
    auto allocation = Json::object(), search = Json::object();
    allocation.set("schema", Json::string("gadget.assay.allocations.v1"));
    allocation.set("operation", Json::string("allocation_bounds"));
    allocation.set("candidates", canonical_candidates);
    allocation.set("classes", shared_classes);
    allocation.set("samples", shared_samples);
    allocation.set("targets", request.at("targets"));
    allocation.set("contrasts", request.at("contrasts"));
    allocation.set("class_partition", Json::string("disjoint_observations"));
    allocation.set("units", Json::string(mode == "equivalence_classes" ? "assigned_fragment_count"
                                                                       : "assigned_segment_count"));
    allocation.set("evidence", evidence_hashes);
    if (mode == "transcript_sequences") {
        search.set("method",
                   Json::string("Complete native mismatch-bounded transcript placement search"));
        search.set("completeness",
                   Json::string("Every declared strand, transcript and contiguous placement under "
                                "the retained IUPAC rule and mismatch budget"));
        search.set("reference_sha256", local[0].at("reference_sha256"));
    } else
        search = source.at("candidate_search");
    allocation.set("candidate_search", search);
    if (request.has("kernel_certificate"))
        allocation.set("kernel_certificate", request.at("kernel_certificate"));
    (void)allocation_analysis(allocation, resources, retained, &shared_counts, true);
    if (retained) {
        require(retained->has_value(), "native expression allocation missing");
        auto &native = retained->value();
        N::Array native_evidence, originals;
        std::optional<N> coupling;
        if (request.has("coupling"))
            coupling.emplace(coupled_expression(allocation, request.at("coupling"), shared_counts,
                                                native, resources, previous_coupling));
        for (std::size_t s = 0; s < local.size(); ++s) {
            auto protocol = Json::parse(local[s].dump());
            json_object_object_del(protocol.raw(), "counts");
            if (protocol.has("allocation_request"))
                json_object_object_del(protocol.at("allocation_request").at("samples").at(0).raw(),
                                       "counts");
            native_evidence.elements.emplace_back(
                N::Object{{{"sample", import_record_view(samples.at(s))},
                           {"evidence", import_record_view(protocol)}}});
        }
        std::set<std::string> used;
        for (std::size_t s = 0; s < samples.size(); ++s)
            used.insert(label(samples.at(s).at("resource")));
        for (const auto *name : {"annotation", "reference_resource"})
            if (source.has(name))
                used.insert(label(source.at(name)));
        for (const auto &id : used) {
            const auto &image = images.at(id);
            const auto key = static_cast<std::uint64_t>(native.sources.size() + 1);
            native.sources.emplace_back(
                key, std::as_bytes(std::span(image.bytes.data(), image.bytes.size())));
            originals.elements.emplace_back(
                N::Object{{{"id", N(id)},
                           {"role", N(std::string("source"))},
                           {"declaration", import_record_view(image.declaration)},
                           {"bytes", N(N::ByteSpan{key, 0, image.bytes.size(), false})}}});
            if (!image.reference_bytes.empty()) {
                const auto reference_key = static_cast<std::uint64_t>(native.sources.size() + 1);
                native.sources.emplace_back(reference_key,
                                            std::as_bytes(std::span(image.reference_bytes.data(),
                                                                    image.reference_bytes.size())));
                originals.elements.emplace_back(N::Object{
                    {{"id", N(id)},
                     {"role", N(std::string("reference"))},
                     {"declaration", import_record_view(image.declaration)},
                     {"bytes",
                      N(N::ByteSpan{reference_key, 0, image.reference_bytes.size(), false})}}});
            }
        }
        native.declaration =
            N(N::Object{{{"schema", N(std::string("gadget.expression.fields.v1"))},
                         {"request", import_record_view(request)},
                         {"allocation_fields", std::move(native.declaration)},
                         {"class_maps", N(std::move(native_maps))},
                         {"evidence_records", N(std::move(native_evidence))},
                         {"resources", N(std::move(originals))},
                         {"candidate_search", import_record_view(search)},
                         {"evidence_hashes", import_record_view(evidence_hashes)}}});
        if (coupling) std::get<N::Object>(native.declaration.data).members.emplace_back("coupling", std::move(*coupling));
        return Json{};
    }
    throw std::logic_error("native expression result missing");
}

Json allocation_fields_summary_view(const NativeRecord &allocation_fields) {
    auto result = Json::object();
    result.set("bound_scope", render_record_view(allocation_fields.at("summary").at("bound_scope")));
    for (const auto *name : {"native_relation_count", "native_occurrence_count"})
        result.set(name, render_record_view(allocation_fields.at(name)));
    result.set("native_complex_sha256", render_record_view(allocation_fields.at("summary").at("native_complex_sha256")));
    for (const auto *name : {"constraint_rank", "allocation_nullity", "target_ambiguity_rank"})
        result.set(name, render_record_view(allocation_fields.at(name)));
    result.set("kernel_residual_zero", Json::boolean(true));
    if (allocation_fields.at("kernel_factors").size() != 0) {
        auto certificate = Json::object(), factors = Json::array();
        const auto entries = [](const NativeRecord &record) {
            auto rows = Json::array();
            const auto &tensor = record.tensor();
            const auto width = tensor.domain().extents().at(1);
            for (const auto &entry : tensor.values().entries()) {
                auto row = Json::array();
                row.push(cardinality(entry.index / width));
                row.push(cardinality(entry.index % width));
                row.push(Json::string(entry.value.get_str()));
                rows.push(row);
            }
            return rows;
        };
        const auto native = allocation_fields.at("kernel_factors");
        for (std::size_t c = 0; c < native.size(); ++c) {
            const auto factor = native.at(c);
            auto row = Json::object();
            for (const auto *name : {"class_id", "occurrence_offset", "kernel_offset"})
                row.set(name, render_record_view(factor.at(name)));
            row.set("native_relation_identity",
                    cardinality(rfc::machine::checked_size_add(1000000, c)));
            row.set("candidates",
                    render_record_view(allocation_fields.at("classes").at(c).at("candidates")));
            row.set("kernel_dimension",
                    cardinality(factor.at("kernel_frame").tensor().domain().extents().at(1)));
            row.set("kernel_frame", entries(factor.at("kernel_frame")));
            row.set("target_image", entries(factor.at("target_image")));
            factors.push(row);
        }
        certificate.set("factors", factors);
        certificate.set("target_coordinate_order",
                        render_record_view(allocation_fields.at("targets")));
        certificate.set(
            "coordinate_rule",
            Json::string(
                "Class order then candidate order; occurrence key is occurrence_offset+j+1; target "
                "rows follow target_coordinate_order; global kernel column is kernel_offset+k"));
        certificate.set("bound_map_nonzeros",
                        cardinality(rfc::machine::checked_size_add(
                            allocation_fields.at("lower_map").tensor().values().entries().size(),
                            allocation_fields.at("upper_map").tensor().values().entries().size())));
        certificate.set(
            "target_image_nonzeros",
            cardinality(allocation_fields.at("target_image").tensor().values().entries().size()));
        result.set("kernel_certificate", certificate);
    }
    return result;
}

Json coupled_fields_view(const NativeRecord &record) {
    const auto coupling = record.at("coupling"), allocation = record.at("allocation_fields");
    const auto coordinates = allocation.at("target_coordinates").tensor();
    auto out = Json::object(), certificates = Json::array(), samples = Json::array(), contrasts = Json::array();
    out.set("status", Json::string(coupling.at("status").text()));
    out.set("variable_origins", render_record_view(coupling.at("variable_origins")));
    out.set("declaration", render_record_view(coupling.at("declaration")));
    std::vector<NativeRecord> solvers;
    std::vector<Json> results;
    for (std::size_t i = 0; i < coupling.at("solvers").size(); ++i) {
        solvers.push_back(NativeRecord::decode(coupling.at("solvers").at(i).at("constraint_publication").bytes()));
        results.push_back(constraint_fields_view(solvers.back()));
    }
    for (std::size_t i = 0; i < coupling.at("certificates").size(); ++i) {
        const auto declaration = coupling.at("certificates").at(i);
        const auto index = declaration.at("solver").integer().get_ui();
        auto certificate = Json::object();
        certificate.set("request", constraint_request_view(solvers.at(index)));
        certificate.set("result", results.at(index));
        for (const auto *name : {"kind", "output_index", "target_index"})
            certificate.set(name, render_record_view(declaration.at(name)));
        certificates.push(certificate);
    }
    if (coupling.at("status").text() == "feasible")
        for (const auto *name : {"samples", "contrasts"})
            for (std::size_t i = 0; i < coupling.at(name).size(); ++i) {
                const auto source = coupling.at(name).at(i);
                const auto &lower = retained_tensor(record, source.at("lower")).tensor,
                           &upper = retained_tensor(record, source.at("upper")).tensor;
                auto row = Json::object(), bounds = Json::array();
                row.set("id", render_record_view(source.at("id")));
                const bool contrast = std::string_view(name) == "contrasts";
                if (contrast)
                    for (const auto *key : {"left", "right"}) row.set(key, render_record_view(source.at(key)));
                for (std::size_t t = 0; t < allocation.at("targets").size(); ++t) {
                    const auto coordinate = coordinates.state(t)->get_num().get_ui();
                    bounds.push(bound(render_record_view(allocation.at("targets").at(t).at("id")),
                                      lower.state(coordinate).value(), upper.state(coordinate).value(), contrast));
                }
                row.set("bounds", bounds);
                (contrast ? contrasts : samples).push(row);
            }
    out.set("certificates", certificates); out.set("samples", samples); out.set("contrasts", contrasts);
    if (coupling.at("status").text() == "feasible")
        out.set("constraint_units", Json::string("Unscaled assigned counts; positive sample scales apply only to output observables"));
    return out;
}
Json expression_fields_view(const NativeRecord &record) {
    const auto allocation_fields = record.at("allocation_fields"), request = record.at("request");
    const auto targets = allocation_fields.at("targets"), samples = allocation_fields.at("samples");
    const auto coordinates = allocation_fields.at("target_coordinates").tensor();
    const auto retained_field =
        [&](const NativeRecord &sample,
            const char *name) -> const rfc::numeric::SparseRationalTensor & {
        const auto key = sample.at(name).integer();
        for (const auto &entry : record.endpoint()->tensors)
            if (key == rfc::numeric::Integer(std::to_string(entry.key)))
                return entry.tensor;
        throw std::invalid_argument("expression field missing from original publication");
    };
    const auto bounds = [&](const NativeRecord &sample, Rational scale, bool contrast) {
        auto result = Json::array();
        const auto &lower = retained_field(sample, "lower"),
                   &upper = retained_field(sample, "upper");
        for (std::size_t t = 0; t < targets.size(); ++t) {
            const auto coordinate = std::stoull(coordinates.state(t)->get_str());
            result.push(bound(render_record_view(targets.at(t).at("id")),
                              lower.state(coordinate).value() * scale,
                              upper.state(coordinate).value() * scale, contrast));
        }
        return result;
    };
    auto result = Json::object(), raw_samples = Json::array(), expression_samples = Json::array(),
         shared_samples = Json::array(), contrasts = Json::array();
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto sample = samples.at(s);
        auto row = Json::object();
        row.set("id", render_record_view(sample.at("id")));
        row.set("total", Json::string(sample.at("total").rational().get_str()));
        row.set("bounds", bounds(sample, Rational{1} / sample.at("scale").rational(), false));
        raw_samples.push(row);
        auto scaled = Json::parse(row.dump());
        scaled.set("bounds", bounds(sample, Rational{1}, false));
        expression_samples.push(scaled);
        auto shared = render_record_view(request.at("samples").at(s));
        json_object_object_del(shared.raw(), "resource");
        Vector counts;
        const auto &values = retained_field(sample, "class_count");
        for (std::size_t c = 0; c < allocation_fields.at("classes").size(); ++c)
            counts.push_back(values.state(c).value());
        shared.set("counts", assay_detail::values(counts));
        shared_samples.push(shared);
    }
    for (std::size_t c = 0; c < allocation_fields.at("contrasts").size(); ++c) {
        const auto contrast = allocation_fields.at("contrasts").at(c);
        auto row = Json::object();
        for (const auto *name : {"id", "left", "right"})
            row.set(name, render_record_view(contrast.at(name)));
        row.set("bounds", bounds(contrast, Rational{1}, true));
        for (std::size_t i = 0; i < row.at("bounds").size(); ++i)
            json_object_object_del(row.at("bounds").at(i).raw(), "identifiable");
        contrasts.push(row);
    }
    result.set("samples", raw_samples);
    result.set("contrasts", contrasts);
    const auto summary = allocation_fields_summary_view(allocation_fields);
    json_object_object_foreach(summary.raw(), name, value) {
        (void)value;
        result.set(name, summary.at(name));
    }
    auto allocation = Json::object();
    allocation.set("schema", Json::string("gadget.assay.allocations.v1"));
    allocation.set("operation", Json::string("allocation_bounds"));
    for (const auto *name : {"candidates", "classes", "targets", "units"})
        allocation.set(name, render_record_view(allocation_fields.at(name)));
    allocation.set("samples", shared_samples);
    allocation.set("contrasts", render_record_view(request.at("contrasts")));
    allocation.set("class_partition", Json::string("disjoint_observations"));
    allocation.set("candidate_search", render_record_view(record.at("candidate_search")));
    allocation.set("evidence", render_record_view(record.at("evidence_hashes")));
    if (request.has("kernel_certificate"))
        allocation.set("kernel_certificate", render_record_view(request.at("kernel_certificate")));
    auto certificate = Json::object();
    certificate.set("schema", Json::string("gadget.assay.result.v1"));
    certificate.set("query", allocation);
    certificate.set("source_digest", Json::string("sha256:" + sha256(allocation.dump())));
    certificate.set("result", result);
    certificate.set("arithmetic", Json::string("rational"));
    certificate.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    certificate.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
    certificate = seal(certificate, "certificate_sha256");
    auto expression = Json::parse(result.dump());
    expression.set("samples", expression_samples);
    expression.set("status", Json::string("feasible"));
    if (record.has("coupling")) expression = coupled_fields_view(record);
    auto evidence = render_record_view(record.at("evidence_records")), maps = Json::array();
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto local = record.at("class_maps").at(s);
        Vector counts;
        for (std::size_t c = 0; c < local.at("local_counts").tensor().domain().logical_size(); ++c)
            counts.push_back(local.at("local_counts").tensor().state(c).value());
        evidence.at(s).at("evidence").set("counts", values(counts));
        if (evidence.at(s).at("evidence").has("allocation_request"))
            evidence.at(s)
                .at("evidence")
                .at("allocation_request")
                .at("samples")
                .at(0)
                .set("counts", values(counts));
        auto map = Json::array();
        for (std::size_t c = 0; c < local.at("class_indices").tensor().domain().logical_size(); ++c)
            map.push(
                cardinality(std::stoull(local.at("class_indices").tensor().state(c)->get_str())));
        maps.push(map);
    }
    auto out = Json::object();
    out.set("evidence_records", evidence);
    out.set("class_maps", maps);
    out.set("allocation_request", allocation);
    out.set("allocation_certificate", certificate);
    out.set("expression", expression);
    out.set("targets", render_record_view(targets));
    out.set("map_evidence", render_record_view(request.at("map_evidence")));
    out.set("units", render_record_view(allocation_fields.at("units")));
    out.set(
        "scope",
        Json::string("Assignment and contrast bounds on the declared candidate universe; "
                     "transcript-to-feature memberships and optional coupling are declarations. "
                     "Each read record, retained RG/query/mate segment or reported class count is "
                     "an observation. Placement multiplicity is not count multiplicity. Bound "
                     "endpoints are separate observables, not fitted abundances."));
    return out;
}
} // namespace gadget::assay_detail
