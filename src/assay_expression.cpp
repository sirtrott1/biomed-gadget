#include "genomic_detail.hpp"
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

Json coupled_expression(const Json &allocation, const Json &declaration,
                        const std::shared_ptr<CpuResources> &resources) {
    members(declaration, {"constraints", "witness_representation", "evidence"});
    (void)label(declaration.at("evidence"));
    const auto samples = allocation.at("samples"), classes = allocation.at("classes"),
               targets = allocation.at("targets"), contrasts = allocation.at("contrasts");
    std::map<std::string, std::size_t> sample_ids, target_ids;
    std::vector<std::set<std::string>> target_members;
    std::vector<Rational> scales;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        sample_ids.emplace(label(samples.at(s).at("id")), s);
        scales.push_back(exact_rational(label(samples.at(s).at("scale"))));
    }
    for (std::size_t t = 0; t < targets.size(); ++t) {
        target_ids.emplace(label(targets.at(t).at("id")), t);
        std::set<std::string> members;
        for (std::size_t i = 0; i < targets.at(t).at("candidates").size(); ++i)
            members.insert(label(targets.at(t).at("candidates").at(i)));
        target_members.push_back(std::move(members));
    }
    auto variables = Json::array(), origins = Json::array(), rows = Json::array();
    std::vector<std::vector<std::size_t>> conservation;
    std::vector<Rational> totals;
    std::vector<std::pair<std::size_t, std::string>> coordinates;
    for (std::size_t s = 0; s < samples.size(); ++s)
        for (std::size_t c = 0; c < classes.size(); ++c) {
            std::vector<std::size_t> columns;
            const auto total = samples.at(s).at("counts").at(c);
            const auto members = classes.at(c).at("candidates");
            for (std::size_t i = 0; i < members.size(); ++i) {
                const auto id = Json::string("allocation:" + std::to_string(s) + ":" +
                                             std::to_string(c) + ":" + std::to_string(i));
                auto v = Json::object(), origin = Json::object();
                v.set("id", id);
                v.set("lower", Json::string("0"));
                v.set("upper", total);
                columns.push_back(variables.size());
                variables.push(v);
                coordinates.emplace_back(s, label(members.at(i)));
                origin.set("id", id);
                origin.set("sample", samples.at(s).at("id"));
                origin.set("class", classes.at(c).at("id"));
                origin.set("candidate", members.at(i));
                origins.push(origin);
            }
            conservation.push_back(std::move(columns));
            totals.push_back(exact_rational(label(total)));
        }
    const auto n = variables.size();
    ExactMatrix target_map(rfc::machine::checked_size_mul(samples.size(), targets.size()), n);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t t = 0; t < targets.size(); ++t)
            if (target_members[t].contains(coordinates[j].second))
                target_map.set(coordinates[j].first * targets.size() + t, j, 1);
    const auto target_transpose = target_map.transpose();
    const auto row = [&](std::size_t i) {
        Vector selector(target_map.rows(), 0);
        selector.at(i) = 1;
        return target_transpose.apply(selector);
    };
    for (std::size_t i = 0; i < conservation.size(); ++i) {
        Vector coefficients(n, 0);
        for (auto j : conservation[i])
            coefficients[j] = 1;
        auto r = Json::object();
        r.set("id", Json::string("conservation:" + std::to_string(i)));
        r.set("coefficients", values(coefficients));
        r.set("sense", Json::string("eq"));
        r.set("rhs", Json::string(totals[i].get_str()));
        rows.push(r);
    }
    std::set<std::string> constraint_ids;
    const auto declared = declaration.at("constraints");
    for (std::size_t i = 0; i < declared.size(); ++i) {
        const auto r = declared.at(i);
        members(r, {"id", "terms", "sense", "rhs", "evidence"});
        require(constraint_ids.insert(label(r.at("id"))).second,
                "duplicate expression constraint identity");
        (void)label(r.at("evidence"));
        Vector coefficients(n, 0);
        for (std::size_t k = 0; k < r.at("terms").size(); ++k) {
            const auto term = r.at("terms").at(k);
            members(term, {"sample", "target", "coefficient"});
            require(sample_ids.contains(label(term.at("sample"))) &&
                        target_ids.contains(label(term.at("target"))),
                    "coupling term requires declared sample and target");
            coefficients = rfc::numeric::add(
                coefficients,
                rfc::numeric::scale(exact_rational(label(term.at("coefficient"))),
                                    row(sample_ids.at(label(term.at("sample"))) * targets.size() +
                                        target_ids.at(label(term.at("target"))))));
        }
        auto translated = Json::object();
        translated.set("id", Json::string("declared:" + label(r.at("id"))));
        translated.set("coefficients", values(coefficients));
        translated.set("sense", r.at("sense"));
        translated.set("rhs", r.at("rhs"));
        rows.push(translated);
    }
    std::vector<Measure> measures;
    auto sample_results = Json::array(), contrast_results = Json::array();
    for (std::size_t s = 0; s < samples.size(); ++s) {
        auto item = Json::object();
        item.set("id", samples.at(s).at("id"));
        item.set("bounds", Json::array());
        sample_results.push(item);
        for (std::size_t t = 0; t < targets.size(); ++t)
            measures.push_back({targets.at(t).at("id"), s, t, false,
                                rfc::numeric::scale(scales[s], row(s * targets.size() + t))});
    }
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto item = contrasts.at(c);
        auto result = Json::parse(item.dump());
        result.set("bounds", Json::array());
        contrast_results.push(result);
        const auto l = sample_ids.at(label(item.at("left"))),
                   r = sample_ids.at(label(item.at("right")));
        for (std::size_t t = 0; t < targets.size(); ++t)
            measures.push_back(
                {targets.at(t).at("id"), c, t, true,
                 rfc::numeric::sub(rfc::numeric::scale(scales[l], row(l * targets.size() + t)),
                                   rfc::numeric::scale(scales[r], row(r * targets.size() + t)))});
    }
    auto problem = Json::object();
    problem.set("schema", Json::string("gadget.assays.query.v1"));
    problem.set("operation", Json::string("constrained_allocations"));
    problem.set("domain", Json::string("bounded_integer"));
    problem.set("variables", variables);
    problem.set("constraints", rows);
    problem.set("units", allocation.at("units"));
    problem.set("evidence", declaration.at("evidence"));
    problem.set("witness_representation", declaration.at("witness_representation"));
    std::map<std::string, Json> solved_requests;
    auto certificates = Json::array(), out = Json::object();
    out.set("variable_origins", origins);
    out.set("declaration", declaration);
    for (const auto &m : measures) {
        Rational low{0}, high{0};
        for (bool maximize : {false, true}) {
            auto objective = Json::object();
            objective.set("sense", Json::string(maximize ? "maximize" : "minimize"));
            objective.set("coefficients", values(m.coefficients));
            problem.set("objective", objective);
            const auto request_bytes = problem.dump();
            if (!solved_requests.contains(request_bytes)) {
                auto reservation = resources->reserve(1);
                solved_requests.emplace(request_bytes, constrained_allocations(problem));
            }
            const auto solved = solved_requests.at(request_bytes);
            auto certificate = Json::object();
            certificate.set("request", Json::parse(problem.dump()));
            certificate.set("result", solved);
            certificate.set("kind", Json::string(m.contrast ? "contrast" : "sample"));
            certificate.set("output_index", cardinality(m.output));
            certificate.set("target_index", cardinality(m.target));
            certificates.push(certificate);
            if (solved.at("status").text() == "infeasible") {
                out.set("status", Json::string("infeasible"));
                out.set("certificates", certificates);
                out.set("samples", Json::array());
                out.set("contrasts", Json::array());
                return out;
            }
            (maximize ? high : low) = exact_rational(label(solved.at("minimum_or_maximum")));
        }
        auto item = (m.contrast ? contrast_results : sample_results).at(m.output);
        item.at("bounds").push(bound(m.id, low, high, m.contrast));
    }
    out.set("status", Json::string("feasible"));
    out.set("samples", sample_results);
    out.set("contrasts", contrast_results);
    out.set("certificates", certificates);
    out.set(
        "constraint_units",
        Json::string(
            "Unscaled assigned counts; positive sample scales apply only to output observables"));
    return out;
}
} // namespace

Json expression_allocations(const ResourceImages &images, const Json &request,
                            std::shared_ptr<CpuResources> resources,
                            std::optional<AssayFieldsSource> *retained) {
    using N = NativeRecordValue;
    require(!retained || !request.has("coupling"),
            "coupled expression has no native field producer");
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
    const auto certificate =
        allocation_analysis(allocation, resources, retained, &shared_counts, true);
    if (retained) {
        require(retained->has_value(), "native expression allocation missing");
        auto &native = retained->value();
        N::Array native_evidence, originals;
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
        return Json{};
    }
    auto expression = Json::parse(certificate.at("result").dump());
    expression.set("status", Json::string("feasible"));
    // Allocation samples report raw counts. Output expression samples apply the
    // exact declared scale, as allocation contrasts already do.
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto bounds = expression.at("samples").at(s).at("bounds");
        Vector lo, hi;
        for (std::size_t t = 0; t < bounds.size(); ++t) {
            lo.push_back(exact_rational(label(bounds.at(t).at("lower"))));
            hi.push_back(exact_rational(label(bounds.at(t).at("upper"))));
        }
        const auto scale = exact_rational(label(samples.at(s).at("scale")));
        lo = rfc::numeric::scale(scale, lo);
        hi = rfc::numeric::scale(scale, hi);
        auto scaled = Json::array();
        for (std::size_t t = 0; t < bounds.size(); ++t)
            scaled.push(bound(bounds.at(t).at("id"), lo[t], hi[t]));
        expression.at("samples").at(s).set("bounds", scaled);
    }
    if (request.has("coupling"))
        expression = coupled_expression(allocation, request.at("coupling"), resources);
    auto out = Json::object();
    out.set("evidence_records", evidence_records);
    out.set("class_maps", class_maps);
    out.set("allocation_request", allocation);
    out.set("allocation_certificate", certificate);
    out.set("expression", expression);
    out.set("targets", request.at("targets"));
    out.set("map_evidence", request.at("map_evidence"));
    out.set("units", allocation.at("units"));
    out.set(
        "scope",
        Json::string("Assignment and contrast bounds on the declared candidate universe; "
                     "transcript-to-feature memberships and optional coupling are declarations. "
                     "Each read record, retained RG/query/mate segment or reported class count is "
                     "an observation. Placement multiplicity is not count multiplicity. Bound "
                     "endpoints are separate observables, not fitted abundances."));
    return out;
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
