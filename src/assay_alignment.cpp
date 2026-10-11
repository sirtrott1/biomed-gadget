#include "sequence_route_witnesses.hpp"
#include "source_grade_fields.hpp"

namespace gadget::assay_detail {
Json alignment_analysis(const ResourceImages &images, const Json &request,
                        std::optional<AssayFieldsSource> *retained) {
    members(request, {"schema", "operation", "query_resource", "query_id", "reference_resource",
                      "reference_id", "scoring", "ambiguity", "mode"});
    require(request.at("mode").text() == "global_linear_gap", "alignment mode");
    const auto ambiguity = request.at("ambiguity").text();
    require(ambiguity == "literal" || ambiguity == "possible" || ambiguity == "definite",
            "alignment ambiguity rule");
    const auto sequence = [&](const Json &resource_id, const Json &record_id) {
        const auto id = text(resource_id);
        require(images.contains(id), "alignment sequence resource");
        const auto &r = images.at(id);
        require(r.declaration.at("format").text() == "fasta", "exact pair alignment uses FASTA");
        std::string bases;
        for (const auto &record : sequence_records(r))
            if (record.id == record_id.text())
                bases = record.bases;
        require(!bases.empty(), "alignment sequence identity not found");
        return bases;
    };
    const auto query = sequence(request.at("query_resource"), request.at("query_id")),
               reference = sequence(request.at("reference_resource"), request.at("reference_id"));
    const auto scoring = request.at("scoring");
    members(scoring, {"match", "mismatch", "query_gap", "reference_gap", "units", "evidence"});
    (void)text(scoring.at("units"), source_byte_limit());
    (void)text(scoring.at("evidence"), source_byte_limit());
    const auto match = exact_rational(text(scoring.at("match"), source_byte_limit())),
               mismatch = exact_rational(text(scoring.at("mismatch"), source_byte_limit())),
               qgap = exact_rational(text(scoring.at("query_gap"), source_byte_limit())),
               rgap = exact_rational(text(scoring.at("reference_gap"), source_byte_limit()));
    require(match >= 0 && mismatch >= 0 && qgap >= 0 && rgap >= 0,
            "native exact route requires nonnegative costs");
    const auto height = rfc::machine::checked_size_add(query.size(), 1),
               width = rfc::machine::checked_size_add(reference.size(), 1),
               n = rfc::machine::checked_size_mul(height, width);
    const auto vertex = [width](std::size_t i, std::size_t j) { return i * width + j; };
    struct Step {
        std::size_t source, target;
        char operation;
    };
    std::vector<Step> steps;
    Vector costs;
    sequence::Routes routes(n, "declared.sequence.edit.step");
    const auto add = [&](std::size_t from, std::size_t to, char op, const Rational &cost) {
        routes.add(from, to, cost);
        steps.push_back({from, to, op});
        costs.push_back(cost);
    };
    for (std::size_t i = 0; i <= query.size(); ++i)
        for (std::size_t j = 0; j <= reference.size(); ++j) {
            const auto from = vertex(i, j);
            if (i < query.size())
                add(from, vertex(i + 1, j), 'I', rgap);
            if (j < reference.size())
                add(from, vertex(i, j + 1), 'D', qgap);
            if (i < query.size() && j < reference.size()) {
                const auto a = iupac_mask(query[i]), b = iupac_mask(reference[j]);
                const bool equal = ambiguity == "literal"    ? query[i] == reference[j]
                                   : ambiguity == "possible" ? (a & b) != 0
                                                             : a == b && (a & (a - 1)) == 0;
                add(from, vertex(i + 1, j + 1), equal ? '=' : 'X', equal ? match : mismatch);
            }
        }
    const auto source = routes.freeze();
    const auto graph = rfc::view::project_pairwise_directed(source, costs);
    const auto start = graph.object_index(1), end = graph.object_index(n);
    const sequence::CompleteRouteWitnesses witnesses(graph, start, end);
    const auto &route = witnesses.chosen;
    const auto &potential = witnesses.prefix, &suffix = witnesses.suffix;
    const auto boundary = rfc::field::boundary_matrix(source, rfc::core::Grade{1});
    witnesses.verify(boundary, costs, start, end);
    std::vector<rfc::path::TraversalStep> traversal;
    for (auto e : route.relations) traversal.push_back({e, Rational{1}});
    const auto realization = rfc::path::realize(source, rfc::core::Grade{1}, std::move(traversal));
    Vector demand(n, 0); demand[start] -= 1; demand[end] += 1;
    require(realization.boundary == demand && rfc::numeric::dot(costs, realization.field) == route.distance,
            "native ordered edit route boundary and cost");
    Vector minimizing(steps.size(), 0), prefix_route(n), prefix_position(n), suffix_route(n), suffix_position(n);
    std::string operations;
    for (std::size_t e = 0; e < steps.size(); ++e) {
        const auto &step = steps[e];
        operations += step.operation;
        if (potential[step.source] + costs[e] + suffix[step.target] == route.distance) minimizing[e] = 1;
    }
    for (std::size_t v = 0; v < n; ++v) {
        prefix_route[v] = Rational(std::to_string(witnesses.prefix_origins[v].route));
        prefix_position[v] = Rational(std::to_string(witnesses.prefix_origins[v].position));
        suffix_route[v] = Rational(std::to_string(witnesses.suffix_origins[v].route));
        suffix_position[v] = Rational(std::to_string(witnesses.suffix_origins[v].position));
    }
    AssayFieldsSource native{rfc::storage::TemporalEndpoint(source), {}, {}};
    for (std::size_t k = 0; k < source.grade_count(); ++k)
        native.endpoint.carriers.push_back({k + 1, rfc::storage::TemporalCarrierKind::SourceGrade,
            rfc::core::Grade{static_cast<std::uint16_t>(k)}, {}, 0,
            rfc::numeric::IndexDomain::linear(source.grade({static_cast<std::uint16_t>(k)}).cell_count())});
    using N = NativeRecordValue;
    const auto integer = [](std::size_t i) {return N(rfc::numeric::Integer(std::to_string(i)));};
    SourceGradeFields store(native.endpoint);
    N::Array fields, cell_ids, histories;
    const auto field = [&](std::string id, std::uint16_t grade, const Vector &values,
                           std::string quantity, std::string units) {
        auto key = store.field(id, grade, values);
        fields.elements.emplace_back(N::Object{{{"id", N(id)}, {"key", key},
            {"quantity", N(quantity)}, {"units", N(units)},
            {"method", N(std::string("Declared edit grid and exact native route witnesses"))}}});
        return key;
    };
    const auto cost_key = field("step_cost", 1, costs, "declared_edit_cost", scoring.at("units").text()),
        chosen_key = field("chosen_path", 1, realization.field, "declared_route_membership", "dimensionless"),
        minimizing_key = field("optimal_support", 1, minimizing, "declared_route_membership", "dimensionless"),
        prefix_key = field("minimum_prefix_cost", 0, potential, "declared_minimum_edit_cost", scoring.at("units").text()),
        suffix_key = field("minimum_suffix_cost", 0, suffix, "declared_minimum_edit_cost", scoring.at("units").text());
    Vector ordered_edges;
    std::size_t begin = 0;
    for (const auto &history : witnesses.histories) {
        histories.elements.emplace_back(N::Object{{{"kind", N(history.kind)}, {"from_vertex", integer(history.from)},
            {"to_vertex", integer(history.to)}, {"minimum_cost", N(history.route.distance)},
            {"begin", integer(begin)}, {"size", integer(history.route.relations.size())}}});
        for (auto e : history.route.relations) ordered_edges.emplace_back(std::to_string(e));
        begin = rfc::machine::checked_size_add(begin, history.route.relations.size());
    }
    const auto occurrence_axis = store.explicit_carrier(ordered_edges.size());
    rfc::numeric::ExactMatrix lift(steps.size(), ordered_edges.size());
    for (std::size_t i = 0; i < ordered_edges.size(); ++i) lift.set(ordered_edges[i].get_num().get_ui(), i, 1);
    N::Array vertices, relations;
    for (std::size_t v = 0; v < n; ++v)
        vertices.elements.emplace_back("prefix:" + std::to_string(v / width) + ":" + std::to_string(v % width));
    for (std::size_t e = 0; e < steps.size(); ++e) relations.elements.emplace_back("edit:" + std::to_string(e));
    cell_ids.elements.emplace_back(std::move(vertices)); cell_ids.elements.emplace_back(std::move(relations));
    native.declaration = N(N::Object{{{"schema", N(std::string("gadget.sequence.alignment.fields.v1"))},
        {"protocol", import_record_view(request)}, {"query_sequence", N(query)}, {"reference_sequence", N(reference)},
        {"operations", N(operations)}, {"fields", N(std::move(fields))}, {"cell_ids", N(std::move(cell_ids))},
        {"step_cost", cost_key}, {"chosen_path", chosen_key}, {"optimal_support", minimizing_key},
        {"minimum_prefix_cost", prefix_key}, {"minimum_suffix_cost", suffix_key},
        {"prefix_route", store.field("prefix_route", 0, prefix_route)},
        {"prefix_position", store.field("prefix_position", 0, prefix_position)},
        {"suffix_route", store.field("suffix_route", 0, suffix_route)},
        {"suffix_position", store.field("suffix_position", 0, suffix_position)},
        {"realized_boundary", store.field("realized_boundary", 0, realization.boundary)},
        {"ordered_relations", store.vector("ordered_relations", occurrence_axis, ordered_edges)},
        {"ordered_source_lift", store.map("ordered_source_lift", occurrence_axis, 2, lift)},
        {"route_witnesses", N(std::move(histories))}}});
    if (retained) {retained->emplace(std::move(native)); return Json::object();}
    const auto publication = NativeRecord::create_from_endpoint(std::move(native.declaration), new_native_record_namespace(),
        std::move(native.sources), std::move(native.endpoint));
    return sequence_alignment_fields_view(publication);
}
Json sequence_alignment_fields_view(const NativeRecord &record) {
    const auto &source = record.endpoint()->source;
    const auto request = render_record_view(record.at("protocol"));
    const auto query = record.at("query_sequence").text(), reference = record.at("reference_sequence").text(),
               ambiguity = request.at("ambiguity").text(), operations = record.at("operations").text();
    const auto scoring = request.at("scoring");
    const auto width = reference.size() + 1, n = source.grade({0}).cell_count();
    const auto values = [&](std::string_view key) {
        const auto &tensor = retained_tensor(record, record.at(key)).tensor;
        Vector result(tensor.domain().logical_size());
        tensor.support().for_each_active([&](auto i) {result[i] = *tensor.state(i);});
        return result;
    };
    const auto costs = values("step_cost"), potential = values("minimum_prefix_cost"), suffix = values("minimum_suffix_cost"),
        chosen = values("chosen_path"), minimizing = values("optimal_support"), ordered = values("ordered_relations"),
        prefix_route = values("prefix_route"), prefix_position = values("prefix_position"),
        suffix_route = values("suffix_route"), suffix_position = values("suffix_position");
    struct Step {std::size_t source, target; char operation;};
    std::vector<Step> steps;
    const auto &grade = source.grade({1});
    for (std::size_t e = 0; e < grade.cell_count(); ++e) {
        const auto [begin, end] = grade.incidence_range(e);
        require(end - begin == 2, "retained edit step incidence");
        steps.push_back({grade.lower_index[begin], grade.lower_index[begin + 1], operations.at(e)});
    }
    auto histories = Json::array();
    std::vector<std::uint32_t> chosen_relations;
    for (std::size_t h = 0; h < record.at("route_witnesses").size(); ++h) {
        const auto history = record.at("route_witnesses").at(h);
        auto item = Json::object(), relations = Json::array();
        const auto begin = history.at("begin").integer().get_ui(), size = history.at("size").integer().get_ui();
        for (std::size_t i = begin; i < begin + size; ++i) {
            const auto e = ordered.at(i).get_num().get_ui();
            relations.push(cardinality(e)); if (h == 0) chosen_relations.push_back(e);
        }
        for (const auto *name : {"kind", "from_vertex", "to_vertex", "minimum_cost"}) item.set(name, render_record_view(history.at(name)));
        item.set("ordered_relations", relations); histories.push(item);
    }
    const auto distance = record.at("route_witnesses").at(0).at("minimum_cost").rational();
    auto potentials = Json::array();
    for (std::size_t v = 0; v < n; ++v) {
        auto p = Json::object();
        p.set("query_prefix", cardinality(v / width)); p.set("reference_prefix", cardinality(v % width));
        p.set("minimum_prefix_cost", Json::string(potential[v].get_str()));
        p.set("minimum_suffix_cost", Json::string(suffix[v].get_str()));
        const auto origin = [&](const Rational &route, const Rational &position) {
            auto out = Json::object(); out.set("route", cardinality(route.get_num().get_ui()));
            out.set("position", cardinality(position.get_num().get_ui())); return out;
        };
        p.set("prefix_witness", origin(prefix_route[v], prefix_position[v]));
        p.set("suffix_witness", origin(suffix_route[v], suffix_position[v])); potentials.push(p);
    }
    auto path = Json::array(), all_optimal = Json::array();
    std::set<std::size_t> optimal_indices,
        chosen_indices(chosen_relations.begin(), chosen_relations.end());
    std::string aligned_query, aligned_reference;
    for (auto relation : chosen_relations) {
        const auto &s = steps[relation];
        const auto i = s.source / width, j = s.source % width;
        aligned_query += s.operation == 'D' ? '-' : query[i];
        aligned_reference += s.operation == 'I' ? '-' : reference[j];
        auto p = Json::object();
        p.set("relation_identity", Json::string(std::to_string(1000000 + relation)));
        p.set("operation", Json::string(std::string(1, s.operation)));
        p.set("query_prefix", cardinality(i));
        p.set("reference_prefix", cardinality(j));
        p.set("cost", Json::string(costs[relation].get_str()));
        path.push(p);
    }
    for (std::size_t k = 0; k < steps.size(); ++k) {
        const auto &s = steps[k];
        if (minimizing[k] == 0)
            continue;
        optimal_indices.insert(k);
        auto p = Json::object();
        p.set("relation_identity", Json::string(std::to_string(1000000 + k)));
        p.set("operation", Json::string(std::string(1, s.operation)));
        p.set("query_prefix", cardinality(s.source / width));
        p.set("reference_prefix", cardinality(s.source % width));
        p.set("cost", Json::string(costs[k].get_str()));
        all_optimal.push(p);
    }
    const auto bytes = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(source), 1, source_codec_limits(source_byte_limit())));
    auto out = Json::object();
    out.set("query_sequence", Json::string(query));
    out.set("reference_sequence", Json::string(reference));
    out.set("aligned_query", Json::string(aligned_query));
    out.set("aligned_reference", Json::string(aligned_reference));
    out.set("minimum_cost", Json::string(distance.get_str()));
    out.set("compatibility_rule",
            Json::string(ambiguity == "literal"    ? "equal_tokens"
                         : ambiguity == "possible" ? "nonempty_IUPAC_intersection"
                                                   : "equal_singleton_IUPAC_sets"));
    out.set("path", path);
    out.set("all_optimal_relations", all_optimal);
    out.set("prefix_suffix_witnesses", potentials);
    out.set("primal_dual_verified", Json::boolean(true));
    out.set("suffix_dual_verified", Json::boolean(true));
    out.set("route_witnesses", histories);
    out.set("native_route_calls", cardinality(histories.size()));
    auto realized = Json::object(), realized_field = Json::array(), realized_order = Json::array();
    for (auto e : chosen_relations) realized_order.push(cardinality(e));
    for (std::size_t e = 0; e < chosen.size(); ++e) if (chosen[e] != 0) {
        auto value = Json::object(); value.set("relation", cardinality(e));
        value.set("coefficient", Json::string(chosen[e].get_str())); realized_field.push(value);
    }
    realized.set("ordered_relations", realized_order); realized.set("realized_field", realized_field);
    realized.set("boundary_verified", Json::boolean(true)); realized.set("cost_verified", Json::boolean(true));
    out.set("realized_path", realized);
    out.set("vertices", cardinality(n));
    out.set("relations", cardinality(steps.size()));
    out.set("native_complex_sha256",
            Json::string(
                "sha256:" +
                sha256(std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    out.set("scope",
            Json::string("Global linear-gap alignment with exact prefix and suffix cost witnesses, "
                         "the complete minimizing step support, and one ordered minimizing path."));
    auto declarations = Json::object(), relations = Json::array(), fields = Json::array();
    for (std::size_t k = 0; k < steps.size(); ++k) {
        const auto &s = steps[k];
        auto r = Json::object(), participants = Json::array();
        r.set("id", Json::string("edit:" + std::to_string(k)));
        r.set("grade", cardinality(1));
        r.set("type", Json::string("declared.sequence.edit.step"));
        r.set("meaning", Json::string("Declared global edit step " + std::string(1, s.operation)));
        for (auto v : {s.source, s.target}) {
            auto p = Json::object();
            p.set("cell", Json::string("prefix:" + std::to_string(v / width) + ":" +
                                       std::to_string(v % width)));
            p.set("exists", Json::boolean(true));
            p.set("negative", Json::boolean(v == s.source));
            p.set("magnitude", Json::string("1"));
            participants.push(p);
        }
        r.set("participants", participants);
        relations.push(r);
    }
    for (const auto &id :
         {std::string("step_cost"), std::string("chosen_path"), std::string("optimal_support")}) {
        auto f = Json::object(), values = Json::array();
        f.set("id", Json::string(id));
        f.set("grade", cardinality(1));
        f.set("channel", Json::string("observation"));
        f.set("quantity",
              Json::string(id == "step_cost" ? "declared_edit_cost" : "declared_route_membership"));
        f.set("units", id == "step_cost" ? scoring.at("units") : Json::string("dimensionless"));
        f.set("method", Json::string("Native exact routes on the complete declared edit grid"));
        for (std::size_t k = 0; k < steps.size(); ++k) {
            auto p = Json::object();
            p.set("relation", Json::string("edit:" + std::to_string(k)));
            p.set("value",
                  Json::string(id == "step_cost"     ? costs[k].get_str()
                               : id == "chosen_path" ? (chosen_indices.contains(k) ? "1" : "0")
                                                     : (optimal_indices.contains(k) ? "1" : "0")));
            values.push(p);
        }
        f.set("values", values);
        fields.push(f);
    }
    declarations.set("relations", relations);
    declarations.set("fields", fields);
    out.set("relation_declaration", declarations);
    return out;
}
} // namespace gadget::assay_detail
