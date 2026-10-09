#include "sequence_route_witnesses.hpp"

namespace gadget::assay_detail {
Json alignment_analysis(const ResourceImages &images, const Json &request) {
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
    auto potentials = Json::array();
    for (std::size_t v = 0; v < n; ++v) {
        const auto index = graph.object_index(v + 1);
        require(index == v, "native edit-grid source coordinate order");
        auto p = Json::object();
        p.set("query_prefix", cardinality(v / width));
        p.set("reference_prefix", cardinality(v % width));
        p.set("minimum_prefix_cost", Json::string(potential[index].get_str()));
        p.set("minimum_suffix_cost", Json::string(suffix[index].get_str()));
        const auto origin = [](const sequence::RouteWitnessOrigin &o) {
            auto item = Json::object();
            item.set("route", cardinality(o.route));
            item.set("position", cardinality(o.position));
            return item;
        };
        p.set("prefix_witness", origin(witnesses.prefix_origins[index]));
        p.set("suffix_witness", origin(witnesses.suffix_origins[index]));
        potentials.push(p);
    }
    const auto boundary = rfc::field::boundary_matrix(source, rfc::core::Grade{1});
    witnesses.verify(boundary, costs, start, end);
    auto path = Json::array(), all_optimal = Json::array();
    std::set<std::size_t> optimal_indices,
        chosen_indices(route.relations.begin(), route.relations.end());
    std::string aligned_query, aligned_reference;
    for (auto relation : route.relations) {
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
        const auto from = graph.object_index(s.source + 1), to = graph.object_index(s.target + 1);
        if (potential[from] + costs[k] + suffix[to] != route.distance)
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
    out.set("minimum_cost", Json::string(route.distance.get_str()));
    out.set("compatibility_rule",
            Json::string(ambiguity == "literal"    ? "equal_tokens"
                         : ambiguity == "possible" ? "nonempty_IUPAC_intersection"
                                                   : "equal_singleton_IUPAC_sets"));
    out.set("path", path);
    out.set("all_optimal_relations", all_optimal);
    out.set("prefix_suffix_witnesses", potentials);
    out.set("primal_dual_verified", Json::boolean(true));
    out.set("suffix_dual_verified", Json::boolean(true));
    out.set("route_witnesses", witnesses.histories);
    out.set("native_route_calls", cardinality(witnesses.histories.size()));
    out.set("realized_path",
            sequence::realized_receipt(source, routes, route.relations, 0, n - 1, route.distance));
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
