#pragma once
#include "sequence_detail.hpp"

namespace gadget::assay_detail::sequence {
struct RouteWitnessOrigin {
    std::size_t route{}, position{};
};
struct CompleteRouteWitnesses {
    Vector prefix, suffix;
    std::vector<RouteWitnessOrigin> prefix_origins, suffix_origins;
    Json histories = Json::array();
    rfc::path::RouteResult chosen;
    std::vector<bool> suffix_reachable;
    std::vector<std::size_t> unreachable_suffix_vertices;

    // All objects must be reachable from start. Optional dead-end suffixes are
    // explicitly absent. Reuse is
    // confined to subpaths of actual native minimum routes, not relaxations or
    // a replacement shortest-path implementation.
    CompleteRouteWitnesses(const rfc::view::DirectedGraphProjection &graph, std::uint32_t start,
                           std::uint32_t end, bool permit_dead_ends = false,
                           std::optional<rfc::path::RouteResult> initial = {}) {
        const auto n = graph.object_count();
        prefix.resize(n);
        suffix.resize(n);
        prefix_origins.resize(n);
        suffix_origins.resize(n);
        std::vector<bool> prefix_known(n, false), suffix_known(n, false);
        const auto record = [&](const rfc::path::RouteResult &route, std::uint32_t from,
                                std::uint32_t to, const std::string &kind) {
            require(route.reachable, "complete native route witness reachability");
            const auto witness = histories.size();
            auto history = Json::object(), relations = Json::array();
            std::vector<std::uint32_t> points{from};
            Vector costs{Rational{0}};
            require(route.projected_edges.size() == route.relations.size(),
                    "native route source lift length");
            for (std::size_t k = 0; k < route.projected_edges.size(); ++k) {
                const auto &edge = graph.edge(route.projected_edges[k]);
                require(edge.source_object == points.back() && edge.measure >= 0 &&
                            edge.source_relation == route.relations[k],
                        "native route ordered source continuity");
                points.push_back(edge.target_object);
                // Scalars are canonical SDK rationals. The route, source lift
                // and dual action remain native; this records its contraction.
                costs.push_back(costs.back() + edge.measure);
                relations.push(cardinality(edge.source_relation));
            }
            require(points.back() == to && costs.back() == route.distance,
                    "native route endpoint and cost witness");
            for (std::size_t k = 0; k < points.size(); ++k) {
                const auto v = points[k];
                if (from == start) {
                    if (prefix_known[v])
                        require(prefix[v] == costs[k], "native optimal prefix consistency");
                    else {
                        prefix_known[v] = true;
                        prefix[v] = costs[k];
                        prefix_origins[v] = {witness, k};
                    }
                }
                if (to == end) {
                    const Rational cost = route.distance - costs[k];
                    if (suffix_known[v])
                        require(suffix[v] == cost, "native optimal suffix consistency");
                    else {
                        suffix_known[v] = true;
                        suffix[v] = cost;
                        suffix_origins[v] = {witness, k};
                    }
                }
            }
            history.set("kind", Json::string(kind));
            history.set("from_vertex", cardinality(from));
            history.set("to_vertex", cardinality(to));
            history.set("minimum_cost", Json::string(route.distance.get_str()));
            history.set("ordered_relations", relations);
            histories.push(history);
        };
        chosen = initial ? std::move(*initial) : rfc::path::dijkstra_exact(graph, start, end);
        record(chosen, start, end, "pair");
        // Farthest declared coordinates first offer longer reusable prefixes.
        // The ordering changes work, never the admitted source or output set.
        for (std::size_t v = n; v-- > 0;)
            if (!prefix_known[v]) {
                const auto route = rfc::path::dijkstra_exact(graph, start, v);
                record(route, start, v, "prefix");
            }
        for (std::size_t v = 0; v < n; ++v)
            if (!suffix_known[v]) {
                const auto route = rfc::path::dijkstra_exact(graph, v, end);
                if (permit_dead_ends && !route.reachable)
                    unreachable_suffix_vertices.push_back(v);
                else
                    record(route, v, end, "suffix");
            }
        require(std::all_of(prefix_known.begin(), prefix_known.end(), [](bool v) { return v; }) &&
                    (permit_dead_ends || std::all_of(suffix_known.begin(), suffix_known.end(),
                                                     [](bool v) { return v; })),
                "complete native prefix and suffix witnesses");
        suffix_reachable = std::move(suffix_known);
    }
    void verify(const rfc::numeric::ExactMatrix &boundary, const Vector &costs, std::uint32_t start,
                std::uint32_t end) const {
        Vector demand(prefix.size(), Rational{0});
        demand[start] -= 1;
        demand[end] += 1;
        const auto field = chosen.materialize_source_field();
        Vector negative_suffix;
        for (const auto &q : suffix)
            negative_suffix.push_back(-q);
        require(rfc::path::shortest_path_certificate(boundary, costs, field, demand, prefix),
                "native complete prefix primal/dual certificate");
        if (unreachable_suffix_vertices.empty()) {
            require(rfc::path::shortest_path_certificate(boundary, costs, field, demand,
                                                         negative_suffix),
                    "native complete suffix primal/dual certificate");
        } else {
            // Dead-end suffix values are absent. The viable vertex/edge
            // restriction is boundary-closed and retains original source maps.
            std::vector<std::size_t> vertices, edges;
            std::vector<std::optional<std::size_t>> vertex_map(prefix.size());
            for (std::size_t v = 0; v < prefix.size(); ++v)
                if (suffix_reachable[v]) {
                    vertex_map[v] = vertices.size();
                    vertices.push_back(v);
                }
            std::vector<bool> viable(costs.size(), true);
            for (const auto &[coordinate, coefficient] : boundary.entries()) {
                (void)coefficient;
                if (!suffix_reachable[coordinate.first])
                    viable[coordinate.second] = false;
            }
            std::vector<std::optional<std::size_t>> edge_map(costs.size());
            for (std::size_t e = 0; e < costs.size(); ++e)
                if (viable[e]) {
                    edge_map[e] = edges.size();
                    edges.push_back(e);
                }
            rfc::numeric::ExactMatrix restricted(vertices.size(), edges.size());
            for (const auto &[coordinate, coefficient] : boundary.entries())
                if (edge_map[coordinate.second])
                    restricted.set(*vertex_map[coordinate.first], *edge_map[coordinate.second],
                                   coefficient);
            Vector restricted_cost, restricted_field, restricted_demand, restricted_potential;
            for (auto e : edges) {
                restricted_cost.push_back(costs[e]);
                restricted_field.push_back(field[e]);
            }
            for (auto v : vertices) {
                restricted_demand.push_back(demand[v]);
                restricted_potential.push_back(negative_suffix[v]);
            }
            require(rfc::path::shortest_path_certificate(restricted, restricted_cost,
                                                         restricted_field, restricted_demand,
                                                         restricted_potential),
                    "native viable suffix primal/dual certificate");
        }
        const auto prefix_action = boundary.transpose().apply(prefix),
                   suffix_action = boundary.transpose().apply(negative_suffix);
        for (std::size_t i = 0; i < histories.size(); ++i) {
            const auto history = histories.at(i);
            for (std::size_t k = 0; k < history.at("ordered_relations").size(); ++k) {
                const auto e = history.at("ordered_relations").at(k).integer();
                if (history.at("from_vertex").integer() == start)
                    require(prefix_action.at(e) == costs.at(e), "tight native prefix witness");
                if (history.at("to_vertex").integer() == end)
                    require(suffix_action.at(e) == costs.at(e), "tight native suffix witness");
            }
        }
    }
};
} // namespace gadget::assay_detail::sequence
