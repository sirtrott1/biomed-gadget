#include "genomic_detail.hpp"
#include <rfc/core/builder.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/path/objectives.hpp>

namespace gadget::assay_detail {
Json assembly_graph(const ResourceImage &resource, const Json &request) {
    using namespace rfc;
    using numeric::ExactMatrix;
    members(request, {"schema", "operation", "resource", "k", "ambiguity", "evidence"}, {"route"});
    require(text(request.at("ambiguity")) == "literal_tokens",
            "assembly retains literal IUPAC tokens");
    (void)text(request.at("evidence"), source_byte_limit());
    const auto k = count(request.at("k"), std::numeric_limits<std::size_t>::max(),
                         "assembly k representation");
    require(k >= 2, "assembly k must be at least two");
    const auto reads = sequence_records(resource);
    require(!reads.empty(), "assembly read count limit");
    std::map<std::string, std::size_t> node_index, edge_index;
    std::vector<std::string> nodes, edges;
    struct Occurrence {
        std::size_t read, offset, edge, source, target;
    };
    std::vector<Occurrence> occurrences;
    std::vector<std::pair<std::size_t, std::size_t>> endpoints;
    const auto admit = [](auto &index, auto &ids, const std::string &id) {
        if (!index.contains(id)) {
            require(ids.size() < std::numeric_limits<std::uint32_t>::max(),
                    "assembly native cell coordinate representation");
            index.emplace(id, ids.size());
            ids.push_back(id);
        }
        return index.at(id);
    };
    for (std::size_t i = 0; i < reads.size(); ++i) {
        require(reads[i].bases.size() >= k, "assembly read is shorter than k");
        for (std::size_t offset = 0; offset <= reads[i].bases.size() - k; ++offset) {
            const auto token = reads[i].bases.substr(offset, k);
            const auto edge = admit(edge_index, edges, token),
                       from = admit(node_index, nodes, token.substr(0, k - 1)),
                       to = admit(node_index, nodes, token.substr(1));
            if (edge == endpoints.size())
                endpoints.emplace_back(from, to);
            require(endpoints.at(edge) == std::pair{from, to},
                    "literal k-mer source endpoint consistency");
            occurrences.push_back({i, offset, edge, from, to});
        }
    }
    ExactMatrix count_map(edges.size(), occurrences.size());
    for (std::size_t i = 0; i < occurrences.size(); ++i)
        count_map.set(occurrences[i].edge, i, Rational{1});
    const auto counts = count_map.apply(Vector(occurrences.size(), Rational{1}));
    core::RelationSignature signature;
    const auto type = signature.add_type("sequence.literal_kmer.transition");
    core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < nodes.size(); ++i)
        builder.declare_witness_identity(i + 1);
    auto relations = Json::array(), values = Json::array(), edge_records = Json::array();
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [from, to] = endpoints[e];
        const auto identity = machine::checked_size_add(1000000, e),
                   occurrence = machine::checked_size_mul(2, e);
        builder.add_grade_one(core::RelationKey{identity}, type,
                              {{from + 1,
                                true,
                                true,
                                Rational{1},
                                {},
                                core::OccurrenceKey{machine::checked_size_add(occurrence, 1)}},
                               {to + 1,
                                true,
                                false,
                                Rational{1},
                                {},
                                core::OccurrenceKey{machine::checked_size_add(occurrence, 2)}}});
        auto r = Json::object(), participants = Json::array(), v = Json::object(),
             record = Json::object();
        const auto id = "kmer:" + std::to_string(e);
        r.set("id", Json::string(id));
        r.set("meaning", Json::string("Literal retained k-mer " + edges[e]));
        r.set("type", Json::string("sequence.literal_kmer.transition"));
        for (const auto &[node, negative] :
             std::vector<std::pair<std::size_t, bool>>{{from, true}, {to, false}}) {
            auto p = Json::object();
            p.set("cell", Json::string("prefix:" + nodes[node]));
            p.set("exists", Json::boolean(true));
            p.set("negative", Json::boolean(negative));
            p.set("magnitude", Json::string("1"));
            participants.push(p);
        }
        r.set("participants", participants);
        relations.push(r);
        v.set("relation", Json::string(id));
        v.set("value", Json::string(counts[e].get_str()));
        values.push(v);
        record.set("id", Json::string(id));
        record.set("sequence", Json::string(edges[e]));
        record.set("source", cardinality(from));
        record.set("target", cardinality(to));
        record.set("count", Json::string(counts[e].get_str()));
        edge_records.push(record);
    }
    const auto source = builder.freeze();
    const auto graph = view::project_pairwise_directed(source);
    std::vector<std::vector<std::size_t>> incoming_edges(nodes.size()),
        outgoing_edges(nodes.size()), read_occurrences(reads.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {
        outgoing_edges[endpoints[e].first].push_back(e);
        incoming_edges[endpoints[e].second].push_back(e);
    }
    for (std::size_t i = 0; i < occurrences.size(); ++i)
        read_occurrences[occurrences[i].read].push_back(i);
    auto read_walks = Json::array(), node_records = Json::array(), node_ids = Json::array(),
         node_counts = Json::array();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        auto n = Json::object(), incoming = Json::array(), outgoing = Json::array();
        for (auto e : incoming_edges[i])
            incoming.push(cardinality(e));
        for (auto e : outgoing_edges[i])
            outgoing.push(cardinality(e));
        n.set("id", Json::string("prefix:" + nodes[i]));
        n.set("sequence", Json::string(nodes[i]));
        n.set("incoming", incoming);
        n.set("outgoing", outgoing);
        n.set("branch", Json::boolean(incoming.size() > 1 || outgoing.size() > 1));
        node_records.push(n);
        node_ids.push(n.at("id"));
        node_counts.push(Json::string(std::to_string(outgoing.size())));
    }
    for (std::size_t i = 0; i < reads.size(); ++i) {
        auto walk = Json::object(), steps = Json::array();
        std::string spelling = reads[i].bases.substr(0, k - 1);
        for (auto occurrence_index : read_occurrences[i]) {
            const auto &occurrence = occurrences[occurrence_index];
            auto e = Json::object();
            e.set("edge", cardinality(occurrence.edge));
            e.set("sequence_offset", cardinality(occurrence.offset));
            steps.push(e);
            spelling += edges[occurrence.edge].back();
        }
        require(spelling == reads[i].bases, "assembly read-walk reconstruction");
        walk.set("read_id", Json::string(reads[i].id));
        walk.set("evidence", reads[i].evidence);
        walk.set("steps", steps);
        walk.set("reconstructed_sequence", Json::string(spelling));
        read_walks.push(walk);
    }
    auto out = Json::object(), declaration = Json::object(), fields = Json::array(),
         f = Json::object();
    f.set("id", Json::string("kmer_occurrence_count"));
    f.set("grade", cardinality(1));
    f.set("channel", Json::string("observation"));
    f.set("quantity", Json::string("observed_kmer_occurrences"));
    f.set("units", Json::string("occurrences"));
    f.set("method", Json::string("Native exact occurrence-to-literal-kmer count map"));
    f.set("values", values);
    fields.push(f);
    declaration.set("relations", relations);
    declaration.set("fields", fields);
    out.set("nodes", node_records);
    out.set("edges", edge_records);
    out.set("read_walks", read_walks);
    out.set("relation_declaration", declaration);
    out.set("source_occurrences", cardinality(occurrences.size()));
    out.set("count_operator_nonzeros", cardinality(count_map.nonzeros()));
    out.set("read_reconstruction_verified", Json::boolean(true));
    out.set("cycle_frame_dimension",
            cardinality(edges.size() -
                        field::materialize_boundary_matrix(source, core::Grade{1}).rank()));
    if (request.has("route")) {
        const auto route = request.at("route");
        members(route, {"start", "end", "objective"});
        require(text(route.at("objective")) == "fewest_steps", "assembly route objective");
        const auto start = text(route.at("start"), source_byte_limit()),
                   end = text(route.at("end"), source_byte_limit());
        require(node_index.contains(start) && node_index.contains(end),
                "assembly route node identity");
        const auto result = path::bfs_min_hops(graph, graph.object_index(node_index.at(start) + 1),
                                               graph.object_index(node_index.at(end) + 1));
        auto receipt = Json::object(), steps = Json::array();
        std::string spelling = start;
        for (auto e : result.relations) {
            steps.push(cardinality(e));
            spelling += edges[e].back();
        }
        receipt.set("found", Json::boolean(result.reachable));
        receipt.set("edges", steps);
        if (result.reachable)
            receipt.set("sequence", Json::string(spelling));
        receipt.set("meaning",
                    Json::string("One native minimum-hop walk between declared endpoints; overlaps "
                                 "alone do not identify a biological contig"));
        out.set("route", receipt);
    }
    out.set("feature_ids", Json::parse("[\"distinct_outgoing_kmers\"]"));
    out.set("observation_ids", node_ids);
    auto rows = Json::array();
    rows.push(node_counts);
    out.set("rows", rows);
    out.set("quantity", Json::string("distinct_outgoing_kmers"));
    out.set("units", Json::string("literal_sequences"));
    out.set("scope",
            Json::string(
                "All observed literal k-mers and read walks, with repeated source occurrences, "
                "self-loops and ambiguous continuations retained; no inferred faces or consensus"));
    return out;
}
} // namespace gadget::assay_detail
