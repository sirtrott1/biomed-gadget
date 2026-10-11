#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
#include <rfc/core/builder.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/path/objectives.hpp>

namespace gadget::assay_detail {
Json assembly_graph(const ResourceImage &resource, const Json &request,
                    std::optional<AssayFieldsSource> *retained) {
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
    builder.declare_grade(core::Grade{1});
    for (std::size_t i = 0; i < nodes.size(); ++i)
        builder.declare_witness_identity(i + 1);
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
    for (std::size_t i = 0; i < reads.size(); ++i) {
        std::string spelling = reads[i].bases.substr(0, k - 1);
        for (auto occurrence : read_occurrences[i]) spelling += edges[occurrences[occurrence].edge].back();
        require(spelling == reads[i].bases, "assembly read-walk reconstruction");
    }
    AssayFieldsSource native{storage::TemporalEndpoint(source), {}, {}};
    for (std::size_t grade = 0; grade < source.grade_count(); ++grade)
        native.endpoint.carriers.push_back({grade + 1, storage::TemporalCarrierKind::SourceGrade,
            core::Grade{static_cast<std::uint16_t>(grade)}, {}, 0,
            numeric::IndexDomain::linear(source.grade({static_cast<std::uint16_t>(grade)}).cell_count())});
    using N = NativeRecordValue;
    const auto integer = [](std::size_t i) {return N(numeric::Integer(std::to_string(i)));};
    SourceGradeFields store(native.endpoint);
    const auto occurrence_axis = store.explicit_carrier(occurrences.size());
    Vector read_indices, offsets, edge_indices, outgoing(nodes.size());
    for (const auto &occurrence : occurrences) {
        read_indices.emplace_back(std::to_string(occurrence.read)); offsets.emplace_back(std::to_string(occurrence.offset));
        edge_indices.emplace_back(std::to_string(occurrence.edge));
    }
    for (std::size_t i = 0; i < nodes.size(); ++i) outgoing[i] = Rational(std::to_string(outgoing_edges[i].size()));
    const auto counts_key = store.field("kmer_occurrence_count", 1, counts),
               outgoing_key = store.field("distinct_outgoing_kmers", 0, outgoing);
    N::Array native_reads, node_ids, edge_ids, node_tokens, edge_tokens;
    for (const auto &read : reads)
        native_reads.elements.emplace_back(N::Object{{{"id", N(read.id)}, {"evidence", import_record_view(read.evidence)}}});
    for (const auto &node : nodes) {node_tokens.elements.emplace_back(node); node_ids.elements.emplace_back("prefix:" + node);}
    for (std::size_t e = 0; e < edges.size(); ++e) {edge_tokens.elements.emplace_back(edges[e]); edge_ids.elements.emplace_back("kmer:" + std::to_string(e));}
    N::Object declaration{{{"schema", N(std::string("gadget.sequence.assembly.fields.v1"))},
        {"protocol", import_record_view(request)}, {"node_tokens", N(std::move(node_tokens))},
        {"edge_tokens", N(std::move(edge_tokens))}, {"reads", N(std::move(native_reads))},
        {"cell_ids", N(N::Array{{N(std::move(node_ids)), N(std::move(edge_ids))}})},
        {"fields", N(N::Array{{N(N::Object{{{"id", N(std::string("kmer_occurrence_count"))}, {"key", counts_key},
            {"quantity", N(std::string("observed_kmer_occurrences"))}, {"units", N(std::string("occurrences"))},
            {"method", N(std::string("Exact source occurrence to literal kmer map"))}}}),
            N(N::Object{{{"id", N(std::string("distinct_outgoing_kmers"))}, {"key", outgoing_key},
            {"quantity", N(std::string("distinct_outgoing_kmers"))}, {"units", N(std::string("literal_sequences"))},
            {"method", N(std::string("Declared outgoing kmer relations"))}}})}})},
        {"kmer_occurrence_count", counts_key}, {"distinct_outgoing_kmers", outgoing_key},
        {"count_map", store.map("count_map", occurrence_axis, 2, count_map)},
        {"source_read", store.vector("source_read", occurrence_axis, read_indices)},
        {"source_offset", store.vector("source_offset", occurrence_axis, offsets)},
        {"source_edge", store.vector("source_edge", occurrence_axis, edge_indices)},
        {"cycle_frame_dimension", integer(edges.size() - field::materialize_boundary_matrix(source, core::Grade{1}).rank())}}};
    if (request.has("route")) {
        const auto route = request.at("route"); members(route, {"start", "end", "objective"});
        require(text(route.at("objective")) == "fewest_steps", "assembly route objective");
        const auto from = text(route.at("start"), source_byte_limit()), to = text(route.at("end"), source_byte_limit());
        require(node_index.contains(from) && node_index.contains(to), "assembly route node identity");
        const auto result = path::bfs_min_hops(graph, graph.object_index(node_index.at(from) + 1), graph.object_index(node_index.at(to) + 1));
        Vector ordered;
        for (auto e : result.relations) ordered.emplace_back(std::to_string(e));
        const auto axis = store.explicit_carrier(ordered.size());
        declaration.members.emplace_back("route", N(N::Object{{{"found", N(result.reachable)},
            {"ordered_relations", store.vector("route.ordered_relations", axis, ordered)}}}));
    }
    native.declaration = N(std::move(declaration));
    if (retained) {retained->emplace(std::move(native)); return Json::object();}
    const auto publication = NativeRecord::create_from_endpoint(std::move(native.declaration), new_native_record_namespace(),
        std::move(native.sources), std::move(native.endpoint));
    return sequence_assembly_fields_view(publication);
}
Json sequence_assembly_fields_view(const NativeRecord &record) {
    using namespace rfc;
    const auto &source = record.endpoint()->source;
    const auto request = render_record_view(record.at("protocol"));
    const auto k = request.at("k").integer();
    std::vector<std::string> nodes, edges;
    for (std::size_t i = 0; i < record.at("node_tokens").size(); ++i) nodes.push_back(record.at("node_tokens").at(i).text());
    for (std::size_t i = 0; i < record.at("edge_tokens").size(); ++i) edges.push_back(record.at("edge_tokens").at(i).text());
    const auto values_for = [&](std::string_view key) {
        const auto &tensor = retained_tensor(record, record.at(key)).tensor;
        Vector values(tensor.domain().logical_size());
        tensor.support().for_each_active([&](auto i) {values[i] = *tensor.state(i);}); return values;
    };
    const auto counts = values_for("kmer_occurrence_count"), source_read = values_for("source_read"),
               source_offset = values_for("source_offset"), source_edge = values_for("source_edge");
    std::vector<std::pair<std::size_t, std::size_t>> endpoints;
    const auto &grade = source.grade({1});
    for (std::size_t e = 0; e < grade.cell_count(); ++e) {
        const auto [begin, end] = grade.incidence_range(e);
        require(end - begin == 2, "retained kmer incidence");
        endpoints.emplace_back(grade.lower_index[begin], grade.lower_index[begin + 1]);
    }
    struct Occurrence {std::size_t read, offset, edge, source, target;};
    std::vector<Occurrence> occurrences;
    std::vector<SequenceRecord> reads;
    for (std::size_t i = 0; i < record.at("reads").size(); ++i) {
        const auto read = record.at("reads").at(i);
        reads.push_back({read.at("id").text(), {}, render_record_view(read.at("evidence"))});
    }
    for (std::size_t i = 0; i < source_read.size(); ++i) {
        const auto e = source_edge[i].get_num().get_ui();
        occurrences.push_back({source_read[i].get_num().get_ui(), source_offset[i].get_num().get_ui(), e, endpoints[e].first, endpoints[e].second});
    }
    std::vector<std::vector<std::size_t>> incoming_edges(nodes.size()), outgoing_edges(nodes.size()), read_occurrences(reads.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {outgoing_edges[endpoints[e].first].push_back(e); incoming_edges[endpoints[e].second].push_back(e);}
    for (std::size_t i = 0; i < occurrences.size(); ++i) read_occurrences[occurrences[i].read].push_back(i);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        require(!read_occurrences[i].empty(), "retained read occurrence family");
        auto &spelling = reads[i].bases;
        spelling = nodes[occurrences[read_occurrences[i][0]].source];
        for (auto occurrence : read_occurrences[i]) spelling += edges[occurrences[occurrence].edge].back();
    }
    auto relations = Json::array(), values = Json::array(), edge_records = Json::array();
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [from, to] = endpoints[e];
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
    out.set("count_operator_nonzeros", cardinality(retained_tensor(record, record.at("count_map")).tensor.values().entries().size()));
    out.set("read_reconstruction_verified", Json::boolean(true));
    out.set("cycle_frame_dimension", render_record_view(record.at("cycle_frame_dimension")));
    if (record.has("route")) {
        const auto route = record.at("route");
        const auto &tensor = retained_tensor(record, route.at("ordered_relations")).tensor;
        auto receipt = Json::object(), steps = Json::array();
        std::string spelling = request.at("route").at("start").text();
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i) {
            const auto e = tensor.state(i)->get_num().get_ui(); steps.push(cardinality(e)); spelling += edges[e].back();
        }
        receipt.set("found", Json::boolean(route.at("found").boolean())); receipt.set("edges", steps);
        if (route.at("found").boolean()) receipt.set("sequence", Json::string(spelling));
        receipt.set("meaning", Json::string("One native minimum-hop walk between declared endpoints; overlaps alone do not identify a biological contig"));
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
