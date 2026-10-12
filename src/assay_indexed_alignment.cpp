#include "sequence_detail.hpp"

namespace gadget::assay_detail {
Json indexed_alignment(const ResourceImages &images, const Json &request) {
    using namespace sequence;
    members(request,
            {"schema", "operation", "query_resource", "reference_resource", "strands",
             "maximum_edits", "scoring", "evidence"},
            {"execution"});
    require(text(request.at("scoring")) == "unit_edit_cost",
            "indexed search requires unit edit costs");
    (void)text(request.at("evidence"), 16384);
    const auto errors = count(request.at("maximum_edits"), std::numeric_limits<std::size_t>::max(),
                              "indexed edit budget");
    const auto orientations = strands(request);
    const SequenceSource queries(images.at(text(request.at("query_resource"))));
    const auto &image = images.at(text(request.at("reference_resource")));
    const SequenceSource refs(image);
    reference_frame(image, refs);
    require(queries.size() > 0, "indexed query source is empty");
    std::size_t partition_records = 256, block_bases = 65536;
    if (request.has("execution")) {
        const auto execution = request.at("execution");
        members(execution, {"read_partition_records", "reference_block_bases"});
        partition_records = count(execution.at("read_partition_records"),
                                  std::numeric_limits<std::size_t>::max(), "read partition width");
        block_bases = count(execution.at("reference_block_bases"),
                            std::numeric_limits<std::size_t>::max(), "reference block width");
        require(partition_records > 0 && block_bases > 0,
                "indexed execution partition and block widths must be positive");
    }
    struct Query {
        std::size_t read;
        std::string strand, bases;
        std::vector<std::pair<std::size_t, std::string>> seeds;
    };
    auto matches = Json::array(), carriers = Json::array(), summaries = Json::array();
    std::vector<std::size_t> query_matches(queries.size());
    std::size_t candidate_count = 0, vertex_count = 0, relation_count = 0, native_calls = 0;
    std::size_t scans = 0, seed_hits = 0, blocks = 0, partitions = 0, peak_block = 0;
    std::size_t peak_queries = 0, peak_query_bases = 0, peak_window = 0;
    std::set<std::size_t> seed_lengths;
    for (std::size_t begin = 0; begin < queries.size();) {
        const auto end = begin + std::min(partition_records, queries.size() - begin);
        std::vector<Query> oriented;
        SeedIndex index;
        std::size_t query_bases = 0;
        for (std::size_t q = begin; q < end; ++q) {
            const auto read = queries.record(q);
            require(read.bases.size() > errors, "indexed query length/edit bound");
            for (auto c : read.bases)
                require(dna(c) >= 0, "indexed query requires literal ACGT");
            for (const auto &strand : orientations) {
                Query v{q, strand, strand == "+" ? read.bases : reverse_complement(read.bases), {}};
                const auto parts = rfc::machine::checked_size_add(errors, 1);
                const auto length = std::min(token_bases, v.bases.size() / parts);
                for (std::size_t s = 0; s <= errors; ++s) {
                    const auto offset = rfc::machine::checked_size_mul(s, v.bases.size()) / parts;
                    const auto seed = v.bases.substr(offset, length);
                    index.request(seed);
                    seed_lengths.insert(length);
                    v.seeds.emplace_back(offset, seed);
                }
                query_bases = rfc::machine::checked_size_add(query_bases, v.bases.size());
                oriented.push_back(std::move(v));
            }
        }
        partitions = rfc::machine::checked_size_add(partitions, 1);
        peak_queries = std::max(peak_queries, end - begin);
        peak_query_bases = std::max(peak_query_bases, query_bases);
        index.build(refs, block_bases);
        const auto counts = index.receipt();
        scans = rfc::machine::checked_size_add(scans, count(counts.at("reference_base_visits"),
                                                            std::numeric_limits<std::size_t>::max(),
                                                            "reference visits"));
        seed_hits = rfc::machine::checked_size_add(
            seed_hits, count(counts.at("retained_seed_hits"),
                             std::numeric_limits<std::size_t>::max(), "seed hits"));
        blocks = rfc::machine::checked_size_add(
            blocks, count(counts.at("reference_blocks"), std::numeric_limits<std::size_t>::max(),
                          "reference blocks"));
        peak_block = std::max(peak_block, count(counts.at("peak_normalized_reference_block_bases"),
                                                std::numeric_limits<std::size_t>::max(),
                                                "normalized reference block"));
        for (const auto &q : oriented) {
            std::set<std::pair<std::size_t, std::size_t>> candidates;
            for (const auto &[offset, seed] : q.seeds)
                for (const auto &hit : index.find(seed)) {
                    const auto last = refs.length(hit.reference) - 1;
                    std::size_t lo = 0, hi = 0;
                    if (hit.offset >= offset) {
                        const auto center = hit.offset - offset;
                        lo = center > errors ? center - errors : 0;
                        hi = center + std::min(errors, last - center);
                    } else {
                        const auto deficit = offset - hit.offset;
                        if (deficit > errors)
                            continue;
                        hi = std::min(last, errors - deficit);
                    }
                    for (std::size_t start = lo;; ++start) {
                        candidates.emplace(hit.reference, start);
                        if (start == hi)
                            break;
                    }
                }
            for (const auto &[r, start] : candidates) {
                candidate_count = rfc::machine::checked_size_add(candidate_count, 1);
                const auto remaining = refs.length(r) - start;
                const auto width =
                    std::min(q.bases.size(), remaining) +
                    std::min(errors, remaining - std::min(q.bases.size(), remaining));
                if (width < q.bases.size() - errors)
                    continue;
                const auto reference = refs.bases(r, start, width);
                peak_window = std::max(peak_window, reference.size());
                std::map<std::pair<std::size_t, std::size_t>, std::size_t> coordinates;
                std::vector<std::pair<std::size_t, std::size_t>> positions;
                for (std::size_t i = 0; i <= q.bases.size(); ++i) {
                    const auto lo = i > errors ? i - errors : 0;
                    const auto hi =
                        std::min(i, width) + std::min(errors, width - std::min(i, width));
                    for (std::size_t j = lo; j <= hi; ++j) {
                        require(positions.size() < std::numeric_limits<std::uint32_t>::max(),
                                "indexed route native vertex coordinate representation");
                        coordinates.emplace(std::pair{i, j}, positions.size());
                        positions.emplace_back(i, j);
                    }
                }
                vertex_count = rfc::machine::checked_size_add(vertex_count, positions.size());
                Routes routes(positions.size(), "sequence.unit_edit.step");
                std::vector<char> operations;
                const auto add = [&](std::size_t from, std::size_t i, std::size_t j, char op,
                                     int cost) {
                    const auto it = coordinates.find({i, j});
                    if (it == coordinates.end() || (errors == 0 && cost != 0))
                        return;
                    routes.add(from, it->second, Rational{cost});
                    operations.push_back(op);
                };
                for (std::size_t v = 0; v < positions.size(); ++v) {
                    const auto [i, j] = positions[v];
                    if (i < q.bases.size())
                        add(v, i + 1, j, 'I', 1);
                    if (j < width)
                        add(v, i, j + 1, 'D', 1);
                    if (i < q.bases.size() && j < width) {
                        const bool equal = q.bases[i] == reference[j];
                        add(v, i + 1, j + 1, equal ? '=' : 'X', equal ? 0 : 1);
                    }
                }
                relation_count =
                    rfc::machine::checked_size_add(relation_count, routes.steps.size());
                // An all-mismatch zero-budget carrier has no grade-one relations.
                if (routes.steps.empty())
                    continue;
                const auto source = routes.freeze();
                const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
                const auto origin = coordinates.at({0, 0});
                auto accepted = Json::array();
                for (std::size_t length = q.bases.size() - errors; length <= width; ++length) {
                    const auto target = coordinates.at({q.bases.size(), length});
                    native_calls = rfc::machine::checked_size_add(native_calls, 1);
                    const auto route = rfc::path::dijkstra_exact(
                        graph, graph.object_index(origin + 1), graph.object_index(target + 1));
                    if (!route.reachable ||
                        route.distance > Rational{static_cast<unsigned long>(errors)})
                        continue;
                    auto m = Json::object();
                    std::string aligned_query, aligned_reference;
                    for (auto e : route.relations) {
                        const auto [i, j] = positions[routes.steps[e].from];
                        aligned_query += operations[e] == 'D' ? '-' : q.bases[i];
                        aligned_reference += operations[e] == 'I' ? '-' : reference[j];
                    }
                    m.set("query_id", Json::string(queries.id(q.read)));
                    m.set("query_record_index", cardinality(q.read));
                    m.set("reference_id", Json::string(refs.id(r)));
                    m.set("strand", Json::string(q.strand));
                    m.set("start", cardinality(start));
                    m.set("end", cardinality(start + length));
                    m.set("minimum_edits", Json::string(route.distance.get_str()));
                    m.set("aligned_query", Json::string(aligned_query));
                    m.set("aligned_reference", Json::string(aligned_reference));
                    m.set("query_evidence", queries.evidence(q.read));
                    m.set("reference_evidence", refs.evidence(r));
                    m.set("route", realized_receipt(source, routes, route.relations, origin, target,
                                                    route.distance));
                    accepted.push(m);
                    query_matches[q.read] =
                        rfc::machine::checked_size_add(query_matches[q.read], 1);
                }
                if (accepted.size() == 0)
                    continue;
                auto c = carrier(source, routes), coords = Json::array();
                for (const auto &[i, j] : positions) {
                    auto row = Json::array();
                    row.push(cardinality(i));
                    row.push(cardinality(j));
                    coords.push(row);
                }
                c.set("prefix_coordinates", coords);
                c.set("reference_start", cardinality(start));
                const auto carrier_index = carriers.size();
                carriers.push(c);
                for (std::size_t i = 0; i < accepted.size(); ++i) {
                    auto m = accepted.at(i);
                    m.set("carrier", cardinality(carrier_index));
                    matches.push(m);
                }
            }
        }
        begin = end;
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        auto row = Json::object();
        row.set("query_id", Json::string(queries.id(q)));
        row.set("query_record_index", cardinality(q));
        row.set("matches", cardinality(query_matches[q]));
        summaries.push(row);
    }
    auto out = Json::object(), work = Json::object();
    work.set("reference_base_visits", cardinality(scans));
    work.set("retained_seed_hits", cardinality(seed_hits));
    work.set("seed_lengths", cardinality(seed_lengths.size()));
    work.set("key_encoding", Json::string("injective_two_bit_DNA"));
    work.set("reference_blocks", cardinality(blocks));
    work.set("read_partitions", cardinality(partitions));
    work.set("read_partition_records", cardinality(partition_records));
    work.set("reference_block_bases", cardinality(block_bases));
    work.set("peak_partition_reads", cardinality(peak_queries));
    work.set("peak_oriented_query_bases", cardinality(peak_query_bases));
    work.set("peak_normalized_reference_block_bases", cardinality(peak_block));
    work.set("peak_normalized_candidate_window_bases", cardinality(peak_window));
    work.set("source_layout",
             Json::string("Original byte spans and sequence coordinates; seed state "
                          "continues across blocks and source lines, and resets at record "
                          "boundaries and ambiguity."));
    work.set("candidate_starts", cardinality(candidate_count));
    work.set("route_vertices", cardinality(vertex_count));
    work.set("route_relations", cardinality(relation_count));
    work.set("native_route_calls", cardinality(native_calls));
    out.set("matches", matches);
    out.set("queries", summaries);
    out.set("carriers", carriers);
    out.set("work", work);
    out.set("complete", Json::boolean(true));
    out.set("coordinate_convention", Json::string("zero_based_half_open"));
    out.set("optimality",
            Json::string("rfc::path::dijkstra_exact on every admitted candidate endpoint"));
    out.set("completeness_argument",
            Json::string(
                "At most e edits disrupt at most e of e+1 disjoint nonempty seeds. A surviving "
                "seed locates the start within e bases of its query offset. Every path of cost at "
                "most e has |query_prefix-reference_prefix|<=e; omitted vertices cannot belong to "
                "an admitted path. Every nonempty endpoint length from m-e to m+e is tested."));
    out.set(
        "search_universe",
        Json::string("All nonempty contiguous reference substrings on the specified strands, with "
                     "a full literal ACGT query and unit-cost substitution, insertion and deletion "
                     "budget. One exact minimizing edit history is retained per placement."));
    return out;
}
} // namespace gadget::assay_detail
