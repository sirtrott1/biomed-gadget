#include "sequence_detail.hpp"

namespace gadget::assay_detail {
AssayFieldsSource indexed_alignment_fields_source(const ResourceImages &images, const Json &request) {
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
    const auto &query_image = images.at(text(request.at("query_resource")));
    const SequenceSource queries(query_image);
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
    N::Array matches, children, summaries, contexts;
    std::vector<PlacementSource> placements;
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(query_image.bytes))),
                                     RawByteSource(2, std::as_bytes(std::span(image.bytes)))};
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
                auto native = route_source(source, routes);
                SourceGradeFields fields(native.endpoint);
                N::Array accepted, native_routes;
                for (std::size_t length = q.bases.size() - errors; length <= width; ++length) {
                    const auto target = coordinates.at({q.bases.size(), length});
                    native_calls = rfc::machine::checked_size_add(native_calls, 1);
                    const auto route = rfc::path::dijkstra_exact(
                        graph, graph.object_index(origin + 1), graph.object_index(target + 1));
                    if (!route.reachable ||
                        route.distance > Rational{static_cast<unsigned long>(errors)})
                        continue;
                    auto m = N::Object{};
                    std::string aligned_query, aligned_reference;
                    for (auto e : route.relations) {
                        const auto [i, j] = positions[routes.steps[e].from];
                        aligned_query += operations[e] == 'D' ? '-' : q.bases[i];
                        aligned_reference += operations[e] == 'I' ? '-' : reference[j];
                    }
                    member(m, "query_id", N(queries.id(q.read)));
                    member(m, "query_record_index", integer(q.read));
                    member(m, "reference_id", N(refs.id(r)));
                    member(m, "strand", N(q.strand));
                    member(m, "start", integer(start));
                    member(m, "end", integer(start + length));
                    member(m, "minimum_edits", N(route.distance));
                    member(m, "aligned_query", N(aligned_query));
                    member(m, "aligned_reference", N(aligned_reference));
                    member(m, "query_evidence", import_record_view(queries.evidence(q.read)));
                    member(m, "reference_evidence", import_record_view(refs.evidence(r)));
                    member(m, "route_index", integer(native_routes.elements.size()));
                    element(native_routes, retain_route(fields, source, routes, route.relations, origin, target, route.distance));
                    element(accepted, m);
                    query_matches[q.read] =
                        rfc::machine::checked_size_add(query_matches[q.read], 1);
                }
                if (accepted.elements.empty())
                    continue;
                N::Object c;
                N::Array coords;
                for (const auto &[i, j] : positions) {
                    auto row = N::Array{};
                    element(row, integer(i));
                    element(row, integer(j));
                    element(coords, row);
                }
                member(c, "prefix_coordinates", coords);
                member(c, "reference_start", integer(start));
                const auto carrier_index = children.elements.size();
                auto &declaration = std::get<N::Object>(native.declaration.data);
                member(declaration, "view", std::move(c));
                member(declaration, "routes", std::move(native_routes));
                member(declaration, "query_prefix", fields.column("query.prefix", 0, positions.size(), [&](auto v) { return Rational(std::to_string(positions[v].first)); }));
                member(declaration, "reference_prefix", fields.column("reference.prefix", 0, positions.size(), [&](auto v) { return Rational(std::to_string(positions[v].second)); }));
                N::Array codes; for (auto code : operations) element(codes, std::string(1, code));
                member(declaration, "operations", std::move(codes));
                const auto child = source_record(std::move(native));
                element(contexts, child_context(child));
                element(children, archive_child(child, sources));
                for (std::size_t i = 0; i < accepted.elements.size(); ++i) {
                    auto m = std::get<N::Object>(std::move(accepted.elements[i].data));
                    const auto end = std::get<rfc::numeric::Integer>(
                        std::find_if(m.members.begin(), m.members.end(),
                            [](const auto &entry) { return entry.first == "end"; })->second.data).get_ui();
                    placements.push_back({q.read, r, {}, carrier_index, i, {{start, end}}});
                    member(m, "carrier", integer(carrier_index));
                    element(matches, m);
                }
            }
        }
        begin = end;
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        auto row = N::Object{};
        member(row, "query_id", N(queries.id(q)));
        member(row, "query_record_index", integer(q));
        member(row, "matches", integer(query_matches[q]));
        element(summaries, row);
    }
    auto out = N::Object{}, work = N::Object{};
    member(work, "reference_base_visits", integer(scans));
    member(work, "retained_seed_hits", integer(seed_hits));
    member(work, "seed_lengths", integer(seed_lengths.size()));
    member(work, "key_encoding", N("injective_two_bit_DNA"));
    member(work, "reference_blocks", integer(blocks));
    member(work, "read_partitions", integer(partitions));
    member(work, "read_partition_records", integer(partition_records));
    member(work, "reference_block_bases", integer(block_bases));
    member(work, "peak_partition_reads", integer(peak_queries));
    member(work, "peak_oriented_query_bases", integer(peak_query_bases));
    member(work, "peak_normalized_reference_block_bases", integer(peak_block));
    member(work, "peak_normalized_candidate_window_bases", integer(peak_window));
    member(work, "source_layout",
             N("Original byte spans and sequence coordinates; seed state "
                          "continues across blocks and source lines, and resets at record "
                          "boundaries and ambiguity."));
    member(work, "candidate_starts", integer(candidate_count));
    member(work, "route_vertices", integer(vertex_count));
    member(work, "route_relations", integer(relation_count));
    member(work, "native_route_calls", integer(native_calls));
    member(out, "matches", matches);
    member(out, "queries", summaries);
    member(out, "work", work);
    member(out, "complete", N(true));
    member(out, "coordinate_convention", N("zero_based_half_open"));
    member(out, "optimality",
            N("rfc::path::dijkstra_exact on every admitted candidate endpoint"));
    member(out, "completeness_argument",
            N(
                "At most e edits disrupt at most e of e+1 disjoint nonempty seeds. A surviving "
                "seed locates the start within e bases of its query offset. Every path of cost at "
                "most e has |query_prefix-reference_prefix|<=e; omitted vertices cannot belong to "
                "an admitted path. Every nonempty endpoint length from m-e to m+e is tested."));
    member(out,
        "search_universe",
        N("All nonempty contiguous reference substrings on the specified strands, with "
                     "a full literal ACGT query and unit-cost substitution, insertion and deletion "
                     "budget. One exact minimizing edit history is retained per placement."));
    return collection_source("gadget.sequence.indexed.fields.v1", request, query_image, image, queries, refs,
                             placements, std::move(out), std::move(children), std::move(contexts), std::move(sources));
}
Json indexed_alignment_fields_view(const NativeRecord &record) {
    return sequence::routed_collection_view(record);
}
Json indexed_alignment(const ResourceImages &images, const Json &request) {
    return indexed_alignment_fields_view(sequence::source_record(indexed_alignment_fields_source(images, request)));
}
} // namespace gadget::assay_detail
