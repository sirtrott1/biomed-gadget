#include "sequence_detail.hpp"

namespace gadget::assay_detail {
AssayFieldsSource splice_search_fields_source(const ResourceImages &images, const Json &request) {
    using namespace sequence;
    members(request, {"schema", "operation", "query_resource", "reference_resource", "strands",
                      "minimum_exon_bases", "minimum_intron_bases", "maximum_intron_bases",
                      "maximum_junctions", "motifs", "evidence"});
    (void)text(request.at("evidence"), source_byte_limit());
    const auto overhang =
        count(request.at("minimum_exon_bases"), std::numeric_limits<std::size_t>::max(),
              "minimum exon representation");
    const auto min_gap =
        count(request.at("minimum_intron_bases"), std::numeric_limits<std::size_t>::max(),
              "minimum intron representation");
    const auto max_gap =
        count(request.at("maximum_intron_bases"), std::numeric_limits<std::size_t>::max(),
              "maximum intron representation");
    const auto max_junctions =
        count(request.at("maximum_junctions"), std::numeric_limits<std::size_t>::max(),
              "splice junction representation");
    require(overhang >= 2 && min_gap >= 4 && min_gap <= max_gap, "splice length declarations");
    std::set<std::string> motifs;
    for (std::size_t i = 0; i < request.at("motifs").size(); ++i) {
        const auto motif = text(request.at("motifs").at(i));
        require(motif.size() == 5 && motif[2] == '-',
                "splice motif requires two DNA pairs, e.g. GT-AG");
        (void)token(motif.substr(0, 2));
        (void)token(motif.substr(3));
        require(motifs.insert(motif).second, "duplicate splice motif");
    }
    require(!motifs.empty(), "splice motif count");
    const auto &query_image = images.at(text(request.at("query_resource")));
    const SequenceSource query_source(query_image);
    const auto queries = records(query_source);
    const auto &reference_image = images.at(text(request.at("reference_resource")));
    const SequenceSource reference_source(reference_image);
    const auto refs = records(reference_source);
    reference_frame(reference_image, refs);
    require(!queries.empty(), "splice query count");
    struct Query {
        std::size_t read;
        std::string strand, bases;
    };
    std::vector<Query> oriented;
    SeedIndex index;
    const auto seed_length = std::min(overhang, token_bases);
    for (std::size_t q = 0; q < queries.size(); ++q) {
        require(queries[q].bases.size() >= overhang, "splice query length");
        for (auto c : queries[q].bases)
            require(dna(c) >= 0, "splice query requires literal ACGT");
        for (const auto &strand : strands(request)) {
            Query v{q, strand,
                    strand == "+" ? queries[q].bases : reverse_complement(queries[q].bases)};
            for (std::size_t i = 0; i + overhang <= v.bases.size(); ++i)
                index.request(v.bases.substr(i, seed_length));
            oriented.push_back(std::move(v));
        }
    }
    index.build(refs);
    struct Exon {
        std::size_t reference, begin, end, start;
    };
    struct Link {
        std::size_t relation, next;
        std::string motif;
    };
    N::Array matches, children, summaries, contexts;
    std::vector<PlacementSource> placements;
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(query_image.bytes))),
                                     RawByteSource(2, std::as_bytes(std::span(reference_image.bytes)))};
    std::vector<std::size_t> counts(queries.size());
    rfc::numeric::Integer total_exons{0}, total_links{0}, visits{0};
    for (const auto &q : oriented) {
        std::vector<Exon> exons;
        std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> beginning;
        for (std::size_t i = 0; i + overhang <= q.bases.size(); ++i)
            for (const auto &hit : index.find(q.bases.substr(i, seed_length))) {
                const auto &r = refs[hit.reference].bases;
                auto length = seed_length;
                while (i + length < q.bases.size() && hit.offset + length < r.size() &&
                       q.bases[i + length] == r[hit.offset + length])
                    ++length;
                for (std::size_t n = overhang; n <= length; ++n) {
                    ++total_exons;
                    beginning[{hit.reference, i}].push_back(exons.size());
                    exons.push_back({hit.reference, i, i + n, hit.offset});
                }
            }
        if (exons.empty())
            continue;
        const auto sink = rfc::machine::checked_size_add(exons.size(), 1);
        Routes routes(rfc::machine::checked_size_add(sink, 1), "sequence.exon_chain.step");
        std::vector<std::vector<Link>> next(exons.size());
        std::vector<std::optional<std::size_t>> entry(exons.size()), exit(exons.size());
        for (std::size_t a = 0; a < exons.size(); ++a) {
            const auto &left = exons[a];
            if (left.begin == 0)
                entry[a] = routes.add(0, a + 1, Rational{0});
            if (left.end == q.bases.size())
                exit[a] = routes.add(a + 1, sink, Rational{0});
            const auto it = beginning.find({left.reference, left.end});
            if (it == beginning.end())
                continue;
            const auto donor = left.start + left.end - left.begin;
            const auto &reference = refs[left.reference].bases;
            if (min_gap > reference.size() - donor)
                continue;
            const auto lower = donor + min_gap,
                       upper = donor + std::min(max_gap, reference.size() - donor);
            // Seed hits are retained in reference-offset order. Exact interval
            // lookup excludes only starts outside the declared intron range.
            const auto begin = std::lower_bound(
                it->second.begin(), it->second.end(), lower,
                [&](std::size_t b, std::size_t start) { return exons[b].start < start; });
            for (auto cursor = begin; cursor != it->second.end() && exons[*cursor].start <= upper;
                 ++cursor) {
                const auto b = *cursor;
                ++visits;
                const auto acceptor = exons[b].start;
                if (acceptor < donor || acceptor - donor < min_gap || acceptor - donor > max_gap)
                    continue;
                // A literal DNA motif cannot match an unresolved or RNA token.
                // Reject this candidate before attempting strand transport;
                // other admissible exon chains must still be evaluated.
                if (dna(reference[donor]) < 0 || dna(reference[donor + 1]) < 0 ||
                    dna(reference[acceptor - 2]) < 0 || dna(reference[acceptor - 1]) < 0)
                    continue;
                const auto motif =
                    q.strand == "+"
                        ? reference.substr(donor, 2) + "-" + reference.substr(acceptor - 2, 2)
                        : reverse_complement(reference.substr(acceptor - 2, 2)) + "-" +
                              reverse_complement(reference.substr(donor, 2));
                if (!motifs.contains(motif))
                    continue;
                ++total_links;
                const auto relation = routes.add(a + 1, b + 1, Rational{1});
                next[a].push_back({relation, b, motif});
            }
        }
        const auto source = routes.freeze();
        const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
        const auto optimum =
            rfc::path::dijkstra_exact(graph, graph.object_index(1), graph.object_index(sink + 1));
        auto native = route_source(source, routes);
        SourceGradeFields fields(native.endpoint);
        N::Object c;
        N::Array exon_records, native_routes;
        for (const auto &exon : exons) {
            auto row = N::Object{};
            member(row, "reference_id", N(refs[exon.reference].id));
            member(row, "query_begin", integer(exon.begin));
            member(row, "query_end", integer(exon.end));
            member(row, "start", integer(exon.start));
            member(row, "end", integer(exon.start + exon.end - exon.begin));
            element(exon_records, row);
        }
        member(c, "exons", exon_records);
        member(c, "source_vertex", integer(0));
        member(c, "sink_vertex", integer(sink));
        member(c, "minimum_junctions",
              optimum.reachable ? N(optimum.distance) : N{});
        const auto carrier_index = children.elements.size();
        std::vector<std::size_t> chain;
        std::vector<std::uint32_t> ordered;
        std::vector<std::string> chain_motifs;
        const auto emit = [&]() {
            auto m = N::Object{}; auto blocks = N::Array{}, junctions = N::Array{};
            std::string spelling;
            for (std::size_t i = 0; i < chain.size(); ++i) {
                const auto &e = exons[chain[i]];
                auto b = std::get<N::Object>(exon_records.elements.at(chain[i]).data);
                member(b, "exon_identity", integer(chain[i]));
                element(blocks, b);
                spelling += refs[e.reference].bases.substr(e.start, e.end - e.begin);
                if (i) {
                    const auto &left = exons[chain[i - 1]];
                    auto j = N::Object{};
                    member(j, "donor", integer(left.start + left.end - left.begin));
                    member(j, "acceptor", integer(e.start));
                    member(j, "motif", N(chain_motifs[i - 1]));
                    member(j, "origin", N("reference_sequence_and_exact_exon_match"));
                    element(junctions, j);
                }
            }
            require(spelling == q.bases, "splice exon gluing must reconstruct the entire query");
            member(m, "query_id", N(queries[q.read].id));
            member(m, "query_record_index", integer(q.read));
            member(m, "reference_id", N(refs[exons[chain.front()].reference].id));
            member(m, "strand", N(q.strand));
            member(m, "blocks", blocks);
            member(m, "junctions", junctions);
            member(m, "query_evidence", import_record_view(queries[q.read].evidence));
            member(m, "reference_evidence", import_record_view(refs[exons[chain.front()].reference].evidence));
            member(m, "carrier", integer(carrier_index));
            member(m, "route_index", integer(native_routes.elements.size()));
            element(native_routes, retain_route(fields, source, routes, ordered, 0, sink,
                                                Rational{std::to_string(junctions.elements.size())}));
            PlacementSource placement{q.read, exons[chain.front()].reference, {},
                carrier_index, native_routes.elements.size() - 1, {}};
            for (const auto e : chain) placement.blocks.emplace_back(
                exons[e].start, exons[e].start + exons[e].end - exons[e].begin);
            placements.push_back(std::move(placement));
            member(m, "query_reconstruction_verified", N(true));
            element(matches, m);
            ++counts[q.read];
        };
        // Enumerate biological exon-chain declarations, preserving all ties.
        // Native routing supplies the independent global minimum on this DAG.
        struct Frame {
            std::size_t exon, next_link;
            bool entered;
        };
        for (std::size_t a = 0; a < exons.size(); ++a) {
            if (!entry[a])
                continue;
            ordered.push_back(static_cast<std::uint32_t>(*entry[a]));
            std::vector<Frame> stack{{a, 0, false}};
            while (!stack.empty()) {
                auto &frame = stack.back();
                if (!frame.entered) {
                    frame.entered = true;
                    ++visits;
                    chain.push_back(frame.exon);
                    if (exit[frame.exon]) {
                        ordered.push_back(static_cast<std::uint32_t>(*exit[frame.exon]));
                        emit();
                        ordered.pop_back();
                    }
                }
                if (chain.size() <= max_junctions && frame.next_link < next[frame.exon].size()) {
                    const auto link = next[frame.exon][frame.next_link++];
                    ordered.push_back(static_cast<std::uint32_t>(link.relation));
                    chain_motifs.push_back(link.motif);
                    stack.push_back({link.next, 0, false});
                } else {
                    stack.pop_back();
                    chain.pop_back();
                    if (!stack.empty()) {
                        ordered.pop_back();
                        chain_motifs.pop_back();
                    }
                }
            }
            ordered.pop_back();
        }
        auto &declaration = std::get<N::Object>(native.declaration.data);
        member(declaration, "view", std::move(c)); member(declaration, "routes", std::move(native_routes));
        member(declaration, "query_record", integer(q.read)); member(declaration, "strand", q.strand);
        std::vector<std::size_t> support;
        for (std::size_t e = 0; e < exons.size(); ++e) support.push_back(e + 1);
        const auto reported = rfc::numeric::StructuralSupport::sparse(source.grade({0}).cell_count(), support);
        member(declaration, "query_begin", fields.column("exon.query_begin", 0, source.grade({0}).cell_count(),
            [&](auto i) { return Rational(std::to_string(exons.at(i - 1).begin)); }, reported));
        member(declaration, "query_end", fields.column("exon.query_end", 0, source.grade({0}).cell_count(),
            [&](auto i) { return Rational(std::to_string(exons.at(i - 1).end)); }, reported));
        member(declaration, "reference_record", fields.column("exon.reference_record", 0, source.grade({0}).cell_count(),
            [&](auto i) { return Rational(std::to_string(exons.at(i - 1).reference)); }, reported));
        member(declaration, "reference_start", fields.column("exon.reference_start", 0, source.grade({0}).cell_count(),
            [&](auto i) { return Rational(std::to_string(exons.at(i - 1).start)); }, reported));
        const auto child = source_record(std::move(native));
        element(contexts, child_context(child));
        element(children, archive_child(child, sources));
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        auto s = N::Object{};
        member(s, "query_id", N(queries[q].id));
        member(s, "query_record_index", integer(q));
        member(s, "matches", integer(counts[q]));
        element(summaries, s);
    }
    N::Object out;
    const auto index_view = index.receipt();
    auto work = std::get<N::Object>(import_record_view(index_view).data);
    member(work, "exon_candidates", N(total_exons));
    member(work, "splice_transitions", N(total_links));
    member(work, "candidate_and_enumeration_visits", N(visits));
    member(out, "matches", matches);
    member(out, "queries", summaries);
    member(out, "work", work);
    member(out, "complete", N(true));
    member(out, "coordinate_convention", N("zero_based_half_open_reference_orientation"));
    member(out,
        "search_universe",
        N("All exact literal ACGT exon chains in the retained reference on the "
                     "specified strands, including unspliced placements, with declared exon/intron "
                     "lengths, transcript-oriented boundary motifs and junction-count limit. "
                     "Junctions are discovered from sequence; no annotation is required."));
    member(out,
        "completeness_argument",
        N(
            "Every admitted exon starts with an indexed exact seed within the declared overhang. "
            "The full overhang and all matching extensions of "
            "that overhang are retained. Every ordered exon pair meeting the query-coordinate, "
            "reference-gap and motif constraints is declared. All complete chains through at most "
            "the declared number of junctions are enumerated, without partial truncation."));
    return collection_source("gadget.sequence.splice.fields.v1", request, query_image, reference_image, query_source, reference_source,
                             placements, std::move(out), std::move(children), std::move(contexts), std::move(sources));
}
Json splice_search_fields_view(const NativeRecord &record) {
    auto result = sequence::routed_collection_view(record);
    for (std::size_t i = 0; i < result.at("matches").size(); ++i) {
        const auto match = result.at("matches").at(i);
        const auto exons = result.at("carriers").at(match.at("carrier").integer()).at("exons");
        for (std::size_t b = 0; b < match.at("blocks").size(); ++b) {
            const auto identity = match.at("blocks").at(b).at("exon_identity");
            exons.at(identity.integer()).set("exon_identity", identity);
        }
    }
    return result;
}
Json splice_search(const ResourceImages &images, const Json &request) {
    return splice_search_fields_view(sequence::source_record(splice_search_fields_source(images, request)));
}
} // namespace gadget::assay_detail
