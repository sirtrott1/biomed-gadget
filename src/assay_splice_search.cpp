#include "sequence_detail.hpp"

namespace gadget::assay_detail {
Json splice_search(const ResourceImages &images, const Json &request) {
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
    const auto queries = sequence_records(images.at(text(request.at("query_resource"))));
    const auto &reference_image = images.at(text(request.at("reference_resource")));
    const auto refs = sequence_records(reference_image);
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
    auto matches = Json::array(), carriers = Json::array(), summaries = Json::array();
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
        auto c = carrier(source, routes), exon_records = Json::array();
        for (const auto &exon : exons) {
            auto row = Json::object();
            row.set("reference_id", Json::string(refs[exon.reference].id));
            row.set("query_begin", cardinality(exon.begin));
            row.set("query_end", cardinality(exon.end));
            row.set("start", cardinality(exon.start));
            row.set("end", cardinality(exon.start + exon.end - exon.begin));
            exon_records.push(row);
        }
        c.set("exons", exon_records);
        c.set("source_vertex", cardinality(0));
        c.set("sink_vertex", cardinality(sink));
        c.set("minimum_junctions",
              optimum.reachable ? Json::string(optimum.distance.get_str()) : Json{});
        const auto carrier_index = carriers.size();
        carriers.push(c);
        std::vector<std::size_t> chain;
        std::vector<std::uint32_t> ordered;
        std::vector<std::string> chain_motifs;
        const auto emit = [&]() {
            auto m = Json::object(), blocks = Json::array(), junctions = Json::array();
            std::string spelling;
            for (std::size_t i = 0; i < chain.size(); ++i) {
                const auto &e = exons[chain[i]];
                auto b = exon_records.at(chain[i]);
                b.set("exon_identity", cardinality(chain[i]));
                blocks.push(b);
                spelling += refs[e.reference].bases.substr(e.start, e.end - e.begin);
                if (i) {
                    const auto &left = exons[chain[i - 1]];
                    auto j = Json::object();
                    j.set("donor", cardinality(left.start + left.end - left.begin));
                    j.set("acceptor", cardinality(e.start));
                    j.set("motif", Json::string(chain_motifs[i - 1]));
                    j.set("origin", Json::string("reference_sequence_and_exact_exon_match"));
                    junctions.push(j);
                }
            }
            require(spelling == q.bases, "splice exon gluing must reconstruct the entire query");
            m.set("query_id", Json::string(queries[q.read].id));
            m.set("query_record_index", cardinality(q.read));
            m.set("reference_id", Json::string(refs[exons[chain.front()].reference].id));
            m.set("strand", Json::string(q.strand));
            m.set("blocks", blocks);
            m.set("junctions", junctions);
            m.set("query_evidence", queries[q.read].evidence);
            m.set("reference_evidence", refs[exons[chain.front()].reference].evidence);
            m.set("carrier", cardinality(carrier_index));
            m.set("route", realized_receipt(source, routes, ordered, 0, sink,
                                            Rational{std::to_string(junctions.size())}));
            m.set("query_reconstruction_verified", Json::boolean(true));
            matches.push(m);
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
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        auto s = Json::object();
        s.set("query_id", Json::string(queries[q].id));
        s.set("query_record_index", cardinality(q));
        s.set("matches", cardinality(counts[q]));
        summaries.push(s);
    }
    auto out = Json::object(), work = index.receipt();
    work.set("exon_candidates", exact_count(total_exons));
    work.set("splice_transitions", exact_count(total_links));
    work.set("candidate_and_enumeration_visits", exact_count(visits));
    out.set("matches", matches);
    out.set("queries", summaries);
    out.set("carriers", carriers);
    out.set("work", work);
    out.set("complete", Json::boolean(true));
    out.set("coordinate_convention", Json::string("zero_based_half_open_reference_orientation"));
    out.set(
        "search_universe",
        Json::string("All exact literal ACGT exon chains in the retained reference on the "
                     "specified strands, including unspliced placements, with declared exon/intron "
                     "lengths, transcript-oriented boundary motifs and junction-count limit. "
                     "Junctions are discovered from sequence; no annotation is required."));
    out.set(
        "completeness_argument",
        Json::string(
            "Every admitted exon starts with an indexed exact seed within the declared overhang. "
            "The full overhang and all matching extensions of "
            "that overhang are retained. Every ordered exon pair meeting the query-coordinate, "
            "reference-gap and motif constraints is declared. All complete chains through at most "
            "the declared number of junctions are enumerated, without partial truncation."));
    return out;
}
} // namespace gadget::assay_detail
