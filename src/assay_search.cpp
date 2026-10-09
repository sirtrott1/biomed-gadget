#include "sequence_detail.hpp"
#include <algorithm>
#include <cctype>
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail {
std::vector<SequenceRecord> sequence_records(const ResourceImage &resource) {
    const auto format = resource.declaration.at("format").text();
    require(format == "fasta" || format == "fastq", "sequence resource requires FASTA or FASTQ");
    const auto records = inspect_resource(resource.bytes, resource.declaration).at("records");
    const auto resource_sha256 = "sha256:" + sha256(resource.bytes);
    std::vector<SequenceRecord> out;
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto record = records.at(i);
        const auto start = record.at("source_start").integer(),
                   end = record.at("source_end").integer();
        const auto bytes = resource.bytes.substr(start, end - start);
        const auto ls = lines(bytes);
        std::string bases;
        for (std::size_t l = 1; l < (format == "fastq" ? 2u : ls.size()); ++l)
            for (auto c : bytes.substr(ls[l].first, ls[l].second - ls[l].first))
                bases += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        auto evidence = Json::object();
        evidence.set("source_sha256", Json::string(resource_sha256));
        evidence.set("record_index", cardinality(i));
        evidence.set("source_start", record.at("source_start"));
        evidence.set("source_end", record.at("source_end"));
        out.push_back({record.at("id").text(), std::move(bases), evidence});
    }
    return out;
}
Json reference_matches(const ResourceImages &images, const Json &request) {
    using rfc::numeric::ExactMatrix;
    members(request,
            {"schema", "operation", "query_resource", "reference_resource", "strands", "ambiguity",
             "maximum_mismatches", "junctions", "minimum_overhang", "evidence"});
    (void)text(request.at("evidence"), source_byte_limit());
    const auto &qimage = images.at(text(request.at("query_resource"))),
               &rimage = images.at(text(request.at("reference_resource")));
    require(rimage.declaration.at("format").text() == "fasta", "reference matching requires FASTA");
    const auto queries = sequence_records(qimage), references = sequence_records(rimage);
    sequence::reference_frame(rimage, references);
    require(!queries.empty(), "reference matching requires query records");
    std::size_t bases = 0;
    for (const auto &r : references)
        bases = rfc::machine::checked_size_add(bases, r.bases.size());
    const auto ambiguity = text(request.at("ambiguity"));
    require(ambiguity == "literal" || ambiguity == "possible" || ambiguity == "definite",
            "reference matching ambiguity rule");
    const auto max_mismatches =
                   count(request.at("maximum_mismatches"), std::numeric_limits<std::size_t>::max(),
                         "mismatch representation"),
               overhang =
                   count(request.at("minimum_overhang"), std::numeric_limits<std::size_t>::max(),
                         "splice overhang representation");
    require(overhang > 0, "splice overhang must be positive");
    const auto strands = sequence::strands(request);
    const auto compatible = [&](char q, char r) {
        if (ambiguity == "literal")
            return q == r;
        const auto a = iupac_mask(q), b = iupac_mask(r);
        return ambiguity == "possible" ? (a & b) != 0 : (a & ~b) == 0;
    };
    struct Junction {
        Json declaration;
        std::string contig, strand;
        std::size_t donor, acceptor;
    };
    std::vector<Junction> junctions;
    std::set<std::string> junction_ids;
    for (std::size_t i = 0; i < request.at("junctions").size(); ++i) {
        const auto j = request.at("junctions").at(i);
        members(j, {"id", "contig", "donor", "acceptor", "strand", "evidence"});
        require(junction_ids.insert(text(j.at("id"), source_byte_limit())).second,
                "duplicate splice identity");
        (void)text(j.at("evidence"), source_byte_limit());
        const auto contig = text(j.at("contig")), strand = text(j.at("strand"));
        require(strand == "+" || strand == "-", "declared splice strand");
        const auto donor = count(j.at("donor"), bases, "splice donor"),
                   acceptor = count(j.at("acceptor"), bases, "splice acceptor");
        const auto r = std::find_if(references.begin(), references.end(),
                                    [&](const auto &v) { return v.id == contig; });
        require(r != references.end() && donor < acceptor && acceptor < r->bases.size(),
                "declared splice interval outside retained contig");
        junctions.push_back({j, contig, strand, donor, acceptor});
    }
    struct Oriented {
        std::size_t query;
        std::string strand, bases;
        bool indexed;
        std::vector<std::pair<std::size_t, std::string>> seeds;
    };
    sequence::SeedIndex index;
    std::vector<Oriented> oriented;
    for (std::size_t q = 0; q < queries.size(); ++q) {
        require(!queries[q].bases.empty(), "reference query is empty");
        for (const auto &strand : strands) {
            Oriented v{q,
                       strand,
                       strand == "+" ? queries[q].bases
                                     : sequence::reverse_complement(queries[q].bases),
                       false,
                       {}};
            v.indexed =
                ambiguity == "literal" && max_mismatches < v.bases.size() &&
                std::all_of(
                    v.bases.begin(), v.bases.end(), [](char c) { return sequence::dna(c) >= 0; });
            if (v.indexed) {
                // At most k mismatches leave one of k+1 disjoint seeds intact.
                // Shortening a seed to the intrinsic uint64 capacity preserves
                // that necessary condition. Every candidate is checked in full.
                const auto parts = rfc::machine::checked_size_add(max_mismatches, 1);
                const auto length = std::min(sequence::token_bases, v.bases.size() / parts);
                for (std::size_t part = 0; part < parts; ++part) {
                    const auto offset =
                        rfc::machine::checked_size_mul(part, v.bases.size()) / parts;
                    const auto seed = v.bases.substr(offset, length);
                    index.request(seed);
                    v.seeds.emplace_back(offset, seed);
                }
            }
            oriented.push_back(std::move(v));
        }
    }
    index.build(references);
    using rfc::numeric::Integer;
    Integer comparisons{0}, candidates{0}, contiguous_universe{0}, junction_universe{0};
    const Rational mismatch_bound{std::to_string(max_mismatches)};
    auto matches = Json::array(), summaries = Json::array();
    std::size_t oriented_position = 0;
    for (const auto &q : queries) {
        std::size_t query_matches = 0;
        ExactMatrix count_map(1, q.bases.size());
        for (std::size_t i = 0; i < q.bases.size(); ++i)
            count_map.set(0, i, Rational{1});
        Vector disagreement(q.bases.size(), Rational{0});
        for (std::size_t orientation = 0; orientation < strands.size(); ++orientation) {
            const auto &v = oriented.at(oriented_position++);
            const auto &strand = v.strand, &sequence = v.bases;
            std::map<std::size_t, std::set<std::size_t>> starts;
            if (v.indexed)
                for (const auto &[offset, seed] : v.seeds)
                    for (const auto &hit : index.find(seed)) {
                        const auto length = references[hit.reference].bases.size();
                        if (hit.offset >= offset && sequence.size() <= length &&
                            hit.offset - offset <= length - sequence.size())
                            starts[hit.reference].insert(hit.offset - offset);
                    }
            const auto check = [&](const SequenceRecord &reference, std::size_t start,
                                   std::optional<std::size_t> split, const Junction *junction) {
                comparisons += Integer(std::to_string(q.bases.size()));
                ++candidates;
                for (std::size_t i = 0; i < sequence.size(); ++i) {
                    const auto position =
                        split && i >= *split ? junction->acceptor + i - *split : start + i;
                    disagreement[i] = compatible(sequence[i], reference.bases[position]) ? 0 : 1;
                }
                const auto mismatch = count_map.apply(disagreement).at(0);
                if (mismatch > mismatch_bound)
                    return;
                ++query_matches;
                auto m = Json::object(), blocks = Json::array();
                m.set("query_id", Json::string(q.id));
                m.set("reference_id", Json::string(reference.id));
                m.set("strand", Json::string(strand));
                m.set("mismatches", Json::string(mismatch.get_str()));
                m.set("query_evidence", q.evidence);
                m.set("reference_evidence", reference.evidence);
                const auto block = [&](std::size_t a, std::size_t b) {
                    auto p = Json::array();
                    p.push(cardinality(a));
                    p.push(cardinality(b));
                    blocks.push(p);
                };
                if (split) {
                    block(start, junction->donor);
                    block(junction->acceptor, junction->acceptor + sequence.size() - *split);
                    m.set("junction", junction->declaration);
                    m.set("query_split_in_reference_orientation", cardinality(*split));
                } else
                    block(start, start + sequence.size());
                m.set("blocks", blocks);
                matches.push(m);
            };
            for (std::size_t r_index = 0; r_index < references.size(); ++r_index) {
                const auto &r = references[r_index];
                if (sequence.size() <= r.bases.size()) {
                    contiguous_universe +=
                        Integer(std::to_string(r.bases.size() - sequence.size() + 1));
                    if (v.indexed) {
                        const auto found = starts.find(r_index);
                        if (found != starts.end())
                            for (auto start : found->second)
                                check(r, start, {}, nullptr);
                    } else
                        for (std::size_t start = 0; start <= r.bases.size() - sequence.size();
                             ++start)
                            check(r, start, {}, nullptr);
                }
                for (const auto &j : junctions) {
                    if (j.contig != r.id || j.strand != strand || overhang > sequence.size() / 2)
                        continue;
                    for (std::size_t split = overhang; split <= sequence.size() - overhang; ++split)
                        if (j.donor >= split &&
                            r.bases.size() - j.acceptor >= sequence.size() - split) {
                            ++junction_universe;
                            check(r, j.donor - split, split, &j);
                        }
                }
            }
        }
        auto s = Json::object();
        s.set("query_id", Json::string(q.id));
        s.set("matches", cardinality(query_matches));
        s.set("status", Json::string(query_matches ? "matched" : "no_match_in_declared_universe"));
        summaries.push(s);
    }
    auto out = Json::object();
    out.set("matches", matches);
    out.set("queries", summaries);
    out.set("tested_candidates", exact_count(candidates));
    out.set("base_comparisons", exact_count(comparisons));
    out.set("contiguous_candidate_universe", exact_count(contiguous_universe));
    out.set("junction_candidate_universe", exact_count(junction_universe));
    out.set("candidate_universe", exact_count(Integer(contiguous_universe + junction_universe)));
    out.set("index", index.receipt());
    out.set("compatibility_rule",
            Json::string(ambiguity == "literal"    ? "equal_tokens"
                         : ambiguity == "possible" ? "nonempty_IUPAC_intersection"
                                                   : "query_IUPAC_set_contained_in_reference"));
    out.set("complete", Json::boolean(true));
    out.set("coordinate_convention", Json::string("zero_based_half_open"));
    out.set(
        "search_universe",
        Json::string(
            "All contiguous placements and all single declared-junction placements, with specified "
            "strands, overhang, IUPAC rule and mismatch bound; no indels or undeclared junctions"));
    return out;
}
} // namespace gadget::assay_detail
