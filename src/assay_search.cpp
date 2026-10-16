#include "sequence_detail.hpp"
#include <algorithm>
#include <cctype>
#include <rfc/numeric/exact_matrix.hpp>

namespace gadget::assay_detail {
SequenceSource::SequenceSource(const ResourceImage &resource)
    : resource_(resource),
      spans_(read_sequence_spans(resource.bytes, resource.declaration.at("format").text())),
      digest_("sha256:" + sha256(resource.bytes)) {
    const std::string_view bytes(resource.bytes);
    segments_.reserve(spans_.size());
    for (const auto &span : spans_) {
        std::vector<Segment> segments;
        auto position = bytes.find('\n', span.start);
        require(position != std::string_view::npos, "sequence header has no bases");
        ++position;
        std::size_t begin = 0;
        while (begin < span.length) {
            auto end = bytes.find('\n', position);
            if (end == std::string_view::npos || end > span.end)
                end = span.end;
            auto stop = end;
            if (stop > position && bytes[stop - 1] == '\r')
                --stop;
            require(stop > position && stop - position <= span.length - begin,
                    "sequence byte span length");
            segments.push_back({position, begin, stop - position});
            begin = rfc::machine::checked_size_add(begin, stop - position);
            position = rfc::machine::checked_size_add(end, 1);
        }
        require(begin == span.length, "sequence byte span coverage");
        segments_.push_back(std::move(segments));
    }
}
std::string SequenceSource::id(std::size_t index) const {
    const auto &span = spans_.at(index);
    return resource_.bytes.substr(span.id_start, span.id_length);
}
Json SequenceSource::evidence(std::size_t index) const {
    const auto &span = spans_.at(index);
    auto result = Json::object();
    result.set("source_sha256", Json::string(digest_));
    result.set("record_index", cardinality(index));
    result.set("source_start", cardinality(span.start));
    result.set("source_end", cardinality(span.end));
    return result;
}
std::string SequenceSource::bases(std::size_t index, std::size_t begin, std::size_t length) const {
    require(begin <= this->length(index) && length <= this->length(index) - begin,
            "sequence window outside original record");
    std::string result;
    result.reserve(length);
    if (length == 0)
        return result;
    const auto &segments = segments_.at(index);
    auto segment =
        std::upper_bound(segments.begin(), segments.end(), begin,
                         [](auto position, const auto &part) { return position < part.begin; });
    --segment;
    while (result.size() < length) {
        const auto local = begin - segment->begin;
        const auto count = std::min(length - result.size(), segment->length - local);
        for (std::size_t i = 0; i < count; ++i)
            result += static_cast<char>(std::toupper(
                static_cast<unsigned char>(resource_.bytes[segment->source + local + i])));
        begin += count;
        ++segment;
    }
    return result;
}
SequenceRecord SequenceSource::record(std::size_t index) const {
    return {id(index), bases(index, 0, length(index)), evidence(index)};
}
std::vector<std::pair<std::size_t, std::size_t>> SequenceSource::byte_spans(
    std::size_t index, std::size_t begin, std::size_t length) const {
    require(begin <= this->length(index) && length <= this->length(index) - begin,
            "sequence byte window outside original record");
    std::vector<std::pair<std::size_t, std::size_t>> result;
    if (!length) return result;
    const auto &segments = segments_.at(index);
    auto segment = std::upper_bound(segments.begin(), segments.end(), begin,
        [](auto position, const auto &part) { return position < part.begin; });
    --segment;
    while (length) {
        const auto local = begin - segment->begin;
        const auto count = std::min(length, segment->length - local);
        result.emplace_back(segment->source + local, count);
        begin += count; length -= count; ++segment;
    }
    return result;
}
void SequenceSource::blocks(std::size_t index, std::size_t width,
                            const std::function<void(std::size_t, std::string_view)> &visit) const {
    require(width > 0, "sequence block width must be positive");
    std::string block;
    block.reserve(std::min(width, length(index)));
    std::size_t begin = 0;
    for (const auto &segment : segments_.at(index))
        for (std::size_t i = 0; i < segment.length; ++i) {
            block += static_cast<char>(
                std::toupper(static_cast<unsigned char>(resource_.bytes[segment.source + i])));
            if (block.size() == width) {
                visit(begin, block);
                begin += block.size();
                block.clear();
            }
        }
    if (!block.empty())
        visit(begin, block);
}
std::vector<SequenceRecord> sequence_records(const ResourceImage &resource) {
    const SequenceSource source(resource);
    std::vector<SequenceRecord> out;
    out.reserve(source.size());
    for (std::size_t i = 0; i < source.size(); ++i)
        out.push_back(source.record(i));
    return out;
}
AssayFieldsSource reference_match_fields_source(const ResourceImages &images, const Json &request) {
    using namespace sequence;
    using rfc::numeric::ExactMatrix;
    members(request,
            {"schema", "operation", "query_resource", "reference_resource", "strands", "ambiguity",
             "maximum_mismatches", "junctions", "minimum_overhang", "evidence"});
    (void)text(request.at("evidence"), source_byte_limit());
    const auto &qimage = images.at(text(request.at("query_resource"))),
               &rimage = images.at(text(request.at("reference_resource")));
    require(rimage.declaration.at("format").text() == "fasta", "reference matching requires FASTA");
    const SequenceSource query_source(qimage), reference_source(rimage);
    const auto queries = records(query_source), references = records(reference_source);
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
    struct Splice {
        const Junction &junction;
        std::size_t split;
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
    const auto dna_only = [](const auto &bases) {
        return std::all_of(bases.begin(), bases.end(),
                           [](char c) { return sequence::dna(c) >= 0; });
    };
    const bool dna_reference = ambiguity == "literal" &&
        std::all_of(references.begin(), references.end(),
                    [&](const auto &reference) { return dna_only(reference.bases); });
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
                (dna_reference || dna_only(v.bases));
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
                    // On a DNA reference, a literal outside ACGT necessarily mismatches.
                    // Such a seed cannot be the intact seed of an admitted placement.
                    if (!dna_only(seed)) continue;
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
    N::Array matches, summaries, children, contexts;
    std::vector<PlacementSource> placements; Vector collection_mismatches;
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(qimage.bytes))),
                                     RawByteSource(2, std::as_bytes(std::span(rimage.bytes)))};
    std::size_t oriented_position = 0;
    for (std::size_t query_index = 0; query_index < queries.size(); ++query_index) {
        const auto &q = queries[query_index];
        std::size_t query_matches = 0;
        ExactMatrix count_map(1, q.bases.size());
        for (std::size_t i = 0; i < q.bases.size(); ++i)
            count_map.set(0, i, Rational{1});
        rfc::core::RelationalComplexBuilder coordinate_builder(rfc::core::RelationSignature{});
        coordinate_builder.declare_grade({1});
        for (std::size_t i = 0; i < q.bases.size(); ++i) coordinate_builder.declare_witness_identity(i + 1);
        rfc::storage::TemporalEndpoint comparison_endpoint(coordinate_builder.freeze());
        comparison_endpoint.carriers.push_back({1, rfc::storage::TemporalCarrierKind::SourceGrade,
            {0}, {}, 0, rfc::numeric::IndexDomain::linear(q.bases.size())});
        SourceGradeFields fields(comparison_endpoint);
        const auto scalar = fields.explicit_carrier(1);
        comparison_endpoint.carriers.push_back({3, rfc::storage::TemporalCarrierKind::SourceGrade,
            {1}, {}, 0, rfc::numeric::IndexDomain::linear(0)});
        N::Object comparison_declaration;
        member(comparison_declaration, "query_record", integer(query_index));
        member(comparison_declaration, "scalar_carrier", integer(scalar));
        member(comparison_declaration, "scalar_coordinate", std::string("mismatch_count"));
        member(comparison_declaration, "count_map", fields.map("query.count_map", 1, scalar, count_map));
        N::Array native_comparisons;
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
                                   std::optional<Splice> splice) {
                comparisons += Integer(std::to_string(q.bases.size()));
                ++candidates;
                for (std::size_t i = 0; i < sequence.size(); ++i) {
                    const auto position =
                        splice && i >= splice->split
                            ? splice->junction.acceptor + i - splice->split : start + i;
                    disagreement[i] = compatible(sequence[i], reference.bases[position]) ? 0 : 1;
                }
                const auto mismatch = count_map.apply(disagreement).at(0);
                if (mismatch > mismatch_bound)
                    return;
                ++query_matches;
                auto m = N::Object{}; auto blocks = N::Array{};
                member(m, "query_id", N(q.id));
                member(m, "query_record_index", integer(query_index));
                member(m, "reference_id", N(reference.id));
                member(m, "strand", N(strand));
                member(m, "mismatches", N(mismatch));
                member(m, "query_evidence", import_record_view(q.evidence));
                member(m, "reference_evidence", import_record_view(reference.evidence));
                const auto block = [&](std::size_t a, std::size_t b) {
                    auto p = N::Array{};
                    element(p, integer(a));
                    element(p, integer(b));
                    element(blocks, p);
                };
                if (splice) {
                    const auto &junction = splice->junction;
                    block(start, junction.donor);
                    block(junction.acceptor, junction.acceptor + sequence.size() - splice->split);
                    member(m, "junction", import_record_view(junction.declaration));
                    member(m, "query_split_in_reference_orientation", integer(splice->split));
                } else
                    block(start, start + sequence.size());
                member(m, "blocks", blocks);
                const auto reference_index = static_cast<std::size_t>(&reference - references.data());
                const auto junction_index = splice ? std::optional<std::size_t>(
                    static_cast<std::size_t>(&splice->junction - junctions.data())) : std::nullopt;
                PlacementSource placement{query_index, reference_index, junction_index,
                    query_index, native_comparisons.elements.size(), {}};
                if (splice) {
                    const auto &junction = splice->junction;
                    placement.blocks.emplace_back(start, junction.donor);
                    placement.blocks.emplace_back(junction.acceptor, junction.acceptor + sequence.size() - splice->split);
                } else placement.blocks.emplace_back(start, start + sequence.size());
                placements.push_back(std::move(placement));
                collection_mismatches.push_back(mismatch);
                Vector original_disagreement(q.bases.size()), positions(q.bases.size()), permutation(q.bases.size());
                for (std::size_t i = 0; i < q.bases.size(); ++i) {
                    const auto original = strand == "+" ? i : q.bases.size() - 1 - i;
                    original_disagreement[original] = disagreement[i];
                    permutation[original] = Rational(std::to_string(i));
                    positions[original] = Rational(std::to_string(splice && i >= splice->split
                        ? splice->junction.acceptor + i - splice->split : start + i));
                }
                require(count_map.apply(original_disagreement).at(0) == mismatch,
                        "original coordinate mismatch-count action changed");
                N::Object operand;
                member(operand, "disagreement", fields.field("placement.disagreement", 0, original_disagreement));
                member(operand, "reference_positions", fields.field("placement.reference_positions", 0, positions));
                member(operand, "oriented_positions", fields.field("placement.oriented_positions", 0, permutation));
                member(operand, "reference_record", integer(reference_index));
                member(operand, "count", fields.vector("placement.mismatch_count", scalar, Vector{mismatch}));
                member(operand, "reverse", strand == "-");
                member(m, "comparison_index", integer(native_comparisons.elements.size()));
                element(native_comparisons, std::move(operand));
                element(matches, m);
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
                                check(r, start, std::nullopt);
                    } else
                        for (std::size_t start = 0; start <= r.bases.size() - sequence.size();
                             ++start)
                            check(r, start, std::nullopt);
                }
                for (const auto &j : junctions) {
                    if (j.contig != r.id || j.strand != strand || overhang > sequence.size() / 2)
                        continue;
                    for (std::size_t split = overhang; split <= sequence.size() - overhang; ++split)
                        if (j.donor >= split &&
                            r.bases.size() - j.acceptor >= sequence.size() - split) {
                            ++junction_universe;
                            check(r, j.donor - split, Splice{j, split});
                        }
                }
            }
        }
        member(comparison_declaration, "comparisons", std::move(native_comparisons));
        auto comparison_record = NativeRecord::create_from_endpoint(N(std::move(comparison_declaration)),
            new_native_record_namespace(), {}, std::move(comparison_endpoint));
        element(contexts, child_context(comparison_record));
        element(children, archive_child(comparison_record, sources));
        auto s = N::Object{};
        member(s, "query_id", N(q.id));
        member(s, "query_record_index", integer(query_index));
        member(s, "matches", integer(query_matches));
        member(s, "status", N(query_matches ? "matched" : "no_match_in_declared_universe"));
        element(summaries, s);
    }
    auto out = N::Object{};
    member(out, "matches", matches);
    member(out, "queries", summaries);
    member(out, "tested_candidates", N(candidates));
    member(out, "base_comparisons", N(comparisons));
    member(out, "contiguous_candidate_universe", N(contiguous_universe));
    member(out, "junction_candidate_universe", N(junction_universe));
    member(out, "candidate_universe", N(Integer(contiguous_universe + junction_universe)));
    member(out, "index", import_record_view(index.receipt()));
    member(out, "compatibility_rule",
            N(ambiguity == "literal"    ? "equal_tokens"
                         : ambiguity == "possible" ? "nonempty_IUPAC_intersection"
                                                   : "query_IUPAC_set_contained_in_reference"));
    member(out, "complete", N(true));
    member(out, "coordinate_convention", N("zero_based_half_open"));
    member(out,
        "search_universe",
        N(
            "All contiguous placements and all single declared-junction placements, with specified "
            "strands, overhang, IUPAC rule and mismatch bound; no indels or undeclared junctions"));
    auto result = collection_source("gadget.sequence.reference.fields.v1", request, qimage, rimage,
        query_source, reference_source, placements, std::move(out), std::move(children), std::move(contexts), std::move(sources));
    SourceGradeFields collection_fields(result.endpoint);
    auto &declaration = std::get<N::Object>(result.declaration.data);
    member(declaration, "mismatches", collection_fields.field("placement.mismatches", 1, collection_mismatches));
    return result;
}
Json reference_match_fields_view(const NativeRecord &record) {
    auto out = render_record_view(record.at("view"));
    const auto &values = retained_tensor(record, record.at("mismatches")).tensor;
    std::vector<NativeRecord> children;
    children.reserve(record.at("children").size());
    for (std::size_t i = 0; i < record.at("children").size(); ++i)
        children.push_back(NativeRecord::decode(record.at("children").at(i).bytes()));
    for (std::size_t i = 0; i < out.at("matches").size(); ++i) {
        auto match = out.at("matches").at(i);
        const auto coupling = record.at("placements").at(i);
        const auto child = children.at(coupling.at("query_record").integer().get_ui());
        const auto comparison = child.at("comparisons").at(coupling.at("operation_index").integer().get_ui());
        const auto count = retained_tensor(child, comparison.at("count")).tensor.state(0).value();
        require(count == values.state(i).value(), "placement count disagrees with its retained comparison");
        json_object_object_del(match.raw(), "comparison_index");
        match.set("mismatches", Json::string(count.get_str()));
    }
    return out;
}
Json reference_matches(const ResourceImages &images, const Json &request) {
    return reference_match_fields_view(sequence::source_record(reference_match_fields_source(images, request)));
}
} // namespace gadget::assay_detail
