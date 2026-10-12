#include "sequence_detail.hpp"
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <sstream>

namespace gadget::assay_detail {
NativeRecordValue read_protocol_fields(ResourceImages images, const Json &request) {
    using N = NativeRecordValue;
    using Bits = rfc::machine::PackedBits;
    const auto integer = [](std::size_t value) {
        return N(rfc::numeric::Integer(std::to_string(value)));
    };
    members(request, {"schema", "operation", "resource", "resource_sha256", "protocol"});
    require(request.at("schema").text() == "gadget.assays.query.v1" &&
                request.at("operation").text() == "read_protocol_evidence",
            "read protocol query schema/operation");
    const auto owners = std::make_shared<const ResourceImages>(std::move(images));
    const auto resource_id = text(request.at("resource"), source_byte_limit());
    require(owners->contains(resource_id), "read protocol resource absent");
    const auto reads = std::make_shared<const SequenceSource>(owners->at(resource_id));
    require(request.at("resource_sha256").text() == reads->digest(),
            "read protocol original source hash mismatch");
    const auto protocol = request.at("protocol");
    members(protocol, {"id", "evidence", "count_unit", "segments"}, {"barcode"});
    (void)text(protocol.at("id"), source_byte_limit());
    (void)text(protocol.at("evidence"), source_byte_limit());
    require(protocol.at("count_unit").text() == "read_record_occurrence",
            "read protocol count unit requires original read record occurrences");
    struct Segment {
        std::string id, role, strand;
        std::size_t start, end;
    };
    std::vector<Segment> segments;
    std::set<std::string> ids;
    require(protocol.at("segments").size() > 0, "read protocol requires segments");
    for (std::size_t i = 0; i < protocol.at("segments").size(); ++i) {
        const auto declaration = protocol.at("segments").at(i);
        members(declaration, {"id", "role", "start", "end", "strand"});
        Segment segment{text(declaration.at("id"), source_byte_limit()),
                        text(declaration.at("role")), text(declaration.at("strand")),
                        count(declaration.at("start"), std::numeric_limits<std::size_t>::max(),
                              "read protocol start representation"),
                        count(declaration.at("end"), std::numeric_limits<std::size_t>::max(),
                              "read protocol end representation")};
        require(ids.insert(segment.id).second, "duplicate read protocol segment identity");
        require(segment.role == "barcode" || segment.role == "umi" || segment.role == "sequence",
                "read protocol segment role");
        require(segment.start < segment.end, "read protocol segment must be nonempty");
        require(segment.strand == "+" || segment.strand == "-", "read protocol segment strand");
        segments.push_back(std::move(segment));
    }
    std::optional<std::size_t> barcode;
    std::size_t budget = 0;
    struct WhitelistEntry {
        std::string id, bases;
        std::size_t start, end;
    };
    std::vector<WhitelistEntry> whitelist;
    std::string whitelist_digest;
    if (protocol.has("barcode")) {
        const auto declaration = protocol.at("barcode");
        members(declaration, {"segment", "whitelist_resource", "whitelist_sha256", "objective",
                              "maximum_substitutions", "symbol_rule"});
        require(declaration.at("objective").text() == "minimum_substitutions" &&
                    declaration.at("symbol_rule").text() == "literal_iupac_tokens",
                "read barcode requires fixed length literal token substitution objective");
        budget = count(declaration.at("maximum_substitutions"),
                       std::numeric_limits<std::size_t>::max(),
                       "read barcode substitution budget representation");
        for (std::size_t i = 0; i < segments.size(); ++i)
            if (segments[i].id == declaration.at("segment").text()) barcode = i;
        require(barcode.has_value() && segments[*barcode].role == "barcode",
                "read barcode segment is not declared with barcode role");
        const auto id = text(declaration.at("whitelist_resource"), source_byte_limit());
        require(owners->contains(id) && owners->at(id).declaration.at("format").text() == "fasta",
                "read barcode whitelist requires an original FASTA resource");
        const SequenceSource source(owners->at(id));
        require(source.digest() == declaration.at("whitelist_sha256").text(),
                "read barcode original whitelist hash mismatch");
        require(source.size() > 0, "read barcode whitelist is empty");
        whitelist_digest = source.digest();
        for (std::size_t i = 0; i < source.size(); ++i) {
            require(source.length(i) == segments[*barcode].end - segments[*barcode].start,
                    "read barcode whitelist length differs from the declared segment");
            const auto &span = source.original_span(i);
            whitelist.push_back({source.id(i), source.bases(i, 0, source.length(i)),
                                 span.start, span.end});
        }
    }
    N::PackedArray records;
    records.count = reads->size();
    records.emit = [owners, reads, segments = std::move(segments), barcode, budget,
                    whitelist = std::move(whitelist), whitelist_digest, integer](const auto &emit) {
        (void)owners; // Own the original byte storage referenced by SequenceSource.
        // Keep every source occurrence. Token matches do not identify molecules.
        for (std::size_t i = 0; i < reads->size(); ++i) {
            N::Array extracted, candidates;
            std::optional<std::string> barcode_tokens;
            bool complete = true;
            for (std::size_t j = 0; j < segments.size(); ++j) {
                const auto &segment = segments[j];
                const bool present = segment.end <= reads->length(i);
                complete &= present;
                N tokens;
                if (present) {
                    auto bases = reads->bases(i, segment.start, segment.end - segment.start);
                    if (segment.strand == "-") bases = sequence::reverse_complement(bases);
                    if (barcode && *barcode == j) barcode_tokens = bases;
                    tokens = N(std::move(bases));
                }
                extracted.elements.emplace_back(N::Object{{
                    {"id", N(segment.id)}, {"role", N(segment.role)},
                    {"start", integer(segment.start)}, {"end", integer(segment.end)},
                    {"strand", N(segment.strand)}, {"tokens", std::move(tokens)},
                    {"status", N(std::string(present ? "present" : "outside_original_read"))}
                }});
            }
            std::optional<std::size_t> minimum;
            if (barcode_tokens) {
                std::vector<std::size_t> costs;
                costs.reserve(whitelist.size());
                Bits mismatches(barcode_tokens->size(), false,
                                Bits::words_for_size(barcode_tokens->size()));
                for (const auto &entry : whitelist) {
                    for (std::size_t p = 0; p < barcode_tokens->size(); ++p)
                        mismatches.set(p, (*barcode_tokens)[p] != entry.bases[p]);
                    const auto cost = mismatches.count();
                    costs.push_back(cost);
                    minimum = minimum ? std::min(*minimum, cost) : cost;
                }
                if (*minimum <= budget)
                    for (std::size_t w = 0; w < whitelist.size(); ++w)
                        if (costs[w] == *minimum)
                            candidates.elements.emplace_back(N::Object{{
                                {"whitelist_record_index", integer(w)},
                                {"id", N(whitelist[w].id)}, {"tokens", N(whitelist[w].bases)},
                                {"substitutions", integer(costs[w])},
                                {"source_sha256", N(whitelist_digest)},
                                {"source_start", integer(whitelist[w].start)},
                                {"source_end", integer(whitelist[w].end)}
                            }});
            }
            const auto &source = reads->original_span(i);
            emit(N(N::Object{{
                {"record_index", integer(i)}, {"read_id", N(reads->id(i))},
                {"source_start", integer(source.start)}, {"source_end", integer(source.end)},
                {"read_length", integer(reads->length(i))},
                {"status", N(std::string(complete ? "complete" : "segment_outside_original_read"))},
                {"segments", N(std::move(extracted))},
                {"minimum_substitutions", minimum ? integer(*minimum) : N{}},
                {"barcode_status", N(std::string(!barcode ? "not_requested" :
                    !barcode_tokens ? "segment_outside_original_read" :
                    candidates.elements.empty() ? "outside_substitution_budget" :
                    candidates.elements.size() == 1 ? "unique_minimum" : "ambiguous_minimum"))},
                {"barcode_candidates", N(std::move(candidates))}
            }}));
        }
    };
    return N(N::Object{{
        {"schema", N(std::string("gadget.read.protocol.evidence.v1"))},
        {"resource", N(resource_id)}, {"source_sha256", N(reads->digest())},
        {"protocol", import_record_view(protocol)}, {"record_count", integer(reads->size())},
        {"records", N(std::move(records))},
        {"coordinate_convention", N(std::string("zero_based_half_open_in_original_read"))},
        {"measurement", N(std::string("Original read occurrences and declared literal segments; "
            "all minimum fixed length whitelist token substitution candidates within the declared "
            "budget. UMI tokens remain observed evidence; no molecule deduplication is performed."))}
    }});
}
Json salmon_classes(const std::string &bytes, Vector *native_counts) {
    const auto ls = lines(bytes);
    require(ls.size() >= 3, "truncated Salmon equivalence-class file");
    const auto line = [&](std::size_t i) {
        return bytes.substr(ls.at(i).first, ls.at(i).second - ls.at(i).first);
    };
    const auto n = rfc::machine::checked_size(unsigned_token(line(0))),
               m = rfc::machine::checked_size(unsigned_token(line(1)));
    require(n > 0 && m > 0 &&
                ls.size() ==
                    rfc::machine::checked_size_add(2, rfc::machine::checked_size_add(n, m)),
            "Salmon class header/dimensions");
    auto candidates = Json::array(), classes = Json::array(), counts = Json::array();
    if (native_counts) native_counts->clear();
    std::vector<std::string> ids;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < n; ++i) {
        const auto id = line(i + 2);
        require(!id.empty() && seen.insert(id).second, "Salmon transcript identity");
        ids.push_back(id);
        candidates.push(Json::string(id));
    }
    for (std::size_t i = 0; i < m; ++i) {
        std::istringstream input(line(2 + n + i));
        std::vector<std::string> fields;
        std::string s;
        while (input >> s)
            fields.push_back(s);
        require(!fields.empty(), "empty Salmon class");
        const auto k = rfc::machine::checked_size(unsigned_token(fields[0]));
        require(k > 0 && fields.size() == rfc::machine::checked_size_add(k, 2),
                "Salmon profile admits --dumpEq classes; conditional-weight classes require a "
                "separate measurement declaration");
        auto c = Json::object(), p = Json::array();
        std::set<std::size_t> members;
        for (std::size_t j = 0; j < k; ++j) {
            const auto index = unsigned_token(fields[j + 1]);
            require(index < n && members.insert(index).second,
                    "Salmon candidate index range/duplicate");
            p.push(Json::string(ids[index]));
        }
        const auto count = exact_rational(fields.back());
        require(count >= 0 && count.get_den() == 1, "Salmon class fragment count");
        c.set("id", Json::string("class:" + std::to_string(i)));
        c.set("candidates", p);
        classes.push(c);
        counts.push(Json::string(count.get_str()));
        if (native_counts) native_counts->push_back(count);
    }
    auto result = Json::object();
    result.set("candidates", candidates);
    result.set("classes", classes);
    result.set("counts", counts);
    result.set("record_count", cardinality(m));
    result.set("measurement", Json::string("Reported unweighted equivalence-class fragment counts; "
                                           "transcript candidates retain their original order"));
    return result;
}
unsigned iupac_mask(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
    case 'A':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 4;
    case 'T':
    case 'U':
        return 8;
    case 'R':
        return 5;
    case 'Y':
        return 10;
    case 'S':
        return 6;
    case 'W':
        return 9;
    case 'K':
        return 12;
    case 'M':
        return 3;
    case 'B':
        return 14;
    case 'D':
        return 13;
    case 'H':
        return 11;
    case 'V':
        return 7;
    case 'N':
        return 15;
    default:
        throw std::invalid_argument("invalid IUPAC motif symbol");
    }
}
namespace {
unsigned complement(unsigned m) {
    return ((m & 1) << 3) | ((m & 2) << 1) | ((m & 4) >> 1) | ((m & 8) >> 3);
}
} // namespace
Json motif_analysis(const std::string &bytes, const Json &resource, const Json &request) {
    require(resource.at("format").text() == "fasta", "motif occurrence analysis requires FASTA");
    const auto inspected = inspect_resource(bytes, resource);
    const auto patterns = request.at("motifs");
    const auto strands = text(request.at("strands"));
    require(strands == "forward" || strands == "both", "motif strand policy");
    require(patterns.size() > 0, "motifs must be nonempty");
    std::vector<std::vector<unsigned>> motifs;
    std::set<std::string> ids;
    for (std::size_t i = 0; i < patterns.size(); ++i) {
        const auto p = patterns.at(i);
        members(p, {"id", "pattern"});
        require(ids.insert(text(p.at("id"), source_byte_limit())).second, "duplicate motif identity");
        const auto s = text(p.at("pattern"), source_byte_limit());
        std::vector<unsigned> v;
        for (auto c : s)
            v.push_back(iupac_mask(c));
        motifs.push_back(v);
    }
    auto hits = Json::array(), summaries = Json::array();
    for (std::size_t r = 0; r < inspected.at("records").size(); ++r) {
        const auto record = inspected.at("records").at(r);
        const auto p = count(record.at("source_start"), bytes.size(), "FASTA offset"),
                   q = count(record.at("source_end"), bytes.size(), "FASTA offset");
        const auto part = bytes.substr(p, q - p);
        const auto ls = lines(part);
        std::vector<unsigned> sequence;
        for (std::size_t i = 1; i < ls.size(); ++i)
            for (auto ch : part.substr(ls[i].first, ls[i].second - ls[i].first))
                sequence.push_back(iupac_mask(ch));
        for (std::size_t t = 0; t < motifs.size(); ++t) {
            std::size_t definite = 0, possible = 0;
            for (int strand = 0; strand < (strands == "both" ? 2 : 1); ++strand) {
                auto motif = motifs[t];
                if (strand) {
                    std::reverse(motif.begin(), motif.end());
                    for (auto &m : motif)
                        m = complement(m);
                }
                if (motif.size() > sequence.size())
                    continue;
                for (std::size_t i = 0; i <= sequence.size() - motif.size(); ++i) {
                    bool certain = true, compatible = true;
                    for (std::size_t j = 0; j < motif.size(); ++j) {
                        certain &= (sequence[i + j] & ~motif[j]) == 0;
                        compatible &= (sequence[i + j] & motif[j]) != 0;
                    }
                    if (!compatible)
                        continue;
                    possible = rfc::machine::checked_size_add(possible, 1);
                    if (certain)
                        definite = rfc::machine::checked_size_add(definite, 1);
                    auto hit = Json::object();
                    hit.set("sequence", record.at("id"));
                    hit.set("motif", patterns.at(t).at("id"));
                    hit.set("start", cardinality(i));
                    hit.set("end", cardinality(i + motif.size()));
                    hit.set("strand", Json::string(strand ? "-" : "+"));
                    hit.set("definite", Json::boolean(certain));
                    hits.push(hit);
                }
            }
            auto summary = Json::object();
            summary.set("sequence", record.at("id"));
            summary.set("motif", patterns.at(t).at("id"));
            summary.set("definite_occurrences", cardinality(definite));
            summary.set("possible_occurrences", cardinality(possible));
            summaries.push(summary);
        }
    }
    auto result = Json::object();
    result.set("hits", hits);
    result.set("summaries", summaries);
    result.set("coordinate_convention", Json::string("zero_based_half_open_within_sequence"));
    result.set(
        "ambiguity",
        Json::string("Definite means every retained IUPAC base choice matches; possible means at "
                     "least one choice matches that window. Possible counts are marginal window "
                     "bounds, not joint extrema for overlapping ambiguous windows."));
    result.set("measurement", Json::string("Sequence occurrences under declared IUPAC patterns; no "
                                           "background enrichment or regulatory interpretation"));
    return result;
}
} // namespace gadget::assay_detail
