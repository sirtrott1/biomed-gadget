#include "assay_detail.hpp"
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <sstream>

namespace gadget::assay_detail {
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
