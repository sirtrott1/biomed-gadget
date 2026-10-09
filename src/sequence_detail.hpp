#pragma once
#include "genomic_detail.hpp"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <map>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/path/objectives.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <set>
#include <tuple>
#include <unordered_map>

namespace gadget::assay_detail::sequence {
// Source indexing and hypothesis declarations belong to the assay adapter.
// Exact path optimization, field realization and gluing remain native APIs.
inline std::string reverse_complement(const std::string &s) {
    const std::string alphabet = "ACGTRYSWKMBDHVN", complement = "TGCAYRSWMKVHDBN";
    std::string out;
    out.reserve(s.size());
    for (auto i = s.rbegin(); i != s.rend(); ++i) {
        const auto p = alphabet.find(*i);
        require(p != std::string::npos, "reverse search requires DNA IUPAC tokens");
        out += complement[p];
    }
    return out;
}
inline std::vector<std::string> strands(const Json &request) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < request.at("strands").size(); ++i) {
        const auto s = text(request.at("strands").at(i));
        require((s == "+" || s == "-") && std::find(out.begin(), out.end(), s) == out.end(),
                "search requires distinct DNA strands");
        out.push_back(s);
    }
    require(!out.empty(), "search requires a strand");
    return out;
}
inline void reference_frame(const ResourceImage &image, const std::vector<SequenceRecord> &refs) {
    require(image.declaration.at("format").text() == "fasta", "reference requires FASTA");
    if (!image.declaration.has("reference"))
        return;
    const auto frame = reference_record(image.declaration.at("reference"));
    require(frame.at("sequence_sha256").text() == "sha256:" + sha256(image.bytes) &&
                frame.at("contigs").size() == refs.size(),
            "reference hash/dictionary mismatch");
    std::map<std::string, std::size_t> lengths;
    for (const auto &r : refs)
        lengths.emplace(r.id, r.bases.size());
    for (std::size_t i = 0; i < frame.at("contigs").size(); ++i) {
        const auto c = frame.at("contigs").at(i);
        require(lengths.contains(c.at("id").text()) &&
                    lengths.at(c.at("id").text()) ==
                        static_cast<std::size_t>(c.at("length").integer()),
                "reference contig length mismatch");
    }
}
inline int dna(char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return -1;
    }
}
inline constexpr std::size_t token_bases = std::numeric_limits<std::uint64_t>::digits / 2;
inline std::uint64_t token(const std::string &s) {
    require(!s.empty() && s.size() <= token_bases, "DNA index token representation");
    std::uint64_t value = 0;
    for (const auto c : s) {
        const auto d = dna(c);
        require(d >= 0,
                "indexed query requires A, C, G or T; ambiguity must be resolved explicitly");
        value = (value << 2u) | static_cast<unsigned>(d);
    }
    return value;
}
struct Hit {
    std::size_t reference, offset;
};
// Injective two-bit keys, not a probabilistic hash or a sampled reference index.
class SeedIndex {
    using Table = std::unordered_map<std::uint64_t, std::vector<Hit>>;
    std::map<std::size_t, Table> tables_;
    std::size_t hits_{}, scans_{};

  public:
    void request(const std::string &seed) { tables_[seed.size()].try_emplace(token(seed)); }
    void build(const std::vector<SequenceRecord> &refs) {
        for (auto &[length, table] : tables_) {
            table.reserve(table.size());
            const auto mask = length == token_bases ? std::numeric_limits<std::uint64_t>::max()
                                                    : (std::uint64_t{1} << (2u * length)) - 1u;
            for (std::size_t r = 0; r < refs.size(); ++r) {
                std::uint64_t code = 0;
                std::size_t run = 0;
                for (std::size_t i = 0; i < refs[r].bases.size(); ++i) {
                    scans_ = rfc::machine::checked_size_add(scans_, 1);
                    const auto d = dna(refs[r].bases[i]);
                    if (d < 0) {
                        code = 0;
                        run = 0;
                        continue;
                    }
                    code = ((code << 2u) | static_cast<unsigned>(d)) & mask;
                    if (++run < length)
                        continue;
                    const auto it = table.find(code);
                    if (it != table.end()) {
                        it->second.push_back({r, i + 1 - length});
                        hits_ = rfc::machine::checked_size_add(hits_, 1);
                    }
                }
            }
        }
    }
    const std::vector<Hit> &find(const std::string &seed) const {
        return tables_.at(seed.size()).at(token(seed));
    }
    Json receipt() const {
        auto out = Json::object();
        out.set("reference_base_visits", cardinality(scans_));
        out.set("retained_seed_hits", cardinality(hits_));
        out.set("seed_lengths", cardinality(tables_.size()));
        out.set("key_encoding", Json::string("injective_two_bit_DNA"));
        return out;
    }
};
struct Step {
    std::size_t from, to;
    Rational cost;
};
class Routes {
    rfc::core::RelationSignature signature_;
    rfc::core::RelationTypeId type_;
    rfc::core::RelationalComplexBuilder builder_;
    std::size_t vertices_;
    mutable std::optional<Vector> prepared_costs_;

  public:
    std::vector<Step> steps;
    explicit Routes(std::size_t vertices, const std::string &type)
        : type_(signature_.add_type(type)), builder_(signature_), vertices_(vertices) {
        builder_.declare_grade(rfc::core::Grade{1});
        require(vertices <= std::numeric_limits<std::uint32_t>::max(),
                "declared route native vertex coordinate representation");
        for (std::size_t i = 0; i < vertices; ++i)
            builder_.declare_witness_identity(i + 1);
    }
    std::size_t add(std::size_t from, std::size_t to, const Rational &cost) {
        prepared_costs_.reset();
        const auto index = steps.size();
        require(index < std::numeric_limits<std::uint32_t>::max() && from < vertices_ &&
                    to < vertices_ && cost >= 0,
                "declared route native relation coordinates/cost");
        const auto occurrence = rfc::machine::checked_size_mul(index, 2);
        steps.push_back({from, to, cost});
        builder_.add_grade_one(
            rfc::core::RelationKey{rfc::machine::checked_size_add(1000000, index)}, type_,
            {{from + 1,
              true,
              true,
              Rational{1},
              {},
              rfc::core::OccurrenceKey{rfc::machine::checked_size_add(occurrence, 1)}},
             {to + 1,
              true,
              false,
              Rational{1},
              {},
              rfc::core::OccurrenceKey{rfc::machine::checked_size_add(occurrence, 2)}}});
        return index;
    }
    rfc::core::RelationalComplex freeze() { return builder_.freeze(); }
    const Vector &costs() const {
        if (!prepared_costs_) {
            Vector costs;
            costs.reserve(steps.size());
            for (const auto &step : steps)
                costs.push_back(step.cost);
            prepared_costs_ = std::move(costs);
        }
        return *prepared_costs_;
    }
};
inline Json realized_receipt(const rfc::core::RelationalComplex &source, const Routes &routes,
                             const std::vector<std::uint32_t> &relations, std::size_t from,
                             std::size_t to, const Rational &cost,
                             rfc::path::RealizedPath *retained = nullptr) {
    std::vector<rfc::path::TraversalStep> traversal;
    for (auto e : relations)
        traversal.push_back({e, Rational{1}});
    auto p = rfc::path::realize(source, rfc::core::Grade{1}, std::move(traversal));
    Vector demand(source.grade(rfc::core::Grade{0}).cell_count(), 0);
    require(source.grade(rfc::core::Grade{0}).identities.at(from) == from + 1 &&
                source.grade(rfc::core::Grade{0}).identities.at(to) == to + 1,
            "native route coordinate declaration");
    demand[from] -= 1;
    demand[to] += 1;
    require(p.boundary == demand && rfc::numeric::dot(routes.costs(), p.field) == cost,
            "native ordered route realization/cost mismatch");
    auto out = Json::object(), ordered = Json::array(), field = Json::array();
    for (auto e : relations)
        ordered.push(cardinality(e));
    for (std::size_t e = 0; e < p.field.size(); ++e)
        if (p.field[e] != 0) {
            auto v = Json::object();
            v.set("relation", cardinality(e));
            v.set("coefficient", Json::string(p.field[e].get_str()));
            field.push(v);
        }
    out.set("ordered_relations", ordered);
    out.set("realized_field", field);
    out.set("boundary_verified", Json::boolean(true));
    out.set("cost_verified", Json::boolean(true));
    if (retained)
        *retained = std::move(p);
    return out;
}
inline Json carrier(const rfc::core::RelationalComplex &source, const Routes &routes) {
    auto out = Json::object(), steps = Json::array();
    const auto bytes = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(source), 1, source_codec_limits(source_byte_limit())));
    for (const auto &s : routes.steps) {
        auto row = Json::array();
        row.push(cardinality(s.from));
        row.push(cardinality(s.to));
        row.push(Json::string(s.cost.get_str()));
        steps.push(row);
    }
    out.set("steps", steps);
    out.set("native_complex_sha256",
            Json::string(
                "sha256:" +
                sha256(std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    out.set("vertices", cardinality(source.grade(rfc::core::Grade{0}).cell_count()));
    out.set("coordinate_rule", Json::string("vertex v has native identity v+1; step e has relation "
                                            "identity 1000000+e and occurrences 2e+1,2e+2"));
    return out;
}
} // namespace gadget::assay_detail::sequence
