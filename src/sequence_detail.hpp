#pragma once
#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
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
inline void reference_frame(const ResourceImage &image, const SequenceSource &refs) {
    require(image.declaration.at("format").text() == "fasta", "reference requires FASTA");
    if (!image.declaration.has("reference"))
        return;
    const auto frame = reference_record(image.declaration.at("reference"));
    require(frame.at("sequence_sha256").text() == refs.digest() &&
                frame.at("contigs").size() == refs.size(),
            "reference hash/dictionary mismatch");
    std::map<std::string, std::size_t> lengths;
    for (std::size_t i = 0; i < refs.size(); ++i)
        lengths.emplace(refs.id(i), refs.length(i));
    for (std::size_t i = 0; i < frame.at("contigs").size(); ++i) {
        const auto contig = frame.at("contigs").at(i);
        require(lengths.contains(contig.at("id").text()) &&
                    lengths.at(contig.at("id").text()) ==
                        static_cast<std::size_t>(contig.at("length").integer()),
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
// Injective two bit keys, not a probabilistic hash or a sampled reference index.
class SeedIndex {
    using Table = std::unordered_map<std::uint64_t, std::vector<Hit>>;
    std::map<std::size_t, Table> tables_;
    std::size_t hits_{}, scans_{}, blocks_{}, peak_block_{};

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
    void build(const SequenceSource &refs, std::size_t block_bases) {
        require(block_bases > 0, "reference block width must be positive");
        for (auto &[length, table] : tables_) {
            table.reserve(table.size());
            const auto mask = length == token_bases ? std::numeric_limits<std::uint64_t>::max()
                                                    : (std::uint64_t{1} << (2u * length)) - 1u;
            for (std::size_t r = 0; r < refs.size(); ++r) {
                std::uint64_t code = 0;
                std::size_t run = 0;
                refs.blocks(r, block_bases, [&](std::size_t begin, std::string_view bases) {
                    blocks_ = rfc::machine::checked_size_add(blocks_, 1);
                    peak_block_ = std::max(peak_block_, bases.size());
                    for (std::size_t i = 0; i < bases.size(); ++i) {
                        scans_ = rfc::machine::checked_size_add(scans_, 1);
                        const auto d = dna(bases[i]);
                        if (d < 0) {
                            code = 0;
                            run = 0;
                            continue;
                        }
                        code = ((code << 2u) | static_cast<unsigned>(d)) & mask;
                        if (++run < length)
                            continue;
                        const auto found = table.find(code);
                        if (found != table.end()) {
                            found->second.push_back({r, begin + i + 1 - length});
                            hits_ = rfc::machine::checked_size_add(hits_, 1);
                        }
                    }
                });
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
        out.set("reference_blocks", cardinality(blocks_));
        out.set("peak_normalized_reference_block_bases", cardinality(peak_block_));
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
using N = NativeRecordValue;
inline std::vector<SequenceRecord> records(const SequenceSource &source) {
    std::vector<SequenceRecord> result;
    result.reserve(source.size());
    for (std::size_t i = 0; i < source.size(); ++i) result.push_back(source.record(i));
    return result;
}
inline N integer(std::size_t value) { return N(rfc::numeric::Integer(std::to_string(value))); }
template<class T> inline void member(N::Object &object, std::string key, T &&value) {
    object.members.emplace_back(std::move(key), N(std::forward<T>(value)));
}
template<class T> inline void element(N::Array &array, T &&value) {
    array.elements.emplace_back(std::forward<T>(value));
}
inline rfc::storage::TemporalEndpoint source_endpoint(const rfc::core::RelationalComplex &source) {
    rfc::storage::TemporalEndpoint endpoint(source);
    for (std::size_t k = 0; k < source.grade_count(); ++k)
        endpoint.carriers.push_back({k + 1, rfc::storage::TemporalCarrierKind::SourceGrade,
            rfc::core::Grade{static_cast<std::uint16_t>(k)}, {}, 0,
            rfc::numeric::IndexDomain::linear(source.grade({static_cast<std::uint16_t>(k)}).cell_count())});
    return endpoint;
}
inline AssayFieldsSource route_source(const rfc::core::RelationalComplex &source, const Routes &routes) {
    auto endpoint = source_endpoint(source);
    SourceGradeFields fields(endpoint);
    N::Object declaration;
    member(declaration, "costs", fields.field("route.cost", 1, routes.costs()));
    return {std::move(endpoint), N(std::move(declaration)), {}};
}
inline N retain_route(SourceGradeFields &fields, const rfc::core::RelationalComplex &source,
                      const Routes &routes, const std::vector<std::uint32_t> &relations,
                      std::size_t from, std::size_t to, const Rational &cost,
                      rfc::path::RealizedPath *retained = nullptr) {
    std::vector<rfc::path::TraversalStep> traversal;
    traversal.reserve(relations.size());
    for (const auto e : relations) traversal.push_back({e, Rational{1}});
    auto path = rfc::path::realize(source, rfc::core::Grade{1}, std::move(traversal));
    Vector demand(source.grade({0}).cell_count(), Rational{0});
    demand.at(from) -= 1; demand.at(to) += 1;
    require(path.boundary == demand && rfc::numeric::dot(routes.costs(), path.field) == cost,
            "native ordered route realization/cost mismatch");
    N::Array ordered;
    ordered.elements.reserve(relations.size());
    for (const auto e : relations) element(ordered, integer(e));
    N::Object result;
    member(result, "ordered_relations", std::move(ordered));
    member(result, "realized_field", fields.field("route.realization", 1, path.field));
    member(result, "boundary", fields.field("route.boundary", 0, path.boundary));
    member(result, "from_vertex", integer(from)); member(result, "to_vertex", integer(to));
    member(result, "minimum_cost", cost);
    if (retained) *retained = std::move(path);
    return N(std::move(result));
}
inline Json retained_route_view(const NativeRecord &record) {
    auto result = Json::object(), values = Json::array();
    result.set("ordered_relations", render_record_view(record.at("ordered_relations")));
    const auto &field = retained_tensor(record, record.at("realized_field")).tensor;
    for (std::size_t e = 0; e < field.domain().logical_size(); ++e) {
        const auto value = field.state(e);
        if (value && *value != 0) {
            auto row = Json::object(); row.set("relation", cardinality(e));
            row.set("coefficient", Json::string(value->get_str())); values.push(row);
        }
    }
    result.set("realized_field", values);
    result.set("boundary_verified", Json::boolean(true));
    result.set("cost_verified", Json::boolean(true));
    return result;
}
inline Json carrier_fields_view(const NativeRecord &record, bool extras = true) {
    const auto &source = record.endpoint()->source;
    const auto &arena = source.grade({1});
    const auto &costs = retained_tensor(record, record.at("costs")).tensor;
    auto result = Json::object(), steps = Json::array();
    for (std::size_t e = 0; e < arena.cell_count(); ++e) {
        const auto range = arena.incidence_range(e);
        require(range.second - range.first == 2, "native route occurrence arity");
        auto step = Json::array();
        step.push(cardinality(arena.lower_index.at(range.first)));
        step.push(cardinality(arena.lower_index.at(range.first + 1)));
        step.push(Json::string(costs.state(e).value().get_str())); steps.push(step);
    }
    const auto bytes = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(source), 1, source_codec_limits(source_byte_limit())));
    result.set("steps", steps);
    result.set("native_complex_sha256", Json::string("sha256:" + sha256(
        std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    result.set("vertices", cardinality(source.grade({0}).cell_count()));
    result.set("coordinate_rule", Json::string("vertex v has native identity v+1; step e has relation "
                                             "identity 1000000+e and occurrences 2e+1,2e+2"));
    if (extras && record.has("view")) {
        const auto extras = render_record_view(record.at("view"));
        for (const auto &name : record.at("view").names()) result.set(name, extras.at(name));
    }
    return result;
}
inline NativeRecord source_record(AssayFieldsSource fields) {
    return NativeRecord::create_from_endpoint(std::move(fields.declaration), new_native_record_namespace(),
                                              std::move(fields.sources), std::move(fields.endpoint));
}
// Keep the complete parent publication for original byte context when selecting child fields.
inline N archive_child(const NativeRecord &child, std::vector<RawByteSource> &sources) {
    const auto key = rfc::machine::checked_size_add(sources.size(), 1);
    auto bytes = child.archive_source(key); const auto size = bytes.bytes().size();
    sources.push_back(std::move(bytes));
    return N(N::ByteSpan{key, 0, size, false});
}
inline N child_context(const NativeRecord &child) {
    N::Object result;
    member(result, "source_namespace", integer(child.source_namespace()));
    member(result, "cutoff", integer(child.endpoint()->sequence));
    N::Array carriers;
    for (const auto &carrier : child.endpoint()->carriers) {
        N::Object row;
        member(row, "key", integer(carrier.key)); member(row, "coordinate_system", integer(carrier.coordinate_system));
        member(row, "grade", integer(carrier.grade.value));
        element(carriers, std::move(row));
    }
    member(result, "carriers", std::move(carriers));
    for (const auto name : {"costs", "count_map"})
        if (child.has(name)) member(result, name, child.at(name).materialize_value());
    const auto name = child.has("comparisons") ? "comparisons" : "routes";
    member(result, "operation_family", std::string(name));
    N::Array operations;
    if (child.has(name)) for (std::size_t i = 0; i < child.at(name).size(); ++i) {
        N::Object row;
        const auto operation = child.at(name).at(i);
        for (const auto key : {"realized_field", "boundary", "disagreement", "reference_positions", "oriented_positions", "count"})
            if (operation.has(key)) member(row, key, operation.at(key).materialize_value());
        element(operations, std::move(row));
    }
    member(result, "operations", std::move(operations));
    return N(std::move(result));
}
inline N original_records(const SequenceSource &source, std::uint64_t key,
                          std::optional<std::size_t> cell_offset = {}) {
    N::Array rows;
    rows.elements.reserve(source.size());
    for (std::size_t i = 0; i < source.size(); ++i) {
        const auto &span = source.original_span(i);
        N row(N::Object{{{"id", N(N::ByteSpan{key, span.id_start, span.id_length})},
            {"record_index", integer(i)}, {"length", integer(source.length(i))},
            {"bytes", N(N::ByteSpan{key, span.start, span.end - span.start, false})},
            {"evidence", import_record_view(source.evidence(i))}}});
        if (cell_offset) row = N(N::Owned{0, rfc::storage::StableCellKey{{0}, *cell_offset + i + 1}, {},
            std::make_shared<const N>(std::move(row))});
        element(rows, std::move(row));
    }
    return N(std::move(rows));
}
inline N original_bases(const SequenceSource &source, std::uint64_t key,
                        std::size_t record, std::size_t begin, std::size_t length) {
    N::Array spans;
    for (const auto &[offset, count] : source.byte_spans(record, begin, length))
        element(spans, N(N::ByteSpan{key, offset, count}));
    return N(std::move(spans));
}
struct PlacementSource {
    std::size_t query, reference;
    std::optional<std::size_t> junction;
    std::size_t child, operation;
    std::vector<std::pair<std::size_t, std::size_t>> blocks;
};
inline AssayFieldsSource collection_source(const std::string &schema, const Json &request,
    const ResourceImage &queries, const ResourceImage &references,
    const SequenceSource &query_source, const SequenceSource &reference_source,
    const std::vector<PlacementSource> &placements, N::Object view, N::Array children, N::Array contexts,
    std::vector<RawByteSource> sources) {
    rfc::core::RelationSignature signature;
    const auto qrole = signature.aspects().add("query"), rrole = signature.aspects().add("reference"),
               jrole = signature.aspects().add("junction");
    const auto type = signature.add_type("sequence.returned_placement", rfc::core::BoundaryPolicy::DeclaredExact, 2, 3);
    const auto junctions = request.has("junctions") ? request.at("junctions").size() : 0;
    rfc::core::RelationalComplexBuilder builder(signature); builder.declare_grade({1});
    const auto records = rfc::machine::checked_size_add(query_source.size(), reference_source.size());
    const auto total = rfc::machine::checked_size_add(records, junctions);
    for (std::size_t i = 0; i < total; ++i) builder.declare_witness_identity(i + 1);
    std::uint64_t occurrence = 1;
    for (std::size_t i = 0; i < placements.size(); ++i) {
        const auto &p = placements[i];
        std::vector<rfc::core::BoundaryOccurrenceInput> participants;
        participants.push_back({p.query + 1, true, false, Rational{1}, {{qrole, true}}, rfc::core::OccurrenceKey{occurrence++}});
        participants.push_back({query_source.size() + p.reference + 1, true, false, Rational{1}, {{rrole, true}}, rfc::core::OccurrenceKey{occurrence++}});
        if (p.junction) participants.push_back({records + *p.junction + 1, true, false, Rational{1}, {{jrole, true}}, rfc::core::OccurrenceKey{occurrence++}});
        builder.add_grade_one(rfc::core::RelationKey{i + 1}, type, std::move(participants));
    }
    N::Object declaration;
    N::Array couplings;
    couplings.elements.reserve(placements.size());
    const auto &matches = std::get<N::Array>(std::find_if(view.members.begin(), view.members.end(),
        [](const auto &item) { return item.first == "matches"; })->second.data);
    const bool edited_splice = text(request.at("operation")) == "sequence_splice_discovery" &&
                               request.has("maximum_edits");
    for (std::size_t i = 0; i < placements.size(); ++i) {
        const auto &p = placements[i];
        N::Object coupling;
        member(coupling, "query_record", integer(p.query));
        member(coupling, "reference_record", integer(p.reference));
        member(coupling, "query_bases", original_bases(query_source, 1, p.query, 0, query_source.length(p.query)));
        N::Array blocks;
        for (const auto &[begin, end] : p.blocks) {
            N::Object block;
            member(block, "begin", integer(begin)); member(block, "end", integer(end));
            member(block, "bytes", original_bases(reference_source, 2, p.reference, begin, end - begin));
            element(blocks, std::move(block));
        }
        member(coupling, "reference_blocks", std::move(blocks));
        if (edited_splice) {
            N::Array introns;
            for (std::size_t b = 1; b < p.blocks.size(); ++b) {
                const auto begin = p.blocks[b - 1].second, end = p.blocks[b].first;
                require(begin <= end, "splice original intron span order");
                N::Object intron;
                member(intron, "begin", integer(begin)); member(intron, "end", integer(end));
                member(intron, "bytes", original_bases(reference_source, 2, p.reference, begin, end - begin));
                element(introns, std::move(intron));
            }
            member(coupling, "reference_introns", std::move(introns));
            const auto &match = std::get<N::Object>(matches.elements.at(i).data);
            member(coupling, "minimum_history_group", std::find_if(match.members.begin(), match.members.end(),
                [](const auto &item) { return item.first == "minimum_history_group"; })->second);
        }
        member(coupling, "child_publication", children.elements.at(p.child));
        member(coupling, "child_context", integer(p.child));
        member(coupling, "operation_index", integer(p.operation));
        const auto &match = std::get<N::Object>(matches.elements.at(i).data);
        const auto family = std::find_if(match.members.begin(), match.members.end(),
            [](const auto &item) { return item.first == "minimum_edit_histories"; });
        if (family != match.members.end()) {
            N::Array operations;
            for (const auto &history : std::get<N::Array>(family->second.data).elements) {
                const auto &members = std::get<N::Object>(history.data).members;
                element(operations, std::find_if(members.begin(), members.end(),
                    [](const auto &item) { return item.first == "route_index"; })->second);
            }
            member(coupling, "minimum_history_operations", std::move(operations));
        }
        const auto history_evidence = std::find_if(match.members.begin(), match.members.end(),
            [](const auto &item) { return item.first == "minimum_history_evidence"; });
        if (history_evidence != match.members.end()) member(coupling, "minimum_history_evidence", history_evidence->second);
        element(couplings, N(N::Owned{0, rfc::storage::StableCellKey{{1}, i + 1}, {},
            std::make_shared<const N>(std::move(coupling))}));
    }
    member(declaration, "schema", schema); member(declaration, "protocol", import_record_view(request));
    member(declaration, "queries", original_records(query_source, 1, 0));
    member(declaration, "references", original_records(reference_source, 2, query_source.size()));
    member(declaration, "query_bytes", N(N::ByteSpan{1, 0, queries.bytes.size(), false}));
    member(declaration, "reference_bytes", N(N::ByteSpan{2, 0, references.bytes.size(), false}));
    if (request.has("junctions")) member(declaration, "junctions", import_record_view(request.at("junctions")));
    member(declaration, "children", std::move(children)); member(declaration, "view", std::move(view));
    member(declaration, "placements", std::move(couplings));
    member(declaration, "child_contexts", std::move(contexts));
    return {source_endpoint(builder.freeze()), N(std::move(declaration)), std::move(sources)};
}
inline Json routed_collection_view(const NativeRecord &record) {
    auto result = render_record_view(record.at("view")), carriers = Json::array();
    const auto matches = record.at("view").at("matches"), children = record.at("children");
    std::vector<std::vector<std::size_t>> indices(children.size());
    for (std::size_t i = 0; i < matches.size(); ++i)
        indices.at(matches.at(i).at("carrier").integer().get_ui()).push_back(i);
    for (std::size_t c = 0; c < children.size(); ++c) {
        const auto child = NativeRecord::decode(children.at(c).bytes());
        carriers.push(carrier_fields_view(child));
        for (const auto i : indices[c]) {
            const auto route = matches.at(i).at("route_index").integer().get_ui();
            auto match = result.at("matches").at(i);
            json_object_object_del(match.raw(), "route_index");
            match.set("route", retained_route_view(child.at("routes").at(route)));
            if (matches.at(i).has("minimum_edit_histories")) {
                const auto histories = matches.at(i).at("minimum_edit_histories");
                for (std::size_t h = 0; h < histories.size(); ++h) {
                    auto history = match.at("minimum_edit_histories").at(h);
                    const auto index = histories.at(h).at("route_index").integer().get_ui();
                    json_object_object_del(history.raw(), "route_index");
                    history.set("route", retained_route_view(child.at("routes").at(index)));
                }
            }
            json_object_object_del(match.raw(), "minimum_history_evidence");
        }
    }
    result.set("carriers", carriers);
    if (result.has("minimum_history_groups"))
        for (std::size_t g = 0; g < result.at("minimum_history_groups").size(); ++g)
            json_object_object_del(result.at("minimum_history_groups").at(g).raw(), "minimum_history_evidence");
    return result;
}
} // namespace gadget::assay_detail::sequence
