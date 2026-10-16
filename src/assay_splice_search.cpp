#include "sequence_route_witnesses.hpp"
#include <array>
#include <numeric>

namespace gadget::assay_detail {
namespace {
struct EditedSpliceState {
    std::size_t query, reference, edits, junctions, query_span, reference_span;
    auto operator<=>(const EditedSpliceState &) const = default;
};
struct EditedSpliceStep {
    std::size_t from, to;
    char operation;
    std::string motif;
};
AssayFieldsSource edited_splice_fields_source(
    const Json &request, const ResourceImage &query_image, const ResourceImage &reference_image,
    const SequenceSource &query_source, const SequenceSource &reference_source,
    const std::vector<SequenceRecord> &queries, const std::vector<SequenceRecord> &refs,
    std::size_t overhang, std::size_t min_gap, std::size_t max_gap,
    std::size_t max_junctions, const std::set<std::string> &motifs) {
    using namespace sequence;
    require(text(request.at("scoring")) == "unit_edit_cost", "edited splice requires unit edit costs");
    const auto errors = count(request.at("maximum_edits"), std::numeric_limits<std::size_t>::max(),
                              "splice edit budget representation");
    const auto orientations = strands(request);
    for (const auto &query : queries) {
        require(query.bases.size() >= overhang, "splice query length");
        for (const auto base : query.bases)
            require(dna(base) >= 0, "splice query requires literal ACGT");
    }
    N::Array matches, children, summaries, contexts, groups;
    std::vector<PlacementSource> placements;
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(query_image.bytes))),
        RawByteSource(2, std::as_bytes(std::span(reference_image.bytes)))};
    std::vector<rfc::numeric::Integer> counts(queries.size(), rfc::numeric::Integer{0});
    rfc::numeric::Integer starts{0}, vertices{0}, relations{0}, histories{0}, native_calls{0}, witness_calls{0};
    for (std::size_t r = 0; r < refs.size(); ++r) {
        const auto &reference = refs[r].bases;
        std::array<std::vector<std::size_t>, 16> acceptors;
        for (std::size_t a = 2; a <= reference.size(); ++a)
            if (dna(reference[a - 2]) >= 0 && dna(reference[a - 1]) >= 0)
                acceptors.at(token(reference.substr(a - 2, 2))).push_back(a);
        for (const auto &strand : orientations) {
            struct Motif { std::string donor, name; std::size_t acceptor; };
            std::vector<Motif> boundaries;
            for (const auto &motif : motifs) {
                const auto donor = strand == "+" ? motif.substr(0, 2) : reverse_complement(motif.substr(3));
                const auto acceptor = strand == "+" ? motif.substr(3) : reverse_complement(motif.substr(0, 2));
                boundaries.push_back({donor, motif, token(acceptor)});
            }
            for (std::size_t q = 0; q < queries.size(); ++q) {
                const auto query = strand == "+" ? queries[q].bases : reverse_complement(queries[q].bases);
                // Each exon consumes at least overhang query bases. This bound
                // removes only junction counts with no possible complete query.
                const auto junctions = std::min(max_junctions, query.size() / overhang - 1);
                for (std::size_t start = 0; start < reference.size(); ++start) {
                    ++starts;
                    const auto remaining = reference.size() - start;
                    if (remaining < overhang) continue;
                    auto width = std::min(query.size(), remaining);
                    width += std::min(errors, remaining - width);
                    if (junctions != 0) {
                        const auto capacity = remaining - width;
                        width += max_gap > capacity / junctions ? capacity :
                            rfc::machine::checked_size_mul(max_gap, junctions);
                    }
                    const auto end = start + width;
                    std::vector<EditedSpliceState> states{{0, start, 0, 0, 0, 0}};
                    std::map<EditedSpliceState, std::size_t> coordinates{{states.front(), 0}};
                    std::vector<EditedSpliceStep> steps;
                    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> endpoints;
                    const auto add = [&](std::size_t from, EditedSpliceState state, char op,
                                         std::string motif = {}) {
                        const auto [it, inserted] = coordinates.emplace(state, states.size());
                        if (inserted) {
                            require(states.size() < std::numeric_limits<std::uint32_t>::max(),
                                    "edited splice native vertex coordinate representation");
                            states.push_back(state);
                        }
                        require(steps.size() < std::numeric_limits<std::uint32_t>::max(),
                                "edited splice native relation coordinate representation");
                        steps.push_back({from, it->second, op, std::move(motif)});
                    };
                    for (std::size_t v = 0; v < states.size(); ++v) {
                        const auto state = states[v];
                        if (state.query == query.size() && state.query_span == overhang &&
                            state.reference_span == overhang)
                            endpoints[{state.reference, state.junctions}].push_back(v);
                        if (state.query < query.size() && state.edits < errors) {
                            auto next = state; ++next.query; ++next.edits;
                            next.query_span += next.query_span < overhang;
                            add(v, next, 'I');
                        }
                        if (state.reference < end && state.edits < errors) {
                            auto next = state; ++next.reference; ++next.edits;
                            next.reference_span += next.reference_span < overhang;
                            add(v, next, 'D');
                        }
                        if (state.query < query.size() && state.reference < end) {
                            const bool equal = query[state.query] == reference[state.reference];
                            if (equal || state.edits < errors) {
                                auto next = state; ++next.query; ++next.reference;
                                next.edits += !equal;
                                next.query_span += next.query_span < overhang;
                                next.reference_span += next.reference_span < overhang;
                                add(v, next, equal ? '=' : 'X');
                            }
                        }
                        if (state.query_span != overhang || state.reference_span != overhang ||
                            state.junctions == junctions || query.size() - state.query < overhang ||
                            end - state.reference < min_gap)
                            continue;
                        const auto lo = state.reference + min_gap;
                        const auto hi = state.reference + std::min(max_gap, end - state.reference);
                        for (const auto &motif : boundaries) {
                            if (reference.compare(state.reference, 2, motif.donor) != 0) continue;
                            const auto &offsets = acceptors.at(motif.acceptor);
                            for (auto it = std::lower_bound(offsets.begin(), offsets.end(), lo);
                                 it != offsets.end() && *it <= hi; ++it) {
                                if (end - *it < overhang) continue;
                                auto next = state; next.reference = *it; ++next.junctions;
                                next.query_span = 0; next.reference_span = 0;
                                add(v, next, 'N', motif.name);
                            }
                        }
                    }
                    if (endpoints.empty()) continue;
                    const auto biological_vertices = states.size();
                    const auto vertex_count = rfc::machine::checked_size_add(states.size(), endpoints.size());
                    Routes routes(vertex_count, "sequence.edited_splice.step");
                    for (const auto &step : steps) {
                        const auto &from = states.at(step.from), &to = states.at(step.to);
                        require(to.query >= from.query && to.reference >= from.reference &&
                                    (to.query > from.query || to.reference > from.reference),
                                "edited splice biological declarations must strictly advance prefixes");
                        routes.add(step.from, step.to,
                            Rational(step.operation == 'I' || step.operation == 'D' || step.operation == 'X' ? 1 : 0));
                    }
                    std::vector<std::size_t> terminals;
                    for (const auto &[key, admitted] : endpoints) {
                        (void)key;
                        const auto sink = biological_vertices + terminals.size();
                        terminals.push_back(sink);
                        for (const auto v : admitted) {
                            routes.add(v, sink, Rational{0});
                            steps.push_back({v, sink, 'T', {}});
                        }
                    }
                    const auto source = routes.freeze();
                    const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
                    vertices += rfc::numeric::Integer(std::to_string(vertex_count));
                    relations += rfc::numeric::Integer(std::to_string(routes.steps.size()));
                    auto native = route_source(source, routes);
                    SourceGradeFields fields(native.endpoint);
                    N::Array native_routes, operations, terminal_values;
                    for (const auto &step : steps) element(operations, std::string(1, step.operation));
                    for (const auto terminal : terminals) element(terminal_values, integer(terminal));
                    N::Object carrier;
                    member(carrier, "source_vertex", integer(0));
                    member(carrier, "operations", std::move(operations));
                    member(carrier, "terminal_vertices", std::move(terminal_values));
                    member(carrier, "reference_record_index", integer(r));
                    member(carrier, "reference_start", integer(start));
                    member(carrier, "query_record_index", integer(q));
                    member(carrier, "strand", strand);
                    const auto child_index = children.elements.size();
                    std::size_t endpoint_index = 0;
                    for (const auto &[key, admitted] : endpoints) {
                        (void)admitted;
                        const auto sink = terminals.at(endpoint_index++);
                        const auto route = rfc::path::dijkstra_exact(graph, graph.object_index(1),
                                                                    graph.object_index(sink + 1));
                        ++native_calls;
                        require(route.reachable && route.distance <= Rational(std::to_string(errors)),
                                "declared edited splice endpoint must have an admitted native route");
                        const CompleteRouteWitnesses witnesses(graph, graph.object_index(1),
                            graph.object_index(sink + 1), true, route);
                        witnesses.verify(rfc::field::boundary_matrix(source, {1}), routes.costs(), 0, sink);
                        witness_calls += rfc::numeric::Integer(std::to_string(witnesses.histories.size() - 1));
                        witness_calls += rfc::numeric::Integer(std::to_string(witnesses.unreachable_suffix_vertices.size()));
                        const auto primary = native_routes.elements.size();
                        N::Array witness_indices;
                        for (const auto &history : witnesses.histories) {
                            element(witness_indices, integer(native_routes.elements.size()));
                            element(native_routes, retain_route(fields, source, routes, history.route.relations,
                                history.from, history.to, history.route.distance));
                        }
                        N::Object evidence;
                        member(evidence, "prefix", fields.field("minimum_history.prefix", 0, witnesses.prefix));
                        std::vector<std::size_t> viable;
                        for (std::size_t v = 0; v < witnesses.suffix.size(); ++v)
                            if (witnesses.suffix_reachable[v]) viable.push_back(v);
                        member(evidence, "suffix", fields.column("minimum_history.suffix", 0, vertex_count,
                            [&](auto v) { return witnesses.suffix[v]; },
                            rfc::numeric::StructuralSupport::sparse(vertex_count, viable)));
                        member(evidence, "witness_routes", std::move(witness_indices));
                        const auto tight = [&](std::size_t e) {
                            const auto &step = routes.steps.at(e);
                            return witnesses.suffix_reachable[step.to] &&
                                witnesses.prefix[step.from] + step.cost + witnesses.suffix[step.to] == route.distance;
                        };
                        member(evidence, "optimal_support", fields.column("minimum_history.optimal", 1,
                            routes.steps.size(), [&](auto e) { return Rational(tight(e) ? 1 : 0); }));
                        const auto group_index = groups.elements.size();
                        rfc::numeric::Integer family_count{0};
                        std::vector<std::uint32_t> ordered;
                        struct Frame { std::uint32_t vertex; std::size_t next; };
                        std::vector<Frame> stack{{0, 0}};
                        const auto path_width = rfc::machine::checked_size_add(query.size(), width);
                        ordered.reserve(rfc::machine::checked_size_add(path_width, 1));
                        stack.reserve(rfc::machine::checked_size_add(path_width, 2));
                        const auto emit = [&]() {
                            N::Object match;
                            N::Array blocks, junction_values;
                            PlacementSource placement{q, r, {}, child_index,
                                ordered == route.relations ? primary : native_routes.elements.size(), {}};
                            if (placement.operation != primary) element(native_routes, retain_route(fields,
                                source, routes, ordered, 0, sink, route.distance));
                            auto query_begin = std::size_t{0}, reference_begin = start;
                            std::string aligned_query, aligned_reference, reconstructed_query;
                            const auto block = [&](std::size_t qe, std::size_t re) {
                                require(qe - query_begin >= overhang && re - reference_begin >= overhang &&
                                    aligned_query.size() == aligned_reference.size(), "edited splice exon span witness");
                                N::Object value;
                                member(value, "query_begin", integer(query_begin)); member(value, "query_end", integer(qe));
                                member(value, "start", integer(reference_begin)); member(value, "end", integer(re));
                                member(value, "aligned_query", aligned_query); member(value, "aligned_reference", aligned_reference);
                                element(blocks, std::move(value)); placement.blocks.emplace_back(reference_begin, re);
                                reconstructed_query += query.substr(query_begin, qe - query_begin);
                                aligned_query.clear(); aligned_reference.clear();
                            };
                            rfc::numeric::Integer observed_edits{0};
                            for (const auto e : ordered) {
                                const auto &step = steps.at(e);
                                const auto &from = states.at(step.from);
                                if (step.operation == 'T') { block(from.query, from.reference); continue; }
                                const auto &to = states.at(step.to);
                                if (step.operation == 'N') {
                                    block(from.query, from.reference);
                                    N::Object junction;
                                    member(junction, "donor", integer(from.reference));
                                    member(junction, "acceptor", integer(to.reference));
                                    member(junction, "motif", step.motif);
                                    member(junction, "origin", "reference_sequence_and_declared_edit_history");
                                    element(junction_values, std::move(junction));
                                    query_begin = to.query; reference_begin = to.reference;
                                } else {
                                    aligned_query += step.operation == 'D' ? '-' : query.at(from.query);
                                    aligned_reference += step.operation == 'I' ? '-' : reference.at(from.reference);
                                    observed_edits += step.operation != '=';
                                }
                            }
                            require(reconstructed_query == query && Rational(observed_edits) == route.distance &&
                                        junction_values.elements.size() == key.second,
                                    "edited splice complete query, edit count and junction witnesses");
                            member(match, "query_id", queries[q].id); member(match, "query_record_index", integer(q));
                            member(match, "reference_id", refs[r].id); member(match, "reference_record_index", integer(r));
                            member(match, "strand", strand); member(match, "start", integer(start));
                            member(match, "end", integer(key.first)); member(match, "minimum_edits", route.distance);
                            member(match, "blocks", std::move(blocks)); member(match, "junctions", std::move(junction_values));
                            member(match, "query_evidence", import_record_view(queries[q].evidence));
                            member(match, "reference_evidence", import_record_view(refs[r].evidence));
                            member(match, "carrier", integer(child_index)); member(match, "route_index", integer(placement.operation));
                            member(match, "minimum_history_group", integer(group_index));
                            member(match, "history_index", family_count);
                            member(match, "query_reconstruction_verified", true);
                            element(matches, std::move(match)); placements.push_back(std::move(placement));
                            ++family_count; ++histories; ++counts[q];
                        };
                        while (!stack.empty()) {
                            auto &frame = stack.back();
                            const auto &outgoing = graph.outgoing(frame.vertex);
                            if (frame.vertex == sink) {
                                emit(); stack.pop_back(); if (!stack.empty()) ordered.pop_back();
                            } else if (frame.next == outgoing.size()) {
                                stack.pop_back(); if (!stack.empty()) ordered.pop_back();
                            } else {
                                const auto &edge = graph.edge(outgoing[frame.next++]);
                                if (!tight(edge.source_relation)) continue;
                                ordered.push_back(edge.source_relation); stack.push_back({edge.target_object, 0});
                            }
                        }
                        require(family_count > 0, "edited splice minimum history family is empty");
                        N::Object group;
                        member(group, "query_record_index", integer(q)); member(group, "reference_record_index", integer(r));
                        member(group, "strand", strand); member(group, "start", integer(start));
                        member(group, "end", integer(key.first)); member(group, "junctions", integer(key.second));
                        member(group, "minimum_edits", route.distance); member(group, "histories", std::move(family_count));
                        member(group, "carrier", integer(child_index)); member(group, "sink_vertex", integer(sink));
                        member(group, "minimum_history_evidence", std::move(evidence)); element(groups, std::move(group));
                    }
                    auto &declaration = std::get<N::Object>(native.declaration.data);
                    member(declaration, "view", std::move(carrier)); member(declaration, "routes", std::move(native_routes));
                    std::vector<std::size_t> coordinate_support(biological_vertices);
                    std::iota(coordinate_support.begin(), coordinate_support.end(), 0);
                    const auto support = rfc::numeric::StructuralSupport::sparse(vertex_count, coordinate_support);
                    const auto coordinate = [&](std::string name, auto select) {
                        member(declaration, name, fields.column("splice." + name, 0, vertex_count,
                            [&](auto v) { return Rational(std::to_string(select(states.at(v)))); }, support));
                    };
                    coordinate("query_prefix", [](const auto &s) { return s.query; });
                    coordinate("reference_prefix", [](const auto &s) { return s.reference; });
                    coordinate("edit_count", [](const auto &s) { return s.edits; });
                    coordinate("junction_count", [](const auto &s) { return s.junctions; });
                    coordinate("query_exon_span_capped", [](const auto &s) { return s.query_span; });
                    coordinate("reference_exon_span_capped", [](const auto &s) { return s.reference_span; });
                    member(declaration, "exon_span_cap", integer(overhang));
                    member(declaration, "operation_codes", fields.column("splice.operation", 1, steps.size(),
                        [&](auto e) { return Rational(steps[e].operation == '=' ? 0 : steps[e].operation == 'X' ? 1 :
                            steps[e].operation == 'I' ? 2 : steps[e].operation == 'D' ? 3 : steps[e].operation == 'N' ? 4 : 5); }));
                    N::Array annotations;
                    for (std::size_t e = 0; e < steps.size(); ++e) if (steps[e].operation == 'N') {
                        N::Object annotation; member(annotation, "motif", steps[e].motif);
                        element(annotations, N(N::Owned{0, rfc::storage::StableCellKey{{1}, 1000000 + e}, {},
                            std::make_shared<const N>(std::move(annotation))}));
                    }
                    member(declaration, "junction_annotations", std::move(annotations));
                    const auto child = source_record(std::move(native));
                    element(contexts, child_context(child)); element(children, archive_child(child, sources));
                }
            }
        }
    }
    for (std::size_t q = 0; q < queries.size(); ++q) {
        N::Object summary; member(summary, "query_id", queries[q].id);
        member(summary, "query_record_index", integer(q)); member(summary, "matches", counts[q]);
        element(summaries, std::move(summary));
    }
    N::Object work, out;
    member(work, "reference_starts", std::move(starts)); member(work, "vertices", std::move(vertices));
    member(work, "relations", std::move(relations)); member(work, "native_route_calls", std::move(native_calls));
    member(work, "native_witness_route_calls", std::move(witness_calls)); member(work, "minimum_edit_histories", std::move(histories));
    member(out, "matches", std::move(matches)); member(out, "queries", std::move(summaries));
    member(out, "minimum_history_groups", std::move(groups)); member(out, "work", std::move(work));
    member(out, "complete", true); member(out, "coordinate_convention", "zero_based_half_open_reference_orientation");
    member(out, "search_universe", "Every minimum unit-edit ordered history per original query, strand, reference start/end and exact junction count, with declared edit/intron/junction bounds and minimum consumed query and reference span per exon. Introns have zero edit cost and transcript-oriented declared sequence motifs.");
    member(out, "completeness_argument", "Every source start is considered. Reachable biological prefix states retain accrued edits, junctions and both exon spans. Every admitted match, substitution, insertion, deletion and motif jump is a native source relation. Exon spans are capped only after meeting their declared minimum. Every history strictly advances query/reference prefixes before one terminal exit. Canonical minimum route certificates select all tight histories per endpoint/junction count; iterative enumeration retains each ordered history without deduplication or truncation.");
    return collection_source("gadget.sequence.splice.fields.v1", request, query_image, reference_image,
        query_source, reference_source, placements, std::move(out), std::move(children), std::move(contexts), std::move(sources));
}
} // namespace
AssayFieldsSource splice_search_fields_source(const ResourceImages &images, const Json &request) {
    using namespace sequence;
    members(request, {"schema", "operation", "query_resource", "reference_resource", "strands",
                      "minimum_exon_bases", "minimum_intron_bases", "maximum_intron_bases",
                      "maximum_junctions", "motifs", "evidence"}, {"maximum_edits", "scoring"});
    require(request.has("maximum_edits") == request.has("scoring"),
            "edited splice requires both maximum_edits and scoring");
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
    if (request.has("maximum_edits"))
        return edited_splice_fields_source(request, query_image, reference_image, query_source,
            reference_source, queries, refs, overhang, min_gap, max_gap, max_junctions, motifs);
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
            // Seed hits are retained in reference offset order. Exact interval
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
        // Enumerate biological exon chain declarations, preserving all ties.
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
        if (!match.at("blocks").at(0).has("exon_identity")) continue;
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
