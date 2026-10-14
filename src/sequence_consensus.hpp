#pragma once
#include "sequence_route_witnesses.hpp"

namespace gadget::assay_detail::consensus {
using namespace sequence;
struct ReadState {
    std::size_t prefix{}, paired{};
    unsigned char phase{};
    auto operator<=>(const ReadState &) const = default;
};
struct State {
    std::vector<ReadState> reads;
    std::size_t length{}, edits{};
    auto operator<=>(const State &) const = default;
};
struct Action {
    ReadState next;
    std::size_t insertions{}, edits{};
    char operation{'O'};
};
struct Transition {
    std::size_t from{}, to{};
    char base{};
    std::vector<Action> actions;
    std::size_t edits{};
};
struct Hypotheses {
    std::vector<State> states;
    std::vector<Transition> transitions;
    std::size_t sink{};
    Hypotheses(const std::vector<std::string> &reads, std::size_t budget, std::size_t paired) {
        State initial; initial.reads.resize(reads.size()); states.push_back(initial);
        std::map<State, std::size_t> index{{initial, 0}};
        std::vector<Transition> tails;
        const auto declare = [&](std::size_t from, State next, char base,
                                 std::vector<Action> actions, std::size_t edits) {
            const auto [it, inserted] = index.emplace(next, states.size());
            if (inserted) {
                require(states.size() < std::numeric_limits<std::uint32_t>::max(),
                        "consensus native vertex coordinate representation");
                states.push_back(std::move(next));
            }
            require(transitions.size() < std::numeric_limits<std::uint32_t>::max(),
                    "consensus native relation coordinate representation");
            transitions.push_back({from, it->second, base, std::move(actions), edits});
        };
        for (std::size_t v = 0; v < states.size(); ++v) {
            const auto state = states[v];
            const auto available = budget - state.edits;
            bool finish = true;
            std::size_t tail_cost = 0;
            std::vector<Action> tail;
            tail.reserve(reads.size());
            for (std::size_t r = 0; r < reads.size(); ++r) {
                const auto &s = state.reads[r];
                const auto count = reads[r].size() - s.prefix;
                if (s.phase == 0 || s.paired < paired || count > available - tail_cost) {
                    finish = false; break;
                }
                tail_cost += count;
                tail.push_back({{reads[r].size(), s.paired, 2}, count, count, 'O'});
            }
            if (finish) tails.push_back({v, 0, 0, std::move(tail), tail_cost});
            for (const char base : std::string_view("ACGT")) {
                std::vector<std::vector<Action>> options(reads.size());
                for (std::size_t r = 0; r < reads.size(); ++r) {
                    const auto s = state.reads[r];
                    auto &choices = options[r];
                    if (s.phase != 1) choices.push_back({s, 0, 0, 'O'});
                    if (s.phase == 2) continue;
                    const auto remaining = reads[r].size() - s.prefix;
                    for (std::size_t n = 0; n <= std::min(available, remaining); ++n) {
                        const auto after = s.prefix + n;
                        if (after == reads[r].size()) {
                            if (s.phase == 1 && s.paired >= paired)
                                choices.push_back({{after, s.paired, 2}, n, n, 'O'});
                            continue;
                        }
                        const bool equal = reads[r][after] == base;
                        const auto aligned = s.paired + (s.paired < paired);
                        const auto prefix = after + 1;
                        // Paired counts are capped only after meeting the
                        // declared minimum. No remaining token can be reused.
                        if ((equal || n < available) &&
                            aligned + (reads[r].size() - prefix) >= paired) {
                            const auto phase = static_cast<unsigned char>(prefix == reads[r].size() ? 2 : 1);
                            choices.push_back({{prefix, aligned, phase}, n, n + !equal, equal ? '=' : 'X'});
                        }
                        if (s.phase == 1 && n < available)
                            choices.push_back({{after, s.paired, 1}, n, n + 1, 'D'});
                    }
                }
                if (std::ranges::any_of(options, [](const auto &v) { return v.empty(); })) continue;
                // Joint choices have one declaration in original record order.
                // Independent insertion interleavings create no extra histories.
                std::vector<std::size_t> selected(reads.size());
                for (;;) {
                    State next = state;
                    std::vector<Action> actions; actions.reserve(reads.size());
                    std::size_t edits = 0;
                    bool fits = true, supported = false;
                    for (std::size_t r = 0; r < reads.size(); ++r) {
                        const auto action = options[r][selected[r]];
                        if (action.edits > available - edits) { fits = false; break; }
                        edits += action.edits; supported |= action.operation == '=';
                        next.reads[r] = action.next; actions.push_back(action);
                    }
                    if (fits && supported) {
                        next.length = rfc::machine::checked_size_add(state.length, 1);
                        next.edits += edits;
                        require(std::ranges::any_of(actions, [&](const auto &a) { return a.operation == '='; }),
                                "consensus emitted token must have exact original support");
                        declare(v, std::move(next), base, std::move(actions), edits);
                    }
                    auto r = selected.size();
                    while (r != 0) {
                        --r;
                        if (++selected[r] < options[r].size()) break;
                        selected[r] = 0;
                    }
                    if (r == 0 && selected[0] == 0) break;
                }
            }
        }
        sink = states.size();
        for (const auto &step : transitions) {
            const auto &from = states.at(step.from), &to = states.at(step.to);
            std::size_t before = 0, after = 0;
            for (std::size_t r = 0; r < reads.size(); ++r) {
                require(to.reads[r].prefix >= from.reads[r].prefix &&
                            to.reads[r].phase >= from.reads[r].phase,
                        "consensus original read prefixes and placement phases must be monotone");
                before = rfc::machine::checked_size_add(before, from.reads[r].prefix);
                after = rfc::machine::checked_size_add(after, to.reads[r].prefix);
            }
            require(to.length == from.length + 1 && after > before &&
                        to.edits == from.edits + step.edits,
                    "consensus biological emission must advance source tokens and output prefix");
        }
        for (auto &step : tails) {
            step.to = sink;
            transitions.push_back(std::move(step));
        }
    }
};
struct History {
    std::size_t route{};
    std::vector<std::string> strands;
    std::string spelling;
    std::vector<std::uint32_t> relations;
    std::size_t edits{};
};
} // namespace gadget::assay_detail::consensus

namespace gadget::assay_detail {
AssayFieldsSource read_consensus_fields_source(const ResourceImages &images, const Json &request) {
    using namespace sequence;
    members(request, {"schema", "operation", "resource", "strands", "maximum_edits",
                      "minimum_paired_bases", "objective", "ambiguity", "evidence"});
    require(text(request.at("objective")) == "minimum_length_then_total_unit_edits" &&
                text(request.at("ambiguity")) == "literal_ACGT", "read consensus profile");
    (void)text(request.at("evidence"), source_byte_limit());
    const auto budget = count(request.at("maximum_edits"), std::numeric_limits<std::size_t>::max(),
                              "consensus total edit budget representation");
    const auto minimum_paired = count(request.at("minimum_paired_bases"),
        std::numeric_limits<std::size_t>::max(), "consensus paired-base declaration");
    require(minimum_paired > 0, "consensus minimum paired bases must be positive");
    const auto allowed = strands(request);
    const auto &resource = images.at(text(request.at("resource")));
    const SequenceSource original(resource);
    const auto reads = records(original);
    require(!reads.empty(), "consensus requires a nonempty read collection");
    std::size_t total_tokens = 0;
    bool possible = true;
    for (const auto &read : reads) {
        require(!read.bases.empty(), "consensus refuses empty reads");
        for (const auto base : read.bases) require(dna(base) >= 0, "consensus requires literal ACGT reads");
        total_tokens = rfc::machine::checked_size_add(total_tokens, read.bases.size());
        possible &= read.bases.size() >= minimum_paired;
    }
    const Rational weight(rfc::numeric::Integer(std::to_string(budget)) + 1);
    struct IndexedState { consensus::State value; std::size_t orientation{}; };
    std::vector<std::optional<IndexedState>> states(1);
    std::vector<consensus::Transition> transitions;
    std::vector<std::vector<std::string>> orientations;
    std::vector<std::size_t> choices(reads.size());
    if (possible) for (;;) {
        std::vector<std::string> bases, orientation;
        for (std::size_t r = 0; r < reads.size(); ++r) {
            const auto &strand = allowed.at(choices[r]);
            orientation.push_back(strand);
            bases.push_back(strand == "+" ? reads[r].bases : reverse_complement(reads[r].bases));
        }
        consensus::Hypotheses declared(bases, budget, minimum_paired);
        const auto offset = states.size(), orientation_index = orientations.size();
        require(rfc::machine::checked_size_add(offset, declared.states.size()) <
                    std::numeric_limits<std::uint32_t>::max(), "consensus native joint vertex representation");
        orientations.push_back(std::move(orientation));
        transitions.push_back({0, offset, 0, {}, 0});
        for (auto &state : declared.states) states.push_back(IndexedState{std::move(state), orientation_index});
        for (auto &step : declared.transitions) {
            step.from += offset;
            step.to = step.to == declared.sink ? std::numeric_limits<std::size_t>::max() : step.to + offset;
            transitions.push_back(std::move(step));
        }
        auto r = choices.size();
        while (r != 0) { --r; if (++choices[r] < allowed.size()) break; choices[r] = 0; }
        if (r == 0 && choices[0] == 0) break;
    }
    const auto sink = states.size(); states.push_back({});
    Routes routes(states.size(), "sequence.consensus.optimization_step");
    for (auto &step : transitions) {
        if (step.to == std::numeric_limits<std::size_t>::max()) step.to = sink;
        routes.add(step.from, step.to, (step.base ? weight : Rational{0}) + Rational(std::to_string(step.edits)));
    }
    const auto source = routes.freeze();
    const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
    const auto optimum = rfc::path::dijkstra_exact(graph, graph.object_index(1), graph.object_index(sink + 1));
    auto optimization = route_source(source, routes);
    SourceGradeFields optimization_fields(optimization.endpoint);
    N::Array native_routes;
    std::vector<consensus::History> histories;
    N::Object certificates;
    if (optimum.reachable) {
        const CompleteRouteWitnesses witnesses(graph, graph.object_index(1), graph.object_index(sink + 1), true, optimum);
        witnesses.verify(rfc::field::boundary_matrix(source, {1}), routes.costs(), 0, sink);
        N::Array indices;
        for (const auto &history : witnesses.histories) {
            element(indices, integer(native_routes.elements.size()));
            element(native_routes, retain_route(optimization_fields, source, routes, history.route.relations,
                history.from, history.to, history.route.distance));
        }
        member(certificates, "prefix", optimization_fields.field("consensus.prefix", 0, witnesses.prefix));
        std::vector<std::size_t> viable;
        for (std::size_t v = 0; v < states.size(); ++v) if (witnesses.suffix_reachable[v]) viable.push_back(v);
        member(certificates, "suffix", optimization_fields.column("consensus.suffix", 0, states.size(),
            [&](auto v) { return witnesses.suffix[v]; }, rfc::numeric::StructuralSupport::sparse(states.size(), viable)));
        member(certificates, "witness_routes", std::move(indices));
        const auto tight = [&](std::size_t e) {
            const auto &step = routes.steps.at(e);
            return witnesses.suffix_reachable[step.to] &&
                witnesses.prefix[step.from] + step.cost + witnesses.suffix[step.to] == optimum.distance;
        };
        member(certificates, "optimal_support", optimization_fields.column("consensus.optimal", 1,
            routes.steps.size(), [&](auto e) { return Rational(tight(e) ? 1 : 0); }));
        std::vector<std::uint32_t> ordered;
        struct Frame { std::uint32_t vertex; std::size_t next; };
        std::vector<Frame> stack{{0, 0}};
        ordered.reserve(rfc::machine::checked_size_add(total_tokens, 2));
        stack.reserve(rfc::machine::checked_size_add(total_tokens, 3));
        while (!stack.empty()) {
            auto &frame = stack.back(); const auto &outgoing = graph.outgoing(frame.vertex);
            if (frame.vertex == sink) {
                consensus::History history;
                history.relations = ordered;
                history.route = ordered == optimum.relations ? 0 : native_routes.elements.size();
                if (history.route != 0) element(native_routes, retain_route(optimization_fields, source, routes,
                    ordered, 0, sink, optimum.distance));
                for (const auto e : ordered) {
                    const auto &step = transitions.at(e);
                    if (step.from == 0) {
                        history.strands = orientations.at(states.at(step.to)->orientation); continue;
                    }
                    if (step.base) history.spelling += step.base;
                    history.edits = rfc::machine::checked_size_add(history.edits, step.edits);
                }
                require(history.spelling.size() <= total_tokens && history.edits <= budget &&
                    weight * Rational(std::to_string(history.spelling.size())) + Rational(std::to_string(history.edits)) == optimum.distance,
                    "consensus complete length and unit-edit witness");
                histories.push_back(std::move(history));
                stack.pop_back(); if (!stack.empty()) ordered.pop_back();
            } else if (frame.next == outgoing.size()) {
                stack.pop_back(); if (!stack.empty()) ordered.pop_back();
            } else {
                const auto &edge = graph.edge(outgoing[frame.next++]);
                if (!tight(edge.source_relation)) continue;
                ordered.push_back(edge.source_relation); stack.push_back({edge.target_object, 0});
            }
        }
    }
    auto &optimization_declaration = std::get<N::Object>(optimization.declaration.data);
    member(optimization_declaration, "routes", std::move(native_routes));
    member(optimization_declaration, "certificates", std::move(certificates));
    N::Array tuple_values;
    for (const auto &tuple : orientations) { N::Array values; for (const auto &s : tuple) element(values, s); element(tuple_values, std::move(values)); }
    N::Object carrier_view; member(carrier_view, "source_vertex", integer(0));
    member(carrier_view, "sink_vertex", integer(sink)); member(carrier_view, "orientations", std::move(tuple_values));
    member(carrier_view, "optimization_contract", "Binary hypothesis optimization edges; biological joint token evidence is retained separately in the parent complex.");
    member(optimization_declaration, "view", std::move(carrier_view));
    std::vector<std::size_t> supported;
    for (std::size_t v = 1; v < sink; ++v) supported.push_back(v);
    const auto support = rfc::numeric::StructuralSupport::sparse(states.size(), supported);
    member(optimization_declaration, "length", optimization_fields.column("consensus.length", 0, states.size(),
        [&](auto v) { return Rational(std::to_string(states[v]->value.length)); }, support));
    member(optimization_declaration, "edit_count", optimization_fields.column("consensus.edits", 0, states.size(),
        [&](auto v) { return Rational(std::to_string(states[v]->value.edits)); }, support));
    member(optimization_declaration, "orientation_tuple", optimization_fields.column("consensus.orientation", 0, states.size(),
        [&](auto v) { return Rational(std::to_string(states[v]->orientation)); }, support));
    N::Array read_coordinates;
    for (std::size_t r = 0; r < reads.size(); ++r) {
        N::Object value;
        member(value, "prefix", optimization_fields.column("consensus.read_prefix", 0, states.size(),
            [&](auto v) { return Rational(std::to_string(states[v]->value.reads[r].prefix)); }, support));
        member(value, "paired_capped", optimization_fields.column("consensus.read_paired", 0, states.size(),
            [&](auto v) { return Rational(std::to_string(states[v]->value.reads[r].paired)); }, support));
        member(value, "phase", optimization_fields.column("consensus.read_phase", 0, states.size(),
            [&](auto v) { return Rational(states[v]->value.reads[r].phase); }, support));
        element(read_coordinates, std::move(value));
    }
    member(optimization_declaration, "read_coordinates", std::move(read_coordinates));
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(resource.bytes)))};
    const auto optimization_record = source_record(std::move(optimization));
    N::Object declaration;
    member(declaration, "optimization_context", child_context(optimization_record));
    member(declaration, "optimization_publication", archive_child(optimization_record, sources));
    // Biological evidence is a primary multiway relation. It is not the
    // binary source used by the separate exact routing objective.
    rfc::core::RelationSignature signature;
    const auto output_role = signature.aspects().add("output"), input_role = signature.aspects().add("original_token"),
               context_role = signature.aspects().add("original_read"), insertion_role = signature.aspects().add("insertion"),
               family_role = signature.aspects().add("family");
    const auto emission_type = signature.add_type("sequence.consensus.joint_emission");
    const auto tail_type = signature.add_type("sequence.consensus.input_tail");
    rfc::core::RelationalComplexBuilder builder(signature); builder.declare_grade({1});
    std::size_t cells = 0, relations = 0;
    const auto vertex = [&]() {
        require(cells < std::numeric_limits<std::uint32_t>::max(), "consensus biological cell representation");
        builder.declare_witness_identity(++cells); return cells;
    };
    for (std::size_t r = 0; r < reads.size(); ++r) (void)vertex();
    std::vector<std::vector<std::size_t>> token_cells(reads.size());
    N::Array original_tokens;
    for (std::size_t r = 0; r < reads.size(); ++r) for (std::size_t q = 0; q < reads[r].bases.size(); ++q) {
        const auto cell = vertex(); token_cells[r].push_back(cell);
        N::Object token_value; member(token_value, "read_record", integer(r)); member(token_value, "offset", integer(q));
        member(token_value, "symbol", std::string(1, reads[r].bases[q]));
        member(token_value, "bytes", original_bases(original, 1, r, q, 1));
        element(original_tokens, N(N::Owned{0, rfc::storage::StableCellKey{{0}, cell}, {},
            std::make_shared<const N>(std::move(token_value))}));
    }
    N::Array families;
    Vector biological_costs;
    std::vector<std::tuple<std::size_t, std::size_t, Rational>> pairs;
    rfc::numeric::Integer family_count{0};
    std::set<std::string> spellings;
    for (const auto &history : histories) {
        const auto family_cell = vertex();
        N::Object family; N::Array tokens, alignments;
        std::vector<N::Array> read_steps(reads.size());
        std::vector<std::optional<std::size_t>> starts(reads.size()), ends(reads.size());
        std::vector<std::size_t> paired_counts(reads.size());
        for (const auto e : history.relations) {
            const auto &step = transitions.at(e);
            if (step.from == 0) continue;
            const auto &from = states.at(step.from)->value;
            const auto output_cell = step.base ? vertex() : family_cell;
            std::vector<rfc::core::BoundaryOccurrenceInput> participants;
            participants.push_back({family_cell, true, false, Rational{1}, {{family_role, true}}, {}});
            if (step.base) participants.push_back({output_cell, true, false, Rational{1}, {{output_role, true}}, {}});
            for (std::size_t r = 0; r < reads.size(); ++r) {
                const auto &action = step.actions[r];
                participants.push_back({r + 1, true, false, Rational{1}, {{context_role, true}}, {}});
                const auto raw_offset = [&](std::size_t oriented) { return history.strands[r] == "+" ? oriented : reads[r].bases.size() - 1 - oriented; };
                const auto append = [&](char op, std::size_t qb, std::size_t qe, std::size_t rb, std::size_t re) {
                    N::Object value; member(value, "operation", std::string(1, op));
                    member(value, "query_begin", integer(qb)); member(value, "query_end", integer(qe));
                    member(value, "output_begin", integer(rb)); member(value, "output_end", integer(re));
                    if (op != 'D') {
                        const auto offset = raw_offset(qb), cell = token_cells[r][offset];
                        member(value, "original_offset", integer(offset)); member(value, "original_token", integer(cell));
                        participants.push_back({cell, true, false, Rational{1}, {{input_role, true}, {insertion_role, op == 'I'}}, {}});
                        if (op != 'I') pairs.emplace_back(output_cell - 1, cell - 1, Rational{1});
                    }
                    element(read_steps[r], std::move(value));
                };
                auto prefix = from.reads[r].prefix;
                if (action.insertions || action.operation == '=' || action.operation == 'X' || action.operation == 'D') {
                    if (!starts[r]) starts[r] = from.length;
                    for (std::size_t i = 0; i < action.insertions; ++i) { append('I', prefix, prefix + 1, from.length, from.length); ++prefix; }
                    if (action.operation != 'O') {
                        const bool pair = action.operation != 'D';
                        append(action.operation, prefix, prefix + pair, from.length, from.length + 1);
                        if (pair) ++paired_counts[r];
                        ends[r] = from.length + 1;
                    } else ends[r] = from.length;
                }
            }
            require(relations < std::numeric_limits<std::uint32_t>::max(), "consensus biological relation representation");
            const auto relation = ++relations;
            builder.add_grade_one(rfc::core::RelationKey{relation}, step.base ? emission_type : tail_type, std::move(participants));
            biological_costs.emplace_back(std::to_string(step.edits));
            if (step.base) {
                N::Object token_value; member(token_value, "symbol", std::string(1, step.base));
                member(token_value, "position", integer(from.length)); member(token_value, "relation", integer(relation));
                element(tokens, N(N::Owned{0, rfc::storage::StableCellKey{{0}, output_cell}, {},
                    std::make_shared<const N>(std::move(token_value))}));
            }
        }
        for (std::size_t r = 0; r < reads.size(); ++r) {
            require(starts[r] && ends[r] && paired_counts[r] >= minimum_paired, "complete consensus read placement witness");
            N::Object value; member(value, "read_record", integer(r)); member(value, "read_id", reads[r].id);
            member(value, "strand", history.strands[r]); member(value, "start", integer(*starts[r]));
            member(value, "end", integer(*ends[r])); member(value, "paired_bases", integer(paired_counts[r]));
            member(value, "steps", std::move(read_steps[r])); element(alignments, std::move(value));
        }
        member(family, "spelling", history.spelling); member(family, "edits", integer(history.edits));
        member(family, "tokens", std::move(tokens)); member(family, "read_alignments", std::move(alignments));
        member(family, "route_index", integer(history.route)); member(family, "family_index", N(family_count));
        element(families, N(N::Owned{0, rfc::storage::StableCellKey{{0}, family_cell}, {},
            std::make_shared<const N>(std::move(family))}));
        ++family_count; spellings.insert(history.spelling);
    }
    auto endpoint = source_endpoint(builder.freeze()); SourceGradeFields fields(endpoint);
    rfc::numeric::ExactMatrix association(cells, cells);
    for (const auto &[output, input, coefficient] : pairs) association.set(output, input, coefficient);
    member(declaration, "paired_occurrence_map", fields.map("consensus.paired_original_occurrences", 1, 1, association));
    member(declaration, "edit_costs", fields.field("consensus.biological_unit_edits", 1, biological_costs));
    member(declaration, "schema", "gadget.sequence.consensus.fields.v1");
    member(declaration, "protocol", import_record_view(request));
    member(declaration, "original_reads", original_records(original, 1, 0));
    member(declaration, "original_tokens", std::move(original_tokens));
    member(declaration, "source_bytes", N(N::ByteSpan{1, 0, resource.bytes.size(), false}));
    N::Object view;
    member(view, "families", std::move(families)); member(view, "family_count", std::move(family_count));
    member(view, "minimum_length", histories.empty() ? N{} : integer(histories.front().spelling.size()));
    member(view, "minimum_edits", histories.empty() ? N{} : integer(histories.front().edits));
    member(view, "minimum_cost", optimum.reachable ? N(optimum.distance) : N{});
    member(view, "identifiable_sequence", spellings.size() == 1); member(view, "complete", true);
    member(view, "status", histories.empty() ? "no_admissible_consensus" : spellings.size() == 1 ? "one_sequence_under_declared_constraints" : "multiple_minimum_sequences");
    member(view, "objective", "minimum_length_then_total_unit_edits");
    member(view, "search_universe", "All complete original-read substring alignments on independently allowed strands, with an explicit total unit-edit budget and minimum paired bases per read. Every output position has an exact paired original token; insertions consume all original tails with explicit cost. All minimizing joint histories remain distinct.");
    member(view, "completeness_argument", "Each emission declares simultaneous original-indexed insertion blocks and paired/deletion/absent actions. Source prefixes advance and every output position consumes an exact original token. Finite output length follows the full original token count. Orientation choices are distinct zero-cost declarations; there are no administrative phase interleavings. Canonical minimum routing with weight budget+1 selects minimum length then total edits; all tight ordered histories are retained.");
    member(declaration, "view", std::move(view));
    return {std::move(endpoint), N(std::move(declaration)), std::move(sources)};
}
Json read_consensus_fields_view(const NativeRecord &record) {
    auto view = render_record_view(record.at("view"));
    const auto optimization = NativeRecord::decode(record.at("optimization_publication").bytes());
    view.set("carrier", sequence::carrier_fields_view(optimization));
    for (std::size_t i = 0; i < view.at("families").size(); ++i) {
        auto family = view.at("families").at(i);
        const auto index = record.at("view").at("families").at(i).at("route_index").integer().get_ui();
        json_object_object_del(family.raw(), "route_index");
        family.set("route", sequence::retained_route_view(optimization.at("routes").at(index)));
    }
    return view;
}
Json read_consensus(const ResourceImages &images, const Json &request) {
    return read_consensus_fields_view(sequence::source_record(read_consensus_fields_source(images, request)));
}
} // namespace gadget::assay_detail
