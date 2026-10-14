#include "sequence_route_witnesses.hpp"
#include "sequence_consensus.hpp"
#include <rfc/local/section.hpp>

namespace gadget::assay_detail {
namespace {
NativeRecord glue_reads(const std::vector<SequenceRecord> &reads, const std::vector<std::size_t> &order,
                const std::vector<std::size_t> &overlaps, const std::string &spelling) {
    using namespace rfc;
    using namespace sequence;
    require(spelling.size() <= std::numeric_limits<std::uint32_t>::max(),
            "glued sequence native cell representation");
    core::RelationSignature signature;
    core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(core::Grade{1});
    for (std::size_t i = 0; i < spelling.size(); ++i)
        builder.declare_witness_identity(i + 1);
    const auto source = builder.freeze();
    const local::RelationalLocalCarrier carrier(source);
    local::BasicSection<char> global;
    algebra::SubcomplexMask global_patch(source);
    auto placements = N::Array{};
    std::size_t start = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto &read = reads[order[i]];
        if (i)
            start += reads[order[i - 1]].bases.size() - overlaps[i - 1];
        local::BasicSection<char> section;
        algebra::SubcomplexMask patch(source);
        for (std::size_t p = 0; p < read.bases.size(); ++p) {
            const core::CellId cell{core::Grade{0}, static_cast<std::uint32_t>(start + p)};
            patch.select(cell);
            section.emplace(cell, read.bases[p]);
        }
        global = local::glue_sections(carrier, global, global_patch, section, patch,
                                      local::OverlapCoverage::RequireComplete);
        for (std::size_t p = 0; p < read.bases.size(); ++p)
            global_patch.select(
                core::CellId{core::Grade{0}, static_cast<std::uint32_t>(start + p)});
        require(local::restrict_section(carrier, global, patch) == section,
                "native read restriction reconstruction");
        auto row = N::Object{};
        member(row, "read_id", N(read.id));
        member(row, "source_read", integer(order[i]));
        member(row, "start", integer(start));
        member(row, "end", integer(start + read.bases.size()));
        member(row, "evidence", import_record_view(read.evidence));
        element(placements, row);
    }
    require(global.size() == spelling.size(),
            "read sections must cover the reconstructed sequence");
    for (std::size_t p = 0; p < spelling.size(); ++p)
        require(global.at(core::CellId{core::Grade{0}, static_cast<std::uint32_t>(p)}) ==
                    spelling[p],
                "native categorical gluing disagrees with spelling");
    auto out = N::Object{};
    member(out, "placements", placements);
    member(out, "complete_overlap_gluing_verified", N(true));
    member(out, "every_read_restriction_verified", N(true));
    member(out, "value_domain",
            N("literal DNA tokens with equality; no numerical encoding"));
    N::Array tokens;
    for (std::size_t i = 0; i < spelling.size(); ++i)
        element(tokens, N(N::Owned{0, storage::StableCellKey{{0}, i + 1}, {},
            std::make_shared<const N>(std::string(1, global.at(core::CellId{{0}, static_cast<std::uint32_t>(i)})))}));
    N::Object declaration; member(declaration, "view", std::move(out)); member(declaration, "tokens", std::move(tokens));
    return NativeRecord::create_from_endpoint(N(std::move(declaration)), new_native_record_namespace(), {}, source_endpoint(source));
}
} // namespace
AssayFieldsSource read_reconstruction_fields_source(const ResourceImage &resource, const Json &request) {
    using namespace sequence;
    members(request,
            {"schema", "operation", "resource", "minimum_overlap", "objective", "evidence"});
    require(text(request.at("objective")) == "minimum_suffix_prefix_read_cover",
            "read reconstruction objective");
    (void)text(request.at("evidence"), source_byte_limit());
    const auto minimum =
        count(request.at("minimum_overlap"), std::numeric_limits<std::size_t>::max(),
              "minimum overlap representation");
    require(minimum > 0, "read overlap must be positive");
    const SequenceSource original_source(resource);
    const auto reads = records(original_source);
    require(!reads.empty(), "read cover requires source reads");
    const auto n = reads.size();
    std::vector<std::vector<std::size_t>> overlap(n, std::vector<std::size_t>(n));
    auto overlap_records = N::Array{}, source_reads = N::Array{};
    Routes overlap_routes(n, "sequence.read_overlap.compatibility");
    for (std::size_t a = 0; a < n; ++a) {
        require(reads[a].bases.size() >= minimum, "read is shorter than declared overlap");
        auto row = N::Object{};
        member(row, "read_id", N(reads[a].id));
        member(row, "bases", N(reads[a].bases));
        member(row, "evidence", import_record_view(reads[a].evidence));
        element(source_reads, row);
        for (std::size_t b = 0; b < n; ++b) {
            if (a == b)
                continue;
            auto lengths = N::Array{};
            for (std::size_t length = minimum;
                 length <= std::min(reads[a].bases.size(), reads[b].bases.size()); ++length)
                if (reads[a].bases.compare(reads[a].bases.size() - length, length, reads[b].bases,
                                           0, length) == 0) {
                    overlap[a][b] = length;
                    element(lengths, integer(length));
                }
            if (lengths.elements.size()) {
                auto r = N::Object{};
                member(r, "left_read", integer(a));
                member(r, "right_read", integer(b));
                member(r, "admissible_lengths", lengths);
                member(r, "retained_length", integer(overlap[a][b]));
                element(overlap_records, r);
                overlap_routes.add(a, b,
                                   Rational{std::to_string(reads[b].bases.size() - overlap[a][b])});
            }
        }
    }
    const auto overlap_source = overlap_routes.freeze();
    std::vector<RawByteSource> sources{RawByteSource(1, std::as_bytes(std::span(resource.bytes)))};
    auto native_overlap = route_source(overlap_source, overlap_routes);
    SourceGradeFields overlap_fields(native_overlap.endpoint);
    N::Array reachability_routes, gluing_publications, gluing_contexts;
    const auto overlap_graph =
        rfc::view::project_pairwise_directed(overlap_source, overlap_routes.costs());
    for (std::size_t a = 0; a < n; ++a)
        require(overlap_graph.object_index(a + 1) == a, "native source read coordinate order");
    auto start_receipts = N::Array{};
    std::vector<std::size_t> starts;
    rfc::numeric::Integer reachability_calls{0};
    for (std::size_t a = 0; a < n; ++a) {
        std::vector<bool> reached(n, false);
        reached[a] = true;
        auto receipt = N::Object{}; auto histories = N::Array{};
        member(receipt, "source_read", integer(a));
        bool eligible = true;
        for (std::size_t b = n; b-- > 0;) {
            if (reached[b])
                continue;
            ++reachability_calls;
            const auto route =
                rfc::path::bfs_min_hops(overlap_graph, overlap_graph.object_index(a + 1),
                                        overlap_graph.object_index(b + 1));
            if (!route.reachable) {
                eligible = false;
                member(receipt, "unreachable_source_read", integer(b));
                break;
            }
            auto history = N::Object{}; auto relations = N::Array{};
            for (auto e : route.projected_edges) {
                const auto &edge = overlap_graph.edge(e);
                reached[edge.target_object] = true;
                element(relations, integer(edge.source_relation));
            }
            member(history, "target_read", integer(b));
            member(history, "ordered_relations", relations);
            std::vector<std::uint32_t> native_relations;
            for (const auto edge : route.projected_edges) native_relations.push_back(overlap_graph.edge(edge).source_relation);
            Rational native_cost{0};
            for (const auto e : native_relations) native_cost += overlap_routes.costs().at(e);
            element(reachability_routes, retain_route(overlap_fields, overlap_source, overlap_routes, native_relations, a, b, native_cost));
            element(histories, history);
        }
        member(receipt, "eligible_start", N(eligible));
        member(receipt, "native_reachability_histories", histories);
        element(start_receipts, receipt);
        if (eligible)
            starts.push_back(a);
    }
    using Bits = rfc::machine::PackedBits;
    struct State {
        Bits seen;
        std::size_t last;
    };
    using StateKey = std::pair<std::vector<Bits::word_type>, std::size_t>;
    struct Transition {
        std::size_t next, read, relation;
    };
    std::map<StateKey, std::size_t> state_index;
    std::vector<State> states;
    std::vector<std::vector<Transition>> next;
    const auto admit = [&](Bits seen, std::size_t last) {
        StateKey key{seen.words(), last};
        if (const auto it = state_index.find(key); it != state_index.end())
            return it->second;
        require(rfc::machine::checked_size_add(states.size(), 3) <=
                    std::numeric_limits<std::uint32_t>::max(),
                "read-cover native state coordinate representation");
        const auto index = states.size();
        state_index.emplace(std::move(key), index);
        states.push_back({std::move(seen), last});
        next.emplace_back();
        return index;
    };
    std::vector<std::pair<std::size_t, std::size_t>> entries;
    for (auto a : starts) {
        Bits seen(n, false, Bits::words_for_size(n));
        seen.set(a);
        entries.emplace_back(a, admit(std::move(seen), a));
    }
    for (std::size_t v = 0; v < states.size(); ++v) {
        const auto seen = states[v].seen;
        const auto a = states[v].last;
        for (std::size_t b = 0; b < n; ++b)
            if (!seen.get(b) && overlap[a][b]) {
                auto added = seen;
                added.set(b);
                const auto target = admit(std::move(added), b);
                next[v].push_back({target, b, 0});
            }
    }
    const auto sink = rfc::machine::checked_size_add(states.size(), 1),
               vertices = rfc::machine::checked_size_add(sink, 1);
    Routes routes(vertices, "sequence.read_cover.step");
    std::vector<std::pair<std::size_t, std::size_t>> entry_relations;
    std::vector<std::optional<std::size_t>> exits(states.size());
    for (const auto &[a, v] : entries)
        entry_relations.emplace_back(
            v, routes.add(0, v + 1, Rational{std::to_string(reads[a].bases.size())}));
    for (std::size_t v = 0; v < states.size(); ++v) {
        if (states[v].seen.count() == n)
            exits[v] = routes.add(v + 1, sink, Rational{0});
        for (auto &link : next[v])
            link.relation =
                routes.add(v + 1, link.next + 1,
                           Rational{std::to_string(reads[link.read].bases.size() -
                                                   overlap[states[v].last][link.read])});
    }
    const auto source = routes.freeze();
    const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
    const auto optimum =
        rfc::path::dijkstra_exact(graph, graph.object_index(1), graph.object_index(sink + 1));
    std::optional<CompleteRouteWitnesses> witnesses;
    if (optimum.reachable) {
        witnesses.emplace(graph, graph.object_index(1), graph.object_index(sink + 1), true,
                          optimum);
        witnesses->verify(rfc::field::boundary_matrix(source, rfc::core::Grade{1}), routes.costs(),
                          0, sink);
    }
    auto native = route_source(source, routes);
    SourceGradeFields fields(native.endpoint);
    N::Object out, c;
    N::Array state_records, layouts, sequences, native_routes, native_histories;
    for (std::size_t v = 0; v < states.size(); ++v) {
        const auto &state = states[v];
        auto row = N::Object{}; auto members = N::Array{};
        rfc::numeric::Integer mask{0};
        for (std::size_t w = 0; w < state.seen.word_count(); ++w)
            mask += rfc::numeric::Integer(std::to_string(state.seen.word(w)))
                    << rfc::machine::checked_size_mul(w, Bits::bits_per_word);
        state.seen.for_each_set([&](std::size_t a) { element(members, integer(a)); });
        member(row, "vertex", integer(v + 1));
        member(row, "source_read_mask", N(mask));
        member(row, "source_reads", members);
        member(row, "last_read", integer(state.last));
        if (witnesses) {
            member(row, "minimum_prefix_length", N(witnesses->prefix[v + 1]));
            member(row, "minimum_suffix_length", witnesses->suffix_reachable[v + 1]
                                                 ? N(witnesses->suffix[v + 1])
                                                 : N{});
        }
        element(state_records, row);
    }
    member(c, "states", state_records);
    member(c, "source_vertex", integer(0));
    member(c, "sink_vertex", integer(sink));
    member(c, "state_domain", N("Reachable prefix states from every start not excluded by "
                                       "native directed reachability; exact subset and last-read "
                                       "coordinates, with no source identity merged"));
    member(out, "source_reads", source_reads);
    member(out, "overlaps", overlap_records);
    member(out, "start_reachability", start_receipts);
    member(out, "native_reachability_calls", N(reachability_calls));
    member(out, "minimum_length",
            optimum.reachable ? N(optimum.distance) : N{});
    if (witnesses) {
        auto absent = N::Array{};
        for (auto v : witnesses->unreachable_suffix_vertices)
            element(absent, integer(v));
        for (const auto &history : witnesses->histories) {
            N::Object row; member(row, "kind", history.kind); member(row, "from_vertex", integer(history.from));
            member(row, "to_vertex", integer(history.to)); member(row, "minimum_cost", history.route.distance);
            N::Array ordered_history; for (auto e : history.route.relations) element(ordered_history, integer(e));
            member(row, "ordered_relations", std::move(ordered_history)); element(native_histories, std::move(row));
            element(native_routes, retain_route(fields, source, routes, history.route.relations, history.from, history.to, history.route.distance));
        }
        member(out, "route_witnesses", native_histories);
        member(out, "unreachable_suffix_vertices", absent);
        member(out, "native_route_calls", integer(witnesses->histories.size() +
                                                  witnesses->unreachable_suffix_vertices.size()));
        member(out, "prefix_dual_verified", N(true));
        member(out, "viable_suffix_dual_verified", N(true));
    }
    std::map<std::string, std::vector<std::size_t>> spellings;
    std::vector<rfc::path::RealizedPath> representatives;
    std::vector<std::size_t> representative_layouts, order, overlaps;
    std::vector<std::uint32_t> relations;
    std::string spelling;
    rfc::numeric::Integer visits{0};
    const auto emit = [&](std::size_t v) {
        relations.push_back(static_cast<std::uint32_t>(*exits[v]));
        auto layout = N::Object{}; auto ids = N::Array{}, lengths = N::Array{};
        for (auto a : order)
            element(ids, integer(a));
        for (auto length : overlaps)
            element(lengths, integer(length));
        member(layout, "read_order", ids);
        member(layout, "overlap_lengths", lengths);
        member(layout, "spelling", N(spelling));
        const bool new_spelling = !spellings.contains(spelling);
        rfc::path::RealizedPath realized;
        member(layout, "route_index", integer(native_routes.elements.size()));
        element(native_routes, retain_route(fields, source, routes, relations, 0, sink, optimum.distance,
                                           new_spelling ? &realized : nullptr));
        member(layout, "gluing_index", integer(gluing_publications.elements.size()));
        const auto glued = glue_reads(reads, order, overlaps, spelling);
        element(gluing_publications, archive_child(glued, sources));
        N::Object context;
        member(context, "publication", gluing_publications.elements.back());
        member(context, "native_context", child_context(glued));
        N::Array inputs;
        for (const auto a : order) {
            N::Object input;
            member(input, "source_read", integer(a));
            member(input, "bases", original_bases(original_source, 1, a, 0, original_source.length(a)));
            element(inputs, std::move(input));
        }
        member(context, "inputs", std::move(inputs));
        element(gluing_contexts, std::move(context));
        if (new_spelling) {
            representatives.push_back(std::move(realized));
            representative_layouts.push_back(layouts.elements.size());
        }
        spellings[spelling].push_back(layouts.elements.size());
        element(layouts, layout);
        relations.pop_back();
    };
    if (witnesses) {
        struct Frame {
            std::size_t state, next_link, previous_length;
            bool entered;
        };
        for (const auto &[v, relation] : entry_relations) {
            const auto a = states[v].last;
            if (!witnesses->suffix_reachable[v + 1] ||
                Rational{std::to_string(reads[a].bases.size())} + witnesses->suffix[v + 1] !=
                    optimum.distance)
                continue;
            order.push_back(a);
            relations.push_back(static_cast<std::uint32_t>(relation));
            spelling = reads[a].bases;
            std::vector<Frame> stack{{v, 0, 0, false}};
            while (!stack.empty()) {
                auto &frame = stack.back();
                if (!frame.entered) {
                    frame.entered = true;
                    ++visits;
                    if (exits[frame.state])
                        emit(frame.state);
                }
                if (frame.next_link < next[frame.state].size()) {
                    const auto link = next[frame.state][frame.next_link++];
                    const auto length = overlap[states[frame.state].last][link.read],
                               extra = reads[link.read].bases.size() - length;
                    if (!witnesses->suffix_reachable[link.next + 1] ||
                        Rational{std::to_string(
                            rfc::machine::checked_size_add(spelling.size(), extra))} +
                                witnesses->suffix[link.next + 1] !=
                            optimum.distance)
                        continue;
                    const auto previous_length = spelling.size();
                    spelling += reads[link.read].bases.substr(length);
                    order.push_back(link.read);
                    overlaps.push_back(length);
                    relations.push_back(static_cast<std::uint32_t>(link.relation));
                    stack.push_back({link.next, 0, previous_length, false});
                } else {
                    const auto previous_length = frame.previous_length;
                    stack.pop_back();
                    if (!stack.empty()) {
                        spelling.resize(previous_length);
                        order.pop_back();
                        overlaps.pop_back();
                        relations.pop_back();
                    }
                }
            }
            order.pop_back();
            relations.pop_back();
        }
    }
    require(!optimum.reachable || layouts.elements.size() > 0,
            "native minimum must have an enumerated cover");
    for (const auto &[text, indices] : spellings) {
        auto row = N::Object{}; auto ids = N::Array{};
        for (auto i : indices)
            element(ids, integer(i));
        member(row, "sequence", N(text));
        member(row, "layouts", ids);
        element(sequences, row);
    }
    N::Object comparison;
    N::Array representative_ids;
    const auto representative_carrier = fields.explicit_carrier(representatives.size());
    rfc::numeric::ExactMatrix G(representatives.size(), representatives.size());
    if (!representatives.empty()) {
        G = rfc::path::inherited_gram(
            representatives, rfc::field::ExactMetric::identity(routes.steps.size()));
        for (std::size_t a = 0; a < representatives.size(); ++a) {
            element(representative_ids, integer(representative_layouts[a]));
        }
    }
    member(comparison, "representative_layouts", representative_ids);
    member(comparison, "all_endpoint_columns_equal", N(true));
    member(comparison, "metric",
                   N("identity on the declared read-cover transition carrier"));
    member(out, "path_comparison", comparison);
    member(out, "layouts", layouts);
    member(out, "sequences", sequences);
    member(out, "identifiable_sequence", N(spellings.size() == 1));
    member(out, "status",
            N(!optimum.reachable      ? "no_complete_read_cover"
                         : spellings.size() == 1 ? "one_sequence_under_declared_constraints"
                                                 : "multiple_minimum_sequences"));
    member(out, "enumeration_visits", N(visits));
    member(out, "complete", N(true));
    member(out, "dominance_argument",
            N(
                "For fixed visited source reads, last read and next read, the longest exact "
                "overlap reaches the same next state at strictly smaller added length than any "
                "shorter overlap. Every valid full cover starts at a read from which all other "
                "source reads are reachable in the complete overlap graph. Native reachability "
                "excludes only impossible starts; exact native suffix witnesses prune only "
                "nonminimizing prefixes. All original overlap alternatives remain retained."));
    member(out,
        "search_universe",
        N(
            "All suffix-prefix covers using every distinct source read exactly once, in its "
            "retained orientation, with at least the declared overlap. All minimum-length layouts "
            "and distinct spellings are retained; source duplicates remain distinct reads. "
            "Interior containment and orientation changes require a different declaration."));
    auto &declaration = std::get<N::Object>(native.declaration.data);
    auto &overlap_declaration = std::get<N::Object>(native_overlap.declaration.data);
    member(overlap_declaration, "reachability_routes", std::move(reachability_routes));
    const auto overlap_record = source_record(std::move(native_overlap));
    member(declaration, "overlap_context", child_context(overlap_record));
    member(declaration, "overlap_publication", archive_child(overlap_record, sources));
    member(declaration, "gluing_publications", std::move(gluing_publications));
    member(declaration, "gluing_contexts", std::move(gluing_contexts));
    member(declaration, "routes", std::move(native_routes));
    member(declaration, "view", std::move(out)); member(declaration, "carrier_view", std::move(c));
    member(declaration, "schema", std::string("gadget.sequence.read_cover.fields.v1"));
    member(declaration, "protocol", import_record_view(request));
    member(declaration, "original_reads", original_records(original_source, 1));
    member(declaration, "source_bytes", N(N::ByteSpan{1, 0, resource.bytes.size(), false}));
    member(declaration, "representative_gram", fields.map("representative.field_gram", representative_carrier,
        representative_carrier, G, rfc::storage::TemporalTensorRole::CrossForm));
    if (witnesses) {
        member(declaration, "prefix", fields.field("route.prefix", 0, witnesses->prefix));
        std::vector<std::size_t> viable;
        for (std::size_t i = 0; i < witnesses->suffix.size(); ++i) if (witnesses->suffix_reachable[i]) viable.push_back(i);
        member(declaration, "suffix", fields.column("route.suffix", 0, witnesses->suffix.size(),
            [&](auto i) { return witnesses->suffix[i]; }, rfc::numeric::StructuralSupport::sparse(witnesses->suffix.size(), viable)));
        member(declaration, "optimal_support", fields.column("route.optimal", 1, routes.steps.size(), [&](auto i) {
            const auto &step = routes.steps[i];
            return Rational(witnesses->suffix_reachable[step.to] &&
                witnesses->prefix[step.from] + step.cost + witnesses->suffix[step.to] == optimum.distance ? 1 : 0);
        }));
    }
    native.sources = std::move(sources);
    return native;
}
Json read_reconstruction_fields_view(const NativeRecord &record) {
    using namespace sequence;
    auto out = render_record_view(record.at("view"));
    auto main = carrier_fields_view(record, false);
    const auto extra = render_record_view(record.at("carrier_view"));
    for (const auto &key : record.at("carrier_view").names()) main.set(key, extra.at(key));
    out.set("carrier", main);
    out.set("overlap_carrier", carrier_fields_view(NativeRecord::decode(record.at("overlap_publication").bytes())));
    for (std::size_t i = 0; i < out.at("layouts").size(); ++i) {
        auto layout = out.at("layouts").at(i);
        const auto native = record.at("view").at("layouts").at(i);
        const auto route = native.at("route_index").integer().get_ui(), glue = native.at("gluing_index").integer().get_ui();
        json_object_object_del(layout.raw(), "route_index"); json_object_object_del(layout.raw(), "gluing_index");
        layout.set("route", retained_route_view(record.at("routes").at(route)));
        layout.set("gluing", render_record_view(NativeRecord::decode(record.at("gluing_publications").at(glue).bytes()).at("view")));
    }
    const auto &gram = retained_tensor(record, record.at("representative_gram")).tensor;
    auto rows = Json::array(); const auto n = gram.domain().extents().at(0);
    for (std::size_t i = 0; i < n; ++i) { auto row = Json::array();
        for (std::size_t j = 0; j < n; ++j) row.push(Json::string(gram.state(i * n + j).value().get_str()));
        rows.push(row);
    }
    out.at("path_comparison").set("shared_field_gram", rows);
    return out;
}
Json read_reconstruction(const ResourceImage &resource, const Json &request) {
    return read_reconstruction_fields_view(sequence::source_record(read_reconstruction_fields_source(resource, request)));
}
} // namespace gadget::assay_detail
