#include "sequence_route_witnesses.hpp"
#include <rfc/local/section.hpp>

namespace gadget::assay_detail {
namespace {
Json glue_reads(const std::vector<SequenceRecord> &reads, const std::vector<std::size_t> &order,
                const std::vector<std::size_t> &overlaps, const std::string &spelling) {
    using namespace rfc;
    require(spelling.size() <= std::numeric_limits<std::uint32_t>::max(),
            "glued sequence native cell representation");
    core::RelationSignature signature;
    core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < spelling.size(); ++i)
        builder.declare_witness_identity(i + 1);
    const auto source = builder.freeze();
    const local::RelationalLocalCarrier carrier(source);
    local::BasicSection<char> global;
    algebra::SubcomplexMask global_patch(source);
    auto placements = Json::array();
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
        auto row = Json::object();
        row.set("read_id", Json::string(read.id));
        row.set("source_read", cardinality(order[i]));
        row.set("start", cardinality(start));
        row.set("end", cardinality(start + read.bases.size()));
        row.set("evidence", read.evidence);
        placements.push(row);
    }
    require(global.size() == spelling.size(),
            "read sections must cover the reconstructed sequence");
    for (std::size_t p = 0; p < spelling.size(); ++p)
        require(global.at(core::CellId{core::Grade{0}, static_cast<std::uint32_t>(p)}) ==
                    spelling[p],
                "native categorical gluing disagrees with spelling");
    auto out = Json::object();
    out.set("placements", placements);
    out.set("complete_overlap_gluing_verified", Json::boolean(true));
    out.set("every_read_restriction_verified", Json::boolean(true));
    out.set("value_domain",
            Json::string("literal DNA tokens with equality; no numerical encoding"));
    return out;
}
} // namespace
Json read_reconstruction(const ResourceImage &resource, const Json &request) {
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
    const auto reads = sequence_records(resource);
    require(!reads.empty(), "read cover requires source reads");
    const auto n = reads.size();
    std::vector<std::vector<std::size_t>> overlap(n, std::vector<std::size_t>(n));
    auto overlap_records = Json::array(), source_reads = Json::array();
    Routes overlap_routes(n, "sequence.read_overlap.compatibility");
    for (std::size_t a = 0; a < n; ++a) {
        require(reads[a].bases.size() >= minimum, "read is shorter than declared overlap");
        auto row = Json::object();
        row.set("read_id", Json::string(reads[a].id));
        row.set("bases", Json::string(reads[a].bases));
        row.set("evidence", reads[a].evidence);
        source_reads.push(row);
        for (std::size_t b = 0; b < n; ++b) {
            if (a == b)
                continue;
            auto lengths = Json::array();
            for (std::size_t length = minimum;
                 length <= std::min(reads[a].bases.size(), reads[b].bases.size()); ++length)
                if (reads[a].bases.compare(reads[a].bases.size() - length, length, reads[b].bases,
                                           0, length) == 0) {
                    overlap[a][b] = length;
                    lengths.push(cardinality(length));
                }
            if (lengths.size()) {
                auto r = Json::object();
                r.set("left_read", cardinality(a));
                r.set("right_read", cardinality(b));
                r.set("admissible_lengths", lengths);
                r.set("retained_length", cardinality(overlap[a][b]));
                overlap_records.push(r);
                overlap_routes.add(a, b,
                                   Rational{std::to_string(reads[b].bases.size() - overlap[a][b])});
            }
        }
    }
    const auto overlap_source = overlap_routes.freeze();
    const auto overlap_graph =
        rfc::view::project_pairwise_directed(overlap_source, overlap_routes.costs());
    for (std::size_t a = 0; a < n; ++a)
        require(overlap_graph.object_index(a + 1) == a, "native source read coordinate order");
    auto start_receipts = Json::array();
    std::vector<std::size_t> starts;
    rfc::numeric::Integer reachability_calls{0};
    for (std::size_t a = 0; a < n; ++a) {
        std::vector<bool> reached(n, false);
        reached[a] = true;
        auto receipt = Json::object(), histories = Json::array();
        receipt.set("source_read", cardinality(a));
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
                receipt.set("unreachable_source_read", cardinality(b));
                break;
            }
            auto history = Json::object(), relations = Json::array();
            for (auto e : route.projected_edges) {
                const auto &edge = overlap_graph.edge(e);
                reached[edge.target_object] = true;
                relations.push(cardinality(edge.source_relation));
            }
            history.set("target_read", cardinality(b));
            history.set("ordered_relations", relations);
            histories.push(history);
        }
        receipt.set("eligible_start", Json::boolean(eligible));
        receipt.set("native_reachability_histories", histories);
        start_receipts.push(receipt);
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
    auto out = Json::object(), c = carrier(source, routes), state_records = Json::array(),
         layouts = Json::array(), sequences = Json::array();
    for (std::size_t v = 0; v < states.size(); ++v) {
        const auto &state = states[v];
        auto row = Json::object(), members = Json::array();
        rfc::numeric::Integer mask{0};
        for (std::size_t w = 0; w < state.seen.word_count(); ++w)
            mask += rfc::numeric::Integer(std::to_string(state.seen.word(w)))
                    << rfc::machine::checked_size_mul(w, Bits::bits_per_word);
        state.seen.for_each_set([&](std::size_t a) { members.push(cardinality(a)); });
        row.set("vertex", cardinality(v + 1));
        row.set("source_read_mask", exact_count(mask));
        row.set("source_reads", members);
        row.set("last_read", cardinality(state.last));
        if (witnesses) {
            row.set("minimum_prefix_length", Json::string(witnesses->prefix[v + 1].get_str()));
            row.set("minimum_suffix_length", witnesses->suffix_reachable[v + 1]
                                                 ? Json::string(witnesses->suffix[v + 1].get_str())
                                                 : Json{});
        }
        state_records.push(row);
    }
    c.set("states", state_records);
    c.set("source_vertex", cardinality(0));
    c.set("sink_vertex", cardinality(sink));
    c.set("state_domain", Json::string("Reachable prefix states from every start not excluded by "
                                       "native directed reachability; exact subset and last-read "
                                       "coordinates, with no source identity merged"));
    out.set("carrier", c);
    out.set("source_reads", source_reads);
    out.set("overlaps", overlap_records);
    out.set("overlap_carrier", carrier(overlap_source, overlap_routes));
    out.set("start_reachability", start_receipts);
    out.set("native_reachability_calls", exact_count(reachability_calls));
    out.set("minimum_length",
            optimum.reachable ? Json::string(optimum.distance.get_str()) : Json{});
    if (witnesses) {
        auto absent = Json::array();
        for (auto v : witnesses->unreachable_suffix_vertices)
            absent.push(cardinality(v));
        out.set("route_witnesses", witnesses->histories);
        out.set("unreachable_suffix_vertices", absent);
        out.set("native_route_calls", cardinality(witnesses->histories.size() +
                                                  witnesses->unreachable_suffix_vertices.size()));
        out.set("prefix_dual_verified", Json::boolean(true));
        out.set("viable_suffix_dual_verified", Json::boolean(true));
    }
    std::map<std::string, std::vector<std::size_t>> spellings;
    std::vector<rfc::path::RealizedPath> representatives;
    std::vector<std::size_t> representative_layouts, order, overlaps;
    std::vector<std::uint32_t> relations;
    std::string spelling;
    rfc::numeric::Integer visits{0};
    const auto emit = [&](std::size_t v) {
        relations.push_back(static_cast<std::uint32_t>(*exits[v]));
        auto layout = Json::object(), ids = Json::array(), lengths = Json::array();
        for (auto a : order)
            ids.push(cardinality(a));
        for (auto length : overlaps)
            lengths.push(cardinality(length));
        layout.set("read_order", ids);
        layout.set("overlap_lengths", lengths);
        layout.set("spelling", Json::string(spelling));
        const bool new_spelling = !spellings.contains(spelling);
        rfc::path::RealizedPath realized;
        layout.set("route", realized_receipt(source, routes, relations, 0, sink, optimum.distance,
                                             new_spelling ? &realized : nullptr));
        layout.set("gluing", glue_reads(reads, order, overlaps, spelling));
        if (new_spelling) {
            representatives.push_back(std::move(realized));
            representative_layouts.push_back(layouts.size());
        }
        spellings[spelling].push_back(layouts.size());
        layouts.push(layout);
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
    require(!optimum.reachable || layouts.size() > 0,
            "native minimum must have an enumerated cover");
    for (const auto &[text, indices] : spellings) {
        auto row = Json::object(), ids = Json::array();
        for (auto i : indices)
            ids.push(cardinality(i));
        row.set("sequence", Json::string(text));
        row.set("layouts", ids);
        sequences.push(row);
    }
    auto comparison = Json::object(), gram = Json::array(), representative_ids = Json::array();
    if (!representatives.empty()) {
        const auto G = rfc::path::inherited_gram(
            representatives, rfc::field::ExactMetric::identity(routes.steps.size()));
        for (std::size_t a = 0; a < representatives.size(); ++a) {
            auto row = Json::array();
            representative_ids.push(cardinality(representative_layouts[a]));
            for (std::size_t b = 0; b < representatives.size(); ++b)
                row.push(Json::string(G.get(a, b).get_str()));
            gram.push(row);
        }
    }
    comparison.set("representative_layouts", representative_ids);
    comparison.set("shared_field_gram", gram);
    comparison.set("all_endpoint_columns_equal", Json::boolean(true));
    comparison.set("metric",
                   Json::string("identity on the declared read-cover transition carrier"));
    out.set("path_comparison", comparison);
    out.set("layouts", layouts);
    out.set("sequences", sequences);
    out.set("identifiable_sequence", Json::boolean(spellings.size() == 1));
    out.set("status",
            Json::string(!optimum.reachable      ? "no_complete_read_cover"
                         : spellings.size() == 1 ? "one_sequence_under_declared_constraints"
                                                 : "multiple_minimum_sequences"));
    out.set("enumeration_visits", exact_count(visits));
    out.set("complete", Json::boolean(true));
    out.set("dominance_argument",
            Json::string(
                "For fixed visited source reads, last read and next read, the longest exact "
                "overlap reaches the same next state at strictly smaller added length than any "
                "shorter overlap. Every valid full cover starts at a read from which all other "
                "source reads are reachable in the complete overlap graph. Native reachability "
                "excludes only impossible starts; exact native suffix witnesses prune only "
                "nonminimizing prefixes. All original overlap alternatives remain retained."));
    out.set(
        "search_universe",
        Json::string(
            "All suffix-prefix covers using every distinct source read exactly once, in its "
            "retained orientation, with at least the declared overlap. All minimum-length layouts "
            "and distinct spellings are retained; source duplicates remain distinct reads. "
            "Interior containment and orientation changes require a different declaration."));
    return out;
}
} // namespace gadget::assay_detail
