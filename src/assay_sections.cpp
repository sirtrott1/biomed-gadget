#include "assay_field_detail.hpp"
#include "scientific_projection.hpp"
#include <gadget/resources.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/local/cellular_sheaf.hpp>
#include <rfc/runtime/cpu_tiles.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::assay_detail {
namespace {
using M = rfc::numeric::ExactMatrix;
using Metric = rfc::field::ExactMetric;
Json sparse_matrix(const M &matrix) {
    auto out = Json::object(), entries = Json::array();
    for (const auto &[ij, value] : matrix.entries()) {
        auto e = Json::object();
        e.set("row", cardinality(ij.first));
        e.set("column", cardinality(ij.second));
        e.set("value", Json::string(value.get_str()));
        entries.push(e);
    }
    out.set("rows", cardinality(matrix.rows()));
    out.set("columns", cardinality(matrix.cols()));
    out.set("entries", entries);
    return out;
}
std::string label(const Json &value) { return text(value, source_byte_limit()); }
struct View {
    std::string id;
    std::vector<AssayField> fields;
    Vector diagonal;
    Json declaration;
};
struct Flag {
    std::size_t view{}, relation{};
    M map;
};
struct Case {
    Vector input, residual, compatible, defect;
    std::vector<unsigned char> support, assessed;
    Rational quadrance{}, orthogonality{};
    std::size_t group{};
};
struct Group {
    std::vector<std::size_t> coordinates, rows, cases;
    M inclusion, restriction, action;
    std::unique_ptr<PreparedProjection> projection;
};
} // namespace
Json assay_sections(const AssayWorkspace &workspace, const Json &request,
                    std::shared_ptr<CpuResources> resources) {
    members(request,
            {"schema", "operation", "views", "sections", "interfaces", "metric", "support_policy"},
            {"workers"});
    require(request.at("support_policy").text() == "require_all_nonzero_occurrences",
            "assay section support policy");
    const auto sections = request.at("sections");
    require(sections.size() > 0, "assay sections must be nonempty");
    std::map<std::string, std::size_t> section_ids, view_ids;
    for (std::size_t i = 0; i < sections.size(); ++i)
        require(section_ids.emplace(label(sections.at(i)), i).second, "duplicate section identity");
    const auto penalty = request.at("metric");
    members(penalty, {"quantity", "units", "method", "evidence"});
    for (const auto *name : {"quantity", "units", "method", "evidence"})
        (void)label(penalty.at(name));
    const auto wanted = request.has("workers")
                            ? count(request.at("workers"), std::numeric_limits<std::size_t>::max(),
                                    "assay section workers")
                            : 1;
    require(wanted > 0, "assay section workers must be positive");
    auto reservation = resources->reserve(wanted);
    const auto workers = reservation.workers();
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("declared_assay_interface",
                                         rfc::core::BoundaryPolicy::DeclaredExact, 1, 0);
    rfc::core::RelationalComplexBuilder builder(signature);
    std::vector<View> views;
    std::vector<Case> cases(sections.size());
    Vector diagonal;
    const auto records = request.at("views");
    require(records.size() > 0, "assay section views must be nonempty");
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto record = records.at(i);
        members(record, {"id", "assay", "partition", "fields", "metric"});
        View view;
        view.id = label(record.at("id"));
        require(view_ids.emplace(view.id, i).second, "duplicate assay view");
        view.declaration = record;
        const auto owner =
            workspace.study(label(record.at("assay")),
                            count(record.at("partition"), std::numeric_limits<std::size_t>::max(),
                                  "assay view partition"));
        const auto fields = record.at("fields");
        require(fields.size() == sections.size(), "every section requires one field per view");
        std::vector<std::optional<std::string>> ordered(sections.size());
        for (std::size_t j = 0; j < fields.size(); ++j) {
            const auto field = fields.at(j);
            members(field, {"section", "field"});
            const auto id = label(field.at("section"));
            require(section_ids.contains(id) && !ordered[section_ids.at(id)],
                    "unknown or repeated view section");
            ordered[section_ids.at(id)] = label(field.at("field"));
        }
        std::vector<std::string> names;
        for (const auto &field : ordered) {
            require(field.has_value(), "section field not declared");
            names.push_back(*field);
        }
        view.fields = load_assay_fields(owner, names);
        require(!view.fields.front().values.empty(), "empty assay section stalk");
        const auto &first = view.fields.front();
        for (const auto &field : view.fields) {
            require(field.grade == first.grade,
                    "section fields must share the actual view carrier");
            for (const auto *name : {"quantity", "units", "channel"})
                require(field.meaning.at(name).text() == first.meaning.at(name).text(),
                        "section field meanings differ within a view");
        }
        const auto form = record.at("metric");
        members(form, {"diagonal", "domain_quantity", "domain_units", "quadrance_quantity",
                       "quadrance_units", "method", "evidence"});
        require(form.at("domain_quantity").text() == first.meaning.at("quantity").text() &&
                    form.at("domain_units").text() == first.meaning.at("units").text() &&
                    form.at("quadrance_quantity").text() == penalty.at("quantity").text() &&
                    form.at("quadrance_units").text() == penalty.at("units").text(),
                "assay section metric meanings disagree");
        (void)label(form.at("method"));
        (void)label(form.at("evidence"));
        require(form.at("diagonal").size() == first.values.size(), "assay view metric dimension");
        for (std::size_t j = 0; j < first.values.size(); ++j) {
            const auto value = exact_rational(label(form.at("diagonal").at(j)));
            require(value > 0, "assay view metric must be positive");
            view.diagonal.push_back(value);
            diagonal.push_back(value);
        }
        for (std::size_t s = 0; s < cases.size(); ++s) {
            const auto &field = view.fields[s];
            cases[s].input.insert(cases[s].input.end(), field.values.begin(), field.values.end());
            cases[s].support.insert(cases[s].support.end(), field.supported.begin(),
                                    field.supported.end());
        }
        builder.declare_witness_identity(i);
        views.push_back(std::move(view));
    }
    const auto interfaces = request.at("interfaces");
    require(interfaces.size() > 0, "assay sections require declared interfaces");
    std::set<std::string> interface_ids;
    std::vector<std::size_t> dimensions;
    std::vector<Flag> flags;
    for (std::size_t i = 0; i < interfaces.size(); ++i) {
        const auto interface = interfaces.at(i);
        members(interface, {"id", "axes", "quantity", "units", "channel", "method", "evidence",
                            "participants"});
        require(interface_ids.insert(label(interface.at("id"))).second, "duplicate interface");
        for (const auto *name : {"quantity", "units", "channel", "method", "evidence"})
            (void)label(interface.at(name));
        const auto channel = interface.at("channel").text();
        require(channel == "G" || channel == "E" || channel == "T" || channel == "observation",
                "interface observable channel");
        const auto axes = interface.at("axes");
        require(axes.size() > 0, "interface observable stalk must be nonempty");
        std::map<std::string, std::size_t> axis_ids;
        for (std::size_t j = 0; j < axes.size(); ++j)
            require(axis_ids.emplace(label(axes.at(j)), j).second, "duplicate interface axis");
        dimensions.push_back(axes.size());
        std::vector<rfc::core::BoundaryOccurrenceInput> occurrences;
        const auto participants = interface.at("participants");
        require(participants.size() > 0, "interface must retain participants");
        for (std::size_t j = 0; j < participants.size(); ++j) {
            const auto participant = participants.at(j);
            members(participant, {"key", "view", "exists", "orientation", "magnitude", "map"});
            const auto id = label(participant.at("view"));
            require(view_ids.contains(id), "interface view not declared");
            const auto v = view_ids.at(id);
            const auto &field = views[v].fields.front();
            const auto map = participant.at("map");
            members(map, {"domain_quantity", "domain_units", "domain_channel", "codomain_quantity",
                          "codomain_units", "codomain_channel", "method", "evidence", "entries"});
            for (const auto *name : {"quantity", "units", "channel"}) {
                require(map.at(std::string("domain_") + name).text() ==
                            field.meaning.at(name).text(),
                        "restriction domain meaning differs from the admitted field");
                require(map.at(std::string("codomain_") + name).text() == interface.at(name).text(),
                        "restriction codomain meaning differs from the declared observable");
            }
            (void)label(map.at("method"));
            (void)label(map.at("evidence"));
            std::map<std::string, std::size_t> coordinates;
            for (std::size_t c = 0; c < field.ids.size(); ++c)
                coordinates.emplace(field.ids.at(c).text(), c);
            M action(axes.size(), field.values.size());
            std::set<M::Key> seen;
            for (std::size_t e = 0; e < map.at("entries").size(); ++e) {
                const auto entry = map.at("entries").at(e);
                members(entry, {"observable_id", "source_id", "coefficient"});
                const auto a = label(entry.at("observable_id")), c = label(entry.at("source_id"));
                require(axis_ids.contains(a) && coordinates.contains(c),
                        "restriction coordinate identity");
                const M::Key k{axis_ids.at(a), coordinates.at(c)};
                require(seen.insert(k).second, "duplicate restriction matrix entry");
                action.set(k.first, k.second, exact_rational(label(entry.at("coefficient"))));
            }
            const auto orientation = label(participant.at("orientation"));
            require(orientation == "positive" || orientation == "negative",
                    "interface orientation");
            const auto magnitude = exact_rational(label(participant.at("magnitude")));
            require(magnitude >= 0, "interface magnitude must be nonnegative");
            occurrences.emplace_back(v, participant.at("exists").boolean(),
                                     orientation == "negative", magnitude,
                                     std::vector<rfc::core::AspectValue>{},
                                     rfc::core::OccurrenceKey{key(participant.at("key"))});
            flags.push_back(Flag{v, i, std::move(action)});
        }
        builder.add_grade_one(rfc::core::RelationKey{i}, type, occurrences);
    }
    const auto source = builder.freeze();
    rfc::local::LinearStalkBundle stalks(source);
    for (std::size_t i = 0; i < views.size(); ++i)
        stalks.set_dimension({rfc::core::Grade{0}, static_cast<std::uint32_t>(i)},
                             views[i].fields.front().values.size());
    for (std::size_t i = 0; i < dimensions.size(); ++i)
        stalks.set_dimension({rfc::core::Grade{1}, static_cast<std::uint32_t>(i)}, dimensions[i]);
    rfc::local::LinearCellularSheaf sheaf(stalks);
    const auto offsets = stalks.offsets(rfc::core::Grade{0});
    const auto row_offsets = stalks.offsets(rfc::core::Grade{1});
    std::vector<std::set<std::size_t>> dependencies(row_offsets.back());
    for (std::size_t slot = 0; slot < flags.size(); ++slot) {
        const auto &flag = flags[slot];
        sheaf.set_restriction(rfc::core::Grade{1}, slot, flag.map);
        if (source.incidence_coefficient(rfc::core::Grade{1}, slot) != 0)
            for (const auto &[ij, value] : flag.map.entries()) {
                (void)value;
                dependencies[row_offsets[flag.relation] + ij.first].insert(offsets[flag.view] +
                                                                           ij.second);
            }
    }
    const auto delta = sheaf.coboundary(rfc::core::Grade{0});
    require(sheaf.is_cochain_complex(), "assay interface sheaf chain law");
    std::map<std::string, std::size_t> patterns;
    std::vector<Group> groups;
    for (std::size_t s = 0; s < cases.size(); ++s) {
        auto &c = cases[s];
        c.assessed.assign(delta.rows(), true);
        std::string pattern;
        for (auto supported : c.support)
            pattern += supported ? '1' : '0';
        for (std::size_t row = 0; row < delta.rows(); ++row) {
            for (auto coordinate : dependencies[row])
                if (!c.support[coordinate])
                    c.assessed[row] = false;
            pattern += c.assessed[row] ? '1' : '0';
        }
        if (!patterns.contains(pattern)) {
            Group group;
            for (std::size_t coordinate = 0; coordinate < c.support.size(); ++coordinate)
                if (c.support[coordinate])
                    group.coordinates.push_back(coordinate);
            for (std::size_t row = 0; row < c.assessed.size(); ++row)
                if (c.assessed[row])
                    group.rows.push_back(row);
            group.inclusion = M(delta.cols(), group.coordinates.size());
            Vector weights;
            for (std::size_t j = 0; j < group.coordinates.size(); ++j) {
                group.inclusion.set(group.coordinates[j], j, 1);
                weights.push_back(diagonal[group.coordinates[j]]);
            }
            group.restriction = M(group.rows.size(), delta.rows());
            for (std::size_t j = 0; j < group.rows.size(); ++j)
                group.restriction.set(j, group.rows[j], 1);
            group.action = group.restriction.multiply(delta).multiply(group.inclusion);
            auto metric = Metric::diagonal(std::move(weights));
            // The metric-orthogonal complement of ker(D) is im(M^-1 D^T).
            // Prepare the smaller native image frame, never an n by n projector.
            const auto image = group.action.transpose().image_frame();
            M frame(image.rows(), image.cols());
            std::vector<Vector> columns(image.cols(), Vector(image.rows()));
            for (const auto &[ij, value] : image.entries())
                columns[ij.second][ij.first] = value;
            for (std::size_t j = 0; j < columns.size(); ++j) {
                const auto column = metric.solve(columns[j]);
                for (std::size_t i = 0; i < column.size(); ++i)
                    if (column[i] != 0)
                        frame.set(i, j, column[i]);
            }
            group.projection =
                std::make_unique<PreparedProjection>(std::move(frame), std::move(metric));
            patterns.emplace(pattern, groups.size());
            groups.push_back(std::move(group));
        }
        c.group = patterns.at(pattern);
        groups[c.group].cases.push_back(s);
    }
    rfc::runtime::CpuTileExecutor executor({workers, workers, cases.size()});
    const auto execution = executor.run_tiles(
        cases.size(), [&](std::size_t s, std::size_t, const rfc::runtime::CpuTileStop &stop) {
            auto &c = cases[s];
            const auto &g = groups[c.group];
            const auto input = g.inclusion.transpose().apply(c.input);
            const auto residual = g.projection->apply(input, [&] { stop.check(); });
            const auto compatible = rfc::numeric::sub(input, residual);
            const auto defect = g.action.apply(compatible);
            require(std::all_of(defect.begin(), defect.end(), [](const auto &q) { return q == 0; }),
                    "compatible section has a nonzero assessed defect");
            c.orthogonality = g.projection->metric.inner(compatible, residual);
            require(c.orthogonality == 0 && rfc::numeric::add(compatible, residual) == input,
                    "assay section projection certificate failed");
            c.quadrance = g.projection->metric.inner(residual, residual);
            c.residual = g.inclusion.apply(residual);
            c.compatible = g.inclusion.apply(compatible);
            c.defect = g.restriction.transpose().apply(g.action.apply(input));
            stop.check();
        });
    auto out = Json::object(), result_sections = Json::array(), prepared = Json::array();
    for (std::size_t g = 0; g < groups.size(); ++g) {
        const auto &group = groups[g];
        auto item = Json::object(), coordinates = Json::array(), rows = Json::array();
        for (auto i : group.coordinates)
            coordinates.push(cardinality(i));
        for (auto i : group.rows)
            rows.push(cardinality(i));
        item.set("coordinates", coordinates);
        item.set("assessed_rows", rows);
        item.set("action", sparse_matrix(group.action));
        item.set("residual_frame", sparse_matrix(group.projection->frame));
        item.set("gram",
                 sparse_matrix(rfc::field::metric_gram(
                     group.projection->frame, group.projection->metric, group.projection->frame)));
        item.set("constraint_rank", cardinality(group.projection->frame.cols()));
        item.set("compatible_dimension",
                 cardinality(group.coordinates.size() - group.projection->frame.cols()));
        prepared.push(item);
    }
    for (std::size_t s = 0; s < cases.size(); ++s) {
        const auto &c = cases[s];
        auto item = Json::object(), fields = Json::array(), constraints = Json::array();
        item.set("id", sections.at(s));
        item.set("prepared_group", cardinality(c.group));
        for (std::size_t v = 0; v < views.size(); ++v) {
            const auto &view = views[v];
            const auto &field = view.fields[s];
            auto values = Json::object(), input = Json::array(), compatible = Json::array(),
                 residual = Json::array();
            for (std::size_t j = offsets[v]; j < offsets[v + 1]; ++j) {
                input.push(c.support[j] ? Json::string(c.input[j].get_str()) : Json{});
                compatible.push(c.support[j] ? Json::string(c.compatible[j].get_str()) : Json{});
                residual.push(c.support[j] ? Json::string(c.residual[j].get_str()) : Json{});
            }
            values.set("view", Json::string(view.id));
            values.set("ids", field.ids);
            values.set("keys", field.keys);
            values.set("source_grade", Json::integer(field.grade));
            values.set("source_digest", Json::string(field.owner->digest()));
            auto meaning = Json::object();
            // Original admission documents remain bound by source_digest and are
            // retained in derivation lineage. Report path-independent field meaning.
            for (const auto *name :
                 {"id", "key", "context", "channel", "quantity", "units", "method", "value_type"})
                meaning.set(name, field.meaning.at(name));
            values.set("meaning", meaning);
            values.set("input", input);
            values.set("compatible", compatible);
            values.set("residual", residual);
            fields.push(values);
        }
        for (std::size_t r = 0; r < delta.rows(); ++r) {
            auto constraint = Json::object(), missing = Json::array();
            for (auto i : dependencies[r])
                if (!c.support[i])
                    missing.push(cardinality(i));
            constraint.set("row", cardinality(r));
            constraint.set("assessed", Json::boolean(c.assessed[r]));
            constraint.set("input_defect",
                           c.assessed[r] ? Json::string(c.defect[r].get_str()) : Json{});
            constraint.set("compatible_defect", c.assessed[r] ? Json::string("0") : Json{});
            constraint.set("missing_coordinates", missing);
            constraints.push(constraint);
        }
        item.set("fields", fields);
        item.set("constraints", constraints);
        item.set("residual_quadrance", Json::string(c.quadrance.get_str()));
        item.set("orthogonality", Json::string(c.orthogonality.get_str()));
        item.set("reconstruction_verified", Json::boolean(true));
        item.set("complete",
                 Json::boolean(
                     std::all_of(c.assessed.begin(), c.assessed.end(), [](auto x) { return x; }) &&
                     std::all_of(c.support.begin(), c.support.end(), [](auto x) { return x; })));
        result_sections.push(item);
    }
    rfc::storage::TemporalEndpoint endpoint(source);
    const auto encoded = rfc::storage::encode_native_series(
        rfc::storage::NativeTemporalSeries(endpoint, 1, source_codec_limits(source_byte_limit())));
    out.set(
        "interface_source_sha256",
        Json::string("sha256:" + sha256(std::string(reinterpret_cast<const char *>(encoded.data()),
                                                    encoded.size()))));
    out.set("interface_grade", Json::integer(1));
    out.set("interface_chain_verified", Json::boolean(sheaf.is_cochain_complex()));
    out.set("views", request.at("views"));
    out.set("interfaces", interfaces);
    out.set("metric", penalty);
    out.set("coboundary", sparse_matrix(delta));
    auto coordinate_map = Json::array(), row_map = Json::array();
    for (std::size_t v = 0; v < views.size(); ++v)
        for (std::size_t c = 0; c < views[v].fields.front().ids.size(); ++c) {
            auto coordinate = Json::object();
            coordinate.set("view", Json::string(views[v].id));
            coordinate.set("source_id", views[v].fields.front().ids.at(c));
            coordinate.set("source_key", views[v].fields.front().keys.at(c));
            coordinate.set("source_grade", Json::integer(views[v].fields.front().grade));
            coordinate_map.push(coordinate);
        }
    for (std::size_t i = 0; i < interfaces.size(); ++i)
        for (std::size_t a = 0; a < interfaces.at(i).at("axes").size(); ++a) {
            auto coordinate = Json::object();
            coordinate.set("interface", interfaces.at(i).at("id"));
            coordinate.set("observable", interfaces.at(i).at("axes").at(a));
            row_map.push(coordinate);
        }
    out.set("coordinates", coordinate_map);
    out.set("constraint_coordinates", row_map);
    out.set("prepared_groups", prepared);
    out.set("sections", result_sections);
    auto receipt = Json::object();
    receipt.set("owner", Json::string("rfc::local::LinearCellularSheaf"));
    receipt.set("projection_owner", Json::string("rfc::numeric::PreparedExactSolve"));
    receipt.set("worker_slots", cardinality(execution.worker_slots));
    receipt.set("completed_tiles", cardinality(execution.completed_tiles));
    receipt.set("prepared_solves", cardinality(groups.size()));
    receipt.set("source_coordinates", cardinality(diagonal.size()));
    receipt.set("constraint_coordinates", cardinality(delta.rows()));
    out.set("execution", receipt);
    return out;
}
} // namespace gadget::assay_detail
