#include "assay_field_detail.hpp"
#include "genomic_detail.hpp"
#include "scientific_projection.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/local/cellular_sheaf.hpp>
#include <rfc/runtime/cpu_tiles.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <rfc/storage/temporal_diagnostics.hpp>

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
                    std::shared_ptr<CpuResources> resources,
                    std::optional<AssayFieldsSource> *retained) {
    members(request,
            {"schema", "operation", "views", "sections", "interfaces", "metric", "support_policy"},
            {"workers", "complete_families"});
    require(!request.has("complete_families") ||
                json_object_is_type(request.at("complete_families").raw(), json_type_boolean),
            "complete assay section families requires a boolean output selection");
    const bool complete_families = request.has("complete_families") &&
                                  request.at("complete_families").boolean();
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
    {
        using N = NativeRecordValue;
        using namespace rfc::storage;
        const auto integer = [](std::size_t x) {
            return N(rfc::numeric::Integer(std::to_string(x)));
        };
        TemporalEndpoint endpoint(source);
        std::vector<std::uint64_t> view_carriers, interface_carriers, stalk_carriers,
                                   active_carriers, assessed_carriers;
        const auto carrier = [&](TemporalCarrierKind kind, std::size_t size,
                                 std::optional<StableCellKey> cell = {}) {
            const auto key = endpoint.carriers.size() + 1;
            endpoint.carriers.push_back({key, kind, cell ? cell->grade : rfc::core::Grade{}, cell,
                                         key, rfc::numeric::IndexDomain::linear(size)});
            return key;
        };
        for (std::size_t v = 0; v < views.size(); ++v)
            view_carriers.push_back(carrier(TemporalCarrierKind::CellStalk,
                                            views[v].fields.front().values.size(),
                                            StableCellKey{rfc::core::Grade{0}, v}));
        for (std::size_t i = 0; i < dimensions.size(); ++i)
            interface_carriers.push_back(carrier(TemporalCarrierKind::CellStalk, dimensions[i],
                                                 StableCellKey{rfc::core::Grade{1}, i}));
        stalk_carriers = view_carriers;
        stalk_carriers.insert(stalk_carriers.end(), interface_carriers.begin(),
                              interface_carriers.end());
        const auto domain = carrier(TemporalCarrierKind::Explicit, delta.cols());
        const auto codomain = carrier(TemporalCarrierKind::Explicit, delta.rows());
        SourceGradeFields store(endpoint);
        N::Array native_stalks, restrictions, prepared, result_sections, coordinates, row_map;
        for (auto key : stalk_carriers)
            native_stalks.elements.push_back(integer(key));
        for (std::size_t slot = 0; slot < flags.size(); ++slot) {
            const auto &flag = flags[slot];
            TemporalTensorOwner owner;
            owner.kind = TemporalOwnerKind::Occurrence;
            owner.occurrence = {rfc::core::Grade{1},
                                source.grade(rfc::core::Grade{1}).occurrence_identities.at(slot)};
            restrictions.elements.push_back(store.map("restriction:" + std::to_string(slot),
                                                      view_carriers[flag.view],
                                                      interface_carriers[flag.relation], flag.map,
                                                      TemporalTensorRole::SheafRestriction, owner));
        }
        const auto delta_key = store.map("coboundary", domain, codomain, delta);
        const auto metric_key =
            store.on_carrier("source.metric_diagonal", domain,
                             SourceGradeFields::tensor(
                                 diagonal.size(), [&](auto i) { return diagonal[i]; },
                                 rfc::numeric::StructuralSupport::full(diagonal.size())));
        for (std::size_t g = 0; g < groups.size(); ++g) {
            const auto &group = groups[g];
            const auto active = carrier(TemporalCarrierKind::Explicit, group.coordinates.size());
            const auto assessed = carrier(TemporalCarrierKind::Explicit, group.rows.size());
            active_carriers.push_back(active);
            assessed_carriers.push_back(assessed);
            const auto frame =
                carrier(TemporalCarrierKind::Explicit, group.projection->frame.cols());
            N::Array coordinate_ids, row_ids;
            for (auto i : group.coordinates)
                coordinate_ids.elements.push_back(integer(i));
            for (auto i : group.rows)
                row_ids.elements.push_back(integer(i));
            prepared.elements.emplace_back(
                N::Object{{{"coordinates", N(std::move(coordinate_ids))},
                           {"assessed_rows", N(std::move(row_ids))},
                           {"inclusion", store.map("group:" + std::to_string(g) + ".inclusion",
                                                   active, domain, group.inclusion)},
                           {"restriction", store.map("group:" + std::to_string(g) + ".restriction",
                                                     codomain, assessed, group.restriction)},
                           {"action", store.map("group:" + std::to_string(g) + ".action", active,
                                                assessed, group.action)},
                           {"residual_frame", store.map("group:" + std::to_string(g) + ".frame",
                                                        frame, active, group.projection->frame)},
                           {"gram", store.map("group:" + std::to_string(g) + ".gram", frame, frame,
                                              rfc::field::metric_gram(group.projection->frame,
                                                                      group.projection->metric,
                                                                      group.projection->frame))},
                           {"constraint_rank", integer(group.projection->frame.cols())},
                           {"compatible_dimension",
                            integer(group.coordinates.size() - group.projection->frame.cols())}}});
        }
        const auto component = [&](std::string name, std::uint64_t target, const Vector &values,
                                   const std::vector<unsigned char> &support, std::size_t begin,
                                   std::size_t end) {
            std::vector<std::size_t> active;
            for (std::size_t j = begin; j < end; ++j)
                if (support[j])
                    active.push_back(j - begin);
            return store.on_carrier(
                std::move(name), target,
                SourceGradeFields::tensor(
                    end - begin, [&](auto i) { return values[begin + i]; },
                    rfc::numeric::StructuralSupport::packed(end - begin, active)));
        };
        for (std::size_t s = 0; s < cases.size(); ++s) {
            const auto &c = cases[s];
            const auto name = label(sections.at(s));
            N::Array fields, constraints;
            for (std::size_t v = 0; v < views.size(); ++v) {
                const auto &field = views[v].fields[s];
                auto meaning = Json::object();
                for (const auto *key : {"id", "key", "context", "channel", "quantity", "units",
                                        "method", "value_type"})
                    meaning.set(key, field.meaning.at(key));
                N::Object record{{{"view", N(views[v].id)},
                                  {"ids", import_record_view(field.ids)},
                                  {"keys", import_record_view(field.keys)},
                                  {"source_grade", integer(field.grade)},
                                  {"source_digest", N(field.owner->digest())},
                                  {"meaning", import_record_view(meaning)}}};
                for (const auto &[part, values] :
                     std::vector<std::pair<std::string, const Vector *>>{
                         {"input", &c.input},
                         {"compatible", &c.compatible},
                         {"residual", &c.residual}})
                    record.members.emplace_back(
                        part, component(name + "." + views[v].id + "." + part, view_carriers[v],
                                        *values, c.support, offsets[v], offsets[v + 1]));
                fields.elements.emplace_back(std::move(record));
            }
            for (std::size_t i = 0; i < dimensions.size(); ++i) {
                const auto defect_key = component(
                    name + ".interface:" + std::to_string(i) + ".defect", interface_carriers[i],
                    c.defect, c.assessed, row_offsets[i], row_offsets[i + 1]);
                constraints.elements.emplace_back(
                    N::Object{{{"field", defect_key},
                               {"begin", integer(row_offsets[i])},
                               {"end", integer(row_offsets[i + 1])}}});
            }
            result_sections.elements.emplace_back(N::Object{
                {{"id", N(name)},
                 {"prepared_group", integer(c.group)},
                 {"fields", N(std::move(fields))},
                 {"constraints", N(std::move(constraints))},
                 {"residual_quadrance", N(c.quadrance)},
                 {"orthogonality", N(c.orthogonality)},
                 {"reconstruction_verified", N(true)},
                 {"complete",
                  N(std::all_of(c.assessed.begin(), c.assessed.end(), [](auto x) { return x; }) &&
                    std::all_of(c.support.begin(), c.support.end(), [](auto x) { return x; }))}}});
        }
        for (std::size_t v = 0; v < views.size(); ++v)
            for (std::size_t c = 0; c < views[v].fields.front().ids.size(); ++c)
                coordinates.elements.emplace_back(
                    N::Object{{{"view", N(views[v].id)},
                               {"source_id", N(views[v].fields.front().ids.at(c).text())},
                               {"source_key", N(views[v].fields.front().keys.at(c).text())},
                               {"source_grade", integer(views[v].fields.front().grade)}}});
        // Add complete families after the existing fields, preserving their keys.
        // Directions are correlated across the original assessed coordinates.
        for (std::size_t g = 0; complete_families && g < groups.size(); ++g) {
            const auto &group = groups[g];
            const auto kernel = group.action.kernel_frame();
            const auto row_dependencies = group.action.transpose().kernel_frame();
            const auto rank = group.projection->frame.cols();
            require(kernel.cols() == group.coordinates.size() - rank &&
                        row_dependencies.cols() == group.rows.size() - rank &&
                        group.action.multiply(kernel).entries().empty() &&
                        group.action.transpose().multiply(row_dependencies).entries().empty(),
                    "assay section complete family certificate failed");
            const auto parameters = carrier(TemporalCarrierKind::Explicit, kernel.cols());
            const auto dependencies = carrier(TemporalCarrierKind::Explicit,
                                              row_dependencies.cols());
            auto &members = std::get<N::Object>(prepared.elements.at(g).data).members;
            members.emplace_back("kernel_frame", store.map("group:" + std::to_string(g) +
                ".kernel", parameters, active_carriers.at(g), kernel));
            members.emplace_back("row_dependencies", store.map("group:" + std::to_string(g) +
                ".row_dependencies", dependencies, assessed_carriers.at(g), row_dependencies));
        }
        for (std::size_t i = 0; i < interfaces.size(); ++i)
            for (std::size_t a = 0; a < interfaces.at(i).at("axes").size(); ++a)
                row_map.elements.emplace_back(
                    N::Object{{{"interface", N(interfaces.at(i).at("id").text())},
                               {"observable", N(interfaces.at(i).at("axes").at(a).text())}}});
        N::Array dependency_rows;
        for (const auto &row : dependencies) {
            N::Array values;
            for (auto j : row)
                values.elements.push_back(integer(j));
            dependency_rows.elements.emplace_back(std::move(values));
        }
        const auto reopened = temporal_linear_cellular_sheaf(endpoint, stalk_carriers);
        require(reopened.sheaf().coboundary(rfc::core::Grade{0}) == delta,
                "retained sheaf coboundary differs from its generating maps");
        const auto encoded = encode_native_series(NativeTemporalSeries(
            TemporalEndpoint(source), 1, source_codec_limits(source_byte_limit())));
        N declaration(N::Object{
            {{"schema", N(std::string("gadget.sections.fields.v1"))},
             {"stalk_carriers", N(std::move(native_stalks))},
             {"restrictions", N(std::move(restrictions))},
             {"coboundary", delta_key},
             {"metric_diagonal", metric_key},
             {"request", import_record_view(request)},
             {"prepared_groups", N(std::move(prepared))},
             {"sections", N(std::move(result_sections))},
             {"coordinates", N(std::move(coordinates))},
             {"constraint_coordinates", N(std::move(row_map))},
             {"dependencies", N(std::move(dependency_rows))},
             {"interface_source_sha256",
              N("sha256:" + sha256(std::string_view(reinterpret_cast<const char *>(encoded.data()),
                                                    encoded.size())))},
             {"execution",
              N(N::Object{{{"owner", N(std::string("rfc::local::LinearCellularSheaf"))},
                           {"projection_owner", N(std::string("rfc::numeric::PreparedExactSolve"))},
                           {"worker_slots", integer(execution.worker_slots)},
                           {"completed_tiles", integer(execution.completed_tiles)},
                           {"prepared_solves", integer(groups.size())},
                           {"source_coordinates", integer(diagonal.size())},
                           {"constraint_coordinates", integer(delta.rows())}}})}}});
        if (retained) {
            retained->emplace(AssayFieldsSource{std::move(endpoint), std::move(declaration), {}});
            return Json::object();
        }
        return section_fields_view(NativeRecord::create_from_endpoint(
            std::move(declaration), new_native_record_namespace(), {}, std::move(endpoint)));
    }
}

Json section_fields_view(const NativeRecord &record) {
    if (record.at("schema").text() == "gadget.sections.selected.v1") {
        const auto complete = NativeRecord::decode(record.at("sections_publication").bytes());
        return section_fields_view(complete.at("result").at("result"));
    }
    const auto request = record.at("request");
    const auto matrix = [&](const NativeRecord &key) {
        return sparse_matrix(rfc::storage::temporal_tensor_matrix(retained_tensor(record, key)));
    };
    const auto vector = [&](const NativeRecord &key) {
        auto result = Json::array();
        const auto &tensor = retained_tensor(record, key).tensor;
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i) {
            const auto value = tensor.state(i);
            result.push(value ? Json::string(value->get_str()) : Json{});
        }
        return result;
    };
    auto out = Json::object(), groups = Json::array(), sections = Json::array();
    for (std::size_t g = 0; g < record.at("prepared_groups").size(); ++g) {
        const auto group = record.at("prepared_groups").at(g);
        auto item = Json::object();
        for (const auto *name :
             {"coordinates", "assessed_rows", "constraint_rank", "compatible_dimension"})
            item.set(name, render_record_view(group.at(name)));
        for (const auto *name : {"action", "residual_frame", "gram"})
            item.set(name, matrix(group.at(name)));
        for (const auto *name : {"kernel_frame", "row_dependencies"})
            if (group.has(name))
                item.set(name, matrix(group.at(name)));
        groups.push(item);
    }
    for (std::size_t s = 0; s < record.at("sections").size(); ++s) {
        const auto section = record.at("sections").at(s);
        auto item = Json::object(), fields = Json::array(), constraints = Json::array();
        for (const auto *name : {"id", "prepared_group", "residual_quadrance", "orthogonality",
                                 "reconstruction_verified", "complete"})
            item.set(name, render_record_view(section.at(name)));
        std::vector<unsigned char> support;
        for (std::size_t v = 0; v < section.at("fields").size(); ++v) {
            const auto field = section.at("fields").at(v);
            auto values = Json::object();
            for (const auto *name :
                 {"view", "ids", "keys", "source_grade", "source_digest", "meaning"})
                values.set(name, render_record_view(field.at(name)));
            for (const auto *part : {"input", "compatible", "residual"})
                values.set(part, vector(field.at(part)));
            const auto &tensor = retained_tensor(record, field.at("input")).tensor;
            for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i)
                support.push_back(tensor.state(i).has_value());
            fields.push(values);
        }
        std::size_t row = 0;
        for (std::size_t i = 0; i < section.at("constraints").size(); ++i) {
            const auto &tensor =
                retained_tensor(record, section.at("constraints").at(i).at("field")).tensor;
            for (std::size_t j = 0; j < tensor.domain().logical_size(); ++j, ++row) {
                const auto state = tensor.state(j);
                auto constraint = Json::object(), missing = Json::array();
                const auto dependencies = record.at("dependencies").at(row);
                for (std::size_t k = 0; k < dependencies.size(); ++k) {
                    const auto coordinate = dependencies.at(k).integer().get_ui();
                    if (!support.at(coordinate))
                        missing.push(cardinality(coordinate));
                }
                constraint.set("row", cardinality(row));
                constraint.set("assessed", Json::boolean(state.has_value()));
                constraint.set("input_defect", state ? Json::string(state->get_str()) : Json{});
                constraint.set("compatible_defect", state ? Json::string("0") : Json{});
                constraint.set("missing_coordinates", missing);
                constraints.push(constraint);
            }
        }
        item.set("fields", fields);
        item.set("constraints", constraints);
        sections.push(item);
    }
    out.set("interface_source_sha256", Json::string(record.at("interface_source_sha256").text()));
    out.set("interface_grade", Json::integer(1));
    out.set("interface_chain_verified", Json::boolean(true));
    for (const auto *name : {"views", "interfaces", "metric"})
        out.set(name, render_record_view(request.at(name)));
    out.set("coboundary", matrix(record.at("coboundary")));
    for (const auto *name : {"coordinates", "constraint_coordinates", "execution"})
        out.set(name, render_record_view(record.at(name)));
    out.set("prepared_groups", groups);
    out.set("sections", sections);
    return out;
}
} // namespace gadget::assay_detail
