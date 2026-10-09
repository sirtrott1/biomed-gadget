#include "research_ablation.hpp"
#include "research_ensemble.hpp"
#include "research_factored.hpp"
#include "research_prepared.hpp"
#include "research_response.hpp"
#include "research_structural_calibration.hpp"
#include "research_temporal.hpp"
#include "scientific_detail.hpp"
#include "scientific_native_hodge.hpp"
#include "scientific_native_novelty.hpp"
#include "scientific_native_observation.hpp"
#include "scientific_native_response.hpp"
#include <gadget/research.hpp>
#include <gadget/resources.hpp>
#include <gadget/sections.hpp>
#include <gadget/structure.hpp>
#include <gadget/study.hpp>
#include <mutex>
#include <rfc/algebra/materialize.hpp>
#include <rfc/database/native_field_cache.hpp>
#include <rfc/database/native_tensor_window.hpp>
#include <rfc/database/parallel_native_field.hpp>
#include <rfc/database/resource_native_field.hpp>
#include <rfc/local/interior_interface_solve.hpp>
#include <rfc/temporal/moment_change.hpp>
#include <rfc/temporal/operation_word.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc::storage;
const std::size_t coordinate_limit = std::numeric_limits<std::size_t>::max();
std::size_t field_memory_headroom() {
    const auto room=rfc::runtime::detect_hardware().advisory_memory_headroom();
    require(room && *room,"native field services require observed memory headroom");
    return rfc::machine::checked_size(*room);
}
rfc::database::NativeFieldCacheLimits field_cache_limits() {
    const auto maximum=std::numeric_limits<std::size_t>::max();
    const auto bytes=field_memory_headroom();
    rfc::database::NativeFieldCacheLimits limits;
    limits.plans=limits.results=limits.result_slots=limits.in_flight=maximum;
    limits.idle_workspaces=std::max(std::size_t{1},rfc::runtime::detect_hardware().allowed_cpus.size());
    limits.idle_workspace_slots=maximum;
    limits.descriptor_bytes=limits.pending_descriptor_bytes=bytes;
    limits.retained_result_bytes=limits.idle_workspace_bytes=bytes;
    return limits;
}
std::uint64_t publication_cutoff(const Json &value) { return value.text() == "0" ? 0 : key(value); }
void check_values(const Vector &values) {
    require(values.size() <= coordinate_limit, "result coordinate budget");
}
Json rational(const Rational &value) {

    auto result = Json::object();
    result.set("rational", Json::string(value.get_str()));
    return result;
}
Json vector_json(const Vector &values) {
    check_values(values);
    auto result = Json::array();
    for (const auto &value : values)
        result.push(rational(value));
    return result;
}
struct PreparedContext {
    std::shared_ptr<const rfc::core::RelationalComplex> source;
    std::vector<std::vector<std::uint32_t>> to_ambient;
    std::shared_ptr<const NativeStudy> study;
    std::unique_ptr<PreparedTemporalStudy> temporal;
};
} // namespace
struct ScientificSession::State {
    std::shared_ptr<const ScientificStudy> owner;
    const ScientificStudy::State &source;
    std::map<std::string, PreparedContext> contexts;
    std::map<std::string, std::unique_ptr<State>> publication_sessions;
    std::map<std::string, rfc::local::PreparedInteriorInterfaceSolve> green_solves;
    std::map<std::uint64_t, std::string> green_spaces;
    std::map<std::string, rfc::database::NativeTensorWindow> windows;
    rfc::database::NativeFieldCache field_cache;
    scientific_detail::DeclaredResponsePlans declared_responses;
    scientific_detail::NativeHodgePlans native_hodge;
    scientific_detail::NativeObservationPlans native_observations;
    scientific_detail::NativeNoveltyPlans native_novelty;
    std::shared_ptr<CpuResources> resources;
    ComputeOptions compute;
    std::map<std::string, std::unique_ptr<PackedOperation>> packed;
    struct Scheduled {
        rfc::database::PreparedParallelNativeField operation;
        rfc::database::ParallelNativeFieldWorkspace workspace;
        std::optional<rfc::database::CpuFieldCalibration> calibration;
        explicit Scheduled(rfc::database::PreparedParallelNativeField prepared)
            : operation(std::move(prepared)),
              workspace(rfc::database::make_parallel_workspace(operation)) {}
    };
    std::map<std::pair<std::string, std::size_t>, std::unique_ptr<Scheduled>> scheduled;
    mutable std::mutex mutex;
    explicit State(std::shared_ptr<const ScientificStudy> study, std::shared_ptr<CpuResources> cpu,
                   ComputeOptions options)
        : owner(std::move(study)), source(*owner->state_), field_cache(field_cache_limits()), resources(std::move(cpu)),
          compute(options) {}
    CpuResources &cpu() {
        if (!resources)
            resources = CpuResources::local();
        return *resources;
    }
    PreparedTemporalStudy &temporal(PreparedContext &context) {
        if (!context.temporal)
            context.temporal = std::make_unique<PreparedTemporalStudy>(context.study);
        return *context.temporal;
    }
    PreparedContext &context(const std::string &name) {
        require(source.contexts.contains(name), "unknown context: " + name);
        if (contexts.contains(name))
            return contexts.at(name);
        const auto endpoint = source.directory->pin_snapshot();
        const auto &declared = source.contexts.at(name);
        rfc::algebra::SubcomplexMask mask(endpoint->source);
        if (declared.block) {
            const auto block =
                std::find_if(endpoint->local_blocks.begin(), endpoint->local_blocks.end(),
                             [&](const auto &b) { return b.block_key == *declared.block; });
            const StableSourceIndex index(endpoint->source);
            for (const auto &cell : block->cells)
                mask.select(index.cell(cell).value());
        } else {
            for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k) {
                const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(k)};
                for (std::size_t i = 0; i < endpoint->source.grade(grade).cell_count(); ++i)
                    mask.select({grade, static_cast<std::uint32_t>(i)});
            }
        }
        // Explicit common ambient grade-zero carrier, including isolated cells.
        for (std::size_t i = 0; i < endpoint->source.grade(rfc::core::Grade{0}).cell_count(); ++i)
            mask.select({rfc::core::Grade{0}, static_cast<std::uint32_t>(i)});
        for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k)
            require(mask.grade(rfc::core::Grade{static_cast<std::uint16_t>(k)}).count() <=
                        coordinate_limit,
                    "context coordinate budget");
        auto materialized = rfc::algebra::materialize_subcomplex(endpoint->source, mask);
        std::vector<rfc::field::ExactMetric> metrics;
        const auto declared_metrics = declared.declaration.at("metrics");
        for (std::size_t k = 0; k < materialized.to_ambient.size(); ++k) {
            const auto &coordinates = materialized.to_ambient[k];
            require(coordinates.size() <= coordinate_limit, "context coordinate budget");
            const auto metric = declared_metrics.at(k);
            if (metric.at("kind").text() == "identity")
                metrics.push_back(rfc::field::ExactMetric::identity(coordinates.size()));
            else {
                Vector diagonal;
                for (const auto coordinate : coordinates)
                    diagonal.push_back(source.metric_value(metric,coordinate));
                metrics.push_back(rfc::field::ExactMetric::diagonal(std::move(diagonal)));
            }
        }
        PreparedContext prepared;
        prepared.to_ambient = std::move(materialized.to_ambient);
        prepared.source =
            std::make_shared<const rfc::core::RelationalComplex>(std::move(materialized.complex));
        prepared.study = std::make_shared<const NativeStudy>(prepared.source, std::move(metrics));
        return contexts.emplace(name, std::move(prepared)).first->second;
    }
    const TemporalTensor &tensor(const std::string &field) const {
        require(source.fields.contains(field), "unknown field: " + field);
        return *temporal_endpoint_detail::tensors(*source.directory->pin_snapshot())
                    .at(source.fields.at(field).tensor);
    }
    std::pair<Vector, std::uint16_t> input(const std::string &name, const std::string &field,
                                           bool legacy_budget = true) {
        const auto &t = tensor(field);
        const auto &meaning = source.fields.at(field).declaration;
        require(meaning.at("value_type").text() != "categorical",
                "categorical codes cannot be linear field amplitudes");
        require(meaning.at("context").text() == name, "field/context mismatch");
        const auto carrier =
            temporal_endpoint_detail::carriers(*source.directory->pin_snapshot()).at(*t.output);
        require(carrier->kind == TemporalCarrierKind::SourceGrade &&
                    carrier->domain.extents().size() == 1,
                "scalar tower field requires a rank-one SourceGrade carrier");
        auto &prepared = context(name);
        const auto grade = carrier->grade.value;
        Vector values;
        for (const auto coordinate : prepared.to_ambient.at(grade)) {
            require(
                t.tensor.support().contains(coordinate),
                "partial observations require an explicit extension before mathematical execution");
            values.push_back(*t.tensor.state(coordinate));
        }
        if (legacy_budget)
            check_values(values);
        return {std::move(values), grade};
    }
    Json ensemble_query(const Json &query, bool validate) {
        const auto name = query.at("context").text();
        const auto ids = query.at("fields");
        source_array(ids, "ensemble fields");
        require(ids.size() > 0, "ensemble requires observed fields");
        std::vector<Vector> values;
        std::optional<std::uint16_t> grade;
        for (std::size_t s = 0; s < ids.size(); ++s) {
            const auto id = ids.at(s).text();
            comparable(ids.at(0).text(), id);
            auto field = input(name, id, false);
            require(!grade || *grade == field.second, "ensemble requires one native source grade");
            grade = field.second;
            values.push_back(std::move(field.first));
        }
        (void)source.evidence(query.at("interpretation").at("evidence"));
        auto &prepared = context(name);
        auto result = research_detail::relational_field_ensemble(
            *prepared.study, query, rfc::core::Grade{*grade}, values, validate, resources);
        if (!validate) {
            result.set("native_carrier",
                       Json::string(std::to_string(*tensor(ids.at(0).text()).output)));
            auto grades = Json::array();
            const auto endpoint = source.directory->pin_snapshot();
            for (std::size_t k = 0; k < prepared.to_ambient.size(); ++k) {
                auto coordinates_ = Json::array();
                for (const auto coordinate : prepared.to_ambient[k]) {
                    const auto cell = stable_cell_key(
                        endpoint->source,
                        {rfc::core::Grade{static_cast<std::uint16_t>(k)}, coordinate});
                    auto value = Json::object();
                    value.set("ambient_coordinate", cardinality(coordinate));
                    value.set("native_key", Json::string(std::to_string(cell.key)));
                    value.set("entity", source.cell_record(cell).at("id"));
                    coordinates_.push(value);
                }
                grades.push(coordinates_);
            }
            result.set("coordinate_carriers", grades);
        }
        return result;
    }
    Json structural_calibration_query(const Json &query, bool validate) {
        const auto name = query.at("context").text();
        const auto axes = query.at("fields");
        source_array(axes, "structural fields");
        require(axes.size(), "structural calibration requires field axes");
        const auto first = axes.at(0).text();
        auto origin = input(name, first, false);
        const rfc::core::Grade grade{origin.second};
        std::map<std::string, Vector> values;
        for (const auto &id : research_detail::structural_calibration_field_ids(query)) {
            auto field = input(name, id, false);
            require(field.second == grade.value && tensor(first).output == tensor(id).output,
                    "structural fields and probes require the identical actual source carrier");
            values.emplace(id, std::move(field.first));
        }
        for (std::size_t a = 0; a < axes.size(); ++a)
            comparable(first, axes.at(a).text());
        if (query.has("uncertainty")) {
            const auto uncertainty = query.at("uncertainty");
            (void)source.evidence(uncertainty.at("evidence"));
            const auto scenarios = uncertainty.at("scenarios");
            for (std::size_t s = 0; s < scenarios.size(); ++s) {
                const auto fields = scenarios.at(s).at("fields");
                require(fields.size() == axes.size(), "uncertainty field axes");
                for (std::size_t a = 0; a < axes.size(); ++a)
                    comparable(axes.at(a).text(), fields.at(a).text());
            }
        }
        const auto hypotheses = query.at("hypotheses");
        for (std::size_t h = 0; h < hypotheses.size(); ++h)
            if (hypotheses.at(h).has("probes")) {
                const auto probes = hypotheses.at(h).at("probes");
                for (std::size_t p = 0; p < probes.size(); ++p)
                    comparable(probes.at(0).text(), probes.at(p).text());
            }
        const auto law = query.at("law");
        (void)source.evidence(law.at("evidence"));
        const auto declarations = law.at("blocks");
        source_array(declarations, "source exchange blocks");
        auto &prepared = context(name);
        const auto &complex = *prepared.source;
        const StableSourceIndex index(complex);
        std::set<std::string> block_ids;
        std::set<std::size_t> used;
        std::vector<research_detail::ExchangeBlock> blocks;
        for (std::size_t b = 0; b < declarations.size(); ++b) {
            const auto declaration = declarations.at(b);
            members(declaration, {"id", "relation", "cells"});
            require(block_ids.insert(text(declaration.at("id"), source_byte_limit())).second,
                    "duplicate exchange block identity");
            const auto relation = declaration.at("relation").text();
            require(source.has_cell(relation), "unknown exchange relation");
            const auto retained = index.cell(source.resolve_cell(relation));
            require(retained && retained->grade.value == grade.value + 1u,
                    "exchange block must name a retained primary relation directly above its field "
                    "grade");
            const auto &arena = complex.grade(retained->grade);
            const auto [begin, end] = arena.incidence_range(retained->index);
            std::set<std::size_t> participants;
            auto occurrences = Json::array();
            for (std::size_t slot = begin; slot < end; ++slot) {
                if (arena.aspects.get(complex.signature().existence_aspect(), slot))
                    participants.insert(arena.lower_index[slot]);
                const StableOccurrenceKey occurrence{retained->grade,
                                                     arena.occurrence_identities[slot]};
                auto witness = Json::parse(source.occurrence_record(occurrence).dump_compact());
                witness.set("coordinate", cardinality(arena.lower_index[slot]));
                witness.set("exists", Json::boolean(arena.aspects.get(
                                          complex.signature().existence_aspect(), slot)));
                witness.set(
                    "coefficient",
                    Json::string(complex.incidence_coefficient(retained->grade, slot).get_str()));
                occurrences.push(witness);
            }
            const auto cells = declaration.at("cells");
            source_array(cells, "exchange participant cells");
            std::set<std::size_t> declared;
            research_detail::ExchangeBlock block;
            for (std::size_t c = 0; c < cells.size(); ++c) {
                const auto id = cells.at(c).text();
                require(source.has_cell(id), "unknown exchange participant");
                const auto cell = index.cell(source.resolve_cell(id));
                require(cell && cell->grade == grade && declared.insert(cell->index).second &&
                            used.insert(cell->index).second,
                        "exchange participants must be distinct retained coordinates in disjoint "
                        "blocks");
                block.coordinates.push_back(cell->index);
            }
            require(!declared.empty() && declared == participants,
                    "exchange block must retain every distinct existing occurrence participant, "
                    "including cancelled and supported-zero participants");
            block.witness = Json::parse(declaration.dump_compact());
            block.witness.set("coordinates", Json::array());
            for (const auto coordinate : block.coordinates)
                block.witness.at("coordinates").push(cardinality(coordinate));
            block.witness.set("source_relation", source.cell_record(source.resolve_cell(relation)));
            block.witness.set("source_occurrences", occurrences);
            blocks.push_back(std::move(block));
        }
        auto result = research_detail::structural_reference_calibration(
            *prepared.study, query, grade, values, blocks, validate, resources);
        if (!validate) {
            result.set("native_carrier", Json::string(std::to_string(*tensor(first).output)));
            auto coordinates = Json::array();
            for (std::size_t i = 0; i < complex.grade(grade).cell_count(); ++i) {
                const auto cell = stable_cell_key(complex, {grade, static_cast<std::uint32_t>(i)});
                auto record = Json::object();
                record.set("coordinate", cardinality(i));
                record.set("ambient_coordinate",
                           cardinality(prepared.to_ambient.at(grade.value).at(i)));
                record.set("source_cell", source.cell_record(cell));
                coordinates.push(record);
            }
            result.set("coordinate_carrier", coordinates);
        }
        return result;
    }
    void comparable(const std::string &old_name, const std::string &new_name) const {
        const auto &a = source.fields.at(old_name).declaration;
        const auto &b = source.fields.at(new_name).declaration;
        for (const auto &meaning : {"channel", "quantity", "units", "method", "value_type"})
            require(a.at(meaning).text() == b.at(meaning).text(),
                    "longitudinal field meaning changed: " + std::string(meaning));
        require(tensor(old_name).output == tensor(new_name).output,
                "longitudinal fields require the identical native carrier");
    }
    rfc::database::NativeFieldAnalysis native_analysis(const std::string &name) const {
        require(source.field_sets.contains(name), "unknown native field set");
        const auto &set = source.field_sets.at(name);
        const auto declaration = set.declaration;
        NativeFieldRequest request;
        request.upper_grade =
            rfc::core::Grade{static_cast<std::uint16_t>(declaration.at("grade").integer())};
        for (std::size_t i = 0; i < declaration.at("carriers").size(); ++i)
            request.stalk_carriers.push_back(key(declaration.at("carriers").at(i)));
        for (std::size_t i = 0; i < declaration.at("fields").size(); ++i) {
            request.fields.push_back(source.field_pins.at(
                source.fields.at(declaration.at("fields").at(i).text()).tensor));
        }
        rfc::runtime::ExecutionHostPolicy host;
        host.name="gadget-native-field";
        const auto maximum=std::numeric_limits<std::size_t>::max();
        host.limits={maximum,maximum,maximum,maximum,maximum};
        auto analysis = rfc::database::NativeFieldAnalysis::inspect(
            source.grade_pins.at(request.upper_grade.value), request, set.realization,
            host,{field_memory_headroom()});
        require(analysis.layout().target_coordinates <= coordinate_limit &&
                    analysis.layout().retained_occurrences <= coordinate_limit,
                "native result/provenance coordinate budget");
        return analysis;
    }
    Json local_declaration(const Json &query, bool legacy_budget = true) {
        auto resolved = Json::parse(query.dump()), sections = Json::array();
        const auto name = query.at("context").text();
        auto &prepared = context(name);
        for (std::size_t s = 0; s < query.at("sections").size(); ++s) {
            auto item = Json::parse(query.at("sections").at(s).dump());
            require(item.has("values") != item.has("fields"),
                    "section must declare values or source fields");
            if (item.has("fields")) {
                const auto ids = item.at("fields");
                require(ids.size(), "section fields must be nonempty");
                const auto grade = static_cast<std::uint16_t>(count(
                    item.at("grade"), std::numeric_limits<std::uint16_t>::max(), "section grade"));
                std::vector<Vector> inputs;
                std::set<std::string> seen;
                for (std::size_t f = 0; f < ids.size(); ++f) {
                    const auto id = text(ids.at(f));
                    require(seen.insert(id).second, "duplicate section source field");
                    const auto [v, g] = input(name, id, legacy_budget);
                    require(g == grade, "section source field grade");
                    const auto meaning = source.fields.at(id).declaration;
                    require(meaning.at("units").text() == item.at("units").text() &&
                                meaning.at("quantity").text() == item.at("quantity").text(),
                            "section field meanings differ");
                    inputs.push_back(v);
                }
                for (std::size_t c = 0; c < query.at("stalks").size(); ++c) {
                    const auto stalk = query.at("stalks").at(c);
                    if (stalk.at("grade").integer() == grade)
                        require(stalk.at("axes").dump() == ids.dump(),
                                "source section axes must name their exact fields in order");
                }
                auto values = Json::array();
                for (std::size_t c = 0; c < prepared.to_ambient.at(grade).size(); ++c)
                    for (const auto &v : inputs)
                        values.push(Json::string(v[c].get_str()));
                item.set("values", values);
            }
            sections.push(item);
        }
        resolved.set("sections", sections);
        return resolved;
    }
    Json ablation_query(const Json &query, bool validate) {
        const auto name = query.at("context").text();
        auto &prepared = context(name);
        const auto ids = query.at("fields");
        source_array(ids, "ablation fields");
        require(ids.size(), "ablation requires observed fields");
        const auto first = ids.at(0).text();
        const auto origin = input(name, first, false);
        const rfc::core::Grade grade{origin.second};
        std::map<std::string, Vector> fields;
        for (std::size_t f = 0; f < ids.size(); ++f) {
            const auto id = ids.at(f).text();
            comparable(first, id);
            require(fields.emplace(id, input(name, id, false).first).second,
                    "duplicate measured ablation field");
        }
        (void)source.evidence(query.at("interpretation").at("evidence"));
        const StableSourceIndex index(*prepared.source);
        std::map<std::string, research_detail::SourceAblationRole> roles;
        source_array(query.at("roles"), "scientific ablation roles");
        require(query.at("roles").size(), "scientific roles must be declared");
        for (std::size_t r = 0; r < query.at("roles").size(); ++r) {
            const auto declaration = query.at("roles").at(r);
            members(declaration, {"id", "kind", "meaning", "evidence"}, {"cells", "fields"});
            (void)text(declaration.at("meaning"), source_byte_limit());
            (void)source.evidence(declaration.at("evidence"));
            const auto kind = declaration.at("kind").text();
            research_detail::SourceAblationRole role;
            role.witness = Json::parse(declaration.dump());
            auto binding = Json::array();
            std::set<std::string> seen;
            if (kind == "source_cells") {
                require(declaration.has("cells") && !declaration.has("fields"),
                        "source role requires cells only");
                source_array(declaration.at("cells"), "source role cells");
                require(declaration.at("cells").size(), "source role must retain a cell binding");
                for (std::size_t c = 0; c < declaration.at("cells").size(); ++c) {
                    const auto id = declaration.at("cells").at(c).text();
                    require(source.has_cell(id) && seen.insert(id).second,
                            "unknown or duplicate role cell");
                    const auto cell = index.cell(source.resolve_cell(id));
                    require(cell.has_value(), "role cell is outside the actual context source");
                    role.cells.push_back(*cell);
                    auto record = Json::parse(source.cell_record(source.resolve_cell(id)).dump());
                    record.set("context_coordinate", cardinality(cell->index));
                    binding.push(record);
                }
            } else {
                require(kind == "field_span" && declaration.has("fields") &&
                            !declaration.has("cells"),
                        "field-span role requires fields only");
                source_array(declaration.at("fields"), "role field span");
                require(declaration.at("fields").size(), "field-span role requires a field");
                for (std::size_t f = 0; f < declaration.at("fields").size(); ++f) {
                    const auto id = declaration.at("fields").at(f).text();
                    require(seen.insert(id).second, "duplicate span field");
                    comparable(declaration.at("fields").at(0).text(), id);
                    const auto value = input(name, id, false);
                    require(value.second == grade.value &&
                                tensor(id).output == tensor(first).output,
                            "span requires the identical actual observed field carrier");
                    fields.emplace(id, value.first);
                    role.fields.push_back(id);
                    binding.push(source.fields.at(id).declaration);
                }
            }
            role.witness.set("source_binding", binding);
            require(roles.emplace(text(declaration.at("id"), source_byte_limit()), std::move(role))
                        .second,
                    "duplicate scientific role");
        }
        auto witness = Json::object(), cells = Json::array(), occurrences = Json::array();
        for (std::size_t k = 0; k < prepared.source->grade_count(); ++k) {
            auto grade_cells = Json::array();
            const auto g = rfc::core::Grade{static_cast<std::uint16_t>(k)};
            const auto &arena = prepared.source->grade(g);
            for (std::size_t i = 0; i < arena.cell_count(); ++i) {
                auto record = Json::parse(
                    source.cell_record(stable_cell_key(*prepared.source, {g, static_cast<std::uint32_t>(i)}))
                        .dump());
                record.set("context_coordinate", cardinality(i));
                record.set("ambient_coordinate", cardinality(prepared.to_ambient[k][i]));
                grade_cells.push(record);
            }
            cells.push(grade_cells);
            for (std::size_t i = 0; i < arena.incidence_count(); ++i) {
                auto record =
                    Json::parse(source.occurrence_record({g, arena.occurrence_identities[i]}).dump());
                record.set("exists", Json::boolean(arena.aspects.get(
                                         prepared.source->signature().existence_aspect(), i)));
                record.set("coefficient",
                           Json::string(prepared.source->incidence_coefficient(g, i).get_str()));
                occurrences.push(record);
            }
        }
        witness.set("cells", cells);
        witness.set("occurrences", occurrences);
        Json sheaf;
        if (query.has("sheaf")) {
            sheaf = Json::parse(query.at("sheaf").dump());
            members(sheaf, {"stalks", "maps", "sections", "metric"});
            sheaf.set("operation", Json::string("local_sections"));
            sheaf.set("context", Json::string(name));
            sheaf.set("construction", Json::string("cellular_sheaf"));
            sheaf = local_declaration(sheaf, false);
        }
        return research_detail::source_structure_ablation(
            *prepared.study, query, grade, fields, roles, witness, sheaf, validate, resources);
    }
    Json response_query(const Json &query, bool validate) {
        const auto name = query.at("context").text();
        auto &prepared = context(name);
        const auto ids = research_detail::response_model_field_ids(query);
        require(!ids.empty(), "response query needs declared fields");
        std::map<std::string, Vector> fields;
        std::map<std::string, Json> meanings;
        std::optional<std::uint64_t> output;
        std::optional<rfc::core::Grade> grade;
        const auto endpoint = source.directory->pin_snapshot();
        const auto carriers = temporal_endpoint_detail::carriers(*endpoint);
        for (const auto &id : ids) {
            const auto &t = tensor(id);
            const auto meaning = source.fields.at(id).declaration;
            require(meaning.at("context").text() == name && t.output.has_value(),
                    "response field/context mismatch");
            const auto *carrier = carriers.at(*t.output);
            require(carrier->kind == TemporalCarrierKind::SourceGrade &&
                        carrier->domain.extents().size() == 1 && (!output || *output == *t.output),
                    "responses and factors require the identical actual SourceGrade carrier");
            output = t.output;
            grade = carrier->grade;
            Vector values;
            for (auto coordinate : prepared.to_ambient.at(grade->value)) {
                require(t.tensor.support().contains(coordinate),
                        "response modeling requires complete observed support");
                values.push_back(*t.tensor.state(coordinate));
            }
            fields.emplace(id, std::move(values));
            meanings.emplace(id, meaning);
        }
        (void)source.evidence(query.at("interpretation").at("evidence"));
        for (std::size_t r = 0; r < query.at("responses").size(); ++r)
            (void)source.evidence(query.at("responses").at(r).at("evidence"));
        if (query.has("physical_shape"))
            (void)source.evidence(query.at("physical_shape").at("evidence"));
        if (query.has("dose_comparisons"))
            for (std::size_t d = 0; d < query.at("dose_comparisons").size(); ++d)
                (void)source.evidence(query.at("dose_comparisons").at(d).at("evidence"));
        if (query.has("reference_law"))
            (void)source.evidence(query.at("reference_law").at("evidence"));
        auto coordinates = Json::array();
        for (std::size_t i = 0; i < prepared.to_ambient.at(grade->value).size(); ++i) {
            auto c = Json::parse(source.cell_record(source.coordinate(grade->value,prepared.to_ambient.at(grade->value)[i]))
                                     .dump());
            c.set("context_coordinate", cardinality(i));
            coordinates.push(c);
        }
        return research_detail::observed_response_model(*prepared.study, query, *grade, fields,
                                                        meanings, coordinates, validate, resources);
    }
    Json temporal_coordinates(const PreparedContext &p, std::uint16_t grade) const {
        auto out = Json::array();
        for (auto i : p.to_ambient.at(grade))
            out.push(source.cell_record(source.coordinate(grade,i)));
        return out;
    }
    void temporal_times(const std::string &old_name, const std::string &new_name) const {
        require(old_name != new_name, "temporal contexts must be distinct");
        const auto a = source.contexts.at(old_name).declaration.at("time"),
                   b = source.contexts.at(new_name).declaration.at("time");
        require(a.native().kind()!=NativeRecordKind::Null &&
                    b.native().kind()!=NativeRecordKind::Null,
                "temporal contexts require declared times");
        for (const auto *key : {"axis", "unit", "population"})
            require(a.at(key).text() == b.at(key).text(),
                    "temporal time/population meaning mismatch");
        require(exact_decimal(b.at("position").text()) > exact_decimal(a.at("position").text()),
                "temporal positions must increase");
    }
    void temporal_meanings(const std::string &a, const std::string &b) const {
        for (const auto *key : {"channel", "quantity", "units", "method", "value_type"})
            require(source.fields.at(a).declaration.at(key).text() ==
                        source.fields.at(b).declaration.at(key).text(),
                    "temporal measurement meaning mismatch");
    }
    State &publication(const std::string &id) {
        if (id == "self")
            return *this;
        require(source.publications.contains(id), "unknown study publication: " + id);
        if (!publication_sessions.contains(id))
            publication_sessions.emplace(
                id, std::make_unique<State>(source.publications.at(id), resources, compute));
        return *publication_sessions.at(id);
    }
    scientific_detail::NativeObservationBinding observation_binding(const Json &query) {
        auto &s = publication(text(query.at("publication"))).source;
        const auto operation = text(query.at("operator")), accession = text(query.at("accession"));
        require(s.operator_pins.contains(operation) && s.operator_pins.contains(accession),
                "observation operator or accession declaration unknown");
        source_array(query.at("fields"), "observation fields");
        std::vector<rfc::storage::PinnedNativeAccession> fields;
        auto declarations = Json::object(), meanings = Json::array(), coordinates = Json::array();
        for (std::size_t i = 0; i < query.at("fields").size(); ++i) {
            const auto id = text(query.at("fields").at(i));
            require(s.fields.contains(id), "observation field declaration unknown");
            const auto &f = s.fields.at(id);
            require(f.declaration.at("value_type").text() != "categorical",
                    "linear observations require quantitative or indicator fields");
            fields.push_back(s.field_pins.at(f.tensor));
            meanings.push(s.referenced(f.declaration));
        }
        const auto op = s.operator_pins.at(operation);
        const auto grade = op.carrier()->grade.value;
        require(grade < s.grade_count(), "observation source grade absent");
        for (const auto &cell : s.coordinate_keys(grade))
            coordinates.push(s.referenced(s.cell_record(cell)));
        declarations.set("source_digest", Json::string(s.digest));
        declarations.set("native_sha256", Json::string(s.native_digest));
        auto declaration = Json::object(), documents = Json::array();
        for (const auto *name : {"schema", "id"})
            declaration.set(name, s.declaration.at(name));
        declaration.set("sha256", Json::string(s.declaration_digest));
        for (const auto &[id, document] : s.documents) {
            (void)id;
            documents.push(document.metadata);
        }
        declaration.set("documents", documents);
        declarations.set("declaration", declaration);
        declarations.set("evidence_profile", Json::string("document-bound-byte-spans-v1"));
        declarations.set("operator", s.referenced(s.operators.at(operation)));
        declarations.set("accession", s.referenced(s.operators.at(accession)));
        declarations.set("fields", meanings);
        declarations.set("source_coordinates", coordinates);
        return {op, s.operator_pins.at(accession), std::move(fields), declarations};
    }
    scientific_detail::NativeNoveltyBinding novelty_binding(const Json &query,
                                                            bool metadata = false) {
        const auto &old_decl = query.at("old"), &new_decl = query.at("new");
        auto &old_source = publication(text(old_decl.at("publication"))).source;
        auto &new_source = publication(text(new_decl.at("publication"))).source;
        const auto field = [](const ScientificStudy::State &s, const Json &id) {
            const auto name = text(id);
            require(s.fields.contains(name), "novelty field declaration unknown");
            const auto &f = s.fields.at(name);
            require(f.declaration.at("value_type").text() != "categorical",
                    "novelty requires quantitative or indicator Field amplitudes");
            return s.field_pins.at(f.tensor);
        };
        const auto old = field(old_source, old_decl.at("field")),
                   current = field(new_source, new_decl.at("field"));
        const auto &old_meaning = old_source.fields.at(text(old_decl.at("field"))).declaration;
        const auto &new_meaning = new_source.fields.at(text(new_decl.at("field"))).declaration;
        const auto same_meaning = [&](const Json &meaning, bool change_field = false) {
            for (const auto *name : {"quantity", "units", "channel", "value_type"}) {
                // A signed change of an indicator can be quantitative. Its
                // declared origin/quantity remains distinct; units and GET
                // channel still have to agree with the target innovation.
                if (change_field && (std::string_view(name) == "quantity" ||
                                     std::string_view(name) == "value_type"))
                    continue;
                require(meaning.at(name).text() == new_meaning.at(name).text(),
                        "novelty field meanings differ");
            }
        };
        same_meaning(old_meaning);
        const auto metric_id = text(new_decl.at("metric"));
        require(new_source.operator_pins.contains(metric_id), "novelty metric declaration unknown");
        std::vector<rfc::storage::PinnedNativeAccession> causes;
        auto meanings = Json::object();
        if (query.has("causes")) {
            source_array(query.at("causes"), "novelty cause word");
            for (std::size_t i = 0; i < query.at("causes").size(); ++i) {
                const auto id = text(query.at("causes").at(i));
                causes.push_back(field(new_source, query.at("causes").at(i)));
                const auto &meaning = new_source.fields.at(id).declaration;
                same_meaning(meaning, true);
                if (metadata && !meanings.has(id))
                    meanings.set(id, new_source.referenced(meaning));
            }
        }
        (void)source.evidence(query.at("comparison").at("evidence"));
        auto declarations = Json::object();
        if (metadata) {
            const auto describe = [](const ScientificStudy::State &s,
                                     const rfc::storage::PinnedNativeAccession &pin,
                                     const Json &meaning) {
                auto out = Json::object(), documents = Json::array(), coordinates = Json::array();
                out.set("source_digest", Json::string(s.digest));
                out.set("native_sha256", Json::string(s.native_digest));
                out.set("declaration_sha256", Json::string(s.declaration_digest));
                out.set("study_id", s.declaration.at("id"));
                out.set("field", s.referenced(meaning));
                out.set("carrier", s.carrier_metadata(*pin.carrier()));
                for (const auto &[id, doc] : s.documents) {
                    (void)id;
                    documents.push(doc.metadata);
                }
                out.set("documents", documents);
                for (const auto &cell : s.coordinate_keys(pin.carrier()->grade.value))
                    coordinates.push(s.referenced(s.cell_record(cell)));
                out.set("source_coordinates", coordinates);
                return out;
            };
            declarations.set("old", describe(old_source, old, old_meaning));
            auto current_declaration = describe(new_source, current, new_meaning);
            current_declaration.set("metric",
                                    new_source.referenced(new_source.operators.at(metric_id)));
            declarations.set("new", current_declaration);
            declarations.set("cause_fields", meanings);
            declarations.set("evidence_profile", Json::string("document-bound-byte-spans-v1"));
        }
        return {old, current, new_source.operator_pins.at(metric_id), std::move(causes),
                declarations};
    }
    scientific_detail::NativeHodgeBinding hodge_binding(const Json &query) {
        auto &s = publication(text(query.at("publication"))).source;
        source_array(query.at("fields"), "Hodge fields");
        require(query.at("fields").size() > 0, "Hodge requires at least one Field");
        std::vector<rfc::storage::PinnedNativeAccession> fields, metrics;
        auto declarations = Json::object(), field_meanings = Json::array(),
             metric_meanings = Json::array(), coordinates = Json::array();
        for (std::size_t i = 0; i < query.at("fields").size(); ++i) {
            const auto id = text(query.at("fields").at(i));
            require(s.fields.contains(id), "Hodge field declaration unknown");
            const auto &f = s.fields.at(id);
            require(f.declaration.at("value_type").text() != "categorical",
                    "Hodge requires quantitative or indicator Field amplitudes");
            fields.push_back(s.field_pins.at(f.tensor));
            field_meanings.push(f.declaration);
        }
        const auto selected = query.at("metrics");
        if (selected.has("tensors")) {
            source_array(selected.at("tensors"), "Hodge metric tensors");
            for (std::size_t k = 0; k < selected.at("tensors").size(); ++k) {
                const auto id = text(selected.at("tensors").at(k));
                require(s.operator_pins.contains(id), "Hodge metric declaration unknown");
                metrics.push_back(s.operator_pins.at(id));
                metric_meanings.push(s.operators.at(id));
            }
        }
        for (const auto &grade : s.coordinate_grades()) {
            auto cells = Json::array();
            for (const auto &cell : grade)
                cells.push(s.cell_record(cell));
            coordinates.push(cells);
        }
        declarations.set("source_digest", Json::string(s.digest));
        declarations.set("native_sha256", Json::string(s.native_digest));
        declarations.set("declaration", s.declaration);
        declarations.set("fields", field_meanings);
        declarations.set("metric_tensors", metric_meanings);
        declarations.set("grade_coordinates", coordinates);
        return {fields.front(), std::move(fields), std::move(metrics), declarations};
    }
    scientific_detail::DeclaredResponseBinding response_binding(const Json &query) {
        auto declarations = Json::object();
        struct EndpointPins {
            rfc::storage::PinnedNativeAccession L, I, u;
            std::optional<rfc::storage::PinnedNativeAccession> A;
        };
        auto bind = [&](const char *name) -> EndpointPins {
            const auto d = query.at(name);
            auto &s = publication(text(d.at("publication"))).source;
            const auto resolving = text(d.at("resolving")), realization = text(d.at("realization")),
                       field = text(d.at("field"));
            require(s.operator_pins.contains(resolving) && s.operator_pins.contains(realization),
                    "response operator declaration unknown");
            require(s.fields.contains(field), "response field declaration unknown");
            require(s.fields.at(field).declaration.at("value_type").text() != "categorical",
                    "response fields require quantitative or indicator coordinates");
            const auto L = s.operator_pins.at(resolving), I = s.operator_pins.at(realization),
                       u = s.field_pins.at(s.fields.at(field).tensor);
            std::optional<rfc::storage::PinnedNativeAccession> A;
            auto record = Json::object(), document_records = Json::array();
            record.set("source_digest", Json::string(s.digest));
            record.set("native_sha256", Json::string(s.native_digest));
            record.set("declaration", s.declaration);
            for (const auto &[id, doc] : s.documents) {
                (void)id;
                document_records.push(doc.metadata);
            }
            record.set("documents", document_records);
            record.set("field", s.fields.at(field).declaration);
            record.set("resolving", s.operators.at(resolving));
            record.set("realization", s.operators.at(realization));
            if (query.has("observation")) {
                const auto accession = text(query.at("observation").at(name));
                require(s.operator_pins.contains(accession),
                        "response accession declaration unknown");
                A = s.operator_pins.at(accession);
                record.set("accession", s.operators.at(accession));
            }
            auto input_values = Vector(u.tensor()->domain().logical_size());
            for (const auto &e : u.tensor()->values().entries())
                input_values[e.index] = e.value;
            auto values = Json::array();
            for (const auto &v : input_values)
                values.push(Json::string(v.get_str()));
            record.set("input_values", values);
            auto carriers = Json::array();
            for (const auto &c : u.endpoint()->carriers) {
                auto item = Json::object(), coordinates = Json::array();
                item.set("key", Json::string(std::to_string(c.key)));
                item.set("kind", Json::integer(static_cast<std::int64_t>(c.kind)));
                item.set("coordinate_system", Json::string(std::to_string(c.coordinate_system)));
                if (c.kind == TemporalCarrierKind::SourceGrade) {
                    item.set("grade", Json::integer(c.grade.value));
                    for (std::size_t i = 0; i < c.domain.logical_size(); ++i) {
                        auto cell = stable_cell_key(u.endpoint()->source,
                                                    {c.grade, static_cast<std::uint32_t>(i)});
                        coordinates.push(s.cell_record(cell));
                    }
                } else {
                    for (std::size_t i = 0; i < c.domain.logical_size(); ++i)
                        coordinates.push(cardinality(i));
                }
                item.set("coordinates", coordinates);
                carriers.push(item);
            }
            record.set("carriers", carriers);
            declarations.set(name, record);
            return {L, I, u, A};
        };
        const auto old = bind("old"), current = bind("new");
        scientific_detail::DeclaredResponseBinding b{
            old.L, current.L, old.I, current.I, old.u, current.u, old.A, current.A, {}, Json{}};
        for (const auto *name : {"quantity", "units", "channel", "value_type"})
            require(declarations.at("old").at("field").at(name).text() ==
                        declarations.at("new").at("field").at(name).text(),
                    "response field meanings differ across publications");
        for (const auto *name : {"ambient", "source"})
            (void)source.evidence(query.at("comparisons").at(name).at("evidence"));
        if (query.has("observation"))
            (void)source.evidence(query.at("observation").at("comparison").at("evidence"));
        if (query.has("pairing")) {
            auto &s = publication(text(query.at("new").at("publication"))).source;
            const auto id = text(query.at("pairing").at("metric"));
            require(s.operator_pins.contains(id), "response metric declaration unknown");
            b.metric = s.operator_pins.at(id);
            auto record = declarations.at("new");
            record.set("metric", s.operators.at(id));
            declarations.set("new", record);
        }
        b.publications = declarations;
        return b;
    }
    static rfc::numeric::ExactMatrix
    temporal_map_between(State &old_owner, const std::string &old_name, State &new_owner,
                         const std::string &new_name, std::uint16_t grade, const Json &declaration,
                         const State &evidence_owner, bool shared_identity) {
        using Matrix = rfc::numeric::ExactMatrix;
        auto &a = old_owner.context(old_name), &b = new_owner.context(new_name);
        const auto &ac = a.to_ambient.at(grade), &bc = b.to_ambient.at(grade);
        members(declaration, {"grade", "kind", "meaning", "evidence"}, {"entries"});
        require(count(declaration.at("grade"), a.to_ambient.size() - 1, "map grade") == grade,
                "map grade differs from carrier");
        (void)text(declaration.at("meaning"), source_byte_limit());
        (void)evidence_owner.source.evidence(declaration.at("evidence"));
        Matrix j(bc.size(), ac.size());
        const auto kind = declaration.at("kind").text();
        if (kind == "retained_identity") {
            require(shared_identity, "retained identity requires declared shared source lineage");
            require(!declaration.has("entries"), "retained identity cannot contain matrix entries");
            std::map<StableCellKey, std::size_t> target;
            for (std::size_t i = 0; i < bc.size(); ++i)
                target.emplace(new_owner.source.coordinate(grade,bc[i]), i);
            for (std::size_t i = 0; i < ac.size(); ++i) {
                auto key_ = old_owner.source.coordinate(grade,ac[i]);
                if (target.contains(key_))
                    j.set(target.at(key_), i, Rational{1});
            }
        } else {
            require(kind == "explicit" && declaration.has("entries"), "temporal map kind/entries");
            std::map<std::string, std::size_t> rows, cols;
            for (std::size_t i = 0; i < bc.size(); ++i)
                rows.emplace(
                    new_owner.source.cell_record(new_owner.source.coordinate(grade,bc[i]))
                        .at("id")
                        .text(),
                    i);
            for (std::size_t i = 0; i < ac.size(); ++i)
                cols.emplace(
                    old_owner.source.cell_record(old_owner.source.coordinate(grade,ac[i]))
                        .at("id")
                        .text(),
                    i);
            source_array(declaration.at("entries"), "temporal map entries");
            std::set<std::pair<std::size_t, std::size_t>> seen;
            for (std::size_t i = 0; i < declaration.at("entries").size(); ++i) {
                auto e = declaration.at("entries").at(i);
                members(e, {"row", "column", "value"});
                require(rows.contains(e.at("row").text()) && cols.contains(e.at("column").text()),
                        "map entry does not name its actual source carriers");
                auto r = rows.at(e.at("row").text()), c = cols.at(e.at("column").text());
                require(seen.emplace(r, c).second, "duplicate temporal map entry");
                j.set(r, c, exact_rational(e.at("value").text()));
            }
        }
        return j;
    }
    rfc::numeric::ExactMatrix temporal_map(const std::string &old_name, const std::string &new_name,
                                           std::uint16_t grade, const Json &declaration) {
        return temporal_map_between(*this, old_name, *this, new_name, grade, declaration, *this,
                                    true);
    }
    Json temporal_binding(const std::string &id) const {
        auto out = Json::object();
        out.set("id", Json::string(id));
        out.set("source_digest", Json::string(source.digest));
        out.set("declaration_sha256", Json::string(source.declaration_digest));
        out.set("native_sha256", Json::string(source.native_digest));
        out.set("source_namespace", source.declaration.at("native").at("source_namespace"));
        out.set("source_sequence", Json::string(std::to_string(source.directory->sequence())));
        if (source.declaration.at("native").has("cutoff"))
            out.set("native_history", source.catalog().at("native_history"));
        out.set("design", source.declaration.at("design"));
        return out;
    }
    std::pair<Vector, std::uint16_t> temporal_input(const std::string &name, const std::string &id,
                                                    const Json &e, Json &observations) {
        const auto &t = tensor(id);
        const auto meaning = source.fields.at(id).declaration;
        require(meaning.at("context").text() == name &&
                    meaning.at("value_type").text() != "categorical",
                "temporal field/context or measurement type mismatch");
        auto carrier =
            temporal_endpoint_detail::carriers(*source.directory->pin_snapshot()).at(*t.output);
        require(carrier->kind == TemporalCarrierKind::SourceGrade &&
                    carrier->domain.extents().size() == 1,
                "history requires a scalar SourceGrade carrier");
        auto &p = context(name);
        const auto grade = carrier->grade.value;
        std::map<std::string, Rational> supplied;
        if (e.has("extension")) {
            auto ext = e.at("extension");
            members(ext, {"meaning", "evidence", "values"});
            (void)text(ext.at("meaning"), source_byte_limit());
            (void)source.evidence(ext.at("evidence"));
            source_array(ext.at("values"), "declared extension values");
            for (std::size_t i = 0; i < ext.at("values").size(); ++i) {
                auto v = ext.at("values").at(i);
                members(v, {"entity", "value"});
                require(
                    supplied.emplace(v.at("entity").text(), exact_rational(v.at("value").text()))
                        .second,
                    "duplicate declared extension entity");
            }
        }
        Vector x;
        observations = Json::array();
        for (auto coordinate : p.to_ambient.at(grade)) {
            const auto entity =
                source.cell_record(source.coordinate(grade,coordinate)).at("id").text();
            const auto present = t.tensor.support().contains(coordinate);
            require(!present || !supplied.contains(entity),
                    "declared extension cannot overwrite an observed coordinate");
            require(present || supplied.contains(entity),
                    "partial temporal observations require an explicit declared extension");
            auto record = Json::object();
            record.set("entity", Json::string(entity));
            record.set("ambient_coordinate", cardinality(coordinate));
            record.set("reported", Json::boolean(present));
            record.set("observed_value", present ? rational(*t.tensor.state(coordinate)) : Json{});
            const auto value = present ? *t.tensor.state(coordinate) : supplied.at(entity);
            if (meaning.at("value_type").text() == "indicator")
                require(value == 0 || value == 1, "indicator extension requires zero or one");
            record.set("execution_value", rational(value));
            observations.push(record);
            x.push_back(value);
            supplied.erase(entity);
        }
        require(supplied.empty(), "extension entity is outside the actual temporal carrier");
        return {std::move(x), grade};
    }
    Json history_query(const Json &query, bool validate) {
        members(query, {"operation", "endpoints", "correspondences", "word", "lambda"},
                {"lineages", "map_output"});
        source_array(query.at("endpoints"), "history endpoints");
        source_array(query.at("correspondences"), "history correspondences");
        require(query.at("endpoints").size() >= 2 &&
                    query.at("correspondences").size() == query.at("endpoints").size() - 1,
                "history endpoint/map count mismatch");
        if (query.has("lineages")) {
            source_array(query.at("lineages"), "history source lineages");
            require(query.at("lineages").size() == query.at("correspondences").size(),
                    "history lineage/step count mismatch");
        }
        std::vector<research_detail::TemporalReading> readings;
        std::vector<research_detail::TemporalMaps> maps;
        std::vector<State *> owners;
        std::vector<std::string> names;
        for (std::size_t t = 0; t < query.at("endpoints").size(); ++t) {
            auto e = query.at("endpoints").at(t);
            members(e, {"context", "field"}, {"publication", "extension"});
            auto pub = e.has("publication") ? e.at("publication").text() : "self";
            auto &owner_ = publication(pub);
            auto name = e.at("context").text(), id = e.at("field").text();
            Json observations;
            auto [x, k] = owner_.temporal_input(name, id, e, observations);
            auto &p = owner_.context(name);
            if (!validate)
                (void)owner_.temporal(p);
            auto d = Json::parse(e.dump());
            d.set("publication_binding", owner_.temporal_binding(pub));
            d.set("time", owner_.source.contexts.at(name).declaration.at("time"));
            require(!json_object_is_type(d.at("time").raw(), json_type_null),
                    "temporal endpoints require declared times");
            auto source_grades = Json::array();
            for (std::size_t g = 0; g < p.to_ambient.size(); ++g)
                source_grades.push(owner_.temporal_coordinates(p, static_cast<std::uint16_t>(g)));
            d.set("source_grades", source_grades);
            d.set("metrics", owner_.source.contexts.at(name).declaration.at("metrics"));
            d.set("field_semantics", owner_.source.fields.at(id).declaration);
            d.set("observations", observations);
            auto occurrences = Json::array();
            for (std::size_t g = 1; g < p.source->grade_count(); ++g) {
                auto grade = rfc::core::Grade{static_cast<std::uint16_t>(g)};
                const auto &arena = p.source->grade(grade);
                for (std::size_t i = 0; i < arena.incidence_count(); ++i) {
                    auto record = Json::parse(
                        owner_.source.occurrence_record({grade, arena.occurrence_identities[i]})
                            .dump());
                    record.set("parent_key",
                               Json::string(std::to_string(
                                   arena.relation_identities[arena.upper_index[i]].value)));
                    record.set(
                        "lower_key",
                        Json::string(std::to_string(
                            stable_cell_key(*p.source,
                                            {rfc::core::Grade{static_cast<std::uint16_t>(g - 1)},
                                             arena.lower_index[i]})
                                .key)));
                    record.set("coefficient",
                               Json::string(p.source->incidence_coefficient(grade, i).get_str()));
                    record.set("lattice_coefficient",
                               Json::string(arena.lattice_coefficient[i].get_str()));
                    record.set("parent_scale", rational(arena.cell_scales[arena.upper_index[i]]));
                    record.set(
                        "lower_scale",
                        rational(
                            p.source->grade(rfc::core::Grade{static_cast<std::uint16_t>(g - 1)})
                                .cell_scales[arena.lower_index[i]]));
                    auto aspects = Json::array();
                    for (std::size_t a = 0; a < p.source->signature().aspects().size(); ++a)
                        aspects.push(Json::boolean(
                            arena.aspects.get(static_cast<rfc::composite::AspectId>(a), i)));
                    record.set("aspects", aspects);
                    record.set("evidence", owner_.source.evidence(record.at("evidence")));
                    occurrences.push(record);
                }
            }
            d.set("source_occurrences", occurrences);
            std::map<std::size_t, Vector> diagonal;
            auto diagonal_meanings = Json::array();
            auto current = k;
            for (std::size_t i = 0; i < query.at("word").size(); ++i) {
                auto op = query.at("word").at(i);
                const auto kind = op.at("kind").text();
                if (kind == "diagonal_field") {
                    source_array(op.at("fields"), "word diagonal fields");
                    require(op.at("fields").size() == query.at("endpoints").size(),
                            "word factor field/endpoint count mismatch");
                    auto id_ = op.at("fields").at(t).text();
                    auto [v, g] = owner_.input(name, id_, false);
                    require(
                        g == current &&
                            owner_.source.fields.at(id_).declaration.at("units").text() ==
                                "dimensionless",
                        "word diagonal field requires current grade and dimensionless multipliers");
                    diagonal.emplace(i, std::move(v));
                    diagonal_meanings.push(owner_.source.fields.at(id_).declaration);
                } else if (kind == "boundary") {
                    require(current > 0, "word boundary at grade zero");
                    --current;
                } else if (kind == "adjoint") {
                    require(current < p.source->max_grade().value,
                            "word adjoint has no upper grade");
                    ++current;
                }
            }
            d.set("diagonal_field_semantics", diagonal_meanings);
            bool shared = true;
            if (t) {
                auto &old_owner = *owners.back();
                const auto old_time = readings.back().declaration.at("time"),
                           new_time = d.at("time");
                for (const auto *key_ : {"axis", "unit"})
                    require(old_time.at(key_).text() == new_time.at(key_).text(),
                            "temporal time-axis/unit mismatch");
                require(exact_decimal(new_time.at("position").text()) >
                            exact_decimal(old_time.at("position").text()),
                        "temporal positions must increase");
                require(&old_owner != &owner_ || names.back() != name,
                        "temporal contexts must be distinct within a publication");
                const auto old_id = query.at("endpoints").at(t - 1).at("field").text();
                for (const auto *key_ : {"channel", "quantity", "units", "method", "value_type"})
                    require(old_owner.source.fields.at(old_id).declaration.at(key_).text() ==
                                owner_.source.fields.at(id).declaration.at(key_).text(),
                            "temporal measurement meaning mismatch");
                require(&old_owner == &owner_ || query.has("lineages"),
                        "cross-publication history requires declared source lineage");
                if (query.has("lineages")) {
                    auto lineage = query.at("lineages").at(t - 1);
                    members(lineage, {"kind", "population_relation", "meaning", "evidence"});
                    (void)text(lineage.at("meaning"), source_byte_limit());
                    (void)source.evidence(lineage.at("evidence"));
                    const auto kind = lineage.at("kind").text();
                    require(kind == "shared_native_namespace" || kind == "independent_sources",
                            "unknown source lineage kind");
                    shared = kind == "shared_native_namespace";
                    if (shared)
                        require(old_owner.source.declaration.at("native")
                                        .at("source_namespace")
                                        .text() == owner_.source.declaration.at("native")
                                                       .at("source_namespace")
                                                       .text(),
                                "shared source lineage requires the same native namespace");
                    const auto population = lineage.at("population_relation").text();
                    require(population == "same_study_population" ||
                                population == "declared_population_change",
                            "unknown population relation");
                    if (population == "same_study_population")
                        require(old_time.at("population").text() ==
                                    new_time.at("population").text(),
                                "same population lineage has different population declarations");
                    lineage = Json::parse(lineage.dump());
                    lineage.set("evidence", source.evidence(lineage.at("evidence")));
                    d.set("lineage", lineage);
                } else
                    require(old_time.at("population").text() == new_time.at("population").text(),
                            "temporal population meaning mismatch");
                require(readings.back().grade == rfc::core::Grade{k},
                        "history field grade changed");
                require(old_owner.context(names.back()).source->grade_count() ==
                            p.source->grade_count(),
                        "history requires the same declared grade family");
                auto declarations = query.at("correspondences").at(t - 1);
                source_array(declarations, "step source-grade maps");
                require(declarations.size() == p.source->grade_count(),
                        "history requires each source-grade correspondence");
                research_detail::TemporalMaps step;
                for (std::size_t grade = 0; grade < declarations.size(); ++grade)
                    step.push_back(temporal_map_between(old_owner, names.back(), owner_, name,
                                                        static_cast<std::uint16_t>(grade),
                                                        declarations.at(grade), *this, shared));
                maps.push_back(std::move(step));
            }
            readings.push_back({p.study.get(), rfc::core::Grade{k}, std::move(x), d,
                                owner_.temporal_coordinates(p, k), std::move(diagonal), shared,
                                p.temporal.get()});
            names.push_back(name);
            owners.push_back(&owner_);
        }
        return research_detail::temporal_field_history(query, readings, maps, validate);
    }
    Json censor_query(const Json &query, bool validate) {
        auto name = query.at("context").text();
        auto &p = context(name);
        std::set<std::string> ids;
        for (const auto *key : {"lower", "upper", "lower_closed", "upper_closed", "upper_finite"})
            ids.insert(query.at(key).text());
        for (const auto *key : {"selection", "strata"})
            for (std::size_t i = 0; i < query.at(key).size(); ++i)
                ids.insert(query.at(key).at(i).text());
        std::map<std::string, Vector> fields;
        std::map<std::string, Json> meanings;
        std::vector<bool> upper_support;
        std::optional<std::uint64_t> carrier;
        for (const auto &id : ids) {
            const auto &t = tensor(id);
            const auto c =
                temporal_endpoint_detail::carriers(*source.directory->pin_snapshot()).at(*t.output);
            require(source.fields.at(id).declaration.at("context").text() == name &&
                        c->kind == TemporalCarrierKind::SourceGrade && c->grade.value == 0 &&
                        c->domain.rank() == 1,
                    "censored observations require the actual context grade-zero carrier");
            require(!carrier || *carrier == *t.output,
                    "censored fields require identical native carriers");
            carrier = *t.output;
            Vector x;
            for (auto coordinate : p.to_ambient[0]) {
                bool present = t.tensor.support().contains(coordinate);
                if (id == query.at("upper").text())
                    upper_support.push_back(present);
                if (id != query.at("upper").text() || id == query.at("lower").text())
                    require(present, "censored endpoint metadata has missing support");
                x.push_back(present ? *t.tensor.state(coordinate) : Rational{0});
            }
            fields.emplace(id, std::move(x));
            meanings.emplace(id, source.fields.at(id).declaration);
        }
        require(query.at("population").text() ==
                    source.declaration.at("design").at("population").text(),
                "censored population differs from declared study population");
        (void)source.evidence(query.at("evidence"));
        if (query.has("joint_family")) {
            const auto family = query.at("joint_family");
            (void)source.evidence(family.at("evidence"));
            if (family.has("constraints"))
                for (std::size_t i = 0; i < family.at("constraints").size(); ++i)
                    (void)source.evidence(family.at("constraints").at(i).at("evidence"));
        }
        return research_detail::censored_endpoint_bounds(
            *p.study, query, fields, upper_support, meanings, temporal_coordinates(p, 0), validate);
    }
    void admit(const Json &query) {
        const auto operation = query.at("operation").text();
        if (operation == "metric_relative_novelty") {
            scientific_detail::NativeNoveltyPlans::validate(query, novelty_binding(query));
            return;
        }
        if (operation == "linear_response_observation") {
            scientific_detail::NativeObservationPlans::validate(query, observation_binding(query));
            return;
        }
        if (operation == "hodge_frame_decomposition") {
            scientific_detail::NativeHodgePlans::validate(query, hodge_binding(query));
            return;
        }
        if (operation == "declared_response_transition") {
            scientific_detail::DeclaredResponsePlans::validate(query, response_binding(query));
            return;
        }
        if (operation == "temporal_field_history") {
            (void)history_query(query, true);
            return;
        }
        if (operation == "censored_endpoint_bounds") {
            (void)censor_query(query, true);
            return;
        }
        if (operation == "observed_response_model" || operation == "constrained_response_model") {
            (void)response_query(query, true);
            return;
        }
        if (operation == "source_structure_ablation") {
            (void)ablation_query(query, true);
            return;
        }
        if (operation == "chain_certificate") {
            (void)analyze_chain(*context(query.at("context").text()).study,
                                chain_declaration(query), true);
            return;
        }
        if (operation == "realized_path_comparison") {
            (void)compare_realized_paths(*context(query.at("context").text()).study, query, true);
            return;
        }
        if (operation == "local_sections") {
            (void)analyze_local_sections(*context(query.at("context").text()).study,
                                         local_declaration(query), true);
            return;
        }
        if (operation == "tensor_window") {
            members(query, {"operation", "field", "lower", "upper", "available_sequence"},
                    {"preparation"});
            if (query.has("preparation"))
                require(query.at("preparation").text() == "reuse" ||
                            query.at("preparation").text() == "one_shot", "window preparation policy");
            const auto &t = tensor(text(query.at("field")));
            require(query.at("lower").size() == t.tensor.domain().rank() &&
                        query.at("upper").size() == t.tensor.domain().rank(),
                    "window rank");
            for (std::size_t a = 0; a < t.tensor.domain().rank(); ++a) {
                const auto lo =
                    count(query.at("lower").at(a), t.tensor.domain().extents()[a], "window lower");
                const auto hi =
                    count(query.at("upper").at(a), t.tensor.domain().extents()[a], "window upper");
                require(lo <= hi, "window bounds");
            }
            const auto pin =
                source.field_pins.at(source.fields.at(query.at("field").text()).tensor);
            rfc::runtime::admit_work_availability(
                pin.source_namespace(), pin.sequence(),
                rfc::runtime::WorkAvailability{pin.source_namespace(),
                                               publication_cutoff(query.at("available_sequence"))});
            std::size_t size = 1;
            for (std::size_t a = 0; a < t.tensor.domain().rank(); ++a)
                size = rfc::machine::checked_size_mul(
                    size, static_cast<std::size_t>(query.at("upper").at(a).integer() -
                                                   query.at("lower").at(a).integer()));
            require(size <= coordinate_limit, "window output budget");
            return;
        }
        if (operation == "green_response") {
            members(query, {"operation", "context", "field", "lambda", "interior", "coordinates"});
            const auto [x, grade] = input(query.at("context").text(), query.at("field").text());
            require(x.size() <= 128, "partitioned Green coordinate budget");
            const auto declaration = query.at("coordinates");
            members(declaration, {"id", "cells"});
            const auto coordinate_id = key(declaration.at("id"));
            require(declaration.at("cells").size() == x.size(), "Green coordinate witness size");
            const auto &prepared = context(query.at("context").text());
            for (std::size_t i = 0; i < x.size(); ++i)
                require(key(declaration.at("cells").at(i)) ==
                            source.coordinate(grade,prepared.to_ambient.at(grade).at(i))
                                .key,
                        "Green coordinate witness differs from source order");
            const auto binding = std::to_string(grade) + "|" + declaration.at("cells").dump();
            require(!green_spaces.contains(coordinate_id) ||
                        green_spaces.at(coordinate_id) == binding,
                    "Green coordinate identity reused with a different basis");
            require(green_spaces.contains(coordinate_id) || green_spaces.size() < 128,
                    "Green coordinate declaration budget");
            green_spaces.insert_or_assign(coordinate_id, binding);
            require(exact_decimal(query.at("lambda").text()) > 0,
                    "positive Green parameter required");
            std::set<std::size_t> indices;
            require(query.at("interior").size() <= 128, "interior partition budget");
            for (std::size_t i = 0; i < query.at("interior").size(); ++i) {
                const auto c =
                    count(query.at("interior").at(i), x.size() - 1, "interior coordinate");
                require(indices.insert(c).second, "duplicate interior coordinate");
            }
            return;
        }
        if (operation == "structural_reference_calibration") {
            (void)cpu();
            (void)structural_calibration_query(query, true);
            return;
        }
        if (operation == "relational_field_ensemble") {
            (void)cpu();
            (void)ensemble_query(query, true);
            return;
        }
        if (is_research_operation(operation)) {
            if (operation == "paired_field_contrasts" || operation == "controlled_pairings" ||
                operation == "linear_model_comparison" ||
                operation == "finite_reference_calibration")
                (void)cpu();
            if (operation == "paired_field_contrasts")
                (void)source.evidence(query.at("normalization").at("evidence"));
            if (operation == "finite_reference_calibration")
                (void)source.evidence(query.at("law").at("evidence"));
            (void)execute_research_query(*owner, query, true, resources);
            return;
        }
        if (operation == "support_mask") {
            members(query, {"operation", "fields", "mode", "backend"},
                    {"workers", "tile_words", "gpu_words"});
            const auto fields = query.at("fields");
            require(json_object_is_type(fields.raw(), json_type_array) && fields.size() &&
                        fields.size() <= source.fields.size(),
                    "support mask requires retained fields");
            std::set<std::string> ids;
            const auto &first = tensor(fields.at(0).text());
            const auto carrier =
                temporal_endpoint_detail::carriers(*source.directory->pin_snapshot())
                    .at(*first.output);
            require(carrier->kind == TemporalCarrierKind::SourceGrade &&
                        carrier->domain.extents().size() == 1,
                    "support mask requires a rank-one SourceGrade carrier");
            require(first.tensor.domain().logical_size() <= coordinate_limit,
                    "support mask coordinate budget");
            for (std::size_t i = 0; i < fields.size(); ++i) {
                const auto id = fields.at(i).text();
                require(ids.insert(id).second, "duplicate support field");
                require(tensor(id).output == first.output,
                        "support fields require the identical native carrier");
            }
            const auto mode = query.at("mode").text();
            require(mode == "all_reported" || mode == "any_reported" || mode == "none_reported",
                    "unknown support mode");
            const auto backend = query.at("backend").text();
            require(backend == "cpu_tiles" || backend == "vulkan" || backend == "cpu_vulkan" ||
                        backend == "measured",
                    "unknown packed backend");
#if !GADGET_HAS_VULKAN
            require(backend == "cpu_tiles", "Vulkan compute is disabled in this build");
#endif
            for (const auto &option : {"workers", "tile_words"})
                if (query.has(option))
                    require(count(query.at(option),
                                  std::string(option) == "workers" ? cpu().coordinator().capacity()
                                                                   : coordinate_limit,
                                  option) > 0,
                            "packed worker/tile count must be positive");
            require(!query.has("workers") || backend == "cpu_tiles" || backend == "cpu_vulkan",
                    "workers only apply to explicit CPU work");
            require(!query.has("gpu_words") || backend == "cpu_vulkan",
                    "gpu_words only applies to cpu_vulkan");
            if (query.has("gpu_words"))
                (void)count(
                    query.at("gpu_words"),
                    rfc::machine::PackedBits::words_for_size(first.tensor.domain().logical_size()),
                    "GPU words");
            return;
        }
        if (operation == "inspect") {
            members(query, {"operation"});
            return;
        }
        if (operation == "coverage") {
            members(query, {"operation", "field"});
            require(tensor(query.at("field").text()).tensor.domain().logical_size() <=
                        coordinate_limit,
                    "coverage coordinate budget");
            return;
        }
        if (operation == "field_summary") {
            members(query, {"operation", "field", "support"});
            require(query.at("support").text() == "reported_only",
                    "summary requires explicit reported_only support");
            (void)tensor(query.at("field").text());
            const auto &field = source.fields.at(query.at("field").text());
            require(field.declaration.at("value_type").text() != "categorical",
                    "categorical codes cannot be linear summary amplitudes");
            require(tensor(query.at("field").text()).tensor.domain().logical_size() <=
                        coordinate_limit,
                    "summary coordinate budget");
            return;
        }
        if (operation == "provenance") {
            members(query, {"operation", "entity"});
            require(source.has_cell(query.at("entity").text()),
                    "unknown declared entity identity");
            return;
        }
        if (operation == "structure_transition") {
            members(query,
                    {"operation", "context", "other_context", "field", "new_field", "lambda"},
                    {"map_output"});
            (void)research_detail::FactoredTemporalMaps::requested(query);
            const auto a = query.at("context").text(), b = query.at("other_context").text();
            const auto old = input(a, query.at("field").text()),
                       next = input(b, query.at("new_field").text());
            require(old.second == next.second, "matching declared structure grade");
            require(
                source.fields.at(query.at("field").text()).declaration.at("units").text() ==
                    source.fields.at(query.at("new_field").text()).declaration.at("units").text(),
                "structure field units");
            require(source.fields.at(query.at("field").text()).declaration.at("quantity").text() ==
                        source.fields.at(query.at("new_field").text())
                            .declaration.at("quantity")
                            .text(),
                    "structure field quantity");
            require(exact_decimal(query.at("lambda").text()) > 0,
                    "structure lambda must be positive");
            const auto t0 = source.contexts.at(a).declaration.at("time"),
                       t1 = source.contexts.at(b).declaration.at("time");
            require(t0.has("axis") && t1.has("axis") &&
                        t0.at("axis").text() == t1.at("axis").text() &&
                        t0.at("unit").text() == t1.at("unit").text() &&
                        exact_decimal(t1.at("position").text()) >
                            exact_decimal(t0.at("position").text()),
                    "structure comparison requires declared increasing time on one axis");
            require(t0.at("population").text() == t1.at("population").text(),
                    "structure comparison requires the same declared population");
            return;
        }
        if (operation == "native_boundary") {
            members(query, {"operation", "field_set", "backend"}, {"workers", "tile_rows"});
            const auto backend = query.at("backend").text();
            require(backend == "cpu_cached" || backend == "cpu_cached_tiles" ||
                        backend == "cpu_tiles" || backend == "cpu_scheduled",
                    "native field backend unavailable in this qualified profile");
            if (backend == "cpu_cached")
                require(!query.has("workers") && !query.has("tile_rows"),
                        "serial cache action has no tile/worker parameters");
            else {
                require(backend != "cpu_scheduled" || !query.has("workers"),
                        "scheduled workers are selected from canonical measurements");
                if (query.has("workers"))
                    require(count(query.at("workers"), cpu().coordinator().capacity(), "workers") >
                                0,
                            "workers must be positive");
                if (query.has("tile_rows"))
                    require(count(query.at("tile_rows"), coordinate_limit, "tile rows") > 0,
                            "tile rows must be positive");
            }
            (void)native_analysis(query.at("field_set").text());
            return;
        }
        const auto name = query.at("context").text();
        auto &prepared = context(name);
        if (operation == "typed_moment") {
            members(query, {"operation", "context", "field", "other_field", "order"});
            const auto a = input(name, query.at("field").text());
            const auto b = input(name, query.at("other_field").text());
            require(a.second == b.second, "moment fields have different grades");
            (void)count(query.at("order"), std::numeric_limits<std::size_t>::max(), "moment order");
        } else if (operation == "temporal_moment_change") {
            members(query,
                    {"operation", "context", "other_context", "field", "other_field", "new_field",
                     "new_other_field", "order"},
                    {"correspondence"});
            const auto other = query.at("other_context").text();
            temporal_times(name, other);
            const auto k = input(name, query.at("field").text(), false).second;
            require(input(name, query.at("other_field").text(), false).second == k &&
                        input(other, query.at("new_field").text(), false).second == k &&
                        input(other, query.at("new_other_field").text(), false).second == k,
                    "moment fields require one declared source grade");
            temporal_meanings(query.at("field").text(), query.at("new_field").text());
            temporal_meanings(query.at("other_field").text(), query.at("new_other_field").text());
            if (query.has("correspondence"))
                (void)temporal_map(name, other, k, query.at("correspondence"));
            else
                require(prepared.to_ambient.at(k) == context(other).to_ambient.at(k),
                        "different temporal carriers require an explicit correspondence");
            (void)count(query.at("order"), std::numeric_limits<std::size_t>::max(), "moment order");
        } else {
            require(operation == "boundary" || operation == "adjoint" || operation == "hodge" ||
                        operation == "green_action",
                    "unsupported scientific operation: " + operation);
            if (operation == "green_action") {
                members(query, {"operation", "context", "field", "lambda"});
                require(exact_rational(query.at("lambda").text()) > 0,
                        "Green action requires positive lambda");
            } else
                members(query, {"operation", "context", "field"});
            const auto grade = input(name, query.at("field").text()).second;
            if (operation == "boundary")
                require(grade > 0, "grade-zero field has no boundary");
            if (operation == "adjoint")
                require(grade + 1u < prepared.source->grade_count(),
                        "no declared upper grade for adjoint");
        }
    }
    Json provenance(StableCellKey cell) const {
        const auto endpoint = source.directory->pin_snapshot();
        auto result = Json::object();
        result.set("entity", source.cell_record(cell));
        result.set("evidence", source.evidence(source.cell_record(cell).at("evidence")));
        auto retained = Json::array();
        for (std::size_t k = 1; k < endpoint->source.grade_count(); ++k) {
            const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(k)};
            const auto &arena = endpoint->source.grade(grade);
            for (std::size_t i = 0; i < arena.cell_count(); ++i) {
                const auto upper =
                    stable_cell_key(endpoint->source, {grade, static_cast<std::uint32_t>(i)});
                const auto [begin, end] = arena.incidence_range(i);
                for (std::size_t slot = begin; slot < end; ++slot) {
                    const auto lower = stable_cell_key(
                        endpoint->source, {rfc::core::Grade{static_cast<std::uint16_t>(k - 1)},
                                           arena.lower_index[slot]});
                    if (cell != upper && cell != lower)
                        continue;
                    auto item = Json::object();
                    const auto occurrence =
                        StableOccurrenceKey{grade, arena.occurrence_identities[slot]};
                    item.set("key", Json::string(std::to_string(occurrence.key.value)));
                    item.set("grade", cardinality(k));
                    item.set("upper", source.cell_record(upper).at("id"));
                    item.set("lower", source.cell_record(lower).at("id"));
                    item.set("coefficient",
                             rational(endpoint->source.incidence_coefficient(grade, slot)));
                    item.set("magnitude",
                             rational(abs(endpoint->source.incidence_coefficient(grade, slot))));
                    auto aspects = Json::object();
                    for (std::size_t a = 0; a < endpoint->source.signature().aspects().size(); ++a)
                        aspects.set(endpoint->source.signature()
                                        .aspects()
                                        .at(static_cast<rfc::composite::AspectId>(a))
                                        .name,
                                    Json::boolean(arena.aspects.get(
                                        static_cast<rfc::composite::AspectId>(a), slot)));
                    item.set("aspects", aspects);
                    item.set("evidence",
                             source.evidence(source.occurrence_record(occurrence).at("evidence")));
                    retained.push(item);
                    require(retained.size() <= coordinate_limit, "provenance occurrence budget");
                }
            }
        }
        result.set("occurrences", retained);
        return result;
    }
    Json coordinates(const PreparedContext &prepared, std::uint16_t grade,
                     const Vector &values) const {
        check_values(values);
        auto result = Json::array();
        const auto endpoint = source.directory->pin_snapshot();
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto cell = stable_cell_key(
                endpoint->source, {rfc::core::Grade{grade}, prepared.to_ambient.at(grade).at(i)});
            auto item = Json::object();
            item.set("coordinate", cardinality(i));
            item.set("native_key", Json::string(std::to_string(cell.key)));
            item.set("entity", source.cell_record(cell).at("id"));
            item.set("value", rational(values[i]));
            result.push(item);
        }
        return result;
    }
    Json certificate(const Json &query) const {
        auto result = Json::object();
        result.set("schema", Json::string("gadget.scientific.certificate.v1"));
        result.set("source_digest", Json::string(source.digest));
        result.set("native_sha256", Json::string(source.native_digest));
        result.set("source_namespace", source.declaration.at("native").at("source_namespace"));
        result.set("source_sequence", Json::string(std::to_string(source.directory->sequence())));
        if (source.declaration.at("native").has("cutoff"))
            result.set("native_history", source.history_metadata());
        result.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        result.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
        result.set("query", query);
        result.set("design", source.declaration.at("design"));
        result.set(
            "arithmetic",
            Json::string("GMP exact rational; no reduction, imputation or numerical fallback"));
        auto meanings = Json::object();
        for (const auto &name : {"field", "other_field", "new_field", "new_other_field"})
            if (query.has(name)) {
                const auto id = query.at(name).text();
                meanings.set(id, source.fields.at(id).declaration);
            }
        if (query.has("field_set")) {
            const auto fields =
                source.field_sets.at(query.at("field_set").text()).declaration.at("fields");
            for (std::size_t i = 0; i < fields.size(); ++i) {
                const auto id = fields.at(i).text();
                meanings.set(id, source.fields.at(id).declaration);
            }
        }
        if (query.has("fields") && query.at("operation").text() != "hodge_frame_decomposition" &&
            query.at("operation").text() != "linear_response_observation")
            for (std::size_t i = 0; i < query.at("fields").size(); ++i) {
                const auto id = query.at("fields").at(i).text();
                meanings.set(id, source.fields.at(id).declaration);
            }
        if (query.at("operation").text() == "local_sections") {
            const auto sections = query.at("sections");
            for (std::size_t i = 0; i < sections.size(); ++i)
                if (sections.at(i).has("fields"))
                    for (std::size_t f = 0; f < sections.at(i).at("fields").size(); ++f) {
                        const auto id = sections.at(i).at("fields").at(f).text();
                        meanings.set(id, source.fields.at(id).declaration);
                    }
        }
        if (query.at("operation").text() == "structural_reference_calibration")
            for (const auto &id : research_detail::structural_calibration_field_ids(query))
                meanings.set(id, source.fields.at(id).declaration);
        if (query.at("operation").text() == "source_structure_ablation") {
            for (std::size_t r = 0; r < query.at("roles").size(); ++r) {
                const auto role = query.at("roles").at(r);
                if (role.has("fields"))
                    for (std::size_t f = 0; f < role.at("fields").size(); ++f) {
                        const auto id = role.at("fields").at(f).text();
                        meanings.set(id, source.fields.at(id).declaration);
                    }
            }
            if (query.has("sheaf"))
                for (std::size_t s = 0; s < query.at("sheaf").at("sections").size(); ++s) {
                    const auto section = query.at("sheaf").at("sections").at(s);
                    if (section.has("fields"))
                        for (std::size_t f = 0; f < section.at("fields").size(); ++f) {
                            const auto id = section.at("fields").at(f).text();
                            meanings.set(id, source.fields.at(id).declaration);
                        }
                }
        }
        if (query.at("operation").text() == "observed_response_model" ||
            query.at("operation").text() == "constrained_response_model")
            for (const auto &id : research_detail::response_model_field_ids(query))
                meanings.set(id, source.fields.at(id).declaration);
        if (query.at("operation").text() == "temporal_field_history") {
            result.set("arithmetic",
                       Json::string("GMP exact rational; original observation support and declared "
                                    "conditional completions retained; no numerical fallback"));
            auto bindings = Json::array();
            for (std::size_t t = 0; t < query.at("endpoints").size(); ++t) {
                const auto e = query.at("endpoints").at(t);
                const auto pub = e.has("publication") ? e.at("publication").text() : "self";
                const auto &owner_ = pub == "self" ? source : *source.publications.at(pub)->state_;
                auto binding = Json::object(), local_meanings = Json::object();
                binding.set("publication", Json::string(pub));
                binding.set("source_digest", Json::string(owner_.digest));
                binding.set("declaration_sha256", Json::string(owner_.declaration_digest));
                if (pub != "self")
                    binding.set("admission", source.publication_declarations.at(pub));
                const auto id = e.at("field").text();
                local_meanings.set(id, owner_.fields.at(id).declaration);
                for (std::size_t i = 0; i < query.at("word").size(); ++i)
                    if (query.at("word").at(i).has("fields")) {
                        const auto factor = query.at("word").at(i).at("fields").at(t).text();
                        local_meanings.set(factor, owner_.fields.at(factor).declaration);
                    }
                binding.set("field_semantics", local_meanings);
                bindings.push(binding);
                // Preserve legacy one-source keys; publication-local meanings
                // are retained without merging equal labels from other owners.
                if (pub == "self") {
                    meanings.set(id, owner_.fields.at(id).declaration);
                    for (std::size_t i = 0; i < query.at("word").size(); ++i)
                        if (query.at("word").at(i).has("fields")) {
                            auto factor = query.at("word").at(i).at("fields").at(t).text();
                            meanings.set(factor, owner_.fields.at(factor).declaration);
                        }
                }
            }
            result.set("publication_bindings", bindings);
        }
        if (query.at("operation").text() == "censored_endpoint_bounds") {
            for (const auto *key :
                 {"lower", "upper", "lower_closed", "upper_closed", "upper_finite"}) {
                auto id = query.at(key).text();
                meanings.set(id, source.fields.at(id).declaration);
            }
            for (const auto *key : {"selection", "strata"})
                for (std::size_t i = 0; i < query.at(key).size(); ++i) {
                    auto id = query.at(key).at(i).text();
                    meanings.set(id, source.fields.at(id).declaration);
                }
        }
        if (query.at("operation").text() == "linear_response_observation" ||
            query.at("operation").text() == "hodge_frame_decomposition") {
            const auto id = text(query.at("publication"));
            const auto &bound = id == "self" ? source : *source.publications.at(id)->state_;
            for (std::size_t i = 0; i < query.at("fields").size(); ++i) {
                const auto field = text(query.at("fields").at(i));
                meanings.set(field, bound.fields.at(field).declaration);
            }
        }
        result.set("field_semantics", meanings);
        if (is_research_operation(query.at("operation").text())) {
            for (const auto &id : research_field_ids(query))
                meanings.set(id, source.fields.at(id).declaration);
        }
        if (query.has("context"))
            result.set("context", source.contexts.at(query.at("context").text()).declaration);
        if (query.has("other_context"))
            result.set("other_context",
                       source.contexts.at(query.at("other_context").text()).declaration);
        return result;
    }
    Json coverage(const std::string &field) const {
        const auto publication = source.directory->pin_snapshot();
        const auto &t = tensor(field);
        auto coverage = Json::object();
        coverage.set("coordinates", cardinality(t.tensor.domain().logical_size()));
        coverage.set("reported", cardinality(t.tensor.support().count()));
        coverage.set("missing",
                     cardinality(t.tensor.domain().logical_size() - t.tensor.support().count()));
        coverage.set("native_tensor", Json::string(std::to_string(t.key)));
        coverage.set("native_carrier", Json::string(std::to_string(*t.output)));
        auto axes = Json::array();
        for (const auto extent : t.tensor.domain().extents())
            axes.push(cardinality(extent));
        coverage.set("extents", axes);
        const auto carriers = temporal_endpoint_detail::carriers(*publication);
        const auto &carrier = *carriers.at(*t.output);
        auto observations = Json::array();
        for (std::size_t i = 0; i < t.tensor.domain().logical_size(); ++i) {
            auto row = Json::object();
            row.set("coordinate", cardinality(i));
            if (carrier.kind == TemporalCarrierKind::SourceGrade) {
                const auto cell = source.coordinate(carrier.grade.value,i);
                row.set("entity", source.cell_record(cell).at("id"));
            }
            const bool present = t.tensor.support().contains(i);
            row.set("reported", Json::boolean(present));
            row.set("value", present ? rational(*t.tensor.state(i)) : Json{});
            observations.push(row);
        }
        coverage.set("observations", observations);
        return coverage;
    }
    Json chain_declaration(const Json &query) {
        members(query, {"operation", "context", "grade", "coefficient_domain"},
                {"field", "values"});
        require(query.has("field") != query.has("values"),
                "chain requires values or a source field");
        if (query.has("values"))
            return query;
        require(query.at("coefficient_domain").text() == "rational",
                "source fields have rational coordinates");
        const auto [v, g] = input(text(query.at("context")), text(query.at("field")));
        require(
            g == count(query.at("grade"), std::numeric_limits<std::uint16_t>::max(), "chain grade"),
            "chain field grade mismatch");
        auto resolved = Json::object(), values = Json::array();
        for (const auto *name : {"operation", "context", "grade", "coefficient_domain"})
            resolved.set(name, query.at(name));
        for (const auto &x : v)
            values.push(Json::string(x.get_str()));
        resolved.set("values", values);
        return resolved;
    }
    Json run(const Json &query, std::optional<rfc::runtime::ExecutionLease> native_lease = {}) {
        auto result = certificate(query);
        const auto operation = query.at("operation").text();
        std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> serial;
        if (resources && operation != "native_boundary" && operation != "support_mask" &&
            operation != "paired_field_contrasts" && operation != "controlled_pairings" &&
            operation != "linear_model_comparison" && operation != "finite_reference_calibration" &&
            operation != "relational_field_ensemble" &&
            operation != "structural_reference_calibration" &&
            operation != "source_structure_ablation" && operation != "observed_response_model" &&
            operation != "constrained_response_model" && operation != "hodge_frame_decomposition" &&
            operation != "linear_response_observation" && operation != "metric_relative_novelty")
            serial.emplace(resources->reserve(1));
        if (operation == "metric_relative_novelty") {
            result.set("result",
                       native_novelty.execute(
                           query, novelty_binding(query, true),
                           native_lease.value_or(scientific_detail::source_execution_lease(
                               "gadget-native-novelty-machine-addressable-v1")),
                           cpu()));
            return seal(result, "certificate_sha256");
        }
        if (operation == "declared_response_transition") {
            result.set("result",
                       declared_responses.execute(
                           query, response_binding(query),
                           native_lease.value_or(scientific_detail::declared_response_lease())));
            return seal(result, "certificate_sha256");
        }
        if (operation == "linear_response_observation") {
            result.set("result",
                       native_observations.execute(
                           query, observation_binding(query),
                           native_lease.value_or(scientific_detail::source_execution_lease(
                               "gadget-native-observation-machine-addressable-v1")),
                           cpu()));
            return seal(result, "certificate_sha256");
        }
        if (operation == "hodge_frame_decomposition") {
            result.set("result",
                       native_hodge.execute(
                           query, hodge_binding(query),
                           native_lease.value_or(scientific_detail::source_execution_lease(
                               "gadget-native-hodge-machine-addressable-v1")),
                           cpu()));
            return seal(result, "certificate_sha256");
        }
        if (operation == "temporal_field_history") {
            result.set("result", history_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (operation == "censored_endpoint_bounds") {
            result.set("result", censor_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (operation == "observed_response_model" || operation == "constrained_response_model") {
            result.set("result", response_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (operation == "source_structure_ablation") {
            result.set("result", ablation_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (operation == "structural_reference_calibration") {
            result.set("result", structural_calibration_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (operation == "relational_field_ensemble") {
            result.set("result", ensemble_query(query, false));
            return seal(result, "certificate_sha256");
        }
        if (is_research_operation(query.at("operation").text())) {
            result.set("result", execute_research_query(*owner, query, false, resources));
            return seal(result, "certificate_sha256");
        }
        if (operation == "chain_certificate" || operation == "realized_path_comparison") {
            auto &prepared = context(query.at("context").text());
            auto out = operation == "chain_certificate"
                           ? analyze_chain(*prepared.study, chain_declaration(query))
                           : compare_realized_paths(*prepared.study, query);
            const auto grade = static_cast<std::uint16_t>(count(
                query.at("grade"), std::numeric_limits<std::uint16_t>::max(), "structural grade"));
            auto coordinates = out.at("coordinates");
            for (std::size_t i = 0; i < coordinates.size(); ++i) {
                const auto stable =
                    source.coordinate(grade,prepared.to_ambient.at(grade).at(i));
                coordinates.at(i).set("entity", source.cell_record(stable).at("id"));
                coordinates.at(i).set("evidence", source.cell_record(stable).at("evidence"));
            }
            if (out.has("occurrences")) {
                auto occurrences = out.at("occurrences");
                for (std::size_t i = 0; i < occurrences.size(); ++i) {
                    const StableOccurrenceKey stable{
                        rfc::core::Grade{grade},
                        rfc::core::OccurrenceKey{key(occurrences.at(i).at("key"))}};
                    occurrences.at(i).set("evidence", source.occurrence_record(stable).at("evidence"));
                }
            }
            result.set("result", out);
        } else if (operation == "local_sections") {
            result.set("result", analyze_local_sections(*context(query.at("context").text()).study,
                                                        local_declaration(query)));
        } else if (operation == "tensor_window") {
            const auto id = query.at("field").text();
            const auto pin = source.field_pins.at(source.fields.at(id).tensor);
            const rfc::runtime::WorkAvailability availability{
                pin.source_namespace(), publication_cutoff(query.at("available_sequence"))};
            auto key_query = Json::parse(query.dump());
            key_query.set("preparation", Json::string("reuse"));
            const auto cache_key = key_query.dump();
            const bool retain = !query.has("preparation") ||
                                query.at("preparation").text() == "reuse";
            const bool reused = retain && windows.contains(cache_key);
            const auto lease = native_lease.value_or(source_execution_lease("scientific-tensor-window"));
            std::optional<rfc::database::NativeTensorWindow> transient;
            if (!reused) {
                std::vector<std::size_t> lo, hi;
                for (std::size_t a = 0; a < query.at("lower").size(); ++a) {
                    lo.push_back(static_cast<std::size_t>(query.at("lower").at(a).integer()));
                    hi.push_back(static_cast<std::size_t>(query.at("upper").at(a).integer()));
                }
                auto prepared = rfc::database::NativeTensorWindow::prepare(
                    pin, lo, hi, availability, lease);
                if (retain)
                    windows.emplace(cache_key, std::move(prepared));
                else
                    transient.emplace(std::move(prepared));
            }
            const auto &window = retain ? windows.at(cache_key) : *transient;
            const auto summary = window.summarize(availability, lease);
            auto out = Json::object(), observations = Json::array();
            require(window.domain().logical_size() <= coordinate_limit, "window output budget");
            for (std::size_t i = 0; i < window.domain().logical_size(); ++i) {
                lease.check_cancelled();
                auto row = Json::object();
                row.set("source_coordinate", cardinality(window.source_index(i)));
                const auto v = window.tensor_pin()->state(window.source_index(i));
                row.set("reported", Json::boolean(v.has_value()));
                row.set("value", v ? rational(*v) : Json{});
                observations.push(row);
            }
            out.set("observations", observations);
            out.set("reported", cardinality(window.supported_count()));
            out.set("missing",
                    cardinality(window.domain().logical_size() - window.supported_count()));
            out.set("sum", summary.sum() ? rational(*summary.sum()) : Json{});
            out.set("lower", query.at("lower"));
            out.set("upper", query.at("upper"));
            out.set("source_namespace", Json::string(std::to_string(pin.source_namespace())));
            out.set("source_sequence", Json::string(std::to_string(pin.sequence())));
            out.set("preparation_reused", Json::boolean(reused));
            out.set("preparation_retained", Json::boolean(retain));
            result.set("result", out);
        } else if (operation == "green_response") {
            const auto name = query.at("context").text();
            const auto [x, grade] = input(name, query.at("field").text());
            auto &prepared = context(name);
            const auto lambda = exact_decimal(query.at("lambda").text());
            auto operation = query;
            operation = Json::parse(operation.dump());
            operation.set("field", Json::string(""));
            const auto cache_key =
                operation.dump() + "|grade:" + std::to_string(grade) +
                "|carrier:" + std::to_string(*tensor(query.at("field").text()).output);
            const bool reused = green_solves.contains(cache_key);
            if (!reused) {
                require(green_solves.size() < 128, "prepared Green count budget");
                auto h = prepared.study->hodge(rfc::core::Grade{grade}).materialize();
                auto a = rfc::numeric::ExactMatrix::identity(x.size()).add(h.scale(lambda));
                std::vector<std::size_t> interior;
                for (std::size_t i = 0; i < query.at("interior").size(); ++i)
                    interior.push_back(
                        static_cast<std::size_t>(query.at("interior").at(i).integer()));
                const auto carrier = *tensor(query.at("field").text()).output;
                const auto coordinates = key(query.at("coordinates").at("id"));
                rfc::local::LinearSolveDeclaration declaration{
                    carrier,
                    carrier,
                    coordinates,
                    coordinates,
                    rfc::numeric::IndexDomain::linear(x.size()),
                    rfc::numeric::IndexDomain::linear(x.size())};
                // Context, grade and carrier are part of the prepared operator's meaning.
                green_solves.emplace(
                    cache_key, rfc::local::PreparedInteriorInterfaceSolve(
                                   std::move(declaration), std::move(a), std::move(interior)));
            }
            const auto &solve = green_solves.at(cache_key);
            const auto y = solve.apply(x);
            check_values(y);
            auto out = Json::object();
            out.set("field", coordinates(prepared, grade, y));
            out.set("interior", query.at("interior"));
            out.set("interface_coordinates", Json::array());
            auto interface = Json::array();
            for (const auto i : solve.interface())
                interface.push(cardinality(i));
            out.set("interface_coordinates", interface);
            out.set("factorization_reused", Json::boolean(reused));
            out.set("residual_zero", Json::boolean(true));
            result.set("result", out);
        } else if (operation == "inspect")
            result.set("result", source.catalog());
        else if (operation == "provenance")
            result.set("result", provenance(source.resolve_cell(query.at("entity").text())));
        else if (operation == "coverage") {
            result.set("result", coverage(query.at("field").text()));
        } else if (operation == "field_summary") {
            const auto &t = tensor(query.at("field").text());
            Vector values;
            auto coordinates = Json::array();
            t.tensor.support().for_each_active([&](std::size_t i) {
                values.push_back(*t.tensor.state(i));
                coordinates.push(cardinality(i));
            });
            check_values(values);
            const Vector ones(values.size(), Rational{1});
            const Rational sum = rfc::numeric::dot(ones, values);
            const Rational square_sum = rfc::numeric::dot(values, values);
            auto summary = Json::object();
            summary.set("reported", cardinality(values.size()));
            summary.set("missing", cardinality(t.tensor.domain().logical_size() - values.size()));
            summary.set("coordinates", coordinates);
            summary.set("sum", rational(sum));
            summary.set("sum_of_squares", rational(square_sum));
            summary.set("mean", values.empty() ? Json{}
                                               : rational(Rational(sum / Rational(values.size()))));
            summary.set("units",
                        source.fields.at(query.at("field").text()).declaration.at("units"));
            result.set(
                "mathematical_definition",
                Json::string("exact canonical dot(1,x), dot(x,x), and dot(1,x)/reported_count on "
                             "explicit native reported support; descriptive observed rows, no "
                             "population estimator or imputation"));
            result.set("result", summary);
        } else if (operation == "support_mask") {
            using namespace rfc::composite;
            const auto fields = query.at("fields");
            const auto &first = tensor(fields.at(0).text());
            const auto identity = fields.dump() + ":" + query.at("mode").text();
            const bool reused = packed.contains(identity);
            if (!reused) {
                AspectSchema schema;
                for (std::size_t i = 0; i < fields.size(); ++i)
                    schema.add(fields.at(i).text());
                AspectPlanes planes(schema, first.tensor.domain().logical_size());
                for (std::size_t i = 0; i < fields.size(); ++i)
                    tensor(fields.at(i).text())
                        .tensor.support()
                        .for_each_active([&](std::size_t coordinate) {
                            planes.set(static_cast<AspectId>(i), coordinate);
                        });
                BinaryProgram program(2);
                program.append({BinaryOpcode::LoadAspect, 0, 0, 0, 0, 0});
                for (std::size_t i = 1; i < fields.size(); ++i) {
                    program.append(
                        {BinaryOpcode::LoadAspect, 1, 0, 0, 0, static_cast<AspectId>(i)});
                    program.append({query.at("mode").text() == "all_reported" ? BinaryOpcode::And
                                                                              : BinaryOpcode::Or,
                                    0, 0, 1});
                }
                if (query.at("mode").text() == "none_reported")
                    program.append({BinaryOpcode::Not, 0, 0});
                packed.emplace(
                    identity,
                    std::make_unique<PackedOperation>(
                        rfc::backend::PreparedPackedBinary(program, planes, {0}), compute));
            }
            auto derived = packed.at(identity)->execute(cpu(), query);
            derived.execution.set("prepared_operation_reused", Json::boolean(reused));
            auto selection = Json::object();
            selection.set("coordinates", cardinality(first.tensor.domain().logical_size()));
            selection.set("selected", cardinality(derived.outputs.at(0).count()));
            selection.set("native_carrier", Json::string(std::to_string(*first.output)));
            auto rows = Json::array();
            const auto carrier =
                temporal_endpoint_detail::carriers(*source.directory->pin_snapshot())
                    .at(*first.output);
            derived.outputs.at(0).for_each_set([&](std::size_t i) {
                auto row = Json::object();
                row.set("coordinate", cardinality(i));
                const auto cell = source.coordinate(carrier->grade.value,i);
                row.set("entity", source.cell_record(cell).at("id"));
                row.set("native_key", Json::string(std::to_string(cell.key)));
                rows.push(row);
            });
            selection.set("selection", rows);
            result.set("arithmetic", Json::string("exact packed Boolean program over native "
                                                  "reported support; supported zero is reported"));
            result.set(
                "mathematical_definition",
                Json::string(
                    "declared ordered AND, OR, or NOT(OR) of reported-support planes on one "
                    "identical native carrier; derived eligibility only, source unmodified"));
            result.set("execution", derived.execution);
            result.set("result", selection);
        } else if (operation == "native_boundary") {
            const auto name = query.at("field_set").text();
            const auto analysis = native_analysis(name);
            Vector values;
            Json execution = Json::object();
            rfc::database::NativeExecutionReport report;
            if (query.at("backend").text() == "cpu_cached" ||
                query.at("backend").text() == "cpu_cached_tiles") {
                const auto derived = [&]() -> rfc::database::CachedNativeField {
                    if (query.at("backend").text() == "cpu_cached") {
                        if (resources)
                            serial.emplace(resources->reserve(1));
                        return field_cache.evaluate(analysis);
                    }
                    const auto rows = query.has("tile_rows") ? count(query.at("tile_rows"),
                                                                     coordinate_limit, "tile rows")
                                                             : 64;
                    const auto workers =
                        query.has("workers")
                            ? count(query.at("workers"), cpu().coordinator().capacity(), "workers")
                            : 1;
                    auto &host = cpu();
                    (void)host.readmit();
                    auto &coordinator = host.coordinator();
                    auto reservation = coordinator.reserve(std::min(
                        workers, std::max(std::size_t{1},
                                          rfc::machine::ceil_div(
                                              analysis.binding().target_dimension(), rows))));
                    rfc::runtime::CpuTileExecutor executor(
                        {reservation.workers(), coordinator.capacity(), coordinate_limit});
                    execution.set("requested_workers", cardinality(workers));
                    execution.set("admitted_workers", cardinality(reservation.workers()));
                    return field_cache.evaluate_parallel(analysis, executor, rows);
                }();
                values = derived.values();
                report = derived.provenance();
                execution.set("plan_reused", Json::boolean(derived.reuse().plan));
                execution.set("result_reused", Json::boolean(derived.reuse().result));
                execution.set("in_flight_reused", Json::boolean(derived.reuse().in_flight));
                auto producer = Json::object();
                if (derived.producer_execution()) {
                    const auto &facts = *derived.producer_execution();
                    producer.set("worker_slots", cardinality(facts.worker_slots));
                    producer.set("launched_threads", cardinality(facts.launched_threads));
                    producer.set("tiles", cardinality(facts.tiles));
                    producer.set("completed_tiles", cardinality(facts.completed_tiles));
                }
                execution.set("producer_execution", producer);
                execution.set(
                    "producer_scope",
                    Json::string(
                        "retained producer facts; a result hit performs no new field action"));
            } else {
                const auto rows = query.has("tile_rows")
                                      ? count(query.at("tile_rows"), coordinate_limit, "tile rows")
                                      : 64;
                const auto workers =
                    query.has("workers")
                        ? count(query.at("workers"), cpu().coordinator().capacity(), "workers")
                        : 1;
                const bool automatic = query.at("backend").text() == "cpu_scheduled";
                auto &host = cpu();
                const auto current = host.readmit();
                auto &coordinator = host.coordinator();
                const auto derived = [&]() -> rfc::database::DerivedParallelNativeField {
                    if (!automatic) {
                        const auto prepared = rfc::database::prepare_parallel_field(analysis, rows);
                        auto workspace = rfc::database::make_parallel_workspace(prepared);
                        auto reservation = coordinator.reserve(std::min(
                            workers, std::max<std::size_t>(1, prepared.ranges().tiles().size())));
                        rfc::runtime::CpuTileExecutor executor(
                            {reservation.workers(), coordinator.capacity(), coordinate_limit});
                        return rfc::database::execute_parallel(prepared, workspace, executor);
                    }
                    const auto identity = std::pair{name, rows};
                    const bool reused = scheduled.contains(identity) &&
                                        scheduled.at(identity)->calibration.has_value();
                    if (!scheduled.contains(identity))
                        scheduled.emplace(
                            identity, std::make_unique<Scheduled>(
                                          rfc::database::prepare_parallel_field(analysis, rows)));
                    auto &entry = *scheduled.at(identity);
                    if (!entry.calibration) {
                        rfc::database::CpuFieldCalibrationPolicy policy;
                        policy.worker_candidates.clear();
                        const auto ceiling = std::max<std::size_t>(
                            1, std::min(coordinator.capacity(),
                                        entry.operation.ranges().tiles().size()));
                        for (std::size_t candidate = 1; candidate < ceiling;) {
                            policy.worker_candidates.push_back(candidate);
                            if (candidate > ceiling / 2)
                                break;
                            candidate *= 2;
                        }
                        policy.worker_candidates.push_back(ceiling);
                        policy.call_limit = rfc::machine::checked_size_add(
                            1, rfc::machine::checked_size_mul(
                                   policy.worker_candidates.size(),
                                   rfc::machine::checked_size_add(policy.samples, policy.warmups)));
                        entry.calibration = rfc::database::calibrate_cpu_field(
                            entry.operation, entry.workspace, coordinator, current, policy);
                    }
                    const auto chosen = rfc::database::execute_scheduled_field(
                        entry.operation, entry.workspace, coordinator, current, *entry.calibration);
                    execution.set("calibration_reused", Json::boolean(reused));
                    execution.set("selected_workers", cardinality(chosen.selected_workers));
                    auto samples = Json::array();
                    for (const auto &timing : entry.calibration->timings()) {
                        auto sample = Json::object();
                        sample.set("workers", cardinality(timing.workers));
                        sample.set("median_nanoseconds",
                                   Json::string(std::to_string(timing.median_nanoseconds)));
                        auto times = Json::array();
                        for (const auto ns : timing.nanosecond_samples)
                            times.push(Json::string(std::to_string(ns)));
                        sample.set("nanosecond_samples", times);
                        samples.push(sample);
                    }
                    execution.set("calibration", samples);
                    execution.set(
                        "selection",
                        Json::string("canonical lowest measured median among currently fitting "
                                     "candidates; owning prepared operation and CPU facts"));
                    return chosen.field;
                }();
                execution.set("cpu_resources", host.report());
                values = derived.values();
                report = derived.provenance();
                execution.set("requested_workers",
                              cardinality(derived.execution().requested_workers));
                execution.set("worker_slots", cardinality(derived.execution().worker_slots));
                execution.set("launched_threads",
                              cardinality(derived.execution().launched_threads));
                execution.set("tiles", cardinality(derived.execution().tiles));
                execution.set("completed_tiles", cardinality(derived.execution().completed_tiles));
            }
            check_values(values);
            auto output = Json::array();
            for (const auto &binding : analysis.binding().outputs()) {
                auto row = Json::object();
                row.set("entity", source.cell_record(binding.cell).at("id"));
                row.set("carrier", Json::string(std::to_string(binding.carrier)));
                row.set("coordinate_system", Json::string(std::to_string(binding.coordinates)));
                row.set("values",
                        vector_json(Vector(values.data() + binding.offset,
                                           values.data() + binding.offset + binding.dimension)));
                row.set("support",
                        Json::string("complete; supported zero and cancellation retained"));
                output.push(row);
            }
            execution.set("profile", Json::string(report.profile));
            execution.set("implementation", Json::string(report.implementation));
            execution.set("backend", Json::string("exact_cpu"));
            execution.set("fallback", Json::boolean(report.fallback));
            execution.set("retained_occurrences",
                          cardinality(analysis.layout().retained_occurrences));
            execution.set("active_occurrences", cardinality(analysis.layout().active_occurrences));
            execution.set("source_coordinates", cardinality(analysis.layout().source_coordinates));
            execution.set("target_coordinates", cardinality(analysis.layout().target_coordinates));
            auto decisions = Json::array();
            for (const auto &decision : report.decisions) {
                auto row = Json::object();
                row.set("subject", Json::string(decision.choice));
                row.set("origin", Json::string(decision.origin));
                row.set("value", Json::string(decision.value));
                row.set("constraints", Json::string(decision.constraints));
                decisions.push(row);
            }
            execution.set("decisions", decisions);
            result.set("execution", execution);
            result.set("field_set", source.field_sets.at(name).declaration);
            auto retained = Json::array();
            for (const auto &binding : analysis.binding().occurrences()) {
                auto occurrence = Json::object();
                occurrence.set("grade", cardinality(binding.occurrence.upper_grade.value));
                occurrence.set("key", Json::string(std::to_string(binding.occurrence.key.value)));
                occurrence.set("upper", source.cell_record(binding.upper).at("id"));
                occurrence.set("lower", source.cell_record(binding.lower).at("id"));
                occurrence.set("active", Json::boolean(binding.active));
                occurrence.set("coefficient", rational(binding.coefficient));
                occurrence.set("map",
                               binding.map ? Json::string(std::to_string(*binding.map)) : Json{});
                occurrence.set(
                    "map_role",
                    Json::string(!binding.map
                                     ? "declared_constant_identity"
                                     : (binding.role == TemporalTensorRole::RealizedIncidenceMap
                                            ? "realized_incidence"
                                            : "coefficient_independent")));
                occurrence.set(
                    "evidence",
                    source.evidence(source.occurrence_record(binding.occurrence).at("evidence")));
                retained.push(occurrence);
            }
            result.set("retained_occurrences", retained);
            result.set(
                "mathematical_definition",
                Json::string("native rational B_k field action; complete declared fibers and "
                             "occurrence maps, with retained inactive and repeated occurrences"));
            result.set("result", output);
        } else {
            const auto name = query.at("context").text();
            auto &prepared = context(name);
            const auto [values, grade] = input(name, query.at("field").text());
            if (operation == "structure_transition") {
                const auto other_name = query.at("other_context").text();
                auto &target = context(other_name);
                const auto next = input(other_name, query.at("new_field").text()).first;
                result.set("result",
                           research_detail::compare_prepared_structures(
                               temporal(prepared), temporal(target), rfc::core::Grade{grade},
                               prepared.to_ambient, target.to_ambient, values, next,
                               exact_decimal(query.at("lambda").text()), nullptr, false,
                               research_detail::FactoredTemporalMaps::requested(query)));
            } else if (operation == "typed_moment") {
                const auto other = input(name, query.at("other_field").text()).first;
                result.set("mathematical_definition",
                           Json::string("<x,L_k^n y>_M; channels remain independent"));
                result.set("result",
                           rational(prepared.study->moment(
                               rfc::core::Grade{grade},
                               count(query.at("order"), std::numeric_limits<std::size_t>::max(),
                                     "moment order"),
                               values, other)));
            } else if (operation == "temporal_moment_change") {
                const auto other_name = query.at("other_context").text();
                auto &target = context(other_name);
                const auto old_y = input(name, query.at("other_field").text()).first;
                const auto new_x = input(other_name, query.at("new_field").text()).first;
                const auto new_y = input(other_name, query.at("new_other_field").text()).first;
                const auto order = count(query.at("order"), std::numeric_limits<std::size_t>::max(),
                                         "moment order");
                const auto form = [order, grade](const NativeStudy &study) {
                    const auto hodge = study.hodge(rfc::core::Grade{grade});
                    const auto power =
                        rfc::field::materialize_action_matrix(hodge.metric(), [&](const Vector &v) {
                            return rfc::statistics::hodge_power_apply(hodge, order, v);
                        });
                    return study.tower().metrics()[grade].materialize().multiply(power);
                };
                const auto old_form = form(*prepared.study);
                const auto new_form = form(*target.study);
                const auto identity =
                    query.has("correspondence")
                        ? temporal_map(name, other_name, grade, query.at("correspondence"))
                        : rfc::numeric::ExactMatrix::identity(values.size());
                const auto change = rfc::temporal::changing_dimension_cross_moment(
                    values, old_y, old_form, new_x, new_y, new_form, identity, identity);
                const Rational direct =
                    target.study->moment(rfc::core::Grade{grade}, order, new_x, new_y) -
                    prepared.study->moment(rfc::core::Grade{grade}, order, values, old_y);
                require(change.total == direct, "native finite cross-moment identity failed");
                auto terms = Json::object();
                terms.set("pullback_form", rational(change.pullback_form));
                terms.set("left_innovation", rational(change.left_innovation));
                terms.set("right_innovation", rational(change.right_innovation));
                terms.set("innovation_interaction", rational(change.innovation_interaction));
                terms.set("total", rational(change.total));
                terms.set("identity_verified", Json::boolean(true));
                const auto t0 = source.contexts.at(name).declaration.at("time");
                const auto t1 = source.contexts.at(other_name).declaration.at("time");
                const Rational interval = exact_decimal(t1.at("position").text()) -
                                          exact_decimal(t0.at("position").text());
                terms.set("interval", rational(interval));
                terms.set("time_unit", t0.at("unit"));
                terms.set("finite_rate",
                          rational(rfc::temporal::finite_rate(change.total, interval)));
                result.set(
                    "mathematical_definition",
                    Json::string(
                        "complete finite cross-moment change with form, left/right innovation and "
                        "interaction; declared or established native correspondence"));
                result.set("result", terms);
            } else {
                rfc::field::OperatorRequest request;
                request.grade = rfc::core::Grade{grade};
                request.kind = rfc::field::OperatorKind::Hodge;
                auto output_grade = grade;
                if (operation == "boundary") {
                    request.kind = rfc::field::OperatorKind::Boundary;
                    --output_grade;
                }
                if (operation == "adjoint") {
                    request.kind = rfc::field::OperatorKind::Adjoint;
                    ++request.grade.value;
                    ++output_grade;
                }
                Vector output;
                if (operation == "green_action") {
                    const auto lambda = exact_rational(query.at("lambda").text());
                    const auto answer =
                        prepared.study->green_action(rfc::core::Grade{grade}, lambda, values);
                    require(answer.exact &&
                                rfc::numeric::is_zero(rfc::numeric::sub(
                                    rfc::numeric::add(
                                        answer.value,
                                        rfc::numeric::scale(
                                            lambda, prepared.study->hodge(rfc::core::Grade{grade})
                                                        .apply(answer.value))),
                                    values)),
                            "native Green action reconstruction residual");
                    output = answer.value;
                    result.set("green_iterations", cardinality(answer.iterations));
                    result.set(
                        "green_execution",
                        Json::string(
                            "rfc::field::exact_green_apply; factored boundary and metric actions"));
                    result.set("residual_verified", Json::boolean(true));
                } else
                    output = prepared.study->execute(request, values);
                result.set(
                    "mathematical_definition",
                    Json::string(
                        operation == "green_action" ? "(I + lambda L_k) x = j; exact operator "
                                                      "action with complete output field"
                        : operation == "hodge"
                            ? "L_k x = (B_k^dagger B_k + B_(k+1) B_(k+1)^dagger) x"
                            : (operation == "boundary" ? "B_k x" : "M_(k+1)^-1 B_(k+1)^T M_k x")));
                result.set("output_grade", cardinality(output_grade));
                result.set("result", coordinates(prepared, output_grade, output));
            }
        }
        return seal(result, "certificate_sha256");
    }
};
ScientificSession::ScientificSession(std::shared_ptr<const ScientificStudy> study,
                                     std::shared_ptr<CpuResources> resources,
                                     ComputeOptions compute) {
    require(static_cast<bool>(study), "session requires an immutable study owner");
    state_ = std::make_unique<State>(std::move(study), std::move(resources), compute);
}
ScientificSession::~ScientificSession() = default;
void ScientificSession::validate(const Json &query) {
    auto admitted = Json::parse(query.dump());
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
}
NativeRecord ScientificSession::execute_record(const NativeRecord &query) {
    const auto result=execute(render_record_view(query));
    return retain_operation(state_->owner->native_publication(),query,import_record_view(result));
}
Json ScientificSession::execute(const Json &query) {
    auto admitted = Json::parse(query.dump());
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted);
}
Json ScientificSession::execute_tensor_window(const Json &query,
                                               rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "tensor_window", "tensor window lease operation");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::tensor_window_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    std::size_t count = 0;
    const auto visit = [&](auto &&self, const State &state) -> void {
        count = rfc::machine::checked_size_add(count, state.windows.size());
        for (const auto &[id, child] : state.publication_sessions) {
            (void)id;
            self(self, *child);
        }
    };
    visit(visit, *state_);
    auto out = Json::object();
    out.set("retained_preparations", cardinality(count));
    out.set("scope", Json::string("this session and its admitted publication sessions"));
    return out;
}
void ScientificSession::clear_tensor_window_preparations() {
    std::lock_guard lock(state_->mutex);
    const auto visit = [&](auto &&self, State &state) -> void {
        state.windows.clear();
        for (auto &[id, child] : state.publication_sessions) {
            (void)id;
            self(self, *child);
        }
    };
    visit(visit, *state_);
}
Json ScientificSession::execute_native_response(const Json &query,
                                                rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "declared_response_transition",
            "native response lease requires declared_response_transition");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::native_response_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    return state_->declared_responses.statistics();
}
Json ScientificSession::execute_native_hodge(const Json &query,
                                             rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "hodge_frame_decomposition",
            "native Hodge lease requires hodge_frame_decomposition");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::native_hodge_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    return state_->native_hodge.statistics();
}
Json ScientificSession::execute_native_observation(const Json &query,
                                                   rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "linear_response_observation",
            "native observation lease requires linear_response_observation");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::native_observation_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    return state_->native_observations.statistics();
}
Json ScientificSession::execute_native_novelty(const Json &query,
                                               rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "metric_relative_novelty",
            "native novelty lease requires metric_relative_novelty");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::native_novelty_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    return state_->native_novelty.statistics();
}
void ScientificSession::clear_native_novelty_preparations() {
    std::lock_guard lock(state_->mutex);
    state_->native_novelty.clear();
}
NativeRecord ScientificSession::execute_plan_record(const NativeRecord &plan) {
    const auto result=execute_plan(render_record_view(plan));
    return retain_operation(state_->owner->native_publication(),plan,import_record_view(result));
}
Json ScientificSession::execute_plan(const Json &plan) {
    auto admitted = Json::parse(plan.dump());
    members(admitted, {"schema", "id", "question", "source_digest", "steps"});
    require(admitted.at("schema").text() == "gadget.scientific.plan.v1", "unsupported plan schema");
    (void)text(admitted.at("id"));
    (void)text(admitted.at("question"));
    require(admitted.at("source_digest").text() == state_->source.digest,
            "plan source binding mismatch");
    auto steps = admitted.at("steps");
    source_array(steps, "plan steps");
    require(steps.size() > 0, "plan requires at least one step");
    std::lock_guard lock(state_->mutex);
    std::set<std::string> identities;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        auto step = steps.at(i);
        members(step, {"id", "query"});
        require(identities.insert(text(step.at("id"))).second, "duplicate plan step identity");
        state_->admit(step.at("query"));
    }
    auto certificates = Json::array();
    for (std::size_t i = 0; i < steps.size(); ++i) {
        auto step = Json::object();
        step.set("id", steps.at(i).at("id"));
        step.set("certificate", state_->run(steps.at(i).at("query")));
        certificates.push(step);
    }
    auto result = Json::object();
    result.set("schema", Json::string("gadget.scientific.plan.certificate.v1"));
    result.set("status", Json::string("complete"));
    result.set("plan", admitted);
    result.set("plan_sha256", Json::string("sha256:" + sha256(admitted.dump())));
    result.set("source_digest", Json::string(state_->source.digest));
    result.set("steps", certificates);
    return seal(result, "certificate_sha256");
}
Json ScientificSession::propose_annotations(const Json &annotations) {
    auto admitted = Json::parse(annotations.dump());
    members(admitted, {"schema", "source_digest", "mentions"});
    require(admitted.at("schema").text() == "gadget.semantic.proposal.v1",
            "unsupported semantic proposal schema");
    require(admitted.at("source_digest").text() == state_->source.digest,
            "semantic source binding mismatch");
    auto mentions = admitted.at("mentions");
    source_array(mentions, "semantic mentions");
    require(mentions.size() > 0, "semantic proposal requires mentions");
    auto candidates = Json::array();
    for (std::size_t i = 0; i < mentions.size(); ++i) {
        auto item = Json::parse(mentions.at(i).dump());
        members(item, {"entity", "kind", "polarity", "evidence"});
        const auto name = item.at("entity").text();
        require(state_->source.has_cell(name), "semantic proposal entity is undeclared");
        const auto &cell = state_->source.cell_record(state_->source.resolve_cell(name));
        require(item.at("kind").text() == cell.at("kind").text(), "semantic entity kind mismatch");
        const auto polarity = item.at("polarity").text();
        require(polarity == "asserted" || polarity == "negated" || polarity == "uncertain",
                "semantic polarity must be explicit");
        auto evidence = state_->source.evidence(item.at("evidence"));
        for (std::size_t j = 0; j < evidence.size(); ++j)
            require(state_->source.documents.at(evidence.at(j).at("document").text())
                            .metadata.at("format")
                            .text() == "utf8",
                    "NER requires a declared UTF-8 document");
        item.set("evidence", evidence);
        auto native_cell = Json::object();
        native_cell.set("grade", cell.at("grade"));
        native_cell.set("key", cell.at("key"));
        native_cell.set("source_namespace",
                        state_->source.declaration.at("native").at("source_namespace"));
        item.set("native_cell", native_cell);
        item.set("status",
                 Json::string(
                     "prospective; evidence identity verified; semantic correctness unverified"));
        candidates.push(item);
    }
    auto result = Json::object();
    result.set("schema", Json::string("gadget.semantic.receipt.v1"));
    result.set("source_digest", Json::string(state_->source.digest));
    result.set("proposal", admitted);
    result.set("mentions", candidates);
    return seal(result, "receipt_sha256");
}
Json ScientificSession::native_cache_statistics() const {
    const auto stats = state_->field_cache.stats();
    auto result = Json::object();
    result.set("plan_hits", cardinality(stats.plan_hits));
    result.set("result_hits", cardinality(stats.result_hits));
    result.set("compiled_plans", cardinality(stats.compiled_plans));
    result.set("applications", cardinality(stats.applications));
    result.set("workspace_creations", cardinality(stats.workspace_creations));
    result.set("workspace_hits", cardinality(stats.workspace_hits));
    result.set("idle_workspaces", cardinality(stats.idle_workspaces));
    result.set("idle_workspace_slots", cardinality(stats.idle_workspace_slots));
    return result;
}
std::size_t ScientificSession::preparation_count() const {
    std::lock_guard lock(state_->mutex);
    auto total = state_->contexts.size();
    for (const auto &[id, child] : state_->publication_sessions) {
        (void)id;
        total = rfc::machine::checked_size_add(total, child->contexts.size());
    }
    return total;
}
Json ScientificSession::temporal_preparation_statistics() const {
    std::lock_guard lock(state_->mutex);
    std::array<std::size_t, 6> total{};
    const auto collect = [&](const State &state, const auto &self) -> void {
        for (const auto &[name, context] : state.contexts) {
            (void)name;
            if (!context.temporal)
                continue;
            const auto counts = context.temporal->statistics();
            for (std::size_t i = 0; i < total.size(); ++i)
                total[i] = rfc::machine::checked_size_add(total[i], counts[i]);
        }
        for (const auto &[name, child] : state.publication_sessions) {
            (void)name;
            self(*child, self);
        }
    };
    collect(*state_, collect);
    auto result = Json::object();
    const std::array<const char *, 6> names{"sector_preparations", "hodge_matrices",
                                            "green_preparations",  "woodbury_preparations",
                                            "green_matrices",      "green_factor_entries"};
    for (std::size_t i = 0; i < total.size(); ++i)
        result.set(names[i], cardinality(total[i]));
    return result;
}
Json ScientificSession::read_source(const Json &request) const {
    std::lock_guard lock(state_->mutex);
    auto worker = state_->cpu().reserve(1);
    return state_->owner->read_source(request);
}
Json apply_scientific_proposal(ScientificSession &session, const Json &proposal) {
    auto result = Json::object();
    result.set("schema", Json::string("gadget.agent.response.v1"));
    result.set("proposal", Json::parse(proposal.dump()));
    try {
        require(proposal.at("schema").text() == "gadget.scientific.agent.proposal.v1",
                "unsupported agent proposal schema");
        const auto action = proposal.at("action").text();
        if (action == "clarify") {
            members(proposal, {"schema", "action", "reason"});
            result.set("reason", Json::string(text(proposal.at("reason"))));
            result.set("status", Json::string("clarification"));
        } else if (action == "retrieve") {
            members(proposal, {"schema", "action", "read"});
            result.set("certificate", session.read_source(proposal.at("read")));
            result.set("status", Json::string("retrieved"));
        } else if (action == "query" || action == "plan" || action == "annotate") {
            const std::string member = action == "annotate" ? "annotations" : action;
            members(proposal, {"schema", "action", member});
            result.set("certificate",
                       action == "query"
                           ? session.execute(proposal.at(member))
                           : (action == "plan" ? session.execute_plan(proposal.at(member))
                                               : session.propose_annotations(proposal.at(member))));
            result.set("status", Json::string(action == "annotate" ? "prospective" : "executed"));
        } else
            throw std::invalid_argument("unsupported scientific model action");
    } catch (const rfc::runtime::CpuResourceBusy &error) {
        result.set("status", Json::string("busy"));
        result.set("error", Json::string(error.what()));
    } catch (const rfc::runtime::ExecutionLimitExceeded &error) {
        result.set("status", Json::string("rejected"));
        result.set("error", Json::string(error.what()));
    } catch (const std::invalid_argument &error) {
        result.set("status", Json::string("rejected"));
        result.set("error", Json::string(error.what()));
    }
    return result;
}
} // namespace gadget
