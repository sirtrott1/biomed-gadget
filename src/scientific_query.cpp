#include "research_ablation.hpp"
#include "scientific_ablation_fields.hpp"
#include "scientific_source_snapshot.hpp"
#include "scientific_structure.hpp"
#include "research_ensemble.hpp"
#include "research_factored.hpp"
#include "research_prepared.hpp"
#include "research_linear.hpp"
#include "research_comparisons.hpp"
#include "research_calibration.hpp"
#include "scientific_result_fields.hpp"
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
#include "scientific_sections.hpp"
#include <gadget/structure.hpp>
#include <gadget/study.hpp>
#include <mutex>
#include <rfc/algebra/materialize.hpp>
#include <rfc/database/native_field_cache.hpp>
#include <rfc/database/native_tensor_window.hpp>
#include <rfc/database/parallel_native_field.hpp>
#include <rfc/database/resource_native_field.hpp>
#include <rfc/local/interior_interface_solve.hpp>
#include <rfc/runtime/staged_resources.hpp>
#include <rfc/temporal/moment_change.hpp>
#include <rfc/temporal/factor_attribution.hpp>
#include <rfc/temporal/operation_word.hpp>

namespace gadget {
bool has_native_scientific_publisher(std::string_view operation) noexcept {
    return operation == "inspect" || operation == "coverage" || operation == "provenance" ||
           operation == "trajectory_structure" || operation == "temporal_field_history" || operation == "structure_transition" || operation == "paired_field_contrasts" || operation == "controlled_pairings" ||
           operation == "linear_model_comparison" || operation == "finite_reference_calibration" ||
           operation == "hodge_frame_decomposition" || operation == "linear_response_observation" ||
           operation == "declared_response_transition" || operation == "metric_relative_novelty" ||
           operation == "observed_response_model" || operation == "constrained_response_model" ||
           operation == "relational_field_ensemble" || operation == "structural_reference_calibration" ||
           operation == "local_sections" || operation == "source_structure_ablation" ||
           operation == "field_summary" || operation == "boundary" || operation == "adjoint" ||
           operation == "hodge" || operation == "green_action" || operation == "censored_endpoint_bounds" || operation == "typed_moment" || operation == "temporal_moment_change" || operation == "tensor_window" || operation == "support_mask" || operation == "native_boundary" || operation == "green_response" || operation == "chain_certificate" || operation == "realized_path_comparison";
}
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
struct PreparedContext {
    std::shared_ptr<const rfc::core::RelationalComplex> source;
    std::vector<std::vector<std::uint32_t>> to_ambient;
    std::shared_ptr<const NativeStudy> study;
    std::unique_ptr<PreparedTemporalStudy> temporal;
};
} // namespace
struct ScientificSession::State {
    std::shared_ptr<const ScientificStudy> owner;
    std::optional<RawByteSource> original_publication;
    const ScientificStudy::State &source;
    std::map<std::string, PreparedContext> contexts;
    std::map<std::string, std::unique_ptr<State>> publication_sessions;
    std::map<std::string, std::shared_ptr<const rfc::local::PreparedInteriorInterfaceSolve>> green_solves;
    struct GreenCoordinates {
        std::string binding;
        rfc::runtime::CpuResourceCoordinator::Reservation charge;
    };
    std::map<std::uint64_t, GreenCoordinates> green_spaces;
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
        // Legacy contexts share ambient witnesses. Explicit block contexts keep
        // precisely the carrier certified by their strict subcomplex.
        if (!declared.declaration.has("witness_scope"))
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
        return *source.field_pins.at(source.fields.at(field).tensor).tensor_declaration();
    }
    std::pair<Vector, std::uint16_t> input(const std::string &name, const std::string &field,
                                           bool legacy_budget = true) {
        const auto &t = tensor(field);
        const auto &meaning = source.fields.at(field).declaration;
        require(meaning.at("value_type").text() != "categorical",
                "categorical codes cannot be linear field amplitudes");
        const auto carrier = source.field_pins.at(t.key).carrier();
        require(carrier->kind == TemporalCarrierKind::SourceGrade &&
                    carrier->domain.extents().size() == 1,
                "scalar tower field requires a rank-one SourceGrade carrier");
        auto &prepared = context(name);
        const auto grade = carrier->grade.value;
        if (meaning.at("context").text() != name) {
            require(source.contexts.at(name).declaration.has("witness_scope"), "field/context mismatch");
            const auto &parent = context(meaning.at("context").text()).to_ambient.at(grade);
            const auto &child = prepared.to_ambient.at(grade);
            require(std::includes(parent.begin(), parent.end(), child.begin(), child.end()),
                    "field carrier does not contain the declared block context");
        }
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
    Json ensemble_query(const Json &query, bool validate, std::optional<NativeRecord> *native = nullptr) {
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
        if (validate) {
            (void)research_detail::relational_field_ensemble(*prepared.study, query,
                rfc::core::Grade{*grade}, values, true, resources, nullptr);
            return Json::object();
        }
        scientific_detail::ResultFields publication(owner->native_publication());
        scientific_detail::ResultFieldContext binding(publication, owner->native_publication(),
            prepared.to_ambient, std::pair<std::size_t, std::uint64_t>{*grade, *tensor(ids.at(0).text()).output});
        auto result = research_detail::relational_field_ensemble(
            *prepared.study, query, rfc::core::Grade{*grade}, values, false, resources, &binding);
        {
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
        const auto record = publication.finish(std::move(result));
        if (native) { *native = record; return Json::object(); }
        return scientific_detail::result_fields_view(record);
    }
    Json structural_calibration_query(const Json &query, bool validate, std::optional<NativeRecord> *native = nullptr) {
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
            auto occurrences = scientific_detail::ResultValue::array();
            for (std::size_t slot = begin; slot < end; ++slot) {
                if (arena.aspects.get(complex.signature().existence_aspect(), slot))
                    participants.insert(arena.lower_index[slot]);
                const StableOccurrenceKey occurrence{retained->grade,
                                                     arena.occurrence_identities[slot]};
                auto witness = scientific_detail::ResultValue(source.occurrence_record(occurrence));
                witness.set("coordinate", cardinality(arena.lower_index[slot]));
                witness.set("exists", Json::boolean(arena.aspects.get(
                                          complex.signature().existence_aspect(), slot)));
                witness.set(
                    "coefficient",
                    scientific_detail::ResultValue(NativeRecordValue(complex.incidence_coefficient(retained->grade, slot))));
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
            block.witness = scientific_detail::ResultValue(declaration);
            auto coordinates = scientific_detail::ResultValue::array();
            for (const auto coordinate : block.coordinates)
                coordinates.push(scientific_detail::ResultValue::count(coordinate));
            block.witness.set("coordinates", std::move(coordinates));
            block.witness.set("source_relation", source.cell_record(source.resolve_cell(relation)));
            block.witness.set("source_occurrences", occurrences);
            blocks.push_back(std::move(block));
        }
        if (validate) {
            (void)research_detail::structural_reference_calibration(
                *prepared.study, query, grade, values, blocks, true, resources, nullptr);
            return Json::object();
        }
        scientific_detail::ResultFields publication(owner->native_publication());
        scientific_detail::ResultFieldContext binding(publication, owner->native_publication(),
            prepared.to_ambient, std::pair<std::size_t, std::uint64_t>{grade.value, *tensor(first).output});
        auto result = research_detail::structural_reference_calibration(
            *prepared.study, query, grade, values, blocks, false, resources, &binding);
        {
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
        const auto record = publication.finish(std::move(result));
        if (native) {*native = record; return Json::object();}
        return scientific_detail::result_fields_view(record);
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
    Json local_declaration(const Json &query, bool legacy_budget = true,
        std::map<std::string, Vector> *native_operands = nullptr) {
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
                const auto axes = item.has("axes") ? item.at("axes") : ids;
                require(axes.size() == ids.size(), "section axis binding dimension");
                std::set<std::string> axis_names;
                for (std::size_t axis = 0; axis < axes.size(); ++axis)
                    require(axis_names.insert(text(axes.at(axis))).second, "duplicate section axis label");
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
                        require(stalk.at("axes").dump() == axes.dump(),
                                "source section axes must follow their declared field bindings");
                }
                Vector values;
                for (std::size_t c = 0; c < prepared.to_ambient.at(grade).size(); ++c)
                    for (const auto &v : inputs) values.push_back(v[c]);
                if (native_operands) {
                    require(native_operands->emplace(item.at("id").text(), std::move(values)).second,
                        "duplicate native section operand");
                } else {
                    auto rendered = Json::array();
                    for (const auto &value : values) rendered.push(Json::string(value.get_str()));
                    item.set("values", rendered);
                }
            }
            sections.push(item);
        }
        resolved.set("sections", sections);
        return resolved;
    }
    NativeRecord local_section_fields(const Json &query, bool validate = false) {
        const auto name = query.at("context").text();
        auto &prepared = context(name);
        std::map<std::string, Vector> operands;
        const auto declaration = local_declaration(query, false, &operands);
        if (validate) {
            (void)scientific_detail::analyze_local_section_fields(*prepared.study, declaration,
                nullptr, &operands, true);
            return owner->native_publication();
        }
        scientific_detail::ResultFields publication(owner->native_publication());
        auto analysis = scientific_detail::analyze_local_section_fields(*prepared.study, declaration,
            &publication, &operands);
        auto source_binding = scientific_detail::ResultValue::object();
        source_binding.set("kind", scientific_detail::ResultValue::string("local_section_source_fields"));
        auto sections = scientific_detail::ResultValue::array();
        for (std::size_t s = 0; s < query.at("sections").size(); ++s) {
            const auto item = query.at("sections").at(s);
            if (!item.has("fields")) continue;
            const auto grade = static_cast<std::size_t>(item.at("grade").integer());
            auto section = scientific_detail::ResultValue::object();
            section.set("id", item.at("id"));
            auto fields = scientific_detail::ResultValue::array();
            for (std::size_t a = 0; a < item.at("fields").size(); ++a) {
                const auto id = item.at("fields").at(a).text();
                const auto &original = tensor(id);
                const auto carrier = source.field_pins.at(original.key).carrier();
                auto field = scientific_detail::ResultValue::object();
                field.set("id", scientific_detail::ResultValue::string(id));
                field.set("field", publication.stored(original.key, "vector"));
                rfc::numeric::ExactMatrix realization(operands.at(item.at("id").text()).size(),
                    carrier->domain.logical_size());
                for (std::size_t c = 0; c < prepared.to_ambient.at(grade).size(); ++c)
                    realization.set(c * item.at("fields").size() + a,
                        prepared.to_ambient.at(grade)[c], Rational{1});
                field.set("realization", scientific_detail::ResultFields::stage_matrix(realization,
                    *original.output, analysis.grade_carriers.at(grade)));
                fields.push(std::move(field));
            }
            section.set("fields", std::move(fields));
            sections.push(std::move(section));
        }
        source_binding.set("sections", std::move(sections));
        publication.bind(std::move(source_binding));
        return publication.finish(std::move(analysis.view));
    }
    Json ablation_query(const Json &query, bool validate, std::optional<NativeRecord> *native = nullptr) {
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
            role.witness = scientific_detail::ResultValue(declaration);
            auto binding = scientific_detail::ResultValue::array();
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
                    auto record = scientific_detail::ResultValue(source.cell_record(source.resolve_cell(id)));
                    record.set("context_coordinate", scientific_detail::ResultValue::count(cell->index));
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
                    binding.push(scientific_detail::ResultValue(
                        source.fields.at(id).declaration.native().unbound_value()));
                }
            }
            role.witness.set("source_binding", binding);
            require(roles.emplace(text(declaration.at("id"), source_byte_limit()), std::move(role))
                        .second,
                    "duplicate scientific role");
        }
        auto witness = scientific_detail::ResultValue::object(), cells = scientific_detail::ResultValue::array(), occurrences = scientific_detail::ResultValue::array();
        for (std::size_t k = 0; k < prepared.source->grade_count(); ++k) {
            auto grade_cells = scientific_detail::ResultValue::array();
            const auto g = rfc::core::Grade{static_cast<std::uint16_t>(k)};
            const auto &arena = prepared.source->grade(g);
            for (std::size_t i = 0; i < arena.cell_count(); ++i) {
                auto record = scientific_detail::ResultValue(
                    source.cell_record(stable_cell_key(*prepared.source, {g, static_cast<std::uint32_t>(i)})));
                record.set("context_coordinate", scientific_detail::ResultValue::count(i));
                record.set("ambient_coordinate", scientific_detail::ResultValue::count(prepared.to_ambient[k][i]));
                grade_cells.push(record);
            }
            cells.push(grade_cells);
            for (std::size_t i = 0; i < arena.incidence_count(); ++i) {
                auto record =
                    scientific_detail::ResultValue(source.occurrence_record({g, arena.occurrence_identities[i]}));
                record.set("exists", scientific_detail::ResultValue::boolean(arena.aspects.get(
                                         prepared.source->signature().existence_aspect(), i)));
                record.set("coefficient",
                           scientific_detail::ResultValue(NativeRecordValue(prepared.source->incidence_coefficient(g, i))));
                occurrences.push(record);
            }
        }
        witness.set("cells", cells);
        witness.set("occurrences", occurrences);
        Json sheaf;
        std::map<std::string, Vector> sheaf_operands;
        if (query.has("sheaf")) {
            sheaf = Json::parse(query.at("sheaf").dump());
            members(sheaf, {"stalks", "maps", "sections", "metric"});
            sheaf.set("operation", Json::string("local_sections"));
            sheaf.set("context", Json::string(name));
            sheaf.set("construction", Json::string("cellular_sheaf"));
            sheaf = local_declaration(sheaf, false, &sheaf_operands);
        }
        if (validate) {
            (void)research_detail::source_structure_ablation(*prepared.study, query, grade, fields,
                roles, witness, sheaf, true, resources, owner->native_publication(), sheaf_operands, nullptr);
            return Json::object();
        }
        scientific_detail::ResultFields publication(owner->native_publication());
        research_detail::AblationFields axes(publication, owner->native_publication());
        const auto &source_axes = axes.source(*prepared.study);
        auto inputs = scientific_detail::ResultValue::array();
        for (const auto &[id, value] : fields) {
            const auto &original = tensor(id);
            const auto pin = source.field_pins.at(original.key).carrier();
            rfc::numeric::ExactMatrix realization(value.size(), pin->domain.logical_size());
            for (std::size_t i = 0; i < value.size(); ++i)
                realization.set(i, prepared.to_ambient.at(grade.value).at(i), Rational{1});
            auto input = scientific_detail::ResultValue::object();
            input.set("id", scientific_detail::ResultValue::string(id));
            input.set("source_field", publication.stored(original.key, "vector"));
            input.set("realization", scientific_detail::ResultFields::stage_matrix(realization,
                *original.output, source_axes.rows.at(grade.value)));
            input.set("context_input", scientific_detail::ResultFields::stage_vector(value,
                source_axes.rows.at(grade.value)));
            inputs.push(std::move(input));
        }
        auto input_binding = scientific_detail::ResultValue::object();
        input_binding.set("kind", scientific_detail::ResultValue::string("actual_ablation_source_fields"));
        input_binding.set("fields", std::move(inputs)); publication.bind(std::move(input_binding));
        auto view = research_detail::source_structure_ablation(*prepared.study, query, grade, fields,
            roles, witness, sheaf, false, resources, owner->native_publication(), sheaf_operands, &axes);
        const auto result = publication.finish(std::move(view));
        if (native) { *native = result; return Json{}; }
        return scientific_detail::result_fields_view(result);
    }
    NativeRecord field_summary_fields(const Json &query) {
        using namespace scientific_detail;
        auto reservation = cpu().reserve(1);
        const auto id = query.at("field").text();
        const auto &input_field = tensor(id);
        ResultFields publication(owner->native_publication());
        Vector values;
        std::vector<std::size_t> coordinates;
        input_field.tensor.support().for_each_active([&](std::size_t i) {
            coordinates.push_back(i); values.push_back(*input_field.tensor.state(i));
        });
        check_values(values);
        const Vector ones(values.size(), Rational{1});
        const auto sum = rfc::numeric::dot(ones, values);
        const auto squares = rfc::numeric::dot(values, values);
        const auto row_axis = publication.axis(coordinates.size());
        const auto scalar_view = [&](const Rational &value) {
            auto field = publication.stage_scalar(value); field.set("display", ResultValue::string("rational_object")); return field;
        };
        auto view = ResultValue::object();
        view.set("reported", ResultValue::count(coordinates.size()));
        view.set("missing", ResultValue::count(input_field.tensor.domain().logical_size() - coordinates.size()));
        view.set("coordinates", publication.indices(coordinates, row_axis));
        view.set("sum", scalar_view(sum));
        view.set("sum_of_squares", scalar_view(squares));
        if (values.empty()) {
            rfc::numeric::SparseRationalTensor absent(rfc::numeric::IndexDomain::linear(1),
                rfc::numeric::StructuralSupport::sparse(1, {}), rfc::numeric::SparseRationalVector(1));
            view.set("mean", ResultFields::staged(std::move(absent), {}, publication.axis(1), "rational_object"));
        } else view.set("mean", scalar_view(sum / Rational(values.size())));
        view.set("units", ResultValue(source.fields.at(id).declaration.at("units").native().unbound_value()));
        auto binding = ResultValue::object();
        binding.set("kind", ResultValue::string("actual_reported_field_summary"));
        binding.set("id", ResultValue::string(id));
        binding.set("source_field", publication.stored(input_field.key, "vector"));
        binding.set("row_carrier", ResultValue::count(row_axis));
        binding.set("row_coordinates", view.at("coordinates"));
        binding.set("reported_values", ResultFields::stage_vector(values, row_axis));
        binding.set("weights", ResultFields::stage_vector(ones, row_axis));
        publication.bind(std::move(binding));
        return publication.finish(std::move(view));
    }
    NativeRecord direct_action_fields(const Json &query) {
        using namespace scientific_detail;
        const auto operation = query.at("operation").text();
        const auto name = query.at("context").text(), id = query.at("field").text();
        auto &prepared = context(name);
        const auto [values, grade] = input(name, id, false);
        auto reservation = cpu().reserve(1);
        const auto lease = source_execution_lease("scientific-direct-source-action");
        lease.check_cancelled();
        const auto snapshot = source_snapshot_publication(*prepared.study, owner->native_publication(), prepared.to_ambient);
        ResultFields publication(snapshot.publication);
        auto output_grade = grade;
        rfc::field::OperatorRequest request;
        request.grade = rfc::core::Grade{grade}; request.kind = rfc::field::OperatorKind::Hodge;
        if (operation == "boundary") {request.kind = rfc::field::OperatorKind::Boundary; --output_grade;}
        if (operation == "adjoint") {request.kind = rfc::field::OperatorKind::Adjoint; ++request.grade.value; ++output_grade;}
        Vector output;
        std::optional<rfc::field::ExactGreenResult> green;
        if (operation == "green_action") {
            const auto lambda = exact_rational(query.at("lambda").text());
            green = prepared.study->green_action(rfc::core::Grade{grade}, lambda, values);
            require(green->exact && rfc::numeric::is_zero(rfc::numeric::sub(
                rfc::numeric::add(green->value, rfc::numeric::scale(lambda,
                    prepared.study->hodge(rfc::core::Grade{grade}).apply(green->value))), values)),
                "native direct Green reconstruction residual");
            output = green->value;
        } else output = prepared.study->execute(request, values);
        lease.check_cancelled();
        const auto output_field = publication.vector(output, snapshot.grade_axes.at(output_grade));
        auto view = ResultValue::object(), result = ResultValue::array();
        for (std::size_t i = 0; i < output.size(); ++i) {
            const auto ambient = prepared.to_ambient.at(output_grade).at(i);
            const auto cell = source.coordinate(output_grade, ambient);
            auto record = ResultValue::object();
            record.set("coordinate", ResultValue::count(i));
            record.set("native_key", ResultValue::string(std::to_string(cell.key)));
            record.set("entity", ResultValue(source.cell_record(cell).at("id")));
            auto value = output_field;
            value.set("display", ResultValue::string("rational_object"));
            value.set("index", ResultValue::count(i));
            record.set("value", std::move(value)); result.push(std::move(record));
        }
        view.set("schema", ResultValue::string("gadget.scientific.direct.result.v1"));
        view.set("result", std::move(result)); view.set("output_grade", ResultValue::count(output_grade));
        if (green) {
            view.set("green_iterations", ResultValue::count(green->iterations));
            view.set("green_execution", ResultValue::string("rfc::field::exact_green_apply; factored boundary and metric actions"));
            view.set("residual_verified", ResultValue::boolean(true));
        }
        auto binding = ResultValue::object();
        const auto &input_tensor = tensor(id);
        binding.set("kind", ResultValue::string("actual_direct_source_action"));
        binding.set("operation", ResultValue::string(operation));
        binding.set("field", ResultValue::string(id));
        binding.set("input_grade", ResultValue::count(grade));
        binding.set("output_grade", ResultValue::count(output_grade));
        binding.set("original_study_publication", ResultValue(snapshot.publication.at("bindings").at(0)
            .at("original_study_publication").unbound_value()));
        binding.set("input_source_namespace", ResultValue::string(std::to_string(source.field_pins.at(input_tensor.key).source_namespace())));
        binding.set("input_cutoff", ResultValue::count(source.field_pins.at(input_tensor.key).sequence()));
        binding.set("input_tensor", ResultValue::count(input_tensor.key));
        binding.set("input_carrier", ResultValue::count(*input_tensor.output));
        const auto input_axis = publication.axis(input_tensor.tensor.domain().logical_size());
        binding.set("input_coordinate_carrier", ResultValue::count(input_axis));
        rfc::numeric::ExactMatrix realization(values.size(), input_tensor.tensor.domain().logical_size());
        for (std::size_t i = 0; i < values.size(); ++i) realization.set(i, prepared.to_ambient.at(grade).at(i), Rational{1});
        binding.set("input_realization", ResultFields::stage_matrix(realization, input_axis, snapshot.grade_axes.at(grade)));
        binding.set("context_input", ResultFields::stage_vector(values, snapshot.grade_axes.at(grade)));
        binding.set("output", output_field);
        binding.set("execution_scope", ResultValue::string("Canonical one shot action; CPU reservation and cancellation before and after the action"));
        publication.bind(std::move(binding));
        return publication.finish(std::move(view));
    }
    NativeRecord support_fields(const Json &query) {
        using namespace scientific_detail;
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
        const auto prepared = packed.at(identity)->prepared();
        auto derived = packed.at(identity)->execute(cpu(), query);
        derived.execution.set("prepared_operation_reused", Json::boolean(reused));
        ResultFields publication(owner->native_publication());
        const auto carrier = *first.output, size = first.tensor.domain().logical_size();
        const auto bits = [&](const rfc::machine::PackedBits &mask) {
            require(mask.size() == size, "packed native field domain");
            std::vector<rfc::numeric::SparseRationalEntry> entries;
            mask.for_each_set([&](auto i) {entries.push_back({i,Rational{1}});});
            return ResultFields::staged(rfc::numeric::SparseRationalTensor(
                rfc::numeric::IndexDomain::linear(size),rfc::numeric::StructuralSupport::full(size),
                rfc::numeric::SparseRationalVector(size,std::move(entries))),{},carrier,"vector");
        };
        const auto words = [](const rfc::machine::PackedBits &mask) {
            auto values = ResultValue::array();
            for (const auto word : mask.words()) values.push(ResultValue::count(word));
            return values;
        };
        auto binding = ResultValue::object(), inputs = ResultValue::array();
        binding.set("kind",ResultValue::string("actual_packed_support_operation"));
        binding.set("source_carrier",ResultValue::count(carrier));
        binding.set("mode",ResultValue::string(query.at("mode").text()));
        binding.set("slots",ResultValue::count(size));
        binding.set("register_count",ResultValue::count(prepared.program().register_count()));
        for (std::size_t i = 0; i < fields.size(); ++i) {
            const auto id = fields.at(i).text();
            const auto &field = tensor(id); const auto &pin = source.field_pins.at(field.key);
            auto input = ResultValue::object();
            input.set("id",ResultValue::string(id)); input.set("tensor",ResultValue::count(field.key));
            input.set("source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));
            input.set("cutoff",ResultValue::count(pin.sequence())); input.set("carrier",ResultValue::count(*field.output));
            input.set("source_field",publication.stored(field.key,"vector"));
            input.set("support_plane",bits(prepared.planes().plane(static_cast<AspectId>(i))));
            input.set("packed_words",words(prepared.planes().plane(static_cast<AspectId>(i))));
            input.set("meaning",ResultValue(source.fields.at(id).declaration)); inputs.push(std::move(input));
        }
        binding.set("inputs",std::move(inputs));
        auto instructions = ResultValue::array(), outputs = ResultValue::array(), masks = ResultValue::array();
        auto output_words = ResultValue::array();
        for (const auto &instruction : prepared.program().instructions()) {
            auto row = ResultValue::object();
            row.set("opcode",ResultValue::count(static_cast<std::size_t>(instruction.opcode)));
            row.set("dst",ResultValue::count(instruction.dst)); row.set("a",ResultValue::count(instruction.a));
            row.set("b",ResultValue::count(instruction.b)); row.set("c",ResultValue::count(instruction.c));
            row.set("aspect",ResultValue::count(instruction.aspect)); instructions.push(std::move(row));
        }
        for (std::size_t i = 0; i < prepared.outputs().size(); ++i) {
            outputs.push(ResultValue::count(prepared.outputs()[i])); masks.push(bits(derived.outputs.at(i)));
            output_words.push(words(derived.outputs.at(i)));
        }
        binding.set("instructions",std::move(instructions)); binding.set("outputs",std::move(outputs));
        binding.set("output_masks",std::move(masks)); binding.set("output_words",std::move(output_words));
        auto view = ResultValue::object(), rows = ResultValue::array();
        view.set("coordinates",ResultValue::count(size));
        view.set("selected",ResultValue::count(derived.outputs.at(0).count()));
        view.set("native_carrier",ResultValue::string(std::to_string(carrier)));
        const auto source_carrier = temporal_endpoint_detail::carriers(*source.directory->pin_snapshot()).at(carrier);
        derived.outputs.at(0).for_each_set([&](auto i) {
            const auto cell = source.coordinate(source_carrier->grade.value,i);
            auto row = ResultValue::object(); row.set("coordinate",ResultValue::count(i));
            row.set("entity",ResultValue(source.cell_record(cell).at("id")));
            row.set("native_key",ResultValue::string(std::to_string(cell.key))); rows.push(std::move(row));
        });
        view.set("selection",std::move(rows)); publication.bind(std::move(binding));
        const auto native = publication.finish(std::move(view));
        auto retained = native.unbound_value();
        std::get<NativeRecordValue::Object>(retained.data).members.emplace_back(
            "execution",import_record_view(derived.execution));
        return native.derive(std::move(retained));
    }
    NativeRecord window_fields(const Json &query,
                              std::optional<rfc::runtime::ExecutionLease> native_lease = {}) {
        using namespace scientific_detail;
        const auto id = query.at("field").text();
        const auto &input_tensor = tensor(id);
        const auto pin = source.field_pins.at(input_tensor.key);
        const rfc::runtime::WorkAvailability availability{pin.source_namespace(),
            publication_cutoff(query.at("available_sequence"))};
        auto key_query = Json::parse(query.dump()); key_query.set("preparation", Json::string("reuse"));
        const auto cache_key = key_query.dump();
        const bool retain = !query.has("preparation") || query.at("preparation").text() == "reuse";
        const bool reused = retain && windows.contains(cache_key);
        const auto lease = native_lease.value_or(source_execution_lease("scientific-tensor-window"));
        std::optional<rfc::database::NativeTensorWindow> transient;
        if (!reused) {
            std::vector<std::size_t> lo, hi;
            for (std::size_t a = 0; a < query.at("lower").size(); ++a) {
                lo.push_back(static_cast<std::size_t>(query.at("lower").at(a).integer()));
                hi.push_back(static_cast<std::size_t>(query.at("upper").at(a).integer()));
            }
            auto window = rfc::database::NativeTensorWindow::prepare(pin,lo,hi,availability,lease);
            if (retain) windows.emplace(cache_key,std::move(window)); else transient.emplace(std::move(window));
        }
        const auto &window = retain ? windows.at(cache_key) : *transient;
        const auto summary = window.summarize(availability,lease);
        ResultFields publication(owner->native_publication());
        auto view = ResultValue::object(), observations = ResultValue::array();
        std::vector<std::size_t> source_coordinates;
        const auto size = window.domain().logical_size(), window_axis = publication.axis(size);
        const auto original = publication.stored(input_tensor.key,"vector");
        for (std::size_t i = 0; i < size; ++i) {
            lease.check_cancelled(); const auto coordinate = window.source_index(i);
            source_coordinates.push_back(coordinate);
            auto row = ResultValue::object(), value = original;
            row.set("source_coordinate",ResultValue::count(coordinate));
            row.set("reported",ResultValue::boolean(window.tensor_pin()->state(coordinate).has_value()));
            value.set("display",ResultValue::string("rational_object")); value.set("index",ResultValue::count(coordinate));
            row.set("value",std::move(value)); observations.push(std::move(row));
        }
        view.set("observations",std::move(observations));
        view.set("reported",ResultValue::count(window.supported_count()));
        view.set("missing",ResultValue::count(size-window.supported_count()));
        ResultValue sum;
        if (summary.sum()) sum=publication.stage_scalar(*summary.sum());
        else sum=ResultFields::staged(rfc::numeric::SparseRationalTensor(
            rfc::numeric::IndexDomain::linear(1),rfc::numeric::StructuralSupport::sparse(1,{}),
            rfc::numeric::SparseRationalVector(1)),{},publication.axis(1),"rational_object");
        sum.set("display",ResultValue::string("rational_object")); view.set("sum",std::move(sum));
        view.set("lower",ResultValue(query.at("lower"))); view.set("upper",ResultValue(query.at("upper")));
        view.set("source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));
        view.set("source_sequence",ResultValue::string(std::to_string(pin.sequence())));
        auto binding=ResultValue::object(); binding.set("kind",ResultValue::string("actual_native_tensor_window"));
        binding.set("field",ResultValue::string(id)); binding.set("input_tensor",ResultValue::count(input_tensor.key));
        binding.set("input_carrier",ResultValue::count(*input_tensor.output));
        binding.set("source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));
        binding.set("cutoff",ResultValue::count(pin.sequence())); binding.set("availability",ResultValue::count(availability.max_sequence));
        binding.set("source_field",original); binding.set("window_carrier",ResultValue::count(window_axis));
        binding.set("source_coordinates",publication.indices(source_coordinates,window_axis));
        auto domain=ResultValue::array(); for (const auto extent:window.domain().extents()) domain.push(ResultValue::count(extent));
        binding.set("window_extents",std::move(domain)); publication.bind(std::move(binding));
        const auto native = publication.finish(std::move(view));
        auto retained = native.unbound_value();
        NativeRecordValue metadata(NativeRecordValue::Object{{
            {"preparation_reused",NativeRecordValue(reused)},
            {"preparation_retained",NativeRecordValue(retain)}}});
        std::get<NativeRecordValue::Object>(retained.data).members.emplace_back("view_execution",std::move(metadata));
        return native.derive(std::move(retained));
    }
    NativeRecord structure_transition_fields(const Json &query,
        std::optional<rfc::runtime::ExecutionLease> requested_lease = {}) {
        using namespace scientific_detail;
        const auto lease = requested_lease ? *requested_lease : source_execution_lease("scientific-native-structure-transition");
        lease.check_cancelled();
        const auto old_name=query.at("context").text(), new_name=query.at("other_context").text();
        auto &old_context=context(old_name); auto &new_context=context(new_name);
        const auto [old_values,grade]=input(old_name,query.at("field").text(),false);
        const auto new_values=input(new_name,query.at("new_field").text(),false).first;
        const auto lambda=exact_rational(query.at("lambda").text());
        std::size_t slots=0,entries=0,bits=0,dimension=0;
        for (const auto &id : {query.at("field").text(),query.at("new_field").text()}) {
            const auto &field=tensor(id);
            slots=rfc::machine::checked_size_add(slots,field.tensor.support().count());
            entries=rfc::machine::checked_size_add(entries,field.tensor.stored_nonzeros());
            for (const auto &v : field.tensor.values().entries())
                bits=rfc::machine::checked_size_add(bits,rfc::storage::native_execution_detail::integer_bits(v.value));
        }
        for (const auto *context : {&old_context,&new_context})
            for (const auto &coordinates : context->to_ambient) dimension=std::max(dimension,coordinates.size());
        rfc::runtime::admit_execution({2,slots,entries,1,bits},lease.policy());
        const auto &policy=lease.policy().limits;
        const rfc::numeric::ExactActionLimits limits{dimension,policy.sparse_entries,policy.scalar_work,policy.integer_bits};
        const auto checkpoint=[lease] {lease.check_cancelled();};
        rfc::numeric::ExactActionGuard guard(limits,checkpoint);
        guard.vector(old_values);guard.vector(new_values);guard.value(lambda);
        auto reservation=cpu().reserve(1);
        const auto admit_matrix=[&](const rfc::numeric::ExactMatrix &matrix) {
            research_detail::admit_structure_matrix(matrix,lease);
        };
        const auto old_snapshot=source_snapshot_publication(*old_context.study,owner->native_publication(),old_context.to_ambient,admit_matrix);
        const auto new_snapshot=source_snapshot_publication(*new_context.study,owner->native_publication(),new_context.to_ambient,admit_matrix);
        ResultFields publication(new_snapshot.publication);
        research_detail::StructuralResultFields admission{publication,{},new_snapshot.grade_axes,lease};
        auto binding=ResultValue::object(), old_grades=ResultValue::array();
        binding.set("kind",ResultValue::string("actual_structure_transition_endpoints"));
        binding.set("old_source",publication.publication(old_snapshot.publication));
        binding.set("new_source",publication.publication(new_snapshot.publication));
        binding.set("original_study_publication",publication.publication(owner->native_publication()));
        std::vector<std::uint64_t> old_axes;
        for (std::size_t k=0;k<old_context.to_ambient.size();++k) {
            const auto axis=publication.axis(old_context.to_ambient[k].size());old_axes.push_back(axis);
            auto record=ResultValue::object();record.set("grade",ResultValue::count(k));
            record.set("source_carrier",ResultValue::count(old_snapshot.grade_axes[k]));
            record.set("coordinate_carrier",ResultValue::count(axis));
            record.set("original_source_namespace",ResultValue::string(std::to_string(source.grade_pins.at(k).source_namespace())));
            record.set("original_source_cutoff",ResultValue::count(source.grade_pins.at(k).sequence()));
            std::vector<std::size_t> coordinates(old_context.to_ambient[k].begin(),old_context.to_ambient[k].end());
            record.set("original_coordinates",publication.indices(coordinates,axis));old_grades.push(std::move(record));
        }
        binding.set("old_grades",std::move(old_grades));
        const auto retain_input=[&](const std::string &id,const Vector &values,std::uint64_t axis,
            const std::vector<std::uint32_t> &coordinates) {
            const auto &field=tensor(id);const auto &pin=source.field_pins.at(field.key);
            auto input=ResultValue::object();input.set("id",ResultValue::string(id));
            input.set("source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));
            input.set("source_cutoff",ResultValue::count(pin.sequence()));
            input.set("source_tensor",ResultValue::count(field.key));input.set("source_carrier",ResultValue::count(*field.output));
            input.set("coordinate_carrier",ResultValue::count(axis));
            const auto original_axis=publication.axis(field.tensor.domain().logical_size());
            input.set("original_coordinate_carrier",ResultValue::count(original_axis));
            rfc::numeric::ExactMatrix lift(values.size(),field.tensor.domain().logical_size());
            for (std::size_t i=0;i<values.size();++i)lift.set(i,coordinates.at(i),Rational{1});
            admission.admit(lift);admission.admit(values);
            input.set("realization",ResultFields::stage_matrix(lift,original_axis,axis));
            input.set("values",ResultFields::stage_vector(values,axis));
            input.set("meaning",ResultValue(source.fields.at(id).declaration));return input;
        };
        binding.set("old_input",retain_input(query.at("field").text(),old_values,old_axes.at(grade),old_context.to_ambient.at(grade)));
        binding.set("new_input",retain_input(query.at("new_field").text(),new_values,new_snapshot.grade_axes.at(grade),new_context.to_ambient.at(grade)));
        binding.set("old_time",ResultValue(source.contexts.at(old_name).declaration.at("time")));
        binding.set("new_time",ResultValue(source.contexts.at(new_name).declaration.at("time")));
        admission.admit(lambda);binding.set("lambda",publication.stage_scalar(lambda));publication.bind(std::move(binding));
        auto definitions=ResultValue::object();
        definitions.set("kind",ResultValue::string("declared_structure_defect_meanings"));
        definitions.set("green_operator_defect",ResultValue::string("L_new J - J L_old"));
        definitions.set("green_response_defect",ResultValue::string("G_new J - J G_old, with G=(I+lambda L)^-1"));
        definitions.set("operator_term",ResultValue::string("(G_new J - J G_old) applied to the original old field"));
        publication.bind(std::move(definitions));
        std::unique_ptr<PreparedTemporalStudy> leased_old,leased_new;
        const auto prepare=[&](PreparedContext &context,std::unique_ptr<PreparedTemporalStudy> &leased) -> PreparedTemporalStudy & {
            if (requested_lease) {
                leased=std::make_unique<PreparedTemporalStudy>(*context.study,limits,checkpoint);
                return *leased;
            }
            if (!context.temporal)
                context.temporal=std::make_unique<PreparedTemporalStudy>(*context.study,limits,checkpoint);
            return *context.temporal;
        };
        auto &old_prepared=prepare(old_context,leased_old),&new_prepared=prepare(new_context,leased_new);
        research_detail::StructuralResultFields fields{publication,old_axes,new_snapshot.grade_axes,lease};
        auto result=research_detail::compare_prepared_structure_fields(old_prepared,new_prepared,rfc::core::Grade{grade},
            old_context.to_ambient,new_context.to_ambient,old_values,new_values,lambda,nullptr,
            research_detail::FactoredTemporalMaps::requested(query),fields);
        lease.check_cancelled();return publication.finish(std::move(result));
    }
    NativeRecord moment_fields(const Json &query) {
        using namespace scientific_detail;
        const auto name = query.at("context").text();
        auto &old_context = context(name);
        const auto [old_x, grade] = input(name, query.at("field").text(), false);
        const auto old_y = input(name, query.at("other_field").text(), false).first;
        const auto order = count(query.at("order"), std::numeric_limits<std::size_t>::max(), "moment order");
        auto reservation = cpu().reserve(1);
        const auto lease = source_execution_lease("scientific-native-moment");
        lease.check_cancelled();
        const auto old_snapshot = source_snapshot_publication(*old_context.study,
            owner->native_publication(), old_context.to_ambient);
        const auto retain_input = [&](ResultFields &publication, const std::string &id,
                                      const Vector &values, std::uint64_t context_axis,
                                      const std::vector<std::uint32_t> &coordinates) {
            const auto &source_tensor = tensor(id);
            const auto &pin = source.field_pins.at(source_tensor.key);
            auto binding = ResultValue::object();
            binding.set("id", ResultValue::string(id));
            binding.set("tensor", ResultValue::count(source_tensor.key));
            binding.set("source_namespace", ResultValue::string(std::to_string(pin.source_namespace())));
            binding.set("cutoff", ResultValue::count(pin.sequence()));
            binding.set("carrier", ResultValue::count(*source_tensor.output));
            const auto axis = publication.axis(source_tensor.tensor.domain().logical_size());
            binding.set("original_coordinate_carrier", ResultValue::count(axis));
            binding.set("context_carrier", ResultValue::count(context_axis));
            rfc::numeric::ExactMatrix realization(values.size(), source_tensor.tensor.domain().logical_size());
            for (std::size_t i = 0; i < values.size(); ++i) realization.set(i, coordinates.at(i), Rational{1});
            binding.set("realization", ResultFields::stage_matrix(realization, axis, context_axis));
            binding.set("values", ResultFields::stage_vector(values, context_axis));
            binding.set("meaning", ResultValue(source.fields.at(id).declaration));
            return binding;
        };
        const auto scalar = [](ResultFields &publication, const Rational &value) {
            auto field = publication.stage_scalar(value);
            field.set("display", ResultValue::string("rational_object")); return field;
        };
        if (query.at("operation").text() == "typed_moment") {
            ResultFields publication(old_snapshot.publication);
            const auto axis = old_snapshot.grade_axes.at(grade);
            const auto hodge = old_context.study->hodge(rfc::core::Grade{grade});
            const auto powered = rfc::statistics::hodge_power_apply(hodge, order, old_y);
            const auto value = hodge.metric().inner(old_x, powered);
            auto binding = ResultValue::object(), inputs = ResultValue::array();
            binding.set("kind", ResultValue::string("actual_typed_moment"));
            binding.set("original_study_publication", ResultValue(old_snapshot.publication.at("bindings").at(0)
                .at("original_study_publication").unbound_value()));
            binding.set("grade", ResultValue::count(grade)); binding.set("order", ResultValue::count(order));
            inputs.push(retain_input(publication, query.at("field").text(), old_x, axis, old_context.to_ambient.at(grade)));
            inputs.push(retain_input(publication, query.at("other_field").text(), old_y, axis, old_context.to_ambient.at(grade)));
            binding.set("inputs", std::move(inputs));
            binding.set("powered_right_field", ResultFields::stage_vector(powered, axis));
            publication.bind(std::move(binding));
            lease.check_cancelled();
            auto view = ResultValue::object();
            view.set("rational", publication.stage_scalar(value));
            return publication.finish(std::move(view));
        }
        const auto other_name = query.at("other_context").text();
        auto &new_context = context(other_name);
        const auto new_x = input(other_name, query.at("new_field").text(), false).first;
        const auto new_y = input(other_name, query.at("new_other_field").text(), false).first;
        const auto new_snapshot = source_snapshot_publication(*new_context.study,
            owner->native_publication(), new_context.to_ambient);
        ResultFields publication(new_snapshot.publication);
        const auto old_axis = publication.axis(old_x.size()), new_axis = new_snapshot.grade_axes.at(grade);
        const auto factored = research_detail::NativeFactoredTemporalMaps::requested(query);
        const auto form = [&](const NativeStudy &study) {
            const auto hodge = study.hodge(rfc::core::Grade{grade});
            return rfc::field::materialize_action_matrix(hodge.metric(), [&](const Vector &v) {
                lease.check_cancelled();
                return hodge.metric().apply(rfc::statistics::hodge_power_apply(hodge, order, v));
            });
        };
        const auto comparison = query.has("correspondence")
            ? temporal_map(name, other_name, grade, query.at("correspondence"))
            : rfc::numeric::ExactMatrix::identity(old_x.size());
        const auto transported_x = comparison.apply(old_x), transported_y = comparison.apply(old_y);
        const auto innovation_x = rfc::numeric::sub(new_x, transported_x),
                   innovation_y = rfc::numeric::sub(new_y, transported_y);
        const auto old_value = old_context.study->moment(rfc::core::Grade{grade}, order, old_x, old_y);
        const auto new_value = new_context.study->moment(rfc::core::Grade{grade}, order, new_x, new_y);
        // These two columns label the actual transported field and innovation.
        // They need not be independent and do not replace either source carrier.
        const auto restrict_form = [&](const NativeStudy &study, const std::vector<Vector> &left,
                                       const std::vector<Vector> &right) {
            const auto hodge = study.hodge(rfc::core::Grade{grade});
            rfc::numeric::ExactMatrix restricted(left.size(), right.size());
            for (std::size_t j = 0; j < right.size(); ++j) {
                lease.check_cancelled();
                const auto powered = rfc::statistics::hodge_power_apply(hodge, order, right[j]);
                for (std::size_t i = 0; i < left.size(); ++i)
                    restricted.set(i, j, hodge.metric().inner(left[i], powered));
            }
            return restricted;
        };
        rfc::numeric::ExactMatrix old_contraction(1, 1), diagnostic_map(2, 1);
        old_contraction.set(0, 0, old_value); diagnostic_map.set(0, 0, Rational{1});
        const auto new_contraction = restrict_form(*new_context.study,
            {transported_x, innovation_x}, {transported_y, innovation_y});
        const auto change = rfc::temporal::changing_dimension_cross_moment(
            Vector{Rational{1}}, Vector{Rational{1}}, old_contraction,
            Vector{Rational{1}, Rational{1}}, Vector{Rational{1}, Rational{1}},
            new_contraction, diagnostic_map, diagnostic_map);
        require(change.total == new_value - old_value, "native finite cross moment identity failed");
        const auto old_time = source.contexts.at(name).declaration.at("time"),
                   new_time = source.contexts.at(other_name).declaration.at("time");
        const Rational interval = exact_decimal(new_time.at("position").text()) - exact_decimal(old_time.at("position").text());
        auto view = ResultValue::object();
        view.set("pullback_form", scalar(publication, change.pullback_form));
        view.set("left_innovation", scalar(publication, change.left_innovation));
        view.set("right_innovation", scalar(publication, change.right_innovation));
        view.set("innovation_interaction", scalar(publication, change.innovation_interaction));
        view.set("total", scalar(publication, change.total));
        view.set("identity_verified", ResultValue::boolean(true));
        view.set("interval", scalar(publication, interval));
        view.set("time_unit", ResultValue(old_time.at("unit")));
        view.set("finite_rate", scalar(publication, rfc::temporal::finite_rate(change.total, interval)));
        auto binding = ResultValue::object(), inputs = ResultValue::array();
        binding.set("kind", ResultValue::string("actual_temporal_moment_change"));
        binding.set("old_context_publication", publication.publication(old_snapshot.publication));
        binding.set("original_study_publication", ResultValue(new_snapshot.publication.at("bindings").at(0)
            .at("original_study_publication").unbound_value()));
        binding.set("old_grade_carrier", ResultValue::count(old_snapshot.grade_axes.at(grade)));
        binding.set("old_coordinate_carrier", ResultValue::count(old_axis));
        binding.set("new_grade_carrier", ResultValue::count(new_axis));
        binding.set("grade", ResultValue::count(grade)); binding.set("order", ResultValue::count(order));
        binding.set("old_time", ResultValue(old_time)); binding.set("new_time", ResultValue(new_time));
        inputs.push(retain_input(publication, query.at("field").text(), old_x, old_axis, old_context.to_ambient.at(grade)));
        inputs.push(retain_input(publication, query.at("other_field").text(), old_y, old_axis, old_context.to_ambient.at(grade)));
        inputs.push(retain_input(publication, query.at("new_field").text(), new_x, new_axis, new_context.to_ambient.at(grade)));
        inputs.push(retain_input(publication, query.at("new_other_field").text(), new_y, new_axis, new_context.to_ambient.at(grade)));
        binding.set("inputs", std::move(inputs));
        if (factored) {
            std::vector<std::uint64_t> old_axes;
            auto axes = ResultValue::array();
            for (std::size_t k = 0; k < old_context.to_ambient.size(); ++k) {
                const auto axis = k == grade ? old_axis : publication.axis(old_context.to_ambient[k].size());
                old_axes.push_back(axis);
                auto item = ResultValue::object();
                item.set("grade", ResultValue::count(k));
                item.set("source_carrier", ResultValue::count(old_snapshot.grade_axes[k]));
                item.set("coordinate_carrier", ResultValue::count(axis));
                std::vector<std::size_t> coordinates(old_context.to_ambient[k].begin(), old_context.to_ambient[k].end());
                item.set("original_coordinates", publication.indices(coordinates, axis));
                axes.push(std::move(item));
            }
            binding.set("old_endpoint_axes", std::move(axes));
            research_detail::NativeFactoredTemporalMaps registry(
                {&temporal(old_context), &temporal(new_context)}, &publication,
                {old_axes, new_snapshot.grade_axes});
            const auto expression = [&](std::size_t endpoint, std::uint64_t axis, std::size_t size) {
                auto power = registry.identity(size, axis);
                if (order) {
                    auto factor = registry.hodge(endpoint, rfc::core::Grade{grade});
                    for (auto remaining = order; remaining; remaining >>= 1u) {
                        if (remaining & 1u) power = registry.compose({power, factor});
                        if (remaining > 1u) factor = registry.compose({factor, factor});
                    }
                }
                return registry.compose({power, registry.metric(endpoint, rfc::core::Grade{grade})});
            };
            const auto old_form = expression(0, old_axis, old_x.size()),
                       new_form = expression(1, new_axis, new_x.size());
            binding.set("old_form", old_form); binding.set("new_form", new_form);
            binding.set("factored_maps", registry.archive());
            auto forms = ResultValue::object();
            forms.set("old_form", old_form); forms.set("new_form", new_form);
            forms.set("factored_maps", registry.archive());
            view.set("forms", std::move(forms));
        } else {
            binding.set("old_form", ResultFields::stage_matrix(form(*old_context.study), old_axis, old_axis));
            binding.set("new_form", ResultFields::stage_matrix(form(*new_context.study), new_axis, new_axis));
        }
        binding.set("comparison", ResultFields::stage_matrix(comparison, old_axis, new_axis));
        binding.set("transported_left", ResultFields::stage_vector(transported_x, new_axis));
        binding.set("transported_right", ResultFields::stage_vector(transported_y, new_axis));
        binding.set("left_innovation", ResultFields::stage_vector(innovation_x, new_axis));
        binding.set("right_innovation", ResultFields::stage_vector(innovation_y, new_axis));
        std::optional<rfc::numeric::ExactMatrix> reverse;
        if (query.has("reverse_correspondence"))
            reverse = temporal_map(other_name, name, grade, query.at("reverse_correspondence"));
        else if (old_context.to_ambient.at(grade) == new_context.to_ambient.at(grade) &&
                 comparison == rfc::numeric::ExactMatrix::identity(old_x.size()))
            reverse = rfc::numeric::ExactMatrix::identity(old_x.size());
        if (reverse) {
            require(reverse->multiply(comparison) == rfc::numeric::ExactMatrix::identity(old_x.size()) &&
                    comparison.multiply(*reverse) == rfc::numeric::ExactMatrix::identity(new_x.size()),
                    "moment common carrier requires two sided inverse correspondences");
            const auto dx = rfc::numeric::sub(reverse->apply(new_x), old_x),
                       dy = rfc::numeric::sub(reverse->apply(new_y), old_y);
            const auto old_restricted = restrict_form(*old_context.study, {old_x, dx}, {old_y, dy});
            const auto new_restricted = restrict_form(*new_context.study,
                {transported_x, comparison.apply(dx)}, {transported_y, comparison.apply(dy)});
            rfc::numeric::ExactMatrix old_left(1, 2), new_left(1, 2);
            old_left.set(0, 0, Rational{1}); new_left.set(0, 0, Rational{1}); new_left.set(0, 1, Rational{1});
            const auto terms = rfc::temporal::finite_ordered_factor_attribution(
                {old_left, old_restricted, old_left.transpose()},
                {new_left, new_restricted, new_left.transpose()}, rfc::temporal::FactorOrder::WrittenProduct);
            require(rfc::temporal::sum_factor_interactions(terms, 1, 1).get(0, 0) == change.total,
                    "native seven term moment identity failed");
            auto interactions = ResultValue::array();
            for (const auto &term : terms) {
                auto item = ResultValue::object(), factors = ResultValue::array();
                for (auto position : term.changed_positions)
                    factors.push(ResultValue::string(position == 0 ? "left_field" : position == 1 ? "form" : "right_field"));
                item.set("changed_factors", std::move(factors));
                item.set("contribution", scalar(publication, term.contribution.get(0, 0)));
                interactions.push(std::move(item));
            }
            view.set("factor_interactions", std::move(interactions));
            auto common = ResultValue::object();
            common.set("reverse_correspondence", ResultFields::stage_matrix(*reverse, new_axis, old_axis));
            common.set("left_change", ResultFields::stage_vector(dx, old_axis));
            common.set("right_change", ResultFields::stage_vector(dy, old_axis));
            const auto left_axis = publication.axis(2), right_axis = publication.axis(2);
            const auto lift = [](const Vector &base, const Vector &delta) {
                rfc::numeric::ExactMatrix result(base.size(), 2);
                for (std::size_t i = 0; i < base.size(); ++i) {
                    result.set(i, 0, base[i]); result.set(i, 1, delta[i]);
                }
                return result;
            };
            common.set("left_realization", ResultFields::stage_matrix(lift(old_x, dx), left_axis, old_axis));
            common.set("right_realization", ResultFields::stage_matrix(lift(old_y, dy), right_axis, old_axis));
            common.set("old_restricted_form", ResultFields::stage_matrix(old_restricted, right_axis, left_axis));
            common.set("new_restricted_form", ResultFields::stage_matrix(new_restricted, right_axis, left_axis));
            common.set("meaning", ResultValue::string("Exact contractions on the two declared field columns, which may be dependent; complete endpoint operators and original carriers retained separately"));
            binding.set("common_carrier", std::move(common));
        }
        binding.set("old_moment", scalar(publication, old_value)); binding.set("new_moment", scalar(publication, new_value));
        publication.bind(std::move(binding));
        lease.check_cancelled();
        return publication.finish(std::move(view));
    }
    Json response_query(const Json &query, bool validate, std::optional<NativeRecord> *native = nullptr) {
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
        if (validate) {
            (void)research_detail::observed_response_model(*prepared.study, query, *grade, fields,
                meanings, coordinates, true, resources, nullptr);
            return Json::object();
        }
        scientific_detail::ResultFields publication(owner->native_publication());
        scientific_detail::ResultFieldContext binding(publication, owner->native_publication(),
            prepared.to_ambient, std::pair<std::size_t, std::uint64_t>{grade->value, *output});
        auto view = research_detail::observed_response_model(*prepared.study, query, *grade, fields,
            meanings, coordinates, false, resources, &binding);
        const auto record = publication.finish(std::move(view));
        if (native) { *native = record; return Json::object(); }
        return scientific_detail::result_fields_view(record);
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
        return {*s.native_publication, op, s.operator_pins.at(accession), std::move(fields), declarations};
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
        return {*old_source.native_publication, *new_source.native_publication, old, current, new_source.operator_pins.at(metric_id), std::move(causes),
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
        return {*s.native_publication, fields.front(), std::move(fields), std::move(metrics), declarations};
    }
    scientific_detail::DeclaredResponseBinding response_binding(const Json &query) {
        auto declarations = Json::object();
        struct EndpointPins {
            NativeRecord source;
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
            auto carriers = Json::array();
            for (const auto &c : u.endpoint()->carriers) {
                auto item = Json::object(), coordinates = Json::array();
                item.set("key", Json::string(std::to_string(c.key)));
                item.set("kind", Json::integer(static_cast<std::int64_t>(c.kind)));
                item.set("coordinate_system", Json::string(std::to_string(c.coordinate_system)));
                if (c.kind == TemporalCarrierKind::SourceGrade) {
                    item.set("grade", Json::integer(c.grade.value));
                    for (std::size_t i = 0; i < c.domain.logical_size(); ++i) {
                        const auto cell = s.coordinate(c.grade.value, i);
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
            return {*s.native_publication, L, I, u, A};
        };
        const auto old = bind("old"), current = bind("new");
        scientific_detail::DeclaredResponseBinding b{
            old.source, current.source, old.L, current.L, old.I, current.I, old.u, current.u, old.A, current.A, {}, Json{}};
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
    Json history_query(const Json &query, bool validate, std::optional<NativeRecord> *native_result = nullptr,
                       std::optional<rfc::runtime::ExecutionLease> requested_lease = {}) {
        const auto lease=requested_lease ? *requested_lease : source_execution_lease("scientific-native-temporal-history");
        lease.check_cancelled();
        std::vector<std::unique_ptr<PreparedTemporalStudy>> leased_preparations;
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
            PreparedTemporalStudy *prepared=nullptr;
            if (!validate) {
                if (requested_lease) {
                    std::size_t dimension=0;for(const auto &coordinates:p.to_ambient)dimension=std::max(dimension,coordinates.size());
                    const auto &policy=lease.policy().limits;
                    leased_preparations.push_back(std::make_unique<PreparedTemporalStudy>(*p.study,
                        rfc::numeric::ExactActionLimits{dimension,policy.sparse_entries,policy.scalar_work,policy.integer_bits},
                        [lease] {lease.check_cancelled();}));prepared=leased_preparations.back().get();
                } else prepared=&owner_.temporal(p);
            }
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
                                prepared});
            names.push_back(name);
            owners.push_back(&owner_);
        }
        if (validate) return research_detail::temporal_field_history(query,readings,maps,true);
        std::size_t slots=0,entries=0,bits=0;
        for(const auto &reading:readings) {
            slots=rfc::machine::checked_size_add(slots,reading.values.size());
            for(const auto &value:reading.values) {entries+=value!=0;bits=rfc::machine::checked_size_add(bits,rfc::storage::native_execution_detail::integer_bits(value));}
            for(const auto &[position,values]:reading.diagonal_fields) {
                (void)position;slots=rfc::machine::checked_size_add(slots,values.size());
                for(const auto &value:values) {entries+=value!=0;bits=rfc::machine::checked_size_add(bits,rfc::storage::native_execution_detail::integer_bits(value));}
            }
        }
        rfc::runtime::admit_execution({readings.size(),slots,entries,query.at("word").size(),bits},lease.policy());
        auto reservation=cpu().reserve(1);
        const auto admit_matrix=[&](const rfc::numeric::ExactMatrix &m) {research_detail::admit_structure_matrix(m,lease);};
        std::vector<SourceSnapshotPublication> snapshots;
        for(std::size_t t=0;t<readings.size();++t) {
            auto &context=owners[t]->context(names[t]);
            snapshots.push_back(source_snapshot_publication(*context.study,owners[t]->owner->native_publication(),context.to_ambient,admit_matrix));
        }
        ResultFields fields(snapshots.back().publication);
        research_detail::StructuralResultFields admission{fields,{},{},lease};
        research_detail::TemporalResultFields binding{fields,{},{},{},lease};
        auto endpoints=ResultValue::array();
        for(std::size_t t=0;t<readings.size();++t) {
            auto &endpoint_owner=*owners[t];auto &context=endpoint_owner.context(names[t]);const auto &reading=readings[t];
            std::vector<std::uint64_t> axes;
            for(const auto &coordinates:context.to_ambient)axes.push_back(fields.axis(coordinates.size()));
            binding.axes.push_back(axes);
            auto endpoint=ResultValue::object();endpoint.set("context_source",fields.publication(snapshots[t].publication));
            endpoint.set("original_study_publication",fields.publication(endpoint_owner.owner->native_publication()));
            auto grades=ResultValue::array();
            for(std::size_t k=0;k<axes.size();++k) {
                auto g=ResultValue::object();g.set("grade",ResultValue::count(k));g.set("coordinate_carrier",ResultValue::count(axes[k]));
                g.set("context_source_carrier",ResultValue::count(snapshots[t].grade_axes[k]));
                const auto &pin=endpoint_owner.source.grade_pins.at(k);
                g.set("original_source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));g.set("original_source_cutoff",ResultValue::count(pin.sequence()));
                std::vector<std::size_t> lift(context.to_ambient[k].begin(),context.to_ambient[k].end());g.set("original_coordinates",fields.indices(lift,axes[k]));grades.push(std::move(g));
            }
            endpoint.set("grades",std::move(grades));
            const auto retain_input=[&](const std::string &id,const Vector &values,std::uint16_t k) {
                const auto &source_field=endpoint_owner.tensor(id);const auto &pin=endpoint_owner.source.field_pins.at(source_field.key);
                auto input=ResultValue::object();input.set("id",ResultValue::string(id));input.set("source_namespace",ResultValue::string(std::to_string(pin.source_namespace())));
                input.set("source_cutoff",ResultValue::count(pin.sequence()));input.set("source_tensor",ResultValue::count(source_field.key));input.set("source_carrier",ResultValue::count(*source_field.output));
                input.set("coordinate_carrier",ResultValue::count(axes[k]));const auto original_axis=fields.axis(source_field.tensor.domain().logical_size());input.set("original_coordinate_carrier",ResultValue::count(original_axis));
                rfc::numeric::ExactMatrix lift(values.size(),source_field.tensor.domain().logical_size());
                for(std::size_t i=0;i<values.size();++i)lift.set(i,context.to_ambient.at(k).at(i),Rational{1});
                admission.admit(lift);admission.admit(values);input.set("realization",ResultFields::stage_matrix(lift,original_axis,axes[k]));
                Vector reported_values;std::vector<std::size_t> reported_rows;
                for(std::size_t i=0;i<values.size();++i) {
                    const auto original=context.to_ambient.at(k).at(i);
                    if(source_field.tensor.support().contains(original)) {reported_rows.push_back(i);reported_values.push_back(*source_field.tensor.state(original));}
                }
                admission.admit(reported_values);input.set("reported_values",ResultFields::stage_selected(reported_values,axes[k],values.size(),reported_rows));
                input.set("execution_values",ResultFields::stage_vector(values,axes[k]));input.set("meaning",ResultValue(endpoint_owner.source.fields.at(id).declaration));return input;
            };
            const auto endpoint_query=query.at("endpoints").at(t);endpoint.set("input",retain_input(endpoint_query.at("field").text(),reading.values,reading.grade.value));
            auto diagonal=ResultValue::array();auto current=reading.grade.value;
            for(std::size_t position=0;position<query.at("word").size();++position) {
                const auto factor=query.at("word").at(position);const auto kind=factor.at("kind").text();
                if(kind=="diagonal_field") {auto input=retain_input(factor.at("fields").at(t).text(),reading.diagonal_fields.at(position),current);input.set("word_position",ResultValue::count(position));diagonal.push(std::move(input));}
                else if(kind=="boundary")--current;else if(kind=="adjoint")++current;
            }
            endpoint.set("diagonal_inputs",std::move(diagonal));endpoints.push(endpoint);
            auto declaration=ResultValue::object();
            for(const auto &key:{"context","field","publication","extension"})if(endpoint_query.has(key))declaration.set(key,ResultValue(endpoint_query.at(key)));
            for(const auto &key:{"publication_binding","time","source_grades","metrics","field_semantics","diagonal_field_semantics","lineage"})
                if(reading.declaration.has(key))declaration.set(key,ResultValue(reading.declaration.at(key)));
            auto observations=ResultValue::array();const auto &source_field=endpoint_owner.tensor(endpoint_query.at("field").text());
            for(std::size_t i=0;i<reading.values.size();++i) {
                const auto original=context.to_ambient.at(reading.grade.value).at(i);const auto reported=source_field.tensor.support().contains(original);
                auto observation=ResultValue::object();observation.set("entity",ResultValue(reading.declaration.at("observations").at(i).at("entity")));
                observation.set("ambient_coordinate",ResultValue::count(original));observation.set("reported",ResultValue::boolean(reported));
                if(reported) {admission.admit(*source_field.tensor.state(original));auto value=fields.stage_scalar(*source_field.tensor.state(original));value.set("display",ResultValue::string("rational_object"));observation.set("observed_value",std::move(value));}
                else observation.set("observed_value",ResultValue{});
                admission.admit(reading.values[i]);auto value=fields.stage_scalar(reading.values[i]);value.set("display",ResultValue::string("rational_object"));observation.set("execution_value",std::move(value));observations.push(std::move(observation));
            }
            declaration.set("observations",std::move(observations));
            auto occurrences=ResultValue::array();
            for(std::size_t k=1;k<context.source->grade_count();++k) {
                const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};const auto &arena=context.source->grade(grade);
                for(std::size_t i=0;i<arena.incidence_count();++i) {
                    auto occurrence=ResultValue(endpoint_owner.source.occurrence_record({grade,arena.occurrence_identities[i]}));
                    occurrence.set("parent_key",ResultValue::string(std::to_string(arena.relation_identities[arena.upper_index[i]].value)));
                    occurrence.set("lower_key",ResultValue::string(std::to_string(stable_cell_key(*context.source,{rfc::core::Grade{static_cast<std::uint16_t>(k-1)},arena.lower_index[i]}).key)));
                    const auto coefficient=context.source->incidence_coefficient(grade,i);admission.admit(coefficient);occurrence.set("coefficient",fields.stage_scalar(coefficient));
                    occurrence.set("lattice_coefficient",ResultValue(NativeRecordValue(arena.lattice_coefficient[i])));
                    const auto parent=arena.cell_scales[arena.upper_index[i]],lower=context.source->grade(rfc::core::Grade{static_cast<std::uint16_t>(k-1)}).cell_scales[arena.lower_index[i]];
                    admission.admit(parent);admission.admit(lower);auto parent_value=fields.stage_scalar(parent),lower_value=fields.stage_scalar(lower);
                    parent_value.set("display",ResultValue::string("rational_object"));lower_value.set("display",ResultValue::string("rational_object"));occurrence.set("parent_scale",std::move(parent_value));occurrence.set("lower_scale",std::move(lower_value));
                    auto aspects=ResultValue::array();for(std::size_t a=0;a<context.source->signature().aspects().size();++a)aspects.push(ResultValue::boolean(arena.aspects.get(static_cast<rfc::composite::AspectId>(a),i)));
                    occurrence.set("aspects",std::move(aspects));occurrence.set("evidence",ResultValue(reading.declaration.at("source_occurrences").at(occurrences.size()).at("evidence")));occurrences.push(std::move(occurrence));
                }
            }
            declaration.set("source_occurrences",std::move(occurrences));binding.declarations.push_back(std::move(declaration));
            binding.coordinates.push_back(ResultValue(reading.coordinates));
        }
        auto owner_binding=ResultValue::object();owner_binding.set("kind",ResultValue::string("actual_temporal_history_endpoints"));owner_binding.set("endpoints",std::move(endpoints));fields.bind(std::move(owner_binding));
        const auto result=research_detail::temporal_history_fields(query,readings,maps,binding);lease.check_cancelled();const auto record=fields.finish(result);
        if(native_result) {*native_result=record;return Json::object();}return result_fields_view(record);
    }
    Json censor_query(const Json &query, bool validate, std::optional<NativeRecord> *native = nullptr) {
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
        if (validate) {
            (void)research_detail::censored_endpoint_bounds(
                *p.study, query, fields, upper_support, meanings, temporal_coordinates(p, 0), true);
            return Json::object();
        }
        auto reservation = cpu().reserve(1);
        const auto lease = source_execution_lease("scientific-censored-endpoint-family");
        lease.check_cancelled();
        const auto snapshot = scientific_detail::source_snapshot_publication(
            *p.study, owner->native_publication(), p.to_ambient);
        scientific_detail::ResultFields publication(snapshot.publication);
        using scientific_detail::ResultValue;
        using scientific_detail::ResultFields;
        const auto group_axis = publication.axis(1 + query.at("strata").size());
        const auto contrast_axis = publication.axis(query.has("contrasts") ? query.at("contrasts").size() : 0);
        const auto input_axis = publication.axis(tensor(*ids.begin()).tensor.domain().logical_size());
        auto binding = ResultValue::object(), inputs = ResultValue::array();
        binding.set("kind", ResultValue::string("actual_censored_endpoint_inputs"));
        binding.set("original_study_publication", ResultValue(snapshot.publication.at("bindings").at(0)
            .at("original_study_publication").unbound_value()));
        binding.set("input_carrier", ResultValue::count(*carrier));
        binding.set("input_coordinate_carrier", ResultValue::count(input_axis));
        binding.set("source_carrier", ResultValue::count(snapshot.grade_axes.at(0)));
        rfc::numeric::ExactMatrix realization(p.to_ambient.at(0).size(), tensor(*ids.begin()).tensor.domain().logical_size());
        for (std::size_t i = 0; i < p.to_ambient.at(0).size(); ++i)
            realization.set(i, p.to_ambient.at(0).at(i), Rational{1});
        binding.set("input_realization", ResultFields::stage_matrix(realization, input_axis, snapshot.grade_axes.at(0)));
        for (const auto &id : ids) {
            const auto &field = tensor(id);
            const auto &pin = source.field_pins.at(field.key);
            auto input = ResultValue::object();
            input.set("id", ResultValue::string(id));
            input.set("tensor", ResultValue::count(field.key));
            input.set("source_namespace", ResultValue::string(std::to_string(pin.source_namespace())));
            input.set("cutoff", ResultValue::count(pin.sequence()));
            input.set("carrier", ResultValue::count(*field.output));
            std::vector<std::size_t> supported;
            Vector values;
            for (std::size_t i = 0; i < p.to_ambient.at(0).size(); ++i) {
                const auto coordinate = p.to_ambient.at(0).at(i);
                if (field.tensor.support().contains(coordinate)) {
                    supported.push_back(i); values.push_back(*field.tensor.state(coordinate));
                }
            }
            input.set("context_values", ResultFields::stage_selected(values, snapshot.grade_axes.at(0),
                p.to_ambient.at(0).size(), supported));
            input.set("meaning", ResultValue(source.fields.at(id).declaration));
            inputs.push(std::move(input));
        }
        binding.set("inputs", std::move(inputs));
        auto groups = ResultValue::array(), contrasts = ResultValue::array();
        groups.push(ResultValue::string("selected_population"));
        for (std::size_t i = 0; i < query.at("strata").size(); ++i) groups.push(ResultValue(query.at("strata").at(i)));
        if (query.has("contrasts"))
            for (std::size_t i = 0; i < query.at("contrasts").size(); ++i) contrasts.push(ResultValue(query.at("contrasts").at(i)));
        binding.set("group_carrier", ResultValue::count(group_axis)); binding.set("groups", std::move(groups));
        binding.set("contrast_carrier", ResultValue::count(contrast_axis)); binding.set("contrasts", std::move(contrasts));
        const auto time_axis = publication.axis(query.at("times").size());
        Vector times;
        for (std::size_t i = 0; i < query.at("times").size(); ++i) times.push_back(exact_rational(query.at("times").at(i).text()));
        binding.set("time_carrier", ResultValue::count(time_axis));
        binding.set("evaluation_times", ResultFields::stage_vector(times, time_axis));
        binding.set("origin", publication.stage_scalar(exact_rational(query.at("origin").text())));
        binding.set("execution_scope", ResultValue::string("Canonical family actions; CPU reservation and cancellation before and after execution"));
        publication.bind(std::move(binding));
        research_detail::CensoredResultFields axes{publication, snapshot.grade_axes.at(0), group_axis, contrast_axis};
        auto result = research_detail::censored_endpoint_bounds(
            *p.study, query, fields, upper_support, meanings, temporal_coordinates(p, 0), false, &axes);
        lease.check_cancelled();
        const auto record = publication.finish(std::move(result));
        if (native) { *native = record; return Json::object(); }
        return scientific_detail::result_fields_view(record);
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
            validate_chain_values(*context(query.at("context").text()).study, query, chain_operand(query));
            return;
        }
        if (operation == "realized_path_comparison") {
            (void)compare_realized_paths(*context(query.at("context").text()).study, query, true);
            return;
        }
        if (operation == "local_sections") {
            (void)local_section_fields(query, true);
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
            const auto [x, grade] = input(query.at("context").text(), query.at("field").text(), false);
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
            const auto binding = std::to_string(grade) + "|" +
                std::to_string(*tensor(query.at("field").text()).output) + "|" + declaration.at("cells").dump();
            require(!green_spaces.contains(coordinate_id) ||
                        green_spaces.at(coordinate_id).binding == binding,
                    "Green coordinate identity reused with a different basis");
            if (!green_spaces.contains(coordinate_id)) {
                const auto bytes = rfc::machine::checked_size_add(binding.size(), sizeof(GreenCoordinates));
                auto charge = cpu().coordinator().reserve(1, bytes);
                cpu().coordinator().resize(charge, 0, bytes);
                green_spaces.emplace(coordinate_id, GreenCoordinates{binding, std::move(charge)});
            }
            require(exact_decimal(query.at("lambda").text()) > 0,
                    "positive Green parameter required");
            std::set<std::size_t> indices;
            require(query.at("interior").size() <= x.size(), "interior exceeds the actual coordinate family");
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
            require(exact_rational(query.at("lambda").text()) > 0,
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
                    {"correspondence", "reverse_correspondence", "map_output"});
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
            (void)research_detail::NativeFactoredTemporalMaps::requested(query);
            if (query.has("reverse_correspondence")) {
                const auto forward = query.has("correspondence")
                    ? temporal_map(name, other, k, query.at("correspondence"))
                    : rfc::numeric::ExactMatrix::identity(prepared.to_ambient.at(k).size());
                const auto reverse = temporal_map(other, name, k, query.at("reverse_correspondence"));
                require(reverse.multiply(forward) == rfc::numeric::ExactMatrix::identity(forward.cols()) &&
                        forward.multiply(reverse) == rfc::numeric::ExactMatrix::identity(forward.rows()),
                        "moment common carrier requires two sided inverse correspondences");
            }
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
    NativeRecord source_fields(const Json &query,
                              std::optional<rfc::runtime::ExecutionLease> supplied_lease = {}) const {
        using R = ResultValue;
        auto lease = supplied_lease.value_or(source_execution_lease("scientific-native-source"));
        lease.check_cancelled();
        rfc::runtime::admit_execution({1, 0, 0, 0, 0}, lease.policy());
        ResultFields fields(owner->native_publication());
        const auto endpoint = source.directory->pin_snapshot();
        auto binding = R::object();
        binding.set("kind", R::string("original_source"));
        binding.set("source_namespace", R::string(std::to_string(source.grade_pins.at(0).source_namespace())));
        binding.set("cutoff", R::count(endpoint->sequence));
        const auto operation = query.at("operation").text();
        binding.set("operation", R::string(operation));
        auto view = R::object();
        if (operation == "inspect") {
            auto records = rfc::machine::checked_size_add(source.fields.size(), endpoint->carriers.size());
            records = rfc::machine::checked_size_add(records, source.cells.size());
            records = rfc::machine::checked_size_add(records, source.documents.size());
            records = rfc::machine::checked_size_add(records, source.operators.size());
            records = rfc::machine::checked_size_add(records, source.contexts.size());
            rfc::runtime::admit_execution({rfc::machine::checked_size_add(records, 1), 0, 0, records, 0}, lease.policy());
            auto sources = owner->native_publication().raw_sources();
            auto value = source.catalog_value(sources);
            fields.retain_sources(std::move(sources));
            view = R(std::move(value));
        } else if (operation == "coverage") {
            const auto id = query.at("field").text();
            const auto &t = tensor(id);
            const auto pin = source.field_pins.at(t.key);
            binding.set("source_namespace", R::string(std::to_string(pin.source_namespace())));
            binding.set("cutoff", R::count(pin.sequence()));
            binding.set("tensor", R::count(t.key));
            binding.set("carrier", R::count(*t.output));
            binding.set("field", R::string(id));
            binding.set("input", fields.stored(t.key, "vector"));
            const auto size = t.tensor.domain().logical_size();
            rfc::runtime::admit_execution({rfc::machine::checked_size_add(size, 1), 0, 0, size, 0}, lease.policy());
            view.set("coordinates", R::count(size));
            view.set("reported", R::count(t.tensor.support().count()));
            view.set("missing", R::count(size - t.tensor.support().count()));
            view.set("native_tensor", R::string(std::to_string(t.key)));
            view.set("native_carrier", R::string(std::to_string(*t.output)));
            auto extents = R::array();
            for (const auto extent : t.tensor.domain().extents()) extents.push(R::count(extent));
            view.set("extents", std::move(extents));
            const auto carriers = temporal_endpoint_detail::carriers(*endpoint);
            const auto &carrier = *carriers.at(*t.output);
            auto observations = R::array();
            for (std::size_t i = 0; i < size; ++i) {
                lease.check_cancelled();
                auto row = R::object();
                row.set("coordinate", R::count(i));
                if (carrier.kind == TemporalCarrierKind::SourceGrade)
                    row.set("entity", R::string(source.coordinate_id(carrier.grade.value, i)));
                row.set("reported", R::boolean(t.tensor.support().contains(i)));
                auto value = fields.stored(t.key, "rational_object");
                value.set("index", R::count(i));
                row.set("value", std::move(value));
                observations.push(std::move(row));
            }
            view.set("observations", std::move(observations));
        } else {
            const auto cell = source.resolve_cell(query.at("entity").text());
            const auto entity = source.cell_record(cell);
            view.set("entity", R(entity));
            view.set("evidence", R(source.evidence(entity.at("evidence"))));
            binding.set("cell_grade", R::count(cell.grade.value));
            binding.set("cell_key", R::count(cell.key));
            auto occurrences = R::array(), identities = R::array();
            std::size_t selected{}, scanned{};
            for (std::size_t k = 1; k < endpoint->source.grade_count(); ++k) {
                const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};
                const auto &arena = endpoint->source.grade(grade);
                for (std::size_t slot = 0; slot < arena.incidence_count(); ++slot) {
                    lease.check_cancelled();
                    scanned = rfc::machine::checked_size_add(scanned, 1);
                    rfc::runtime::admit_execution({1, 0, 0, scanned, 0}, lease.policy());
                    const auto upper = stable_cell_key(endpoint->source, {grade, arena.upper_index[slot]});
                    const auto lower = stable_cell_key(endpoint->source,
                        {{static_cast<std::uint16_t>(k - 1)}, arena.lower_index[slot]});
                    if (cell != upper && cell != lower) continue;
                    selected = rfc::machine::checked_size_add(selected, 1);
                    const auto values = rfc::machine::checked_size_mul(selected, 2);
                    rfc::runtime::admit_execution({rfc::machine::checked_size_add(selected, 1), values, values, scanned, 0}, lease.policy());
                    const StableOccurrenceKey key{grade, arena.occurrence_identities[slot]};
                    auto item = R::object(), identity = R::object();
                    item.set("key", R::string(std::to_string(key.key.value)));
                    item.set("grade", R::count(k));
                    item.set("upper", R::string(source.coordinate_id(k, arena.upper_index[slot])));
                    item.set("lower", R::string(source.coordinate_id(k - 1, arena.lower_index[slot])));
                    auto coefficient = fields.scalar(endpoint->source.incidence_coefficient(grade, slot));
                    coefficient.set("display", R::string("rational_object"));
                    item.set("coefficient", std::move(coefficient));
                    auto magnitude = fields.scalar(abs(endpoint->source.incidence_coefficient(grade, slot)));
                    magnitude.set("display", R::string("rational_object"));
                    item.set("magnitude", std::move(magnitude));
                    auto aspects = R::object();
                    for (std::size_t a = 0; a < endpoint->source.signature().aspects().size(); ++a)
                        aspects.set(endpoint->source.signature().aspects().at(
                            static_cast<rfc::composite::AspectId>(a)).name,
                            R::boolean(arena.aspects.get(static_cast<rfc::composite::AspectId>(a), slot)));
                    item.set("aspects", std::move(aspects));
                    item.set("evidence", R(source.evidence(source.occurrence_record(key).at("evidence"))));
                    identity.set("grade", R::count(k)); identity.set("key", R::count(key.key.value));
                    identity.set("slot", R::count(slot)); identity.set("upper", R::count(upper.key));
                    identity.set("lower", R::count(lower.key));
                    identities.push(std::move(identity)); occurrences.push(std::move(item));
                }
            }
            binding.set("occurrences", std::move(identities));
            view.set("occurrences", std::move(occurrences));
        }
        fields.bind(std::move(binding));
        lease.check_cancelled();
        return fields.finish(std::move(view));
    }
    Vector chain_operand(const Json &query) {
        members(query, {"operation", "context", "grade", "coefficient_domain"}, {"field", "values"});
        require(query.has("field") != query.has("values"), "chain requires values or a source field");
        auto &prepared = context(text(query.at("context")));
        const auto grade = count(query.at("grade"), prepared.study->source().grade_count() - 1, "chain grade");
        if (query.has("values")) return declared_chain_values(query, prepared.study->source().grade(
            rfc::core::Grade{static_cast<std::uint16_t>(grade)}).cell_count());
        require(text(query.at("coefficient_domain")) == "rational", "source fields have rational coordinates");
        const auto [values, field_grade] = input(text(query.at("context")), text(query.at("field")));
        require(field_grade == grade, "chain field grade mismatch");
        return values;
    }
    NativeRecord structural_fields(const Json &query,
        std::optional<rfc::runtime::ExecutionLease> native_lease = {}) {
        const auto name = text(query.at("context"));
        auto &prepared = context(name);
        const auto lease = native_lease.value_or(source_execution_lease("scientific-structure"));
        lease.check_cancelled();
        const auto snapshot = source_snapshot_publication(*prepared.study,
            owner->native_publication(), prepared.to_ambient);
        ResultFields publication(snapshot.publication);
        const auto grade = count(query.at("grade"), snapshot.grade_axes.size() - 1, "structural grade");
        const bool chain = text(query.at("operation")) == "chain_certificate";
        auto result = chain ? chain_fields(*prepared.study, query, chain_operand(query),
            publication, snapshot.grade_axes, lease) : path_fields(*prepared.study, query,
            publication, snapshot.grade_axes, lease);
        auto decorated = ResultValue::array();
        const auto original_coordinates = result.at("coordinates");
        for (std::size_t i = 0; i < original_coordinates.size(); ++i) {
            auto coordinate = original_coordinates.at(i);
            const auto cell = source.coordinate(grade, prepared.to_ambient.at(grade).at(i));
            coordinate.set("entity", ResultValue(source.cell_record(cell).at("id")));
            coordinate.set("evidence", ResultValue(source.cell_record(cell).at("evidence")));
            decorated.push(std::move(coordinate));
        }
        result.set("coordinates", std::move(decorated));
        if (chain) {
            auto occurrences = ResultValue::array();
            const auto original = result.at("occurrences");
            for (std::size_t i = 0; i < original.size(); ++i) {
                auto occurrence = original.at(i);
                const StableOccurrenceKey cell{rfc::core::Grade{static_cast<std::uint16_t>(grade)},
                    rfc::core::OccurrenceKey{key(Json::string(occurrence.at("key").text()))}};
                occurrence.set("evidence", ResultValue(source.occurrence_record(cell).at("evidence")));
                occurrences.push(std::move(occurrence));
            }
            result.set("occurrences", std::move(occurrences));
        }
        auto binding = ResultValue::object();
        binding.set("kind", ResultValue::string(chain ? "declared_chain_certificate" : "declared_path_family"));
        binding.set("context", ResultValue::string(name));
        binding.set("grade", ResultValue::count(grade));
        binding.set("original_study_publication", ResultValue(snapshot.publication.at("bindings").at(0)
            .at("original_study_publication").unbound_value()));
        if (chain && query.has("field")) {
            const auto &actual = tensor(text(query.at("field")));
            const auto &pin = source.field_pins.at(actual.key);
            binding.set("source_namespace", ResultValue::string(std::to_string(pin.source_namespace())));
            binding.set("source_cutoff", ResultValue::count(pin.sequence()));
            binding.set("source_tensor", ResultValue::count(actual.key));
            binding.set("source_carrier", ResultValue::count(*actual.output));
            const auto ambient = publication.axis(actual.tensor.domain().logical_size());
            binding.set("source_field", ResultFields::staged(actual.tensor, {}, ambient, "vector"));
            rfc::numeric::ExactMatrix lift(prepared.to_ambient.at(grade).size(), actual.tensor.domain().logical_size());
            for (std::size_t i = 0; i < lift.rows(); ++i) lift.set(i, prepared.to_ambient.at(grade).at(i), 1);
            binding.set("context_realization", ResultFields::stage_matrix(lift, ambient, snapshot.grade_axes.at(grade)));
        }
        binding.set("result", result);
        publication.bind(std::move(binding));
        lease.check_cancelled();
        return publication.finish(std::move(result));
    }
    NativeRecord green_response_fields(const Json &query,
        std::optional<rfc::runtime::ExecutionLease> native_lease = {}) {
        const auto name = query.at("context").text();
        const auto [x, grade] = input(name, query.at("field").text(), false);
        auto &prepared = context(name);
        const auto lease = native_lease.value_or(source_execution_lease("scientific-partitioned-green"));
        lease.check_cancelled();
        (void)cpu().readmit();
        const auto &policy = lease.policy();
        rfc::runtime::admit_execution({1, x.size(), 0, 0, 0}, policy);
        const rfc::numeric::ExactActionLimits action_limits{x.size(), policy.limits.sparse_entries,
            policy.limits.scalar_work, policy.limits.integer_bits};
        const auto checkpoint = [&] {lease.check_cancelled();};
        const auto lambda = exact_decimal(query.at("lambda").text());
        auto operation = query;
        operation = Json::parse(operation.dump());
        operation.set("field", Json::string(""));
        const auto cache_key =
            operation.dump() + "|grade:" + std::to_string(grade) +
            "|carrier:" + std::to_string(*tensor(query.at("field").text()).output);
        const bool reused = green_solves.contains(cache_key);
        if (!reused) {
            // Charge the complete named square coefficient containers before
            // preparation. This envelope excludes GMP limbs and process RSS.
            const auto square = rfc::machine::checked_size_mul(x.size(), x.size());
            const auto slots = rfc::machine::checked_size_add(
                rfc::machine::checked_size_mul(square, 8),
                rfc::machine::checked_size_mul(x.size(), 8));
            const auto slot_bytes = sizeof(std::pair<const std::pair<std::size_t, std::size_t>, Rational>);
            const auto envelope = rfc::machine::checked_size_mul(slots, slot_bytes);
            rfc::runtime::StagedResourceLease stage(cpu().coordinator(),
                {"partitioned-green-preparation", 1, 0, 0, envelope, 0}, lease);
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
            rfc::local::PreparedInteriorInterfaceSolve solve(std::move(declaration),
                std::move(a), std::move(interior), action_limits, checkpoint);
            const auto retained_bytes = rfc::machine::checked_size_add(
                rfc::machine::checked_size_mul(solve.retained_entries(), slot_bytes),
                rfc::machine::checked_size_add(sizeof(solve), cache_key.size()));
            auto retained = stage.finish(std::move(solve), retained_bytes);
            green_solves.emplace(cache_key, std::move(retained));
        }
        const auto &solve = *green_solves.at(cache_key);
        auto execution = cpu().coordinator().reserve(1,
            rfc::machine::checked_size_mul(x.size(), sizeof(Rational) * 8));
        const auto y = solve.apply(x, action_limits, checkpoint);
        check_values(y);
        lease.check_cancelled();
        const auto snapshot = source_snapshot_publication(*prepared.study,
            owner->native_publication(), prepared.to_ambient);
        ResultFields publication(snapshot.publication);
        const auto axis = snapshot.grade_axes.at(grade);
        const auto output = publication.vector(y, axis);
        auto view = ResultValue::object(), field = ResultValue::array();
        for (std::size_t i = 0; i < y.size(); ++i) {
            const auto cell = source.coordinate(grade, prepared.to_ambient.at(grade).at(i));
            auto row = ResultValue::object();
            row.set("coordinate", ResultValue::count(i));
            row.set("native_key", ResultValue::string(std::to_string(cell.key)));
            row.set("entity", ResultValue(source.cell_record(cell).at("id")));
            auto value = output;
            value.set("display", ResultValue::string("rational_object"));
            value.set("index", ResultValue::count(i));
            row.set("value", std::move(value)); field.push(std::move(row));
        }
        view.set("field", std::move(field));
        view.set("interior", ResultValue(query.at("interior")));
        auto interface = ResultValue::array();
        for (const auto i : solve.interface()) interface.push(ResultValue::count(i));
        view.set("interface_coordinates", std::move(interface));
        view.set("residual_zero", ResultValue::boolean(true));
        auto binding = ResultValue::object();
        binding.set("kind", ResultValue::string("declared_partitioned_green_action"));
        binding.set("original_study_publication", ResultValue(snapshot.publication.at("bindings").at(0)
            .at("original_study_publication").unbound_value()));
        const auto &input_tensor = tensor(query.at("field").text());
        const auto &pin = source.field_pins.at(input_tensor.key);
        binding.set("input_source_namespace", ResultValue::string(std::to_string(pin.source_namespace())));
        binding.set("input_cutoff", ResultValue::count(pin.sequence()));
        binding.set("input_tensor", ResultValue::count(input_tensor.key));
        binding.set("input_carrier", ResultValue::count(*input_tensor.output));
        binding.set("context", ResultValue::string(name));
        binding.set("grade", ResultValue::count(grade));
        binding.set("coordinates", ResultValue(query.at("coordinates")));
        binding.set("lambda", publication.stage_scalar(lambda));
        const auto ambient = publication.axis(input_tensor.tensor.domain().logical_size());
        rfc::numeric::ExactMatrix realization(x.size(), input_tensor.tensor.domain().logical_size());
        for (std::size_t i = 0; i < x.size(); ++i)
            realization.set(i, prepared.to_ambient.at(grade).at(i), Rational{1});
        binding.set("input_realization", ResultFields::stage_matrix(realization, ambient, axis));
        binding.set("context_input", ResultFields::stage_vector(x, axis));
        binding.set("output", output);
        const auto interior_axis = publication.axis(solve.interior().size());
        const auto interface_axis = publication.axis(solve.interface().size());
        rfc::numeric::ExactMatrix interior_lift(x.size(), solve.interior().size());
        rfc::numeric::ExactMatrix interface_lift(x.size(), solve.interface().size());
        for (std::size_t i = 0; i < solve.interior().size(); ++i) interior_lift.set(solve.interior()[i], i, 1);
        for (std::size_t i = 0; i < solve.interface().size(); ++i) interface_lift.set(solve.interface()[i], i, 1);
        binding.set("interior_coordinates", publication.indices(solve.interior(), interior_axis));
        binding.set("interface_coordinates", publication.indices(solve.interface(), interface_axis));
        binding.set("interior_lift", ResultFields::stage_matrix(interior_lift, interior_axis, axis));
        binding.set("interface_lift", ResultFields::stage_matrix(interface_lift, interface_axis, axis));
        binding.set("schur", ResultFields::stage_matrix(solve.schur(), interface_axis, interface_axis));
        const auto residual = rfc::numeric::sub(rfc::numeric::add(y, rfc::numeric::scale(lambda,
            prepared.study->hodge(rfc::core::Grade{grade}).apply(y))), x);
        require(rfc::numeric::is_zero(residual), "native partitioned Green reconstruction residual");
        binding.set("residual", ResultFields::stage_vector(residual, axis));
        publication.bind(std::move(binding));
        const auto native = publication.finish(std::move(view));
        auto root = native.materialize_value();
        auto execution_view = ResultValue::object();
        execution_view.set("factorization_reused", ResultValue::boolean(reused));
        execution_view.set("retained_entries", ResultValue::count(solve.retained_entries()));
        execution_view.set("retained_preparations", ResultValue::count(green_solves.size()));
        execution_view.set("declared_coordinate_families", ResultValue::count(green_spaces.size()));
        std::get<NativeRecordValue::Object>(root.data).members.emplace_back("view_execution", std::move(execution_view).release());
        lease.check_cancelled();
        return native.derive(std::move(root));
    }
    NativeRecord native_boundary_fields(const Json &query) {
        const auto name = query.at("field_set").text();
        const auto analysis = native_analysis(name);
        std::shared_ptr<const Vector> values_pin;
        std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> serial;
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
            values_pin = derived.values_pin();
            require(derived.source_publication().get() == analysis.source_publication().get(),
                    "native boundary result lost its source owner");
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
            values_pin = derived.values_pin();
            require(derived.source_publication().get() == analysis.source_publication().get(),
                    "native boundary result lost its source owner");
            report = derived.provenance();
            execution.set("requested_workers",
                          cardinality(derived.execution().requested_workers));
            execution.set("worker_slots", cardinality(derived.execution().worker_slots));
            execution.set("launched_threads",
                          cardinality(derived.execution().launched_threads));
            execution.set("tiles", cardinality(derived.execution().tiles));
            execution.set("completed_tiles", cardinality(derived.execution().completed_tiles));
        }
        const auto &values = *values_pin;
        check_values(values);
        ResultFields publication(owner->native_publication());
        auto view = ResultValue::object();
        view.set("schema", ResultValue::string("gadget.scientific.native.boundary.result.v1"));
        auto action = ResultValue::object();
        action.set("kind", ResultValue::string("actual_native_cell_boundary_action"));
        action.set("source_publication", publication.publication(owner->native_publication()));
        auto inputs = ResultValue::array();
        for (const auto &binding : analysis.binding().inputs()) {
            auto item = ResultValue::object();
            item.set("cell_grade", ResultValue::count(binding.cell.grade.value));
            item.set("cell_key", ResultValue::count(binding.cell.key));
            item.set("carrier", ResultValue::count(binding.carrier));
            item.set("coordinate_system", ResultValue::count(binding.coordinates));
            item.set("offset", ResultValue::count(binding.offset));
            item.set("dimension", ResultValue::count(binding.dimension));
            item.set("field", publication.stored(binding.field, "vector"));
            item.set("source_namespace", ResultValue::string(std::to_string(source.field_pins.at(binding.field).source_namespace())));
            item.set("cutoff", ResultValue::count(source.field_pins.at(binding.field).sequence()));
            inputs.push(std::move(item));
        }
        action.set("inputs", std::move(inputs));
        auto outputs = ResultValue::array();
        auto output = ResultValue::array();
        for (const auto &binding : analysis.binding().outputs()) {
            auto row = ResultValue::object();
            row.set("entity", ResultValue(source.cell_record(binding.cell).at("id")));
            row.set("carrier", ResultValue::string(std::to_string(binding.carrier)));
            row.set("coordinate_system", ResultValue::string(std::to_string(binding.coordinates)));
            const auto field = publication.cell_vector(
                Vector(values.data() + binding.offset, values.data() + binding.offset + binding.dimension),
                binding.carrier, binding.cell, "vector_objects");
            row.set("values", field);
            row.set("support", ResultValue::string("complete; supported zero and cancellation retained"));
            output.push(std::move(row));
            auto item = ResultValue::object();
            item.set("cell_grade", ResultValue::count(binding.cell.grade.value));
            item.set("cell_key", ResultValue::count(binding.cell.key));
            item.set("carrier", ResultValue::count(binding.carrier));
            item.set("coordinate_system", ResultValue::count(binding.coordinates));
            item.set("offset", ResultValue::count(binding.offset));
            item.set("dimension", ResultValue::count(binding.dimension));
            item.set("field", field);
            outputs.push(std::move(item));
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

        view.set("field_set", ResultValue(source.field_sets.at(name).declaration));
        auto retained = ResultValue::array(), occurrences = ResultValue::array();
        for (const auto &binding : analysis.binding().occurrences()) {
            auto occurrence = ResultValue::object();
            occurrence.set("grade", ResultValue::count(binding.occurrence.upper_grade.value));
            occurrence.set("key", ResultValue::string(std::to_string(binding.occurrence.key.value)));
            occurrence.set("upper", ResultValue(source.cell_record(binding.upper).at("id")));
            occurrence.set("lower", ResultValue(source.cell_record(binding.lower).at("id")));
            occurrence.set("active", ResultValue::boolean(binding.active));
            auto coefficient = publication.stage_scalar(binding.coefficient);
            coefficient.set("display", ResultValue::string("rational_object"));
            occurrence.set("coefficient", std::move(coefficient));
            occurrence.set("map", binding.map ? ResultValue::string(std::to_string(*binding.map)) : ResultValue{});
            occurrence.set("map_role", ResultValue::string(!binding.map ? "declared_constant_identity" :
                (binding.role == TemporalTensorRole::RealizedIncidenceMap ? "realized_incidence" : "coefficient_independent")));
            occurrence.set("evidence", ResultValue(source.evidence(source.occurrence_record(binding.occurrence).at("evidence"))));
            retained.push(std::move(occurrence));
            auto item = ResultValue::object();
            item.set("grade", ResultValue::count(binding.occurrence.upper_grade.value));
            item.set("key", ResultValue::count(binding.occurrence.key.value));
            item.set("upper_grade", ResultValue::count(binding.upper.grade.value));
            item.set("upper_key", ResultValue::count(binding.upper.key));
            item.set("lower_grade", ResultValue::count(binding.lower.grade.value));
            item.set("lower_key", ResultValue::count(binding.lower.key));
            item.set("input_carrier", ResultValue::count(binding.input_carrier));
            item.set("output_carrier", ResultValue::count(binding.output_carrier));
            item.set("source_offset", ResultValue::count(binding.source_offset));
            item.set("target_offset", ResultValue::count(binding.target_offset));
            item.set("active", ResultValue::boolean(binding.active));
            item.set("coefficient", publication.stage_scalar(binding.coefficient));
            item.set("map", binding.map ? publication.stored(*binding.map, "declared_map") : ResultValue{});
            item.set("role", ResultValue::count(static_cast<std::size_t>(binding.role)));
            occurrences.push(std::move(item));
        }
        action.set("outputs", std::move(outputs));
        action.set("occurrences", std::move(occurrences));
        publication.bind(std::move(action));
        view.set("retained_occurrences", std::move(retained));
        view.set(
            "mathematical_definition",
            ResultValue::string("native rational B_k field action; complete declared fibers and "
                         "occurrence maps, with retained inactive and repeated occurrences"));
        view.set("result", std::move(output));
        const auto native = publication.finish(std::move(view));
        auto root = native.materialize_value();
        std::get<NativeRecordValue::Object>(root.data).members.emplace_back("execution", import_record_view(execution));
        return native.derive(std::move(root));

    }
    NativeRecord record(const NativeRecord &query, const Json &admitted) {
        using N = NativeRecordValue;
        const auto operation = admitted.at("operation").text();
        require(has_native_scientific_publisher(operation), "scientific operation has no native publisher");
        const auto fields = [&] {
            if (operation == "inspect" || operation == "coverage" || operation == "provenance")
                return source_fields(admitted);
            if (operation == "trajectory_structure") return research_detail::trajectory_fields(*owner,admitted,resources,
                source_execution_lease("scientific-native-trajectory"));
            if (operation == "temporal_field_history") {std::optional<NativeRecord> fields;(void)history_query(admitted,false,&fields);return *fields;}
            if (operation == "structure_transition") return structure_transition_fields(admitted);
            if (operation == "chain_certificate" || operation == "realized_path_comparison") return structural_fields(admitted);
            if (operation == "green_response") return green_response_fields(admitted);
            if (operation == "native_boundary") return native_boundary_fields(admitted);
            if (operation == "support_mask") return support_fields(admitted);
            if (operation == "tensor_window") return window_fields(admitted);
            if (operation == "typed_moment" || operation == "temporal_moment_change") return moment_fields(admitted);
            if (operation == "field_summary") return field_summary_fields(admitted);
            if (operation == "boundary" || operation == "adjoint" || operation == "hodge" || operation == "green_action")
                return direct_action_fields(admitted);
            if (operation == "local_sections") return local_section_fields(admitted);
            if (operation == "censored_endpoint_bounds") {
                std::optional<NativeRecord> result;
                (void)censor_query(admitted, false, &result);
                return *result;
            }
            if (operation == "source_structure_ablation") {
                std::optional<NativeRecord> result;
                (void)ablation_query(admitted, false, &result);
                return *result;
            }
            if (operation == "structural_reference_calibration") {
                std::optional<NativeRecord> result;
                (void)structural_calibration_query(admitted, false, &result);
                return *result;
            }
            if (operation == "relational_field_ensemble") {
                std::optional<NativeRecord> result;
                (void)ensemble_query(admitted, false, &result);
                return *result;
            }
            if (operation == "observed_response_model" || operation == "constrained_response_model") {
                std::optional<NativeRecord> result;
                (void)response_query(admitted, false, &result);
                return *result;
            }
            if (operation == "paired_field_contrasts" || operation == "controlled_pairings")
                return research_detail::comparison_fields(*owner, admitted, resources,
                    source_execution_lease("scientific-declared-comparison"));
            if (operation == "linear_model_comparison")
                return research_detail::linear_model_fields(*owner, admitted, resources).publication;
            if (operation == "finite_reference_calibration")
                return research_detail::finite_reference_fields(*owner, admitted, resources);
            if (operation == "hodge_frame_decomposition")
                return native_hodge.execute_fields(admitted, hodge_binding(admitted),
                    source_execution_lease("gadget-native-hodge-machine-addressable-v1"), cpu());
            if (operation == "metric_relative_novelty")
                return native_novelty.execute_fields(admitted, novelty_binding(admitted, true),
                    source_execution_lease("gadget-native-novelty-machine-addressable-v1"), cpu());
            if (operation == "declared_response_transition") {
                auto reservation = cpu().reserve(1);
                return declared_responses.execute_fields(admitted, response_binding(admitted),
                    scientific_detail::declared_response_lease());
            }
            return native_observations.execute_fields(admitted, observation_binding(admitted),
                source_execution_lease("gadget-native-observation-machine-addressable-v1"), cpu());
        }();
        std::uint64_t source_key = 1;
        for (const auto &bytes : fields.raw_sources())
            source_key = std::max(source_key, rfc::machine::checked_size_add(bytes.key(), 1));
        if (!original_publication)
            original_publication.emplace(owner->native_publication().archive_source(1));
        auto bytes = original_publication->rekey(source_key);
        N value(N::Object{{{"schema", N(std::string("gadget.scientific.result.fields.v1"))},
            {"certificate", import_record_view(certificate(admitted))},
            {"fields", fields.unbound_value()},
            {"source_study_publication", N(N::ByteSpan{source_key, 0, bytes.bytes().size(), false})}}});
        const auto retained = fields.derive(std::move(value), {std::move(bytes)});
        return retain_operation(retained, query, retained.unbound_value());
    }
    Json run(const Json &query, std::optional<rfc::runtime::ExecutionLease> native_lease = {}) {
        auto result = certificate(query);
        const auto operation = query.at("operation").text();
        std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> serial;
        if (resources && operation != "typed_moment" && operation != "temporal_moment_change" && operation != "censored_endpoint_bounds" && operation != "field_summary" && operation != "boundary" &&
            operation != "adjoint" && operation != "hodge" && operation != "green_action" &&
            operation != "native_boundary" && operation != "support_mask" &&
            operation != "paired_field_contrasts" && operation != "controlled_pairings" &&
            operation != "trajectory_structure" &&
            operation != "linear_model_comparison" && operation != "finite_reference_calibration" &&
            operation != "relational_field_ensemble" &&
            operation != "structural_reference_calibration" &&
            operation != "source_structure_ablation" && operation != "observed_response_model" &&
            operation != "constrained_response_model" && operation != "hodge_frame_decomposition" &&
            operation != "linear_response_observation" && operation != "metric_relative_novelty" &&
            operation != "green_response")
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
        if (operation == "field_summary") {
            result.set("result", result_fields_view(field_summary_fields(query)));
            result.set("mathematical_definition", Json::string(
                "exact canonical dot(1,x), dot(x,x), and dot(1,x)/reported_count on explicit native reported support; descriptive observed rows, no population estimator or imputation"));
            return seal(result, "certificate_sha256");
        }
        if (operation == "boundary" || operation == "adjoint" || operation == "hodge" || operation == "green_action") {
            const auto display = result_fields_view(direct_action_fields(query));
            result.set("result", display.at("result"));
            result.set("output_grade", display.at("output_grade"));
            if (operation == "green_action") {
                result.set("green_iterations", display.at("green_iterations"));
                result.set("green_execution", display.at("green_execution"));
                result.set("residual_verified", display.at("residual_verified"));
            }
            result.set("mathematical_definition", Json::string(
                operation == "green_action" ? "(I + lambda L_k) x = j; exact operator action with complete output field"
                : operation == "hodge" ? "L_k x = (B_k^dagger B_k + B_(k+1) B_(k+1)^dagger) x"
                : operation == "boundary" ? "B_k x" : "M_(k+1)^-1 B_(k+1)^T M_k x"));
            return seal(result, "certificate_sha256");
        }
        if (operation == "typed_moment" || operation == "temporal_moment_change") {
            result.set("result", result_fields_view(moment_fields(query)));
            return seal(result, "certificate_sha256");
        }
        if (operation == "temporal_field_history") {
            result.set("result", history_query(query, false,nullptr,native_lease));
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
        if (is_research_operation(operation) && operation != "trajectory_structure" && operation != "paired_field_contrasts" &&
            operation != "controlled_pairings") {
            result.set("result", execute_research_query(*owner, query, false, resources));
            return seal(result, "certificate_sha256");
        }
        if (operation == "trajectory_structure") {
            result.set("result",result_fields_view(research_detail::trajectory_fields(*owner,query,resources,
                native_lease.value_or(source_execution_lease("scientific-native-trajectory")))));
        } else if (operation == "structure_transition") {
            result.set("result",result_fields_view(structure_transition_fields(query,native_lease)));
        } else if (operation == "paired_field_contrasts" || operation == "controlled_pairings") {
            result.set("result", result_fields_view(research_detail::comparison_fields(*owner, query, resources,
                native_lease.value_or(source_execution_lease("scientific-declared-comparison")))));
        } else if (operation == "chain_certificate" || operation == "realized_path_comparison") {
            result.set("result", result_fields_view(structural_fields(query, native_lease)));
        } else if (operation == "local_sections") {
            result.set("result", scientific_detail::result_fields_view(local_section_fields(query)));
        } else if (operation == "tensor_window") {
            result.set("result",result_fields_view(window_fields(query, native_lease)));
        } else if (operation == "green_response") {
            result.set("result", result_fields_view(green_response_fields(query, native_lease)));
        } else if (operation == "inspect" || operation == "provenance" || operation == "coverage") {
            result.set("result", result_fields_view(source_fields(query, native_lease)));
        } else if (operation == "support_mask") {
            const auto native = support_fields(query);
            result.set("result",result_fields_view(native));
            result.set("execution",render_record_view(native.at("execution")));
            result.set("arithmetic",Json::string("exact packed Boolean program over native reported support; supported zero is reported"));
            result.set("mathematical_definition",Json::string("declared ordered AND, OR, or NOT(OR) of reported support planes on one identical native carrier; derived eligibility only, source unmodified"));
        } else if (operation == "native_boundary") {
            const auto native = native_boundary_fields(query);
            const auto view = result_fields_view(native);
            for (const auto &name : {"result", "field_set", "retained_occurrences", "mathematical_definition"})
                result.set(name, view.at(name));
            result.set("execution", render_record_view(native.at("execution")));
        } else throw std::invalid_argument("direct scientific action dispatch mismatch");
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
    const auto admitted = render_record_view(query);
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->record(query, admitted);
}
Json ScientificSession::execute(const Json &query) {
    auto admitted = Json::parse(query.dump());
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted);
}
Json ScientificSession::execute_source(const Json &query, rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    const auto admitted = Json::parse(query.dump());
    const auto operation = admitted.at("operation").text();
    require(operation == "inspect" || operation == "provenance" || operation == "coverage",
            "source lease requires an original source operation");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
}
Json ScientificSession::execute_temporal_history(const Json &query,rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();const auto admitted=Json::parse(query.dump());
    require(admitted.at("operation").text()=="temporal_field_history","temporal lease requires a declared history");
    std::lock_guard lock(state_->mutex);state_->admit(admitted);return state_->run(admitted,std::move(lease));
}
Json ScientificSession::execute_trajectory(const Json &query,rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();const auto admitted=Json::parse(query.dump());
    require(admitted.at("operation").text()=="trajectory_structure","trajectory lease requires declared observed changes");
    std::lock_guard lock(state_->mutex);state_->admit(admitted);return state_->run(admitted,std::move(lease));
}
Json ScientificSession::execute_structure_transition(const Json &query, rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();auto admitted=Json::parse(query.dump());
    require(admitted.at("operation").text()=="structure_transition","structure lease requires a declared transition");
    std::lock_guard lock(state_->mutex);state_->admit(admitted);return state_->run(admitted,std::move(lease));
}
Json ScientificSession::execute_comparison(const Json &query, rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    const auto operation = admitted.at("operation").text();
    require(operation == "paired_field_contrasts" || operation == "controlled_pairings",
            "comparison lease requires paired or controlled fields");
    std::lock_guard lock(state_->mutex);
    state_->admit(admitted);
    return state_->run(admitted, std::move(lease));
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
Json ScientificSession::execute_green_response(const Json &query,
                                               rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    auto admitted = Json::parse(query.dump());
    require(admitted.at("operation").text() == "green_response", "Green response lease operation");
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
    const auto admitted = render_record_view(plan);
    members(admitted, {"schema", "id", "question", "source_digest", "steps"});
    require(admitted.at("schema").text() == "gadget.scientific.plan.v1", "unsupported plan schema");
    (void)text(admitted.at("id"));
    (void)text(admitted.at("question"));
    require(admitted.at("source_digest").text() == state_->source.digest, "plan source binding mismatch");
    const auto steps = admitted.at("steps");
    source_array(steps, "plan steps");
    require(steps.size() > 0, "plan requires at least one step");
    std::lock_guard lock(state_->mutex);
    std::set<std::string> identities;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        members(steps.at(i), {"id", "query"});
        require(identities.insert(text(steps.at(i).at("id"))).second, "duplicate plan step identity");
        state_->admit(steps.at(i).at("query"));
    }
    std::vector<std::pair<std::string, NativeRecord>> results;
    for (std::size_t i = 0; i < steps.size(); ++i)
        results.emplace_back(steps.at(i).at("id").text(),
            state_->record(plan.at("steps").at(i).at("query"), steps.at(i).at("query")));
    return scientific_detail::retain_scientific_plan(state_->owner->native_publication(), plan,
                                                    state_->source.digest, results);
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
NativeRecord ScientificSession::apply_proposal_record(const NativeRecord &proposal) {
    using N = NativeRecordValue;
    const auto admitted = render_record_view(proposal);
    const auto rejected = [&](const char *status, const std::exception &error) {
        N value(N::Object{{{"schema", N(std::string("gadget.agent.response.v1"))},
            {"proposal", import_record_view(admitted)}, {"status", N(std::string(status))},
            {"error", N(std::string(error.what()))}}});
        return retain_operation(state_->owner->native_publication(), proposal, std::move(value));
    };
    try {
        require(admitted.at("schema").text() == "gadget.scientific.agent.proposal.v1",
                "unsupported agent proposal schema");
        const auto action = admitted.at("action").text();
        require(action == "query" || action == "plan", "native proposal requires query or plan");
        members(admitted, {"schema", "action", action});
        const auto execution = action == "query" ? execute_record(proposal.at(action))
                                                   : execute_plan_record(proposal.at(action));
        return scientific_detail::retain_scientific_proposal(execution, proposal);
    } catch (const rfc::runtime::CpuResourceBusy &error) {
        return rejected("busy", error);
    } catch (const rfc::runtime::ExecutionLimitExceeded &error) {
        return rejected("rejected", error);
    } catch (const std::invalid_argument &error) {
        return rejected("rejected", error);
    }
}
Json apply_scientific_proposal(ScientificSession &session, const Json &proposal) {
    if (proposal.has("action") && json_object_is_type(proposal.at("action").raw(), json_type_string) &&
        (proposal.at("action").text() == "query" || proposal.at("action").text() == "plan"))
        return render_record_view(session.apply_proposal_record(admit_record_view(proposal)).at("result"));
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
        } else if (action == "annotate") {
            const std::string member = "annotations";
            members(proposal, {"schema", "action", member});
            result.set("certificate",
                       session.propose_annotations(proposal.at(member)));
            result.set("status", Json::string("prospective"));
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
