#include "research_ensemble.hpp"
#include "research_fields.hpp"
#include "scientific_result_fields.hpp"
#include "scientific_ensemble.hpp"
#include "scientific_relational_field.hpp"
#include <array>
#include <gadget/resources.hpp>
#include <rfc/core/axioms.hpp>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::research_detail {
namespace {
using Matrix = numeric::ExactMatrix;
constexpr std::array<const char *, 4> channels{"full", "gradient", "curl", "harmonic"};
using Value = scientific_detail::ResultValue;
using Publication = scientific_detail::ResultFields;
using Context = scientific_detail::ResultFieldContext;
struct Scope {
    std::string id;
    std::vector<std::size_t> samples;
};
struct Reading {
    std::array<Vector, 4> fields;
    Value record;
};
Value describe(Context &context, std::size_t grade, std::uint64_t samples,
               const scientific_detail::WeightedFieldEnsemble &ensemble) {
    auto &publication = context.publication;
    const auto rows = context.grades.at(grade).rows;
    const auto contraction = ensemble.sample_space ? samples : rows;
    const auto rational = [&](const Rational &value) {return publication.stage_scalar(value);};
    auto out = Value::object();
    out.set("mean", context.field(ensemble.mean, grade));
    out.set("centered_fields", Publication::stage_matrix(ensemble.centered, rows, samples));
    out.set("rank", Value::count(ensemble.rank));
    out.set("contraction_space",
            Value::string(ensemble.sample_space ? "samples" : "source_coordinates"));
    out.set("contraction_metric", Publication::stage_matrix(ensemble.contraction_metric, contraction, contraction, storage::TemporalTensorRole::Metric));
    out.set("second_moment_or_gram", Publication::stage_matrix(ensemble.contraction, contraction, contraction));
    out.set("status", Value::string(ensemble.moments ? "completed" : "constant_fields"));
    out.set("variance", rational(ensemble.moments ? ensemble.moments->variance : Rational{0}));
    out.set("second_contraction",
            rational(ensemble.moments ? ensemble.moments->second_contraction : Rational{0}));
    out.set("effective_dimension_2",
            ensemble.moments ? rational(ensemble.moments->dimension) : Value{});
    out.set("collision_probability",
            ensemble.moments ? rational(ensemble.moments->collision_probability) : Value{});
    out.set("weighted_centering_exact", Value::boolean(true));
    return out;
}
} // namespace

Value relational_field_ensemble(const NativeStudy &study, const Json &query, core::Grade grade,
                               const std::vector<Vector> &fields, bool validate,
                               const std::shared_ptr<CpuResources> &resources, Context *context) {
    members(query,
            {"operation", "context", "fields", "weights", "interpretation", "green_parameters",
             "workers"},
            {"groups", "comparisons"});
    const auto ids = names(query.at("fields"), fields.size());
    require(ids.size() == fields.size() && !fields.empty(), "ensemble field axes");
    const auto interpretation = query.at("interpretation");
    members(interpretation, {"kind", "population", "weighting", "evidence"});
    require(interpretation.at("kind").text() == "observed_fields" ||
                interpretation.at("kind").text() == "declared_scenarios",
            "ensemble interpretation kind");
    (void)text(interpretation.at("population"), source_byte_limit());
    (void)text(interpretation.at("weighting"), source_byte_limit());
    source_array(interpretation.at("evidence"), "ensemble source evidence");
    require(interpretation.at("evidence").size(), "ensemble requires source evidence");
    const auto &metric = study.tower().metrics().at(grade.value);
    require(metric.size() > 0, "ensemble requires a nonempty source grade");
    source_array(query.at("weights"), "ensemble weights");
    require(query.at("weights").size() == fields.size(), "ensemble weight axes");
    Vector weights;
    Rational total{0};
    for (std::size_t s = 0; s < fields.size(); ++s) {
        const auto weight = exact_rational(text(query.at("weights").at(s), source_byte_limit()));
        require(weight > 0 && fields[s].size() == metric.size(),
                "ensemble positive weights and complete carrier");
        weights.push_back(weight);
        total += weight;
    }
    require(total == 1, "ensemble weights must sum exactly to one");
    std::vector<Scope> scopes{{"all", {}}};
    for (std::size_t s = 0; s < fields.size(); ++s)
        scopes.front().samples.push_back(s);
    std::map<std::string, std::size_t> scope_indices{{"all", 0}}, positions;
    for (std::size_t s = 0; s < ids.size(); ++s)
        positions.emplace(ids[s], s);
    if (query.has("groups")) {
        source_array(query.at("groups"), "ensemble groups");
        for (std::size_t g = 0; g < query.at("groups").size(); ++g) {
            const auto declaration = query.at("groups").at(g);
            members(declaration, {"id", "fields"});
            Scope scope{text(declaration.at("id"), source_byte_limit()), {}};
            require(scope_indices.emplace(scope.id, scopes.size()).second,
                    "duplicate ensemble group");
            for (const auto &id : names(declaration.at("fields"), fields.size())) {
                require(positions.contains(id), "ensemble group names an unselected field");
                scope.samples.push_back(positions.at(id));
            }
            scopes.push_back(std::move(scope));
        }
    }
    if (query.has("comparisons")) {
        source_array(query.at("comparisons"), "ensemble comparisons");
        std::set<std::string> comparison_ids;
        for (std::size_t i = 0; i < query.at("comparisons").size(); ++i) {
            const auto comparison = query.at("comparisons").at(i);
            members(comparison, {"id", "left", "right"});
            require(comparison_ids.insert(text(comparison.at("id"), source_byte_limit())).second,
                    "duplicate ensemble comparison");
            require(scope_indices.contains(comparison.at("left").text()) &&
                        scope_indices.contains(comparison.at("right").text()),
                    "unknown ensemble comparison group");
        }
    }
    Vector parameters;
    std::set<Rational> unique;
    source_array(query.at("green_parameters"), "ensemble Green parameters");
    for (std::size_t i = 0; i < query.at("green_parameters").size(); ++i) {
        const auto value =
            exact_rational(text(query.at("green_parameters").at(i), source_byte_limit()));
        require(value > 0 && unique.insert(value).second,
                "ensemble positive unique Green parameters");
        parameters.push_back(value);
    }
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        count(query.at("workers"), cpu->coordinator().capacity(), "ensemble workers");
    require(workers > 0, "ensemble workers must be positive");
    if (validate)
        return Value::object();
    require(context != nullptr, "ensemble native publication owner required");
    auto &publication = context->publication;
    const auto rows = context->grades.at(grade.value).rows;
    const auto sample_axis = publication.axis(fields.size());
    auto sample_binding = Value::object(), sample_names = Value::array();
    for (const auto &id : ids) sample_names.push(Value::string(id));
    sample_binding.set("sample_carrier", Value::count(sample_axis));
    sample_binding.set("fields", std::move(sample_names));
    publication.bind(std::move(sample_binding));
    const auto empty_axis = publication.axis(0);
    const auto rational = [&](const Rational &value) {return publication.stage_scalar(value);};
    const auto vector_json = [&](const Vector &value) {return context->field(value, grade.value);};
    std::optional<runtime::CpuResourceCoordinator::Reservation> reservation;
    reservation.emplace(cpu->reserve(workers));
    const scientific_detail::PreparedRelationalField prepared(study, grade);
    const auto &projections = prepared.projections;
    const auto hodge = study.hodge(grade);
    const auto nonharmonic_rank =
        machine::checked_size_add(projections[0].frame.cols(), projections[1].frame.cols());
    require(nonharmonic_rank <= metric.size(), "Hodge complementary dimension");
    std::vector<Reading> readings(fields.size());
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), fields.size()});
    const auto report = executor.run_tiles(fields.size(), [&](std::size_t s, std::size_t,
                                                              const runtime::CpuTileStop &stop) {
        auto &reading = readings[s];
        reading.fields = prepared.apply(fields[s], [&] { stop.check(); });
        Rational sector_energy{0};
        reading.record = Value::object();
        reading.record.set("field", Value::string(ids[s]));
        auto quadrances = Value::object();
        for (std::size_t sector = 0; sector < channels.size(); ++sector) {
            reading.record.set(channels[sector], vector_json(reading.fields[sector]));
            const auto quadrance = metric.inner(reading.fields[sector], reading.fields[sector]);
            quadrances.set(channels[sector], rational(quadrance));
            if (sector)
                sector_energy += quadrance;
        }
        reading.record.set("quadrances", quadrances);
        for (std::size_t a = 1; a < channels.size(); ++a)
            for (std::size_t b = a + 1; b < channels.size(); ++b)
                require(metric.inner(reading.fields[a], reading.fields[b]) == 0,
                        "ensemble sector orthogonality");
        require(sector_energy == metric.inner(fields[s], fields[s]),
                "ensemble sector quadrance budget");
        const auto down = hodge.apply_down(fields[s]), up = hodge.apply_up(fields[s]);
        reading.record.set("hodge_down", vector_json(down));
        reading.record.set("hodge_up", vector_json(up));
        reading.record.set("hodge", vector_json(numeric::add(down, up)));
        auto tower = field::zero_exact_tower(study.source());
        tower.grades.at(grade.value) = fields[s];
        const auto dirac = study.tower().apply_dirac(tower);
        const auto square = study.tower().apply_dirac(dirac);
        const auto laplacian = study.tower().apply_hodge_tower(tower);
        require(square.grades == laplacian.grades, "ensemble Dirac square identity");
        auto dirac_fields = Value::array();
        for (std::size_t k = 0; k < dirac.grades.size(); ++k)
            dirac_fields.push(context->field(dirac.grades[k], k));
        reading.record.set("dirac", dirac_fields);
        reading.record.set("boundary",
                           grade.value ? context->field(dirac.grades[grade.value - 1u], grade.value - 1u)
                                       : Publication::stage_vector(Vector{}, empty_axis));
        reading.record.set("upper_adjoint",
                           grade.value < study.source().max_grade().value
                               ? context->field(dirac.grades[grade.value + 1u], grade.value + 1u)
                               : Publication::stage_vector(Vector{}, empty_axis));
        auto green = Value::array();
        for (const auto &parameter : parameters) {
            stop.check();
            const auto response = study.green_action(grade, parameter, fields[s]);
            require(response.exact &&
                        numeric::add(response.value,
                                     numeric::scale(parameter, hodge.apply(response.value))) ==
                            fields[s],
                    "ensemble Green residual");
            auto record = Value::object();
            record.set("lambda", rational(parameter));
            record.set("field", vector_json(response.value));
            record.set("iterations", Value::count(response.iterations));
            record.set("residual_exact", Value::boolean(true));
            green.push(record);
        }
        reading.record.set("green", green);
        reading.record.set("reconstruction_exact", Value::boolean(true));
        reading.record.set("orthogonality_exact", Value::boolean(true));
        reading.record.set("dirac_square_exact", Value::boolean(true));
        reading.record.set("harmonic_kernel_exact", Value::boolean(true));
    });
    reservation.reset();
    auto serial = cpu->reserve(1);
    auto out = Value::object(), samples = Value::array(), groups = Value::array(),
         frames = Value::object();
    out.set("grade", Value::count(grade.value));
    out.set("interpretation", interpretation);
    out.set("weights", Publication::stage_vector(weights, sample_axis));
    out.set("metric", Publication::stage_matrix(metric.materialize(), rows, rows, storage::TemporalTensorRole::Metric));
    for (const auto &reading : readings)
        samples.push(reading.record);
    out.set("samples", samples);
    for (std::size_t s = 0; s < projections.size(); ++s) {
        const auto frame_axis = publication.axis(projections[s].frame.cols());
        auto record = Value::object();
        record.set("frame", Publication::stage_matrix(projections[s].frame, frame_axis, rows));
        record.set("gram", Publication::stage_matrix(projections[s].solve.matrix(), frame_axis, frame_axis, storage::TemporalTensorRole::Metric));
        record.set("rank", Value::count(projections[s].frame.cols()));
        record.set("factor_entries", Value::count(projections[s].solve.stored_entries()));
        frames.set(channels[s + 1], record);
    }
    auto harmonic = Value::object();
    harmonic.set("rank", Value::count(metric.size() - nonharmonic_rank));
    harmonic.set("representation", Value::string("identity_minus_gradient_minus_curl"));
    harmonic.set("basis_materialized", Value::boolean(false));
    frames.set("harmonic", harmonic);
    out.set("frames", frames);
    std::vector<std::array<Vector, 4>> means;
    for (const auto &scope : scopes) {
        Vector conditional;
        Rational mass{0};
        for (const auto s : scope.samples)
            mass += weights[s];
        const auto group_axis = publication.axis(scope.samples.size());
        auto binding = Value::object();
        binding.set("group", Value::string(scope.id));
        binding.set("sample_carrier", Value::count(sample_axis));
        binding.set("group_carrier", Value::count(group_axis));
        binding.set("sample_coordinates", publication.indices(scope.samples, group_axis));
        publication.bind(std::move(binding));
        auto members_ = Value::array();
        for (const auto s : scope.samples) {
            conditional.push_back(weights[s] / mass);
            members_.push(Value::string(ids[s]));
        }
        auto record = Value::object(), summaries = Value::object(), couplings = Value::array();
        record.set("id", Value::string(scope.id));
        record.set("fields", members_);
        record.set("original_mass", rational(mass));
        record.set("conditional_weights", Publication::stage_vector(conditional, group_axis));
        std::array<Vector, 4> mean;
        std::vector<scientific_detail::WeightedFieldEnsemble> ensembles;
        for (std::size_t channel = 0; channel < channels.size(); ++channel) {
            std::vector<Vector> values;
            for (const auto s : scope.samples)
                values.push_back(readings[s].fields[channel]);
            ensembles.emplace_back(values, conditional, metric);
            mean[channel] = ensembles.back().mean;
            summaries.set(channels[channel], describe(*context, grade.value, group_axis, ensembles.back()));
        }
        Rational variance{0};
        for (std::size_t s = 1; s < channels.size(); ++s)
            if (ensembles[s].moments)
                variance += ensembles[s].moments->variance;
        require(variance == (ensembles[0].moments ? ensembles[0].moments->variance : Rational{0}),
                "ensemble variance budget");
        for (std::size_t a = 1; a < channels.size(); ++a)
            for (std::size_t b = a + 1; b < channels.size(); ++b) {
                auto value = Value::object();
                value.set("left", Value::string(channels[a]));
                value.set("right", Value::string(channels[b]));
                const bool defined = ensembles[a].moments && ensembles[b].moments;
                value.set("status", Value::string(defined ? "completed" : "zero_sector_variance"));
                value.set("value", defined ? rational(scientific_detail::weighted_field_coupling(
                                                 ensembles[a], ensembles[b], conditional, metric))
                                           : Value{});
                if (ensembles[a].sample_space) {
                    Matrix terms(scope.samples.size(), scope.samples.size());
                    for (std::size_t i = 0; i < scope.samples.size(); ++i)
                        for (std::size_t j = 0; j < scope.samples.size(); ++j)
                            terms.set(i, j, conditional[i] * conditional[j] *
                                ensembles[a].contraction.get(i, j) * ensembles[b].contraction.get(i, j));
                    value.set("ordered_contributions", Publication::stage_matrix(terms, group_axis, group_axis));
                } else {
                    const auto cross = field::metric_gram(ensembles[a].centered,
                        field::ExactMetric::diagonal(conditional), ensembles[b].centered);
                    value.set("cross_covariance", Publication::stage_matrix(cross, rows, rows));
                }
                couplings.push(value);
            }
        record.set("ensembles", summaries);
        record.set("sector_coupling", couplings);
        record.set("variance_budget_exact", Value::boolean(true));
        groups.push(record);
        means.push_back(std::move(mean));
    }
    out.set("groups", groups);
    auto comparisons = Value::array();
    if (query.has("comparisons"))
        for (std::size_t i = 0; i < query.at("comparisons").size(); ++i) {
            const auto declaration = query.at("comparisons").at(i);
            const auto left = scope_indices.at(declaration.at("left").text()),
                       right = scope_indices.at(declaration.at("right").text());
            auto record = Value::object(), changes = Value::object(), overlap = Value::array();
            record.set("declaration", declaration);
            for (const auto s : scopes[left].samples)
                if (std::find(scopes[right].samples.begin(), scopes[right].samples.end(), s) !=
                    scopes[right].samples.end())
                    overlap.push(Value::string(ids[s]));
            record.set("shared_fields", overlap);
            Vector reconstructed(metric.size(), Rational{0});
            Rational sector_quadrance{0};
            for (std::size_t channel = 0; channel < channels.size(); ++channel) {
                const auto delta = numeric::sub(means[right][channel], means[left][channel]);
                auto change = Value::object();
                change.set("field", vector_json(delta));
                change.set("quadrance", rational(metric.inner(delta, delta)));
                changes.set(channels[channel], change);
                if (channel) {
                    reconstructed = numeric::add(reconstructed, delta);
                    sector_quadrance += metric.inner(delta, delta);
                }
            }
            const auto delta = numeric::sub(means[right][0], means[left][0]);
            require(reconstructed == delta && sector_quadrance == metric.inner(delta, delta),
                    "ensemble mean contrast reconstruction");
            record.set("right_minus_left", changes);
            record.set("reconstruction_exact", Value::boolean(true));
            comparisons.push(record);
        }
    out.set("comparisons", comparisons);
    auto chain = Value::array();
    for (const auto &check : prepared.bianchi.interfaces) {
        auto value = Value::object();
        value.set("lower_boundary_grade", Value::count(check.lower_boundary_grade.value));
        value.set("satisfied", Value::boolean(check.satisfied));
        chain.push(value);
    }
    out.set("chain_interfaces", chain);
    auto execution = Value::object();
    execution.set("projector_preparations", Value::count(projections.size()));
    execution.set("applications_per_projector", Value::count(fields.size()));
    execution.set("requested_workers", Value::count(report.requested_workers));
    execution.set("worker_slots", Value::count(report.worker_slots));
    execution.set("completed_tiles", Value::count(report.completed_tiles));
    execution.set("preparation",
                  Value::string("native boundary frames and prepared exact Gram solves"));
    execution.set("green", Value::string("native exact operator-first actions"));
    out.set("execution", execution);
    return out;
}
} // namespace gadget::research_detail
