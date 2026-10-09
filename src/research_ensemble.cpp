#include "research_ensemble.hpp"
#include "research_fields.hpp"
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
Json rational(Rational value) {
    value.canonicalize();
    return Json::string(value.get_str());
}
Json vector_json(const Vector &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(rational(value));
    return out;
}
Json matrix_json(const Matrix &matrix) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(matrix.rows()));
    out.set("columns", cardinality(matrix.cols()));
    for (const auto &[position, value] : matrix.entries()) {
        auto entry = Json::object();
        entry.set("row", cardinality(position.first));
        entry.set("column", cardinality(position.second));
        entry.set("value", rational(value));
        entries.push(entry);
    }
    out.set("entries", entries);
    return out;
}
struct Scope {
    std::string id;
    std::vector<std::size_t> samples;
};
struct Reading {
    std::array<Vector, 4> fields;
    Json record;
};
Json describe(const scientific_detail::WeightedFieldEnsemble &ensemble) {
    auto out = Json::object();
    out.set("mean", vector_json(ensemble.mean));
    out.set("centered_fields", matrix_json(ensemble.centered));
    out.set("rank", cardinality(ensemble.rank));
    out.set("contraction_space",
            Json::string(ensemble.sample_space ? "samples" : "source_coordinates"));
    out.set("contraction_metric", matrix_json(ensemble.contraction_metric));
    out.set("second_moment_or_gram", matrix_json(ensemble.contraction));
    out.set("status", Json::string(ensemble.moments ? "completed" : "constant_fields"));
    out.set("variance", rational(ensemble.moments ? ensemble.moments->variance : Rational{0}));
    out.set("second_contraction",
            rational(ensemble.moments ? ensemble.moments->second_contraction : Rational{0}));
    out.set("effective_dimension_2",
            ensemble.moments ? rational(ensemble.moments->dimension) : Json{});
    out.set("collision_probability",
            ensemble.moments ? rational(ensemble.moments->collision_probability) : Json{});
    out.set("weighted_centering_exact", Json::boolean(true));
    return out;
}
} // namespace

Json relational_field_ensemble(const NativeStudy &study, const Json &query, core::Grade grade,
                               const std::vector<Vector> &fields, bool validate,
                               const std::shared_ptr<CpuResources> &resources) {
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
        return Json::object();
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
        reading.record = Json::object();
        reading.record.set("field", Json::string(ids[s]));
        auto quadrances = Json::object();
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
        auto dirac_fields = Json::array();
        for (const auto &values : dirac.grades)
            dirac_fields.push(vector_json(values));
        reading.record.set("dirac", dirac_fields);
        reading.record.set("boundary",
                           vector_json(grade.value ? dirac.grades[grade.value - 1u] : Vector{}));
        reading.record.set("upper_adjoint",
                           vector_json(grade.value < study.source().max_grade().value
                                           ? dirac.grades[grade.value + 1u]
                                           : Vector{}));
        auto green = Json::array();
        for (const auto &parameter : parameters) {
            stop.check();
            const auto response = study.green_action(grade, parameter, fields[s]);
            require(response.exact &&
                        numeric::add(response.value,
                                     numeric::scale(parameter, hodge.apply(response.value))) ==
                            fields[s],
                    "ensemble Green residual");
            auto record = Json::object();
            record.set("lambda", rational(parameter));
            record.set("field", vector_json(response.value));
            record.set("iterations", cardinality(response.iterations));
            record.set("residual_exact", Json::boolean(true));
            green.push(record);
        }
        reading.record.set("green", green);
        reading.record.set("reconstruction_exact", Json::boolean(true));
        reading.record.set("orthogonality_exact", Json::boolean(true));
        reading.record.set("dirac_square_exact", Json::boolean(true));
        reading.record.set("harmonic_kernel_exact", Json::boolean(true));
    });
    reservation.reset();
    auto serial = cpu->reserve(1);
    auto out = Json::object(), samples = Json::array(), groups = Json::array(),
         frames = Json::object();
    out.set("grade", cardinality(grade.value));
    out.set("interpretation", interpretation);
    out.set("weights", vector_json(weights));
    out.set("metric", matrix_json(metric.materialize()));
    for (const auto &reading : readings)
        samples.push(reading.record);
    out.set("samples", samples);
    for (std::size_t s = 0; s < projections.size(); ++s) {
        auto record = Json::object();
        record.set("frame", matrix_json(projections[s].frame));
        record.set("gram", matrix_json(projections[s].solve.matrix()));
        record.set("rank", cardinality(projections[s].frame.cols()));
        record.set("factor_entries", cardinality(projections[s].solve.stored_entries()));
        frames.set(channels[s + 1], record);
    }
    auto harmonic = Json::object();
    harmonic.set("rank", cardinality(metric.size() - nonharmonic_rank));
    harmonic.set("representation", Json::string("identity_minus_gradient_minus_curl"));
    harmonic.set("basis_materialized", Json::boolean(false));
    frames.set("harmonic", harmonic);
    out.set("frames", frames);
    std::vector<std::array<Vector, 4>> means;
    for (const auto &scope : scopes) {
        Vector conditional;
        Rational mass{0};
        for (const auto s : scope.samples)
            mass += weights[s];
        auto members_ = Json::array();
        for (const auto s : scope.samples) {
            conditional.push_back(weights[s] / mass);
            members_.push(Json::string(ids[s]));
        }
        auto record = Json::object(), summaries = Json::object(), couplings = Json::array();
        record.set("id", Json::string(scope.id));
        record.set("fields", members_);
        record.set("original_mass", rational(mass));
        record.set("conditional_weights", vector_json(conditional));
        std::array<Vector, 4> mean;
        std::vector<scientific_detail::WeightedFieldEnsemble> ensembles;
        for (std::size_t channel = 0; channel < channels.size(); ++channel) {
            std::vector<Vector> values;
            for (const auto s : scope.samples)
                values.push_back(readings[s].fields[channel]);
            ensembles.emplace_back(values, conditional, metric);
            mean[channel] = ensembles.back().mean;
            summaries.set(channels[channel], describe(ensembles.back()));
        }
        Rational variance{0};
        for (std::size_t s = 1; s < channels.size(); ++s)
            if (ensembles[s].moments)
                variance += ensembles[s].moments->variance;
        require(variance == (ensembles[0].moments ? ensembles[0].moments->variance : Rational{0}),
                "ensemble variance budget");
        for (std::size_t a = 1; a < channels.size(); ++a)
            for (std::size_t b = a + 1; b < channels.size(); ++b) {
                auto value = Json::object();
                value.set("left", Json::string(channels[a]));
                value.set("right", Json::string(channels[b]));
                const bool defined = ensembles[a].moments && ensembles[b].moments;
                value.set("status", Json::string(defined ? "completed" : "zero_sector_variance"));
                value.set("value", defined ? rational(scientific_detail::weighted_field_coupling(
                                                 ensembles[a], ensembles[b], conditional, metric))
                                           : Json{});
                couplings.push(value);
            }
        record.set("ensembles", summaries);
        record.set("sector_coupling", couplings);
        record.set("variance_budget_exact", Json::boolean(true));
        groups.push(record);
        means.push_back(std::move(mean));
    }
    out.set("groups", groups);
    auto comparisons = Json::array();
    if (query.has("comparisons"))
        for (std::size_t i = 0; i < query.at("comparisons").size(); ++i) {
            const auto declaration = query.at("comparisons").at(i);
            const auto left = scope_indices.at(declaration.at("left").text()),
                       right = scope_indices.at(declaration.at("right").text());
            auto record = Json::object(), changes = Json::object(), overlap = Json::array();
            record.set("declaration", declaration);
            for (const auto s : scopes[left].samples)
                if (std::find(scopes[right].samples.begin(), scopes[right].samples.end(), s) !=
                    scopes[right].samples.end())
                    overlap.push(Json::string(ids[s]));
            record.set("shared_fields", overlap);
            Vector reconstructed(metric.size(), Rational{0});
            Rational sector_quadrance{0};
            for (std::size_t channel = 0; channel < channels.size(); ++channel) {
                const auto delta = numeric::sub(means[right][channel], means[left][channel]);
                auto change = Json::object();
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
            record.set("reconstruction_exact", Json::boolean(true));
            comparisons.push(record);
        }
    out.set("comparisons", comparisons);
    auto chain = Json::array();
    for (const auto &check : prepared.bianchi.interfaces) {
        auto value = Json::object();
        value.set("lower_boundary_grade", cardinality(check.lower_boundary_grade.value));
        value.set("satisfied", Json::boolean(check.satisfied));
        chain.push(value);
    }
    out.set("chain_interfaces", chain);
    auto execution = Json::object();
    execution.set("projector_preparations", cardinality(projections.size()));
    execution.set("applications_per_projector", cardinality(fields.size()));
    execution.set("requested_workers", cardinality(report.requested_workers));
    execution.set("worker_slots", cardinality(report.worker_slots));
    execution.set("completed_tiles", cardinality(report.completed_tiles));
    execution.set("preparation",
                  Json::string("native boundary frames and prepared exact Gram solves"));
    execution.set("green", Json::string("native exact operator-first actions"));
    out.set("execution", execution);
    return out;
}
} // namespace gadget::research_detail
