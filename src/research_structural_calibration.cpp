#include "research_structural_calibration.hpp"
#include "research_fields.hpp"
#include "scientific_ensemble.hpp"
#include "scientific_reference.hpp"
#include "scientific_relational_field.hpp"
#include <algorithm>
#include <gadget/resources.hpp>
#include <limits>
#include <numeric>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::research_detail {
namespace {
using Matrix = numeric::ExactMatrix;
constexpr std::array<const char *, 4> sectors{"full", "gradient", "curl", "harmonic"};
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
Json indices_json(const std::vector<std::size_t> &indices) {
    auto out = Json::array();
    for (const auto index : indices)
        out.push(cardinality(index));
    return out;
}
struct Hypothesis {
    std::string statistic;
    std::size_t axis, sector{0};
    Rational lambda{0};
    std::vector<std::string> probes;
    std::shared_ptr<const scientific_detail::PreparedProjection> projection;
};
struct Scenario {
    std::string id;
    Rational weight;
    std::vector<std::string> fields;
};
void interpretation(const Json &value, const std::set<std::string> &keys) {
    for (const auto &key_ : keys)
        (void)text(value.at(key_), source_byte_limit());
    (void)names(value.at("assumptions"), value.at("assumptions").size());
    source_array(value.at("evidence"), "structural calibration evidence");
    require(value.at("evidence").size(), "structural calibration requires source evidence");
}
std::size_t sector_index(const Json &value) {
    const auto name = value.text();
    const auto it = std::find(sectors.begin(), sectors.end(), name);
    require(it != sectors.end(), "unknown structural sector");
    return static_cast<std::size_t>(it - sectors.begin());
}
// Factorial unranking retains coordinate multiplicity even when values coincide.
// This constructs declared experimental assignments, not mathematical kernels.
std::vector<std::size_t> permutation(std::vector<std::size_t> available, std::size_t rank) {
    std::vector<std::size_t> factors(available.size() + 1, 1), out;
    for (std::size_t i = 1; i <= available.size(); ++i)
        factors[i] = machine::checked_size_mul(factors[i - 1], i);
    while (!available.empty()) {
        const auto factor = factors[available.size() - 1], index = rank / factor;
        rank %= factor;
        out.push_back(available.at(index));
        available.erase(available.begin() + static_cast<std::ptrdiff_t>(index));
    }
    return out;
}
} // namespace
std::set<std::string> structural_calibration_field_ids(const Json &query) {
    const auto axes = names(query.at("fields"), query.at("fields").size());
    std::set<std::string> ids(axes.begin(), axes.end());
    const auto hypotheses = query.at("hypotheses");
    source_array(hypotheses, "structural hypotheses");
    for (std::size_t h = 0; h < hypotheses.size(); ++h)
        if (hypotheses.at(h).has("probes"))
            for (const auto &id :
                 names(hypotheses.at(h).at("probes"), hypotheses.at(h).at("probes").size()))
                ids.insert(id);
    if (query.has("uncertainty")) {
        const auto scenarios = query.at("uncertainty").at("scenarios");
        source_array(scenarios, "joint uncertainty scenarios");
        for (std::size_t s = 0; s < scenarios.size(); ++s)
            for (const auto &id :
                 names(scenarios.at(s).at("fields"), scenarios.at(s).at("fields").size()))
                ids.insert(id);
    }
    return ids;
}
Json structural_reference_calibration(const NativeStudy &study, const Json &query,
                                      core::Grade grade,
                                      const std::map<std::string, Vector> &fields,
                                      const std::vector<ExchangeBlock> &blocks, bool validate,
                                      const std::shared_ptr<CpuResources> &resources) {
    members(query,
            {"operation", "context", "fields", "law", "hypotheses", "workers", "execution",
             "reference_output"},
            {"alpha", "uncertainty"});
    const auto axes = names(query.at("fields"), fields.size());
    const auto &metric = study.tower().metrics().at(grade.value);
    const auto law = query.at("law");
    members(law, {"kind", "blocks", "null_hypothesis", "population", "assumptions", "derivation",
                  "evidence"});
    require(law.at("kind").text() == "blocked_permutations", "structural law kind");
    interpretation(law, {"null_hypothesis", "population", "derivation"});
    require(!blocks.empty() && blocks.size() == law.at("blocks").size(),
            "exchange blocks required");
    numeric::Integer state_count{1}, assignment_coordinates{0};
    std::vector<numeric::Integer> radices;
    for (const auto &block : blocks) {
        require(!block.coordinates.empty(), "empty exchange block");
        numeric::Integer radix{1};
        for (std::size_t i = 1; i <= block.coordinates.size(); ++i)
            radix *= numeric::Integer{std::to_string(i)};
        radices.push_back(radix);
        state_count *= radix;
        assignment_coordinates += numeric::Integer{std::to_string(block.coordinates.size())};
    }
    const auto family = query.at("hypotheses");
    source_array(family, "structural hypothesis family");
    require(family.size(), "structural calibration requires hypotheses");
    std::vector<Hypothesis> hypotheses;
    std::set<std::string> hypothesis_ids;
    for (std::size_t h = 0; h < family.size(); ++h) {
        const auto declaration = family.at(h);
        Hypothesis hypothesis;
        hypothesis.statistic = text(declaration.at("statistic"), source_byte_limit());
        if (hypothesis.statistic == "sector_fraction")
            members(declaration, {"id", "field", "statistic", "sector"});
        else if (hypothesis.statistic == "probe_fraction")
            members(declaration, {"id", "field", "statistic", "sector", "probes"});
        else {
            require(hypothesis.statistic == "green_change_fraction",
                    "unknown structural statistic");
            members(declaration, {"id", "field", "statistic", "lambda"});
            hypothesis.lambda = exact_rational(text(declaration.at("lambda"), source_byte_limit()));
            require(hypothesis.lambda > 0, "positive Green parameter required");
        }
        require(hypothesis_ids.insert(text(declaration.at("id"), source_byte_limit())).second,
                "duplicate structural hypothesis identity");
        const auto field = declaration.at("field").text();
        const auto found = std::find(axes.begin(), axes.end(), field);
        require(found != axes.end(), "structural hypothesis field is not a selected axis");
        hypothesis.axis = static_cast<std::size_t>(found - axes.begin());
        if (declaration.has("sector"))
            hypothesis.sector = sector_index(declaration.at("sector"));
        if (declaration.has("probes"))
            hypothesis.probes = names(declaration.at("probes"), fields.size());
        hypotheses.push_back(std::move(hypothesis));
    }
    std::vector<Scenario> scenarios;
    std::size_t observed_scenario{0};
    if (query.has("uncertainty")) {
        const auto uncertainty = query.at("uncertainty");
        members(uncertainty, {"kind", "observed_scenario", "population", "weighting", "derivation",
                              "assumptions", "evidence", "scenarios"});
        require(uncertainty.at("kind").text() == "declared_joint_scenarios",
                "uncertainty interpretation kind");
        interpretation(uncertainty, {"population", "weighting", "derivation"});
        const auto declarations = uncertainty.at("scenarios");
        source_array(declarations, "joint uncertainty scenarios");
        std::set<std::string> ids;
        std::optional<std::size_t> observed;
        Rational total{0};
        for (std::size_t s = 0; s < declarations.size(); ++s) {
            const auto value = declarations.at(s);
            members(value, {"id", "weight", "fields"});
            Scenario scenario{text(value.at("id"), source_byte_limit()),
                              exact_rational(text(value.at("weight"), source_byte_limit())),
                              names(value.at("fields"), fields.size())};
            require(ids.insert(scenario.id).second && scenario.weight > 0 &&
                        scenario.fields.size() == axes.size(),
                    "joint scenario identity, positive weight and complete field axes");
            if (scenario.id == uncertainty.at("observed_scenario").text())
                observed = s;
            total += scenario.weight;
            scenarios.push_back(std::move(scenario));
        }
        require(total == 1 && observed.has_value(),
                "joint scenarios must sum exactly to one and retain the observed scenario");
        observed_scenario = *observed;
        for (std::size_t a = 0; a < axes.size(); ++a)
            require(fields.at(scenarios[observed_scenario].fields[a]) == fields.at(axes[a]),
                    "observed uncertainty scenario differs from original field");
    } else
        scenarios.push_back({"observed", Rational{1}, axes});
    for (const auto &[id, field] : fields) {
        (void)id;
        require(field.size() == metric.size(), "structural calibration complete field carrier");
    }
    std::optional<Rational> alpha;
    if (query.has("alpha")) {
        alpha = exact_rational(text(query.at("alpha"), source_byte_limit()));
        require(*alpha > 0 && *alpha < 1,
                "structural alpha must lie strictly between zero and one");
    }
    const auto execution = query.at("execution").text();
    const auto profile = query.at("reference_output").text();
    require(execution == "plan" || execution == "exhaustive", "structural execution profile");
    require(profile == "statistics" || profile == "full_fields",
            "structural reference output profile");
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        count(query.at("workers"), cpu->coordinator().capacity(), "structural workers");
    require(workers, "positive structural worker budget required");
    if (validate)
        return Json::object();
    auto out = Json::object(), cost = Json::object(), witnesses = Json::array();
    out.set("schema", Json::string("gadget.structural.reference.calibration.v1"));
    out.set("arithmetic", Json::string("exact_rational"));
    out.set("law", law);
    out.set("grade", cardinality(grade.value));
    cost.set("reference_states", exact_count(state_count));
    cost.set("scenarios", cardinality(scenarios.size()));
    const numeric::Integer tasks = state_count * numeric::Integer{std::to_string(scenarios.size())};
    cost.set("field_applications",
             exact_count(tasks * numeric::Integer{std::to_string(axes.size())}));
    cost.set("statistics",
             exact_count(tasks * numeric::Integer{std::to_string(hypotheses.size())}));
    cost.set("assignment_coordinates", exact_count(state_count * assignment_coordinates));
    cost.set("full_field_coordinates",
             exact_count(profile == "full_fields"
                             ? tasks * numeric::Integer{std::to_string(axes.size())} *
                                   numeric::Integer{std::to_string(metric.size())} * 4
                             : numeric::Integer{0}));
    cost.set("reference_output", Json::string(profile));
    cost.set("method",
             Json::string("Complete product of factorial block assignments; native preparations "
                          "and exact bit growth have additional source-dependent costs"));
    for (const auto &block : blocks)
        witnesses.push(block.witness);
    out.set("exchange_blocks", witnesses);
    out.set("cost", cost);
    if (execution == "plan") {
        out.set("status", Json::string("planned"));
        return out;
    }
    require(tasks <= numeric::Integer{std::to_string(std::numeric_limits<std::size_t>::max())},
            "exhaustive reference task count exceeds machine addressability; use execution plan to "
            "inspect the declared law");
    const auto states = static_cast<std::size_t>(state_count.get_ui());
    const auto task_count = machine::checked_size_mul(states, scenarios.size());
    (void)machine::checked_size_mul(task_count, hypotheses.size());
    std::vector<std::size_t> machine_radices;
    for (const auto &radix : radices)
        machine_radices.push_back(static_cast<std::size_t>(radix.get_ui()));
    std::optional<runtime::CpuResourceCoordinator::Reservation> reservation;
    reservation.emplace(cpu->reserve(workers));
    const scientific_detail::PreparedRelationalField prepared(study, grade);
    std::map<std::pair<std::size_t, std::vector<std::string>>,
             std::shared_ptr<const scientific_detail::PreparedProjection>>
        probe_preparations;
    std::map<std::string, std::array<Vector, 4>> probe_readings;
    for (auto &hypothesis : hypotheses)
        if (!hypothesis.probes.empty()) {
            auto names_ = hypothesis.probes;
            std::sort(names_.begin(), names_.end());
            const auto binding = std::make_pair(hypothesis.sector, names_);
            if (!probe_preparations.contains(binding)) {
                Matrix frame(metric.size(), names_.size());
                for (std::size_t p = 0; p < names_.size(); ++p) {
                    if (!probe_readings.contains(names_[p]))
                        probe_readings.emplace(names_[p], prepared.apply(fields.at(names_[p])));
                    const auto &projection = probe_readings.at(names_[p]);
                    for (std::size_t c = 0; c < metric.size(); ++c)
                        frame.set(c, p, projection[hypothesis.sector][c]);
                }
                probe_preparations.emplace(
                    binding, std::make_shared<const scientific_detail::PreparedProjection>(
                                 frame.image_frame(), metric));
            }
            hypothesis.projection = probe_preparations.at(binding);
        }
    std::vector<Vector> scores(task_count, Vector(hypotheses.size(), Rational{0}));
    // uint8_t storage permits independent task writes; vector<bool> packed bits do not.
    std::vector<std::vector<std::uint8_t>> defined(task_count,
                                                   std::vector<std::uint8_t>(hypotheses.size(), 1));
    std::vector<Json> maps(states), details(profile == "full_fields" ? task_count : 0);
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), task_count});
    const auto report = executor.run_tiles(task_count, [&](std::size_t task, std::size_t,
                                                           const runtime::CpuTileStop &stop) {
        const auto scenario = task / states, state = task % states;
        Matrix action = Matrix::identity(metric.size());
        auto map = Json::array();
        auto rank = state;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            const auto &coordinates = blocks[b].coordinates;
            const auto assigned = permutation(coordinates, rank % machine_radices[b]);
            rank /= machine_radices[b];
            for (std::size_t i = 0; i < coordinates.size(); ++i) {
                action.set(coordinates[i], coordinates[i], Rational{0});
                action.set(coordinates[i], assigned[i], Rational{1});
            }
            map.push(indices_json(assigned));
        }
        if (scenario == 0)
            maps[state] = map;
        std::vector<std::array<Vector, 4>> readings;
        for (const auto &field : scenarios[scenario].fields)
            readings.push_back(
                prepared.apply(action.apply(fields.at(field)), [&] { stop.check(); }));
        auto record = Json::object(), field_records = Json::array(),
             hypothesis_records = Json::array();
        if (profile == "full_fields")
            for (std::size_t a = 0; a < axes.size(); ++a) {
                auto field = Json::object();
                field.set("axis", Json::string(axes[a]));
                field.set("source_field", Json::string(scenarios[scenario].fields[a]));
                for (std::size_t c = 0; c < sectors.size(); ++c)
                    field.set(sectors[c], vector_json(readings[a][c]));
                field.set("reconstruction_exact", Json::boolean(true));
                field.set("orthogonality_exact", Json::boolean(true));
                field.set("harmonic_kernel_exact", Json::boolean(true));
                field_records.push(field);
            }
        std::map<std::pair<std::size_t, Rational>, Vector> green_readings;
        for (std::size_t h = 0; h < hypotheses.size(); ++h) {
            stop.check();
            const auto &hypothesis = hypotheses[h];
            const auto &reading = readings[hypothesis.axis];
            const auto &x = reading[0];
            Rational denominator = metric.inner(x, x), numerator{0};
            auto detail = Json::object();
            if (hypothesis.statistic == "sector_fraction")
                numerator = metric.inner(reading[hypothesis.sector], reading[hypothesis.sector]);
            else if (hypothesis.statistic == "probe_fraction") {
                const auto &sector = reading[hypothesis.sector];
                const auto captured = hypothesis.projection->apply(sector, [&] { stop.check(); });
                const auto residual = numeric::sub(sector, captured);
                denominator = metric.inner(sector, sector);
                numerator = metric.inner(captured, captured);
                require(metric.inner(captured, residual) == 0 &&
                            numerator + metric.inner(residual, residual) == denominator,
                        "probe span orthogonality and quadrance budget");
                if (profile == "full_fields") {
                    detail.set("captured", vector_json(captured));
                    detail.set("remaining", vector_json(residual));
                    detail.set("orthogonality_exact", Json::boolean(true));
                }
            } else {
                const auto binding = std::make_pair(hypothesis.axis, hypothesis.lambda);
                if (!green_readings.contains(binding)) {
                    const auto green = study.green_action(grade, hypothesis.lambda, x);
                    require(green.exact &&
                                numeric::add(green.value, numeric::scale(hypothesis.lambda,
                                                                         study.hodge(grade).apply(
                                                                             green.value))) == x,
                            "structural Green residual");
                    green_readings.emplace(binding, green.value);
                }
                const auto &green = green_readings.at(binding);
                const auto change = numeric::sub(x, green);
                numerator = metric.inner(change, change);
                if (profile == "full_fields") {
                    detail.set("green", vector_json(green));
                    detail.set("change", vector_json(change));
                    detail.set("residual_exact", Json::boolean(true));
                }
            }
            detail.set("numerator", rational(numerator));
            detail.set("denominator", rational(denominator));
            if (denominator == 0)
                defined[task][h] = 0;
            else {
                scores[task][h] = numerator / denominator;
                require(scores[task][h] >= 0 && scores[task][h] <= 1,
                        "structural dimensionless statistic bound");
            }
            if (profile == "full_fields")
                hypothesis_records.push(detail);
        }
        if (profile == "full_fields") {
            record.set("fields", field_records);
            record.set("hypotheses", hypothesis_records);
            details[task] = record;
        }
    });
    reservation.reset();
    auto serial = cpu->reserve(1);
    const Rational mass{numeric::Integer{1}, state_count};
    const Vector masses(states, mass);
    auto state_records = Json::array();
    for (std::size_t s = 0; s < states; ++s) {
        auto record = Json::object();
        record.set("index", cardinality(s));
        record.set("mass", rational(mass));
        record.set("block_source_coordinates", maps[s]);
        state_records.push(record);
    }
    out.set("states", state_records);
    out.set("observed_state", cardinality(0));
    out.set("map_contract",
            Json::string(
                "For each block, target coordinates are its declared coordinate list and source "
                "coordinates are the state list; every other coordinate is fixed. One shared "
                "assignment acts on all field axes. Probes and the source complex stay fixed. "
                "These reassignments are not asserted to be chain maps or metric isometries."));
    auto scenario_records = Json::array();
    for (std::size_t s = 0; s < scenarios.size(); ++s) {
        std::vector<Vector> values;
        std::vector<std::vector<bool>> available;
        auto records = Json::array();
        for (std::size_t state = 0; state < states; ++state) {
            const auto task = s * states + state;
            values.push_back(scores[task]);
            available.emplace_back(defined[task].begin(), defined[task].end());
            auto record = Json::object(), statistics = Json::array();
            record.set("state", cardinality(state));
            for (std::size_t h = 0; h < hypotheses.size(); ++h)
                statistics.push(defined[task][h] ? rational(scores[task][h]) : Json{});
            record.set("statistics", statistics);
            if (profile == "full_fields")
                record.set("readings", details[task]);
            records.push(record);
        }
        auto calibration =
            scientific_detail::calibrate_finite_scores(values, available, masses, 0, family, alpha);
        calibration.set("id", Json::string(scenarios[s].id));
        calibration.set("weight", rational(scenarios[s].weight));
        calibration.set("statistics", records);
        scenario_records.push(calibration);
    }
    const auto observed = scenario_records.at(observed_scenario);
    out.set("status", observed.at("status"));
    out.set("hypotheses", observed.at("hypotheses"));
    out.set("joint_family", observed.at("joint_family"));
    out.set("scenarios", scenario_records);
    out.set("observed_scenario", Json::string(scenarios[observed_scenario].id));
    if (alpha)
        out.set("alpha", rational(*alpha));
    if (query.has("uncertainty")) {
        out.set("uncertainty", query.at("uncertainty"));
        auto sensitivity = Json::array();
        for (std::size_t h = 0; h < hypotheses.size(); ++h) {
            auto summary = Json::object();
            summary.set("hypothesis", family.at(h).at("id"));
            bool complete = true;
            Vector values, tails;
            Rational mean{0}, variance{0};
            for (std::size_t s = 0; s < scenarios.size(); ++s) {
                const auto value = scenario_records.at(s).at("hypotheses").at(h);
                if (!value.at("observed").raw() || !value.at("reference_tail_mass").raw()) {
                    complete = false;
                    break;
                }
                values.push_back(exact_rational(value.at("observed").text()));
                tails.push_back(exact_rational(value.at("reference_tail_mass").text()));
            }
            summary.set("status",
                        Json::string(complete ? "completed" : "undefined_scenario_statistic"));
            if (complete) {
                std::vector<Vector> samples;
                Vector weights;
                for (std::size_t s = 0; s < scenarios.size(); ++s) {
                    samples.push_back(Vector{values[s]});
                    weights.push_back(scenarios[s].weight);
                }
                const scientific_detail::WeightedFieldEnsemble moments(
                    samples, weights, field::ExactMetric::identity(1));
                mean = moments.mean[0];
                variance = moments.moments ? moments.moments->variance : Rational{0};
                summary.set("minimum", rational(*std::min_element(values.begin(), values.end())));
                summary.set("maximum", rational(*std::max_element(values.begin(), values.end())));
                summary.set("weighted_mean", rational(mean));
                summary.set("weighted_variance", rational(variance));
                summary.set("minimum_conditional_tail",
                            rational(*std::min_element(tails.begin(), tails.end())));
                summary.set("maximum_conditional_tail",
                            rational(*std::max_element(tails.begin(), tails.end())));
            }
            sensitivity.push(summary);
        }
        out.set("scenario_sensitivity", sensitivity);
        bool complete_observed = true;
        std::vector<Vector> observed_scores;
        Vector scenario_weights;
        for (std::size_t s = 0; s < scenarios.size(); ++s) {
            observed_scores.push_back(scores[s * states]);
            scenario_weights.push_back(scenarios[s].weight);
            complete_observed = complete_observed &&
                                std::all_of(defined[s * states].begin(), defined[s * states].end(),
                                            [](std::uint8_t value) { return value != 0; });
        }
        if (complete_observed)
            out.set(
                "scenario_observable_ensemble",
                scientific_detail::score_ensemble(observed_scores, scenario_weights,
                                                  "declared joint sensitivity weights; separate "
                                                  "from reference probability masses"));
        else {
            auto status = Json::object();
            status.set("status", Json::string("undefined_observed_scenario_statistic"));
            out.set("scenario_observable_ensemble", status);
        }
    }
    auto chain = Json::array();
    for (const auto &check : prepared.bianchi.interfaces) {
        auto value = Json::object();
        value.set("lower_boundary_grade", cardinality(check.lower_boundary_grade.value));
        value.set("satisfied", Json::boolean(check.satisfied));
        chain.push(value);
    }
    out.set("chain_interfaces", chain);
    auto frames = Json::object();
    for (std::size_t p = 0; p < prepared.projections.size(); ++p) {
        auto frame = Json::object();
        frame.set("rank", cardinality(prepared.projections[p].frame.cols()));
        frame.set("factor_entries", cardinality(prepared.projections[p].solve.stored_entries()));
        frames.set(sectors[p + 1], frame);
    }
    frames.set("harmonic_rank", cardinality(metric.size() - prepared.projections[0].frame.cols() -
                                            prepared.projections[1].frame.cols()));
    out.set("frames", frames);
    auto probes = Json::array();
    for (std::size_t h = 0; h < hypotheses.size(); ++h)
        if (hypotheses[h].projection) {
            auto record = Json::object();
            record.set("hypothesis", family.at(h).at("id"));
            record.set("rank", cardinality(hypotheses[h].projection->frame.cols()));
            record.set("factor_entries",
                       cardinality(hypotheses[h].projection->solve.stored_entries()));
            probes.push(record);
        }
    out.set("probe_frames", probes);
    auto execution_record = Json::object();
    execution_record.set("requested_workers", cardinality(report.requested_workers));
    execution_record.set("worker_slots", cardinality(report.worker_slots));
    execution_record.set("completed_tiles", cardinality(report.completed_tiles));
    execution_record.set("projector_preparations", cardinality(2));
    execution_record.set("probe_projector_preparations", cardinality(probe_preparations.size()));
    execution_record.set("probe_field_decompositions", cardinality(probe_readings.size()));
    out.set("execution", execution_record);
    out.set(
        "calibration_scope",
        Json::string(
            "Inclusive conditional tails under the complete declared joint exchange law; source "
            "evidence records the assumptions without proving exchangeability. Scenario "
            "sensitivity keeps the declared measurement weights separate from the reference law."));
    return out;
}
} // namespace gadget::research_detail
