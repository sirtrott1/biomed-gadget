#include "count_analysis_detail.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <numeric>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::statistics {
namespace {
using namespace count_detail;
using namespace scientific_detail;
using rfc::machine::checked_size_add;
using rfc::machine::checked_size_mul;
using Window = rfc::database::NativeTensorWindow;
using Value = ResultValue;
using Integer = rfc::numeric::Integer;
void check(std::stop_token stop) {
    if (stop.stop_requested())
        throw rfc::runtime::ExecutionCancelled{};
}
double finite(long double value) {
    const double result = value;
    if (!std::isfinite(result))
        throw std::invalid_argument("count parameter exceeds binary64 support");
    return result;
}
Vector exact(std::span<const double> values) {
    Vector out;
    out.reserve(values.size());
    for (const auto value : values) {
        require(std::isfinite(value), "nonfinite count output");
        out.emplace_back(value);
    }
    return out;
}
Value labels(const std::vector<std::string> &names) {
    auto out = Value::array();
    for (const auto &name : names)
        out.push(Value::string(name));
    return out;
}
Value scalar(ResultFields &fields, double value) {
    return std::isfinite(value) ? fields.scalar(Rational(value)) : Value{};
}
Value profile_view(ResultFields &fields, const DispersionFit &fit) {
    auto out = Value::object();
    out.set("value", scalar(fields, fit.dispersion));
    out.set("objective", scalar(fields, fit.objective));
    out.set("converged", Value::boolean(fit.converged));
    out.set("lower_boundary", Value::boolean(fit.lower_boundary));
    out.set("upper_boundary", Value::boolean(fit.upper_boundary));
    out.set("evaluations", Value::count(fit.evaluations));
    auto grid = Value::array();
    for (std::size_t j = 0; j < fit.log_grid.size(); ++j) {
        auto item = Value::object();
        item.set("log_dispersion", scalar(fields, fit.log_grid[j]));
        item.set("objective", scalar(fields, fit.objectives[j]));
        grid.push(std::move(item));
    }
    out.set("grid", std::move(grid));
    return out;
}
} // namespace
struct PreparedCountAnalysis::State {
    Configuration c;
    NativeRecord source, request;
    std::string digest;
    std::shared_ptr<CpuResources> resources;
    std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> memory;
    std::map<std::string, Window> fields;
    std::map<std::string, std::uint64_t> field_keys;
    std::optional<PreparedCountDesign> design, reduced;
    std::vector<double> factors, offsets;
    std::vector<std::string> eligible_references;
    std::size_t preparation_work{}, working_bytes{}, source_size{}, source_samples{};
    std::uint64_t source_carrier{};
    mutable std::mutex mutex;
    State(Configuration config, NativeRecord input, NativeRecord query, std::string pin,
          std::shared_ptr<CpuResources> cpu)
        : c(std::move(config)), source(std::move(input)), request(std::move(query)),
          digest(std::move(pin)), resources(std::move(cpu)) {}
    Observation observe(const std::string &field, CountWork &work, std::stop_token stop) const {
        work.charge(c.samples.size(), stop);
        Observation out;
        out.values.resize(c.samples.size());
        const auto &window = fields.at(field);
        for (std::size_t i = 0; i < c.samples.size(); ++i) {
            if ((i & 255) == 0)
                check(stop);
            const auto value = window.coordinate(c.samples[i].coordinate).state();
            if (!value) {
                ++out.missing;
                continue;
            }
            try {
                out.values[i] = exact_count(*value);
            } catch (const std::invalid_argument &) {
                ++out.invalid;
            }
        }
        return out;
    }
    std::vector<double> complete(const std::string &field, CountWork &work,
                                 std::stop_token stop) const {
        auto observation = observe(field, work, stop);
        require(!observation.missing && !observation.invalid,
                "normalization/dispersion reference has missing or invalid counts");
        return std::move(observation.values);
    }
    void normalize(CountWork &work, std::stop_token stop) {
        const auto n = c.samples.size();
        factors = c.factors;
        offsets.resize(n);
        std::vector<double> logs(n);
        if (c.normalization == "median_ratio") {
            std::vector<double> geometric(c.normalization_fields.size());
            std::vector<std::size_t> eligible;
            for (std::size_t f = 0; f < c.normalization_fields.size(); ++f) {
                const auto values = complete(c.normalization_fields[f], work, stop);
                if (std::any_of(values.begin(), values.end(), [](double v) { return v == 0; }))
                    continue;
                long double sum{};
                for (const auto v : values)
                    sum += std::log(v);
                geometric[f] = finite(sum / n);
                eligible.push_back(f);
                eligible_references.push_back(c.normalization_fields[f]);
            }
            require(!eligible.empty(),
                    "median ratio requires a reference positive in every selected sample");
            // Pins are retained, while each sample uses only one reference buffer.
            std::vector<double> ratios;
            ratios.reserve(eligible.size());
            for (std::size_t i = 0; i < n; ++i) {
                ratios.clear();
                work.charge(eligible.size(), stop);
                for (const auto f : eligible) {
                    const auto value = fields.at(c.normalization_fields[f])
                                           .coordinate(c.samples[i].coordinate)
                                           .state();
                    ratios.push_back(std::log(exact_count(*value)) - geometric[f]);
                }
                const auto middle = ratios.size() / 2;
                std::nth_element(ratios.begin(), ratios.begin() + middle, ratios.end());
                logs[i] = ratios[middle];
                if (ratios.size() % 2 == 0) {
                    const auto lower = *std::max_element(ratios.begin(), ratios.begin() + middle);
                    logs[i] += std::log1p(std::exp(lower - logs[i])) - std::log(2.0);
                }
            }
        } else if (!factors.empty())
            for (std::size_t i = 0; i < n; ++i)
                logs[i] = std::log(factors[i]);
        if (c.center) {
            long double sum{};
            for (const auto value : logs)
                sum += value;
            const auto center = sum / n;
            for (auto &value : logs)
                value = finite(value - center);
        }
        factors.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            factors[i] = std::exp(logs[i]);
            require(std::isfinite(factors[i]) && factors[i] > 0,
                    "normalization factor outside binary64 support");
            offsets[i] = finite(static_cast<long double>(logs[i]) +
                                (c.exposure.empty() ? 0 : std::log(c.exposure[i])));
        }
    }
};
PreparedCountAnalysis PreparedCountAnalysis::prepare(const ScientificStudy &study,
                                                     const NativeRecord &request,
                                                     std::shared_ptr<CpuResources> resources,
                                                     std::stop_token stop) {
    check(stop);
    if (!resources)
        resources = CpuResources::local();
    (void)resources->readmit();
    // Native limits admit the request render and every retained output first.
    const auto raw_bytes = request.at("limits").at("declared_buffer_bytes").integer();
    require(raw_bytes > 0 &&
                raw_bytes <= Integer(std::to_string(std::numeric_limits<std::size_t>::max())),
            "count byte admission outside addressability");
    auto memory = resources->coordinator().reserve(1, static_cast<std::size_t>(raw_bytes.get_ui()));
    auto c = admit_count_request(request);
    require(request.at("source_digest").text() == study.digest(),
            "count source digest differs from actual study");
    auto state = std::make_shared<State>(std::move(c), study.native_publication(), request,
                                         study.digest(), resources);
    auto &s = *state;
    const auto n = s.c.samples.size(), p = s.c.columns.size(), f = s.c.features.size();
    require(s.c.workers <= resources->coordinator().capacity(),
            "count requested workers exceed admitted CPU capacity");
    const auto add = [&](std::size_t bytes) {
        s.working_bytes = checked_size_add(s.working_bytes, bytes);
        require(s.working_bytes <= s.c.bytes, "count working and output byte admission");
    };
    add(checked_size_mul(s.c.request_bytes, 32));
    add(checked_size_mul(checked_size_mul(n, p), 1024));
    add(checked_size_mul(checked_size_mul(p, p), checked_size_mul(s.c.workers, 512)));
    add(checked_size_mul(checked_size_mul(n, s.c.workers), 512));
    add(checked_size_mul(checked_size_mul(f, checked_size_add(p, 32)), 1024));
    if (s.c.fitted)
        add(checked_size_mul(checked_size_mul(f, n), 1024));
    add(checked_size_mul(checked_size_add(s.c.normalization_fields.size(), s.c.pool_fields.size()),
                         1024));
    if (s.c.dispersion_method != "fixed")
        add(checked_size_mul(checked_size_mul(s.c.profile.grid_points, checked_size_add(f, 1)),
                             1024));
    rfc::storage::stream_native_series(
        *s.source.publication_history(), [&](std::span<const std::byte> bytes) {
            check(stop);
            s.source_size = checked_size_add(s.source_size, bytes.size());
            require(s.source_size <= s.c.bytes / 16, "count source publication byte admission");
        });
    add(checked_size_mul(s.source_size, 16));
    std::size_t design_bytes{}, max_limb_words{1};
    for (const auto &row : s.c.rows)
        for (const auto &value : row) {
            const auto bytes = rfc::machine::checked_size_add(mpz_size(value.get_num_mpz_t()),
                                                              mpz_size(value.get_den_mpz_t()));
            max_limb_words = std::max(max_limb_words, bytes);
            design_bytes =
                checked_size_add(design_bytes, checked_size_mul(bytes, sizeof(mp_limb_t)));
        }
    add(checked_size_mul(checked_size_mul(design_bytes, checked_size_mul(p, p)), 16));
    CountWork work(s.c.work);
    work.charge(checked_size_add(checked_size_mul(n, checked_size_mul(p, p)),
                                 checked_size_mul(p, checked_size_mul(p, p))),
                stop);
    work.charge(checked_size_mul(checked_size_mul(n, checked_size_mul(p, p)),
                                 checked_size_mul(max_limb_words, max_limb_words)),
                stop);
    s.design.emplace(s.c.rows);
    require(s.design->rank() == p, "count design is exactly rank deficient");
    if (s.c.inference == "likelihood_ratio") {
        std::vector<rfc::numeric::RationalVector> reduced_rows;
        for (const auto &row : s.c.rows) {
            rfc::numeric::RationalVector values;
            for (auto j : s.c.reduced)
                values.push_back(row[j]);
            reduced_rows.push_back(std::move(values));
        }
        s.reduced.emplace(reduced_rows);
    }
    std::set<std::string> needed;
    for (const auto &feature : s.c.features)
        needed.insert(feature.field);
    needed.insert(s.c.normalization_fields.begin(), s.c.normalization_fields.end());
    needed.insert(s.c.pool_fields.begin(), s.c.pool_fields.end());
    add(checked_size_mul(needed.size(), 4096));
    std::optional<rfc::storage::PinnedNativeAccession> first;
    const auto cancellation = std::make_shared<std::atomic_bool>(false);
    std::stop_callback on_stop(stop, [&] { cancellation->store(true, std::memory_order_relaxed); });
    const auto base_lease = source_execution_lease("statistics-count-selection");
    const auto lease = rfc::runtime::ExecutionLease(base_lease.policy(), cancellation);
    for (const auto &id : needed) {
        check(stop);
        const auto declaration = study.field_declaration(id);
        const auto pin = study.field_accession(id);
        require(declaration.at("quantity").text() == s.c.quantity &&
                    declaration.at("units").text() == s.c.units,
                "count quantity or units differ from source declaration");
        require(declaration.at("value_type").text() == "quantitative",
                "count model requires quantitative fields");
        const auto carrier = pin.carrier();
        require(carrier->kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                    carrier->grade == rfc::core::Grade{0},
                "count field requires its original grade zero carrier");
        require(pin.tensor()->domain().rank() == 1, "count scalar field requires one axis");
        if (first)
            require(pin.same_execution_authority(*first) && pin.sequence() == first->sequence() &&
                        pin.source_namespace() == first->source_namespace() &&
                        carrier->key == first->carrier()->key,
                    "count fields do not share their actual carrier and cutoff");
        else {
            first = pin;
            s.source_carrier = carrier->key;
            s.source_samples = pin.tensor()->domain().logical_size();
        }
        s.field_keys.emplace(id, pin.tensor_declaration()->key);
        s.fields.emplace(id, Window::prepare(pin, {0}, {s.source_samples},
                                             {pin.source_namespace(), pin.sequence()}, lease));
    }
    for (auto &sample : s.c.samples) {
        require(sample.coordinate < s.source_samples,
                "count sample coordinate outside original carrier");
        sample.native_id = study.cell_identifier({0}, sample.coordinate);
    }
    s.normalize(work, stop);
    s.preparation_work = work.used();
    resources->coordinator().resize(memory, 0, s.c.bytes);
    s.memory.emplace(std::move(memory));
    check(stop);
    return PreparedCountAnalysis(std::move(state));
}
NativeRecord PreparedCountAnalysis::execute(std::stop_token stop) const {
    auto &s = *state_;
    std::lock_guard lock(s.mutex);
    check(stop);
    auto reservation = s.resources->reserve(s.c.workers);
    CountWork work(s.c.work);
    work.charge(s.preparation_work, stop);
    const auto &c = s.c;
    std::optional<DispersionFit> pool;
    double pooled_alpha = c.alpha;
    if (c.dispersion_method == "pooled") {
        pool = profile_count_dispersion(
            *s.design, c.pool_fields.size(),
            [&](std::size_t i) { return s.complete(c.pool_fields[i], work, stop); }, s.offsets,
            c.profile, c.control, work, stop);
        require(pool->converged, "pooled dispersion profile did not converge");
        pooled_alpha = pool->dispersion;
    }
    std::vector<FeatureResult> results(c.features.size());
    rfc::runtime::CpuTileExecutor executor(
        {c.workers, s.resources->coordinator().capacity(), c.features.size()});
    (void)executor.run_tiles(c.features.size(), [&](std::size_t i, std::size_t,
                                                    const rfc::runtime::CpuTileStop &tile) {
        tile.check();
        auto &result = results[i];
        auto observation = s.observe(c.features[i].field, work, stop);
        result.missing = observation.missing;
        result.invalid = observation.invalid;
        result.alpha = pooled_alpha;
        if (result.missing || result.invalid) {
            result.status = result.invalid ? "invalid_counts" : "missing_counts";
            return;
        }
        const auto support = s.design->support_status(observation.values, work, stop);
        if (support != "finite_coefficients") {
            result.status = support;
            return;
        }
        const auto response = [&](std::size_t) { return observation.values; };
        if (c.dispersion_method == "profile") {
            result.profile = profile_count_dispersion(*s.design, 1, response, s.offsets, c.profile,
                                                      c.control, work, stop);
            result.alpha = result.profile->dispersion;
            if (!result.profile->converged) {
                result.status = "dispersion_not_converged";
                return;
            }
        }
        result.fit = s.design->fit(observation.values, s.offsets, result.alpha, c.contrast,
                                   c.control, work, stop);
        result.status = result.fit.status;
        if (!c.fitted)
            result.fit.means = std::vector<double>{};
        if (!result.fit.converged)
            return;
        const auto effect = *result.fit.effect, variance = *result.fit.variance,
                   se = std::sqrt(variance);
        const auto critical = normal_critical(c.confidence);
        result.low = finite(static_cast<long double>(effect) - critical * se);
        result.high = finite(static_cast<long double>(effect) + critical * se);
        result.z = effect / se;
        if (c.inference == "wald")
            result.log_p = normal_log_two_sided(*result.z);
        else {
            double alpha = result.alpha;
            if (c.dispersion_method == "profile") {
                result.reduced_profile = profile_count_dispersion(
                    *s.reduced, 1, response, s.offsets, c.profile, c.control, work, stop);
                if (!result.reduced_profile->converged) {
                    result.status = "reduced_dispersion_not_converged";
                    return;
                }
                alpha = result.reduced_profile->dispersion;
            }
            result.reduced_fit =
                s.reduced->fit(observation.values, s.offsets, alpha, {}, c.control, work, stop);
            result.reduced_fit->means = std::vector<double>{};
            if (!result.reduced_fit->converged) {
                result.status = "reduced_model_not_converged";
                return;
            }
            const auto difference = result.fit.log_likelihood - result.reduced_fit->log_likelihood;
            const auto tolerance = 1e-8 * (1 + std::abs(result.fit.log_likelihood) +
                                           std::abs(result.reduced_fit->log_likelihood));
            if (difference < -tolerance) {
                result.status = "nested_likelihood_order_failed";
                return;
            }
            result.statistic = 2 * std::max(0.0, difference);
            const auto degrees = c.columns.size() - c.reduced.size();
            work.charge(degrees, stop);
            result.log_p = chi_square_log_survival(*result.statistic, degrees);
        }
        if (c.effect_prior)
            result.posterior = shrink_normal_effect(effect, variance, c.effect_prior->mean,
                                                    c.effect_prior->variance, c.confidence);
        if (!c.fitted)
            result.fit.means = std::vector<double>{};
        if (result.reduced_fit)
            result.reduced_fit->means = std::vector<double>{};
        tile.check();
    });
    std::vector<double> probabilities;
    for (const auto &result : results)
        probabilities.push_back(result.log_p.value_or(0));
    work.charge(checked_size_mul(probabilities.size(),
                                 checked_size_add(std::bit_width(probabilities.size()), 2)),
                stop);
    const auto adjusted = adjust_log_probabilities(probabilities, c.correction);
    for (std::size_t i = 0; i < results.size(); ++i)
        if (results[i].log_p)
            results[i].adjusted_log_p = adjusted[i];
    check(stop);
    ResultFields fields(s.source);
    auto out = Value::object();
    out.set("schema", Value::string("gadget.statistics.count.result.v1"));
    out.set("operation", Value::string("count_model_analysis"));
    out.set("source_digest", Value::string(s.digest));
    out.set("numerical_kind",
            Value::string("binary64 model estimates retained as native rational fields"));
    out.set("family", Value::string(c.family));
    out.set("dispersion_method", Value::string(c.dispersion_method));
    out.set("normalization", Value::string(c.normalization));
    out.set("inference", Value::string(c.inference));
    out.set("correction", Value::string(c.correction));
    out.set("independence",
            Value::string("Samples are declared independent experimental units. Repeated or "
                          "correlated observations require a different model."));
    out.set("inference_scope",
            Value::string("Asymptotic normal or chi square calibration conditional on the declared "
                          "model, offsets and dispersion policy. Shared dispersion and "
                          "normalization uncertainty are not propagated."));
    out.set("multiplicity_scope",
            Value::string(
                "The complete declared feature family includes unavailable tests as probability "
                "one. BH requires valid null tails and independence or suitable positive "
                "dependence. BY permits arbitrary dependence with valid null tails."));
    out.set("parameter_scale", Value::string("natural logarithm"));
    out.set("link",
            Value::string(
                "log mean = design times coefficients + log normalization factor + log exposure"));
    out.set("variance",
            Value::string(c.family == "poisson" ? "mean" : "mean + dispersion times mean squared"));
    auto controls = Value::object();
    controls.set("score_tolerance", scalar(fields, c.control.score_tolerance));
    controls.set("linear_tolerance", scalar(fields, c.control.linear_tolerance));
    controls.set("iterations", Value::count(c.control.iterations));
    controls.set("line_search_steps", Value::count(c.control.line_search_steps));
    controls.set("finite_coefficients", Value::string("exact design recession constraints"));
    out.set("solver_controls", std::move(controls));
    if (c.dispersion_method == "pooled")
        out.set("dispersion_objective",
                Value::string("sum of reference marginal profile log likelihoods with optional log "
                              "dispersion prior"));
    out.set("columns", labels(c.columns));
    const auto observations = fields.axis(c.samples.size()),
               parameters = fields.axis(c.columns.size());
    out.set("normalization_factors", fields.vector(exact(s.factors), observations));
    out.set("offsets", fields.vector(exact(s.offsets), observations));
    out.set("contrast", fields.vector(exact(c.contrast), parameters));
    out.set("normalization_references", labels(s.eligible_references));
    auto samples = Value::array();
    std::vector<std::size_t> coordinates;
    rfc::numeric::ExactMatrix selection(s.source_samples, c.samples.size());
    for (std::size_t i = 0; i < c.samples.size(); ++i) {
        const auto &sample = c.samples[i];
        auto item = Value::object();
        item.set("id", Value::string(sample.id));
        item.set("experimental_unit", Value::string(sample.unit));
        item.set("native_id", Value::string(sample.native_id));
        item.set("coordinate", Value::count(sample.coordinate));
        samples.push(std::move(item));
        coordinates.push_back(sample.coordinate);
        selection.set(sample.coordinate, i, 1);
    }
    out.set("samples", std::move(samples));
    out.set("source_coordinates", fields.indices(coordinates, observations));
    out.set("source_selection_map",
            fields.matrix(selection, observations, fields.source_grade(0, s.source_carrier)));
    rfc::numeric::ExactMatrix design(c.samples.size(), c.columns.size());
    for (std::size_t i = 0; i < c.samples.size(); ++i)
        for (std::size_t j = 0; j < c.columns.size(); ++j)
            design.set(i, j, c.rows[i][j]);
    out.set("design", fields.matrix(design, parameters, observations));
    out.set("design_rank", Value::count(s.design->rank()));
    std::optional<std::uint64_t> reduced_axis;
    if (s.reduced) {
        reduced_axis = fields.axis(c.reduced.size());
        out.set("reduced_columns", fields.indices(c.reduced, *reduced_axis));
        out.set("likelihood_ratio_degrees", Value::count(c.columns.size() - c.reduced.size()));
        rfc::numeric::ExactMatrix reduction(c.columns.size(), c.reduced.size());
        for (std::size_t j = 0; j < c.reduced.size(); ++j)
            reduction.set(c.reduced[j], j, 1);
        out.set("reduced_design_map", fields.matrix(reduction, *reduced_axis, parameters));
    }
    if (pool)
        out.set("pooled_dispersion", profile_view(fields, *pool));
    auto features = Value::array();
    for (std::size_t i = 0; i < results.size(); ++i) {
        check(stop);
        const auto &result = results[i];
        auto item = Value::object();
        item.set("id", Value::string(c.features[i].id));
        item.set("field", Value::string(c.features[i].field));
        item.set("input", fields.stored(s.field_keys.at(c.features[i].field), "vector"));
        // The original declaration and archived study retain the field accession.
        item.set("status", Value::string(result.status));
        item.set("missing", Value::count(result.missing));
        item.set("invalid", Value::count(result.invalid));
        item.set("dispersion", scalar(fields, result.alpha));
        item.set("iterations", Value::count(result.fit.iterations));
        if (!result.fit.coefficients.empty())
            item.set("coefficients", fields.vector(exact(result.fit.coefficients), parameters));
        if (result.fit.converged) {
            item.set("log_likelihood", scalar(fields, result.fit.log_likelihood));
            item.set("relative_score", scalar(fields, result.fit.relative_score));
            item.set("pearson_sum_squares", scalar(fields, result.fit.pearson_sum_squares));
            item.set("score", fields.vector(exact(result.fit.score), parameters));
            item.set("effect", scalar(fields, *result.fit.effect));
            item.set("variance", scalar(fields, *result.fit.variance));
            item.set("confidence_low", scalar(fields, *result.low));
            item.set("confidence_high", scalar(fields, *result.high));
            item.set("wald_z", scalar(fields, *result.z));
            if (c.fitted)
                item.set("fitted_means", fields.vector(exact(result.fit.means), observations));
        }
        if (result.profile)
            item.set("dispersion_profile", profile_view(fields, *result.profile));
        if (result.reduced_profile)
            item.set("reduced_dispersion_profile", profile_view(fields, *result.reduced_profile));
        if (result.reduced_fit && result.reduced_fit->converged) {
            item.set("reduced_log_likelihood", scalar(fields, result.reduced_fit->log_likelihood));
            item.set("reduced_coefficients",
                     fields.vector(exact(result.reduced_fit->coefficients), *reduced_axis));
            item.set("reduced_score",
                     fields.vector(exact(result.reduced_fit->score), *reduced_axis));
            item.set("reduced_relative_score", scalar(fields, result.reduced_fit->relative_score));
            item.set("reduced_iterations", Value::count(result.reduced_fit->iterations));
        }
        if (result.statistic)
            item.set("likelihood_ratio", scalar(fields, *result.statistic));
        for (const auto &[name, value] : std::vector<std::pair<std::string, std::optional<double>>>{
                 {"p_value", result.log_p}, {"adjusted_p_value", result.adjusted_log_p}}) {
            item.set(name, value ? scalar(fields, std::exp(*value)) : Value{});
            item.set("log_" + name, value ? scalar(fields, *value) : Value{});
            item.set(name + "_underflow", Value::boolean(value && std::exp(*value) == 0));
        }
        if (result.posterior) {
            auto posterior = Value::object();
            posterior.set(
                "method",
                Value::string("normal likelihood approximation with declared normal prior"));
            posterior.set("mean", scalar(fields, result.posterior->mean));
            posterior.set("variance", scalar(fields, result.posterior->variance));
            posterior.set("credible_low", scalar(fields, result.posterior->low));
            posterior.set("credible_high", scalar(fields, result.posterior->high));
            item.set("effect_posterior", std::move(posterior));
        }
        features.push(std::move(item));
    }
    out.set("features", std::move(features));
    out.set("tested_family_size", Value::count(results.size()));
    out.set("numerical_work", Value::count(work.used()));
    out.set("declared_buffer_bytes", Value::count(c.bytes));
    out.set("working_byte_estimate", Value::count(s.working_bytes));
    out.set("workers", Value::count(c.workers));
    const auto archive = fields.publication(s.source);
    out.set("source_study_publication", archive.at("publication"));
    check(stop);
    auto result = fields.finish_operation(std::move(out), {}, "", {}, s.request);
    check(stop);
    return result;
}
Json count_result_view(const NativeRecord &result) {
    require(result.at("result").at("view").at("schema").text() ==
                "gadget.statistics.count.result.v1",
            "count result schema");
    return scientific_detail::result_fields_view(result.at("result"));
}
NativeRecord count_model_analysis(const ScientificStudy &study, const NativeRecord &request,
                                  std::shared_ptr<CpuResources> resources, std::stop_token stop) {
    return PreparedCountAnalysis::prepare(study, request, std::move(resources), stop).execute(stop);
}
} // namespace gadget::statistics
