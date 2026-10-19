#include <algorithm>
#include <cmath>
#include <gadget/statistics/counts.hpp>
#include <limits>
#include <numbers>
#include <numeric>
#include <rfc/machine/checked_size.hpp>
#include <rfc/runtime/execution_policy.hpp>

namespace gadget::statistics {
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::invalid_argument(message);
}
void policy(const std::string &name) {
    require(name == "none" || name == "bonferroni" || name == "benjamini_hochberg" ||
                name == "benjamini_yekutieli",
            "unknown probability adjustment policy");
}
double log_add(double a, double b) {
    if (a < b)
        std::swap(a, b);
    if (!std::isfinite(a))
        return a;
    return a + std::log1p(std::exp(b - a));
}
} // namespace
double normal_log_two_sided(double z) {
    require(std::isfinite(z), "normal reference statistic must be finite");
    z = std::abs(z);
    const auto tail = std::erfc(z / std::numbers::sqrt2);
    if (tail > 0)
        return std::log(tail);
    const long double inverse = 1.L / (static_cast<long double>(z) * z);
    long double term = 1, sum = 1;
    for (unsigned k = 1; k < 128; ++k) {
        const auto next = term * (-static_cast<long double>(2 * k - 1)) * inverse;
        if (std::abs(next) > std::abs(term))
            break;
        sum += next;
        term = next;
        if (std::abs(term) < std::numeric_limits<long double>::epsilon() * std::abs(sum))
            break;
    }
    const long double result =
        -0.5L * static_cast<long double>(z) * z - std::log(static_cast<long double>(z)) +
        0.5L * std::log(2.L / std::numbers::pi_v<long double>) + std::log(sum);
    return static_cast<double>(result);
}
double normal_critical(double confidence) {
    require(std::isfinite(confidence) && confidence > 0 && confidence < 1,
            "normal confidence must lie strictly between zero and one in binary64");
    const auto target = std::log1p(-confidence);
    double lower = 0, upper = 1;
    while (normal_log_two_sided(upper) > target)
        upper *= 2;
    for (unsigned i = 0; i < 100; ++i) {
        const auto middle = std::midpoint(lower, upper);
        if (normal_log_two_sided(middle) > target)
            lower = middle;
        else
            upper = middle;
    }
    return std::midpoint(lower, upper);
}
double chi_square_log_survival(double statistic, std::size_t degrees) {
    require(std::isfinite(statistic) && statistic >= 0 && degrees > 0,
            "chi square reference domain");
    if (!statistic)
        return 0;
    const auto x = statistic / 2;
    if (!x)
        return 0;
    const auto log_x = std::log(x);
    double result, shape;
    std::size_t additions;
    if (degrees % 2) {
        result = normal_log_two_sided(std::sqrt(statistic));
        shape = 0.5;
        additions = (degrees - 1) / 2;
    } else {
        result = -x;
        shape = 1;
        additions = degrees / 2 - 1;
    }
    for (std::size_t i = 0; i < additions; ++i, shape += 1)
        result = log_add(result, shape * log_x - x - std::lgamma(shape + 1));
    return std::min(0.0, result);
}
std::vector<double> adjust_log_probabilities(std::span<const double> values,
                                             const std::string &name) {
    policy(name);
    const auto n = values.size();
    std::vector<double> out(values.begin(), values.end());
    for (const auto value : out)
        require(!std::isnan(value) && value <= 0, "invalid log probability");
    if (!n || name == "none")
        return out;
    double multiplier = std::log(static_cast<double>(n));
    if (name == "benjamini_yekutieli") {
        long double harmonic{};
        for (std::size_t i = 1; i <= n; ++i)
            harmonic += 1.L / i;
        multiplier += static_cast<double>(std::log(harmonic));
    }
    if (name == "bonferroni") {
        for (auto &value : out)
            value = std::min(0.0, value + multiplier);
        return out;
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](auto a, auto b) { return values[a] < values[b]; });
    double following = 0;
    for (std::size_t i = n; i > 0; --i) {
        following = std::min(following,
                             values[order[i - 1]] + multiplier - std::log(static_cast<double>(i)));
        out[order[i - 1]] = following;
    }
    return out;
}
std::vector<rfc::numeric::Rational>
adjust_probabilities(std::span<const rfc::numeric::Rational> values, const std::string &name,
                     std::stop_token stop, const std::function<void(std::size_t)> &admit) {
    using R = rfc::numeric::Rational;
    using I = rfc::numeric::Integer;
    policy(name);
    const auto n = values.size();
    std::vector<R> out(values.begin(), values.end());
    for (const auto &value : out)
        require(value >= 0 && value <= 1, "invalid exact probability");
    if (!n || name == "none")
        return out;
    R multiplier(std::to_string(n));
    if (name == "benjamini_yekutieli") {
        R harmonic{0};
        for (std::size_t i = 1; i <= n; ++i) {
            if (stop.stop_requested())
                throw rfc::runtime::ExecutionCancelled{};
            harmonic += R(I{1}, I(std::to_string(i)));
        }
        const auto bits =
            rfc::machine::checked_size_add(mpz_sizeinbase(harmonic.get_num_mpz_t(), 2),
                                           mpz_sizeinbase(harmonic.get_den_mpz_t(), 2));
        if (admit)
            admit(rfc::machine::checked_size_mul(n,
                                                 rfc::machine::checked_size_mul(bits / 8 + 2, 16)));
        multiplier *= harmonic;
    }
    if (name == "bonferroni") {
        for (auto &value : out)
            value = std::min(R{1}, R(value * multiplier));
        return out;
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](auto a, auto b) { return values[a] < values[b]; });
    R following{1};
    for (std::size_t i = n; i > 0; --i) {
        if (stop.stop_requested())
            throw rfc::runtime::ExecutionCancelled{};
        following =
            std::min(following, R(values[order[i - 1]] * multiplier / R(std::to_string(i))));
        out[order[i - 1]] = following;
    }
    return out;
}
NormalEffectPosterior shrink_normal_effect(double effect, double variance, double prior_mean,
                                           double prior_variance, double confidence) {
    require(std::isfinite(effect) && std::isfinite(variance) && variance > 0 &&
                std::isfinite(prior_mean) && std::isfinite(prior_variance) && prior_variance > 0,
            "normal effect prior/likelihood domain");
    const long double ratio = static_cast<long double>(variance) / prior_variance;
    const auto weight = 1.L / (1 + ratio);
    const auto mean = static_cast<double>(weight * effect + (1 - weight) * prior_mean);
    const auto posterior_variance = static_cast<double>(weight * variance);
    const auto width = normal_critical(confidence) * std::sqrt(posterior_variance);
    require(std::isfinite(mean) && posterior_variance > 0 && std::isfinite(width) &&
                std::isfinite(mean - width) && std::isfinite(mean + width),
            "normal effect posterior exceeds binary64 support");
    return {mean, posterior_variance, mean - width, mean + width};
}
DispersionFit
profile_count_dispersion(const PreparedCountDesign &design, std::size_t features,
                         const std::function<std::vector<double>(std::size_t)> &counts,
                         std::span<const double> offsets, const DispersionControl &dispersion,
                         const CountFitControl &control, CountWork &work, std::stop_token stop) {
    require(features > 0 && static_cast<bool>(counts),
            "dispersion profile requires reference features");
    require(std::isfinite(dispersion.minimum) && dispersion.minimum > 0 &&
                std::isfinite(dispersion.maximum) && dispersion.maximum > dispersion.minimum &&
                std::isfinite(dispersion.log_tolerance) && dispersion.log_tolerance > 0 &&
                dispersion.grid_points >= 3 && dispersion.iterations > 0,
            "dispersion profile domain");
    if (dispersion.prior)
        require(std::isfinite(dispersion.prior->mean) &&
                    std::isfinite(dispersion.prior->variance) && dispersion.prior->variance > 0,
                "log dispersion prior domain");
    DispersionFit result;
    double best = -std::numeric_limits<double>::infinity(), best_log{};
    bool failed{};
    const auto evaluate = [&](double log_alpha) {
        work.charge(1, stop);
        ++result.evaluations;
        long double objective{};
        for (std::size_t feature = 0; feature < features; ++feature) {
            const auto response = counts(feature);
            const auto fit =
                design.fit(response, offsets, std::exp(log_alpha), {}, control, work, stop);
            if (!fit.converged) {
                failed = true;
                return -std::numeric_limits<double>::infinity();
            }
            objective += fit.log_likelihood;
        }
        if (dispersion.prior) {
            const auto delta = static_cast<long double>(log_alpha) - dispersion.prior->mean;
            objective -=
                0.5L *
                (delta * delta / dispersion.prior->variance +
                 std::log(2.L * std::numbers::pi_v<long double> * dispersion.prior->variance));
        }
        const auto value = static_cast<double>(objective);
        if (!std::isfinite(value)) {
            failed = true;
            return -std::numeric_limits<double>::infinity();
        }
        if (value > best) {
            best = value;
            best_log = log_alpha;
        }
        return value;
    };
    const auto lower = std::log(dispersion.minimum), upper = std::log(dispersion.maximum);
    std::size_t selected{};
    for (std::size_t i = 0; i < dispersion.grid_points; ++i) {
        const auto point = std::lerp(
            lower, upper, static_cast<double>(i) / static_cast<double>(dispersion.grid_points - 1));
        result.log_grid.push_back(point);
        result.objectives.push_back(evaluate(point));
        if (result.objectives.back() > result.objectives[selected])
            selected = i;
    }
    if (!std::isfinite(best))
        return result;
    double a = result.log_grid[selected ? selected - 1 : 0],
           b = result.log_grid[std::min(selected + 1, dispersion.grid_points - 1)];
    constexpr double fraction = 0.6180339887498948482;
    double c = b - fraction * (b - a), d = a + fraction * (b - a), fc = evaluate(c),
           fd = evaluate(d);
    for (std::size_t i = 0; i < dispersion.iterations && b - a > dispersion.log_tolerance; ++i) {
        if (fc > fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - fraction * (b - a);
            fc = evaluate(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + fraction * (b - a);
            fd = evaluate(d);
        }
    }
    result.dispersion = std::exp(best_log);
    result.objective = best;
    result.lower_boundary = best_log - lower <= dispersion.log_tolerance;
    result.upper_boundary = upper - best_log <= dispersion.log_tolerance;
    result.converged = !failed && b - a <= dispersion.log_tolerance;
    return result;
}
} // namespace gadget::statistics
