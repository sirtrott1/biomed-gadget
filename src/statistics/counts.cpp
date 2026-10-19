#include "count_support.hpp"
#include <algorithm>
#include <cmath>
#include <gadget/statistics/counts.hpp>
#include <limits>
#include <numbers>
#include <numeric>
#include <rfc/approximation/dense.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/runtime/execution_policy.hpp>

namespace gadget::statistics {
namespace {
using Matrix = rfc::approximation::DenseMatrix;
void require(bool value, const char *message) {
    if (!value)
        throw std::invalid_argument(message);
}
double finite(long double value) {
    const auto converted = static_cast<double>(value);
    if (!std::isfinite(converted))
        throw std::runtime_error("count model exceeds binary64 numerical range");
    return converted;
}
long double stirling_correction(long double value) {
    const auto inverse = 1 / value, squared = inverse * inverse;
    return inverse *
           (1.L / 12 + squared * (-1.L / 360 + squared * (1.L / 1260 + squared * (-1.L / 1680))));
}
long double log_remainder(long double value) {
    if (std::abs(value) >= .125L)
        return std::log1p(value) - value;
    long double power = -value * value, result = power / 2;
    for (unsigned k = 3; k < 256; ++k) {
        power *= -value;
        const auto term = power / k;
        result += term;
        if (std::abs(term) <= std::numeric_limits<long double>::epsilon() * std::abs(result))
            break;
    }
    return result;
}
long double nb_rising_log(double count, double dispersion) {
    const long double y = count, r = 1.L / dispersion;
    if (y < 8) {
        long double out{};
        for (unsigned k = 0; k < static_cast<unsigned>(count); ++k)
            out += std::log1p(static_cast<long double>(dispersion) * k);
        return out;
    }
    if (r < 32)
        return std::lgamma(r + y) - std::lgamma(r) - y * std::log(r);
    const auto t = y / r;
    long double main;
    if (t < 0.125L) {
        long double power = t * t, series = power / 2;
        for (unsigned k = 3; k < 128; ++k) {
            power *= -t;
            const auto term = power / (static_cast<long double>(k) * (k - 1));
            series += term;
            if (std::abs(term) <= std::numeric_limits<long double>::epsilon() * std::abs(series))
                break;
        }
        main = r * series - 0.5L * std::log1p(t);
    } else
        main = (r + y - 0.5L) * std::log1p(t) - y;
    return main + stirling_correction(r + y) - stirling_correction(r);
}
std::vector<double> solve_information(const Matrix &matrix, std::span<const double> rhs,
                                      double tolerance, CountWork &work, std::stop_token stop) {
    const auto p = matrix.rows();
    work.charge(rfc::machine::checked_size_mul(rfc::machine::checked_size_mul(p, p), p), stop);
    Matrix normalized(p, p);
    std::vector<double> diagonal(p), scaled(p);
    for (std::size_t j = 0; j < p; ++j) {
        if (!(matrix(j, j) > 0) || !std::isfinite(matrix(j, j)))
            throw std::runtime_error("count information has a zero direction");
        diagonal[j] = std::sqrt(matrix(j, j));
        scaled[j] = rhs[j] / diagonal[j];
        if (!std::isfinite(scaled[j]))
            throw std::runtime_error("count information scaling is not finite");
    }
    for (std::size_t j = 0; j < p; ++j)
        for (std::size_t k = 0; k < p; ++k)
            normalized(j, k) =
                finite(static_cast<long double>(matrix(j, k)) / diagonal[j] / diagonal[k]);
    auto result = normalized.solve(scaled);
    const auto reconstructed = normalized.apply(result);
    double error{}, magnitude = 1;
    for (std::size_t j = 0; j < p; ++j) {
        error = std::max(error, std::abs(reconstructed[j] - scaled[j]));
        magnitude = std::max(magnitude, std::abs(scaled[j]));
        result[j] = finite(static_cast<long double>(result[j]) / diagonal[j]);
    }
    if (!std::isfinite(error) || error > tolerance * magnitude)
        throw std::runtime_error("count information solve residual exceeds tolerance");
    return result;
}
} // namespace
CountWork::CountWork(std::size_t limit) : limit_(limit) {
    require(limit > 0, "count numerical work budget must be positive");
}
void CountWork::charge(std::size_t work, std::stop_token stop) {
    if (stop.stop_requested())
        throw rfc::runtime::ExecutionCancelled{};
    auto used = used_.load(std::memory_order_relaxed);
    for (;;) {
        if (work > limit_ - used)
            throw CountWorkLimit{};
        if (used_.compare_exchange_weak(used, used + work, std::memory_order_relaxed))
            return;
    }
}
double count_log_mass(double count, double mean, double dispersion) {
    require(std::isfinite(count) && count >= 0 && std::floor(count) == count,
            "count likelihood requires finite integer observations");
    require(std::isfinite(mean) && mean > 0 && std::isfinite(dispersion) && dispersion >= 0,
            "count likelihood mean/dispersion domain");
    const long double y = count, mu = mean;
    if (y >= 32) {
        const auto log_two_pi = std::log(2.L * std::numbers::pi_v<long double>);
        if (!dispersion) {
            const auto delta = (mu - y) / y;
            const auto deviance = delta > -.9L ? y * log_remainder(delta)
                                               : y * (std::log(mu) - std::log(y)) - (mu - y);
            return finite(deviance - .5L * (log_two_pi + std::log(y)) - stirling_correction(y));
        }
        const long double r = 1.L / dispersion;
        if (r >= 32) {
            const auto delta = mu - y, total = r + y;
            long double deviance;
            if (delta / y > -.9L && delta / total > -.9L)
                deviance = y * log_remainder(delta / y) - total * log_remainder(delta / total);
            else
                deviance =
                    y * (std::log(mu) - std::log(y)) - total * (std::log(r + mu) - std::log(total));
            return finite(
                deviance + .5L * (std::log(r) - std::log(total) - std::log(y) - log_two_pi) +
                stirling_correction(total) - stirling_correction(r) - stirling_correction(y));
        }
        const auto gamma_ratio = (r - 1) * std::log(y) + (y + r - .5L) * std::log1p(r / y) - r +
                                 stirling_correction(y + r) - stirling_correction(y);
        return finite(gamma_ratio - std::lgamma(r) + r * (std::log(r) - std::log(r + mu)) -
                      y * std::log1p(r / mu));
    }
    if (!dispersion)
        return finite(y * std::log(mu) - mu - std::lgamma(y + 1));
    const long double a = dispersion;
    return finite(y * std::log(mu) - std::lgamma(y + 1) + nb_rising_log(count, dispersion) -
                  (y + 1 / a) * std::log1p(a * mu));
}
struct PreparedCountDesign::State {
    std::size_t n, p, exact_rank;
    Matrix x, gram;
    rfc::numeric::ExactMatrix exact;
    std::vector<double> scales;
    explicit State(const std::vector<rfc::numeric::RationalVector> &rows)
        : n(rows.size()), p(n ? rows.front().size() : 0), exact_rank(0), x(), gram(), scales(p) {
        require(n > 0, "count design requires observations");
        (void)rfc::machine::checked_size_mul(rfc::machine::checked_size_mul(n, p), sizeof(double));
        (void)rfc::machine::checked_size_mul(rfc::machine::checked_size_mul(p, p), sizeof(double));
        x = Matrix(n, p);
        gram = Matrix(p, p);
        exact = rfc::numeric::ExactMatrix(n, p);
        for (std::size_t i = 0; i < n; ++i) {
            require(rows[i].size() == p, "count design rows have different shapes");
            for (std::size_t j = 0; j < p; ++j) {
                exact.set(i, j, rows[i][j]);
                const auto value = rows[i][j].get_d();
                require(std::isfinite(value), "count design exceeds binary64 input range");
                x(i, j) = value;
                scales[j] = std::hypot(scales[j], value);
            }
        }
        exact_rank = exact.rank();
        for (std::size_t j = 0; j < p; ++j) {
            if (scales[j] == 0)
                scales[j] = 1;
            require(std::isfinite(scales[j]), "count design column norm exceeds binary64 range");
            for (std::size_t i = 0; i < n; ++i)
                x(i, j) /= scales[j];
        }
        for (std::size_t j = 0; j < p; ++j)
            for (std::size_t k = 0; k <= j; ++k) {
                long double sum{};
                for (std::size_t i = 0; i < n; ++i)
                    sum += static_cast<long double>(x(i, j)) * x(i, k);
                gram(j, k) = gram(k, j) = finite(sum);
            }
    }
};
PreparedCountDesign::PreparedCountDesign(const std::vector<rfc::numeric::RationalVector> &rows)
    : state_(std::make_shared<State>(rows)) {}
std::size_t PreparedCountDesign::samples() const noexcept { return state_->n; }
std::size_t PreparedCountDesign::columns() const noexcept { return state_->p; }
std::size_t PreparedCountDesign::rank() const noexcept { return state_->exact_rank; }
std::string PreparedCountDesign::support_status(std::span<const double> counts, CountWork &work,
                                                std::stop_token stop) const {
    require(counts.size() == state_->n, "count support shape differs");
    work.charge(counts.size(), stop);
    for (const auto value : counts)
        require(std::isfinite(value) && value >= 0 && std::floor(value) == value,
                "count support requires nonnegative integers");
    if (state_->exact_rank != state_->p)
        return "rank_deficient";
    if (finite_count_coefficients(state_->exact, counts, work, stop))
        return "finite_coefficients";
    return std::all_of(counts.begin(), counts.end(), [](double value) { return value == 0; })
               ? "all_zero_counts"
               : "boundary_likelihood";
}
CountFit PreparedCountDesign::fit(std::span<const double> counts, std::span<const double> offsets,
                                  double alpha, std::span<const double> contrast,
                                  const CountFitControl &control, CountWork &work,
                                  std::stop_token stop, std::span<const double> initial) const {
    const auto &s = *state_;
    const auto n = s.n, p = s.p;
    require(counts.size() == n && offsets.size() == n && (contrast.empty() || contrast.size() == p),
            "count fit input shape differs");
    require(std::isfinite(alpha) && alpha >= 0, "count dispersion domain");
    require(control.iterations > 0 && control.line_search_steps > 0 &&
                control.score_tolerance > 0 && std::isfinite(control.score_tolerance) &&
                control.linear_tolerance > 0 && std::isfinite(control.linear_tolerance),
            "count fit control domain");
    require(initial.empty() || initial.size() == p, "count initial coefficient shape differs");
    for (std::size_t i = 0; i < n; ++i)
        require(std::isfinite(counts[i]) && counts[i] >= 0 && std::floor(counts[i]) == counts[i] &&
                    std::isfinite(offsets[i]),
                "count observations/offsets outside numerical support");
    for (const auto value : contrast)
        require(std::isfinite(value), "count contrast is not finite");
    CountFit fit;
    fit.status = support_status(counts, work, stop);
    if (fit.status != "finite_coefficients")
        return fit;
    std::vector<double> beta(p), gradient(p), means(n);
    Matrix information(p, p);
    if (!initial.empty())
        for (std::size_t j = 0; j < p; ++j)
            beta[j] = finite(static_cast<long double>(initial[j]) * s.scales[j]);
    else if (p) {
        std::vector<double> rhs(p);
        work.charge(rfc::machine::checked_size_mul(n, p), stop);
        for (std::size_t j = 0; j < p; ++j) {
            long double sum{};
            for (std::size_t i = 0; i < n; ++i)
                sum +=
                    static_cast<long double>(s.x(i, j)) * (std::log(counts[i] + 0.5) - offsets[i]);
            rhs[j] = finite(sum);
        }
        try {
            beta = solve_information(s.gram, rhs, control.linear_tolerance, work, stop);
        } catch (const CountWorkLimit &) {
            throw;
        } catch (const rfc::runtime::ExecutionCancelled &) {
            throw;
        } catch (const std::runtime_error &) {
            fit.status = "numerically_singular_design";
            return fit;
        }
    }
    const auto evaluate = [&](std::span<const double> coefficients, bool derivatives) {
        work.charge(rfc::machine::checked_size_mul(
                        n, rfc::machine::checked_size_add(
                               1, derivatives ? rfc::machine::checked_size_mul(p, p) : p)),
                    stop);
        std::vector<long double> scores(derivatives ? p : 0),
            entries(derivatives ? rfc::machine::checked_size_mul(p, p) : 0);
        long double likelihood{}, pearson{};
        for (std::size_t i = 0; i < n; ++i) {
            if ((i & 255) == 0)
                work.charge(0, stop);
            long double eta = offsets[i];
            for (std::size_t j = 0; j < p; ++j)
                eta += static_cast<long double>(s.x(i, j)) * coefficients[j];
            const auto mu = static_cast<double>(std::exp(eta));
            if (!(mu > 0) || !std::isfinite(mu))
                throw std::runtime_error("count log mean exceeds binary64 support");
            if (derivatives)
                means[i] = mu;
            likelihood += count_log_mass(counts[i], mu, alpha);
            if (derivatives) {
                const auto denominator = 1 + static_cast<long double>(alpha) * mu;
                const auto residual = (static_cast<long double>(counts[i]) - mu) / denominator;
                const auto weight = static_cast<long double>(mu) / denominator;
                pearson += (static_cast<long double>(counts[i]) - mu) * residual / mu;
                for (std::size_t j = 0; j < p; ++j) {
                    scores[j] += s.x(i, j) * residual;
                    for (std::size_t k = 0; k <= j; ++k)
                        entries[j * p + k] += s.x(i, j) * s.x(i, k) * weight;
                }
            }
        }
        if (derivatives) {
            for (std::size_t j = 0; j < p; ++j) {
                gradient[j] = finite(scores[j]);
                for (std::size_t k = 0; k <= j; ++k)
                    information(j, k) = information(k, j) = finite(entries[j * p + k]);
            }
            fit.pearson_sum_squares = finite(pearson);
        }
        return finite(likelihood);
    };
    try {
        for (std::size_t iteration = 0; iteration <= control.iterations; ++iteration) {
            const auto likelihood = evaluate(beta, true);
            fit.log_likelihood = likelihood;
            fit.iterations = iteration;
            double error{};
            for (std::size_t j = 0; j < p; ++j)
                error = std::max(error, std::abs(gradient[j]) /
                                            std::max(1.0, std::sqrt(information(j, j))));
            fit.relative_score = error;
            if (error <= control.score_tolerance) {
                fit.converged = true;
                fit.status = "converged";
                break;
            }
            if (iteration == control.iterations) {
                fit.status = "iteration_limit";
                break;
            }
            const auto step =
                solve_information(information, gradient, control.linear_tolerance, work, stop);
            bool accepted{};
            double length = 1;
            for (std::size_t trial = 0; trial < control.line_search_steps; ++trial, length *= 0.5) {
                std::vector<double> candidate(p);
                for (std::size_t j = 0; j < p; ++j)
                    candidate[j] = finite(static_cast<long double>(beta[j]) + length * step[j]);
                try {
                    if (evaluate(candidate, false) >= likelihood) {
                        beta = std::move(candidate);
                        accepted = true;
                        break;
                    }
                } catch (const CountWorkLimit &) {
                    throw;
                } catch (const rfc::runtime::ExecutionCancelled &) {
                    throw;
                } catch (const std::runtime_error &) {
                }
            }
            if (!accepted) {
                fit.status = "line_search_failed";
                break;
            }
        }
        fit.coefficients.resize(p);
        fit.score.resize(p);
        for (std::size_t j = 0; j < p; ++j) {
            fit.coefficients[j] = finite(static_cast<long double>(beta[j]) / s.scales[j]);
            fit.score[j] = finite(static_cast<long double>(gradient[j]) * s.scales[j]);
        }
        fit.means = std::move(means);
        if (fit.converged && !contrast.empty()) {
            std::vector<double> covector(p);
            long double effect{};
            for (std::size_t j = 0; j < p; ++j) {
                covector[j] = contrast[j] / s.scales[j];
                effect += static_cast<long double>(contrast[j]) * fit.coefficients[j];
            }
            const auto response =
                solve_information(information, covector, control.linear_tolerance, work, stop);
            long double variance{};
            for (std::size_t j = 0; j < p; ++j)
                variance += static_cast<long double>(covector[j]) * response[j];
            fit.effect = finite(effect);
            fit.variance = finite(variance);
            if (!(*fit.variance > 0)) {
                fit.converged = false;
                fit.status = "nonpositive_contrast_variance";
            }
        }
    } catch (const CountWorkLimit &) {
        throw;
    } catch (const rfc::runtime::ExecutionCancelled &) {
        throw;
    } catch (const std::runtime_error &error) {
        fit.converged = false;
        fit.status = error.what();
    }
    return fit;
}
} // namespace gadget::statistics
