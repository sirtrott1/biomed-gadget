#pragma once
#include <atomic>
#include <functional>
#include <memory>
#include <optional>
#include <rfc/numeric/exact.hpp>
#include <span>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <vector>

namespace gadget::statistics {
class CountWorkLimit : public std::runtime_error {
  public:
    CountWorkLimit() : std::runtime_error("count model numerical work budget exhausted") {}
};
// One execution shares this admission across feature fits and profile evaluations.
class CountWork {
    std::size_t limit_;
    std::atomic_size_t used_{};

  public:
    explicit CountWork(std::size_t limit);
    void charge(std::size_t work, std::stop_token stop = {});
    std::size_t used() const noexcept { return used_.load(std::memory_order_relaxed); }
};
struct CountFitControl {
    std::size_t iterations{200}, line_search_steps{40};
    double score_tolerance{1e-8}, linear_tolerance{1e-8};
};
struct CountFit {
    bool converged{};
    std::string status;
    std::size_t iterations{};
    std::vector<double> coefficients, means, score;
    double log_likelihood{}, relative_score{}, pearson_sum_squares{};
    std::optional<double> effect, variance;
};
// Exact rank is checked once. Binary64 design columns and their Gram matrix are
// retained for repeated feature fits; no observation or feature is removed.
class PreparedCountDesign {
    struct State;
    std::shared_ptr<const State> state_;

  public:
    explicit PreparedCountDesign(const std::vector<rfc::numeric::RationalVector> &rows);
    std::size_t samples() const noexcept;
    std::size_t columns() const noexcept;
    std::size_t rank() const noexcept;
    std::string support_status(std::span<const double> counts, CountWork &work,
                               std::stop_token stop = {}) const;
    CountFit fit(std::span<const double> counts, std::span<const double> offsets, double dispersion,
                 std::span<const double> contrast, const CountFitControl &control, CountWork &work,
                 std::stop_token stop = {}, std::span<const double> initial = {}) const;
};
struct LogDispersionPrior {
    double mean{}, variance{1};
};
struct DispersionControl {
    double minimum{1e-6}, maximum{10}, log_tolerance{1e-5};
    std::size_t grid_points{12}, iterations{100};
    std::optional<LogDispersionPrior> prior;
};
struct DispersionFit {
    double dispersion{}, objective{};
    std::size_t evaluations{};
    bool lower_boundary{}, upper_boundary{}, converged{};
    std::vector<double> log_grid, objectives;
};
// Optimize the sum of the supplied marginal profile likelihoods. The callback
// supplies one complete feature at a time, so reference counts need no dense cache.
DispersionFit
profile_count_dispersion(const PreparedCountDesign &design, std::size_t features,
                         const std::function<std::vector<double>(std::size_t)> &counts,
                         std::span<const double> offsets, const DispersionControl &dispersion,
                         const CountFitControl &control, CountWork &work,
                         std::stop_token stop = {});
double count_log_mass(double count, double mean, double dispersion);
double normal_log_two_sided(double z);
double normal_critical(double confidence);
double chi_square_log_survival(double statistic, std::size_t degrees);
std::vector<double> adjust_log_probabilities(std::span<const double> values,
                                             const std::string &policy);
std::vector<rfc::numeric::Rational>
adjust_probabilities(std::span<const rfc::numeric::Rational> values, const std::string &policy,
                     std::stop_token stop = {}, const std::function<void(std::size_t)> &admit = {});
struct NormalEffectPosterior {
    double mean{}, variance{}, low{}, high{};
};
NormalEffectPosterior shrink_normal_effect(double effect, double variance, double prior_mean,
                                           double prior_variance, double confidence);
} // namespace gadget::statistics
