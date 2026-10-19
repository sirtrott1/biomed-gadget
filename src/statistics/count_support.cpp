#include "count_support.hpp"
#include "../scientific_shape.hpp"
#include <rfc/runtime/execution_policy.hpp>

namespace gadget::statistics {
bool finite_count_coefficients(const rfc::numeric::ExactMatrix &design,
                               std::span<const double> counts, CountWork &work,
                               std::stop_token stop) {
    using namespace rfc::numeric;
    using rfc::machine::checked_size_mul;
    const auto n = design.rows(), p = design.cols();
    if (!p)
        return true;
    const auto positive =
        std::count_if(counts.begin(), counts.end(), [](double value) { return value > 0; });
    if (static_cast<std::size_t>(positive) == n)
        return true;
    const auto checkpoint = [&] {
        work.charge(checked_size_mul(n, checked_size_mul(p, checked_size_mul(p, p))), stop);
    };
    checkpoint();
    ExactMatrix observed(positive, p);
    std::size_t row{};
    for (std::size_t i = 0; i < n; ++i)
        if (counts[i] > 0) {
            for (std::size_t j = 0; j < p; ++j)
                observed.set(row, j, design.get(i, j));
            ++row;
        }
    const auto kernel = observed.kernel_frame();
    const auto d = kernel.cols();
    if (!d)
        return true;
    // A recession direction is zero on positive counts, nonpositive on zeros,
    // and strictly negative on at least one zero. Normalize its sum to minus one.
    ExactMatrix zeros(n - positive, p);
    row = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (counts[i] == 0) {
            for (std::size_t j = 0; j < p; ++j)
                zeros.set(row, j, design.get(i, j));
            ++row;
        }
    const auto z = zeros.multiply(kernel);
    RationalVector sum(d);
    for (std::size_t i = 0; i < z.rows(); ++i)
        for (std::size_t j = 0; j < d; ++j)
            sum[j] += z.get(i, j);
    if (is_zero(sum))
        return true;
    ExactMatrix normal(1, d);
    for (std::size_t j = 0; j < d; ++j)
        normal.set(0, j, sum[j]);
    const auto transpose = normal.transpose();
    const auto origin = transpose.apply(normal.multiply(transpose).solve({Rational{-1}}));
    const auto free = normal.kernel_frame(), a = z.multiply(free);
    auto rhs = scale(Rational{-1}, z.apply(origin));
    if (!free.cols())
        return std::any_of(rhs.begin(), rhs.end(), [](const auto &value) { return value < 0; });
    try {
        checkpoint();
        const auto representatives = scientific_detail::physical_shape_representatives(a, rhs);
        (void)scientific_detail::physical_shape_fit(
            a, rhs, rfc::field::ExactMetric::identity(free.cols()), RationalVector(free.cols()),
            representatives, checkpoint);
        return false;
    } catch (const CountWorkLimit &) {
        throw;
    } catch (const rfc::runtime::ExecutionCancelled &) {
        throw;
    } catch (const std::invalid_argument &error) {
        if (std::string_view(error.what()) == "constant physical constraint is infeasible")
            return true;
        throw;
    } catch (const std::runtime_error &error) {
        if (std::string_view(error.what()) ==
            "physical shape profile lacks its exact optimality witness")
            return true;
        throw;
    }
}
} // namespace gadget::statistics
