#pragma once
#include <functional>
#include <limits>
#include <optional>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/checked_exact_actions.hpp>

namespace gadget::scientific_detail {
// Owns one declared frame and its canonical factorization. The metric is an
// explicit request input; all numerical actions are canonical SDK operations.
struct PreparedProjection {
    rfc::numeric::ExactMatrix frame, dual;
    rfc::field::ExactMetric metric;
    rfc::numeric::ExactActionLimits limits;
    rfc::numeric::PreparedExactSolve solve;
    explicit PreparedProjection(rfc::numeric::ExactMatrix independent_frame,
                                std::optional<rfc::field::ExactMetric> form = {},
                                std::optional<rfc::numeric::ExactActionLimits> action_limits = {},
                                std::function<void()> checkpoint = {})
        : frame(std::move(independent_frame)), dual(frame.transpose()),
          metric(form ? std::move(*form) : rfc::field::ExactMetric::identity(frame.rows())),
          limits(action_limits ? *action_limits : rfc::numeric::ExactActionLimits{frame.cols(),
                 rfc::machine::checked_size_mul(
                     rfc::machine::checked_size_mul(frame.cols(), frame.cols()), 2),
                 std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()}),
          solve(rfc::field::metric_gram(frame, metric, frame), limits, std::move(checkpoint)) {}
    rfc::numeric::RationalVector coefficients(const rfc::numeric::RationalVector &input,
                                              std::function<void()> checkpoint = {}) const {
        return solve.apply(dual.apply(metric.apply(input)), limits, std::move(checkpoint));
    }
    rfc::numeric::RationalVector apply(const rfc::numeric::RationalVector &input,
                                       std::function<void()> checkpoint = {}) const {
        return frame.apply(coefficients(input, std::move(checkpoint)));
    }
};
} // namespace gadget::scientific_detail
