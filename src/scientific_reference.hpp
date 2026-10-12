#pragma once
#include "scientific_detail.hpp"
#include <gadget/study.hpp>
#include "scientific_result_fields.hpp"
namespace gadget::scientific_detail {
Json score_ensemble(const std::vector<Vector> &scores, const Vector &weights,
                    const std::string &weighting);
struct FiniteScoreValues {
    std::vector<std::optional<Rational>> observed, reference_tail;
};
ResultValue score_ensemble_fields(const std::vector<Vector> &scores, const Vector &weights,
    const std::string &weighting, ResultFields &fields, std::uint64_t state_carrier,
    std::uint64_t hypothesis_carrier);
ResultValue calibrate_finite_scores_fields(const std::vector<Vector> &scores,
                             const std::vector<std::vector<bool>> &defined, const Vector &weights,
                             std::size_t observed_state, const Json &family,
                             const std::optional<Rational> &alpha, ResultFields &fields,
                             std::uint64_t state_carrier, std::uint64_t hypothesis_carrier,
                             FiniteScoreValues *native_values = nullptr);
// One finite-tail owner for model and structural reference laws.
Json calibrate_finite_scores(const std::vector<Vector> &scores,
                             const std::vector<std::vector<bool>> &defined, const Vector &weights,
                             std::size_t observed_state, const Json &family,
                             const std::optional<Rational> &alpha);
} // namespace gadget::scientific_detail
