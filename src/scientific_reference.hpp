#pragma once
#include "scientific_detail.hpp"
#include <gadget/study.hpp>
namespace gadget::scientific_detail {
Json score_ensemble(const std::vector<Vector> &scores, const Vector &weights,
                    const std::string &weighting);
// One finite-tail owner for model and structural reference laws.
Json calibrate_finite_scores(const std::vector<Vector> &scores,
                             const std::vector<std::vector<bool>> &defined, const Vector &weights,
                             std::size_t observed_state, const Json &family,
                             const std::optional<Rational> &alpha);
} // namespace gadget::scientific_detail
