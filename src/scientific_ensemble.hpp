#pragma once
#include "scientific_detail.hpp"
#include <optional>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/statistics/ensemble.hpp>

namespace gadget::scientific_detail {
// Weighted ensemble contractions on a declared metric carrier.
// All maps, ranks, Gram contractions and trace evaluations belong to RexGraph.
// This owner keeps the original coordinate fields; its smaller contraction
// space is an evaluation choice and does not truncate their span.
struct WeightedFieldEnsemble {
    rfc::numeric::RationalVector mean;
    rfc::numeric::ExactMatrix centered, contraction, contraction_metric;
    std::size_t rank{};
    bool sample_space{};
    std::optional<rfc::statistics::EffectiveDimension2> moments;
    WeightedFieldEnsemble(const std::vector<Vector> &fields, const Vector &weights,
                          const rfc::field::ExactMetric &metric) {
        require(!fields.empty() && fields.size() == weights.size(), "ensemble sample shape");
        Rational mass{0};
        mean.assign(metric.size(), Rational{0});
        for (std::size_t s = 0; s < fields.size(); ++s) {
            require(weights[s] > 0 && fields[s].size() == metric.size(),
                    "ensemble requires positive weights and identical coordinate carriers");
            mass += weights[s];
            mean = rfc::numeric::add(mean, rfc::numeric::scale(weights[s], fields[s]));
        }
        require(mass == 1, "ensemble weights must sum exactly to one");
        centered = rfc::numeric::ExactMatrix(fields.size(), metric.size());
        for (std::size_t s = 0; s < fields.size(); ++s) {
            const auto departure = rfc::numeric::sub(fields[s], mean);
            for (std::size_t i = 0; i < departure.size(); ++i)
                centered.set(s, i, departure[i]);
        }
        require(rfc::numeric::is_zero(centered.transpose().apply(weights)),
                "ensemble weighted centering");
        rank = centered.rank();
        sample_space = fields.size() < metric.size();
        if (sample_space) {
            const auto columns = centered.transpose();
            contraction = rfc::field::metric_gram(columns, metric, columns);
            contraction_metric = rfc::numeric::ExactMatrix::diagonal(weights);
        } else {
            contraction = rfc::field::metric_gram(
                centered, rfc::field::ExactMetric::diagonal(weights), centered);
            contraction_metric = metric.materialize();
        }
        if (rank) {
            moments = rfc::statistics::effective_dimension_2(contraction_metric, contraction);
            require(rank < fields.size() && moments->dimension >= 1 &&
                        moments->dimension <= Rational{rfc::numeric::Integer{std::to_string(rank)}},
                    "ensemble effective dimension bound");
        }
    }
};
// Direct sample-Gram evaluation of covariance-trace coupling. The two trace
// contractions are dual to its covariance expression; no eigensystem or
// cross-carrier pairing is introduced. The caller retains both complete fields.
inline Rational weighted_field_coupling(const WeightedFieldEnsemble &left,
                                        const WeightedFieldEnsemble &right, const Vector &weights,
                                        const rfc::field::ExactMetric &metric) {
    if (!left.sample_space) {
        const auto cross = rfc::field::metric_gram(
            left.centered, rfc::field::ExactMetric::diagonal(weights), right.centered);
        const auto form = metric.materialize();
        return rfc::statistics::cross_sector_coupling(form, left.contraction, cross, form,
                                                      right.contraction, cross.transpose());
    }
    const auto d = rfc::numeric::ExactMatrix::diagonal(weights);
    Rational value = rfc::statistics::trace_product(d.multiply(left.contraction).multiply(d),
                                                    right.contraction) /
                     (left.moments->variance * right.moments->variance);
    value.canonicalize();
    require(value >= 0 && value <= 1, "ensemble coupling bound");
    return value;
}
} // namespace gadget::scientific_detail
