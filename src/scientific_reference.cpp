#include "scientific_reference.hpp"
#include "scientific_ensemble.hpp"
#include <algorithm>
namespace gadget::scientific_detail {
using namespace rfc;
namespace {
using Matrix = numeric::ExactMatrix;
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
Json matrix_json(const Matrix &matrix) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(matrix.rows()));
    out.set("columns", cardinality(matrix.cols()));
    for (const auto &[position, value] : matrix.entries()) {
        auto entry = Json::object();
        entry.set("row", cardinality(position.first));
        entry.set("column", cardinality(position.second));
        entry.set("value", rational(value));
        entries.push(entry);
    }
    out.set("entries", entries);
    return out;
}
// Inclusive tails of the declared discrete law; no fitted sampling distribution.
struct Tails {
    std::map<Rational, Rational> masses, inclusive;
    explicit Tails(const Vector &values, const Vector &weights) {
        require(values.size() == weights.size(), "finite law tail shape");
        for (std::size_t s = 0; s < values.size(); ++s)
            masses[values[s]] += weights[s];
        Rational accumulated{0};
        for (auto it = masses.rbegin(); it != masses.rend(); ++it) {
            accumulated += it->second;
            inclusive.emplace(it->first, accumulated);
        }
        // Check P(p(T) <= a) <= a at every attainable tail level.
        std::map<Rational, Rational> tail_distribution;
        for (std::size_t s = 0; s < values.size(); ++s)
            tail_distribution[inclusive.at(values[s])] += weights[s];
        Rational probability{0};
        for (const auto &[alpha, mass] : tail_distribution) {
            probability += mass;
            require(probability <= alpha, "finite law tail calibration");
        }
    }
    Rational at(const Rational &threshold) const {
        const auto found = inclusive.lower_bound(threshold);
        return found == inclusive.end() ? Rational{0} : found->second;
    }
    Json distribution() const {
        auto out = Json::array();
        for (const auto &[value, mass] : masses) {
            auto record = Json::object();
            record.set("value", rational(value));
            record.set("mass", rational(mass));
            record.set("inclusive_tail_mass", rational(inclusive.at(value)));
            out.push(record);
        }
        return out;
    }
};
} // namespace
Json score_ensemble(const std::vector<Vector> &scores, const Vector &weights,
                    const std::string &weighting) {
    const auto states = scores.size(), hypotheses = scores.front().size();
    const scientific_detail::WeightedFieldEnsemble summary(
        scores, weights, field::ExactMetric::identity(hypotheses));
    const auto &mean = summary.mean;
    const auto &centered = summary.centered;
    auto fields = Json::array();
    for (std::size_t s = 0; s < states; ++s) {
        const auto values = numeric::sub(scores[s], mean);
        fields.push(vector_json(values));
    }
    require(numeric::is_zero(centered.transpose().apply(weights)), "finite law weighted centering");
    const auto rank = summary.rank;
    auto out = Json::object();
    out.set("mean", vector_json(mean));
    out.set("centered_statistics", fields);
    out.set("rank", cardinality(rank));
    out.set("metric", Json::string("identity on the declared dimensionless hypothesis axes"));
    out.set("normalization", Json::string(weighting));
    // Select coordinate covariance or sample-Gram space from the dimensions.
    // A positive rational law metric avoids square roots of unequal masses.
    const bool state_space = summary.sample_space;
    const auto &form = summary.contraction;
    const auto &metric = summary.contraction_metric;
    out.set("contraction_space", Json::string(state_space ? "reference_states" : "hypotheses"));
    out.set("contraction_metric", matrix_json(metric));
    out.set("second_moment_or_gram", matrix_json(form));
    if (rank == 0) {
        out.set("status", Json::string("constant_statistics"));
        out.set("variance", rational(0));
        out.set("second_contraction", rational(0));
        out.set("effective_dimension_2", Json{});
        out.set("collision_probability", Json{});
    } else {
        const auto moments = *summary.moments;
        require(rank < states && moments.dimension >= 1 &&
                    moments.dimension <= Rational{numeric::Integer{std::to_string(rank)}},
                "finite reference effective dimension bound");
        out.set("status", Json::string("completed"));
        out.set("variance", rational(moments.variance));
        out.set("second_contraction", rational(moments.second_contraction));
        out.set("effective_dimension_2", rational(moments.dimension));
        out.set("collision_probability", rational(moments.collision_probability));
    }
    out.set("weighted_centering_exact", Json::boolean(true));
    return out;
}
Json calibrate_finite_scores(const std::vector<Vector> &scores,
                             const std::vector<std::vector<bool>> &defined, const Vector &weights,
                             std::size_t observed_state, const Json &family,
                             const std::optional<Rational> &alpha) {
    require(!scores.empty() && scores.size() == defined.size() && scores.size() == weights.size() &&
                observed_state < scores.size() && family.size(),
            "finite score family shape");
    std::vector<bool> available(family.size(), true);
    Rational total{0};
    for (std::size_t s = 0; s < scores.size(); ++s) {
        require(scores[s].size() == family.size() && defined[s].size() == family.size() &&
                    weights[s] > 0,
                "finite score axes and positive mass");
        total += weights[s];
        for (std::size_t h = 0; h < family.size(); ++h)
            available[h] = available[h] && defined[s][h];
    }
    require(total == 1, "finite score mass must sum exactly to one");
    const bool complete = std::all_of(available.begin(), available.end(), [](bool a) { return a; });
    Vector maxima(scores.size(), Rational{0});
    if (complete)
        for (std::size_t s = 0; s < scores.size(); ++s)
            maxima[s] = *std::max_element(scores[s].begin(), scores[s].end());
    const std::optional<Tails> maximum_tails =
        complete ? std::optional<Tails>{Tails{maxima, weights}} : std::nullopt;
    auto results = Json::array();
    for (std::size_t h = 0; h < family.size(); ++h) {
        auto result = Json::parse(family.at(h).dump_compact());
        result.set("status", Json::string(available[h] ? "completed" : "undefined_statistic"));
        result.set("observed", defined.at(observed_state).at(h)
                                   ? rational(scores.at(observed_state).at(h))
                                   : Json{});
        result.set("reference_tail_mass", Json{});
        result.set("max_statistic_tail_mass", Json{});
        if (available[h]) {
            Vector values;
            for (const auto &state_scores : scores)
                values.push_back(state_scores[h]);
            const Tails tails(values, weights);
            const auto threshold = values[observed_state], tail = tails.at(threshold);
            require(tail >= weights[observed_state] && tail <= 1,
                    "inclusive finite tail retains observed state mass");
            result.set("reference_tail_mass", rational(tail));
            result.set("tie_mass", rational(tails.masses.at(threshold)));
            result.set("strict_tail_mass", rational(Rational{tail - tails.masses.at(threshold)}));
            result.set("distribution", tails.distribution());
            if (alpha)
                result.set("reference_tail_at_or_below_alpha", Json::boolean(tail <= *alpha));
            if (maximum_tails) {
                const auto adjusted = maximum_tails->at(threshold);
                require(adjusted >= tail, "joint maximum tail must dominate marginal tail");
                result.set("max_statistic_tail_mass", rational(adjusted));
                if (alpha)
                    result.set("max_tail_at_or_below_alpha", Json::boolean(adjusted <= *alpha));
            }
        }
        results.push(result);
    }
    auto joint = Json::object();
    joint.set("status", Json::string(complete ? "completed" : "undefined_family_statistic"));
    joint.set("global_tail_mass",
              maximum_tails ? rational(maximum_tails->at(maxima[observed_state])) : Json{});
    if (maximum_tails) {
        joint.set("observed_maximum", rational(maxima[observed_state]));
        joint.set("maxima", vector_json(maxima));
        joint.set("distribution", maximum_tails->distribution());
        joint.set(
            "ensemble",
            score_ensemble(scores, weights, "declared probability masses, summing exactly to one"));
    }
    auto out = Json::object();
    out.set("status", Json::string(complete ? "completed" : "undefined_family_statistic"));
    out.set("hypotheses", results);
    out.set("joint_family", joint);
    return out;
}
} // namespace gadget::scientific_detail
