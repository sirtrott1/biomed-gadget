#pragma once
#include "scientific_model.hpp"
#include <algorithm>
#include <numeric>
#include <optional>
#include <rfc/path/objectives.hpp>

namespace gadget::scientific_detail {
inline std::vector<std::size_t>
physical_shape_representatives(const rfc::numeric::ExactMatrix &normals,
                               const rfc::numeric::RationalVector &rhs) {
    using namespace rfc::numeric;
    std::set<std::pair<bool, RationalVector>> seen;
    std::vector<std::size_t> representatives;
    for (std::size_t i = 0; i < normals.rows(); ++i) {
        ExactMatrix row(1, rfc::machine::checked_size_add(normals.cols(), 1));
        for (std::size_t j = 0; j < normals.cols(); ++j)
            row.set(0, j, normals.get(i, j));
        row.set(0, normals.cols(), rhs.at(i));
        std::vector<std::size_t> pivots;
        const auto normalized = row.rref(pivots);
        if (pivots.empty() || pivots.front() == normals.cols()) {
            if (rhs[i] < 0)
                throw std::invalid_argument("constant physical constraint is infeasible");
            continue;
        }
        RationalVector plane;
        for (std::size_t j = 0; j <= normals.cols(); ++j)
            plane.push_back(normalized.get(0, j));
        // The sign retains inequality orientation; only positive proportional
        // rows share a representative. Every original row still tests the fit.
        const bool positive = normals.get(i, pivots.front()) > 0;
        if (seen.emplace(positive, std::move(plane)).second)
            representatives.push_back(i);
    }
    return representatives;
}
// Application face declarations; every affine/quadratic action belongs to RFC.
// The returned exact primal/dual witness proves the global constrained minimum.
struct PhysicalShapeFit {
    rfc::numeric::RationalVector coordinates, multipliers, slack, stationarity;
    std::vector<std::size_t> active;
    rfc::numeric::Integer tested{0}, dependent{0};
};
inline PhysicalShapeFit physical_shape_fit(const rfc::numeric::ExactMatrix &normals,
                                           const rfc::numeric::RationalVector &rhs,
                                           const rfc::field::ExactMetric &gram,
                                           const rfc::numeric::RationalVector &unconstrained,
                                           const std::vector<std::size_t> &representatives,
                                           const std::function<void()> &checkpoint) {
    using namespace rfc::numeric;
    if (normals.rows() != rhs.size() || normals.cols() != gram.size() ||
        unconstrained.size() != gram.size())
        throw std::invalid_argument("physical shape coordinate domains");
    PhysicalShapeFit out;
    const auto admit = [&](const std::vector<std::size_t> &active) {
        checkpoint();
        ExactMatrix h(active.size(), gram.size());
        RationalVector targets;
        for (std::size_t i = 0; i < active.size(); ++i) {
            for (std::size_t j = 0; j < gram.size(); ++j)
                h.set(i, j, normals.get(active[i], j));
            targets.push_back(rhs[active[i]]);
        }
        if (h.rank() != active.size()) {
            ++out.dependent;
            return false;
        }
        ++out.tested;
        const auto ht = h.transpose(), euclidean = h.multiply(ht);
        const auto demand = sub(targets, h.apply(unconstrained));
        const auto feasible = ht.apply(euclidean.solve(demand));
        const auto correction = rfc::path::least_quadrance_feasible(h, gram, feasible);
        if (h.apply(correction) != demand)
            throw std::runtime_error("native physical affine face residual");
        const auto v = add(unconstrained, correction);
        const auto slack = sub(rhs, normals.apply(v));
        if (std::any_of(slack.begin(), slack.end(), [](const auto &q) { return q < 0; }))
            return false;
        const auto gradient = gram.apply(correction);
        const auto multipliers = scale(Rational{-1}, euclidean.solve(h.apply(gradient)));
        if (std::any_of(multipliers.begin(), multipliers.end(),
                        [](const auto &q) { return q < 0; }))
            return false;
        const auto stationarity = add(gradient, ht.apply(multipliers));
        if (!is_zero(stationarity))
            throw std::runtime_error("native physical dual residual");
        out.coordinates = v;
        out.active = active;
        out.multipliers = multipliers;
        out.slack = slack;
        out.stationarity = stationarity;
        return true;
    };
    if (admit({}))
        return out;
    // Observed violated rows are a deterministic first candidate; the certificate
    // is authoritative even when a different face is needed.
    const auto initial_slack = sub(rhs, normals.apply(unconstrained));
    auto order = representatives;
    std::stable_sort(order.begin(), order.end(),
                     [&](auto a, auto b) { return initial_slack[a] < initial_slack[b]; });
    std::vector<std::size_t> seed;
    for (auto row : order) {
        if (initial_slack[row] >= 0)
            break;
        auto candidate = seed;
        candidate.push_back(row);
        ExactMatrix h(candidate.size(), gram.size());
        for (std::size_t i = 0; i < candidate.size(); ++i)
            for (std::size_t j = 0; j < gram.size(); ++j)
                h.set(i, j, normals.get(candidate[i], j));
        if (h.rank() == candidate.size())
            seed = std::move(candidate);
    }
    if (!seed.empty() && admit(seed))
        return out;
    // Complete rank derived independent affine faces. No truncation or policy
    // ceiling. Stop as soon as the native witness proves global optimality.
    for (std::size_t k = 1; k <= std::min(order.size(), gram.size()); ++k) {
        std::vector<std::size_t> chosen(k);
        std::iota(chosen.begin(), chosen.end(), 0);
        for (;;) {
            std::vector<std::size_t> active;
            for (auto i : chosen)
                active.push_back(order[i]);
            if (admit(active))
                return out;
            std::size_t j = k;
            while (j && chosen[j - 1] == order.size() - k + j - 1)
                --j;
            if (!j)
                break;
            ++chosen[j - 1];
            for (std::size_t i = j; i < k; ++i)
                chosen[i] = chosen[i - 1] + 1;
        }
    }
    throw std::runtime_error("physical shape profile lacks its exact optimality witness");
}
} // namespace gadget::scientific_detail
