#include "allocation_detail.hpp"
#include <numeric>

namespace gadget::assay_detail {
namespace {
using namespace allocation_detail;
struct Search {
    Integer tested{0}, nonsingular{0}, feasible{0};
    std::optional<Rational> optimum;
    Vector marginal_min, marginal_max, optimal_min, optimal_max;
    std::vector<Vector> optimizers;
};
Search solve_block(const Problem &p, const Block &local) {
    Search search;
    std::set<Vector> seen;
    const auto evaluate = [&](const Vector &z) {
        const auto values = local.normals.apply(z);
        for (std::size_t i = 0; i < values.size(); ++i)
            if (!satisfies(values[i], local.senses[i], local.rhs[i]))
                return;
        const auto x = rfc::numeric::add(local.embedding.apply(z), local.offset);
        if (p.domain == "bounded_integer")
            for (const auto &q : x)
                if (q.get_den() != 1)
                    return;
        if (p.domain == "bounded_rational" && !seen.insert(z).second)
            return;
        ++search.feasible;
        if (search.feasible == 1)
            search.marginal_min = search.marginal_max = x;
        for (std::size_t i = 0; i < x.size(); ++i) {
            search.marginal_min[i] = std::min(search.marginal_min[i], x[i]);
            search.marginal_max[i] = std::max(search.marginal_max[i], x[i]);
        }
        const auto value = local.objective.apply(z).at(0);
        if (!search.optimum || (p.direction == "minimize" ? value < *search.optimum : value > *search.optimum)) {
            search.optimum = value;
            search.optimizers.clear();
            search.optimal_min = search.optimal_max = x;
        }
        if (value == *search.optimum) {
            search.optimizers.push_back(z);
            for (std::size_t i = 0; i < x.size(); ++i) {
                search.optimal_min[i] = std::min(search.optimal_min[i], x[i]);
                search.optimal_max[i] = std::max(search.optimal_max[i], x[i]);
            }
        }
    };
    const auto d = local.embedding.cols();
    require(d > 0, "empty allocation factor");
    if (p.domain == "bounded_integer") {
        // Free parameters are actual source coordinates, not rounded kernel
        // coefficients. Reconstructed dependent coordinates must be integral.
        Integer box{1};
        for (std::size_t j = 0; j < d; ++j) {
            const Rational width = local.upper[j] - local.lower[j] + 1;
            box *= width.get_num();
        }
        Vector z = local.lower;
        for (;;) {
            ++search.tested;
            evaluate(z);
            std::size_t k = d;
            while (k > 0 && z[k - 1] == local.upper[k - 1]) {
                z[k - 1] = local.lower[k - 1];
                --k;
            }
            if (!k)
                break;
            ++z[k - 1];
        }
        require(search.tested == box, "complete local integer box traversal");
    } else {
        // Exact native RREF canonicalizes active hyperplanes. Scaled duplicate
        // rows share one hyperplane; every original row still tests feasibility.
        std::set<Vector> planes;
        for (std::size_t i = 0; i < local.normals.rows(); ++i) {
            ExactMatrix row(1, d + 1);
            for (std::size_t j = 0; j < d; ++j)
                row.set(0, j, local.normals.get(i, j));
            row.set(0, d, local.rhs[i]);
            std::vector<std::size_t> pivots;
            const auto normalized = row.rref(pivots);
            require(!pivots.empty() && pivots.front() < d, "nonconstant local hyperplane");
            Vector plane;
            for (std::size_t j = 0; j <= d; ++j)
                plane.push_back(normalized.get(0, j));
            planes.insert(plane);
        }
        const std::vector<Vector> active_rows(planes.begin(), planes.end());
        const auto rows = active_rows.size();
        require(rows >= d, "bounded affine factor lacks active normals");
        Integer combinations{1};
        for (std::size_t j = 1; j <= d; ++j) {
            combinations *= Integer(std::to_string(rows - j + 1));
            combinations /= Integer(std::to_string(j));
        }
        std::vector<std::size_t> active(d);
        std::iota(active.begin(), active.end(), 0);
        for (;;) {
            ++search.tested;
            ExactMatrix matrix(d, d);
            Vector rhs;
            for (std::size_t i = 0; i < d; ++i) {
                const auto &row = active_rows[active[i]];
                for (std::size_t j = 0; j < d; ++j)
                    matrix.set(i, j, row[j]);
                rhs.push_back(row[d]);
            }
            if (matrix.rank() == d) {
                ++search.nonsingular;
                const auto z = matrix.solve(rhs);
                require(matrix.apply(z) == rhs, "native local active-basis residual");
                evaluate(z);
            }
            std::size_t k = d;
            while (k > 0 && active[k - 1] == rows - d + k - 1)
                --k;
            if (!k)
                break;
            ++active[k - 1];
            for (std::size_t j = k; j < d; ++j)
                active[j] = active[j - 1] + 1;
        }
        require(search.tested == combinations, "complete distinct local active-basis traversal");
    }
    std::sort(search.optimizers.begin(), search.optimizers.end());
    return search;
}
} // namespace

Json constrained_allocations(const Json &request) {
    using namespace allocation_detail;
    const Problem p(request);
    const Factors f(p);
    const bool factorized = request.has("witness_representation") &&
                            text(request.at("witness_representation")) == "factorized";
    const bool certificate = factorized || (request.has("factor_certificate") &&
                                             request.at("factor_certificate").boolean());
    std::vector<Search> searches;
    Integer tested{0}, nonsingular{0}, feasible{f.consistent ? 1 : 0}, optimum_count{f.consistent ? 1 : 0};
    Rational optimum{0};
    Vector marginal_min = f.offset, marginal_max = f.offset,
           optimal_min = f.offset, optimal_max = f.offset;
    if (f.consistent) {
        optimum = p.objective.apply(f.offset).at(0);
        for (const auto &block : f.blocks) {
            auto search = solve_block(p, block);
            tested += search.tested;
            nonsingular += search.nonsingular;
            feasible *= search.feasible;
            optimum_count *= Integer(std::to_string(search.optimizers.size()));
            if (search.optimum) {
                optimum += *search.optimum;
                for (std::size_t i = 0; i < block.variables.size(); ++i) {
                    const auto j = block.variables[i];
                    marginal_min[j] = search.marginal_min[i];
                    marginal_max[j] = search.marginal_max[i];
                    optimal_min[j] = search.optimal_min[i];
                    optimal_max[j] = search.optimal_max[i];
                }
            }
            searches.push_back(std::move(search));
        }
    }
    auto out = Json::object(), bounds = Json::array(), witnesses = Json::array(), blocks = Json::array();
    out.set("status", Json::string(feasible != 0 ? "feasible" : "infeasible"));
    out.set("domain", request.at("domain"));
    out.set("variables", p.names);
    out.set("tested_candidates", exact_count(tested));
    out.set("nonsingular_bases", exact_count(nonsingular));
    out.set("feasible_points_or_vertices", exact_count(feasible));
    out.set("optimal_points_or_vertices", exact_count(optimum_count));
    out.set("equality_rank", cardinality(f.equality_rank));
    if (f.consistent)
        out.set("affine_dimension", cardinality(f.kernel.cols()));
    out.set("independent_blocks", cardinality(f.blocks.size()));
    out.set("constraint_carrier", f.carrier);
    if (p.domain == "bounded_integer")
        out.set("declared_box_points", exact_count(p.declared_box));
    out.set("complete", Json::boolean(true));
    out.set("witness_representation", Json::string(factorized ? "factorized" : "explicit"));
    if (feasible != 0) {
        out.set("minimum_or_maximum", Json::string(optimum.get_str()));
        for (std::size_t i = 0; i < p.lower.size(); ++i) {
            auto b = Json::object();
            b.set("id", p.names.at(i));
            b.set("lower", Json::string(marginal_min[i].get_str()));
            b.set("upper", Json::string(marginal_max[i].get_str()));
            b.set("optimal_lower", Json::string(optimal_min[i].get_str()));
            b.set("optimal_upper", Json::string(optimal_max[i].get_str()));
            bounds.push(b);
        }
        if (!factorized) {
            // The caller explicitly requested every Cartesian witness. Never
            // truncate or substitute factors for this output contract.
            std::vector<std::size_t> selected(searches.size(), 0);
            Integer emitted{0};
            for (;;) {
                Vector z(f.kernel.cols(), Rational{0});
                for (std::size_t i = 0; i < searches.size(); ++i) {
                    const auto &tuple = searches[i].optimizers[selected[i]];
                    for (std::size_t j = 0; j < tuple.size(); ++j)
                        z[f.blocks[i].parameters[j]] = tuple[j];
                }
                const auto x = rfc::numeric::add(f.kernel.apply(z), f.offset);
                for (std::size_t i = 0; i < x.size(); ++i) {
                    require(x[i] >= p.lower[i] && x[i] <= p.upper[i] &&
                                (p.domain != "bounded_integer" || x[i].get_den() == 1),
                            "complete allocation witness bounds/integrality");
                }
                const auto values = p.constraints.apply(x);
                for (std::size_t i = 0; i < values.size(); ++i)
                    require(satisfies(values[i], p.senses[i], p.rhs[i]), "source constraint witness residual");
                require(p.objective.apply(x).at(0) == optimum, "source objective witness residual");
                auto w = Json::object();
                w.set("values", vector_json(x));
                w.set("constraint_values", vector_json(values));
                w.set("objective", Json::string(optimum.get_str()));
                witnesses.push(w);
                ++emitted;
                std::size_t k = selected.size();
                while (k > 0 && selected[k - 1] + 1 == searches[k - 1].optimizers.size()) {
                    selected[k - 1] = 0;
                    --k;
                }
                if (!k)
                    break;
                ++selected[k - 1];
            }
            require(emitted == optimum_count, "complete optimal Cartesian witness expansion");
        }
    }
    for (std::size_t i = 0; i < searches.size(); ++i) {
        const auto &search = searches[i];
        auto item = Json::object();
        item.set("block_index", cardinality(i));
        item.set("tested_candidates", exact_count(search.tested));
        item.set("nonsingular_bases", exact_count(search.nonsingular));
        item.set("feasible_points_or_vertices", exact_count(search.feasible));
        item.set("optimal_points_or_vertices", cardinality(search.optimizers.size()));
        if (search.optimum)
            item.set("objective_contribution", Json::string(search.optimum->get_str()));
        if (certificate) {
            auto tuples = Json::array();
            for (const auto &z : search.optimizers)
                tuples.push(vector_json(z));
            item.set("optimal_parameter_tuples", tuples);
        }
        blocks.push(item);
    }
    out.set("bounds", bounds);
    out.set("optimizers", witnesses);
    out.set("block_results", blocks);
    if (certificate)
        out.set("factor_certificate", f.certificate(p));
    if (!f.consistent)
        out.set("contradiction", f.contradiction);
    out.set("optimizer_semantics", Json::string(p.domain == "bounded_integer"
        ? "All optimal integer points, reconstructed from the Cartesian product of complete block tuples"
        : "All optimal vertices, reconstructed from complete block tuples; their convex hull is the complete optimal face"));
    out.set("certificate_method", Json::string("Native exact affine kernel and support-indicator factors; complete local finite search and original constraint residuals"));
    return out;
}
} // namespace gadget::assay_detail
