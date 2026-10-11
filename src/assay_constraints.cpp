#include "allocation_detail.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
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

allocation_detail::Result allocation_detail::solve_fields(const Problem &p) {
    const auto &request = p.declaration;
    p.validate();
    using namespace allocation_detail;
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

    using N = NativeRecordValue;
    using namespace rfc::storage;
    TemporalEndpoint endpoint(*f.source);
    endpoint.carriers.push_back({1, TemporalCarrierKind::SourceGrade, {0}, {}, 0,
                                rfc::numeric::IndexDomain::linear(p.lower.size())});
    endpoint.carriers.push_back({2, TemporalCarrierKind::SourceGrade, {1}, {}, 0,
                                rfc::numeric::IndexDomain::linear(p.rhs.size())});
    SourceGradeFields store(endpoint);
    N::Object tensors;
    const auto put = [&](std::string name, N key) {
        tensors.members.emplace_back(std::move(name), std::move(key));
    };
    const auto column = [&](std::string name, std::uint64_t carrier, const Vector &v) {
        auto key = store.vector(name, carrier, v);
        put(name, key);
        return key;
    };
    const auto map = [&](std::string name, std::uint64_t input, std::uint64_t output,
                         const ExactMatrix &value) {
        const auto size = [&](auto key) {
            for (const auto &carrier : endpoint.carriers)
                if (carrier.key == key) return carrier.domain.logical_size();
            throw std::invalid_argument("constraint map carrier absent");
        };
        if (!f.consistent && (value.cols() != size(input) || value.rows() != size(output)))
            put(name, N{});
        else
            put(name, store.map(name, input, output, value));
    };
    const auto indices = [&](std::string name, const std::vector<std::size_t> &v) {
        Vector values;
        for (auto index : v) values.emplace_back(std::to_string(index));
        return column(name, store.explicit_carrier(values.size()), values);
    };
    const auto scalar = store.explicit_carrier(1),
               equalities = store.explicit_carrier(f.equalities.rows()),
               parameters = store.explicit_carrier(f.kernel.cols()),
               reduced = store.explicit_carrier(f.inequalities.rows()),
               grouping = store.explicit_carrier(f.grouping.rows()),
               factors = store.explicit_carrier(f.indicator.cols());
    map("constraints", 1, 2, p.constraints);
    map("objective", 1, scalar, p.objective);
    column("rhs", 2, p.rhs);
    column("declared_lower", 1, p.lower);
    column("declared_upper", 1, p.upper);
    map("equality_operator", 1, equalities, f.equalities);
    column("equality_rhs", equalities, f.equality_rhs);
    if (f.offset.size() == p.lower.size()) column("offset", 1, f.offset);
    map("kernel", parameters, 1, f.kernel);
    map("reduced_inequalities", parameters, reduced, f.inequalities);
    column("reduced_rhs", reduced, f.rhs);
    map("grouping_contrasts", parameters, grouping, f.grouping);
    map("group_indicators", factors, parameters, f.indicator);
    map("reduced_objective", parameters, scalar, f.objective);
    indices("free_variable_indices", f.free_variables);
    if (!f.contradiction_left.empty()) column("contradiction_left", equalities, f.contradiction_left);
    N::Array declarations, blocks;
    for (const auto &[name, values] : std::vector<std::pair<std::string, const Vector *>>{
            {"lower", &marginal_min}, {"upper", &marginal_max},
            {"optimal_lower", &optimal_min}, {"optimal_upper", &optimal_max}}) {
        const auto key = store.column(name, 0, p.lower.size(), [&](auto i) { return values->at(i); },
            feasible != 0 ? rfc::numeric::StructuralSupport::full(p.lower.size())
                          : rfc::numeric::StructuralSupport::sparse(p.lower.size(), {}));
        put(name, key);
        declarations.elements.emplace_back(N::Object{{{"id", N(name == "lower" ? "feasible_lower" : name == "upper" ? "feasible_upper" : name)}, {"key", key},
            {"quantity", N(std::string("constrained_allocation_bound"))},
            {"units", N(text(request.at("units")))},
            {"method", N(std::string("Declared constraints and complete exact local search"))}}});
    }
    const auto count = [](const Integer &value) { return N(value); };
    for (std::size_t i = 0; i < f.blocks.size(); ++i) {
        const auto &block = f.blocks[i];
        const auto prefix = "block:" + std::to_string(i) + ":";
        const auto local_parameters = store.explicit_carrier(block.parameters.size()),
                   local_variables = store.explicit_carrier(block.variables.size()),
                   local_rows = store.explicit_carrier(block.inequalities.size());
        N::Object item{{{"parameter_indices", indices(prefix + "parameter_indices", block.parameters)},
                       {"variable_indices", indices(prefix + "variable_indices", block.variables)},
                       {"inequality_indices", indices(prefix + "inequality_indices", block.inequalities)}}};
        const auto block_column = [&](std::string name, std::uint64_t carrier, const Vector &v) {
            item.members.emplace_back(name, column(prefix + name, carrier, v));
        };
        const auto block_map = [&](std::string name, std::uint64_t input, std::uint64_t output,
                                   const ExactMatrix &matrix) {
            auto key = store.map(prefix + name, input, output, matrix);
            item.members.emplace_back(name, key);
        };
        block_map("embedding", local_parameters, local_variables, block.embedding);
        block_map("normals", local_parameters, local_rows, block.normals);
        block_map("objective", local_parameters, scalar, block.objective);
        block_column("offset", local_variables, block.offset);
        block_column("rhs", local_rows, block.rhs);
        block_column("lower", local_parameters, block.lower);
        block_column("upper", local_parameters, block.upper);
        const auto &search = searches.at(i);
        ExactMatrix tuples(block.parameters.size(), search.optimizers.size());
        for (std::size_t j = 0; j < search.optimizers.size(); ++j)
            for (std::size_t k = 0; k < block.parameters.size(); ++k)
                tuples.set(k, j, search.optimizers[j][k]);
        block_map("optimal_parameter_tuples", store.explicit_carrier(search.optimizers.size()),
                  local_parameters, tuples);
        item.members.emplace_back("tested_candidates", count(search.tested));
        item.members.emplace_back("nonsingular_bases", count(search.nonsingular));
        item.members.emplace_back("feasible_points_or_vertices", count(search.feasible));
        item.members.emplace_back("optimal_points_or_vertices", count(Integer(std::to_string(search.optimizers.size()))));
        if (search.optimum) item.members.emplace_back("objective_contribution", N(*search.optimum));
        blocks.elements.emplace_back(std::move(item));
    }
    auto summary = Json::object();
    summary.set("status", Json::string(feasible != 0 ? "feasible" : "infeasible"));
    summary.set("domain", request.at("domain"));
    summary.set("equality_rank", cardinality(f.equality_rank));
    if (f.consistent) summary.set("affine_dimension", cardinality(f.kernel.cols()));
    summary.set("independent_blocks", cardinality(f.blocks.size()));
    summary.set("constraint_carrier", f.carrier);
    summary.set("complete", Json::boolean(true));
    summary.set("witness_representation", Json::string(factorized ? "factorized" : "explicit"));
    summary.set("optimizer_semantics", Json::string(p.domain == "bounded_integer"
        ? "All optimal integer points, reconstructed from the Cartesian product of complete block tuples"
        : "All optimal vertices, reconstructed from complete block tuples; their convex hull is the complete optimal face"));
    summary.set("certificate_method", Json::string("Native exact affine kernel and support-indicator factors; complete local finite search and original constraint residuals"));
    N::Array senses, reduced_senses;
    for (const auto &sense : p.senses) senses.elements.emplace_back(sense);
    for (const auto &sense : f.senses) reduced_senses.elements.emplace_back(sense);
    auto contradiction = Json::parse(f.contradiction.dump());
    for (const auto *name : {"left_null_vector", "rhs_action", "rhs", "value"})
        json_object_object_del(contradiction.raw(), name);
    N declaration(N::Object{{{"schema", N(std::string("gadget.constraints.fields.v1"))},
        {"fields", N(std::move(declarations))}, {"tensors", N(std::move(tensors))},
        {"variables", import_record_view(p.names)}, {"protocol", import_record_view(request)},
        {"summary", import_record_view(summary)}, {"consistent", N(f.consistent)},
        {"factor_certificate", N(certificate)}, {"senses", N(std::move(senses))},
        {"reduced_senses", N(std::move(reduced_senses))},
        {"equality_origins", import_record_view(f.equality_origins)},
        {"inequality_origins", import_record_view(f.inequality_origins)},
        {"blocks", N(std::move(blocks))}, {"contradiction", import_record_view(contradiction)},
        {"tested_candidates", count(tested)}, {"nonsingular_bases", count(nonsingular)},
        {"feasible_points_or_vertices", count(feasible)},
        {"optimal_points_or_vertices", count(optimum_count)},
        {"declared_box_points", count(p.declared_box)},
        {"minimum_or_maximum", feasible != 0 ? N(optimum) : N{}},
        {"cell_ids", N(N::Array{{import_record_view(p.names), N{}}})}}});
    return {{std::move(endpoint), std::move(declaration), {}},
            feasible != 0 ? std::optional<Rational>(optimum) : std::nullopt};
}

Json constrained_allocations(const Json &request, std::optional<AssayFieldsSource> *retained) {
    auto result = allocation_detail::solve_fields(allocation_detail::Problem(request));
    if (retained) {
        retained->emplace(std::move(result.fields));
        return Json::object();
    }
    const auto record = NativeRecord::create_from_endpoint(std::move(result.fields.declaration),
        new_native_record_namespace(), std::move(result.fields.sources), std::move(result.fields.endpoint));
    return constraint_fields_view(record);
}

Json constraint_request_view(const NativeRecord &record) {
    auto protocol = render_record_view(record.at("protocol"));
    if (!protocol.has("typed_source")) return protocol;
    json_object_object_del(protocol.raw(), "typed_source");
    const auto keys = record.at("tensors");
    const auto &lower = retained_tensor(record, keys.at("declared_lower")).tensor,
               &upper = retained_tensor(record, keys.at("declared_upper")).tensor,
               &rhs = retained_tensor(record, keys.at("rhs")).tensor;
    const auto constraints = rfc::storage::temporal_tensor_matrix(retained_tensor(record, keys.at("constraints"))),
               objective = rfc::storage::temporal_tensor_matrix(retained_tensor(record, keys.at("objective")));
    for (std::size_t i = 0; i < protocol.at("variables").size(); ++i) {
        protocol.at("variables").at(i).set("lower", Json::string(lower.state(i)->get_str()));
        protocol.at("variables").at(i).set("upper", Json::string(upper.state(i)->get_str()));
    }
    for (std::size_t i = 0; i < constraints.rows(); ++i) {
        Vector row;
        for (std::size_t j = 0; j < constraints.cols(); ++j) row.push_back(constraints.get(i, j));
        auto item = protocol.at("constraints").at(i);
        item.set("coefficients", allocation_detail::vector_json(row));
        item.set("sense", Json::string(record.at("senses").at(i).text()));
        item.set("rhs", Json::string(rhs.state(i)->get_str()));
    }
    Vector row;
    for (std::size_t j = 0; j < objective.cols(); ++j) row.push_back(objective.get(0, j));
    protocol.at("objective").set("coefficients", allocation_detail::vector_json(row));
    return protocol;
}
Json constraint_fields_view(const NativeRecord &record) {
    using namespace allocation_detail;
    const auto keys = record.at("tensors");
    const auto vector = [&](const NativeRecord &key) {
        const auto &tensor = retained_tensor(record, key).tensor;
        Vector result;
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i)
            result.push_back(tensor.state(i).value());
        return result;
    };
    const auto matrix = [&](const NativeRecord &key) {
        return key.kind() == NativeRecordKind::Null ? ExactMatrix{}
            : rfc::storage::temporal_tensor_matrix(retained_tensor(record, key));
    };
    const auto indices = [&](const NativeRecord &key) {
        auto result = Json::array();
        for (const auto &x : vector(key)) {
            require(x >= 0 && x.get_den() == 1, "retained constraint coordinate");
            result.push(cardinality(x.get_num().get_ui()));
        }
        return result;
    };
    auto out = render_record_view(record.at("summary"));
    out.set("variables", render_record_view(record.at("variables")));
    for (const auto *name : {"tested_candidates", "nonsingular_bases", "feasible_points_or_vertices",
                              "optimal_points_or_vertices"})
        out.set(name, exact_count(record.at(name).integer()));
    if (out.at("domain").text() == "bounded_integer")
        out.set("declared_box_points", exact_count(record.at("declared_box_points").integer()));
    auto bounds = Json::array(), witnesses = Json::array(), blocks = Json::array(), factor_blocks = Json::array();
    const bool feasible = out.at("status").text() == "feasible";
    if (feasible) {
        out.set("minimum_or_maximum", Json::string(record.at("minimum_or_maximum").rational().get_str()));
        const auto lo = vector(keys.at("lower")), hi = vector(keys.at("upper")),
                   olo = vector(keys.at("optimal_lower")), ohi = vector(keys.at("optimal_upper"));
        for (std::size_t i = 0; i < lo.size(); ++i) {
            auto item = Json::object();
            item.set("id", render_record_view(record.at("variables").at(i)));
            item.set("lower", Json::string(lo[i].get_str())); item.set("upper", Json::string(hi[i].get_str()));
            item.set("optimal_lower", Json::string(olo[i].get_str())); item.set("optimal_upper", Json::string(ohi[i].get_str()));
            bounds.push(item);
        }
    }
    std::vector<ExactMatrix> tuples;
    std::vector<std::vector<std::size_t>> parameter_indices;
    for (std::size_t i = 0; i < record.at("blocks").size(); ++i) {
        const auto block = record.at("blocks").at(i);
        auto item = Json::object(), factor = Json::object();
        item.set("block_index", cardinality(i));
        for (const auto *name : {"tested_candidates", "nonsingular_bases", "feasible_points_or_vertices", "optimal_points_or_vertices"})
            item.set(name, exact_count(block.at(name).integer()));
        if (block.has("objective_contribution"))
            item.set("objective_contribution", Json::string(block.at("objective_contribution").rational().get_str()));
        tuples.push_back(matrix(block.at("optimal_parameter_tuples")));
        std::vector<std::size_t> parameters;
        for (const auto &x : vector(block.at("parameter_indices"))) parameters.push_back(x.get_num().get_ui());
        parameter_indices.push_back(std::move(parameters));
        if (record.at("factor_certificate").boolean()) {
            auto values = Json::array();
            for (std::size_t j = 0; j < tuples.back().cols(); ++j) {
                Vector selector(tuples.back().cols(), 0); selector[j] = 1;
                values.push(vector_json(tuples.back().apply(selector)));
            }
            item.set("optimal_parameter_tuples", values);
        }
        for (const auto *name : {"parameter_indices", "variable_indices", "inequality_indices"})
            factor.set(name, indices(block.at(name)));
        blocks.push(item); factor_blocks.push(factor);
    }
    if (feasible && out.at("witness_representation").text() == "explicit") {
        const auto kernel = matrix(keys.at("kernel")), constraints = matrix(keys.at("constraints")),
                   objective = matrix(keys.at("objective"));
        const auto offset = vector(keys.at("offset")), rhs = vector(keys.at("rhs")),
                   lower = vector(keys.at("declared_lower")), upper = vector(keys.at("declared_upper"));
        std::vector<std::size_t> selected(tuples.size(), 0);
        Integer emitted{0};
        for (;;) {
            Vector z(kernel.cols(), 0);
            for (std::size_t b = 0; b < tuples.size(); ++b)
                for (std::size_t j = 0; j < tuples[b].rows(); ++j)
                    z[parameter_indices[b][j]] = tuples[b].get(j, selected[b]);
            const auto x = rfc::numeric::add(kernel.apply(z), offset), values = constraints.apply(x);
            for (std::size_t i = 0; i < x.size(); ++i)
                require(x[i] >= lower[i] && x[i] <= upper[i] &&
                    (out.at("domain").text() != "bounded_integer" || x[i].get_den() == 1),
                    "retained allocation witness bounds/integrality");
            for (std::size_t i = 0; i < values.size(); ++i)
                require(satisfies(values[i], record.at("senses").at(i).text(), rhs[i]),
                        "retained source constraint witness residual");
            const auto optimum = record.at("minimum_or_maximum").rational();
            require(objective.apply(x).at(0) == optimum, "retained objective witness residual");
            auto witness = Json::object();
            witness.set("values", vector_json(x)); witness.set("constraint_values", vector_json(values));
            witness.set("objective", Json::string(optimum.get_str())); witnesses.push(witness); ++emitted;
            std::size_t k = selected.size();
            while (k > 0 && selected[k-1] + 1 == tuples[k-1].cols()) {selected[k-1] = 0; --k;}
            if (!k) break;
            ++selected[k-1];
        }
        require(emitted == record.at("optimal_points_or_vertices").integer(),
                "complete retained optimal Cartesian witness expansion");
    }
    out.set("bounds", bounds); out.set("optimizers", witnesses); out.set("block_results", blocks);
    if (record.at("factor_certificate").boolean()) {
        auto certificate = Json::object();
        certificate.set("declaration", constraint_request_view(record));
        certificate.set("consistent", Json::boolean(record.at("consistent").boolean()));
        for (const auto *name : {"equality_operator", "kernel", "reduced_inequalities", "grouping_contrasts", "group_indicators", "reduced_objective"})
            certificate.set(name, matrix_json(matrix(keys.at(name))));
        for (const auto *name : {"equality_rhs", "offset", "reduced_rhs"})
            certificate.set(name, keys.has(name) ? vector_json(vector(keys.at(name))) : Json::array());
        certificate.set("free_variable_indices", indices(keys.at("free_variable_indices")));
        certificate.set("reduced_senses", render_record_view(record.at("reduced_senses")));
        for (const auto *name : {"equality_origins", "inequality_origins"}) certificate.set(name, render_record_view(record.at(name)));
        certificate.set("blocks", factor_blocks);
        certificate.set("reconstruction", Json::string("x = offset + kernel z; z equals the original free-variable coordinates; each inequality and nonconstant variable belongs to exactly one block"));
        out.set("factor_certificate", certificate);
    }
    if (!record.at("consistent").boolean()) {
        auto contradiction = render_record_view(record.at("contradiction"));
        const auto kind = contradiction.at("kind").text();
        if (kind == "inconsistent_equalities") {
            const auto left = vector(keys.at("contradiction_left"));
            contradiction.set("left_null_vector", vector_json(left));
            contradiction.set("rhs_action", Json::string(rfc::numeric::dot(left, vector(keys.at("equality_rhs"))).get_str()));
        } else if (kind == "constant_inequality")
            contradiction.set("rhs", Json::string(vector(keys.at("reduced_rhs")).at(contradiction.at("reduced_row").integer()).get_str()));
        else if (kind == "nonintegral_fixed_coordinate")
            contradiction.set("value", Json::string(vector(keys.at("offset")).at(contradiction.at("variable_index").integer()).get_str()));
        out.set("contradiction", contradiction);
        if (out.has("factor_certificate")) out.at("factor_certificate").set("contradiction", contradiction);
    }
    return out;
}
} // namespace gadget::assay_detail
