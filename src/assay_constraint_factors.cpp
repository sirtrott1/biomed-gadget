#include "allocation_detail.hpp"
#include <rfc/core/builder.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <limits>

namespace gadget::assay_detail::allocation_detail {
Json vector_json(const Vector &values) {
    auto out = Json::array();
    for (const auto &q : values)
        out.push(Json::string(q.get_str()));
    return out;
}
Json matrix_json(const ExactMatrix &matrix) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(matrix.rows()));
    out.set("columns", cardinality(matrix.cols()));
    for (const auto &[ij, q] : matrix.entries()) {
        auto entry = Json::array();
        entry.push(cardinality(ij.first));
        entry.push(cardinality(ij.second));
        entry.push(Json::string(q.get_str()));
        entries.push(entry);
    }
    out.set("entries", entries);
    return out;
}
bool satisfies(const Rational &value, const std::string &sense, const Rational &rhs) {
    return sense == "eq" ? value == rhs : sense == "le" ? value <= rhs : value >= rhs;
}
Problem::Problem(const Json &request) : declaration(request) {
    members(request, {"schema", "operation", "domain", "variables", "constraints", "objective",
                      "units", "evidence"}, {"factor_certificate", "witness_representation"});
    domain = text(request.at("domain"));
    require(domain == "bounded_integer" || domain == "bounded_rational", "allocation domain");
    (void)text(request.at("units"), source_byte_limit());
    (void)text(request.at("evidence"), source_byte_limit());
    if (request.has("factor_certificate"))
        (void)request.at("factor_certificate").boolean();
    if (request.has("witness_representation")) {
        const auto mode = text(request.at("witness_representation"));
        require(mode == "explicit" || mode == "factorized", "allocation witness representation");
    }
    const auto variables = request.at("variables"), rows = request.at("constraints");
    const auto n = variables.size(), m = rows.size();
    require(n > 0 && n <= std::numeric_limits<std::uint32_t>::max() &&
                m <= std::numeric_limits<std::uint32_t>::max(),
            "allocation native cell coordinate representation");
    std::set<std::string> ids, row_ids;
    for (std::size_t i = 0; i < n; ++i) {
        const auto v = variables.at(i);
        members(v, {"id", "lower", "upper"});
        require(ids.insert(text(v.at("id"), source_byte_limit())).second,
                "duplicate allocation variable");
        names.push(v.at("id"));
        lower.push_back(exact_rational(text(v.at("lower"), source_byte_limit())));
        upper.push_back(exact_rational(text(v.at("upper"), source_byte_limit())));
        require(lower.back() <= upper.back(), "allocation variable bounds");
        if (domain == "bounded_integer") {
            require(lower.back().get_den() == 1 && upper.back().get_den() == 1,
                    "integer allocation requires integer bounds");
            const Rational width = upper.back() - lower.back() + 1;
            declared_box *= width.get_num();
        }
    }
    const auto coefficients = [&](const Json &values) {
        require(values.size() == n, "allocation coefficient width");
        Vector out;
        for (std::size_t i = 0; i < n; ++i)
            out.push_back(exact_rational(text(values.at(i), source_byte_limit())));
        return out;
    };
    constraints = ExactMatrix(m, n);
    for (std::size_t i = 0; i < m; ++i) {
        const auto row = rows.at(i);
        members(row, {"id", "coefficients", "sense", "rhs"});
        require(row_ids.insert(text(row.at("id"), source_byte_limit())).second,
                "duplicate allocation constraint");
        senses.push_back(text(row.at("sense")));
        require(senses.back() == "eq" || senses.back() == "le" || senses.back() == "ge",
                "allocation constraint sense");
        rhs.push_back(exact_rational(text(row.at("rhs"), source_byte_limit())));
        const auto row_values = coefficients(row.at("coefficients"));
        for (std::size_t j = 0; j < n; ++j)
            constraints.set(i, j, row_values[j]);
    }
    const auto obj = request.at("objective");
    members(obj, {"sense", "coefficients"});
    direction = text(obj.at("sense"));
    require(direction == "minimize" || direction == "maximize", "allocation objective sense");
    objective = ExactMatrix(1, n);
    const auto values = coefficients(obj.at("coefficients"));
    for (std::size_t j = 0; j < n; ++j)
        objective.set(0, j, values[j]);
}

void Problem::validate() const {
    require(domain == "bounded_integer" || domain == "bounded_rational", "allocation domain");
    require(direction == "minimize" || direction == "maximize", "allocation objective sense");
    const auto n = lower.size(), m = rhs.size();
    require(n > 0 && n <= std::numeric_limits<std::uint32_t>::max() &&
                m <= std::numeric_limits<std::uint32_t>::max(),
            "allocation native cell coordinate representation");
    require(upper.size() == n && names.size() == n && constraints.cols() == n &&
                constraints.rows() == m && senses.size() == m && objective.rows() == 1 &&
                objective.cols() == n, "typed allocation source shape");
    Integer box{1};
    std::set<std::string> ids;
    for (std::size_t i = 0; i < n; ++i) {
        require(ids.insert(text(names.at(i), source_byte_limit())).second,
                "duplicate allocation variable");
        require(lower[i] <= upper[i], "allocation variable bounds");
        if (domain == "bounded_integer") {
            require(lower[i].get_den() == 1 && upper[i].get_den() == 1,
                    "integer allocation requires integer bounds");
            const Rational width = upper[i] - lower[i] + 1;
            box *= width.get_num();
        }
    }
    require(domain != "bounded_integer" || declared_box == box, "typed allocation box count");
    for (const auto &sense : senses)
        require(sense == "eq" || sense == "le" || sense == "ge", "allocation constraint sense");
    (void)text(declaration.at("units"), source_byte_limit());
    (void)text(declaration.at("evidence"), source_byte_limit());
}
namespace {
Json origin(const std::string &kind, std::size_t index) {
    auto out = Json::object();
    out.set("kind", Json::string(kind));
    out.set("source_index", cardinality(index));
    return out;
}
Json index_json(const std::vector<std::size_t> &indices) {
    auto out = Json::array();
    for (auto i : indices)
        out.push(cardinality(i));
    return out;
}
Json source_carrier(const Problem &p, std::optional<rfc::core::RelationalComplex> &retained) {
    rfc::core::RelationSignature signature;
    std::map<std::string, rfc::core::RelationTypeId> types;
    for (const auto *sense : {"eq", "le", "ge"})
        types.emplace(sense, signature.add_type(std::string("declared.allocation.linear.") + sense));
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(rfc::core::Grade{1});
    const auto n = p.constraints.cols(), m = p.constraints.rows();
    for (std::size_t j = 0; j < n; ++j)
        builder.declare_witness_identity(j + 1);
    std::size_t occurrences = 0;
    for (std::size_t i = 0; i < m; ++i) {
        std::vector<rfc::core::BoundaryOccurrenceInput> slots;
        for (std::size_t j = 0; j < n; ++j) {
            const auto q = p.constraints.get(i, j);
            occurrences = rfc::machine::checked_size_add(occurrences, 1);
            // Dense source rows declare zero coefficients too: they remain
            // existent occurrences, distinct from absent numerical support.
            slots.push_back({j + 1, true, q < 0, q < 0 ? -q : q, {},
                             rfc::core::OccurrenceKey{occurrences}});
        }
        builder.add_grade_one(rfc::core::RelationKey{i + 1}, types.at(p.senses[i]), slots);
    }
    const auto complex = builder.freeze();
    retained.emplace(complex);
    if (m)
        require(rfc::field::materialize_boundary_matrix(complex, {1}) == p.constraints.transpose(),
                "native allocation source boundary residual");
    const auto bytes = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(complex), 1, source_codec_limits(source_byte_limit())));
    auto out = Json::object();
    out.set("native_complex_sha256", Json::string("sha256:" + sha256(std::string(
        reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    out.set("variable_count", cardinality(n));
    out.set("relation_count", cardinality(m));
    out.set("occurrence_count", cardinality(occurrences));
    out.set("coordinates", Json::string("Grade-zero variables and grade-one authored constraints in source order; boundary equals A transpose"));
    return out;
}
} // namespace

Factors::Factors(const Problem &p) {
    carrier = source_carrier(p, source);
    const auto n = p.constraints.cols(), m = p.constraints.rows();
    std::vector<std::size_t> eq_rows, fixed;
    for (std::size_t i = 0; i < m; ++i)
        if (p.senses[i] == "eq") {
            eq_rows.push_back(i);
            equality_origins.push(origin("source_constraint", i));
            equality_rhs.push_back(p.rhs[i]);
        }
    for (std::size_t j = 0; j < n; ++j)
        if (p.lower[j] == p.upper[j]) {
            fixed.push_back(j);
            equality_origins.push(origin("fixed_bound", j));
            equality_rhs.push_back(p.lower[j]);
        }
    equalities = ExactMatrix(equality_rhs.size(), n);
    std::map<std::size_t, std::size_t> eq_index;
    for (std::size_t i = 0; i < eq_rows.size(); ++i)
        eq_index.emplace(eq_rows[i], i);
    for (const auto &[ij, q] : p.constraints.entries())
        if (eq_index.contains(ij.first))
            equalities.set(eq_index.at(ij.first), ij.second, q);
    for (std::size_t i = 0; i < fixed.size(); ++i)
        equalities.set(eq_rows.size() + i, fixed[i], 1);
    ExactMatrix augmented(equalities.rows(), rfc::machine::checked_size_add(n, 1));
    for (const auto &[ij, q] : equalities.entries())
        augmented.set(ij.first, ij.second, q);
    for (std::size_t i = 0; i < equality_rhs.size(); ++i)
        augmented.set(i, n, equality_rhs[i]);
    std::vector<std::size_t> pivots;
    const auto reduced = augmented.rref(pivots);
    if (!pivots.empty() && pivots.back() == n) {
        consistent = false;
        equality_rank = pivots.size() - 1;
        const auto left = equalities.transpose().kernel_frame();
        ExactMatrix b(1, equality_rhs.size());
        for (std::size_t i = 0; i < equality_rhs.size(); ++i)
            b.set(0, i, equality_rhs[i]);
        for (std::size_t j = 0; j < left.cols(); ++j) {
            Vector unit(left.cols(), Rational{0});
            unit[j] = 1;
            const auto lambda = left.apply(unit);
            const auto value = b.apply(lambda).at(0);
            if (value != 0) {
                require(equalities.transpose().apply(lambda) == Vector(n, Rational{0}),
                        "native left-null contradiction residual");
                contradiction.set("kind", Json::string("inconsistent_equalities"));
                contradiction.set("left_null_vector", vector_json(lambda));
                contradiction_left = lambda;
                contradiction.set("rhs_action", Json::string(value.get_str()));
                return;
            }
        }
        require(false, "missing native equality contradiction");
    }
    equality_rank = pivots.size();
    offset = Vector(n, Rational{0});
    std::vector<bool> pivot(n, false);
    for (std::size_t i = 0; i < pivots.size(); ++i) {
        pivot[pivots[i]] = true;
        offset[pivots[i]] = reduced.get(i, n);
    }
    for (std::size_t j = 0; j < n; ++j)
        if (!pivot[j])
            free_variables.push_back(j);
    kernel = equalities.kernel_frame();
    const auto d = kernel.cols();
    require(d == free_variables.size() && equalities.apply(offset) == equality_rhs &&
                equalities.multiply(kernel).nonzeros() == 0,
            "native affine equality residual");
    std::vector<std::size_t> free_index(n, d);
    for (std::size_t j = 0; j < d; ++j) {
        free_index[free_variables[j]] = j;
        require(offset[free_variables[j]] == 0, "native free-coordinate offset");
    }
    std::size_t free_entries = 0;
    for (const auto &[ij, q] : kernel.entries())
        if (free_index[ij.first] != d) {
            require(ij.second == free_index[ij.first] && q == 1,
                    "native source free-coordinate identity");
            ++free_entries;
        }
    require(free_entries == d, "complete native free-coordinate identity");
    const auto action = p.constraints.multiply(kernel);
    const auto shift = p.constraints.apply(offset);
    const auto row_count = rfc::machine::checked_size_add(m - eq_rows.size(),
        rfc::machine::checked_size_mul(n, 2));
    inequalities = ExactMatrix(row_count, d);
    std::map<std::size_t, std::size_t> source_index;
    for (std::size_t i = 0; i < m; ++i)
        if (p.senses[i] != "eq") {
            source_index.emplace(i, rhs.size());
            inequality_origins.push(origin("source_constraint", i));
            senses.push_back(p.senses[i]);
            rhs.push_back(p.rhs[i] - shift[i]);
        }
    for (const auto &[ij, q] : action.entries())
        if (source_index.contains(ij.first))
            inequalities.set(source_index.at(ij.first), ij.second, q);
    const auto bound_begin = rhs.size();
    for (std::size_t j = 0; j < n; ++j) {
        inequality_origins.push(origin("lower_bound", j));
        inequality_origins.push(origin("upper_bound", j));
        senses.push_back("ge");
        senses.push_back("le");
        rhs.push_back(p.lower[j] - offset[j]);
        rhs.push_back(p.upper[j] - offset[j]);
    }
    for (const auto &[ij, q] : kernel.entries()) {
        inequalities.set(bound_begin + 2 * ij.first, ij.second, q);
        inequalities.set(bound_begin + 2 * ij.first + 1, ij.second, q);
    }
    objective = p.objective.multiply(kernel);
    std::vector<std::vector<std::size_t>> supports(rhs.size()), variable_support(n);
    for (const auto &[ij, q] : inequalities.entries()) {
        (void)q;
        supports[ij.first].push_back(ij.second);
    }
    for (const auto &[ij, q] : kernel.entries()) {
        (void)q;
        variable_support[ij.first].push_back(ij.second);
    }
    for (std::size_t i = 0; i < rhs.size(); ++i)
        if (supports[i].empty() && !satisfies(Rational{0}, senses[i], rhs[i])) {
            consistent = false;
            contradiction.set("kind", Json::string("constant_inequality"));
            contradiction.set("reduced_row", cardinality(i));
            contradiction.set("origin", inequality_origins.at(i));
            contradiction.set("sense", Json::string(senses[i]));
            contradiction.set("rhs", Json::string(rhs[i].get_str()));
            return;
        }
    if (p.domain == "bounded_integer")
        for (std::size_t j = 0; j < n; ++j)
            if (variable_support[j].empty() && offset[j].get_den() != 1) {
                consistent = false;
                contradiction.set("kind", Json::string("nonintegral_fixed_coordinate"));
                contradiction.set("variable_index", cardinality(j));
                contradiction.set("value", Json::string(offset[j].get_str()));
                return;
            }
    // Auxiliary pair contrasts encode support connectivity only. Their native
    // kernel has one constant indicator per independent parameter component.
    // These are not additional constraints on x or biological higher faces.
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    for (const auto &support : supports)
        for (std::size_t k = 1; k < support.size(); ++k)
            pairs.emplace(support.front(), support[k]);
    grouping = ExactMatrix(pairs.size(), d);
    std::size_t row = 0;
    for (const auto &[a, b] : pairs) {
        grouping.set(row, a, 1);
        grouping.set(row++, b, -1);
    }
    indicator = grouping.kernel_frame();
    require(grouping.multiply(indicator).nonzeros() == 0, "native grouping kernel residual");
    blocks.resize(indicator.cols());
    std::vector<std::size_t> groups(d, blocks.size());
    for (const auto &[ij, q] : indicator.entries()) {
        require(q == 1 && groups[ij.first] == blocks.size(), "native disjoint parameter indicators");
        groups[ij.first] = ij.second;
        blocks[ij.second].parameters.push_back(ij.first);
    }
    for (auto g : groups)
        require(g < blocks.size(), "native complete parameter indicators");
    for (std::size_t i = 0; i < supports.size(); ++i)
        if (!supports[i].empty()) {
            const auto g = groups[supports[i].front()];
            for (auto j : supports[i])
                require(groups[j] == g, "inequality spans independent factors");
            blocks[g].inequalities.push_back(i);
        }
    for (std::size_t i = 0; i < variable_support.size(); ++i)
        if (!variable_support[i].empty()) {
            const auto g = groups[variable_support[i].front()];
            for (auto j : variable_support[i])
                require(groups[j] == g, "source coordinate spans independent factors");
            blocks[g].variables.push_back(i);
        }
    // Pack declared restrictions once in source order. Repeated selection of
    // the entire sparse operator for each independent block would scan unrelated
    // coordinates and duplicate preparation work.
    std::vector<std::size_t> parameter_local(d), variable_local(n), row_local(rhs.size());
    for (auto &block : blocks) {
        block.embedding = ExactMatrix(block.variables.size(), block.parameters.size());
        block.normals = ExactMatrix(block.inequalities.size(), block.parameters.size());
        block.objective = ExactMatrix(1, block.parameters.size());
        for (std::size_t j = 0; j < block.parameters.size(); ++j) {
            const auto p_index = block.parameters[j];
            parameter_local[p_index] = j;
            block.lower.push_back(p.lower[free_variables[p_index]]);
            block.upper.push_back(p.upper[free_variables[p_index]]);
        }
        for (std::size_t i = 0; i < block.variables.size(); ++i) {
            variable_local[block.variables[i]] = i;
            block.offset.push_back(offset[block.variables[i]]);
        }
        for (std::size_t i = 0; i < block.inequalities.size(); ++i) {
            const auto source_row = block.inequalities[i];
            row_local[source_row] = i;
            block.rhs.push_back(rhs[source_row]);
            block.senses.push_back(senses[source_row]);
        }
    }
    for (const auto &[ij, q] : kernel.entries())
        blocks[groups[ij.second]].embedding.set(variable_local[ij.first], parameter_local[ij.second], q);
    for (const auto &[ij, q] : inequalities.entries())
        blocks[groups[ij.second]].normals.set(row_local[ij.first], parameter_local[ij.second], q);
    for (const auto &[ij, q] : objective.entries())
        blocks[groups[ij.second]].objective.set(0, parameter_local[ij.second], q);
}

Json Factors::certificate(const Problem &p) const {
    auto out = Json::object(), blocks_json = Json::array();
    out.set("declaration", p.declaration);
    out.set("consistent", Json::boolean(consistent));
    out.set("equality_operator", matrix_json(equalities));
    out.set("equality_rhs", vector_json(equality_rhs));
    out.set("equality_origins", equality_origins);
    out.set("offset", vector_json(offset));
    out.set("kernel", matrix_json(kernel));
    out.set("free_variable_indices", index_json(free_variables));
    out.set("reduced_inequalities", matrix_json(inequalities));
    out.set("reduced_rhs", vector_json(rhs));
    auto senses_json = Json::array();
    for (const auto &sense : senses)
        senses_json.push(Json::string(sense));
    out.set("reduced_senses", senses_json);
    out.set("inequality_origins", inequality_origins);
    out.set("grouping_contrasts", matrix_json(grouping));
    out.set("group_indicators", matrix_json(indicator));
    out.set("reduced_objective", matrix_json(objective));
    for (const auto &block : blocks) {
        auto item = Json::object();
        item.set("parameter_indices", index_json(block.parameters));
        item.set("variable_indices", index_json(block.variables));
        item.set("inequality_indices", index_json(block.inequalities));
        blocks_json.push(item);
    }
    out.set("blocks", blocks_json);
    out.set("reconstruction", Json::string("x = offset + kernel z; z equals the original free-variable coordinates; each inequality and nonconstant variable belongs to exactly one block"));
    if (!consistent)
        out.set("contradiction", contradiction);
    return out;
}
} // namespace gadget::assay_detail::allocation_detail
