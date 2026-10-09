#include "research_linear.hpp"
#include "research_fields.hpp"
#include "scientific_model.hpp"
#include <algorithm>
#include <gadget/resources.hpp>
#include <limits>
#include <rfc/machine/checked_size.hpp>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::research_detail {
namespace {
using Matrix = rfc::numeric::ExactMatrix;
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
Json strings(const std::vector<std::string> &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(Json::string(value));
    return out;
}
Json sparse_matrix(const Matrix &matrix) {
    auto out = Json::object(), entries = Json::array();
    for (const auto &[position, value] : matrix.entries()) {
        auto entry = Json::object();
        entry.set("row", cardinality(position.first));
        entry.set("column", cardinality(position.second));
        entry.set("value", rational(value));
        entries.push(entry);
    }
    out.set("rows", cardinality(matrix.rows()));
    out.set("columns", cardinality(matrix.cols()));
    out.set("entries", entries);
    return out;
}
struct Column {
    std::string id, role;
    std::vector<std::string> factors;
};
struct Variant {
    std::string id;
    std::set<std::string> remove;
};
struct Scope {
    std::string id;
    std::vector<std::string> selection;
    Indices rows;
};
struct Contrast {
    std::string id;
    Vector weights;
};
// The original column map and its complete coefficient fibre remain distinct
// from the independent frame used by the canonical prepared action.
struct Model : scientific_detail::PreparedObservedModel {
    std::string id;
    std::vector<std::size_t> columns;
    Model(std::string name, std::vector<std::size_t> selected, Matrix map,
          const rfc::field::ExactMetric &metric)
        : PreparedObservedModel(std::move(map), metric), id(std::move(name)),
          columns(std::move(selected)) {}
};
struct Fit {
    Vector coefficients, represented, residual;
    Rational quadrance;
    Json record;
};
Fit fit(const Model &model, const Vector &input, const rfc::field::ExactMetric &metric,
        const rfc::runtime::CpuTileStop &stop) {
    const auto native = model.apply(input, [&] { stop.check(); });
    Fit out{native.coefficients, native.represented, native.residual, native.quadrance, Json{}};
    const auto normal = model.design.transpose().apply(metric.apply(out.residual));
    out.record = Json::object();
    out.record.set("particular_coefficients", vector_json(out.coefficients));
    out.record.set("represented", vector_json(out.represented));
    out.record.set("residual", vector_json(out.residual));
    out.record.set("residual_quadrance", rational(out.quadrance));
    out.record.set("normal_equations", vector_json(normal));
    out.record.set("orthogonality_exact", Json::boolean(true));
    out.record.set("reconstruction_exact", Json::boolean(true));
    return out;
}
} // namespace

std::set<std::string> linear_model_field_ids(const Json &query) {
    std::set<std::string> out;
    for (const auto *member : {"responses", "controls", "selection"})
        for (const auto &id : names(query.at(member), query.at(member).size(), true))
            out.insert(id);
    const auto terms = query.at("terms");
    require(json_object_is_type(terms.raw(), json_type_array),
            "linear model terms must be an array");
    for (std::size_t j = 0; j < terms.size(); ++j) {
        const auto factors = terms.at(j).at("factors");
        require(json_object_is_type(factors.raw(), json_type_array),
                "term factors must be an array");
        for (std::size_t k = 0; k < factors.size(); ++k)
            out.insert(text(factors.at(k), source_byte_limit()));
    }
    if (query.has("strata"))
        for (std::size_t j = 0; j < query.at("strata").size(); ++j)
            for (const auto &id : names(query.at("strata").at(j).at("selection"),
                                        query.at("strata").at(j).at("selection").size(), true))
                out.insert(id);
    if (json_object_is_type(query.at("metric").raw(), json_type_object))
        out.insert(text(query.at("metric").at("field"), source_byte_limit()));
    require(!out.empty(), "linear model requires admitted fields");
    return out;
}

Json linear_model_comparison(const ScientificStudy &study, const Json &query, bool validate,
                             const std::shared_ptr<CpuResources> &resources) {
    members(query,
            {"operation", "responses", "terms", "controls", "selection", "support", "metric",
             "intercept"},
            {"strata", "ablations", "coefficient_contrasts", "workers"});
    require(query.at("support").text() == "joint_reported", "linear model joint reported support");
    const Inputs input(study, linear_model_field_ids(query), false, std::nullopt);
    const auto responses = names(query.at("responses"), input.fields.size()),
               controls = names(query.at("controls"), input.fields.size(), true),
               selection = names(query.at("selection"), input.fields.size(), true);
    std::vector<Column> columns;
    if (query.at("intercept").boolean())
        columns.push_back({"intercept", "intercept", {}});
    for (const auto &id : controls)
        columns.push_back({"control:" + id, "control", {id}});
    const auto baseline_columns = columns.size();
    std::set<std::string> column_ids, term_ids, design_fields;
    for (const auto &column : columns) {
        column_ids.insert(column.id);
        design_fields.insert(column.factors.begin(), column.factors.end());
    }
    const auto terms = query.at("terms");
    for (std::size_t j = 0; j < terms.size(); ++j) {
        const auto term = terms.at(j);
        members(term, {"id", "factors"});
        Column column{text(term.at("id"), source_byte_limit()), "term", {}};
        require(column_ids.insert(column.id).second, "duplicate linear model column identity");
        term_ids.insert(column.id);
        const auto factors = term.at("factors");
        require(factors.size() > 0, "linear model term needs explicit factors");
        for (std::size_t k = 0; k < factors.size(); ++k) {
            const auto id = text(factors.at(k), source_byte_limit());
            column.factors.push_back(id);
            design_fields.insert(id);
        }
        columns.push_back(std::move(column));
    }
    for (const auto &response : responses)
        require(!design_fields.contains(response), "a response cannot be its own model input");
    std::vector<std::string> complete(design_fields.begin(), design_fields.end());
    complete.insert(complete.end(), responses.begin(), responses.end());
    std::optional<std::string> metric_field;
    const auto metric_request = query.at("metric");
    if (json_object_is_type(metric_request.raw(), json_type_object)) {
        members(metric_request, {"kind", "field", "meaning"});
        require(metric_request.at("kind").text() == "diagonal_field", "linear model metric kind");
        metric_field = text(metric_request.at("field"), source_byte_limit());
        (void)text(metric_request.at("meaning"), source_byte_limit());
        require(input.meaning.at(*metric_field).at("units").text() == "dimensionless",
                "linear model metric weights must be dimensionless");
        complete.push_back(*metric_field);
    } else
        require(metric_request.text() == "identity", "linear model metric must be declared");
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        query.has("workers")
            ? count(query.at("workers"), cpu->coordinator().capacity(), "model workers")
            : 1;
    require(workers > 0, "linear model workers must be positive");
    std::vector<Variant> variants{{"controls", term_ids}, {"full", {}}};
    std::set<std::string> variant_ids{"controls", "full"};
    if (query.has("ablations")) {
        const auto declarations = query.at("ablations");
        require(json_object_is_type(declarations.raw(), json_type_array),
                "model ablations must be an array");
        for (std::size_t j = 0; j < declarations.size(); ++j) {
            const auto declared = declarations.at(j);
            members(declared, {"id", "remove"});
            Variant variant{text(declared.at("id"), source_byte_limit()), {}};
            require(variant_ids.insert(variant.id).second, "duplicate model variant identity");
            for (const auto &id : names(declared.at("remove"), term_ids.size())) {
                require(term_ids.contains(id), "ablation must name a declared term");
                variant.remove.insert(id);
            }
            variants.push_back(std::move(variant));
        }
    }
    std::vector<Contrast> contrasts;
    std::set<std::string> contrast_ids;
    if (query.has("coefficient_contrasts")) {
        const auto declarations = query.at("coefficient_contrasts");
        require(json_object_is_type(declarations.raw(), json_type_array),
                "coefficient contrasts must be an array");
        for (std::size_t j = 0; j < declarations.size(); ++j) {
            const auto declared = declarations.at(j);
            members(declared, {"id", "weights"});
            Contrast contrast{text(declared.at("id"), source_byte_limit()),
                              Vector(columns.size(), Rational{0})};
            require(contrast_ids.insert(contrast.id).second,
                    "duplicate coefficient contrast identity");
            const auto weights = declared.at("weights");
            require(json_object_is_type(weights.raw(), json_type_array) && weights.size(),
                    "coefficient contrast needs explicit weights");
            std::set<std::string> used;
            for (std::size_t k = 0; k < weights.size(); ++k) {
                const auto weight = weights.at(k);
                members(weight, {"column", "value"});
                const auto id = text(weight.at("column"), source_byte_limit());
                require(column_ids.contains(id) && used.insert(id).second,
                        "unknown or repeated contrast column");
                const auto position = std::find_if(columns.begin(), columns.end(),
                                                   [&](const auto &c) { return c.id == id; });
                contrast.weights[static_cast<std::size_t>(position - columns.begin())] =
                    exact_rational(text(weight.at("value"), source_byte_limit()));
            }
            contrasts.push_back(std::move(contrast));
        }
    }
    std::vector<Scope> scopes{{"all", selection, {}}};
    std::set<std::string> scope_ids{"all"};
    if (query.has("strata")) {
        const auto declarations = query.at("strata");
        require(json_object_is_type(declarations.raw(), json_type_array),
                "model strata must be an array");
        for (std::size_t j = 0; j < declarations.size(); ++j) {
            const auto declared = declarations.at(j);
            members(declared, {"id", "selection"});
            Scope scope{text(declared.at("id"), source_byte_limit()), selection, {}};
            require(scope_ids.insert(scope.id).second, "duplicate model stratum identity");
            for (const auto &id : names(declared.at("selection"), input.fields.size(), true))
                if (std::find(scope.selection.begin(), scope.selection.end(), id) ==
                    scope.selection.end())
                    scope.selection.push_back(id);
            scopes.push_back(std::move(scope));
        }
    }
    for (auto &scope : scopes) {
        scope.rows = input.select(complete, scope.selection, true);
        if (metric_field)
            for (const auto row : scope.rows)
                require(input.value(*metric_field, row) > 0,
                        "model metric diagonal must be positive");
    }
    if (validate)
        return Json{};
    auto reservation = cpu->reserve(workers);
    auto result = Json::object(), declarations = Json::array(), strata = Json::array();
    result.set("schema", Json::string("gadget.linear.model.comparison.v1"));
    result.set("native_carrier", Json::string(std::to_string(*input.carrier)));
    result.set("grade", cardinality(input.grade.value));
    result.set("arithmetic", Json::string("exact_rational"));
    result.set("coefficient_rule", Json::string("beta = particular + kernel * free; only "
                                                "kernel-annihilating contrasts are identified"));
    for (const auto &column : columns) {
        auto declaration = Json::object(), units = Json::array();
        declaration.set("id", Json::string(column.id));
        declaration.set("role", Json::string(column.role));
        declaration.set("factors", strings(column.factors));
        for (const auto &id : column.factors)
            units.push(input.meaning.at(id).at("units"));
        declaration.set("unit_factors", units);
        declarations.push(declaration);
    }
    result.set("columns", declarations);
    for (const auto &scope : scopes) {
        auto stratum = Json::object(), coordinates = Json::array(), exclusions = Json::array();
        stratum.set("id", Json::string(scope.id));
        stratum.set("selection", strings(scope.selection));
        stratum.set("entities", input.identity(scope.rows));
        for (const auto row : scope.rows)
            coordinates.push(cardinality(row));
        stratum.set("coordinates", coordinates);
        stratum.set("n", cardinality(scope.rows.size()));
        std::set<std::size_t> retained(scope.rows.begin(), scope.rows.end());
        for (std::size_t row = 0; row < input.entities.size(); ++row) {
            if (retained.contains(row))
                continue;
            auto exclusion = Json::object(), missing = Json::array(), selectors = Json::array();
            for (const auto &id : complete)
                if (!input.fields.at(id)->tensor.support().contains(row))
                    missing.push(Json::string(id));
            for (const auto &id : scope.selection)
                if (!input.fields.at(id)->tensor.support().contains(row) ||
                    input.value(id, row) != 1)
                    selectors.push(Json::string(id));
            exclusion.set("entity", Json::string(input.entities[row]));
            exclusion.set("coordinate", cardinality(row));
            exclusion.set("unobserved_fields", missing);
            exclusion.set("unselected_indicators", selectors);
            exclusions.push(exclusion);
        }
        stratum.set("excluded", exclusions);
        if (scope.rows.empty()) {
            stratum.set("status", Json::string("no_joint_observations"));
            strata.push(stratum);
            continue;
        }
        stratum.set("status", Json::string("completed"));
        Vector weights(scope.rows.size(), Rational{1});
        if (metric_field)
            weights = input.values(*metric_field, scope.rows);
        const auto metric = field::ExactMetric::diagonal(weights);
        stratum.set("metric_diagonal", vector_json(weights));
        Matrix all(scope.rows.size(), columns.size());
        for (std::size_t j = 0; j < columns.size(); ++j)
            for (std::size_t i = 0; i < scope.rows.size(); ++i) {
                Rational value{1};
                for (const auto &id : columns[j].factors)
                    value *= input.value(id, scope.rows[i]);
                all.set(i, j, value);
            }
        auto models = Json::array();
        std::vector<std::vector<Fit>> fits;
        for (const auto &variant : variants) {
            std::vector<std::size_t> selected;
            for (std::size_t j = 0; j < columns.size(); ++j)
                if (j < baseline_columns || !variant.remove.contains(columns[j].id))
                    selected.push_back(j);
            Matrix map(scope.rows.size(), selected.size());
            for (std::size_t j = 0; j < selected.size(); ++j)
                for (std::size_t i = 0; i < scope.rows.size(); ++i)
                    map.set(i, j, all.get(i, selected[j]));
            const Model model(variant.id, selected, std::move(map), metric);
            fits.emplace_back(responses.size());
            auto &computed = fits.back();
            runtime::CpuTileExecutor executor(
                {workers, cpu->coordinator().capacity(), responses.size()});
            const auto report =
                executor.run_tiles(responses.size(), [&](std::size_t j, std::size_t,
                                                         const runtime::CpuTileStop &stop) {
                    computed[j] = fit(model, input.values(responses[j], scope.rows), metric, stop);
                });
            auto record = Json::object(), axis = Json::array(), identified = Json::array(),
                 response_records = Json::array();
            for (std::size_t j = 0; j < selected.size(); ++j) {
                axis.push(Json::string(columns[selected[j]].id));
                bool unique = true;
                for (std::size_t k = 0; k < model.kernel.cols(); ++k)
                    unique = unique && model.kernel.get(j, k) == 0;
                identified.push(Json::boolean(unique));
            }
            record.set("id", Json::string(model.id));
            record.set("columns", axis);
            record.set("rank", cardinality(model.projection.frame.cols()));
            record.set("residual_dimension",
                       cardinality(scope.rows.size() - model.projection.frame.cols()));
            record.set("coefficient_nullity", cardinality(model.kernel.cols()));
            record.set("coefficient_identified", identified);
            record.set("design", sparse_matrix(model.design));
            record.set("independent_frame", sparse_matrix(model.projection.frame));
            record.set("coefficient_kernel", sparse_matrix(model.kernel));
            for (std::size_t j = 0; j < responses.size(); ++j) {
                auto response = computed[j].record;
                response.set("field", Json::string(responses[j]));
                response.set("meaning", input.meaning.at(responses[j]));
                response.set("input", vector_json(input.values(responses[j], scope.rows)));
                if (model.id == "full") {
                    auto readings = Json::array();
                    for (const auto &contrast : contrasts) {
                        const auto ambiguity = model.kernel.transpose().apply(contrast.weights);
                        const bool unique = numeric::is_zero(ambiguity);
                        auto reading = Json::object();
                        reading.set("id", Json::string(contrast.id));
                        reading.set("weights", vector_json(contrast.weights));
                        reading.set("kernel_pairings", vector_json(ambiguity));
                        reading.set("identified", Json::boolean(unique));
                        reading.set("value", unique
                                                 ? rational(numeric::dot(contrast.weights,
                                                                         computed[j].coefficients))
                                                 : Json{});
                        readings.push(reading);
                    }
                    response.set("coefficient_contrasts", readings);
                }
                response_records.push(response);
            }
            record.set("responses", response_records);
            auto execution = Json::object();
            execution.set("owner", Json::string("rfc::numeric::PreparedExactSolve"));
            execution.set("gram_preparations", cardinality(1));
            execution.set("applications", cardinality(responses.size()));
            execution.set("factor_entries", cardinality(model.projection.solve.stored_entries()));
            execution.set("requested_workers", cardinality(report.requested_workers));
            execution.set("worker_slots", cardinality(report.worker_slots));
            execution.set("completed_tiles", cardinality(report.completed_tiles));
            record.set("execution", execution);
            models.push(record);
        }
        stratum.set("models", models);
        auto comparisons = Json::array();
        for (std::size_t m = 0; m < variants.size(); ++m) {
            if (m == 1)
                continue;
            auto comparison = Json::object(), readings = Json::array();
            comparison.set("reduced", Json::string(variants[m].id));
            comparison.set("full", Json::string("full"));
            for (std::size_t j = 0; j < responses.size(); ++j) {
                const auto contribution =
                    numeric::sub(fits[1][j].represented, fits[m][j].represented);
                const Rational gain = fits[m][j].quadrance - fits[1][j].quadrance;
                require(gain >= 0 && gain == metric.inner(contribution, contribution) &&
                            metric.inner(contribution, fits[1][j].residual) == 0,
                        "nested linear model quadrance identity");
                auto reading = Json::object();
                reading.set("field", Json::string(responses[j]));
                reading.set("added_representation", vector_json(contribution));
                reading.set("added_quadrance", rational(gain));
                reading.set("reduced_residual_quadrance", rational(fits[m][j].quadrance));
                reading.set("fraction_of_reduced_residual",
                            fits[m][j].quadrance == 0
                                ? Json{}
                                : rational(Rational{gain / fits[m][j].quadrance}));
                reading.set("nested_identity_exact", Json::boolean(true));
                readings.push(reading);
            }
            comparison.set("responses", readings);
            comparisons.push(comparison);
        }
        stratum.set("comparisons", comparisons);
        strata.push(stratum);
    }
    result.set("strata", strata);
    result.set(
        "definition",
        Json::string("For each declared stratum, minimize (y-X beta)* M (y-X beta) on fixed joint "
                     "observed support. Retain Pi_image(X) y, its orthogonal residual, the "
                     "complete coefficient fibre, and full-versus-reduced field differences. "
                     "Indicator responses describe supplied memberships."));
    return result;
}
} // namespace gadget::research_detail
