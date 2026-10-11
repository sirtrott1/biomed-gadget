#include "research_linear.hpp"
#include "research_fields.hpp"
#include "scientific_result_fields.hpp"
#include "scientific_model.hpp"
#include <algorithm>
#include <gadget/resources.hpp>
#include <limits>
#include <rfc/machine/checked_size.hpp>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::research_detail {
namespace {
using Matrix = rfc::numeric::ExactMatrix;
using Value = scientific_detail::ResultValue;
Value strings(const std::vector<std::string> &values) {
    auto out = Value::array();
    for (const auto &value : values) out.push(Value::string(value));
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
    Vector normal;
};
Fit fit(const Model &model, const Vector &input, const rfc::field::ExactMetric &metric,
        const rfc::runtime::CpuTileStop &stop) {
    const auto native = model.apply(input, [&] { stop.check(); });
    Fit out{native.coefficients, native.represented, native.residual, native.quadrance,
            model.design.transpose().apply(metric.apply(native.residual))};
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

std::optional<LinearModelFields> evaluate_linear_model(const ScientificStudy &study, const Json &query, bool validate,
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
        return std::nullopt;
    auto reservation = cpu->reserve(workers);
    scientific_detail::ResultFields publication(study.native_publication());
    const auto full_columns = publication.axis(columns.size());
    std::map<LinearModelAxis, std::optional<Rational>> comparisons_by_axis;
    const auto rational = [&](const Rational &value) {return publication.scalar(value);};
    auto result = Value::object(), declarations = Value::array(), strata = Value::array();
    result.set("schema", Value::string("gadget.linear.model.comparison.v1"));
    result.set("native_carrier", Value::string(std::to_string(*input.carrier)));
    result.set("grade", Value::count(input.grade.value));
    result.set("arithmetic", Value::string("exact_rational"));
    result.set("coefficient_rule", Value::string("beta = particular + kernel * free; only "
                                                "kernel-annihilating contrasts are identified"));
    for (const auto &column : columns) {
        auto declaration = Value::object(), units = Value::array();
        declaration.set("id", Value::string(column.id));
        declaration.set("role", Value::string(column.role));
        declaration.set("factors", strings(column.factors));
        for (const auto &id : column.factors)
            units.push(input.meaning.at(id).at("units"));
        declaration.set("unit_factors", units);
        declarations.push(declaration);
    }
    result.set("columns", declarations);
    for (const auto &scope : scopes) {
        auto stratum = Value::object(), coordinates = Value::array(), exclusions = Value::array();
        stratum.set("id", Value::string(scope.id));
        stratum.set("selection", strings(scope.selection));
        stratum.set("entities", input.identity(scope.rows));
        for (const auto row : scope.rows)
            coordinates.push(Value::count(row));
        stratum.set("coordinates", coordinates);
        stratum.set("n", Value::count(scope.rows.size()));
        std::set<std::size_t> retained(scope.rows.begin(), scope.rows.end());
        for (std::size_t row = 0; row < input.entities.size(); ++row) {
            if (retained.contains(row))
                continue;
            auto exclusion = Value::object(), missing = Value::array(), selectors = Value::array();
            for (const auto &id : complete)
                if (!input.fields.at(id)->tensor.support().contains(row))
                    missing.push(Value::string(id));
            for (const auto &id : scope.selection)
                if (!input.fields.at(id)->tensor.support().contains(row) ||
                    input.value(id, row) != 1)
                    selectors.push(Value::string(id));
            exclusion.set("entity", Value::string(input.entities[row]));
            exclusion.set("coordinate", Value::count(row));
            exclusion.set("unobserved_fields", missing);
            exclusion.set("unselected_indicators", selectors);
            exclusions.push(exclusion);
        }
        stratum.set("excluded", exclusions);
        if (scope.rows.empty()) {
            stratum.set("status", Value::string("no_joint_observations"));
            strata.push(stratum);
            continue;
        }
        stratum.set("status", Value::string("completed"));
        const auto row_carrier = publication.axis(scope.rows.size());
        const auto row_coordinates = publication.coordinates(scope.rows, row_carrier);
        Matrix row_lift(input.entities.size(), scope.rows.size());
        for (std::size_t i = 0; i < scope.rows.size(); ++i) row_lift.set(scope.rows[i], i, 1);
        auto scope_binding = Value::object();
        scope_binding.set("stratum", Value::string(scope.id));
        scope_binding.set("source_carrier", Value::count(*input.carrier));
        scope_binding.set("row_carrier", Value::count(row_carrier));
        scope_binding.set("row_coordinates", Value(row_coordinates));
        scope_binding.set("source_lift", publication.matrix(row_lift, row_carrier, *input.carrier));
        publication.bind(std::move(scope_binding));
        const auto selected = [&](const Vector &value) {
            return publication.selected(value, *input.carrier, scope.rows, row_coordinates);
        };
        std::map<std::string, Value> observations;
        for (const auto &id : responses) observations.emplace(id, selected(input.values(id, scope.rows)));
        Vector weights(scope.rows.size(), Rational{1});
        if (metric_field)
            weights = input.values(*metric_field, scope.rows);
        const auto metric = field::ExactMetric::diagonal(weights);
        stratum.set("metric_diagonal", selected(weights));
        Matrix all(scope.rows.size(), columns.size());
        for (std::size_t j = 0; j < columns.size(); ++j)
            for (std::size_t i = 0; i < scope.rows.size(); ++i) {
                Rational value{1};
                for (const auto &id : columns[j].factors)
                    value *= input.value(id, scope.rows[i]);
                all.set(i, j, value);
            }
        auto models = Value::array();
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
            const auto coefficient_carrier = publication.axis(selected.size()),
                       frame_carrier = publication.axis(model.projection.frame.cols()),
                       kernel_carrier = publication.axis(model.kernel.cols());
            auto model_binding = Value::object();
            model_binding.set("stratum", Value::string(scope.id));
            model_binding.set("model", Value::string(model.id));
            model_binding.set("coefficient_carrier", Value::count(coefficient_carrier));
            model_binding.set("frame_carrier", Value::count(frame_carrier));
            model_binding.set("kernel_carrier", Value::count(kernel_carrier));
            model_binding.set("column_origins", Value(publication.coordinates(selected, coefficient_carrier)));
            model_binding.set("frame_column_origins", Value(publication.coordinates(model.pivots, frame_carrier)));
            Matrix column_lift(columns.size(), selected.size());
            for (std::size_t i = 0; i < selected.size(); ++i) column_lift.set(selected[i], i, 1);
            model_binding.set("column_lift", publication.matrix(column_lift, coefficient_carrier, full_columns));
            publication.bind(std::move(model_binding));
            auto record = Value::object(), axis = Value::array(), identified = Value::array(),
                 response_records = Value::array();
            for (std::size_t j = 0; j < selected.size(); ++j) {
                axis.push(Value::string(columns[selected[j]].id));
                bool unique = true;
                for (std::size_t k = 0; k < model.kernel.cols(); ++k)
                    unique = unique && model.kernel.get(j, k) == 0;
                identified.push(Value::boolean(unique));
            }
            record.set("id", Value::string(model.id));
            record.set("columns", axis);
            record.set("rank", Value::count(model.projection.frame.cols()));
            record.set("residual_dimension",
                       Value::count(scope.rows.size() - model.projection.frame.cols()));
            record.set("coefficient_nullity", Value::count(model.kernel.cols()));
            record.set("coefficient_identified", identified);
            record.set("design", publication.matrix(model.design, coefficient_carrier, row_carrier));
            record.set("independent_frame", publication.matrix(model.projection.frame, frame_carrier, row_carrier));
            record.set("coefficient_kernel", publication.matrix(model.kernel, kernel_carrier, coefficient_carrier));
            for (std::size_t j = 0; j < responses.size(); ++j) {
                auto response = Value::object();
                response.set("particular_coefficients", publication.vector(computed[j].coefficients, coefficient_carrier));
                response.set("represented", publication.selected(computed[j].represented, *input.carrier, scope.rows, row_coordinates));
                response.set("residual", publication.selected(computed[j].residual, *input.carrier, scope.rows, row_coordinates));
                response.set("residual_quadrance", rational(computed[j].quadrance));
                response.set("normal_equations", publication.vector(computed[j].normal, coefficient_carrier));
                response.set("orthogonality_exact", Value::boolean(true));
                response.set("reconstruction_exact", Value::boolean(true));
                response.set("field", Value::string(responses[j]));
                response.set("meaning", input.meaning.at(responses[j]));
                response.set("input", observations.at(responses[j]));
                if (model.id == "full") {
                    auto readings = Value::array();
                    for (const auto &contrast : contrasts) {
                        const auto ambiguity = model.kernel.transpose().apply(contrast.weights);
                        const bool unique = numeric::is_zero(ambiguity);
                        auto reading = Value::object();
                        reading.set("id", Value::string(contrast.id));
                        reading.set("weights", publication.vector(contrast.weights, coefficient_carrier));
                        reading.set("kernel_pairings", publication.vector(ambiguity, kernel_carrier));
                        reading.set("identified", Value::boolean(unique));
                        reading.set("value", unique
                                                 ? rational(numeric::dot(contrast.weights,
                                                                         computed[j].coefficients))
                                                 : Value{});
                        readings.push(reading);
                    }
                    response.set("coefficient_contrasts", readings);
                }
                response_records.push(response);
            }
            record.set("responses", response_records);
            auto execution = Value::object();
            execution.set("owner", Value::string("rfc::numeric::PreparedExactSolve"));
            execution.set("gram_preparations", Value::count(1));
            execution.set("applications", Value::count(responses.size()));
            execution.set("factor_entries", Value::count(model.projection.solve.stored_entries()));
            execution.set("requested_workers", Value::count(report.requested_workers));
            execution.set("worker_slots", Value::count(report.worker_slots));
            execution.set("completed_tiles", Value::count(report.completed_tiles));
            record.set("execution", execution);
            models.push(record);
        }
        stratum.set("models", models);
        auto comparisons = Value::array();
        for (std::size_t m = 0; m < variants.size(); ++m) {
            if (m == 1)
                continue;
            auto comparison = Value::object(), readings = Value::array();
            comparison.set("reduced", Value::string(variants[m].id));
            comparison.set("full", Value::string("full"));
            for (std::size_t j = 0; j < responses.size(); ++j) {
                const auto contribution =
                    numeric::sub(fits[1][j].represented, fits[m][j].represented);
                const Rational gain = fits[m][j].quadrance - fits[1][j].quadrance;
                require(gain >= 0 && gain == metric.inner(contribution, contribution) &&
                            metric.inner(contribution, fits[1][j].residual) == 0,
                        "nested linear model quadrance identity");
                auto reading = Value::object();
                reading.set("field", Value::string(responses[j]));
                reading.set("added_representation", selected(contribution));
                reading.set("added_quadrance", rational(gain));
                reading.set("reduced_residual_quadrance", rational(fits[m][j].quadrance));
                reading.set("fraction_of_reduced_residual",
                            fits[m][j].quadrance == 0
                                ? Value{}
                                : rational(Rational{gain / fits[m][j].quadrance}));
                comparisons_by_axis.emplace(LinearModelAxis{scope.id, variants[m].id, responses[j]},
                    fits[m][j].quadrance == 0 ? std::optional<Rational>{}
                        : std::optional<Rational>{Rational{gain / fits[m][j].quadrance}});
                reading.set("nested_identity_exact", Value::boolean(true));
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
        Value::string("For each declared stratum, minimize (y-X beta)* M (y-X beta) on fixed joint "
                     "observed support. Retain Pi_image(X) y, its orthogonal residual, the "
                     "complete coefficient fibre, and full-versus-reduced field differences. "
                     "Indicator responses describe supplied memberships."));
    return LinearModelFields{publication.finish(std::move(result)), std::move(comparisons_by_axis)};
}
LinearModelFields linear_model_fields(const ScientificStudy &study, const Json &query,
                                      const std::shared_ptr<CpuResources> &resources) {
    return evaluate_linear_model(study, query, false, resources).value();
}
Json linear_model_comparison(const ScientificStudy &study, const Json &query, bool validate,
                             const std::shared_ptr<CpuResources> &resources) {
    auto result = evaluate_linear_model(study, query, validate, resources);
    return result ? render_record_view(result->publication) : Json{};
}
} // namespace gadget::research_detail
