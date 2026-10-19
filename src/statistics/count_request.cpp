#include "../archive_identity.hpp"
#include "count_analysis_detail.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <set>

namespace gadget::statistics::count_detail {
namespace {
using namespace scientific_detail;
void keys(const Json &value, const std::set<std::string> &required,
          const std::set<std::string> &optional = {}) {
    require(json_object_is_type(value.raw(), json_type_object),
            "count request member requires an object");
    for (const auto &name : required)
        require(value.has(name), "count request member absent: " + name);
    json_object_object_foreach(value.raw(), name, item) {
        (void)item;
        require(required.contains(name) || optional.contains(name),
                "unknown count request member: " + std::string(name));
    }
}
double number(const Json &value, bool positive = false) {
    const auto token = value.text();
    double out{};
    if (token.find('/') != std::string::npos)
        out = exact_rational(token).get_d();
    else {
        auto literal = std::string_view(token);
        if (literal.starts_with('+'))
            literal.remove_prefix(1);
        const auto result = std::from_chars(literal.data(), literal.data() + literal.size(), out);
        require(result.ec == std::errc{} && result.ptr == literal.data() + literal.size(),
                "count decimal parameter outside binary64 support");
    }
    require(std::isfinite(out) && (!positive || out > 0),
            "count numerical parameter outside binary64 support");
    return out;
}
std::vector<std::string> labels(const Json &value) {
    source_array(value, "count label array");
    std::vector<std::string> out;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < value.size(); ++i) {
        const auto id = text(value.at(i));
        require(seen.insert(id).second, "duplicate count declaration identity");
        out.push_back(id);
    }
    return out;
}
std::vector<double> numbers(const Json &value, bool positive = false) {
    source_array(value, "count numerical parameter array");
    std::vector<double> out;
    for (std::size_t i = 0; i < value.size(); ++i)
        out.push_back(number(value.at(i), positive));
    return out;
}
std::optional<LogDispersionPrior> prior(const Json &value) {
    if (!value.raw())
        return {};
    keys(value, {"mean", "variance"});
    return LogDispersionPrior{number(value.at("mean")), number(value.at("variance"), true)};
}
std::size_t size(const Json &value, bool positive = true) {
    const auto out =
        count(value, std::numeric_limits<std::size_t>::max(), "count execution admission");
    require(!positive || out > 0, "count execution admission must be positive");
    return out;
}
} // namespace
double exact_count(const rfc::numeric::Rational &value) {
    if (value < 0 || value.get_den() != 1)
        throw std::invalid_argument("count observations must be nonnegative integers");
    const auto converted = value.get_d();
    if (!std::isfinite(converted) || rfc::numeric::Rational(converted) != value)
        throw std::invalid_argument("exact count is not representable in binary64");
    return converted;
}
Configuration admit_count_request(const NativeRecord &request) {
    using namespace scientific_detail;
    const auto native_limits = request.at("limits");
    const auto raw_budget = native_limits.at("request_bytes").integer();
    require(raw_budget > 0 &&
                raw_budget <= rfc::numeric::Integer(std::to_string(source_byte_limit())),
            "count request archive byte admission");
    (void)archive_identity(request.owning_publication(),
                           static_cast<std::size_t>(raw_budget.get_ui()));
    const auto q = render_record_view(request);
    Configuration c;
    keys(q,
         {"schema", "operation", "source_digest", "samples", "features", "counts", "design",
          "model", "normalization", "exposure", "inference", "effect_prior", "correction",
          "outputs", "limits"},
         {"workers"});
    require(q.at("schema").text() == "gadget.statistics.count.request.v1" &&
                q.at("operation").text() == "count_model_analysis",
            "count request schema/operation");
    const auto limits = q.at("limits");
    keys(limits, {"samples", "features", "design_columns", "iterations", "line_search_steps",
                  "numerical_work", "declared_buffer_bytes", "request_bytes"});
    const auto n = size(limits.at("samples")), f = size(limits.at("features")),
               p = size(limits.at("design_columns"));
    require(rfc::machine::checked_size_mul(c.request_bytes =
                                               static_cast<std::size_t>(raw_budget.get_ui()),
                                           32) <= size(limits.at("declared_buffer_bytes")),
            "count request working byte admission");
    c.work = size(limits.at("numerical_work"));
    c.bytes = size(limits.at("declared_buffer_bytes"));
    c.request_bytes = size(limits.at("request_bytes"));
    c.control.iterations = size(limits.at("iterations"));
    c.control.line_search_steps = size(limits.at("line_search_steps"));
    if (q.has("workers"))
        c.workers = size(q.at("workers"));
    const auto samples = q.at("samples"), features = q.at("features");
    source_array(samples, "count samples");
    source_array(features, "count features");
    require(samples.size() > 0 && samples.size() <= n && features.size() > 0 &&
                features.size() <= f,
            "count sample/feature admission");
    const auto minimum = rfc::machine::checked_size_mul(
        rfc::machine::checked_size_add(samples.size(), features.size()), 1024);
    require(minimum <= c.bytes, "count declaration working byte admission");
    std::set<std::string> ids, units, fields;
    std::set<std::size_t> coordinates;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto sample = samples.at(i);
        keys(sample, {"id", "experimental_unit", "coordinate"});
        const auto id = text(sample.at("id")), unit = text(sample.at("experimental_unit"));
        const auto coordinate = size(sample.at("coordinate"), false);
        require(ids.insert(id).second && units.insert(unit).second &&
                    coordinates.insert(coordinate).second,
                "independent count samples require distinct identities, units and coordinates");
        c.samples.push_back({id, unit, {}, coordinate});
    }
    ids.clear();
    for (std::size_t i = 0; i < features.size(); ++i) {
        const auto feature = features.at(i);
        keys(feature, {"id", "field"});
        const auto id = text(feature.at("id")), field = text(feature.at("field"));
        require(ids.insert(id).second && fields.insert(field).second,
                "count feature identity/field repeated");
        c.features.push_back({id, field});
    }
    keys(q.at("counts"), {"quantity", "units"});
    c.quantity = text(q.at("counts").at("quantity"));
    c.units = text(q.at("counts").at("units"));
    const auto design = q.at("design");
    keys(design, {"columns", "rows", "contrast"});
    c.columns = labels(design.at("columns"));
    require(!c.columns.empty() && c.columns.size() <= p, "count design column admission");
    const auto rows = design.at("rows");
    source_array(rows, "count design rows");
    require(rows.size() == c.samples.size(), "count design observation shape differs");
    const auto cells = rfc::machine::checked_size_mul(rows.size(), c.columns.size());
    require(rfc::machine::checked_size_mul(cells, 512) <= c.bytes,
            "count design working byte admission");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        source_array(rows.at(i), "count design row");
        require(rows.at(i).size() == c.columns.size(), "count design column shape differs");
        rfc::numeric::RationalVector values;
        for (std::size_t j = 0; j < c.columns.size(); ++j)
            values.push_back(exact_rational(rows.at(i).at(j).text()));
        c.rows.push_back(std::move(values));
    }
    c.contrast = numbers(design.at("contrast"));
    require(c.contrast.size() == c.columns.size() &&
                std::any_of(c.contrast.begin(), c.contrast.end(), [](double v) { return v != 0; }),
            "count contrast shape or zero contrast");
    const auto model = q.at("model");
    keys(model, {"family", "dispersion"});
    c.family = text(model.at("family"));
    const auto dispersion = model.at("dispersion");
    if (c.family == "poisson") {
        require(!dispersion.raw(), "Poisson dispersion is fixed at zero");
        c.dispersion_method = "fixed";
    } else {
        require(c.family == "negative_binomial", "unsupported count probability family");
        c.dispersion_method = dispersion.at("method").text();
        if (c.dispersion_method == "fixed") {
            keys(dispersion, {"method", "value"});
            c.alpha = number(dispersion.at("value"), true);
        } else {
            require(c.dispersion_method == "profile" || c.dispersion_method == "pooled",
                    "unknown count dispersion method");
            const std::set<std::string> base{"method",        "minimum",     "maximum",
                                             "log_tolerance", "grid_points", "iterations",
                                             "log_prior"};
            keys(dispersion, base,
                 c.dispersion_method == "pooled" ? std::set<std::string>{"reference_fields"}
                                                 : std::set<std::string>{});
            c.profile.minimum = number(dispersion.at("minimum"), true);
            c.profile.maximum = number(dispersion.at("maximum"), true);
            c.profile.log_tolerance = number(dispersion.at("log_tolerance"), true);
            c.profile.grid_points = size(dispersion.at("grid_points"));
            c.profile.iterations = size(dispersion.at("iterations"));
            c.profile.prior = prior(dispersion.at("log_prior"));
            require(c.profile.maximum > c.profile.minimum && c.profile.grid_points >= 3 &&
                        rfc::machine::checked_size_mul(c.profile.grid_points, 32) <= c.bytes,
                    "dispersion domain or grid working admission");
            if (c.dispersion_method == "pooled") {
                require(dispersion.has("reference_fields"),
                        "pooled dispersion requires explicit reference fields");
                c.pool_fields = labels(dispersion.at("reference_fields"));
                require(!c.pool_fields.empty() && c.pool_fields.size() <= f,
                        "dispersion reference admission");
            }
        }
    }
    const auto normalization = q.at("normalization");
    c.normalization = normalization.at("method").text();
    if (c.normalization == "none")
        keys(normalization, {"method"});
    else if (c.normalization == "declared_factors") {
        keys(normalization, {"method", "factors", "center"});
        c.factors = numbers(normalization.at("factors"), true);
        c.center = normalization.at("center").boolean();
        require(c.factors.size() == c.samples.size(), "normalization factor shape differs");
    } else {
        require(c.normalization == "median_ratio", "unknown count normalization method");
        keys(normalization, {"method", "reference_fields"});
        c.normalization_fields = labels(normalization.at("reference_fields"));
        require(!c.normalization_fields.empty() && c.normalization_fields.size() <= f,
                "normalization reference admission");
        c.center = true;
    }
    if (q.at("exposure").raw()) {
        keys(q.at("exposure"), {"values", "meaning"});
        (void)text(q.at("exposure").at("meaning"));
        c.exposure = numbers(q.at("exposure").at("values"), true);
        require(c.exposure.size() == c.samples.size(), "exposure shape differs");
    }
    const auto inference = q.at("inference");
    keys(inference, {"method", "confidence", "reduced_columns"});
    c.inference = text(inference.at("method"));
    c.confidence = number(inference.at("confidence"), true);
    (void)normal_critical(c.confidence);
    require(c.inference == "wald" || c.inference == "likelihood_ratio",
            "unknown count inference method");
    if (c.inference == "wald")
        require(!inference.at("reduced_columns").raw(), "Wald inference has no reduced design");
    else {
        const auto reduced = inference.at("reduced_columns");
        source_array(reduced, "reduced design columns");
        std::set<std::size_t> selected;
        for (std::size_t i = 0; i < reduced.size(); ++i) {
            const auto j = size(reduced.at(i), false);
            require(j < c.columns.size() && selected.insert(j).second,
                    "reduced design column absent/repeated");
            c.reduced.push_back(j);
        }
        require(c.reduced.size() < c.columns.size() && !c.profile.prior,
                "likelihood ratio requires a nested unpenalized model");
    }
    c.effect_prior = prior(q.at("effect_prior"));
    c.correction = text(q.at("correction"));
    (void)adjust_log_probabilities({}, c.correction);
    keys(q.at("outputs"), {"fitted_means"});
    c.fitted = q.at("outputs").at("fitted_means").boolean();
    return c;
}
} // namespace gadget::statistics::count_detail
