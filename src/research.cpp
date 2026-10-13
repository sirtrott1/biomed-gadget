#include "research_calibration.hpp"
#include "research_factored.hpp"
#include "research_fields.hpp"
#include "research_linear.hpp"
#include "research_comparisons.hpp"
#include "scientific_result_fields.hpp"
#include <rfc/numeric/prepared_exact_frames.hpp>
#include "research_prepared.hpp"
#include "scientific_detail.hpp"
#include "scientific_projection.hpp"
#include "scientific_relational_field.hpp"
#include <algorithm>
#include <array>
#include <gadget/research.hpp>
#include <gadget/resources.hpp>
#include <limits>
#include <map>
#include <rfc/algebra/materialize.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/field/augmented_green.hpp>
#include <rfc/geometry/chain_map.hpp>
#include <rfc/numeric/checked_exact_actions.hpp>
#include <rfc/runtime/cpu_tiles.hpp>
#include <rfc/temporal/sector_transport.hpp>
#include <rfc/temporal/transfer.hpp>
#include <rfc/topology/integral_homology.hpp>
#include <set>

namespace gadget {
namespace {
using namespace rfc;
using namespace scientific_detail;
using Matrix = numeric::ExactMatrix;
using Indices = std::vector<std::size_t>;
Json q(Rational v) {
    v.canonicalize();

    return Json::string(v.get_str());
}
Json vec(const Vector &v) {
    auto a = Json::array();
    for (const auto &x : v)
        a.push(q(x));
    return a;
}
Json matrix(const Matrix &m) {
    auto a = Json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        auto r = Json::array();
        for (std::size_t j = 0; j < m.cols(); ++j)
            r.push(q(m.get(i, j)));
        a.push(r);
    }
    auto out = Json::object();
    out.set("rows", cardinality(m.rows()));
    out.set("columns", cardinality(m.cols()));
    out.set("values", a);
    return out;
}
Json indices(const auto &v) {
    auto a = Json::array();
    for (auto x : v)
        a.push(cardinality(x));
    return a;
}
using research_detail::Inputs;
using research_detail::names;
numeric::ExactActionLimits comparison_limits(const runtime::ExecutionLease &lease, std::size_t dimension) {
    const auto &limits = lease.policy().limits;
    return {dimension, limits.sparse_entries, limits.scalar_work, limits.integer_bits};
}
void admit_comparison(const Inputs &input, const runtime::ExecutionLease &lease) {
    std::size_t slots = 0, bits = 0, entries = 0;
    for (const auto &[id, tensor] : input.fields) {
        (void)id;
        slots = machine::checked_size_add(slots, tensor->tensor.support().count());
        entries = machine::checked_size_add(entries, tensor->tensor.stored_nonzeros());
        for (const auto &entry : tensor->tensor.values().entries())
            bits = machine::checked_size_add(bits, storage::native_execution_detail::integer_bits(entry.value));
    }
    runtime::admit_execution({input.fields.size(), slots, entries, 1, bits}, lease.policy());
}
ResultValue comparison_matrix(ResultFields &fields, const Matrix &value,
                              std::uint64_t input, std::uint64_t output) {
    auto result = fields.matrix(value, input, output);
    result.set("display", ResultValue::string("dense_matrix"));
    return result;
}
void bind_comparison_source(const ScientificStudy &study, const Inputs &input, ResultFields &fields) {
    auto binding = ResultValue::object(), pins = ResultValue::array();
    binding.set("kind", ResultValue::string("declared_comparison_source"));
    binding.set("original_study_publication", fields.publication(study.native_publication()));
    binding.set("grade", ResultValue::count(input.grade.value));
    binding.set("source_carrier", ResultValue::count(*input.carrier));
    for (const auto &[id, tensor] : input.fields) {
        const auto actual = study.field_accession(id);
        require(actual.tensor() && actual.carrier() && actual.carrier()->key == *tensor->output &&
                storage::temporal_endpoint_detail::tensor_equal(*actual.tensor(), tensor->tensor),
                "comparison field differs from its original accession");
        auto pin = ResultValue::object();
        pin.set("field", ResultValue::string(id));
        pin.set("source_namespace", ResultValue::string(std::to_string(actual.source_namespace())));
        pin.set("source_cutoff", ResultValue::count(actual.sequence()));
        pin.set("source_tensor", ResultValue::count(tensor->key));
        pin.set("source_carrier", ResultValue::count(*tensor->output));
        pin.set("source_field", fields.stored(tensor->key, "vector"));
        pins.push(std::move(pin));
    }
    binding.set("fields", std::move(pins));
    fields.bind(std::move(binding));
}
std::uint64_t bind_comparison_rows(const Inputs &input, const Indices &rows, ResultFields &fields) {
    const auto axis = fields.axis(rows.size());
    auto binding = ResultValue::object();
    binding.set("kind", ResultValue::string("original_comparison_rows"));
    binding.set("source_carrier", ResultValue::count(*input.carrier));
    binding.set("row_carrier", ResultValue::count(axis));
    binding.set("row_coordinates", ResultValue(fields.coordinates(rows, axis)));
    Matrix lift(input.entities.size(), rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) lift.set(rows[i], i, 1);
    binding.set("source_lift", comparison_matrix(fields, lift, axis, *input.carrier));
    fields.bind(std::move(binding));
    return axis;
}
// Inputs to one declared scientific projection request. All maps, Gram
// construction, factorization and solve actions belong to canonical RexGraph.
struct FrameRequest : PreparedProjection {
    FrameRequest(Matrix independent_frame, numeric::ExactActionLimits limits,
                 const runtime::ExecutionLease &lease)
        : PreparedProjection(std::move(independent_frame), {}, limits,
            [lease] {lease.check_cancelled();}) {}
    ResultValue receipt(std::size_t applications, const runtime::CpuTileRunReport &report) const {
        auto out = ResultValue::object();
        out.set("owner", ResultValue::string("rfc::numeric::PreparedExactSolve"));
        out.set("gram_preparations", ResultValue::count(1));
        out.set("applications", ResultValue::count(applications));
        out.set("carrier_dimension", ResultValue::count(frame.rows()));
        out.set("frame_rank", ResultValue::count(frame.cols()));
        out.set("factor_entries", ResultValue::count(solve.stored_entries()));
        out.set("requested_workers", ResultValue::count(report.requested_workers));
        out.set("worker_slots", ResultValue::count(report.worker_slots));
        out.set("completed_tiles", ResultValue::count(report.completed_tiles));
        out.set("arithmetic", ResultValue::string("exact rational"));
        return out;
    }
};
// A contrast family is a field on explicitly identified matched units.
// The native constant-span projector retains its complete consensus and residual.
std::optional<NativeRecord> paired(const ScientificStudy &study, const Inputs &in, const Json &query, bool validate,
            const std::shared_ptr<CpuResources> &resources,
            const runtime::ExecutionLease &lease) {
    members(
        query,
        {"operation", "pairs", "selection", "support", "normalization", "output_units", "workers"},
        {"benchmark"});
    require(query.at("support").text() == "joint_reported", "paired contrast support");
    const auto pairs = query.at("pairs");
    require(json_object_is_type(pairs.raw(), json_type_array) && pairs.size(),
            "paired contrast requires matched units");
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        count(query.at("workers"), cpu->coordinator().capacity(), "paired workers");
    require(workers > 0, "paired workers must be positive");
    const auto normalization = query.at("normalization");
    members(normalization, {"method", "evidence"});
    (void)text(normalization.at("method"));
    (void)text(query.at("output_units"));
    auto selection = names(query.at("selection"), in.fields.size(), true);
    std::vector<std::string> complete;
    std::set<std::string> used, identities;
    std::vector<Rational> left_scales, right_scales;
    auto scale = [&lease](const Json &v) {
        const auto value = exact_rational(text(v, source_byte_limit()));
        require(value > 0, "scale must be positive");
        runtime::admit_execution({1, 1, 0, 0, storage::native_execution_detail::integer_bits(value)}, lease.policy());
        return value;
    };
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto pair = pairs.at(i);
        members(pair, {"id", "left", "right", "left_scale", "right_scale"});
        require(identities.insert(text(pair.at("id"))).second, "duplicate matched unit");
        const auto left = text(pair.at("left")), right = text(pair.at("right"));
        require(used.insert(left).second && used.insert(right).second,
                "a sample field cannot occur in two matched units");
        complete.push_back(left);
        complete.push_back(right);
        left_scales.push_back(scale(pair.at("left_scale")));
        right_scales.push_back(scale(pair.at("right_scale")));
    }
    const auto first = in.meaning.at(complete.front());
    require(first.at("channel").text() == "E" && first.at("quantity").text() == "read_count",
            "paired expression requires counts");
    for (const auto &id : complete) {
        const auto meaning = in.meaning.at(id);
        for (const auto *member : {"channel", "quantity", "units", "value_type", "method"})
            require(meaning.at(member).text() == first.at(member).text(),
                    "paired sample measurement meanings differ");
        for (std::size_t row = 0; row < in.entities.size(); ++row)
            if (in.fields.at(id)->tensor.support().contains(row)) {
                const auto value = in.value(id, row);
                require(value >= 0 && value.get_den() == 1,
                        "read counts must be nonnegative integers");
            }
    }
    if (query.has("benchmark")) {
        const auto b = query.at("benchmark");
        members(b, {"effect", "adjusted_p", "alpha"});
        const auto alpha = exact_decimal(b.at("alpha").text());
        require(alpha > 0 && alpha < 1, "benchmark alpha");
        for (const auto *role : {"effect", "adjusted_p"}) {
            const auto id = text(b.at(role));
            const auto meaning = in.meaning.at(id);
            require(meaning.at("channel").text() == "observation" &&
                        meaning.at("units").text() == "dimensionless",
                    "benchmark fields require declared dimensionless upstream observations");
            require(
                meaning.at("quantity").text() ==
                    (std::string(role) == "effect" ? "model_log2_fold_change" : "adjusted_p_value"),
                "benchmark field quantity");
        }
        const auto id = b.at("adjusted_p").text();
        for (std::size_t r = 0; r < in.entities.size(); ++r)
            if (in.fields.at(id)->tensor.support().contains(r))
                require(in.value(id, r) >= 0 && in.value(id, r) <= 1,
                        "adjusted p token outside [0,1]");
    }
    const auto rows = in.select(complete, selection);
    const auto selected_rows = in.select({}, selection);
    lease.check_cancelled();
    admit_comparison(in, lease);
    if (validate)
        return std::nullopt;
    ResultFields publication(study.native_publication());
    bind_comparison_source(study, in, publication);
    const auto n = pairs.size();
    Matrix consensus(n, 1);
    for (std::size_t i = 0; i < n; ++i)
        consensus.set(i, 0, Rational{1});
    auto reservation = cpu->reserve(workers);
    const auto limits = comparison_limits(lease, n);
    const numeric::PreparedExactFrames frames(std::move(consensus), limits, [lease] {lease.check_cancelled();});
    const FrameRequest projection(frames.image(), limits, lease);
    const auto pair_axis = publication.axis(n), coefficient_axis = publication.axis(1),
               frame_axis = publication.axis(frames.rank()), kernel_axis = publication.axis(frames.kernel().cols());
    auto axes = ResultValue::object();
    axes.set("kind", ResultValue::string("ordered_matched_units"));
    axes.set("pairs", pairs);
    axes.set("pair_carrier", ResultValue::count(pair_axis));
    axes.set("coefficient_carrier", ResultValue::count(coefficient_axis));
    axes.set("frame_carrier", ResultValue::count(frame_axis));
    axes.set("kernel_carrier", ResultValue::count(kernel_axis));
    axes.set("design", comparison_matrix(publication, frames.original(), coefficient_axis, pair_axis));
    axes.set("frame", comparison_matrix(publication, frames.image(), frame_axis, pair_axis));
    axes.set("coordinates", comparison_matrix(publication, frames.coordinates(), coefficient_axis, frame_axis));
    axes.set("kernel", comparison_matrix(publication, frames.kernel(), kernel_axis, coefficient_axis));
    axes.set("metric_diagonal", publication.vector(Vector(n, Rational{1}), pair_axis));
    axes.set("left_scales", publication.vector(left_scales, pair_axis));
    axes.set("right_scales", publication.vector(right_scales, pair_axis));
    publication.bind(std::move(axes));
    bind_comparison_rows(in, rows, publication);
    const Vector ones(n, Rational{1});
    const Rational denominator{std::to_string(n)};
    std::vector<ResultValue> records(rows.size());
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), rows.size()});
    const auto report = executor.run_tiles(rows.size(), [&](std::size_t tile, std::size_t,
                                                            const runtime::CpuTileStop &stop) {
        stop.check();
        lease.check_cancelled();
        numeric::ExactActionGuard guard(limits, [&] {stop.check(); lease.check_cancelled();});
        const auto row = rows[tile];
        Vector left, right;
        std::size_t positive = 0, negative = 0, zero = 0;
        for (std::size_t i = 0; i < n; ++i) {
            left.push_back(guard.mul(in.value(complete[2 * i], row), left_scales[i]));
            right.push_back(guard.mul(in.value(complete[2 * i + 1], row), right_scales[i]));
        }
        const auto delta = numeric::sub(left, right);
        const auto consensus_coefficients = projection.coefficients(delta, [&] {stop.check(); lease.check_cancelled();});
        const auto glued = projection.frame.apply(consensus_coefficients);
        const auto residual = numeric::sub(delta, glued);
        guard.vector(delta); guard.vector(glued); guard.vector(residual);
        require(numeric::add(glued, residual) == delta && numeric::dot(glued, residual) == 0,
                "paired consensus reconstruction");
        for (const auto &v : delta) {
            if (v > 0)
                ++positive;
            else if (v < 0)
                ++negative;
            else
                ++zero;
        }
        const auto scalar = [&](const Rational &v) {guard.value(v); return publication.stage_scalar(v);};
        auto record = ResultValue::object();
        record.set("entity", ResultValue::string(in.entities[row]));
        record.set("coordinate", ResultValue::count(row));
        record.set("left", ResultFields::stage_vector(left, pair_axis));
        record.set("right", ResultFields::stage_vector(right, pair_axis));
        record.set("delta", ResultFields::stage_vector(delta, pair_axis));
        record.set("consensus", ResultFields::stage_vector(glued, pair_axis));
        record.set("consensus_coefficients", ResultFields::stage_vector(consensus_coefficients, frame_axis));
        record.set("residual", ResultFields::stage_vector(residual, pair_axis));
        record.set("mean_difference", scalar(glued.front()));
        const Rational lm = numeric::dot(left, ones) / denominator;
        const Rational rm = numeric::dot(right, ones) / denominator;
        record.set("left_mean", scalar(lm));
        record.set("right_mean", scalar(rm));
        record.set("mean_ratio", rm != 0 ? scalar(Rational(lm / rm)) : ResultValue{});
        record.set("consensus_quadrance", scalar(numeric::dot(glued, glued)));
        record.set("heterogeneity_quadrance", scalar(numeric::dot(residual, residual)));
        record.set("positive_pairs", ResultValue::count(positive));
        record.set("negative_pairs", ResultValue::count(negative));
        record.set("zero_pairs", ResultValue::count(zero));
        if (query.has("benchmark")) {
            const auto b = query.at("benchmark");
            const auto effect = b.at("effect").text(), adj = b.at("adjusted_p").text();
            auto comparison = ResultValue::object();
            const bool has_effect = in.fields.at(effect)->tensor.support().contains(row);
            const bool has_p = in.fields.at(adj)->tensor.support().contains(row);
            comparison.set("effect", has_effect ? scalar(in.value(effect, row)) : ResultValue{});
            comparison.set("adjusted_p", has_p ? scalar(in.value(adj, row)) : ResultValue{});
            comparison.set("significant", has_p ? ResultValue::boolean(in.value(adj, row) <
                                                                exact_decimal(b.at("alpha").text()))
                                                : ResultValue{});
            if (has_effect) {
                const auto score = in.value(effect, row), mean = glued.front();
                comparison.set("direction_agrees", ResultValue::boolean(((score > 0) - (score < 0)) ==
                                                                 ((mean > 0) - (mean < 0))));
            } else
                comparison.set("direction_agrees", Json{});
            record.set("benchmark", comparison);
        }
        records[tile] = std::move(record);
    });
    auto out = ResultValue::object(), values = ResultValue::array();
    for (const auto &record : records)
        values.push(record);
    out.set("features", values);
    out.set("pairs", pairs);
    out.set("reported", ResultValue::count(rows.size()));
    out.set("missing", ResultValue::count(selected_rows.size() - rows.size()));
    out.set("selected_rows", ResultValue::count(selected_rows.size()));
    out.set("excluded_by_selection", ResultValue::count(in.entities.size() - selected_rows.size()));
    out.set("units", query.at("output_units"));
    out.set("normalization", normalization);
    out.set("metric", ResultValue::string("identity on the declared matched-unit carrier"));
    out.set("workers", ResultValue::count(report.worker_slots));
    out.set("tiles", ResultValue::count(report.completed_tiles));
    out.set("reconstruction_exact", ResultValue::boolean(true));
    out.set("execution", projection.receipt(rows.size(), report));
    lease.check_cancelled();
    return publication.finish(std::move(out));
}
std::optional<NativeRecord> controlled(const ScientificStudy &study, const Inputs &in, const Json &query, bool validate,
                const std::shared_ptr<CpuResources> &resources,
            const runtime::ExecutionLease &lease) {
    members(query,
            {"operation", "fields", "endpoint", "controls", "selection", "support", "metric",
             "intercept"},
            {"workers", "gram_fields"});
    require(query.at("support").text() == "joint_reported" &&
                query.at("metric").text() == "identity" && query.at("intercept").boolean(),
            "controlled profile requires joint reported support, identity metric and intercept");
    const auto fs = names(query.at("fields"), in.fields.size()),
               controls = names(query.at("controls"), in.fields.size(), true),
               selection = names(query.at("selection"), in.fields.size(), true);
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        query.has("workers")
            ? count(query.at("workers"), cpu->coordinator().capacity(), "control workers")
            : 1;
    require(workers > 0, "control workers must be positive");
    const auto gram_fields =
        query.has("gram_fields") ? names(query.at("gram_fields"), fs.size(), true) : fs;
    std::map<std::string, std::size_t> gram_columns;
    for (std::size_t i = 0; i < gram_fields.size(); ++i) {
        require(std::find(fs.begin(), fs.end(), gram_fields[i]) != fs.end(),
                "Gram selection must name a tested field");
        gram_columns.emplace(gram_fields[i], i);
    }
    const auto endpoint = text(query.at("endpoint"));
    require(std::find(controls.begin(), controls.end(), endpoint) == controls.end(),
            "endpoint cannot be its own control");
    require(std::find(fs.begin(), fs.end(), endpoint) == fs.end(),
            "endpoint cannot be its own tested field");
    auto complete = fs;
    complete.insert(complete.end(), controls.begin(), controls.end());
    complete.push_back(endpoint);
    const auto rows = in.select(complete, selection);
    lease.check_cancelled();
    admit_comparison(in, lease);
    if (validate)
        return std::nullopt;
    ResultFields publication(study.native_publication());
    bind_comparison_source(study, in, publication);
    auto reservation = cpu->reserve(workers);
    Matrix declared(rows.size(), controls.size() + 1);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        declared.set(i, 0, Rational{1});
        for (std::size_t j = 0; j < controls.size(); ++j)
            declared.set(i, j + 1, in.value(controls[j], rows[i]));
    }
    const auto limits = comparison_limits(lease, std::max(rows.size(), declared.cols()));
    const numeric::PreparedExactFrames frames(std::move(declared), limits, [lease] {lease.check_cancelled();});
    const FrameRequest projection(frames.image(), limits, lease);
    const auto row_axis = bind_comparison_rows(in, rows, publication),
               coefficient_axis = publication.axis(frames.original().cols()),
               frame_axis = publication.axis(frames.rank()),
               kernel_axis = publication.axis(frames.kernel().cols()),
               gram_axis = publication.axis(gram_fields.size());
    const auto selected = [&](const Vector &v) {return ResultFields::stage_selected(v, *in.carrier, in.entities.size(), rows);};
    auto axes = ResultValue::object(), columns = ResultValue::array();
    columns.push(ResultValue::string("intercept"));
    for (const auto &id : controls) columns.push(ResultValue::string(id));
    axes.set("kind", ResultValue::string("declared_controls"));
    axes.set("columns", columns);
    axes.set("coefficient_carrier", ResultValue::count(coefficient_axis));
    axes.set("frame_carrier", ResultValue::count(frame_axis));
    axes.set("kernel_carrier", ResultValue::count(kernel_axis));
    axes.set("gram_carrier", ResultValue::count(gram_axis));
    axes.set("metric_diagonal", publication.vector(Vector(rows.size(), Rational{1}), row_axis));
    axes.set("design", comparison_matrix(publication, frames.original(), coefficient_axis, row_axis));
    axes.set("frame", comparison_matrix(publication, frames.image(), frame_axis, row_axis));
    axes.set("coordinates", comparison_matrix(publication, frames.coordinates(), coefficient_axis, frame_axis));
    axes.set("kernel", comparison_matrix(publication, frames.kernel(), kernel_axis, coefficient_axis));
    publication.bind(std::move(axes));
    const auto coefficients = [&](const Vector &c) {
        Vector out(frames.original().cols());
        for (std::size_t i = 0; i < frames.pivots().size(); ++i) out[frames.pivots()[i]] = c[i];
        require(frames.original().apply(out) == projection.frame.apply(c), "control coefficient reconstruction");
        return out;
    };
    const auto &frame = projection.frame;
    const auto y = in.values(endpoint, rows), cy = projection.coefficients(y, [lease] {lease.check_cancelled();}),
               py = projection.frame.apply(cy), ry = numeric::sub(y, py);
    require(numeric::is_zero(projection.dual.apply(ry)) && numeric::dot(py, ry) == 0 &&
                numeric::add(py, ry) == y,
            "endpoint projection reconstruction");
    auto result = ResultValue::object();
    result.set("entities", in.identity(rows));
    result.set("coordinates", publication.indices(rows, row_axis));
    result.set("excluded_rows", ResultValue::count(in.entities.size() - rows.size()));
    result.set("n", ResultValue::count(rows.size()));
    result.set("control_columns", ResultValue::count(frames.original().cols()));
    result.set("control_rank", ResultValue::count(frame.cols()));
    result.set("control_frame", comparison_matrix(publication, frame, frame_axis, row_axis));
    result.set("residual_dimension", ResultValue::count(rows.size() - frame.cols()));
    result.set("endpoint_input", selected(y));
    result.set("endpoint_represented", selected(py));
    result.set("endpoint_residual", selected(ry));
    result.set("control_design", comparison_matrix(publication, frames.original(), coefficient_axis, row_axis));
    result.set("control_coordinates", comparison_matrix(publication, frames.coordinates(), coefficient_axis, frame_axis));
    result.set("control_kernel", comparison_matrix(publication, frames.kernel(), kernel_axis, coefficient_axis));
    result.set("endpoint_control_coefficients", publication.vector(coefficients(cy), coefficient_axis));
    const auto ey = numeric::dot(ry, ry);
    result.set("endpoint_residual_quadrance", publication.stage_scalar(ey));
    result.set("endpoint_identified", ResultValue::boolean(ey != 0));
    std::vector<Vector> gram_vectors(gram_fields.size());
    std::vector<ResultValue> records(fs.size());
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), fs.size()});
    const auto report = executor.run_tiles(
        fs.size(), [&](std::size_t f, std::size_t, const runtime::CpuTileStop &stop) {
            stop.check();
            lease.check_cancelled();
            auto x = in.values(fs[f], rows), cx = projection.coefficients(x, [&] {stop.check(); lease.check_cancelled();}),
                 p = projection.frame.apply(cx), r = numeric::sub(x, p);
            numeric::ExactActionGuard guard(limits, [&] {stop.check(); lease.check_cancelled();});
            guard.vector(x); guard.vector(p); guard.vector(r);
            require(numeric::is_zero(projection.dual.apply(r)) && numeric::dot(p, r) == 0 &&
                        numeric::add(p, r) == x,
                    "field projection reconstruction");
            if (gram_columns.contains(fs[f]))
                gram_vectors[gram_columns.at(fs[f])] = r;
            const auto pairing = numeric::dot(r, ry);
            const auto scalar = [&](const Rational &v) {guard.value(v); return publication.stage_scalar(v);};
            auto record = ResultValue::object();
            record.set("field", ResultValue::string(fs[f]));
            record.set("units", in.meaning.at(fs[f]).at("units"));
            record.set("input", selected(x));
            record.set("control_coefficients", ResultFields::stage_vector(coefficients(cx), coefficient_axis));
            record.set("represented", selected(p));
            record.set("residual", selected(r));
            record.set("pairing", scalar(pairing));
            record.set("residual_quadrance", scalar(numeric::dot(r, r)));
            record.set("endpoint_coefficient", ey == 0 ? ResultValue{} : scalar(Rational{pairing / ey}));
            record.set("orthogonality_exact", ResultValue::boolean(true));
            records[f] = std::move(record);
        });
    auto output = ResultValue::array(), gram_names = ResultValue::array();
    for (const auto &record : records)
        output.push(record);
    Matrix residuals(rows.size(), gram_fields.size());
    for (std::size_t j = 0; j < gram_fields.size(); ++j) {
        gram_names.push(ResultValue::string(gram_fields[j]));
        for (std::size_t i = 0; i < rows.size(); ++i)
            residuals.set(i, j, gram_vectors[j][i]);
    }
    result.set("fields", output);
    result.set("residual_gram",
               comparison_matrix(publication, field::metric_gram(residuals, projection.metric, residuals), gram_axis, gram_axis));
    result.set("gram_fields", gram_names);
    result.set("execution", projection.receipt(machine::checked_size_add(fs.size(), 1), report));
    result.set(
        "definition",
        ResultValue::string(
            "r_x=(I-P_C)x; r_y=(I-P_C)y; pairing=<r_x,r_y>; endpoint coefficient=pairing/<r_y,r_y> "
            "when identified. C contains the intercept and declared controls. Descriptive "
            "conditional linear contrast on the observed carrier."));
    lease.check_cancelled();
    return publication.finish(std::move(result));
}
using Pattern = std::vector<int>;
struct Vocabulary {
    std::vector<Pattern> patterns;
    std::map<Pattern, std::size_t> lookup;
    std::vector<std::array<std::size_t, 4>> faces;
    core::RelationalComplex open;
    core::RelationalComplex complex;
    explicit Vocabulary(const std::set<Pattern> &observed)
        : open(core::RelationalComplexBuilder(core::RelationSignature{}).freeze()), complex(open) {
        patterns.assign(observed.begin(), observed.end());
        const auto width = patterns.at(0).size();
        require(patterns.size() <= 128, "trajectory state budget (128 observed states)");
        core::RelationSignature sig;
        const auto state =
            sig.add_type("declared_joint_change_state", core::BoundaryPolicy::DeclaredExact,
                         static_cast<std::uint16_t>(width), static_cast<std::uint16_t>(width));
        const auto face = sig.add_type("declared_coordinate_exchange_identity",
                                       core::BoundaryPolicy::DeclaredExact, 4, 4);
        core::RelationalComplexBuilder builder(sig);
        for (std::size_t f = 0; f < width; ++f)
            for (std::size_t s = 0; s < 3; ++s)
                builder.declare_witness_identity(core::IdentityKey{1 + f * 3 + s});
        std::uint64_t occurrence = 0;
        for (std::size_t i = 0; i < patterns.size(); ++i) {
            lookup.emplace(patterns[i], i);
            std::vector<core::BoundaryOccurrenceInput> terms;
            for (std::size_t f = 0; f < width; ++f)
                terms.push_back(
                    {core::IdentityKey{1 + f * 3 + static_cast<std::size_t>(patterns[i][f])},
                     true,
                     false,
                     Rational{1},
                     {},
                     core::OccurrenceKey{++occurrence}});
            builder.add_grade_one(core::RelationKey{i + 1}, state, terms);
        }
        open = builder.freeze();
        builder.declare_grade(core::Grade{2});
        for (std::size_t a = 0; a < width; ++a)
            for (std::size_t b = a + 1; b < width; ++b)
                for (std::size_t i = 0; i < patterns.size(); ++i)
                    for (std::size_t j = 0; j < patterns.size(); ++j) {
                        const auto &x = patterns[i];
                        const auto &y = patterns[j];
                        bool fixed = true;
                        for (std::size_t f = 0; f < width; ++f)
                            if (f != a && f != b && x[f] != y[f])
                                fixed = false;
                        if (!fixed || x[a] >= y[a] || x[b] >= y[b])
                            continue;
                        auto u = x, v = x;
                        u[a] = y[a];
                        v[b] = y[b];
                        if (!lookup.contains(u) || !lookup.contains(v))
                            continue;
                        faces.push_back({i, lookup.at(u), lookup.at(v), j});
                        require(faces.size() <= 4096, "trajectory face budget");
                        std::vector<core::HigherIncidenceInput> terms;
                        for (std::size_t k = 0; k < 4; ++k)
                            terms.push_back({static_cast<std::uint32_t>(faces.back()[k]),
                                             true,
                                             k == 1 || k == 2,
                                             Rational{1},
                                             {},
                                             core::OccurrenceKey{++occurrence}});
                        builder.add_higher(core::Grade{2}, core::RelationKey{faces.size()}, face,
                                           terms);
                    }
        complex = builder.freeze();
    }
    algebra::SubcomplexMask mask(const std::vector<std::size_t> &counts,
                                 std::size_t minimum = 1) const {
        algebra::SubcomplexMask out(complex);
        for (std::size_t i = 0; i < complex.grade(core::Grade{0}).cell_count(); ++i)
            out.select({core::Grade{0}, static_cast<std::uint32_t>(i)});
        for (std::size_t i = 0; i < counts.size(); ++i)
            if (counts[i] >= minimum)
                out.select({core::Grade{1}, static_cast<std::uint32_t>(i)});
        for (std::size_t i = 0; i < faces.size(); ++i) {
            bool all = true;
            for (auto s : faces[i])
                all = all && counts[s] >= minimum;
            if (all)
                out.select({core::Grade{2}, static_cast<std::uint32_t>(i)});
        }
        require(out.boundary_closed(complex), "trajectory closure");
        return out;
    }
    algebra::MaterializedSubcomplex realize(const algebra::SubcomplexMask &mask_) const {
        if (mask_.grade(core::Grade{2}).count())
            return algebra::materialize_subcomplex(complex, mask_);
        algebra::SubcomplexMask lower(open);
        for (std::size_t k = 0; k < 2; ++k)
            mask_.grade(core::Grade{static_cast<std::uint16_t>(k)})
                .for_each_set([&](std::size_t i) {
                    lower.select({core::Grade{static_cast<std::uint16_t>(k)},
                                  static_cast<std::uint32_t>(i)});
                });
        auto out = algebra::materialize_subcomplex(open, lower);
        out.to_ambient.resize(3);
        out.from_ambient.resize(3);
        out.from_ambient[2].assign(faces.size(), algebra::RelationalMorphism::absent);
        return out;
    }
};
std::vector<field::ExactMetric> metrics(const core::RelationalComplex &c) {
    std::vector<field::ExactMetric> out;
    for (std::size_t k = 0; k < c.grade_count(); ++k)
        out.push_back(field::ExactMetric::identity(
            c.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count()));
    return out;
}
Matrix inclusion(const algebra::MaterializedSubcomplex &a, const algebra::MaterializedSubcomplex &b,
                 std::size_t grade = 1) {
    std::vector<geometry::LinearGradeMap> maps;
    for (std::size_t k = 0; k < a.complex.grade_count(); ++k) {
        geometry::LinearGradeMap m(a.to_ambient[k].size());
        for (std::size_t i = 0; i < m.size(); ++i) {
            const auto target = b.from_ambient[k].at(a.to_ambient[k][i]);
            require(target != algebra::RelationalMorphism::absent, "invalid chain inclusion");
            m[i].push_back({target, Rational{1}});
        }
        maps.push_back(std::move(m));
    }
    require(geometry::RelationalChainMap(a.complex, b.complex, std::move(maps)).is_chain_map(),
            "inclusion chain certificate");
    Matrix out(b.to_ambient[grade].size(), a.to_ambient[grade].size());
    for (std::size_t i = 0; i < a.to_ambient[grade].size(); ++i) {
        const auto target = b.from_ambient[grade].at(a.to_ambient[grade][i]);
        require(target != algebra::RelationalMorphism::absent, "invalid inclusion");
        out.set(target, i, Rational{1});
    }
    return out;
}
Json topology_record(const core::RelationalComplex &c) {
    auto out = Json::object();
    const auto h = topology::integral_homology(c, core::Grade{1});
    out.set("H1_rank", cardinality(h.free_rank));
    out.set("cycle_rank", cardinality(h.cycle_rank));
    out.set("boundary_rank", cardinality(h.boundary_rank));
    auto torsion = Json::array();
    for (const auto &v : h.torsion_invariants)
        torsion.push(q(Rational{v}));
    out.set("torsion", torsion);
    out.set("B1", matrix(field::materialize_boundary_matrix(c, core::Grade{1})));
    out.set("B2", matrix(c.grade_count() > 2 ? field::materialize_boundary_matrix(c, core::Grade{2})
                                             : Matrix(c.grade(core::Grade{1}).cell_count(), 0)));
    return out;
}
Matrix image_on_harmonic(const Matrix &cycles, const Matrix &map, const field::HodgeProjectors &p) {
    Matrix out(map.rows(), cycles.cols());
    for (std::size_t j = 0; j < cycles.cols(); ++j) {
        auto x = p.project_harmonic(map.apply(field::matrix_column(cycles, j)));
        for (std::size_t i = 0; i < x.size(); ++i)
            out.set(i, j, x[i]);
    }
    return out;
}
Json sectors_from_fields(const field::ExactHodgeOperator &op, const std::array<Vector, 4> &fields,
                         const Rational &lambda, const Vector *prepared_green = nullptr) {
    const auto &x = fields[0];
    auto out = Json::object();
    auto g = fields[1], c = fields[2], h = fields[3];
    require(numeric::add(numeric::add(g, c), h) == x && op.metric().inner(g, c) == 0 &&
                op.metric().inner(g, h) == 0 && op.metric().inner(c, h) == 0 &&
                numeric::is_zero(op.apply(h)),
            "trajectory Hodge identity");
    const auto green =
        prepared_green ? *prepared_green : field::exact_green_apply(op, lambda, x).value;
    require(numeric::add(green, numeric::scale(lambda, op.apply(green))) == x,
            "trajectory Green residual");
    out.set("input", vec(x));
    out.set("gradient", vec(g));
    out.set("curl", vec(c));
    out.set("harmonic", vec(h));
    out.set("green", vec(green));
    out.set("gradient_quadrance", q(op.metric().inner(g, g)));
    out.set("curl_quadrance", q(op.metric().inner(c, c)));
    out.set("harmonic_quadrance", q(op.metric().inner(h, h)));
    return out;
}
Json sectors(const field::ExactHodgeOperator &op, const field::HodgeProjectors &p, const Vector &x,
             const Rational &lambda) {
    return sectors_from_fields(
        op, {x, p.project_gradient(x), p.project_curl(x), p.project_harmonic(x)}, lambda);
}
Json trajectory(const Inputs &in, const Json &query, bool validate) {
    members(query,
            {"operation", "series", "selection", "support", "time", "faces", "semantics", "lambda",
             "minimum_witnesses", "endpoint_fields"},
            {"controls"});
    require(query.at("support").text() == "joint_reported" &&
                query.at("faces").text() == "coordinate_exchange_identities" &&
                query.at("semantics").text() == "observed_feature_change_relations",
            "trajectory declarations");
    auto series = query.at("series");
    require(json_object_is_type(series.raw(), json_type_array) && series.size() >= 2 &&
                series.size() <= 6,
            "trajectory requires 2..6 channels");
    std::vector<std::vector<std::string>> fs;
    std::vector<std::string> complete;
    std::set<std::string> labels;
    for (std::size_t f = 0; f < series.size(); ++f) {
        auto s = series.at(f);
        members(s, {"id", "fields"});
        require(labels.insert(text(s.at("id"))).second, "duplicate trajectory channel");
        fs.push_back(names(s.at("fields"), 8));
        require(fs.back().size() >= 2 && fs.back().size() == fs.front().size(),
                "trajectory visit cardinality");
        const auto unit = in.meaning.at(fs.back()[0]).at("units").text();
        for (const auto &name : fs.back()) {
            require(in.meaning.at(name).at("units").text() == unit,
                    "trajectory units differ across time");
            complete.push_back(name);
        }
    }
    const auto time = query.at("time");
    members(time, {"axis", "unit", "positions"});
    text(time.at("axis"));
    text(time.at("unit"));
    require(time.at("positions").size() == fs.front().size(), "trajectory time cardinality");
    Rational previous = exact_decimal(time.at("positions").at(0).text());
    for (std::size_t i = 1; i < fs.front().size(); ++i) {
        auto t = exact_decimal(time.at("positions").at(i).text());
        require(t > previous, "trajectory times must increase");
        previous = t;
    }
    const auto selection = names(query.at("selection"), 32, true),
               outcomes = names(query.at("endpoint_fields"), 8, true);
    const auto controls =
        query.has("controls") ? names(query.at("controls"), 32, true) : std::vector<std::string>{};
    for (const auto &control : controls) {
        require(std::find(outcomes.begin(), outcomes.end(), control) == outcomes.end(),
                "endpoint cannot control itself");
        complete.push_back(control);
    }
    const auto rows = in.select(complete, selection);
    const auto lambda = exact_decimal(query.at("lambda").text());
    require(lambda > 0, "Green lambda must be positive");
    const auto minimum = count(query.at("minimum_witnesses"), 4096, "minimum witnesses");
    require(minimum >= 2, "repeated support threshold must be at least two");
    if (validate)
        return Json{};
    Matrix declared_controls(rows.size(), controls.size() + 1);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        declared_controls.set(i, 0, Rational{1});
        for (std::size_t c = 0; c < controls.size(); ++c)
            declared_controls.set(i, c + 1, in.value(controls[c], rows[i]));
    }
    const auto control_frame = declared_controls.image_frame();
    const auto patient_metric = field::ExactMetric::identity(rows.size());
    const auto intervals = fs.front().size() - 1;
    std::vector<std::vector<Pattern>> patterns(intervals);
    std::vector<std::vector<Vector>> changes(intervals);
    std::set<Pattern> observed;
    for (std::size_t t = 0; t < intervals; ++t)
        for (auto row : rows) {
            Pattern pattern;
            Vector delta;
            for (const auto &channel : fs) {
                const Rational d = in.value(channel[t + 1], row) - in.value(channel[t], row);
                delta.push_back(d);
                pattern.push_back(d < 0 ? 0 : d == 0 ? 1 : 2);
            }
            observed.insert(pattern);
            patterns[t].push_back(pattern);
            changes[t].push_back(delta);
        }
    Vocabulary vocabulary(observed);
    auto result = Json::object(), states = Json::array(), faces = Json::array(),
         events = Json::array();
    for (const auto &p : vocabulary.patterns)
        states.push(indices(p));
    for (const auto &f : vocabulary.faces)
        faces.push(indices(f));
    result.set("control_frame", matrix(control_frame));
    result.set("control_rank", cardinality(control_frame.cols()));
    result.set("states", states);
    result.set("faces", faces);
    result.set("entities", in.identity(rows));
    result.set("n", cardinality(rows.size()));
    result.set("excluded_rows", cardinality(in.entities.size() - rows.size()));
    std::vector<algebra::MaterializedSubcomplex> blocks;
    std::vector<std::vector<std::size_t>> counts;
    std::vector<std::vector<Vector>> fields;
    auto readings = Json::array();
    for (std::size_t t = 0; t < intervals; ++t) {
        counts.emplace_back(vocabulary.patterns.size(), 0);
        for (const auto &p : patterns[t])
            ++counts.back()[vocabulary.lookup.at(p)];
        blocks.push_back(vocabulary.realize(vocabulary.mask(counts.back())));
        const auto &block = blocks.back();
        const auto ms = metrics(block.complex);
        field::ExactHodgeOperator op(block.complex, core::Grade{1}, ms);
        const auto projectors = field::exact_hodge_projectors(op);
        fields.emplace_back(series.size() + 1, Vector(op.size(), Rational{0}));
        auto sources = Json::array();
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto state = vocabulary.lookup.at(patterns[t][i]);
            const auto local = block.from_ambient[1][state];
            fields.back()[0][local] += 1;
            for (std::size_t f = 0; f < series.size(); ++f)
                fields.back()[f + 1][local] += changes[t][i][f];
            auto event = Json::object();
            event.set("entity", Json::string(in.entities[rows[i]]));
            event.set("interval", cardinality(t));
            event.set("state", cardinality(state));
            event.set("changes", vec(changes[t][i]));
            events.push(event);
            auto witness = Json::array();
            witness.push(Json::string(in.entities[rows[i]]));
            witness.push(cardinality(state));
            sources.push(witness);
        }
        auto record = topology_record(block.complex);
        record.set("state_indices", indices(block.to_ambient[1]));
        record.set("face_indices", indices(block.to_ambient[2]));
        record.set("witnesses", sources);
        record.set("interval", cardinality(t));
        auto stable = vocabulary.realize(vocabulary.mask(counts.back(), minimum));
        const auto cycles =
            field::materialize_boundary_matrix(stable.complex, core::Grade{1}).kernel_frame();
        record.set("repeated_support", topology_record(stable.complex));
        const auto stable_images =
            image_on_harmonic(cycles, inclusion(stable, block), projectors).image_frame();
        record.set("repeated_support_H1_image_rank", cardinality(stable_images.cols()));
        record.set("repeated_support_harmonic_frame", matrix(stable_images));
        auto field_records = Json::array();
        for (std::size_t f = 0; f < fields.back().size(); ++f) {
            auto r = sectors(op, projectors, fields.back()[f], lambda);
            r.set("id", f == 0 ? Json::string("patient_count") : series.at(f - 1).at("id"));
            r.set("units",
                  f == 0 ? Json::string("patients") : in.meaning.at(fs[f - 1][0]).at("units"));
            field_records.push(r);
        }
        for (const auto &outcome : outcomes) {
            Vector y(op.size(), Rational{0});
            bool full = true;
            auto known = Json::array();
            for (std::size_t i = 0; i < rows.size(); ++i) {
                const auto state = vocabulary.lookup.at(patterns[t][i]);
                if (in.fields.at(outcome)->tensor.support().contains(rows[i])) {
                    y[block.from_ambient[1][state]] += in.value(outcome, rows[i]);
                    known.push(Json::string(in.entities[rows[i]]));
                } else
                    full = false;
            }
            auto r = full ? sectors(op, projectors, y, lambda) : Json::object();
            r.set("id", Json::string(outcome));
            r.set("complete_support", Json::boolean(full));
            r.set("observed_entities", known);
            if (full) {
                auto pairings = Json::array();
                auto hy = projectors.project_harmonic(y);
                for (const auto &x : fields.back())
                    pairings.push(q(numeric::dot(projectors.project_harmonic(x), hy)));
                r.set("harmonic_pairings", pairings);
                const auto patient_y = in.values(outcome, rows);
                const auto residual_y = numeric::sub(
                    patient_y, field::projector_apply(control_frame, patient_metric, patient_y));
                require(numeric::is_zero(control_frame.transpose().apply(residual_y)),
                        "controlled topology endpoint orthogonality");
                Vector state_y(op.size(), Rational{0});
                for (std::size_t i = 0; i < rows.size(); ++i)
                    state_y[block.from_ambient[1][vocabulary.lookup.at(patterns[t][i])]] +=
                        residual_y[i];
                auto adjusted = sectors(op, projectors, state_y, lambda);
                adjusted.set("patient_input", vec(patient_y));
                adjusted.set("patient_residual", vec(residual_y));
                const auto harmonic_y = projectors.project_harmonic(state_y);
                Vector evaluations(rows.size(), Rational{0});
                for (std::size_t i = 0; i < rows.size(); ++i)
                    evaluations[i] =
                        harmonic_y[block.from_ambient[1][vocabulary.lookup.at(patterns[t][i])]];
                adjusted.set("patient_harmonic_evaluation", vec(evaluations));
                const auto repeated_y =
                    field::projector_apply(stable_images, op.metric(), harmonic_y);
                adjusted.set("repeated_support_harmonic", vec(repeated_y));
                adjusted.set("repeated_support_quadrance",
                             q(op.metric().inner(repeated_y, repeated_y)));
                auto adjusted_fields = Json::array();
                for (std::size_t f = 0; f < series.size(); ++f) {
                    Vector patient_x(rows.size(), Rational{0}), state_x(op.size(), Rational{0});
                    for (std::size_t i = 0; i < rows.size(); ++i)
                        patient_x[i] = changes[t][i][f];
                    const auto rx =
                        numeric::sub(patient_x, field::projector_apply(control_frame,
                                                                       patient_metric, patient_x));
                    require(numeric::is_zero(control_frame.transpose().apply(rx)),
                            "controlled topology feature orthogonality");
                    for (std::size_t i = 0; i < rows.size(); ++i)
                        state_x[block.from_ambient[1][vocabulary.lookup.at(patterns[t][i])]] +=
                            rx[i];
                    auto feature = Json::object();
                    feature.set("id", series.at(f).at("id"));
                    feature.set("patient_residual", vec(rx));
                    const auto hx = projectors.project_harmonic(state_x);
                    feature.set("harmonic", vec(hx));
                    feature.set("pairing", q(op.metric().inner(hx, harmonic_y)));
                    feature.set(
                        "repeated_support_pairing",
                        q(op.metric().inner(field::projector_apply(stable_images, op.metric(), hx),
                                            repeated_y)));
                    adjusted_fields.push(feature);
                }
                adjusted.set("fields", adjusted_fields);
                r.set("controlled", adjusted);
            }
            field_records.push(r);
        }
        record.set("fields", field_records);
        readings.push(record);
    }
    auto comparisons = Json::array();
    for (std::size_t t = 0; t + 1 < intervals; ++t) {
        const auto &a = blocks[t], &b = blocks[t + 1];
        auto shared =
            vocabulary.realize(meet(vocabulary.mask(counts[t]), vocabulary.mask(counts[t + 1])));
        auto ma = metrics(a.complex), mb = metrics(b.complex);
        field::ExactHodgeOperator oa(a.complex, core::Grade{1}, ma),
            ob(b.complex, core::Grade{1}, mb);
        auto pa = field::exact_hodge_projectors(oa), pb = field::exact_hodge_projectors(ob);
        const auto z =
            field::materialize_boundary_matrix(shared.complex, core::Grade{1}).kernel_frame();
        const auto ia = inclusion(shared, a), ib = inclusion(shared, b),
                   ha = image_on_harmonic(z, ia, pa), hb = image_on_harmonic(z, ib, pb);
        auto r = Json::object();
        r.set("from", cardinality(t));
        r.set("to", cardinality(t + 1));
        r.set("shared", topology_record(shared.complex));
        r.set("shared_states", indices(shared.to_ambient[1]));
        r.set("shared_cycle_frame", matrix(z));
        r.set("left_harmonic_images", matrix(ha));
        r.set("right_harmonic_images", matrix(hb));
        const auto joint = Matrix::vertical_concat(ha, hb).rank();
        r.set("left_image_rank", cardinality(ha.rank()));
        r.set("right_image_rank", cardinality(hb.rank()));
        r.set("joint_image_rank", cardinality(joint));
        r.set("common_quotient_dimension", cardinality(ha.rank() + hb.rank() - joint));
        Matrix j(b.to_ambient[1].size(), a.to_ambient[1].size());
        for (std::size_t i = 0; i < a.to_ambient[1].size(); ++i)
            if (b.from_ambient[1][a.to_ambient[1][i]] != algebra::RelationalMorphism::absent)
                j.set(b.from_ambient[1][a.to_ambient[1][i]], i, Rational{1});
        const auto ba = field::materialize_boundary_matrix(a.complex, core::Grade{1}),
                   bb = field::materialize_boundary_matrix(b.complex, core::Grade{1});
        const auto defect = bb.multiply(j).add(ba.scale(Rational{-1}));
        r.set("retained_state_correspondence", matrix(j));
        r.set("boundary_defect", matrix(defect));
        r.set("boundary_compatible", Json::boolean(defect.entries().empty()));
        auto transport = Json::array();
        for (std::size_t f = 0; f < fields[t].size(); ++f) {
            auto tr = Json::object();
            tr.set("id", f == 0 ? Json::string("patient_count") : series.at(f - 1).at("id"));
            const auto old_h = pa.project_harmonic(fields[t][f]);
            auto destinations = Json::array();
            for (const auto sector : {temporal::HodgeSector::Gradient, temporal::HodgeSector::Curl,
                                      temporal::HodgeSector::Harmonic})
                destinations.push(vec(temporal::sector_transport_apply(
                    pa, pb, j, sector, temporal::HodgeSector::Harmonic, fields[t][f])));
            tr.set("old_harmonic_destination_sectors", destinations);
            tr.set("transported_old_harmonic", vec(j.apply(old_h)));
            const auto change = temporal::green_response_delta(
                oa.materialize(), ob.materialize(), j, Matrix::identity(oa.size()),
                Matrix::identity(ob.size()), j, fields[t][f], fields[t + 1][f], lambda);
            require(change.aggregate ==
                        numeric::sub(change.new_response, j.apply(change.old_response)),
                    "finite Green transport identity");
            tr.set("realization_term", vec(change.realization_term));
            tr.set("innovation_term", vec(change.innovation_term));
            tr.set("operator_term", vec(change.operator_term));
            tr.set("total", vec(change.aggregate));
            tr.set("identity_exact", Json::boolean(true));
            transport.push(tr);
        }
        r.set("fields", transport);
        comparisons.push(r);
    }
    result.set("events", events);
    result.set("intervals", readings);
    result.set("temporal_comparisons", comparisons);
    result.set(
        "interpretation",
        Json::string("The carrier records measured feature-change states and declared exchange "
                     "identities. Shared-state inclusions support homology comparisons; a "
                     "defective full retained-state correspondence is not an induced homology map. "
                     "Neither H1 nor harmonic amplitudes assert a biological feedback mechanism."));
    return result;
}
} // namespace
Json compare_research_structures(const NativeStudy &old_study, const NativeStudy &new_study,
                                 core::Grade grade,
                                 const std::vector<std::vector<std::uint32_t>> &old_coordinates,
                                 const std::vector<std::vector<std::uint32_t>> &new_coordinates,
                                 const Vector &old_field, const Vector &new_field,
                                 const Rational &lambda, const std::vector<Matrix> *declared_maps,
                                 bool sparse_output) {
    scientific_detail::PreparedTemporalStudy old_prepared(old_study), new_prepared(new_study);
    return research_detail::compare_prepared_structures(
        old_prepared, new_prepared, grade, old_coordinates, new_coordinates, old_field, new_field,
        lambda, declared_maps, sparse_output);
}
Json research_detail::compare_prepared_structures(
    scientific_detail::PreparedTemporalStudy &old_prepared,
    scientific_detail::PreparedTemporalStudy &new_prepared, core::Grade grade,
    const std::vector<std::vector<std::uint32_t>> &old_coordinates,
    const std::vector<std::vector<std::uint32_t>> &new_coordinates, const Vector &old_field,
    const Vector &new_field, const Rational &lambda, const std::vector<Matrix> *declared_maps,
    bool sparse_output, bool factored_output) {
    const auto &old_study = old_prepared.study(), &new_study = new_prepared.study();
    require(lambda > 0, "structure Green lambda");
    const auto encode = [&](const Matrix &m) {
        if (!sparse_output)
            return matrix(m);
        auto out = Json::object(), entries = Json::array();
        out.set("rows", cardinality(m.rows()));
        out.set("columns", cardinality(m.cols()));
        for (const auto &[p, v] : m.entries()) {
            auto e = Json::object();
            e.set("row", cardinality(p.first));
            e.set("column", cardinality(p.second));
            e.set("value", Json::string(v.get_str()));
            entries.push(e);
        }
        out.set("entries", entries);
        return out;
    };
    auto correspondence = [&](std::size_t k) {
        if (declared_maps)
            return declared_maps->at(k);
        Matrix j(new_coordinates.at(k).size(), old_coordinates.at(k).size());
        std::map<std::uint32_t, std::size_t> target;
        for (std::size_t i = 0; i < new_coordinates.at(k).size(); ++i)
            target.emplace(new_coordinates[k][i], i);
        for (std::size_t i = 0; i < old_coordinates.at(k).size(); ++i)
            if (target.contains(old_coordinates[k][i]))
                j.set(target.at(old_coordinates[k][i]), i, Rational{1});
        return j;
    };
    const auto j = correspondence(grade.value);
    const auto old_hodge = old_study.hodge(grade), new_hodge = new_study.hodge(grade);
    // Reuse the existing source-bound owners. Their harmonic complement
    // avoids a second full harmonic frame and repeated Gram factorizations.
    const auto &old_p = old_prepared.sectors(grade), &new_p = new_prepared.sectors(grade);
    auto &old_green_plan = old_prepared.green(grade, lambda);
    auto &new_green_plan = new_prepared.green(grade, lambda);
    const auto old_green = factored_output ? old_green_plan.apply(old_field)
                                           : old_green_plan.matrix().apply(old_field);
    const auto new_green = factored_output ? new_green_plan.apply(new_field)
                                           : new_green_plan.matrix().apply(new_field);
    research_detail::FactoredTemporalMaps registry({&old_prepared, &new_prepared});
    auto out = Json::object();
    out.set("correspondence", encode(j));
    if (grade.value) {
        const auto old_b = field::materialize_boundary_matrix(old_study.source(), grade);
        const auto new_b = field::materialize_boundary_matrix(new_study.source(), grade);
        auto defect = new_b.multiply(j).add(
            correspondence(grade.value - 1).multiply(old_b).scale(Rational{-1}));
        out.set("boundary_defect", encode(defect));
        out.set("boundary_compatible", Json::boolean(defect.entries().empty()));
    }
    out.set("metric_pullback_defect",
            encode(temporal::metric_pullback_defect(old_hodge.metric(), new_hodge.metric(), j)));
    const auto old_sectors = old_p.apply(old_field);
    out.set("old", sectors_from_fields(old_hodge, old_sectors, lambda, &old_green));
    out.set("new", sectors_from_fields(new_hodge, new_p.apply(new_field), lambda, &new_green));
    auto transport = Json::array();
    Vector reconstructed(new_field.size(), Rational{0});
    for (std::size_t a = 1; a < old_sectors.size(); ++a) {
        auto row = Json::array();
        const auto transported = new_p.apply(j.apply(old_sectors[a]));
        for (std::size_t b = 1; b < transported.size(); ++b) {
            row.push(vec(transported[b]));
            reconstructed = numeric::add(reconstructed, transported[b]);
        }
        transport.push(row);
    }
    require(reconstructed == j.apply(old_field), "complete nine-sector transport");
    if (factored_output) {
        auto blocks = Json::array();
        const auto j_map = registry.matrix(j);
        for (std::size_t a = 0; a < 3; ++a) {
            auto row = Json::array();
            for (std::size_t b = 0; b < 3; ++b)
                row.push(registry.compose(
                    {registry.sector(0, grade, a), j_map, registry.sector(1, grade, b)}));
            blocks.push(row);
        }
        out.set("sector_maps", blocks);
    } else if (sparse_output) {
        std::vector<std::vector<Matrix>> maps(3,
                                              std::vector<Matrix>(3, Matrix(j.rows(), j.cols())));
        for (std::size_t c = 0; c < j.cols(); ++c) {
            Vector e(j.cols(), Rational{0});
            e[c] = 1;
            const auto projected = old_p.apply(e);
            for (std::size_t a = 1; a < projected.size(); ++a) {
                // A certified zero sector has zero images in every target
                // sector. Its complete declared block remains in the output.
                if (numeric::is_zero(projected[a]))
                    continue;
                const auto transported = new_p.apply(j.apply(projected[a]));
                for (std::size_t b = 1; b < transported.size(); ++b)
                    for (std::size_t i = 0; i < transported[b].size(); ++i)
                        if (transported[b][i] != 0)
                            maps[a - 1][b - 1].set(i, c, transported[b][i]);
            }
        }
        auto blocks = Json::array();
        Matrix total(j.rows(), j.cols());
        for (const auto &source_row : maps) {
            auto row = Json::array();
            for (const auto &block : source_row) {
                row.push(encode(block));
                total = total.add(block);
            }
            blocks.push(row);
        }
        require(total.entries() == j.entries(), "complete nine-sector map reconstruction");
        out.set("sector_maps", blocks);
    }
    out.set("sector_transport", transport);
    out.set("old_transported", vec(j.apply(old_field)));
    out.set("field_innovation", vec(temporal::field_delta(j, old_field, new_field).innovation));
    // Identity injections and complete endpoint Green actions. The ordered
    // native word defect keeps rectangular transport and every new direction.
    const auto innovation = temporal::field_delta(j, old_field, new_field).innovation;
    const Vector realization(new_field.size(), Rational{0});
    Vector innovation_term, operator_term;
    if (factored_output) {
        innovation_term = new_green_plan.apply(innovation);
        operator_term = numeric::sub(new_green_plan.apply(j.apply(old_field)), j.apply(old_green));
        out.set("green_operator_defect",
                registry.defect(registry.hodge(0, grade), registry.hodge(1, grade),
                                registry.matrix(j), registry.matrix(j)));
        out.set("factored_maps", registry.archive());
    } else {
        const auto &old_g = old_green_plan.matrix(), &new_g = new_green_plan.matrix();
        const auto defect = temporal::finite_word_defect(temporal::OperationWord({old_g}),
                                                         temporal::OperationWord({new_g}), {j, j});
        const auto hodge_defect = temporal::finite_word_defect(
            temporal::OperationWord({old_prepared.hodge(grade)}),
            temporal::OperationWord({new_prepared.hodge(grade)}), {j, j});
        require(defect.exact_identity() && hodge_defect.exact_identity(),
                "native complete Green and Hodge word maps");
        innovation_term = new_g.apply(innovation);
        operator_term = defect.direct.apply(old_field);
        out.set("green_operator_defect", encode(hodge_defect.direct));
    }
    const auto aggregate = numeric::add(innovation_term, operator_term);
    require(aggregate == numeric::sub(new_green, j.apply(old_green)),
            "native Green word and field reconstruction");
    out.set("realization_term", vec(realization));
    out.set("innovation_term", vec(innovation_term));
    out.set("operator_term", vec(operator_term));
    out.set("total", vec(aggregate));
    out.set("identity_exact", Json::boolean(true));
    return out;
}
bool is_research_operation(const std::string &op) {
    return op == "controlled_pairings" || op == "trajectory_structure" ||
           op == "paired_field_contrasts" || op == "linear_model_comparison" ||
           op == "finite_reference_calibration";
}
std::set<std::string> research_field_ids(const Json &query) {
    if (query.at("operation").text() == "finite_reference_calibration")
        return research_detail::finite_reference_field_ids(query);
    if (query.at("operation").text() == "linear_model_comparison")
        return research_detail::linear_model_field_ids(query);
    std::set<std::string> out;
    for (const auto *key_ : {"fields", "controls", "selection", "endpoint_fields"})
        if (query.has(key_))
            for (const auto &s : names(query.at(key_),
                                       query.at("operation").text() != "trajectory_structure"
                                           ? query.at(key_).size()
                                           : 64,
                                       true))
                out.insert(s);
    if (query.has("endpoint"))
        out.insert(text(query.at("endpoint")));
    if (query.has("series")) {
        const auto series = query.at("series");
        require(json_object_is_type(series.raw(), json_type_array) && series.size() <= 6,
                "trajectory series budget");
        for (std::size_t i = 0; i < series.size(); ++i)
            for (const auto &s : names(series.at(i).at("fields"), 8))
                out.insert(s);
    }
    if (query.has("benchmark")) {
        out.insert(text(query.at("benchmark").at("effect")));
        out.insert(text(query.at("benchmark").at("adjusted_p")));
    }
    if (query.has("pairs")) {
        const auto pairs = query.at("pairs");
        require(json_object_is_type(pairs.raw(), json_type_array),
                "paired fields must be an array");
        for (std::size_t i = 0; i < pairs.size(); ++i) {
            out.insert(text(pairs.at(i).at("left")));
            out.insert(text(pairs.at(i).at("right")));
        }
    }
    require(!out.empty() &&
                (query.at("operation").text() != "trajectory_structure" || out.size() <= 128),
            "research field budget");
    return out;
}
NativeRecord research_detail::comparison_fields(const ScientificStudy &study, const Json &query,
    const std::shared_ptr<CpuResources> &resources, const runtime::ExecutionLease &lease) {
    const auto operation = text(query.at("operation"));
    require(operation == "paired_field_contrasts" || operation == "controlled_pairings",
            "declared comparison operation");
    Inputs input(study, research_field_ids(query), false);
    auto result = operation == "paired_field_contrasts" ? paired(study, input, query, false, resources, lease)
                                                      : controlled(study, input, query, false, resources, lease);
    return std::move(result).value();
}
Json execute_research_query(const ScientificStudy &study, const Json &query, bool validate_only,
                            std::shared_ptr<CpuResources> resources) {
    const auto op = query.at("operation").text();
    require(is_research_operation(op), "unsupported research operation");
    if (op == "finite_reference_calibration")
        return research_detail::finite_reference_calibration(study, query, validate_only,
                                                             resources);
    if (op == "linear_model_comparison")
        return research_detail::linear_model_comparison(study, query, validate_only, resources);
    Inputs inputs(study, research_field_ids(query), op == "trajectory_structure");
    if (op == "paired_field_contrasts" || op == "controlled_pairings") {
        const auto lease = source_execution_lease("scientific-declared-comparison");
        auto result = op == "paired_field_contrasts" ? paired(study, inputs, query, validate_only, resources, lease)
                                                   : controlled(study, inputs, query, validate_only, resources, lease);
        return result ? result_fields_view(*result) : Json{};
    }
    return trajectory(inputs, query, validate_only);
}
} // namespace gadget
