#include "scientific_native_hodge.hpp"
#include "scientific_detail.hpp"
#include <charconv>
#include <rfc/runtime/cpu_tiles.hpp>
#include <set>

namespace gadget::scientific_detail {
namespace {
using namespace rfc;
using Pin = storage::PinnedNativeAccession;
using Matrix = numeric::ExactMatrix;
Json vector(const numeric::RationalVector &values) {
    auto result = Json::array();
    for (const auto &value : values)
        result.push(Json::string(value.get_str()));
    return result;
}
Json indices(const std::vector<std::size_t> &values) {
    auto result = Json::array();
    for (auto value : values)
        result.push(cardinality(value));
    return result;
}
Json matrix(const Matrix &value) {
    auto result = Json::object(), entries = Json::array();
    result.set("rows", cardinality(value.rows()));
    result.set("columns", cardinality(value.cols()));
    for (const auto &[position, coefficient] : value.entries()) {
        auto entry = Json::object();
        entry.set("row", cardinality(position.first));
        entry.set("column", cardinality(position.second));
        entry.set("value", Json::string(coefficient.get_str()));
        entries.push(entry);
    }
    result.set("entries", entries);
    return result;
}
Json reduction(const numeric::PreparedExactFrames &plan) {
    auto result = Json::object();
    result.set("original", matrix(plan.original()));
    result.set("image", matrix(plan.image()));
    result.set("coefficients", matrix(plan.coordinates()));
    result.set("kernel", matrix(plan.kernel()));
    result.set("pivot_columns", indices(plan.pivots()));
    result.set("free_columns", indices(plan.free_columns()));
    return result;
}
Json requirements(const runtime::ExecutionRequirements &value) {
    auto result = Json::object();
    // Limit values may exceed the signed integer range of the JSON codec.
    result.set("records", Json::string(std::to_string(value.records)));
    result.set("rational_slots", Json::string(std::to_string(value.rational_slots)));
    result.set("sparse_entries", Json::string(std::to_string(value.sparse_entries)));
    result.set("scalar_work", Json::string(std::to_string(value.scalar_work)));
    result.set("integer_bits", Json::string(std::to_string(value.integer_bits)));
    return result;
}
const storage::TemporalCarrier &carrier(const Pin &pin, std::uint64_t id) {
    for (const auto &value : pin.endpoint()->carriers)
        if (value.key == id)
            return value;
    throw std::invalid_argument("Hodge frame carrier absent from actual publication");
}
inference::NativeHodgeFrameRequest request(const Json &query, const NativeHodgeBinding &binding) {
    members(query,
            {"operation", "publication", "source_sequence", "grade", "grade_carriers", "fields",
             "metrics"},
            {"workers"});
    require(query.at("operation").text() == "hodge_frame_decomposition", "wrong Hodge operation");
    (void)text(query.at("publication"));
    if (query.has("workers"))
        require(count(query.at("workers"), std::numeric_limits<std::size_t>::max(),
                      "native Hodge workers") > 0,
                "native Hodge workers must be positive");
    const auto sequence = query.at("source_sequence").text();
    require(sequence == std::to_string(binding.publication.sequence()),
            "Hodge source_sequence differs from actual publication");
    const auto endpoint = binding.publication.endpoint();
    const auto grade =
        count(query.at("grade"), std::numeric_limits<std::uint16_t>::max(), "Hodge grade");
    require(grade < endpoint->source.grade_count(), "Hodge grade outside source");
    source_array(query.at("grade_carriers"), "Hodge grade carriers");
    source_array(query.at("fields"), "Hodge fields");
    require(query.at("grade_carriers").size() == endpoint->source.grade_count() &&
                query.at("fields").size() > 0 && query.at("fields").size() == binding.fields.size(),
            "Hodge requires all grade carriers and nonempty ordered Fields");
    inference::NativeHodgeFrameRequest result;
    result.grade = core::Grade{static_cast<std::uint16_t>(grade)};
    std::vector<std::size_t> sizes;
    std::size_t dimension = 1;
    for (std::size_t k = 0; k < endpoint->source.grade_count(); ++k) {
        const auto id = key(query.at("grade_carriers").at(k));
        const auto &value = carrier(binding.publication, id);
        const auto n =
            endpoint->source.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count();
        require(value.kind == storage::TemporalCarrierKind::SourceGrade && value.grade.value == k &&
                    value.domain == numeric::IndexDomain::linear(n),
                "Hodge requires actual complete SourceGrade carriers in grade order");
        result.grade_carriers.push_back(id);
        sizes.push_back(n);
        dimension = std::max(dimension, n);
    }
    // Native harmonic reduction stacks both adjacent constraints, so its rows
    // can exceed every individual grade. Admission follows the real shape.
    const auto lower = grade ? sizes[grade - 1] : 0;
    const auto upper = grade + 1 < sizes.size() ? sizes[grade + 1] : 0;
    result.dimension_limit = std::max(dimension, machine::checked_size_add(lower, upper));
    const auto metrics = query.at("metrics");
    const auto kind = metrics.at("kind").text();
    if (kind == "unit") {
        members(metrics, {"kind"});
        require(binding.metrics.empty(), "unit Hodge metrics cannot consume stored tensors");
        result.metrics = inference::NativeMetricSelection::NamedUnit;
    } else {
        members(metrics, {"kind", "tensors"});
        require(kind == "diagonal" || kind == "positive", "unknown Hodge metric selection");
        source_array(metrics.at("tensors"), "Hodge metric family");
        require(metrics.at("tensors").size() == sizes.size() &&
                    binding.metrics.size() == sizes.size(),
                "Hodge requires the full declared metric family");
        result.metrics = kind == "diagonal" ? inference::NativeMetricSelection::DeclaredDiagonal
                                            : inference::NativeMetricSelection::DeclaredPositive;
        result.metric_tensors = binding.metrics;
        for (std::size_t k = 0; k < sizes.size(); ++k) {
            (void)text(metrics.at("tensors").at(k));
            const auto &pin = binding.metrics[k];
            require(pin.same_execution_authority(binding.publication) &&
                        pin.endpoint().get() == endpoint.get() &&
                        pin.kind() == storage::NativeAccessionKind::Tensor,
                    "Hodge metric requires the actual publication's Tensor pin");
            const auto t = pin.tensor_declaration();
            require(t->role == storage::TemporalTensorRole::Metric &&
                        t->input == result.grade_carriers[k] && t->output == t->input &&
                        t->tensor.domain().extents() ==
                            std::vector<std::size_t>{sizes[k], sizes[k]} &&
                        t->tensor.support().count() == t->tensor.domain().logical_size(),
                    "Hodge metric requires complete support on its declared grade");
        }
    }
    for (std::size_t i = 0; i < binding.fields.size(); ++i) {
        (void)text(query.at("fields").at(i));
        const auto &pin = binding.fields[i];
        require(pin.same_execution_authority(binding.publication) &&
                    pin.endpoint().get() == endpoint.get() &&
                    pin.kind() == storage::NativeAccessionKind::Field,
                "Hodge Field requires its actual publication");
        const auto t = pin.tensor_declaration();
        require(t->role == storage::TemporalTensorRole::Field &&
                    t->output == result.grade_carriers[grade] &&
                    t->tensor.domain() == numeric::IndexDomain::linear(sizes[grade]) &&
                    t->tensor.support().count() == sizes[grade],
                "Hodge requires complete Field support on the selected actual grade");
    }
    return result;
}
} // namespace
void NativeHodgePlans::validate(const Json &query, const NativeHodgeBinding &binding) {
    (void)request(query, binding);
}
Json NativeHodgePlans::execute(const Json &query, const NativeHodgeBinding &binding,
                               runtime::ExecutionLease lease, CpuResources &cpu) {
    lease.check_cancelled();
    const auto selected = request(query, binding);
    const auto capacity = cpu.coordinator().capacity();
    const auto workers = query.has("workers")
                             ? count(query.at("workers"), capacity, "native Hodge workers")
                             : std::min(capacity, binding.fields.size());
    require(workers > 0, "native Hodge workers must be positive");
    auto reservation = cpu.reserve(workers);
    auto descriptor = Json::object();
    for (const auto *name :
         {"publication", "source_sequence", "grade", "grade_carriers", "metrics"})
        descriptor.set(name, query.at(name));
    const auto key = descriptor.dump_compact();
    const runtime::WorkAvailability available{binding.publication.source_namespace(),
                                              binding.publication.sequence()};
    // Commit a new preparation only after every requested action succeeds.
    // A rejected current lease cannot change the retained preparation statistics.
    auto existing = plans_.find(key);
    auto plan = existing == plans_.end() ? inference::PreparedNativeHodgeFrames::prepare(
                                               binding.publication, selected, available, lease)
                                         : existing->second;
    std::vector<Json> records(binding.fields.size());
    runtime::CpuTileExecutor executor({workers, capacity, binding.fields.size()});
    const auto work = executor.run_tiles(
        binding.fields.size(), [&](std::size_t i, std::size_t, const runtime::CpuTileStop &stop) {
            stop.check();
            const auto derived = plan.execute(binding.fields[i], available, lease);
            const auto &value = derived.values();
            const auto &report = derived.provenance();
            auto field = Json::object(), sectors = Json::object(), coordinates = Json::object(),
                 receipt = Json::object(), limits = Json::object();
            field.set("id", query.at("fields").at(i));
            field.set("input_field_key", Json::string(std::to_string(report.input_field)));
            field.set("meaning", binding.declarations.at("fields").at(i));
            numeric::RationalVector input(binding.fields[i].tensor()->domain().logical_size());
            for (const auto &e : binding.fields[i].tensor()->values().entries())
                input[e.index] = e.value;
            field.set("input", vector(input));
            sectors.set("gradient", vector(value.gradient));
            sectors.set("curl", vector(value.curl));
            sectors.set("harmonic", vector(value.harmonic));
            coordinates.set("gradient", vector(value.gradient_coordinates));
            coordinates.set("curl", vector(value.curl_coordinates));
            coordinates.set("harmonic", vector(value.harmonic_coordinates));
            field.set("sectors", sectors);
            field.set("frame_coordinates", coordinates);
            field.set("reconstructed", vector(value.reconstructed));
            receipt.set("execution_policy", Json::string(report.execution_policy));
            receipt.set("preparation_policy", Json::string(report.preparation_policy));
            receipt.set("preparation_source_sequence",
                        Json::string(std::to_string(report.preparation_availability.max_sequence)));
            receipt.set("execution_source_sequence",
                        Json::string(std::to_string(report.execution_availability.max_sequence)));
            receipt.set("execution_steps", cardinality(report.execution_steps));
            limits.set("dimension", Json::string(std::to_string(report.action_limits.dimension)));
            limits.set("entries", Json::string(std::to_string(report.action_limits.entries)));
            limits.set("steps", Json::string(std::to_string(report.action_limits.steps)));
            limits.set("integer_bits",
                       Json::string(std::to_string(report.action_limits.integer_bits)));
            receipt.set("action_limits", limits);
            receipt.set("requirements", requirements(report.selected_requirements));
            receipt.set("source_namespace", Json::string(std::to_string(report.source_namespace)));
            receipt.set("source_sequence", Json::string(std::to_string(report.sequence)));
            receipt.set("carrier", Json::string(std::to_string(report.carrier)));
            receipt.set("fallback", Json::boolean(report.fallback));
            field.set("native_execution", receipt);
            records[i] = std::move(field);
            stop.check();
        });
    auto fields = Json::array();
    for (const auto &field : records)
        fields.push(field);
    const auto &numeric_plan = plan.numeric_plan();
    auto frames = Json::object(), reductions = Json::object(), result = Json::object(),
         preparation = Json::object();
    frames.set("gradient", matrix(numeric_plan.gradient_frame()));
    frames.set("curl", matrix(numeric_plan.curl_frame()));
    frames.set("harmonic", matrix(numeric_plan.harmonic_frame()));
    reductions.set("gradient", reduction(numeric_plan.gradient_reduction()));
    reductions.set("curl", reduction(numeric_plan.curl_reduction()));
    reductions.set("harmonic", reduction(numeric_plan.harmonic_reduction()));
    preparation.set("steps", cardinality(plan.preparation_steps()));
    preparation.set("retained_requirements", requirements(plan.retained_requirements()));
    preparation.set("dimension_admission", cardinality(selected.dimension_limit));
    preparation.set("chain_interfaces_checked",
                    cardinality(numeric_plan.chain_interfaces_checked()));
    auto execution = Json::object();
    execution.set("reserved_workers", cardinality(workers));
    execution.set("completed_fields", cardinality(records.size()));
    execution.set("executor", Json::string("rfc::runtime::CpuTileExecutor"));
    (void)work; // Wall time remains a separate measurement, outside scientific certificates.
    result.set("native_work", execution);
    result.set("profile", Json::string("source-grade-complete-hodge-frames-v1"));
    result.set("publication", binding.declarations);
    result.set("grade", cardinality(selected.grade.value));
    result.set("grade_carriers", query.at("grade_carriers"));
    result.set("metrics", query.at("metrics"));
    result.set("frames", frames);
    result.set("reductions", reductions);
    result.set("fields", fields);
    result.set("native_preparation", preparation);
    lease.check_cancelled();
    if (existing == plans_.end())
        plans_.emplace(key, std::move(plan));
    applications_ = machine::checked_size_add(applications_, binding.fields.size());
    return Json::parse(result.dump());
}
Json NativeHodgePlans::statistics() const {
    auto result = Json::object();
    result.set("hodge_frame_preparations", cardinality(plans_.size()));
    result.set("field_applications", cardinality(applications_));
    return result;
}
} // namespace gadget::scientific_detail
