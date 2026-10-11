#include "scientific_native_hodge.hpp"
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <charconv>
#include <rfc/runtime/cpu_tiles.hpp>
#include <set>

namespace gadget::scientific_detail {
namespace {
using namespace rfc;
using Pin = storage::PinnedNativeAccession;
using Matrix = numeric::ExactMatrix;
ResultValue reduction(ResultFields &fields, const numeric::PreparedExactFrames &plan,
                      std::uint64_t domain, std::uint64_t range,
                      std::uint64_t image, std::uint64_t kernel) {
    auto result = ResultValue::object();
    result.set("original", fields.matrix(plan.original(), domain, range));
    result.set("image", fields.matrix(plan.image(), image, range));
    result.set("coefficients", fields.matrix(plan.coordinates(), domain, image));
    result.set("kernel", fields.matrix(plan.kernel(), kernel, domain));
    result.set("pivot_columns", fields.indices(plan.pivots(), image));
    result.set("free_columns", fields.indices(plan.free_columns(), kernel));
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
NativeRecord NativeHodgePlans::execute_fields(const Json &query, const NativeHodgeBinding &binding,
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
    std::vector<std::optional<inference::DerivedNativeHodgeDecomposition>> derived_fields(binding.fields.size());
    runtime::CpuTileExecutor executor({workers, capacity, binding.fields.size()});
    (void)executor.run_tiles(
        binding.fields.size(), [&](std::size_t i, std::size_t, const runtime::CpuTileStop &stop) {
            stop.check();
            derived_fields[i].emplace(plan.execute(binding.fields[i], available, lease));
            stop.check();
        });
    ResultFields stored(binding.source);
    const auto source_carrier = plan.carrier();
    const auto &native_frames = plan.numeric_plan();
    const auto gradient_axis = stored.axis(native_frames.gradient_frame().cols());
    const auto curl_axis = stored.axis(native_frames.curl_frame().cols());
    const auto harmonic_axis = stored.axis(native_frames.harmonic_frame().cols());
    auto fields = ResultValue::array();
    for (std::size_t i = 0; i < derived_fields.size(); ++i) {
        const auto &derived = *derived_fields[i];
        const auto &value = derived.values();
        const auto &report = derived.provenance();
        auto field = ResultValue::object(), sectors = ResultValue::object(), coordinates = ResultValue::object();
        auto receipt = Json::object(), limits = Json::object();
        field.set("id", query.at("fields").at(i));
        field.set("input_field_key", ResultValue::string(std::to_string(report.input_field)));
        field.set("meaning", binding.declarations.at("fields").at(i));
        field.set("input", stored.stored(report.input_field, "vector"));
        sectors.set("gradient", stored.vector(value.gradient, source_carrier));
        sectors.set("curl", stored.vector(value.curl, source_carrier));
        sectors.set("harmonic", stored.vector(value.harmonic, source_carrier));
        coordinates.set("gradient", stored.vector(value.gradient_coordinates, gradient_axis));
        coordinates.set("curl", stored.vector(value.curl_coordinates, curl_axis));
        coordinates.set("harmonic", stored.vector(value.harmonic_coordinates, harmonic_axis));
        field.set("sectors", std::move(sectors));
        field.set("frame_coordinates", std::move(coordinates));
        field.set("reconstructed", stored.vector(value.reconstructed, source_carrier));
        receipt.set("execution_policy", Json::string(report.execution_policy));
        receipt.set("preparation_policy", Json::string(report.preparation_policy));
        receipt.set("preparation_source_sequence", Json::string(std::to_string(report.preparation_availability.max_sequence)));
        receipt.set("execution_source_sequence", Json::string(std::to_string(report.execution_availability.max_sequence)));
        receipt.set("execution_steps", cardinality(report.execution_steps));
        limits.set("dimension", Json::string(std::to_string(report.action_limits.dimension)));
        limits.set("entries", Json::string(std::to_string(report.action_limits.entries)));
        limits.set("steps", Json::string(std::to_string(report.action_limits.steps)));
        limits.set("integer_bits", Json::string(std::to_string(report.action_limits.integer_bits)));
        receipt.set("action_limits", limits);
        receipt.set("requirements", requirements(report.selected_requirements));
        receipt.set("source_namespace", Json::string(std::to_string(report.source_namespace)));
        receipt.set("source_sequence", Json::string(std::to_string(report.sequence)));
        receipt.set("carrier", Json::string(std::to_string(report.carrier)));
        receipt.set("fallback", Json::boolean(report.fallback));
        field.set("native_execution", receipt);
        fields.push(std::move(field));
    }
    const auto &numeric_plan = plan.numeric_plan();
    auto frames = ResultValue::object(), reductions = ResultValue::object(), result = ResultValue::object();
    auto preparation = Json::object();
    frames.set("gradient", stored.matrix(numeric_plan.gradient_frame(), gradient_axis, source_carrier));
    frames.set("curl", stored.matrix(numeric_plan.curl_frame(), curl_axis, source_carrier));
    frames.set("harmonic", stored.matrix(numeric_plan.harmonic_frame(), harmonic_axis, source_carrier));
    const auto lower_axis = selected.grade.value ? selected.grade_carriers[selected.grade.value - 1]
                                                 : stored.axis(0);
    const auto upper_axis = std::size_t{selected.grade.value} + 1 < selected.grade_carriers.size()
                             ? selected.grade_carriers[selected.grade.value + 1] : stored.axis(0);
    const auto gradient_kernel = stored.axis(numeric_plan.gradient_reduction().kernel().cols());
    const auto curl_kernel = stored.axis(numeric_plan.curl_reduction().kernel().cols());
    const auto constraint_axis = stored.axis(numeric_plan.harmonic_reduction().original().rows());
    const auto constraint_image = stored.axis(numeric_plan.harmonic_reduction().image().cols());
    reductions.set("gradient", reduction(stored, numeric_plan.gradient_reduction(), lower_axis,
                                          source_carrier, gradient_axis, gradient_kernel));
    reductions.set("curl", reduction(stored, numeric_plan.curl_reduction(), upper_axis,
                                      source_carrier, curl_axis, curl_kernel));
    reductions.set("harmonic", reduction(stored, numeric_plan.harmonic_reduction(), source_carrier,
                                          constraint_axis, constraint_image, harmonic_axis));
    auto binding_record = ResultValue::object();
    binding_record.set("kind", ResultValue::string("source_grade_hodge_frames"));
    binding_record.set("source_namespace", ResultValue::string(std::to_string(binding.publication.source_namespace())));
    binding_record.set("source_sequence", ResultValue::string(std::to_string(binding.publication.sequence())));
    binding_record.set("grade_carriers", query.at("grade_carriers"));
    binding_record.set("metrics", query.at("metrics"));
    binding_record.set("frames", frames);
    binding_record.set("reductions", reductions);
    stored.bind(std::move(binding_record));
    preparation.set("steps", cardinality(plan.preparation_steps()));
    preparation.set("retained_requirements", requirements(plan.retained_requirements()));
    preparation.set("dimension_admission", cardinality(selected.dimension_limit));
    preparation.set("chain_interfaces_checked",
                    cardinality(numeric_plan.chain_interfaces_checked()));
    auto execution = Json::object();
    execution.set("reserved_workers", cardinality(workers));
    execution.set("completed_fields", cardinality(derived_fields.size()));
    execution.set("executor", Json::string("rfc::runtime::CpuTileExecutor"));
    result.set("native_work", execution);
    result.set("profile", Json::string("source-grade-complete-hodge-frames-v1"));
    result.set("publication", binding.declarations);
    result.set("grade", ResultValue::count(selected.grade.value));
    result.set("grade_carriers", query.at("grade_carriers"));
    result.set("metrics", query.at("metrics"));
    result.set("frames", frames);
    result.set("reductions", reductions);
    result.set("fields", fields);
    result.set("native_preparation", preparation);
    const auto publication = stored.finish(std::move(result));
    lease.check_cancelled();
    if (existing == plans_.end())
        plans_.emplace(key, std::move(plan));
    applications_ = machine::checked_size_add(applications_, binding.fields.size());
    return publication;
}
Json NativeHodgePlans::execute(const Json &query, const NativeHodgeBinding &binding,
                               runtime::ExecutionLease lease, CpuResources &cpu) {
    return render_record_view(execute_fields(query, binding, std::move(lease), cpu));
}
Json NativeHodgePlans::statistics() const {
    auto result = Json::object();
    result.set("hodge_frame_preparations", cardinality(plans_.size()));
    result.set("field_applications", cardinality(applications_));
    return result;
}
} // namespace gadget::scientific_detail
