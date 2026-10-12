#include "scientific_detail.hpp"
#include <chrono>
#include <gadget/compute.hpp>
#include <limits>
#if GADGET_HAS_VULKAN
#include <rfc/backend/mixed_binary_calibration.hpp>
#include <rfc/backend/vulkan_workspace_cache.hpp>
#endif

namespace gadget {
using namespace scientific_detail;
namespace {
Json tiles(const rfc::runtime::CpuTileRunReport &r) {
    auto j = Json::object();
    j.set("requested_workers", cardinality(r.requested_workers));
    j.set("worker_slots", cardinality(r.worker_slots));
    j.set("launched_threads", cardinality(r.launched_threads));
    j.set("tiles", cardinality(r.tiles));
    j.set("completed_tiles", cardinality(r.completed_tiles));
    return j;
}
#if GADGET_HAS_VULKAN
Json device_json(const rfc::backend::VulkanDeviceInfo &d) {
    auto j = Json::object();
    j.set("name", Json::string(d.name));
    j.set("physical_gpu", Json::boolean(d.physical_gpu));
    j.set("compute", Json::boolean(d.compute));
    j.set("vendor", cardinality(d.vendor));
    j.set("device", cardinality(d.device));
    j.set("driver_version", cardinality(d.driver_version));
    j.set("api_version", cardinality(d.api_version));
    auto uuid = Json::array();
    for (auto b : d.uuid)
        uuid.push(cardinality(b));
    j.set("uuid", uuid);
    auto heaps = Json::array();
    for (const auto &h : d.heaps) {
        auto v = Json::object();
        v.set("bytes", Json::string(std::to_string(h.bytes)));
        v.set("device_local", Json::boolean(h.device_local));
        v.set("advisory_budget",
              h.advisory_budget ? Json::string(std::to_string(*h.advisory_budget)) : Json{});
        heaps.push(v);
    }
    j.set("heaps", heaps);
    j.set("memory_scope",
          Json::string("heaps may be shared on UMA; these are not additive VRAM capacities"));
    return j;
}
#endif
} // namespace
struct PackedOperation::State {
    rfc::backend::PreparedPackedBinary owner;
    ComputeOptions options;
#if GADGET_HAS_VULKAN
    std::optional<rfc::backend::VulkanBinaryDevice> device;
    std::unique_ptr<rfc::backend::VulkanBinaryWorkspaceCache> workspace;
    std::map<std::size_t, std::unique_ptr<rfc::backend::MixedBinaryCalibration>> calibrations;
#endif
    State(rfc::backend::PreparedPackedBinary p, ComputeOptions o)
        : owner(std::move(p)), options(o) {}
};
PackedOperation::PackedOperation(rfc::backend::PreparedPackedBinary p, ComputeOptions o)
    : state_(std::make_unique<State>(std::move(p), o)) {}
PackedOperation::~PackedOperation() = default;
PackedResult PackedOperation::execute(CpuResources &resources, const Json &query) {
    return execute(resources, query, {});
}
PackedResult PackedOperation::execute(CpuResources &resources, const Json &query, std::stop_token stop) {
    using namespace rfc::backend;
    const auto backend = query.at("backend").text();
    require(backend == "cpu_tiles" || backend == "vulkan" || backend == "cpu_vulkan" ||
                backend == "measured",
            "unknown packed backend");
    auto &coordinator = resources.coordinator();
    const auto workers = query.has("workers")
                             ? count(query.at("workers"), coordinator.capacity(), "workers") : 1;
    const auto tile =
        query.has("tile_words") ? count(query.at("tile_words"),
                                        std::numeric_limits<std::size_t>::max(), "tile words") : 64;
    const auto &owner = state_->owner;
    require(workers && tile, "packed workers and tile words must be positive");
    if (query.has("gpu_words"))
        (void)count(query.at("gpu_words"), owner.words(), "GPU words");
    const auto start = std::chrono::steady_clock::now();
    const auto check_stop = [&] {
        if (stop.stop_requested()) throw rfc::runtime::ExecutionCancelled{};
    };
    check_stop();
    (void)resources.readmit();
    const auto bytes = rfc::machine::checked_size_mul(
        rfc::machine::checked_size_mul(owner.words(), owner.outputs().size()), sizeof(std::uint64_t));
    const auto verified_bytes = rfc::machine::checked_size_add(bytes, bytes);
    std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> reservation;
    PackedResult result;
    result.execution = Json::object();
    result.execution.set("backend", Json::string(backend));
    result.execution.set("declared_output_bytes", Json::string(std::to_string(bytes)));
    result.execution.set("declared_verification_bytes", Json::string(std::to_string(bytes)));
    if (backend == "cpu_tiles") {
        const auto num_tiles = rfc::machine::ceil_div(owner.words(), tile);
        reservation.emplace(
            coordinator.reserve(std::min(workers, std::max(std::size_t{1}, num_tiles)), verified_bytes));
        result.outputs = owner.empty_result(state_->options.allocation_byte_limit);
        rfc::runtime::CpuTileExecutor executor(
            {reservation->workers(), coordinator.capacity(), coordinator.policy().tile_limit});
        const auto run_start = std::chrono::steady_clock::now();
        auto report = executor.run_tiles(num_tiles, [&](auto i, auto, auto peers) {
            peers.check();
            check_stop();
            const auto begin = i * tile;
            owner.cpu_words(begin, begin + std::min(tile, owner.words() - begin), result.outputs, stop);
        });
        result.execution.set("cpu_run_ns", Json::string(std::to_string(
                                               std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - run_start)
                                                   .count())));
        result.execution.set("cpu_tiles", tiles(report));
        result.execution.set("gpu_words", cardinality(0));
        result.execution.set("cpu_words", cardinality(owner.words()));
    } else {
#if GADGET_HAS_VULKAN
        if (!state_->device)
            state_->device.emplace(VulkanBinaryOptions{state_->options.device_index,
                                                       state_->options.allow_software_device,
                                                       state_->options.allocation_byte_limit});
        MixedBinaryResult mixed = [&]() {
            if (backend == "measured") {
                const bool reused = state_->calibrations.contains(tile);
                if (!reused) {
                    std::vector<MixedBinaryCandidate> candidates;
                    for (auto w : {1u, 2u, 4u, 8u})
                        if (w <= coordinator.capacity()) {
                            candidates.push_back({0, w, tile});
                            if (owner.words() > 1)
                                candidates.push_back({owner.words() / 2, w, tile});
                        }
                    candidates.push_back({owner.words(), 1, tile});
                    state_->calibrations.emplace(tile, std::make_unique<MixedBinaryCalibration>(
                                                           owner, *state_->device, coordinator,
                                                           candidates, 3,
                                                           state_->options.allocation_byte_limit));
                }
                auto &calibration = *state_->calibrations.at(tile);
                auto samples = Json::array();
                for (const auto &s : calibration.samples()) {
                    auto j = Json::object();
                    j.set("gpu_words", cardinality(s.candidate.gpu_words));
                    j.set("cpu_workers", cardinality(s.candidate.cpu_workers));
                    j.set("median_complete_call_ns", Json::string(std::to_string(s.median_ns)));
                    auto raw = Json::array();
                    for (auto t : s.complete_call_ns)
                        raw.push(Json::string(std::to_string(t)));
                    j.set("complete_call_ns", raw);
                    samples.push(j);
                }
                result.execution.set("calibration_reused", Json::boolean(reused));
                result.execution.set("calibration", samples);
                check_stop();
                return calibration.execute(coordinator, stop);
            }
            if (!state_->workspace)
                state_->workspace = std::make_unique<VulkanBinaryWorkspaceCache>(
                    *state_->device,
                    VulkanWorkspaceCacheLimits{1, state_->options.allocation_byte_limit});
            const auto gpu = backend == "vulkan" ? owner.words()
                             : query.has("gpu_words")
                                 ? static_cast<std::size_t>(query.at("gpu_words").integer())
                                 : owner.words() / 2;
            auto executed = state_->workspace->execute(
                owner, coordinator,
                {gpu, workers, tile, state_->options.allocation_byte_limit, stop, {}});
            const auto stats = state_->workspace->stats();
            auto residency = Json::object();
            residency.set("creations", cardinality(stats.creations));
            residency.set("hits", cardinality(stats.hits));
            residency.set("discarded", cardinality(stats.discarded));
            residency.set("active", cardinality(stats.active));
            residency.set("idle", cardinality(stats.idle));
            residency.set("allocation_bytes", Json::string(std::to_string(stats.allocation_bytes)));
            residency.set("scope",
                          Json::string("canonical exclusive workspace loans; actual Vulkan buffer "
                                       "allocations, not total driver/process memory"));
            result.execution.set("workspace_cache", residency);
            return executed;
        }();
        result.outputs = std::move(mixed.outputs);
        result.execution.set("device", device_json(mixed.device));
        result.execution.set("gpu_words", cardinality(mixed.gpu_words));
        result.execution.set("cpu_words", cardinality(mixed.cpu_words));
        result.execution.set("cpu_tiles", tiles(mixed.cpu_report));
        auto timing = Json::object();
        timing.set("upload_ns", Json::string(std::to_string(mixed.timings.host_upload_ns)));
        timing.set("submit_ns", Json::string(std::to_string(mixed.timings.record_submit_ns)));
        timing.set("cpu_ns", Json::string(std::to_string(mixed.timings.cpu_ns)));
        timing.set("readback_ns",
                   Json::string(std::to_string(mixed.timings.completion_readback_ns)));
        timing.set("total_ns", Json::string(std::to_string(mixed.timings.total_ns)));
        result.execution.set("native_timings", timing);
        result.execution.set("rexgraph_backend_sources_sha256",
                             Json::string(GADGET_REXGRAPH_BACKEND_SHA256));
#else
        throw std::invalid_argument("Vulkan compute is disabled in this build");
#endif
    }
    // Canonical serial verification, separately timed; never silently substitute.
    check_stop();
    if (reservation)
        coordinator.resize(*reservation, 1, verified_bytes);
    else
        reservation.emplace(coordinator.reserve(1, verified_bytes));
    auto reference = owner.empty_result(state_->options.allocation_byte_limit);
    owner.cpu_words(0, owner.words(), reference, stop);
    require(reference.size() == result.outputs.size(), "packed output count mismatch");
    for (std::size_t i = 0; i < reference.size(); ++i)
        require(reference[i].words() == result.outputs[i].words(),
                "packed exact CPU reference mismatch");
    check_stop();
    result.execution.set("exact_cpu_reference_passed", Json::boolean(true));
    result.execution.set(
        "complete_verified_call_seconds",
        Json::number(
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()));
    return result;
}
Json compute_devices() {
    auto result = Json::array();
#if GADGET_HAS_VULKAN
    for (const auto &d : rfc::backend::enumerate_vulkan_devices())
        result.push(device_json(d));
#else
    throw std::invalid_argument("Vulkan compute is disabled in this build");
#endif
    return result;
}
} // namespace gadget
