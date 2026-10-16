#pragma once
#include <gadget/resources.hpp>
#include <rfc/backend/packed_binary.hpp>
#include <stop_token>

namespace gadget {
// Device policy belongs to the host, never to a model authored query.
struct ComputeOptions {
    std::size_t device_index{};
    bool allow_software_device{false};
    std::size_t allocation_byte_limit{256u * 1024u * 1024u};
};
struct PackedResult {
    std::vector<rfc::machine::PackedBits> outputs;
    Json execution;
};
// Retains the canonical immutable operation and its private calibration owner.
class PackedOperation {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit PackedOperation(rfc::backend::PreparedPackedBinary prepared,
                             ComputeOptions options = {});
    ~PackedOperation();
    rfc::backend::PreparedPackedBinary prepared() const;
    PackedResult execute(CpuResources &resources, const Json &query);
    // Host cancellation is checked at canonical CPU word boundaries. A call
    // reserves every output and its complete verification buffer until joined.
    PackedResult execute(CpuResources &resources, const Json &query, std::stop_token stop);
};
Json compute_devices();
} // namespace gadget
