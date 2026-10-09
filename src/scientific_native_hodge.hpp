#pragma once
#include <gadget/io.hpp>
#include <gadget/resources.hpp>
#include <rfc/inference/native_hodge_frames.hpp>

namespace gadget::scientific_detail {
struct NativeHodgeBinding {
    rfc::storage::PinnedNativeAccession publication;
    std::vector<rfc::storage::PinnedNativeAccession> fields, metrics;
    Json declarations;
};
// Plans own actual native publications; results are always newly executed.
class NativeHodgePlans {
    std::map<std::string, rfc::inference::PreparedNativeHodgeFrames> plans_;
    std::size_t applications_{};

  public:
    static void validate(const Json &query, const NativeHodgeBinding &binding);
    Json execute(const Json &query, const NativeHodgeBinding &binding,
                 rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json statistics() const;
};
} // namespace gadget::scientific_detail
