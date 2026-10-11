#pragma once
#include <gadget/io.hpp>
#include <gadget/resources.hpp>
#include <gadget/native_record.hpp>
#include <rfc/inference/native_metric_novelty.hpp>

namespace gadget::scientific_detail {
struct NativeNoveltyBinding {
    NativeRecord old_source, new_source;
    rfc::storage::PinnedNativeAccession old_field, new_field, metric;
    std::vector<rfc::storage::PinnedNativeAccession> causes;
    Json declarations;
};
// Only canonical immutable plans are retained. Input fields, output vectors,
// quotations and receipts are never cached. Clearing releases our plan owners.
class NativeNoveltyPlans {
    std::map<std::string, rfc::inference::PreparedNativeMetricNovelty> plans_;
    std::size_t applications_{}, preparations_{};

  public:
    static void validate(const Json &query, const NativeNoveltyBinding &binding);
    NativeRecord execute_fields(const Json &query, const NativeNoveltyBinding &binding,
                 rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json execute(const Json &query, const NativeNoveltyBinding &binding,
                 rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json statistics() const;
    void clear();
};
} // namespace gadget::scientific_detail
