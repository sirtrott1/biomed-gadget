#pragma once
#include <gadget/io.hpp>
#include <gadget/native_record.hpp>
#include <gadget/resources.hpp>
#include <rfc/inference/native_linear_observation.hpp>
#include <rfc/inference/native_selected_observation.hpp>
#include <variant>

namespace gadget::scientific_detail {
struct NativeObservationBinding {
    NativeRecord source;
    rfc::storage::PinnedNativeAccession operation, accession;
    std::vector<rfc::storage::PinnedNativeAccession> fields;
    Json declarations;
};
class NativeObservationPlans {
    using Base = rfc::inference::PreparedNativeWoodburySolve::Base;
    using Plan = std::variant<rfc::inference::PreparedNativeLinearObservation,
                              rfc::inference::PreparedNativeSelectedObservation>;
    struct BaseRecord {
        Base owner;
    };
    std::map<std::string, BaseRecord> bases_;
    std::map<std::string, Plan> plans_;
    std::size_t applications_{};

  public:
    static void validate(const Json &query, const NativeObservationBinding &binding);
    NativeRecord execute_fields(const Json &query, const NativeObservationBinding &binding,
                              rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json execute(const Json &query, const NativeObservationBinding &binding,
                 rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json statistics() const;
};
} // namespace gadget::scientific_detail
