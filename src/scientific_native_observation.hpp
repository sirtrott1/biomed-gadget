#pragma once
#include <gadget/io.hpp>
#include <gadget/resources.hpp>
#include <rfc/inference/native_linear_observation.hpp>
#include <rfc/inference/native_selected_observation.hpp>
#include <variant>

namespace gadget::scientific_detail {
struct NativeObservationBinding {
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
        Json ordering;
    };
    std::map<std::string, BaseRecord> bases_;
    std::map<std::string, Plan> plans_;
    std::size_t applications_{};

  public:
    static void validate(const Json &query, const NativeObservationBinding &binding);
    Json execute(const Json &query, const NativeObservationBinding &binding,
                 rfc::runtime::ExecutionLease lease, CpuResources &cpu);
    Json statistics() const;
};
} // namespace gadget::scientific_detail
