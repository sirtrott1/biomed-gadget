#pragma once
#include <gadget/io.hpp>
#include <rfc/inference/native_green_response.hpp>

namespace gadget::scientific_detail {
struct DeclaredResponseBinding {
    rfc::storage::PinnedNativeAccession old_operator, new_operator, old_realization,
        new_realization, old_field, new_field;
    std::optional<rfc::storage::PinnedNativeAccession> old_accession, new_accession, metric;
    Json publications;
};
// Stores native preparations only. Each application binds the current input
// pins and rechecks its host lease; neither numerical results nor evidence are cached.
class DeclaredResponsePlans {
    std::map<std::string, rfc::inference::PreparedNativeGreenResponseDelta> green_;
    std::map<std::string, rfc::inference::PreparedNativeObservedGreenResponseDelta> observed_;
    std::size_t applications_{};

  public:
    static void validate(const Json &query, const DeclaredResponseBinding &binding);
    Json execute(const Json &query, const DeclaredResponseBinding &binding,
                 rfc::runtime::ExecutionLease lease);
    Json statistics() const;
};
rfc::runtime::ExecutionLease declared_response_lease();
} // namespace gadget::scientific_detail
