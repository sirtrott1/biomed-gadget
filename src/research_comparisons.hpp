#pragma once
#include <gadget/research.hpp>
#include <rfc/runtime/execution_policy.hpp>

namespace gadget::research_detail {
// Each canonical action observes the declared host policy and cancellation.
NativeRecord comparison_fields(const ScientificStudy &study, const Json &query,
    const std::shared_ptr<CpuResources> &resources, const rfc::runtime::ExecutionLease &lease);
} // namespace gadget::research_detail
