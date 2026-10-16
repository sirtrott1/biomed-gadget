#pragma once
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <rfc/runtime/execution_policy.hpp>

namespace gadget {
// Each source is its complete original native publication at the named cutoff.
NativeRecord homology_diagram_record(const NativeRecord &request,
    const std::vector<NativeRecord> &sources, std::shared_ptr<CpuResources> resources = {},
    std::optional<rfc::runtime::ExecutionLease> lease = {});
NativeRecord run_homology_diagram(const fs::path &specification,
    std::shared_ptr<CpuResources> resources = {});
NativeRecord replay_homology_diagram(const NativeRecord &record,
    std::shared_ptr<CpuResources> resources = {},
    std::optional<rfc::runtime::ExecutionLease> lease = {});
Json homology_diagram_view(const NativeRecord &record);
} // namespace gadget
