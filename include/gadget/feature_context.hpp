#pragma once
#include <gadget/biomedical.hpp>
#include <gadget/resources.hpp>
#include <stop_token>

namespace gadget {
// Keep the unchanged analysis and the separate mapping/evidence owners. All
// biological correspondence is an explicit released or declared source map.
NativeRecord contextualize_feature_evidence(const NativeRecord &analysis,
    const EvidenceStore &evidence, const NativeRecord &mapping_source,
    const NativeRecord &request, std::shared_ptr<CpuResources> resources = {});
// Host cancellation stays outside the source authored request.
NativeRecord contextualize_feature_evidence(const NativeRecord &analysis,
    const EvidenceStore &evidence, const NativeRecord &mapping_source,
    const NativeRecord &request, std::shared_ptr<CpuResources> resources, std::stop_token stop);
} // namespace gadget
