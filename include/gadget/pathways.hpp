#pragma once
#include <gadget/statistics/pathways.hpp>

namespace gadget {
// Compatibility entry points delegate to gadget::statistics. New consumers
// include gadget/statistics/pathways.hpp and select that namespace explicitly.
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources = {});
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources,
                                        std::stop_token stop);
} // namespace gadget
