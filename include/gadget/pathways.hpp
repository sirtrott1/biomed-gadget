#pragma once
#include <gadget/native_record.hpp>
#include <stop_token>

namespace gadget {
class CpuResources;
// Select exact retained feature occurrences and explicitly mapped identifiers.
// The request pins both owning archives. The operation retains the unchanged
// molecular source/history, the original request and released source archive.
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources = {});
// Cancellation is observed between bounded native/archive and SDK stages;
// an individual native encoding or binomial action runs to completion.
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources,
                                        std::stop_token stop);
} // namespace gadget
