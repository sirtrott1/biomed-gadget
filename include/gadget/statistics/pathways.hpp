#pragma once
#include <gadget/native_record.hpp>
#include <stop_token>

namespace gadget {
class CpuResources;
}

namespace gadget::statistics {
// Test mapped feature selections under an explicitly declared uniform subset
// law. Exact arithmetic evaluates the law; its scientific assumptions remain
// declarations. Both source archives and the original request remain retained.
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources = {});
// Cancellation is observed between bounded native/archive and SDK stages;
// an individual native encoding or binomial action runs to completion.
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources,
                                        std::stop_token stop);
} // namespace gadget::statistics
