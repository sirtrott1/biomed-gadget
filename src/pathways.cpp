#include <gadget/pathways.hpp>

namespace gadget {
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources) {
    return statistics::pathway_overrepresentation(analysis, gene_sets, request, std::move(resources));
}
NativeRecord pathway_overrepresentation(const NativeRecord &analysis, const NativeRecord &gene_sets,
                                        const NativeRecord &request,
                                        std::shared_ptr<CpuResources> resources,
                                        std::stop_token stop) {
    return statistics::pathway_overrepresentation(analysis, gene_sets, request, std::move(resources), stop);
}
} // namespace gadget
