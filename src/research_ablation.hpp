#pragma once
#include <gadget/research.hpp>
#include <map>
namespace gadget::research_detail {
struct SourceAblationRole {
    std::vector<rfc::core::CellId> cells;
    std::vector<std::string> fields;
    Json witness;
};
Json source_structure_ablation(const NativeStudy &study, const Json &query, rfc::core::Grade grade,
                               const std::map<std::string, Vector> &fields,
                               const std::map<std::string, SourceAblationRole> &roles,
                               const Json &source_witness, const Json &sheaf, bool validate,
                               const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
