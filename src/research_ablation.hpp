#pragma once
#include <gadget/research.hpp>
#include <gadget/native_record.hpp>
#include <map>
#include "scientific_result_fields.hpp"
namespace gadget::research_detail {
class AblationFields;
struct SourceAblationRole {
    std::vector<rfc::core::CellId> cells;
    std::vector<std::string> fields;
    scientific_detail::ResultValue witness;
};
scientific_detail::ResultValue source_structure_ablation(const NativeStudy &study, const Json &query, rfc::core::Grade grade,
                               const std::map<std::string, Vector> &fields,
                               const std::map<std::string, SourceAblationRole> &roles,
                               const scientific_detail::ResultValue &source_witness, const Json &sheaf, bool validate,
                               const std::shared_ptr<CpuResources> &resources,
                               const NativeRecord &source_publication,
                               const std::map<std::string, Vector> &sheaf_operands,
                               AblationFields *publication);
} // namespace gadget::research_detail
