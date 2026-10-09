#pragma once
#include <gadget/research.hpp>

namespace gadget::research_detail {
std::set<std::string> linear_model_field_ids(const Json &query);
Json linear_model_comparison(const ScientificStudy &study, const Json &query, bool validate,
                             const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
