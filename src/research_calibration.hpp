#pragma once
#include <gadget/research.hpp>

namespace gadget::research_detail {
std::set<std::string> finite_reference_field_ids(const Json &query);
Json finite_reference_calibration(const ScientificStudy &study, const Json &query, bool validate,
                                  const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
