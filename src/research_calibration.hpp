#pragma once
#include <gadget/research.hpp>

namespace gadget::research_detail {
NativeRecord finite_reference_fields(const ScientificStudy &study, const Json &query,
                                      const std::shared_ptr<CpuResources> &resources);
std::set<std::string> finite_reference_field_ids(const Json &query);
Json finite_reference_calibration(const ScientificStudy &study, const Json &query, bool validate,
                                  const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
