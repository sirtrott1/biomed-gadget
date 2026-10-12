#pragma once
#include <gadget/research.hpp>
#include <map>
#include "scientific_result_fields.hpp"
namespace gadget::research_detail {
struct ExchangeBlock {
    std::vector<std::size_t> coordinates;
    scientific_detail::ResultValue witness;
};
std::set<std::string> structural_calibration_field_ids(const Json &query);
scientific_detail::ResultValue structural_reference_calibration(const NativeStudy &study, const Json &query,
                                      rfc::core::Grade grade,
                                      const std::map<std::string, Vector> &fields,
                                      const std::vector<ExchangeBlock> &blocks, bool validate,
                                      const std::shared_ptr<CpuResources> &resources,
                                      scientific_detail::ResultFieldContext *context);
} // namespace gadget::research_detail
