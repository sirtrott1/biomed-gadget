#pragma once
#include <gadget/research.hpp>
#include <map>
namespace gadget::research_detail {
struct ExchangeBlock {
    std::vector<std::size_t> coordinates;
    Json witness;
};
std::set<std::string> structural_calibration_field_ids(const Json &query);
Json structural_reference_calibration(const NativeStudy &study, const Json &query,
                                      rfc::core::Grade grade,
                                      const std::map<std::string, Vector> &fields,
                                      const std::vector<ExchangeBlock> &blocks, bool validate,
                                      const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
