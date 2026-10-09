#pragma once
#include <gadget/research.hpp>
#include <map>
namespace gadget::research_detail {
std::set<std::string> response_model_field_ids(const Json &query);
Json observed_response_model(const NativeStudy &study, const Json &query, rfc::core::Grade grade,
                             const std::map<std::string, Vector> &fields,
                             const std::map<std::string, Json> &meanings, const Json &coordinates,
                             bool validate, const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
