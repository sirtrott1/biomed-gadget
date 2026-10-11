#pragma once
#include <gadget/research.hpp>
#include <map>
#include <optional>
#include <tuple>

namespace gadget::research_detail {
using LinearModelAxis = std::tuple<std::string, std::string, std::string>;
struct LinearModelFields {
    NativeRecord publication;
    std::map<LinearModelAxis, std::optional<Rational>> comparisons;
};
LinearModelFields linear_model_fields(const ScientificStudy &study, const Json &query,
                                      const std::shared_ptr<CpuResources> &resources);
std::set<std::string> linear_model_field_ids(const Json &query);
Json linear_model_comparison(const ScientificStudy &study, const Json &query, bool validate,
                             const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
