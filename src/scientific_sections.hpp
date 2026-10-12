#pragma once
#include "scientific_result_fields.hpp"
#include <gadget/study.hpp>
#include <map>

namespace gadget::scientific_detail {
struct LocalSectionAnalysis {
    ResultValue view;
    std::vector<rfc::numeric::ExactMatrix> coboundaries;
    std::vector<Vector> glued;
    std::vector<std::uint64_t> grade_carriers;
};
LocalSectionAnalysis analyze_local_section_fields(const NativeStudy &study, const Json &query,
    ResultFields *publication, const std::map<std::string, Vector> *operands = nullptr,
    bool validate_only = false, bool include_cochain_maps = false);
} // namespace gadget::scientific_detail
