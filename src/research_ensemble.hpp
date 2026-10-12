#pragma once
#include <gadget/research.hpp>
#include "scientific_result_fields.hpp"

namespace gadget::research_detail {
// Fields have already been admitted on one actual context/grade/carrier by the
// scientific session. The query owns sample weights and comparison meanings.
scientific_detail::ResultValue relational_field_ensemble(const NativeStudy &study, const Json &query, rfc::core::Grade grade,
                               const std::vector<Vector> &fields, bool validate,
                               const std::shared_ptr<CpuResources> &resources,
                               scientific_detail::ResultFieldContext *publication);
} // namespace gadget::research_detail
