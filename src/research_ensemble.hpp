#pragma once
#include <gadget/research.hpp>

namespace gadget::research_detail {
// Fields have already been admitted on one actual context/grade/carrier by the
// scientific session. The query owns sample weights and comparison meanings.
Json relational_field_ensemble(const NativeStudy &study, const Json &query, rfc::core::Grade grade,
                               const std::vector<Vector> &fields, bool validate,
                               const std::shared_ptr<CpuResources> &resources);
} // namespace gadget::research_detail
