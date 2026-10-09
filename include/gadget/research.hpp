#pragma once
#include <gadget/scientific.hpp>
#include <gadget/study.hpp>
#include <set>

namespace gadget {
// Research queries consume the same admitted immutable native study as other
// scientific operations. Endpoint roles, controls and temporal semantics are
// explicit query declarations; canonical RexGraph owns every algebraic action.
bool is_research_operation(const std::string &operation);
std::set<std::string> research_field_ids(const Json &query);
Json execute_research_query(const ScientificStudy &study, const Json &query,
                            bool validate_only = false,
                            std::shared_ptr<CpuResources> resources = {});
Json compare_research_structures(
    const NativeStudy &old_study, const NativeStudy &new_study, rfc::core::Grade grade,
    const std::vector<std::vector<std::uint32_t>> &old_coordinates,
    const std::vector<std::vector<std::uint32_t>> &new_coordinates, const Vector &old_field,
    const Vector &new_field, const Rational &lambda,
    const std::vector<rfc::numeric::ExactMatrix> *declared_maps = nullptr,
    bool sparse_output = false);
} // namespace gadget
