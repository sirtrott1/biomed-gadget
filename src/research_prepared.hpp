#pragma once
#include "scientific_temporal_prepare.hpp"

namespace gadget::research_detail {
Json compare_prepared_structures(scientific_detail::PreparedTemporalStudy &old_prepared,
                                 scientific_detail::PreparedTemporalStudy &new_prepared,
                                 rfc::core::Grade grade,
                                 const std::vector<std::vector<std::uint32_t>> &old_coordinates,
                                 const std::vector<std::vector<std::uint32_t>> &new_coordinates,
                                 const Vector &old_field, const Vector &new_field,
                                 const Rational &lambda,
                                 const std::vector<rfc::numeric::ExactMatrix> *declared_maps,
                                 bool sparse_output, bool factored_output = false);
} // namespace gadget::research_detail
