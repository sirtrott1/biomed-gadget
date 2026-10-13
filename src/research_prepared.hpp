#pragma once
#include "scientific_temporal_prepare.hpp"
#include "scientific_result_fields.hpp"

namespace gadget::research_detail {
inline void admit_structure_matrix(const rfc::numeric::ExactMatrix &matrix,
                                   const rfc::runtime::ExecutionLease &lease) {
    std::size_t bits=0;
    for (const auto &[position,value] : matrix.entries()) {
        (void)position;
        bits=rfc::machine::checked_size_add(bits,rfc::storage::native_execution_detail::integer_bits(value));
    }
    rfc::runtime::admit_execution({1,rfc::machine::checked_size_mul(matrix.rows(),matrix.cols()),
        matrix.nonzeros(),0,bits},lease.policy());
    lease.check_cancelled();
}
struct StructuralResultFields {
    scientific_detail::ResultFields &publication;
    std::vector<std::uint64_t> old_axes, new_axes;
    rfc::runtime::ExecutionLease lease;
    void admit(const Rational &value) const {
        rfc::runtime::admit_execution({1,1,value!=0,0,
            rfc::storage::native_execution_detail::integer_bits(value)},lease.policy());
        lease.check_cancelled();
    }
    void admit(const rfc::numeric::ExactMatrix &matrix) const {
        admit_structure_matrix(matrix,lease);
    }
    void admit(const Vector &vector) const {
        std::size_t bits=0,entries=0;
        for (const auto &value : vector) {
            bits=rfc::machine::checked_size_add(bits,rfc::storage::native_execution_detail::integer_bits(value));
            entries+=value!=0;
        }
        rfc::runtime::admit_execution({1,vector.size(),entries,0,bits},lease.policy());
        lease.check_cancelled();
    }
};
scientific_detail::ResultValue compare_prepared_structure_fields(
    scientific_detail::PreparedTemporalStudy &old_prepared,
    scientific_detail::PreparedTemporalStudy &new_prepared, rfc::core::Grade grade,
    const std::vector<std::vector<std::uint32_t>> &old_coordinates,
    const std::vector<std::vector<std::uint32_t>> &new_coordinates,
    const Vector &old_field, const Vector &new_field, const Rational &lambda,
    const std::vector<rfc::numeric::ExactMatrix> *declared_maps,
    bool factored_output, StructuralResultFields &fields);
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
