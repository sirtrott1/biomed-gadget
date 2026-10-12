#pragma once
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <rfc/topology/integral_homology.hpp>

namespace gadget::research_detail {
// Snapshot lattice certificates retain the canonical bases and transformations.
// Temporal diagrams require their separately declared modules and maps.
inline rfc::numeric::ExactMatrix lattice_matrix(const rfc::numeric::IntegerMatrix &value) {
    rfc::numeric::ExactMatrix result(value.rows(), value.cols());
    for (std::size_t i = 0; i < value.rows(); ++i)
        for (std::size_t j = 0; j < value.cols(); ++j)
            if (value.get(i, j) != 0)
                result.set(i, j, rfc::numeric::Rational{value.get(i, j)});
    return result;
}
inline scientific_detail::ResultValue lattice_snapshot(
    const rfc::core::RelationalComplex &source, rfc::core::Grade grade,
    scientific_detail::ResultFields &fields,
    const std::vector<std::uint64_t> &lattice_axes) {
    using namespace scientific_detail;
    using N = NativeRecordValue;
    const auto group = rfc::topology::integral_homology(source, grade);
    const auto chain_axis = lattice_axes.at(grade.value);
    const auto upper_axis = grade.value + 1u < lattice_axes.size()
        ? lattice_axes.at(grade.value + 1u) : fields.axis(0);
    const auto cycle_axis = fields.axis(group.cycle_rank);
    const auto smith_row = fields.axis(group.boundary_smith.diagonal.rows());
    const auto smith_column = fields.axis(group.boundary_smith.diagonal.cols());
    auto result = ResultValue::object();
    result.set("grade", ResultValue::count(grade.value));
    result.set("coefficient_domain", ResultValue::string("Z"));
    result.set("chain_rank", ResultValue::count(group.chain_rank));
    result.set("cycle_rank", ResultValue::count(group.cycle_rank));
    result.set("boundary_rank", ResultValue::count(group.boundary_rank));
    result.set("free_rank", ResultValue::count(group.free_rank));
    auto scales = ResultValue::array();
    for (const auto &value : source.grade(grade).cell_scales) scales.push(ResultValue(N(value)));
    result.set("cell_scales", std::move(scales));
    auto torsion = ResultValue::array(), invariants = ResultValue::array();
    for (const auto &value : group.torsion_invariants) torsion.push(ResultValue(N(value)));
    for (const auto &value : group.boundary_smith.invariant_factors)
        invariants.push(ResultValue(N(value)));
    result.set("torsion", std::move(torsion));
    result.set("invariant_factors", std::move(invariants));
    const auto cycles = lattice_matrix(group.cycle_basis);
    const auto boundaries = lattice_matrix(group.boundaries_in_cycle_coordinates);
    const auto left = lattice_matrix(group.boundary_smith.left);
    const auto right = lattice_matrix(group.boundary_smith.right);
    const auto diagonal = lattice_matrix(group.boundary_smith.diagonal);
    require(left.multiply(boundaries).multiply(right).entries() == diagonal.entries(),
            "complete native Smith transformation certificate");
    result.set("cycle_basis", ResultFields::stage_matrix(cycles, cycle_axis, chain_axis));
    result.set("boundaries_in_cycle_coordinates",
               ResultFields::stage_matrix(boundaries, upper_axis, cycle_axis));
    auto smith = ResultValue::object();
    smith.set("rank", ResultValue::count(group.boundary_smith.rank));
    smith.set("diagonal", ResultFields::stage_matrix(diagonal, smith_column, smith_row));
    smith.set("left", ResultFields::stage_matrix(left, cycle_axis, smith_row));
    smith.set("right", ResultFields::stage_matrix(right, smith_column, upper_axis));
    result.set("boundary_smith", std::move(smith));
    auto certificates = ResultValue::array();
    for (std::size_t column = 0; column < group.cycle_basis.cols(); ++column) {
        std::vector<rfc::numeric::Integer> chain;
        for (std::size_t row = 0; row < group.cycle_basis.rows(); ++row)
            chain.push_back(group.cycle_basis.get(row, column));
        const auto proof = rfc::topology::integral_chain_certificate(source, grade, chain);
        require(proof.cycle, "native integral cycle basis certificate");
        auto record = ResultValue::object();
        record.set("column", ResultValue::count(column));
        record.set("cycle", ResultValue::boolean(proof.cycle));
        record.set("boundary", ResultValue::boolean(proof.boundary));
        record.set("nonzero_class", ResultValue::boolean(proof.represents_nonzero_homology_class()));
        if (proof.boundary) {
            Vector preimage;
            for (const auto &value : proof.boundary_preimage) preimage.emplace_back(value);
            record.set("boundary_preimage", ResultFields::stage_vector(preimage, upper_axis));
        } else record.set("boundary_preimage", ResultValue{});
        certificates.push(std::move(record));
    }
    result.set("cycle_certificates", std::move(certificates));
    return result;
}
} // namespace gadget::research_detail
