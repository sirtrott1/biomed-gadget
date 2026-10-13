#pragma once
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"

namespace gadget::scientific_detail {
Vector declared_chain_values(const Json &query, std::size_t size);
void validate_chain_values(const NativeStudy &study, const Json &query, const Vector &input);
ResultValue chain_fields(const NativeStudy &study, const Json &query, const Vector &input,
    ResultFields &fields, const std::vector<std::uint64_t> &axes,
    const rfc::runtime::ExecutionLease &lease, bool validate_only = false);
ResultValue path_fields(const NativeStudy &study, const Json &query,
    ResultFields &fields, const std::vector<std::uint64_t> &axes,
    const rfc::runtime::ExecutionLease &lease, bool validate_only = false);
} // namespace gadget::scientific_detail
