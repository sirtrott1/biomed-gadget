#pragma once
#include <memory>
#include <rfc/inference/native_gluing.hpp>
#include <rfc/runtime/cpu_tiles.hpp>
#include <rfc/runtime/resource_budget.hpp>
namespace gadget {
struct NativeGluingSection {
    rfc::core::Grade grade;
    std::uint64_t field;
    rfc::storage::PinnedNativeAccession left, right;
};
struct NativeGluingResult {
    rfc::inference::PreparedNativeGluing owner;
    std::vector<std::vector<rfc::inference::DerivedNativeGluedChain>> sections;
    rfc::runtime::CpuTileRunReport execution;
};
// Output identities are authored explicitly; source grades and complete field
// support come from the retained native chain, never from presentation rows.
struct NativeGluingTensorBinding {
    std::size_t section{}, lineage{};
    std::uint64_t tensor{}, carrier{};
};
rfc::storage::TemporalEndpoint native_gluing_endpoint(
    const NativeGluingResult &result, const std::vector<NativeGluingTensorBinding> &bindings);
// This application layer coordinates owning canonical actions. It performs no
// arithmetic, boundary construction, catalogue generation or JSON conversion.
NativeGluingResult execute_native_gluing(
    rfc::inference::PreparedNativeGluing plan, std::vector<NativeGluingSection> sections,
    rfc::inference::NativeGluingAvailability availability, bool boundaries,
    rfc::runtime::CpuResourceCoordinator::Reservation reservation,
    rfc::runtime::ExecutionLease lease);
} // namespace gadget
