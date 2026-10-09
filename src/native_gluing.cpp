#include <gadget/native_gluing.hpp>
#include <set>
namespace gadget {
rfc::storage::TemporalEndpoint native_gluing_endpoint(
    const NativeGluingResult &result, const std::vector<NativeGluingTensorBinding> &bindings) {
    using namespace rfc;
    storage::TemporalEndpoint endpoint(result.owner.constructed_source()->source);
    std::set<std::pair<std::size_t, std::size_t>> selected;
    std::set<std::uint64_t> tensors, carriers;
    std::size_t total = 0;
    for (const auto &lineage : result.sections) total = machine::checked_size_add(total, lineage.size());
    if (bindings.size() != total)
        throw std::invalid_argument("native gluing export requires every retained chain");
    for (const auto &binding : bindings) {
        if (!binding.tensor || !binding.carrier || !tensors.insert(binding.tensor).second ||
            !carriers.insert(binding.carrier).second ||
            !selected.emplace(binding.section, binding.lineage).second)
            throw std::invalid_argument("duplicate native gluing output identity");
        const auto &chain = result.sections.at(binding.section).at(binding.lineage).chain();
        const auto size = endpoint.source.grade(chain.grade()).cell_count();
        const auto domain = numeric::IndexDomain::linear(size);
        std::vector<std::size_t> support;
        std::vector<numeric::SparseRationalEntry> entries;
        for (const auto &[index, value] : chain.values()) {
            support.push_back(index);
            if (value != 0) entries.push_back({index, value});
        }
        endpoint.carriers.push_back({binding.carrier, storage::TemporalCarrierKind::SourceGrade,
                                     chain.grade(), {}, 0, domain});
        storage::TemporalTensor tensor;
        tensor.key = binding.tensor;
        tensor.role = storage::TemporalTensorRole::Field;
        tensor.output = binding.carrier;
        tensor.tensor = numeric::SparseRationalTensor(domain,
            numeric::StructuralSupport::sparse(size, std::move(support)),
            numeric::SparseRationalVector(size, std::move(entries)));
        endpoint.tensors.push_back(std::move(tensor));
    }
    storage::validate_temporal_endpoint(endpoint);
    return endpoint;
}
NativeGluingResult execute_native_gluing(
    rfc::inference::PreparedNativeGluing plan, std::vector<NativeGluingSection> inputs,
    rfc::inference::NativeGluingAvailability available, bool boundaries,
    rfc::runtime::CpuResourceCoordinator::Reservation reservation,
    rfc::runtime::ExecutionLease lease) {
    lease.check_cancelled();
    std::vector<std::vector<rfc::inference::DerivedNativeGluedChain>> results(inputs.size());
    rfc::runtime::CpuTileExecutor executor({reservation.workers(), reservation.workers(), inputs.size()});
    const auto execution = executor.run_tiles(inputs.size(),
        [&](std::size_t i, std::size_t, const rfc::runtime::CpuTileStop &stop) {
            lease.check_cancelled(); stop.check();
            const auto &input = inputs[i];
            auto &lineage = results[i];
            lineage.push_back(plan.glue_sections(input.field, input.grade, input.left, input.right,
                                                  available, lease));
            while (boundaries && lineage.back().chain().grade().value > 0) {
                lease.check_cancelled(); stop.check();
                lineage.push_back(plan.boundary(lineage.back(), available, lease));
            }
            lease.check_cancelled(); stop.check();
        });
    lease.check_cancelled();
    return {std::move(plan), std::move(results), execution};
}
} // namespace gadget
