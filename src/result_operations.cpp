#include "pipeline_result_binding.hpp"
#include "archive_identity.hpp"
#include <gadget/feature_context.hpp>
#include <gadget/result_operations.hpp>
#include <gadget/statistics/pathways.hpp>
#include <gadget/record_view.hpp>
#include <rfc/runtime/execution_policy.hpp>

namespace gadget {
NativeRecord read_analysis_record(const fs::path &path, const NativeRecord &request) {
    const auto spec = pipeline_detail::admit_result_operation(render_record_view(request));
    const auto image = read_text(path, spec.analysis_bytes);
    return NativeRecord::decode(std::as_bytes(std::span(image.data(), image.size())));
}
NativeRecord execute_result_record(const NativeRecord &analysis,
                                             const NativeRecord &request,
                                             const EvidenceStore *evidence,
                                             std::shared_ptr<CpuResources> resources,
                                             std::stop_token stop) {
    const auto cancelled = [&] {
        if (stop.stop_requested())
            throw rfc::runtime::ExecutionCancelled{};
    };
    cancelled();
    const auto spec = pipeline_detail::admit_result_operation(render_record_view(request));
    const auto image = pipeline_detail::verified_record_image(spec.source_path, spec.source_sha256,
                                                              spec.source_bytes);
    const auto source = NativeRecord::decode(std::as_bytes(std::span(image.data(), image.size())));
    const auto analysis_digest = archive_identity(analysis.owning_publication(), spec.analysis_bytes, stop);
    cancelled();
    auto query = request.at("query");
    if (spec.analysis_binding == "actual_result") {
        auto bound = query.unbound_value();
        std::get<NativeRecordValue::Object>(bound.data)
            .members.emplace_back(
                "analysis_sha256",
                NativeRecordValue(
                    analysis_digest));
        query = query.derive(std::move(bound));
    }
    NativeRecord result = [&] {
        if (spec.operation == "pathway_overrepresentation")
            return statistics::pathway_overrepresentation(analysis, source, query, resources, stop);
        if (!evidence)
            throw std::invalid_argument(
                "feature evidence context requires an actual EvidenceStore owner");
        return contextualize_feature_evidence(analysis, *evidence, source, query, resources,
                                                 stop);
    }();
    cancelled();
    // Preserve the dispatch and source file declaration in addition to the producer's
    // complete actual analysis, released source and biological request owners.
    return annotate_operation(result, request, result.at("result").unbound_value());
}
} // namespace gadget
