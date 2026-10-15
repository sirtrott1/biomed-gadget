#pragma once
#include <concepts>
#include <functional>
#include <gadget/resources.hpp>
#include <gadget/scientific.hpp>

namespace gadget {
class EvidenceStore;
struct WorkflowModel {
    using NativeInvoke = std::function<NativeRecord(const std::string &, std::shared_ptr<CpuResources>)>;
    std::size_t workers{1};
    std::function<Json(const std::string &, std::shared_ptr<CpuResources>)> invoke;
    NativeInvoke invoke_native;
    template <class F>
        requires std::invocable<F, const std::string &>
    WorkflowModel(F f)
        : invoke([f = std::move(f)](const std::string &s, std::shared_ptr<CpuResources>) {
              return f(s);
          }) {}
    WorkflowModel(std::size_t n, decltype(invoke) f) : workers(n), invoke(std::move(f)) {}
    WorkflowModel(std::size_t n, NativeInvoke f) : workers(n), invoke_native(std::move(f)) {}
};
Json check_certificate_claims(const Json &certificate, const Json &claims);
// Compare control literals with an actual native operation's values. The
// request names its publication namespace and cutoff; tensor coordinates are
// original flat coordinates, independent of an optional display order.
NativeRecord check_certificate_claims_record(const NativeRecord &certificate,
                                            const NativeRecord &request);
// Prospective, typed model claims; numerical agreement is checked separately.
Json admit_certificate_claims(const Json &proposal, const std::string &source_digest);
using WorkflowModels = std::map<std::string, WorkflowModel>;
// Biomedical control flow delegates numerical leaves to owning native providers.
Json execute_workflow(std::shared_ptr<const ScientificStudy> study, const Json &workflow,
                      const WorkflowModels &models = {},
                      std::shared_ptr<CpuResources> resources = {}, ComputeOptions compute = {},
                      const EvidenceStore *evidence = nullptr);
// Compatibility display from the original typed declaration and byte owner.
Json execute_workflow(std::shared_ptr<const ScientificStudy> study, const NativeRecord &workflow,
                      const WorkflowModels &models = {},
                      std::shared_ptr<CpuResources> resources = {}, ComputeOptions compute = {},
                      const EvidenceStore *evidence = nullptr);
NativeRecord execute_workflow_record(std::shared_ptr<const ScientificStudy> study,
    const NativeRecord &workflow,const WorkflowModels &models = {},
    std::shared_ptr<CpuResources> resources = {},ComputeOptions compute = {},
    const EvidenceStore *evidence = nullptr);
// Optional presentation of retained publications; never used to execute a leaf
// or select an actual tensor coefficient.
Json workflow_record_view(const NativeRecord &publication);
} // namespace gadget
