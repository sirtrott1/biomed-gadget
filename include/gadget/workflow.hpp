#pragma once
#include <concepts>
#include <functional>
#include <gadget/resources.hpp>
#include <gadget/scientific.hpp>

namespace gadget {
class EvidenceStore;
struct WorkflowModel {
    std::size_t workers{1};
    std::function<Json(const std::string &, std::shared_ptr<CpuResources>)> invoke;
    template <class F>
        requires std::invocable<F, const std::string &>
    WorkflowModel(F f)
        : invoke([f = std::move(f)](const std::string &s, std::shared_ptr<CpuResources>) {
              return f(s);
          }) {}
    WorkflowModel(std::size_t n, decltype(invoke) f) : workers(n), invoke(std::move(f)) {}
};
Json check_certificate_claims(const Json &certificate, const Json &claims);
// Prospective, typed model claims; numerical agreement is checked separately.
Json admit_certificate_claims(const Json &proposal, const std::string &source_digest);
using WorkflowModels = std::map<std::string, WorkflowModel>;
// Biomedical control flow only. Each numerical leaf delegates to ScientificSession.
Json execute_workflow(std::shared_ptr<const ScientificStudy> study, const Json &workflow,
                      const WorkflowModels &models = {},
                      std::shared_ptr<CpuResources> resources = {}, ComputeOptions compute = {},
                      const EvidenceStore *evidence = nullptr);
NativeRecord execute_workflow_record(std::shared_ptr<const ScientificStudy> study,
    const NativeRecord &workflow,const WorkflowModels &models = {},
    std::shared_ptr<CpuResources> resources = {},ComputeOptions compute = {},
    const EvidenceStore *evidence = nullptr);
} // namespace gadget
