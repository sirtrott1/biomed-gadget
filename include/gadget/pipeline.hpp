#pragma once
#include <gadget/assay.hpp>

namespace gadget {
// Profiles describe assay-specific choices; a plan contains only selected steps.
Json pipeline_profiles();
Json prepare_pipeline(const AssayWorkspace &workspace, const Json &design);
Json default_pipeline(const AssayWorkspace &workspace, const Json &settings);
Json admit_pipeline_proposal(const Json &proposal, const std::string &source_digest);
Json pipeline_model_catalog(const AssayWorkspace &workspace);
// Intent bindings retain researcher-declared parameters. Models select names;
// prerequisite closure and compilation reuse the existing default recipes.
Json recipe_model_catalog(const AssayWorkspace &workspace, const Json &intent);
Json admit_recipe_proposal(const AssayWorkspace &workspace, const Json &intent,
                           const Json &proposal, const std::string &question);
// Every producing artifact is retained separately, with owning native jobs and
// an atomic result index. Derived workspaces become inputs of dependent steps.
Json execute_pipeline(std::shared_ptr<const AssayWorkspace> workspace, const Json &pipeline,
                      const fs::path &destination, const WorkflowModels &models = {},
                      std::shared_ptr<CpuResources> resources = {},
                      const EvidenceStore *evidence = nullptr, const Json *planning = nullptr);
NativeRecord execute_pipeline_record(std::shared_ptr<const AssayWorkspace> workspace,
    const NativeRecord &pipeline,const fs::path &destination,const WorkflowModels &models = {},
    std::shared_ptr<CpuResources> resources = {},const EvidenceStore *evidence = nullptr,
    const NativeRecord *planning = nullptr);
} // namespace gadget
