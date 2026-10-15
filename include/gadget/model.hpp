#pragma once
#include <gadget/execution.hpp>
#include <gadget/scientific.hpp>

namespace gadget {
class AssayWorkspace;
class ScientificFacts;
class EvidenceStore;
struct LocalModelOptions {
    fs::path model;
    std::uint32_t context_tokens{4096}, output_tokens{256}, threads{8};
    std::int32_t gpu_layers{};
    bool cache_prompt{true};
    std::uint32_t prompt_cache_mib{256};
    bool compact_catalog{true};
    bool paged_catalog{false};
    std::string task{"planner"};
    ComputeOptions compute;
    fs::path projector;
};
// Admission stays independent of the model decoder and available without llama.
Json apply_model_proposal(QuerySession &session, const Json &proposal);
#if GADGET_HAS_LOCAL_MODEL
class LocalAgent {
    struct State;
    std::unique_ptr<State> state_;
    Json ask_assigned_impl(const std::string &prompt, const std::vector<fs::path> &images,
                          std::shared_ptr<CpuResources> resources,
                          std::optional<NativeRecord> *native_response);

  public:
    LocalAgent(std::shared_ptr<const Workspace> study, LocalModelOptions options);
    LocalAgent(std::shared_ptr<const ScientificStudy> study, LocalModelOptions options);
    LocalAgent(const Json &pipeline_catalog, LocalModelOptions options);
    // Accepts a pipeline intent or an explicit source-context request. Each
    // profile re-admits proposals against the same immutable owning workspace.
    LocalAgent(std::shared_ptr<const AssayWorkspace> workspace, const Json &intent,
               LocalModelOptions options);
    LocalAgent(std::shared_ptr<const ScientificFacts> facts, LocalModelOptions options);
    LocalAgent(std::shared_ptr<const EvidenceStore> evidence, const NativeRecord &context_request,
               LocalModelOptions options);
    ~LocalAgent();
    LocalAgent(const LocalAgent &) = delete;
    LocalAgent &operator=(const LocalAgent &) = delete;
    // Explicit publication switch; weights and decoding context remain loaded.
    // A failed switch retains the previous complete study and preparation.
    void bind(std::shared_ptr<const Workspace> study);
    void bind_scientific(std::shared_ptr<const ScientificStudy> study);
    // Drop the model's presentation receipt without reloading source, weights
    // or mathematical preparations. The caller still owns earlier transcripts.
    void clear_source_context();
    Json ask(const std::string &prompt);
    Json ask(const std::string &prompt, const std::vector<fs::path> &images);
    Json ask_assigned(const std::string &prompt, const std::vector<fs::path> &images,
                      std::shared_ptr<CpuResources> resources);
    // Scientific query and plan responses retain their actual native execution.
    // Other responses are source bound control records with no numerical import.
    NativeRecord ask_record(const std::string &prompt, const std::vector<fs::path> &images = {});
    NativeRecord ask_assigned_record(const std::string &prompt, const std::vector<fs::path> &images,
                                    std::shared_ptr<CpuResources> resources);
    // Latest admitted evidence proposal, retaining its actual source and request.
    // A rejected generation does not replace a previously admitted publication.
    std::optional<NativeRecord> prospective_publication() const;
};
#endif
} // namespace gadget
