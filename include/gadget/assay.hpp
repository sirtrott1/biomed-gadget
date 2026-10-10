#pragma once
#include <gadget/scientific.hpp>
#include <gadget/workflow.hpp>
#include <functional>
#include <map>
#include <string_view>

namespace gadget {
// A portable collection of independently admitted assays and physical resources.
// Each assay's numerical storage is canonical RFCXS, not an application matrix format.
Json bioinformatics_capabilities();
Json inspect_bioinformatics_file(const fs::path &specification);
// Native physical records retain source spans and typed source-grade columns.
NativeRecord inspect_bioinformatics_record(const fs::path &specification);
void import_assays(const fs::path &specification, const fs::path &destination);
// Report completed import stages. An observer exception aborts the import and
// removes its unpublished destination, as for other admission failures.
using AssayImportObserver = std::function<void(std::string_view)>;
void import_assays(const fs::path &specification, const fs::path &destination,
                   const AssayImportObserver &observer);
class AssayWorkspace;
// Return the admitted workspace and its immutable native publication. Source
// owners used during import remain available for subsequent analyses.
std::shared_ptr<const AssayWorkspace> import_assay_workspace(
    const fs::path &specification, const fs::path &destination,
    const AssayImportObserver &observer = {});
class AssayWorkspace {
    struct State;
    std::unique_ptr<State> state_;
    using AdmittedSources =
        std::map<std::string, std::vector<std::shared_ptr<const ScientificStudy>>>;
    AssayWorkspace(const fs::path &manifest, AdmittedSources sources);
    AssayWorkspace(const fs::path &manifest, AdmittedSources sources, NativeRecord publication);
    friend void import_assays(const fs::path &, const fs::path &);
    friend void import_assays(const fs::path &, const fs::path &, const AssayImportObserver &);
    friend std::shared_ptr<const AssayWorkspace> import_assay_workspace(
        const fs::path &, const fs::path &, const AssayImportObserver &);
    NativeRecord execute_coverage(const NativeRecord &, std::shared_ptr<CpuResources>,
        std::optional<std::uint64_t> source_namespace = {}) const;
    NativeRecord execute_spatial(const NativeRecord &, std::shared_ptr<CpuResources>,
        std::optional<std::uint64_t> source_namespace = {}) const;
    Json materialize_native_fields(const NativeRecord &, const fs::path &,
        std::shared_ptr<CpuResources>,
        std::shared_ptr<const ScientificStudy> *retained = nullptr) const;
    Json replay_native_derivation(const std::string &, const fs::path &,
        std::shared_ptr<CpuResources>) const;

  public:
    explicit AssayWorkspace(const fs::path &manifest);
    static NativeRecord publish_studies(const NativeRecord &request, const std::string &id,
        const std::map<std::string, std::shared_ptr<const ScientificStudy>> &studies);
    Json materialize_record(const NativeRecord &request, const fs::path &destination,
        std::shared_ptr<CpuResources> resources = {}) const;
    ~AssayWorkspace();
    AssayWorkspace(const AssayWorkspace &) = delete;
    AssayWorkspace &operator=(const AssayWorkspace &) = delete;
    Json catalog() const;
    const NativeRecord &native_publication() const;
    void write_native(const fs::path &path) const;
    std::string digest() const;
    // A researcher-selected view of source identities and complete quoted spans.
    Json source_context(const Json &request) const;
    // Read the selected actual native partition through its original owner.
    Json read_source(const Json &request) const;
    // Bind explicitly selected admitted assays/resources without identifying
    // their independent native carriers or rebuilding their numerical packets.
    static Json
    bind_inputs(const std::map<std::string, std::shared_ptr<const AssayWorkspace>> &inputs,
                const Json &request, const fs::path &destination);
    Json execute(const Json &request, std::shared_ptr<CpuResources> resources = {}) const;
    NativeRecord execute_record(const NativeRecord &request,
                                std::shared_ptr<CpuResources> resources = {}) const;
    Json native_novelty_preparation_statistics() const;
    void clear_native_novelty_preparations();
    // Produces a portable derived assay with its complete generating certificate
    // and used inputs. Physical neighborhoods append to existing authored relations.
    Json materialize(const Json &request, const fs::path &destination,
                     std::shared_ptr<CpuResources> resources = {}) const;
    Json replay_derivation(const std::string &assay, const fs::path &destination,
                           std::shared_ptr<CpuResources> resources = {}) const;
    Json derivation_certificate(const std::string &assay) const;
    // Writes bounded per-partition certificates for large assays, then seals a
    // complete result index. Failed runs do not leave a published result directory.
    Json run(const Json &request, const fs::path &destination,
             std::shared_ptr<CpuResources> resources = {}) const;
    NativeRecord run_record(const NativeRecord &request, const fs::path &destination,
                            std::shared_ptr<CpuResources> resources = {}) const;
    // Selects an actual admitted native partition for the existing model/work-graph interface.
    std::shared_ptr<const ScientificStudy> study(const std::string &assay,
                                                 std::size_t partition = 0) const;
};
// Allocation classes are retained multiway observations. This profile solves
// independent nonnegative allocation simplexes; it does not fit abundances.
Json analyze_allocations(const Json &request, std::shared_ptr<CpuResources> resources = {});
// Models select declared context targets and propose interpretations. The
// receipt verifies owners and bytes; interpretations remain prospective.
Json admit_source_annotations(const AssayWorkspace &workspace, const Json &context,
                              const Json &proposal, const std::string &question);
} // namespace gadget
