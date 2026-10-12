#pragma once
#include <gadget/compute.hpp>
#include <gadget/io.hpp>
#include <gadget/native_record.hpp>
#include <functional>
#include <string_view>
#include <rfc/runtime/execution_policy.hpp>
#include <rfc/storage/native_accession.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <rfc/storage/temporal_endpoint.hpp>

namespace rfc::storage { class NativeAtlasSource; }
namespace gadget {
class CpuResources;
// Reports producer availability. Source and query admission remain required.
bool has_native_scientific_publisher(std::string_view operation) noexcept;
// Select actual retained native values through archived scientific views.
// Publication and overlay traversal preserves the selected value's owner.
NativeRecord select_scientific_result(const NativeRecord &value, const NativeRecord &path);
// Recompute a retained model or finite calibration operation from its original
// study bytes and query, then verify its complete fields, carriers and support.
NativeRecord replay_scientific_record(const NativeRecord &publication,
    std::shared_ptr<CpuResources> resources = {});
class AssayWorkspace;
// Reports completed admission stages. Throwing aborts construction; no source
// publication or input file is changed by admission.
using ScientificAdmissionObserver = std::function<void(std::string_view)>;
// Biomedical admission over a canonical archive and explicitly pinned publication
// declarations. Each referenced archive owns its native endpoint. Mathematical source/storage
// and tensor ownership remain RexGraph's. Public JSON views are detached copies.
class ScientificStudy {
    struct State;
    std::unique_ptr<State> state_;
    explicit ScientificStudy(std::unique_ptr<State> state);
    friend class ScientificSession;
    friend void import_assays(const fs::path &, const fs::path &);
    friend void import_assays(const fs::path &, const fs::path &,
                              const std::function<void(std::string_view)> &);
    friend std::shared_ptr<const AssayWorkspace> import_assay_workspace(
        const fs::path &, const fs::path &, const std::function<void(std::string_view)> &);
    static std::shared_ptr<const ScientificStudy> admit_generated(const fs::path &,
        rfc::storage::NativeTemporalSeries history);

  public:
    explicit ScientificStudy(const fs::path &declaration);
    explicit ScientificStudy(NativeRecord publication);
    ScientificStudy(const fs::path &declaration, const ScientificAdmissionObserver &observer);
    ScientificStudy(NativeRecord publication, const ScientificAdmissionObserver &observer);
    // Bind application declarations to an unchanged, complete producer history.
    // Document byte spans must name raw sources already retained by the producer.
    static NativeRecord publish_annotations(const NativeRecord &producer,
        NativeRecordValue study, NativeRecordValue documents);
    NativeRecord producer_publication() const;
    ~ScientificStudy();
    ScientificStudy(const ScientificStudy &) = delete;
    ScientificStudy &operator=(const ScientificStudy &) = delete;
    Json catalog() const;
    // A compact inventory plus complete, independently addressed source pages.
    // Offsets select presentation only; mathematics retains the original owner.
    Json source_index() const;
    Json read_source(const Json &request) const;
    Json read_source(const Json &request, rfc::runtime::ExecutionLease lease) const;
    Json declaration() const;
    std::string digest() const;
    std::shared_ptr<const rfc::storage::TemporalEndpoint> endpoint() const;
    const NativeRecord &native_publication() const;
    void write_native(const fs::path &path) const;
    // Reuse the admitted immutable field pin, including its actual publication.
    rfc::storage::PinnedNativeAccession field_accession(const std::string &id) const;
    // An atlas member reuses this study's actual directory and publication.
    rfc::storage::NativeAtlasSource native_atlas_source(std::uint64_t key) const;
    // Compatibility labels resolve to the already admitted primary identity;
    // gluing never reconstructs its source from presentation catalogues.
    rfc::storage::StableCellKey cell_identity(rfc::core::Grade grade, const std::string &id) const;
    // Read the declared name at an actual native grade coordinate.
    std::string cell_identifier(rfc::core::Grade grade, std::size_t coordinate) const;
    Json field_declaration(const std::string &id) const;
    Json source_binding() const;
    // Owns only the actual base and frames through the admitted cutoff.
    std::shared_ptr<const rfc::storage::NativeSeriesPrefix> native_history() const;
    // Resolve a declared identity and quoted evidence in this immutable owner.
    // Annotation text does not define a chain map or alter the native source.
    Json source_target(const Json &target, const Json &evidence) const;
};

class ScientificSession {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit ScientificSession(std::shared_ptr<const ScientificStudy> study,
                               std::shared_ptr<CpuResources> resources = {},
                               ComputeOptions compute = {});
    ~ScientificSession();
    ScientificSession(const ScientificSession &) = delete;
    ScientificSession &operator=(const ScientificSession &) = delete;
    Json execute(const Json &query);
    NativeRecord execute_record(const NativeRecord &query);
    Json read_source(const Json &request) const;
    // Explicit current host limits and cancellation for declared native responses.
    // Other query families retain their own execution contracts.
    Json execute_native_response(const Json &query, rfc::runtime::ExecutionLease lease);
    Json execute_native_hodge(const Json &query, rfc::runtime::ExecutionLease lease);
    Json execute_native_observation(const Json &query, rfc::runtime::ExecutionLease lease);
    Json execute_native_novelty(const Json &query, rfc::runtime::ExecutionLease lease);
    Json execute_tensor_window(const Json &query, rfc::runtime::ExecutionLease lease);
    Json execute_green_response(const Json &query, rfc::runtime::ExecutionLease lease);
    void validate(const Json &query);
    Json execute_plan(const Json &plan);
    NativeRecord execute_plan_record(const NativeRecord &plan);
    Json propose_annotations(const Json &annotations);
    Json native_cache_statistics() const;
    // Counts source-bound preparations separately from immutable certificates.
    Json temporal_preparation_statistics() const;
    Json native_response_preparation_statistics() const;
    Json native_hodge_preparation_statistics() const;
    Json native_observation_preparation_statistics() const;
    Json native_novelty_preparation_statistics() const;
    void clear_native_novelty_preparations();
    Json tensor_window_preparation_statistics() const;
    void clear_tensor_window_preparations();
    std::size_t preparation_count() const;
};
// Runtime-independent admission. Models may propose queries, bounded plans,
// evidence-anchored annotations or a clarification; none mutate the source.
Json apply_scientific_proposal(ScientificSession &session, const Json &proposal);
} // namespace gadget
