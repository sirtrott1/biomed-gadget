#pragma once
#include <gadget/scientific.hpp>

namespace rfc::database { class InterpretedSnapshot; }

namespace gadget {
// HTTPS retrieval is optional. Admission and queries work entirely offline.
Json plan_database_fetch(const Json &request);
Json fetch_databases(const Json &request, const fs::path &destination);
Json database_capabilities();
// Retain declared original files beside an unchanged admitted study. The
// declaration records provenance; measurement correspondence is the adapter's
// separately checked contract, rather than an inference from matching hashes.
NativeRecord retain_study_sources(const ScientificStudy &study, const NativeRecord &request);
NativeRecord capture_database_snapshot(const fs::path &directory);
void import_evidence(const fs::path &snapshot, const fs::path &destination,
                     std::uint32_t adapter_version = 2);
// Source adapters validate their original format and select through the same
// canonical native interpretation. Applications retain the adapter with its store.
class EvidenceAdapter {
  public:
    virtual ~EvidenceAdapter() = default;
    virtual void verify(const NativeRecord &, const rfc::database::InterpretedSnapshot &) const = 0;
    virtual void validate(const Json &, const std::string &digest) const = 0;
    virtual NativeRecord execute(const NativeRecord &, const rfc::database::InterpretedSnapshot &,
                                 const NativeRecord &) const = 0;
    virtual NativeRecord context(const NativeRecord &, const rfc::database::InterpretedSnapshot &,
                                 const NativeRecord &) const = 0;
};
class EvidenceStore {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit EvidenceStore(const fs::path &directory);
    EvidenceStore(const NativeRecord &publication, std::shared_ptr<const EvidenceAdapter> adapter = {});
    ~EvidenceStore();
    const rfc::database::InterpretedSnapshot &interpreted_source() const;
    std::string digest() const;
    const NativeRecord &native_publication() const;
    // The request selects original assertions and a question on this frozen
    // snapshot. Native context retains those assertions and their byte owners.
    NativeRecord source_context(const NativeRecord &request) const;
    // Check the admitted reported point projection against original assertions.
    // The receipt retains separate study/evidence owners and a native selector.
    NativeRecord bind_reported_field(const ScientificStudy &study,
                                    const NativeRecord &request) const;
    // Admit a prospective query or source annotation without executing it.
    // Snapshot retrieval does not establish historical or biological time.
    NativeRecord admit_proposal(const NativeRecord &context_request,
                               const NativeRecord &proposal,
                               const std::string &question) const;
    void validate(const Json &query) const;
    NativeRecord execute_record(const NativeRecord &query) const;
    Json execute(const Json &query) const;
};
// Non numerical model observations retain input identities and remain prospective.
Json admit_observation(const Json &proposal, const std::string &source_digest, const Json &media);
} // namespace gadget
