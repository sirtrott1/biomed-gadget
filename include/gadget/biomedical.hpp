#pragma once
#include <gadget/scientific.hpp>

namespace gadget {
// HTTPS retrieval is optional. Admission and queries work entirely offline.
Json plan_database_fetch(const Json &request);
Json fetch_databases(const Json &request, const fs::path &destination);
Json database_capabilities();
NativeRecord capture_database_snapshot(const fs::path &directory);
void import_evidence(const fs::path &snapshot, const fs::path &destination,
                     std::uint32_t adapter_version = 2);
class EvidenceStore {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit EvidenceStore(const fs::path &directory);
    ~EvidenceStore();
    std::string digest() const;
    const NativeRecord &native_publication() const;
    void validate(const Json &query) const;
    Json execute(const Json &query) const;
};
// Non-numerical model observations retain input identities and remain prospective.
Json admit_observation(const Json &proposal, const std::string &source_digest, const Json &media);
} // namespace gadget
