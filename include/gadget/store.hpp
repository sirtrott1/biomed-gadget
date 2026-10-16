#pragma once
#include <gadget/snapshot.hpp>
#if GADGET_HAS_DURABLE_STORE
#include <rfc/storage/concurrent_native_history.hpp>
#include <rfc/storage/native_durable_store.hpp>

namespace gadget {
enum class StudyPreparationStage {
    DirectoryCreated,
    FileWritten,
    FileSynced,
    DirectorySynced,
    Renamed,
    ParentSynced
};
// Fault injection seam for the application provenance preparation, before HEAD.
using StudyPreparationHook = std::function<void(StudyPreparationStage)>;
// Exclusive Linux local profile. Every operation opens the canonical coordinator.
// Owning readers returned by load_store survive its closure and later publication.
// Preparation preserves the required native archive and present, validated JSON
// views byte for byte. A prepared transaction requires the same artifacts on retry.
// Snapshot decoding and archive digest use the same opened byte owner.
Json initialize_store(const fs::path &snapshot, const fs::path &destination,
                      rfc::storage::DurableFaultHook hook = {},
                      StudyPreparationHook preparation_hook = {});
Json publish_store(const fs::path &store, const fs::path &snapshot, const Json &request,
                   rfc::storage::DurableFaultHook hook = {},
                   StudyPreparationHook preparation_hook = {});
Json store_status(const fs::path &store, std::optional<std::uint64_t> transaction = {});
Workspace load_store(const fs::path &store);
std::uint64_t store_transaction(const std::string &decimal);
class ConcurrentStudyStore;
class StudySnapshot {
    struct State;
    std::shared_ptr<const State> state_;
    explicit StudySnapshot(std::shared_ptr<const State> state);
    friend class ConcurrentStudyStore;

  public:
    rfc::storage::DurableNativeRoot root() const;
    std::shared_ptr<const Workspace> workspace() const;
};
// One service process: owning readers, serialized writers, whole root conflicts.
// Active calls require service lifetime; StudySnapshot pins outlive it.
class ConcurrentStudyStore {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit ConcurrentStudyStore(const fs::path &directory);
    ~ConcurrentStudyStore();
    ConcurrentStudyStore(const ConcurrentStudyStore &) = delete;
    ConcurrentStudyStore &operator=(const ConcurrentStudyStore &) = delete;
    StudySnapshot snapshot() const;
    Json status(std::optional<std::uint64_t> transaction = {}) const;
    Json publish(const StudySnapshot &expected, const fs::path &snapshot, const Json &request,
                 std::stop_token stop = {}, rfc::storage::DurableFaultHook hook = {},
                 StudyPreparationHook preparation_hook = {});
    void recover();
};
} // namespace gadget
#endif
