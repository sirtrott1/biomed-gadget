#pragma once
#include <gadget/source.hpp>
#include <rfc/storage/native_accession.hpp>

namespace gadget {
// Immutable native views own the publication and survive this adapter's lifetime.
class NativeSnapshot {
    std::shared_ptr<const rfc::storage::NativeAccessionDirectory> directory_;
    std::map<std::string, std::uint64_t> fields_;

  public:
    NativeSnapshot(const rfc::storage::TemporalEndpoint &endpoint, std::uint64_t source_namespace,
                   std::map<std::string, std::uint64_t> fields);
    NativeSnapshot(const NativeRecord &publication,std::map<std::string,std::uint64_t> fields);
    NativeSnapshot(const NativeSnapshot &) = delete;
    NativeSnapshot &operator=(const NativeSnapshot &) = delete;
    std::shared_ptr<const rfc::storage::TemporalEndpoint> endpoint() const;
    rfc::storage::PinnedNativeAccession field(const std::string &name) const;
    rfc::storage::PinnedNativeAccession source_grade(rfc::core::Grade grade) const;
};

// Admit complete source bytes and attach descriptions to the actual native cells.
NativeRecord capture_workspace(const Workspace &source);
void admit_workspace(Workspace &source);

// Version 1 captures the admitted PAGER/GET study, including imported measurements.
// These functions use canonical codecs; publication is a local staged export.
Json save_snapshot(const Workspace &source, const fs::path &destination,
                   bool export_json_views = false);
Workspace load_snapshot(const fs::path &directory);
} // namespace gadget
