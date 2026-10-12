#pragma once
#include "scientific_source_snapshot.hpp"

namespace gadget::research_detail {
using AblationSourcePublication = scientific_detail::SourceSnapshotPublication;
inline AblationSourcePublication ablation_source_publication(
    const NativeStudy &study, const NativeRecord &original,
    const std::vector<std::vector<std::uint32_t>> &to_original) {
    return scientific_detail::source_snapshot_publication(study, original, to_original);
}
} // namespace gadget::research_detail
