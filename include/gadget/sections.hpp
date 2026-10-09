#pragma once
#include <gadget/io.hpp>
#include <gadget/study.hpp>
namespace gadget {
// Declarations identify every local carrier and retained incidence map. RexGraph
// owns restriction assembly, cohomology, projection and relation tensors.
Json analyze_local_sections(const NativeStudy &study, const Json &query, bool validate_only = false,
                            bool include_cochain_maps = false);
} // namespace gadget
