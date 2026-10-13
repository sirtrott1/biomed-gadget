#pragma once
#include <gadget/io.hpp>
#include <gadget/study.hpp>

namespace gadget {
// Exact structural witnesses on the caller's retained source and metrics.
// These functions render results produced from native fields and certificates.
Json analyze_chain(const NativeStudy &study, const Json &query, bool validate_only = false);
Json compare_realized_paths(const NativeStudy &study, const Json &query,
                            bool validate_only = false);
} // namespace gadget
