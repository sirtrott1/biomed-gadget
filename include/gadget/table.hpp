#pragma once
#include <gadget/scientific.hpp>

namespace gadget {
// Strict declared unquoted table admission. The importer preserves rows, fields
// and original artifact bytes and builds the source/archive with canonical APIs.
// evidence_profile=document_spans_v2 avoids repeated quotations; omitted retains
// quoted_v1. Optional intake budgets count source/record bytes and token offsets.
// A private retained image is mapped through RexGraph; native admission stays complete.
void import_table(const fs::path &specification, const fs::path &destination);
} // namespace gadget
