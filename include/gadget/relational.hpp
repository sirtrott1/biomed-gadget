#pragma once
#include <gadget/io.hpp>

namespace gadget {
// Bounded retained complex research profile. Declarations own source meaning;
// canonical RexGraph owns topology, metric fields, Green solves and scheduling.
void run_relational_audit(const fs::path &input, const fs::path &output, std::size_t workers);
} // namespace gadget
