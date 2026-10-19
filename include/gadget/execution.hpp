#pragma once
#include <gadget/statistics/execution.hpp>

// Compatibility names for the contextual numerical analysis interface.
namespace gadget {
using statistics::Execution;
using statistics::QuerySession;
using statistics::analyze;
using statistics::approximate_values;
using statistics::binary64_values;
using statistics::cross_relation_gram;
using statistics::execute_query;
using statistics::norm;
using statistics::union_field;
} // namespace gadget
