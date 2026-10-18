#pragma once
#include <cstddef>
#include <optional>
namespace gadget {
// Bound additional address space before admitting a CLI operation. The default
// reserves half of the currently available memory for other processes.
std::optional<std::size_t> limit_process_memory(std::optional<std::size_t> bytes = {});
}
