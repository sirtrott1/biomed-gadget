#pragma once
#include <gadget/native_record.hpp>

namespace gadget {
// Bind an RFCXS source to its declared axes, units and original record addresses.
// Optional value and support tensors must match the retained source digests and
// logical domain. Invalid declarations throw before a NativeRecord is returned.
[[nodiscard]] NativeRecord import_data_field(const RawByteSource &source,
                                             const NativeRecord &coordinates,
                                             std::uint64_t source_namespace,
                                             const RawByteSource *values = nullptr,
                                             const RawByteSource *support = nullptr);
} // namespace gadget
