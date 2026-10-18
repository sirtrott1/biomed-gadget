#pragma once
#include <gadget/native_record.hpp>

namespace gadget {
// The native source owns coefficients. The declaration contains only source
// identities, axis coordinates, units and references to immutable records.
[[nodiscard]] NativeRecord import_data_field(const RawByteSource &source,
                                             const NativeRecord &coordinates,
                                             std::uint64_t source_namespace,
                                             const RawByteSource *values = nullptr,
                                             const RawByteSource *support = nullptr);
} // namespace gadget
