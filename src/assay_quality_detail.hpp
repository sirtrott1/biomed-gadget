#pragma once
#include <gadget/record_view.hpp>

namespace gadget::assay_detail {
struct ResourceImage;
// The numerical result is an authored native declaration, not an imported view.
NativeRecordValue quality_fields(const ResourceImage &, const Json &);
Json quality_fields_view(const NativeRecord &);
inline constexpr std::string_view quality_fields_schema = "gadget.read.quality.fields.v1";
} // namespace gadget::assay_detail
