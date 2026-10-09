#pragma once
#include "genomic_detail.hpp"
#include <gadget/native_record.hpp>

namespace gadget::assay_detail {
inline constexpr std::string_view coverage_fields_schema = "gadget.coverage.fields.v1";
struct CoverageSource {
    rfc::storage::TemporalEndpoint endpoint;
    NativeRecordValue declaration;
    std::vector<RawByteSource> sources;
};
CoverageSource controlled_coverage_source(const ResourceImages &, const Json &);
NativeRecord controlled_coverage_record(const ResourceImages &, const Json &);
Json coverage_fields_view(const NativeRecord &);
} // namespace gadget::assay_detail
