#pragma once
#include <cstddef>
#include <functional>
#include <map>
#include <span>
#include <stop_token>
#include <string_view>
#include <vector>

namespace gadget::examples {
// This profile reads immutable, decompressed source bytes. It neither creates
// native fields nor infers normalization, exposure, sample selection or design.
enum class GtexLayout { Gct12, PortalTissue13 };
enum class GtexFeature { Gene, Transcript };
enum class GtexMeasurement { IntegerCounts, EstimatedCounts, Tpm };
struct GtexLimits {
    std::size_t source_bytes = 1024u * 1024u * 1024u;
    std::size_t rows = 1000000;
    std::size_t columns = 100000;
    std::size_t record_bytes = 32u * 1024u * 1024u;
    std::size_t field_bytes = 1024u * 1024u;
    std::size_t cells = 100000000;
};
struct GtexDeclaration {
    GtexLayout layout = GtexLayout::Gct12;
    GtexFeature feature = GtexFeature::Gene;
    GtexMeasurement measurement = GtexMeasurement::IntegerCounts;
};
struct GtexSourceRecord {
    std::size_t line{}, byte_offset{}, byte_length{};
    std::string_view original_bytes; // Includes the original line terminator, if any.
};
struct GtexHeader {
    GtexDeclaration declaration;
    std::size_t features{}, declared_columns{};
    std::vector<std::string_view> samples;
    GtexSourceRecord version_record, dimensions_record, labels_record;
};
struct GtexValue {
    std::string_view token;
    bool observed{}; // Only the exact source token NA denotes missingness.
};
struct GtexRowView {
    const GtexHeader &header;
    std::size_t ordinal{};
    GtexSourceRecord source;
    std::string_view feature, description;
    std::span<const GtexValue> values;
};
struct GtexSummary {
    std::size_t rows{}, samples{}, observed{}, missing{}, observed_zeros{};
};
// All views refer to the caller's unchanged source bytes. Row/value spans and
// the header reference are callback only. Callbacks are provisional until the
// whole visit returns successfully (late invalid rows or cancellation can fail).
GtexSummary visit_gtex_gct(std::string_view source, GtexDeclaration declaration,
                           const GtexLimits &limits,
                           const std::function<void(const GtexRowView &)> &visitor,
                           std::stop_token stop = {});

enum class GtexMetadataKind { SampleAttributes, SubjectPhenotypes };
struct GtexMetadataRow {
    GtexSourceRecord source;
    std::vector<std::string_view> fields;
};
struct GtexMetadataTable {
    GtexMetadataKind kind;
    GtexSourceRecord header_record;
    std::vector<std::string_view> columns;
    std::map<std::string_view, GtexMetadataRow, std::less<>> rows;
};
// All original fields, including empty cells, NA and categorical age ranges,
// remain views of their original source; no metadata is silently dropped.
GtexMetadataTable read_gtex_metadata(std::string_view source, GtexMetadataKind kind,
                                     const GtexLimits &limits, std::stop_token stop = {});
struct GtexSampleView {
    std::string_view sample, donor, tissue;
    const GtexMetadataRow *sample_attributes{}, *subject_phenotypes{};
    bool repeated_donor{};
};
struct GtexSampleSet {
    std::vector<GtexSampleView> samples;
    std::size_t distinct_donors{}, repeated_donors{};
    // A count model must explicitly check its experimental units. This guard
    // refuses repeated donors; admission itself preserves all original samples.
    void require_independent_donors() const;
};
// Returned views require both source strings and tables to remain alive and
// unchanged. No library total is inferred from metadata or from matrix values.
GtexSampleSet join_gtex_samples(std::span<const std::string_view> samples,
                                const GtexMetadataTable &sample_attributes,
                                const GtexMetadataTable &subject_phenotypes,
                                std::stop_token stop = {});
} // namespace gadget::examples
