#pragma once
#include "assay_detail.hpp"

namespace gadget::assay_detail {
struct ResourceImage {
    Json declaration;
    std::string bytes, reference_bytes;
};
using ResourceImages = std::map<std::string, ResourceImage>;
struct Alignment {
    std::string query, group, contig, strand, decoded;
    std::uint16_t flags{};
    std::uint8_t mapq{};
    std::uint64_t start{}, end{}, query_length{};
    std::vector<std::pair<std::uint64_t, std::uint64_t>> blocks, junctions;
};
struct Alignments {
    std::string header;
    std::vector<Alignment> records;
};
Alignments decode_alignments(const ResourceImage &resource, bool render_records = true);
Json inspect_hts(const ResourceImage &resource, bool summary = false);
std::string hts_codec_version();
struct AssayFieldsSource {
    rfc::storage::TemporalEndpoint endpoint;
    NativeRecordValue declaration;
    std::vector<RawByteSource> sources;
};
AssayFieldsSource interval_fields_source(const ResourceImage &resource, const Json &request);
Json interval_fields_view(const NativeRecord &record);
AssayFieldsSource annotation_fields_source(const ResourceImage &resource, const Json &request);
Json annotation_fields_view(const NativeRecord &record);
AssayFieldsSource motif_fields_source(const ResourceImage &resource, const Json &request);
Json motif_fields_view(const NativeRecord &record);
AssayFieldsSource alignment_fields_source(const ResourceImage &resource, const Json &request);
AssayFieldsSource feature_count_fields_source(const ResourceImages &images, const Json &request,
                                              std::shared_ptr<CpuResources> resources);
Json allocation_analysis(const Json &request, std::shared_ptr<CpuResources> resources,
                         std::optional<AssayFieldsSource> *retained = nullptr,
                         const std::vector<Vector> *native_counts = nullptr,
                         bool scale_samples = false);
Json alignment_fields_view(const NativeRecord &record);
Json feature_count_fields_view(const NativeRecord &record);
Json expression_fields_view(const NativeRecord &record);
Json aggregation_fields_view(const NativeRecord &record);
Json allocation_fields_summary_view(const NativeRecord &allocation_fields);
Json genotype_analysis(const ResourceImage &resource, const Json &request,
                       std::optional<AssayFieldsSource> *retained = nullptr);
Json genotype_fields_view(const NativeRecord &record);
Json genomic_analysis(const ResourceImages &resources, const Json &request,
                      std::shared_ptr<CpuResources> budget,
                      std::optional<AssayFieldsSource> *retained = nullptr);
Json chromatin_fields_view(const NativeRecord &record);
Json feature_allocation_evidence(const ResourceImages &images, const Json &request,
                                 Vector *native_counts = nullptr);
Json expression_allocations(const ResourceImages &images, const Json &request,
                            std::shared_ptr<CpuResources> resources,
                            std::optional<AssayFieldsSource> *retained = nullptr,
                            const NativeRecord *previous_coupling = nullptr);
Json constraint_request_view(const NativeRecord &record);
Json controlled_coverage(const ResourceImages &images, const Json &request);
using SpatialRelationSink = std::function<void(std::size_t, std::size_t, const Rational &)>;
Json spatial_analysis(const Json &assay, const Json &request, const SpatialRelationSink &sink = {});
AssayFieldsSource spatial_fields_source(const ScientificStudy &source, const Json &assay,
                                        const Json &request);
Json quality_analysis(const ResourceImage &resource, const Json &request);
Json aggregate_analysis(const Matrix &matrix, const Json &declaration, const Json &request,
                        const rfc::runtime::CpuResourceCoordinator::Reservation &reservation);
// Disjoint original-coordinate storage ranges. These correspondences address
// quantitative measurements; they do not glue independent source complexes.
struct NativeAssayPartition {
    std::size_t begin{}, end{};
    std::shared_ptr<const ScientificStudy> owner;
};
// Native field pins supply source coordinates to the canonical map actions.
Json native_aggregate_analysis(const Axis &features, const Axis &observations,
                               const std::vector<NativeAssayPartition> &partitions,
                               const Json &declaration, const Json &request,
                               const rfc::runtime::CpuResourceCoordinator::Reservation &reservation,
                               std::optional<AssayFieldsSource> *retained = nullptr);
Json design_aggregation(const Matrix &matrix, const Json &declaration, const Json &request,
                        const rfc::runtime::CpuResourceCoordinator::Reservation &reservation,
                        const std::vector<NativeAssayPartition> &partitions = {},
                        std::optional<AssayFieldsSource> *retained = nullptr);
Json alignment_analysis(const ResourceImages &images, const Json &request,
                        std::optional<AssayFieldsSource> *retained = nullptr);
Json sequence_alignment_fields_view(const NativeRecord &record);
struct SequenceRecord {
    std::string id, bases;
    Json evidence;
};
// Views address original bytes; only a requested record or block is normalized.
class SequenceSource {
    struct Segment {
        std::size_t source, begin, length;
    };
    const ResourceImage &resource_;
    std::vector<SequenceSpan> spans_;
    std::vector<std::vector<Segment>> segments_;
    std::string digest_;

  public:
    explicit SequenceSource(const ResourceImage &resource);
    SequenceSource(ResourceImage &&) = delete;
    std::size_t size() const { return spans_.size(); }
    std::size_t length(std::size_t record) const { return spans_.at(record).length; }
    const std::string &digest() const { return digest_; }
    const SequenceSpan &original_span(std::size_t record) const { return spans_.at(record); }
    std::string id(std::size_t record) const;
    Json evidence(std::size_t record) const;
    std::string bases(std::size_t record, std::size_t begin, std::size_t length) const;
    std::vector<std::pair<std::size_t, std::size_t>> byte_spans(
        std::size_t record, std::size_t begin, std::size_t length) const;
    SequenceRecord record(std::size_t index) const;
    void blocks(std::size_t record, std::size_t width,
                const std::function<void(std::size_t, std::string_view)> &visit) const;
};
std::vector<SequenceRecord> sequence_records(const ResourceImage &resource);
NativeRecordValue read_protocol_fields(ResourceImages images, const Json &request);
Json reference_matches(const ResourceImages &images, const Json &request);
Json assembly_graph(const ResourceImage &resource, const Json &request,
                    std::optional<AssayFieldsSource> *retained = nullptr);
Json sequence_assembly_fields_view(const NativeRecord &record);
Json indexed_alignment(const ResourceImages &images, const Json &request);
Json splice_search(const ResourceImages &images, const Json &request);
Json read_reconstruction(const ResourceImage &resource, const Json &request);
AssayFieldsSource indexed_alignment_fields_source(const ResourceImages &images, const Json &request);
AssayFieldsSource reference_match_fields_source(const ResourceImages &images, const Json &request);
AssayFieldsSource splice_search_fields_source(const ResourceImages &images, const Json &request);
AssayFieldsSource read_reconstruction_fields_source(const ResourceImage &resource, const Json &request);
Json indexed_alignment_fields_view(const NativeRecord &record);
Json reference_match_fields_view(const NativeRecord &record);
Json splice_search_fields_view(const NativeRecord &record);
Json read_reconstruction_fields_view(const NativeRecord &record);
Json query_batch(const AssayWorkspace &workspace, const Json &request,
                 const std::shared_ptr<CpuResources> &resources);
NativeRecord query_batch_record(const AssayWorkspace &workspace, const NativeRecord &request,
                               const std::shared_ptr<CpuResources> &resources);
Json query_batch_fields_view(const NativeRecord &record);
Json constrained_allocations(const Json &request,
                             std::optional<AssayFieldsSource> *retained = nullptr);
Json constraint_fields_view(const NativeRecord &record);
Json relational_persistence(const ScientificStudy &study, const Json &request);
} // namespace gadget::assay_detail
