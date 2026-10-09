#pragma once
#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/native_execution_binding.hpp>
#include <rfc/storage/native_series.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <set>
#include <tuple>

namespace gadget::scientific_detail {
void require(bool valid, const std::string &message);
void source_array(const Json &value, const std::string &name);
void members(const Json &object, const std::set<std::string> &required,
             const std::set<std::string> &optional = {});
std::size_t count(const Json &value, std::size_t limit, const std::string &name);
Json cardinality(std::size_t value);
Json exact_count(const rfc::numeric::Integer &value);
rfc::numeric::SparseRationalTensor declared_comparison(const Json &map, std::size_t rows,
                                                       std::size_t columns);
std::string local_artifact(const fs::path &root, const Json &record,
                           std::size_t limit = source_byte_limit());
rfc::storage::NativeSeriesLimits source_codec_limits(std::size_t bytes);
std::uint64_t key(const Json &value);
std::string text(const Json &value, std::size_t limit = source_byte_limit());
Json seal(Json result, const std::string &name);
rfc::runtime::ExecutionLease source_execution_lease(const std::string &name);
struct Document {
    RetainedRecord metadata;
    std::string_view bytes;
};
Json source_evidence(const Json &value, const std::map<std::string, Document> &documents,
                     bool require_utf8 = false, bool allow_references = false);
fs::path local_artifact_path(const fs::path &root, const Json &record);
struct Context {
    RetainedRecord declaration;
    std::optional<std::uint64_t> block;
};
struct Field {
    RetainedRecord declaration;
    std::uint64_t tensor;
};
struct FieldSet {
    RetainedRecord declaration;
    rfc::storage::NativeFieldRealization realization;
};
} // namespace gadget::scientific_detail

namespace gadget {
struct ScientificStudy::State {
    fs::path root;
    RetainedRecord declaration;
    std::optional<NativeRecord> native_publication;
    std::string digest, native_digest, declaration_digest;
    std::map<std::string, std::shared_ptr<const ScientificStudy>> publications;
    std::map<std::string, RetainedRecord> publication_declarations;
    std::shared_ptr<const rfc::storage::NativeAccessionDirectory> directory;
    std::shared_ptr<const rfc::storage::NativeSeriesPrefix> history;
    std::map<std::uint64_t, rfc::storage::PinnedNativeAccession> field_pins;
    std::map<std::string, rfc::storage::PinnedNativeAccession> operator_pins;
    std::map<std::string, RetainedRecord> operators;
    std::vector<rfc::storage::PinnedNativeAccession> grade_pins;
    std::map<std::string, scientific_detail::Document> documents;
    std::map<rfc::storage::StableCellKey, RetainedRecord> cells;
    std::vector<std::vector<rfc::storage::StableCellKey>> coordinate_cells;
    std::map<std::string, rfc::storage::StableCellKey> cell_ids;
    std::map<rfc::storage::StableOccurrenceKey, RetainedRecord> occurrences;
    std::map<std::string, scientific_detail::Context> contexts;
    std::map<std::string, scientific_detail::Field> fields;
    std::map<std::string, scientific_detail::FieldSet> field_sets;
    // Presentation indices only; numerical data and declarations retain their
    // existing owners. Late pages do not rescan all earlier catalogue entries.
    std::map<std::string, std::vector<std::string>> catalog_keys;
    std::vector<rfc::storage::StableCellKey> cell_keys;
    std::vector<rfc::storage::StableOccurrenceKey> occurrence_keys;
    mutable std::mutex evidence_reference_mutex;
    mutable std::map<std::tuple<std::string, std::size_t, std::size_t>, std::string>
        evidence_reference_hashes;
    explicit State(const fs::path &path, std::optional<NativeRecord> publication = {},
                   std::optional<rfc::storage::NativeTemporalSeries> producer = {});
    Json evidence(const Json &value) const;
    Json catalog() const;
    Json history_metadata() const;
    Json referenced(const Json &metadata) const;
    Json field_metadata(const scientific_detail::Field &field) const;
    Json operator_metadata(const std::string &id) const;
    Json carrier_metadata(const rfc::storage::TemporalCarrier &carrier) const;
};
} // namespace gadget
