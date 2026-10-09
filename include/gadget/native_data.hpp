#pragma once
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <rfc/storage/native_atlas_types.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <rfc/storage/declaration_codec.hpp>
#include <rfc/storage/prepared_source_program.hpp>
#include <rfc/field/augmented_green.hpp>
#include <rfc/field/execution_tower.hpp>
#include <rfc/statistics/field_moments.hpp>

namespace gadget {
enum class NativeDataFormat { Unknown, RelationalComplex, TemporalSeries };
NativeDataFormat detect_native_data(std::span<const std::byte> bytes) noexcept;
// Host addressability, not a scientific population or dimensionality ceiling.
rfc::storage::NativeSeriesLimits native_data_limits();
rfc::storage::NativeAccessionLimits native_accession_limits();

// Own the unchanged input bytes. Identity is supplied by the source declaration,
// rather than inferred from a filename, hash, display label or equal contents.
class RawByteSource {
    std::uint64_t key_;
    std::shared_ptr<const std::vector<std::byte>> bytes_;
    bool frozen_{true};
    // A public vector transfer can have prior mutable aliases. Admission
    // freezes it; only fresh reader/compiler output may adopt without copying.
    RawByteSource freeze() const;
    static RawByteSource adopt_fresh(std::uint64_t, std::vector<std::byte>);
    friend class NativeStudy;
    friend class NativeRecord;
  public:
    RawByteSource(std::uint64_t key, std::span<const std::byte> bytes);
    static RawByteSource own(std::uint64_t key, std::vector<std::byte> bytes);
    static RawByteSource read(std::uint64_t key, const std::filesystem::path &path);
    // Independent declared identity sharing the same immutable byte owner.
    RawByteSource rekey(std::uint64_t key) const;
    bool shares_owner(const RawByteSource &other) const noexcept { return bytes_ == other.bytes_; }
    std::uint64_t key() const noexcept { return key_; }
    std::span<const std::byte> bytes() const noexcept { return *bytes_; }
};
// The caller's native Product declaration supplies these stable factor names.
// Each factor is an exact Integer: source identity, byte offset and byte length.
struct RawSpanDeclaration {
    rfc::composite::DeclarationId type;
    std::uint64_t source_factor{}, offset_factor{}, length_factor{};
};
class PinnedRawSpan {
    RawByteSource source_;
    rfc::storage::PinnedNativeAccession declaration_;
    std::size_t offset_, length_;
    friend class NativeStudy;
    PinnedRawSpan(RawByteSource source, rfc::storage::PinnedNativeAccession declaration,
                  std::size_t offset, std::size_t length);
  public:
    std::span<const std::byte> bytes() const noexcept;
    std::uint64_t source_key() const noexcept { return source_.key(); }
    rfc::storage::PinnedNativeAccession declaration() const { return declaration_; }
};
struct DeclaredByteValue {
    PinnedRawSpan source;
    rfc::composite::CompiledLayout layout;
    rfc::composite::DeclaredValue value;
};

// The canonical source is the actual relational complex, native tensor fields,
// declaration DAG and retained operation programs. There is no JSON catalogue.
class NativeStudy {
    struct State;
    std::shared_ptr<const State> state_;
    explicit NativeStudy(std::shared_ptr<const State> state);
  public:
    NativeStudy(std::shared_ptr<const rfc::core::RelationalComplex> source,
                std::vector<rfc::field::ExactMetric> metrics);
    NativeStudy(rfc::core::RelationalComplex source,
                std::vector<rfc::field::ExactMetric> metrics);
    NativeStudy(const rfc::storage::TemporalEndpoint &endpoint, std::uint64_t source_namespace,
                std::vector<RawByteSource> sources = {},
                rfc::storage::NativeAccessionLimits limits = native_accession_limits());
    explicit NativeStudy(const rfc::storage::NativeSeriesEndpoint &endpoint,
                std::vector<RawByteSource> sources = {},
                rfc::storage::NativeAccessionLimits limits = native_accession_limits());
    // Retain a canonical publication; prepare SDK accession indexes on demand.
    static NativeStudy from_publication(const rfc::storage::NativeSeriesEndpoint &endpoint,
                                       std::vector<RawByteSource> sources = {});
    // Inspect canonical byte signatures. RFCX needs an explicit namespace;
    // RFCXS owns its namespace and an optional actual history cutoff.
    static NativeStudy read(const std::filesystem::path &path,
                            std::vector<RawByteSource> sources = {},
                            std::optional<std::uint64_t> source_namespace = {},
                            std::optional<std::uint64_t> cutoff = {},
                            rfc::storage::NativeSeriesLimits limits = native_data_limits());
    std::shared_ptr<const rfc::storage::TemporalEndpoint> endpoint() const;
    std::shared_ptr<const rfc::storage::NativeAccessionDirectory> accession_directory() const;
    const rfc::core::RelationalComplex &source() const;
    const rfc::field::ExactExecutionTower &tower() const;
    rfc::field::ExactHodgeOperator hodge(rfc::core::Grade grade) const;
    rfc::numeric::RationalVector execute(const rfc::field::OperatorRequest &request,
                                       const rfc::numeric::RationalVector &input) const;
    rfc::numeric::Rational moment(rfc::core::Grade grade, std::size_t order,
        const rfc::numeric::RationalVector &x, const rfc::numeric::RationalVector &y) const;
    rfc::field::ExactGreenResult green_action(rfc::core::Grade grade,
        const rfc::numeric::Rational &lambda, const rfc::numeric::RationalVector &input) const;
    rfc::storage::NativeAtlasSource atlas_source(std::uint64_t key) const;
    rfc::storage::PinnedNativeAccession access(const rfc::storage::NativeAccessionRequest &request) const;
    PinnedRawSpan raw_span(std::uint64_t value, const RawSpanDeclaration &declaration) const;
    // Decode an authored composite layout directly from its original byte span.
    // The returned exact value retains both the layout and provenance owners.
    DeclaredByteValue decode_raw_span(std::uint64_t value,
        const RawSpanDeclaration &span, rfc::composite::DeclarationId type,
        std::size_t bit_count) const;
    // Executes the authored instruction word over this publication's original
    // packed source aspects, using the canonical prepared program/workspace.
    rfc::storage::PreparedNativeSourceProgram program(
        std::uint64_t key, std::uint64_t layout,
        rfc::storage::StableOccurrenceKey occurrence) const;
};
} // namespace gadget
