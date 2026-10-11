#pragma once
#include <gadget/native_data.hpp>
#include <rfc/storage/native_encoded_tensor.hpp>
#include <functional>
#include <variant>

namespace gadget {
// Application data is compiled into authored native declaration constructors.
// Text addresses byte owners; source bytes are retained without text encoding.
struct NativeRecordValue {
    struct Object { std::vector<std::pair<std::string, NativeRecordValue>> members; };
    struct Array { std::vector<NativeRecordValue> elements; };
    // Emit exactly the declared number of complete values in order. Each value
    // is encoded with its native declaration before the next value is requested.
    struct PackedArray {
        std::size_t count{};
        std::function<void(const std::function<void(NativeRecordValue)> &)> emit;
    };
    struct ByteSpan { std::uint64_t source; std::size_t offset, length; bool text{true}; };
    struct Decimal { rfc::numeric::Rational exact; std::string literal; };
    struct ExactNumber { rfc::numeric::Rational exact; bool approximate_view; };
    struct Owned {
        std::uint64_t block{};
        std::optional<rfc::storage::StableCellKey> cell;
        std::optional<rfc::storage::StableOccurrenceKey> occurrence;
        std::shared_ptr<const NativeRecordValue> payload;
    };
    using Data = std::variant<std::monostate, bool, rfc::numeric::Integer,
        rfc::numeric::Rational, std::string, Object, Array, ByteSpan, Decimal,
        rfc::numeric::SparseRationalTensor, ExactNumber, Owned, PackedArray>;
    Data data;
    NativeRecordValue() = default;
    template<class T> explicit NativeRecordValue(T value) : data(std::move(value)) {}
};
enum class NativeRecordKind { Null, Boolean, Integer, Rational, Text, Object, Array, Decimal, Tensor, ExactNumber, Bytes };
using NativeRecordObserver = std::function<void(std::string_view)>;

class NativeRecord {
    struct State;
    std::shared_ptr<const State> state_;
    rfc::composite::DeclarationId type_;
    const rfc::composite::DeclaredValue *value_{};
    std::shared_ptr<const rfc::composite::DeclaredValue> value_owner_;
    bool root_view_{};
    NativeRecordValue copy_value(const std::map<std::uint64_t,std::uint64_t> &,bool) const;
    NativeRecord(std::shared_ptr<const State>, rfc::composite::DeclarationId,
                 const rfc::composite::DeclaredValue *,
                 std::shared_ptr<const rfc::composite::DeclaredValue> = {});
    static NativeRecord create_owned(NativeRecordValue, std::uint64_t,
        std::vector<RawByteSource>, rfc::storage::TemporalEndpoint,
        const NativeRecord *previous,
        std::optional<rfc::storage::NativeTemporalSeries> history = {},
        const NativeRecordObserver &observer = {});
  public:
    struct PackedEncoding {
        rfc::composite::DeclarationId type;
        std::size_t bit_count;
        std::span<const std::byte> bytes;
    };
    // A new application publication has its own explicit namespace. Supplying
    // an endpoint retains its actual complex, tensor fields and existing codes.
    static NativeRecord create(const NativeRecordValue &value, std::uint64_t source_namespace,
        std::vector<RawByteSource> sources = {},
        const rfc::storage::TemporalEndpoint *endpoint = nullptr);
    // Consume temporary field storage before freezing the native publication.
    static NativeRecord create_from_endpoint(const NativeRecordValue &value,
        std::uint64_t source_namespace, std::vector<RawByteSource> sources,
        rfc::storage::TemporalEndpoint endpoint);
    static NativeRecord create_from_endpoint(NativeRecordValue &&value,
        std::uint64_t source_namespace, std::vector<RawByteSource> sources,
        rfc::storage::TemporalEndpoint endpoint, const NativeRecordObserver &observer = {});
    // Publish added fields on the source's unchanged complex and carriers under
    // a new namespace. Retain its declarations and immutable original bytes.
    static NativeRecord publish_fields(const NativeRecord &source, std::uint64_t source_namespace,
                                       std::vector<rfc::storage::TemporalTensor> fields,
                                       NativeRecordValue value,
                                       std::vector<RawByteSource> added_sources = {},
                                       std::vector<rfc::storage::TemporalCarrier> added_carriers = {});
    // Publish declarations after an actual owning history. Unchanged source
    // fields stay in its original base; the canonical frame retains the update.
    static NativeRecord annotate(const NativeRecordValue &value,
        rfc::storage::NativeTemporalSeries history, std::vector<RawByteSource> sources = {});
    static NativeRecord annotate(const NativeRecordValue &value,
        const NativeRecord &producer, std::vector<RawByteSource> sources = {});
    static NativeRecord from_history(rfc::storage::NativeTemporalSeries history);
    std::optional<rfc::storage::NativeTemporalSeries> publication_history() const;
    // Select a tensor from this publication's retained binary field storage.
    std::optional<rfc::storage::NativeEncodedTensor> encoded_tensor(std::uint64_t key) const;
    NativeRecord derive(const NativeRecordValue &value,std::vector<RawByteSource> added_sources = {}) const;
    NativeRecord derive(NativeRecordValue &&value,std::vector<RawByteSource> added_sources = {},
                        const NativeRecordObserver &observer = {}) const;
    static NativeRecord read(const std::filesystem::path &path);
    static std::optional<NativeRecord> read_if_record(const std::filesystem::path &path);
    static NativeRecord decode(std::span<const std::byte> bytes,
        rfc::storage::NativeSeriesLimits limits = native_data_limits());
    // Inspect retained byte identities without copying or materializing their payloads.
    static std::vector<std::uint64_t> source_keys(const rfc::storage::TemporalEndpoint &endpoint);
    // Historical byte references require from_history; a bare endpoint must be self-contained.
    static NativeRecord from_endpoint(const rfc::storage::TemporalEndpoint &endpoint,
                                      std::uint64_t source_namespace);
    void write(const std::filesystem::path &path) const;
    std::vector<std::byte> encode() const;
    // Encode once into a fresh immutable byte owner for archive composition.
    RawByteSource archive_source(std::uint64_t source_key) const;
    bool is_root() const noexcept;
    std::uint64_t format_version() const noexcept;
    // The root declaration shares this view's actual immutable publication.
    NativeRecord owning_publication() const;
    // Capture an actual ancestor; unchanged byte owners are validated and shared.
    NativeRecord historical_publication(std::uint64_t cutoff) const;
    NativeRecordKind kind() const;
    std::size_t size() const;
    bool has(std::string_view name) const;
    NativeRecord at(std::string_view name) const;
    NativeRecord at(std::size_t index) const;
    // The immutable parent owns this encoding. An ordinary array returns no
    // encoding; selection still checks the index against its actual count.
    std::optional<PackedEncoding> packed_encoding(std::size_t index) const;
    std::vector<std::string> names() const;
    bool boolean() const;
    rfc::numeric::Integer integer() const;
    rfc::numeric::Rational rational() const;
    bool approximate_view() const;
    std::string text() const;
    const rfc::numeric::SparseRationalTensor &tensor() const;
    std::span<const std::byte> bytes() const;
    RawByteSource byte_source() const;
    std::vector<RawByteSource> raw_sources() const;
    NativeRecordValue materialize_value(const std::map<std::uint64_t,std::uint64_t> &source_map = {}) const;
    // Copy a declaration into another complex without transferring cell or
    // occurrence bindings. Keep its original publication for those bindings.
    NativeRecordValue unbound_value(const std::map<std::uint64_t,std::uint64_t> &source_map = {}) const;
    std::optional<NativeRecord> parent() const;
    std::shared_ptr<const rfc::storage::TemporalEndpoint> endpoint() const;
    std::shared_ptr<const rfc::storage::NativeAccessionDirectory> accession_directory() const;
    std::uint64_t source_namespace() const;
    NativeRecord declared_value(std::uint64_t key) const;
    rfc::composite::DeclarationId type() const noexcept { return type_; }
    const rfc::composite::DeclaredValue &value() const noexcept { return *value_; }
};
// Host allocation of a new application publication identity; content equality
// does not identify sources. Persisted records keep their admitted namespace.
std::uint64_t new_native_record_namespace();
NativeRecord retain_operation(const NativeRecord &source,const NativeRecord &request,
                              NativeRecordValue result, const NativeRecordObserver &observer = {});
// Publish a newly constructed mathematical source and result together.
NativeRecord retain_operation(const rfc::storage::TemporalEndpoint &source,
    std::uint64_t source_namespace, std::vector<RawByteSource> sources,
    const NativeRecord &request, NativeRecordValue result);
NativeRecord retain_operation(rfc::storage::TemporalEndpoint &&source,
    std::uint64_t source_namespace, std::vector<RawByteSource> sources,
    const NativeRecord &request, NativeRecordValue result);
} // namespace gadget
