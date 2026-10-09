#include <gadget/native_data.hpp>
#include <algorithm>
#include <limits>
#include <set>
#include <atomic>
#include <mutex>

namespace gadget {
using namespace rfc;
namespace {
constexpr auto addressable = std::numeric_limits<std::size_t>::max();
void require(bool condition, const char *message) {
    if (!condition) throw std::invalid_argument(message);
}
std::size_t machine_integer(const composite::DeclaredValue &value) {
    const auto *integer = std::get_if<numeric::Integer>(&value.data);
    require(integer && *integer >= 0 && *integer <= numeric::Integer(std::to_string(addressable)),
            "raw byte span requires a nonnegative machine-addressable integer");
    return static_cast<std::size_t>(std::stoull(integer->get_str()));
}
} // namespace
NativeDataFormat detect_native_data(std::span<const std::byte> bytes) noexcept {
    const auto matches = [&](const auto &magic) {
        return bytes.size() >= magic.size() && std::equal(magic.begin(), magic.end(), bytes.begin());
    };
    if (matches(storage::native_series_codec_detail::magic)) return NativeDataFormat::TemporalSeries;
    if (matches(storage::rfcx_detail::magic)) return NativeDataFormat::RelationalComplex;
    return NativeDataFormat::Unknown;
}
storage::NativeAccessionLimits native_accession_limits() {
    storage::NativeAccessionLimits limits;
    limits.max_handles = limits.max_bucket_entries = limits.max_index_entries =
        limits.max_selection_cells = limits.max_descriptor_bytes = limits.max_total_descriptor_bytes =
        limits.max_binding_bytes = limits.max_factor_depth = limits.max_action_work = addressable;
    limits.max_integer_bytes = addressable / std::numeric_limits<unsigned char>::digits;
    return limits;
}
storage::NativeSeriesLimits native_data_limits() {
    storage::NativeSeriesLimits limits;
    limits.max_segments = limits.max_checkpoints = limits.max_archive_bytes =
        limits.max_analysis_steps = limits.max_analysis_dimension = limits.max_analysis_scalars =
        limits.max_analysis_work = addressable;
    auto &read = limits.temporal;
    read.max_file_bytes = read.max_records = read.max_dimension = read.max_rank =
        read.max_integer_bytes = addressable;
    read.source_declarations = {addressable, addressable};
    auto &d = read.declarations;
    d.max_definitions = d.max_components = d.max_value_nodes = d.max_encoded_bytes =
        d.max_integer_bits = d.max_materialized_bytes = addressable;
    // Use the native codec's supported stack/rank profile, not an app ceiling.
    d.max_depth = d.max_tensor_rank = 128;
    return limits;
}
RawByteSource::RawByteSource(std::uint64_t key, std::span<const std::byte> bytes)
    : key_(key), bytes_(std::make_shared<const std::vector<std::byte>>(bytes.begin(), bytes.end())) {
    require(key != 0, "raw source requires its declared identity");
}
RawByteSource RawByteSource::read(std::uint64_t key, const std::filesystem::path &path) {
    return adopt_fresh(key, storage::read_file_exact(path, addressable));
}
RawByteSource RawByteSource::own(std::uint64_t key, std::vector<std::byte> bytes) {
    RawByteSource source(key, {});
    source.bytes_ = std::make_shared<const std::vector<std::byte>>(std::move(bytes));
    source.frozen_ = false;
    return source;
}
RawByteSource RawByteSource::adopt_fresh(std::uint64_t key,std::vector<std::byte> bytes) {
    auto source=own(key,std::move(bytes));source.frozen_=true;return source;
}
RawByteSource RawByteSource::freeze() const {
    return frozen_?*this:RawByteSource(key_,bytes());
}
RawByteSource RawByteSource::rekey(std::uint64_t key) const {
    require(key != 0, "raw source requires its declared identity");
    auto source = freeze();
    source.key_ = key;
    return source;
}
PinnedRawSpan::PinnedRawSpan(RawByteSource source, storage::PinnedNativeAccession declaration,
                           std::size_t offset, std::size_t length)
    : source_(std::move(source)), declaration_(std::move(declaration)), offset_(offset), length_(length) {}
std::span<const std::byte> PinnedRawSpan::bytes() const noexcept {
    return source_.bytes().subspan(offset_, length_);
}
struct NativeStudy::State {
    mutable std::unique_ptr<storage::NativeAccessionDirectory> directory;
    std::optional<storage::NativeSeriesEndpoint> publication;
    storage::NativeAccessionLimits accession_limits;
    mutable std::once_flag accession_once;
    mutable std::mutex accession_mutex;
    mutable std::once_flag declaration_once;
    mutable std::unique_ptr<composite::PreparedLayoutRegistry> declarations;
    std::map<std::uint64_t, RawByteSource> sources;
    std::unique_ptr<field::ExactExecutionTower> tower;
    State(const storage::TemporalEndpoint &endpoint, std::uint64_t ns,
          const std::vector<RawByteSource> &raw, storage::NativeAccessionLimits limits)
        : directory(std::make_unique<storage::NativeAccessionDirectory>(endpoint, ns,
              std::vector<storage::RelationalCodeObservationQuery>{},
              [](auto bytes) { return storage::rfcx_detail::fnv1a64(bytes); }, limits)) {
        for (const auto &source : raw)
            require(sources.emplace(source.key(), source.freeze()).second, "duplicate raw source identity");
    }
    State(const core::RelationalComplex &source, std::vector<field::ExactMetric> metrics)
        : State(storage::TemporalEndpoint(source), ephemeral_namespace(), {}, native_accession_limits()) {
        tower = std::make_unique<field::ExactExecutionTower>(field::ExactExecutionTower::weighted(
            directory->pin_snapshot()->source, std::move(metrics)));
    }
    State(const storage::NativeSeriesEndpoint &endpoint,
          const std::vector<RawByteSource> &raw, storage::NativeAccessionLimits limits,
          bool defer_accessions = false)
        : publication(endpoint), accession_limits(limits) {
        (void)endpoint.pin();
        if (!defer_accessions)
            (void)prepare_accessions();
        for (const auto &source : raw)
            require(sources.emplace(source.key(), source.freeze()).second, "duplicate raw source identity");
    }
    const storage::NativeAccessionDirectory &prepare_accessions() const {
        if (publication)
            std::call_once(accession_once, [&] {
                directory = std::make_unique<storage::NativeAccessionDirectory>(
                    *publication, std::vector<storage::RelationalCodeObservationQuery>{},
                    [](auto bytes) { return storage::rfcx_detail::fnv1a64(bytes); }, accession_limits);
            });
        return *directory;
    }
    const composite::PreparedLayoutRegistry &prepare_declarations() const {
        std::call_once(declaration_once, [&] {
            const auto endpoint = snapshot();
            require(endpoint->declarations.has_value(), "native source has no declarations");
            declarations = std::make_unique<composite::PreparedLayoutRegistry>(
                endpoint->declarations->registry, native_data_limits().temporal.declarations);
        });
        return *declarations;
    }
    std::shared_ptr<const storage::TemporalEndpoint> snapshot() const {
        if (publication)
            return publication->pin();
        const std::lock_guard lock(accession_mutex);
        return directory->pin_snapshot();
    }
    static std::uint64_t ephemeral_namespace() {
        static std::atomic<std::uint64_t> next{1};
        auto value=next.load(std::memory_order_relaxed);
        for(;;) {
            if(!value)throw std::overflow_error("ephemeral native source identity exhausted");
            if(next.compare_exchange_weak(value,value==std::numeric_limits<std::uint64_t>::max()?0:value+1,std::memory_order_relaxed))return value;
        }
    }
};
NativeStudy::NativeStudy(std::shared_ptr<const core::RelationalComplex> source,
                         std::vector<field::ExactMetric> metrics) {
    require(static_cast<bool>(source), "native study requires an immutable source owner");
    state_ = std::make_shared<const State>(*source, std::move(metrics));
}
NativeStudy::NativeStudy(core::RelationalComplex source, std::vector<field::ExactMetric> metrics)
    : state_(std::make_shared<const State>(source, std::move(metrics))) {}
NativeStudy::NativeStudy(std::shared_ptr<const State> state) : state_(std::move(state)) {}
NativeStudy::NativeStudy(const storage::TemporalEndpoint &endpoint, std::uint64_t ns,
                         std::vector<RawByteSource> sources, storage::NativeAccessionLimits limits)
    : NativeStudy(std::make_shared<const State>(endpoint, ns, sources, limits)) {}
NativeStudy::NativeStudy(const storage::NativeSeriesEndpoint &endpoint,
                         std::vector<RawByteSource> sources, storage::NativeAccessionLimits limits)
    : NativeStudy(std::make_shared<const State>(endpoint, sources, limits)) {}
NativeStudy NativeStudy::from_publication(const storage::NativeSeriesEndpoint &endpoint,
                                        std::vector<RawByteSource> sources) {
    return NativeStudy(std::make_shared<const State>(endpoint, sources, native_accession_limits(), true));
}
NativeStudy NativeStudy::from_history(const storage::NativeTemporalSeries &history,
                                     std::vector<RawByteSource> sources,
                                     std::optional<std::uint64_t> cutoff) {
    const auto endpoint = !cutoff || *cutoff == history.last_sequence()
        ? history.tail() : history.at(*cutoff);
    return from_publication(endpoint, std::move(sources));
}
NativeStudy NativeStudy::read(const std::filesystem::path &path, std::vector<RawByteSource> sources,
                             std::optional<std::uint64_t> ns, std::optional<std::uint64_t> cutoff,
                             storage::NativeSeriesLimits limits) {
    const storage::ReadOnlyMappedFile image(path);
    require(image.size() <= limits.max_archive_bytes, "native input exceeds caller byte resources");
    switch (detect_native_data(image.bytes())) {
    case NativeDataFormat::TemporalSeries: {
        const auto series = storage::decode_native_series(image.bytes(), limits);
        require(!ns || *ns == series.source_namespace(), "native source namespace mismatch");
        return from_history(series, std::move(sources), cutoff);
    }
    case NativeDataFormat::RelationalComplex: {
        require(ns && *ns, "RFCX input requires an explicit source namespace");
        require(!cutoff || *cutoff == 0, "RFCX input has no temporal history cutoff");
        return NativeStudy(storage::TemporalEndpoint(storage::MappedRfcx(image.bytes()).materialize()),
                           *ns, std::move(sources));
    }
    default: throw std::invalid_argument("unrecognized canonical native data bytes");
    }
}
std::shared_ptr<const storage::TemporalEndpoint> NativeStudy::endpoint() const {
    return state_->snapshot();
}
std::shared_ptr<const storage::NativeAccessionDirectory> NativeStudy::accession_directory() const {
    return {state_,&state_->prepare_accessions()};
}
const core::RelationalComplex &NativeStudy::source() const { return endpoint()->source; }
const field::ExactExecutionTower &NativeStudy::tower() const {
    require(static_cast<bool>(state_->tower), "native study requires explicitly declared exact metrics");
    return *state_->tower;
}
field::ExactHodgeOperator NativeStudy::hodge(core::Grade grade) const {
    return {source(), grade, tower().metrics()};
}
numeric::RationalVector NativeStudy::execute(const field::OperatorRequest &request,
                                            const numeric::RationalVector &input) const {
    return tower().apply(request, input);
}
numeric::Rational NativeStudy::moment(core::Grade grade, std::size_t order,
    const numeric::RationalVector &x, const numeric::RationalVector &y) const {
    auto op = hodge(grade);
    return op.metric().inner(x, statistics::hodge_power_apply(op, order, y));
}
field::ExactGreenResult NativeStudy::green_action(core::Grade grade, const numeric::Rational &lambda,
                                                const numeric::RationalVector &input) const {
    return field::exact_green_apply(hodge(grade), lambda, input);
}
storage::NativeAtlasSource NativeStudy::atlas_source(std::uint64_t key) const {
    const auto &directory = state_->prepare_accessions();
    const std::lock_guard lock(state_->accession_mutex);
    return {key, directory};
}
storage::PinnedNativeAccession NativeStudy::access(const storage::NativeAccessionRequest &request) const {
    const auto &directory = state_->prepare_accessions();
    const std::lock_guard lock(state_->accession_mutex);
    return directory.pin(request);
}
PinnedRawSpan NativeStudy::raw_span(std::uint64_t key, const RawSpanDeclaration &binding) const {
    require(binding.type.valid() && binding.source_factor && binding.offset_factor && binding.length_factor &&
                std::set<std::uint64_t>{binding.source_factor, binding.offset_factor, binding.length_factor}.size() == 3,
            "raw span requires three distinct authored product factors");
    storage::NativeAccessionRequest request;
    request.kind = storage::NativeAccessionKind::DeclaredValue;
    request.subject = key;
    const auto pin = access(request);
    require(pin.declared_owner()->type == binding.type, "raw span declaration type mismatch");
    const auto endpoint = pin.endpoint();
    const auto &root = endpoint->declarations->registry.at(binding.type);
    require(root.kind == composite::DeclarationKind::Product, "raw span declaration must be a Product");
    const auto &parts = std::get<composite::DeclaredValue::Product>(pin.declared_value()->data).elements;
    const auto factor = [&](std::uint64_t id) {
        const auto found = std::find_if(root.factors.begin(), root.factors.end(),
                                       [&](const auto &f) { return f.key == id; });
        require(found != root.factors.end() && endpoint->declarations->registry.at(found->type).kind ==
                    composite::DeclarationKind::Integer, "raw span factor is absent or has the wrong type");
        return machine_integer(parts.at(static_cast<std::size_t>(found - root.factors.begin())));
    };
    const auto source = factor(binding.source_factor), offset = factor(binding.offset_factor),
               length = factor(binding.length_factor);
    require(state_->sources.contains(source), "raw span source bytes are unavailable");
    const auto &raw = state_->sources.at(source);
    require(offset <= raw.bytes().size() && length <= raw.bytes().size() - offset,
            "raw span is outside its original source bytes");
    return {raw, pin, offset, length};
}
storage::PreparedNativeSourceProgram NativeStudy::program(
    std::uint64_t key, std::uint64_t layout, storage::StableOccurrenceKey occurrence) const {
    storage::NativeAccessionRequest request;
    request.kind = storage::NativeAccessionKind::SourceAspects;
    request.subject = layout;
    request.occurrence = occurrence;
    composite::BinaryExecutionLimits limits;
    limits.max_instructions = limits.max_outputs = limits.max_output_words = addressable;
    return storage::PreparedNativeSourceProgram(access(request), key, limits);
}
DeclaredByteValue NativeStudy::decode_raw_span(std::uint64_t value,
    const RawSpanDeclaration &span, composite::DeclarationId type, std::size_t bit_count) const {
    auto source = raw_span(value, span);
    const auto endpoint = source.declaration().endpoint();
    composite::CompiledLayout layout(state_->prepare_declarations(), type,
        endpoint->declarations->profile);
    auto decoded = storage::decode_declared_bits(layout, source.bytes(), bit_count);
    return {std::move(source), std::move(layout), std::move(decoded)};
}
} // namespace gadget
