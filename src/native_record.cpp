#include <cerrno>
#include <cstring>
#include <gadget/native_record.hpp>
#include <iterator>
#include <rfc/core/builder.hpp>
#include <rfc/runtime/hardware_inventory.hpp>
#include <rfc/storage/file_io.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <rfc/storage/native_series_prefix.hpp>
#include <set>
#include <sys/random.h>
#include <tuple>

namespace gadget {
namespace {
using namespace rfc;
using V = composite::DeclaredValue;
using D = composite::Declaration;
using ID = composite::DeclarationId;
constexpr std::uint64_t format_namespace = 0x4741444745540201ULL;
constexpr std::uint64_t header_key = 0x4741444745540201ULL;
constexpr std::uint64_t root_key = header_key + 1, labels_key = header_key + 2,
                        sources_key = header_key + 3;
constexpr ID integer_type{format_namespace, 1, 1}, boolean_type{format_namespace, 2, 1},
    rational_type{format_namespace, 3, 1}, null_type{format_namespace, 4, 1},
    text_type{format_namespace, 5, 1}, source_type{format_namespace, 6, 1},
    sources_type{format_namespace, 7, 1}, label_type{format_namespace, 8, 1},
    labels_type{format_namespace, 9, 1}, header_type{format_namespace, 10, 1},
    decimal_type{format_namespace, 11, 1}, exact_number_type{format_namespace, 12, 1},
    bytes_type{format_namespace, 13, 1}, owned_reference_type{format_namespace, 14, 1},
    historical_source_type{format_namespace, 15, 1}, source_payload_type{format_namespace, 16, 1},
    referenced_source_type{format_namespace, 17, 1},
    referenced_sources_type{format_namespace, 18, 1},
    packed_value_type{format_namespace, 19, 1};
constexpr auto addressable = std::numeric_limits<std::size_t>::max();
void require(bool ok, const char *message) {
    if (!ok)
        throw std::invalid_argument(message);
}
V number(std::uint64_t n) { return V::integer(numeric::Integer(std::to_string(n))); }
std::uint64_t unsigned_value(const V &value) {
    const auto &n = std::get<numeric::Integer>(value.data);
    require(n >= 0 &&
                n <= numeric::Integer(std::to_string(std::numeric_limits<std::uint64_t>::max())),
            "native record integer is outside its declared host identity domain");
    return std::stoull(n.get_str());
}
const auto &parts(const V &v) { return std::get<V::Product>(v.data).elements; }
std::vector<D> definitions(bool historical_sources = true, bool packed_values = false) {
    std::vector<D> out{
        D::integer(integer_type),
        D::boolean(boolean_type),
        D::rational(rational_type),
        D::product(null_type, {}),
        D::product(text_type,
                   {{1, integer_type, {}}, {2, integer_type, {}}, {3, integer_type, {}}}),
        D::product(source_type,
                   {{1, integer_type, {}}, {2, integer_type, {}}, {3, integer_type, {}}}),
        D::repeated(sources_type, source_type, 0, addressable),
        D::product(label_type, {{1, integer_type, {}}, {2, text_type, {}}}),
        D::repeated(labels_type, label_type, 0, addressable),
        D::product(header_type, {{1, integer_type, {}},
                                 {2, integer_type, {}},
                                 {3, integer_type, {}},
                                 {4, integer_type, {}},
                                 {5, integer_type, {}},
                                 {6, integer_type, {}},
                                 {7, integer_type, {}}}),
        D::product(decimal_type, {{1, rational_type, {}}, {2, text_type, {}}}),
        D::product(exact_number_type, {{1, rational_type, {}}, {2, boolean_type, {}}}),
        D::product(bytes_type,
                   {{1, integer_type, {}}, {2, integer_type, {}}, {3, integer_type, {}}}),
        D::product(owned_reference_type, {{1, integer_type, {}},
                                          {2, integer_type, {}},
                                          {3, integer_type, {}},
                                          {4, integer_type, {}}})};
    if (historical_sources) {
        out.push_back(
            D::product(historical_source_type,
                       {{1, integer_type, {}}, {2, integer_type, {}}, {3, integer_type, {}}}));
        out.push_back(D::alternative(source_payload_type,
                                     {{1, integer_type, {}}, {2, historical_source_type, {}}}));
        out.push_back(D::product(
            referenced_source_type,
            {{1, integer_type, {}}, {2, integer_type, {}}, {3, source_payload_type, {}}}));
        out.push_back(D::repeated(referenced_sources_type, referenced_source_type, 0, addressable));
    }
    if (packed_values)
        out.push_back(D::product(packed_value_type,
            {{1, integer_type, {}}, {2, integer_type, {}}, {3, integer_type, {}},
             {4, integer_type, {}}, {5, integer_type, {}}, {6, integer_type, {}},
             {7, integer_type, {}}}));
    return out;
}
struct Compiler {
    std::uint64_t ns, next{1}, pool;
    std::vector<D> defs = definitions();
    std::map<std::string,std::uint64_t,std::less<>> name_keys;
    std::vector<V> labels;
    std::vector<std::byte> text_bytes;
    std::map<std::string, std::pair<std::size_t, std::size_t>, std::less<>> text_spans;
    std::map<std::string,ID> shapes;
    std::map<ID, composite::CompiledLayout> layouts;
    bool packed_values{};
    std::set<std::uint64_t> source_keys,value_keys;
    std::vector<storage::TemporalDeclaredValue> owned_values;
    std::uint64_t next_value{1};
    Compiler(std::uint64_t n, const std::vector<RawByteSource> &sources) : ns(n), pool(1) {
        require(ns && ns != format_namespace, "native record requires an independent declaration namespace");
        for (const auto &s:sources)
            require(source_keys.insert(s.key()).second,"duplicate native record raw source identity");
        while (source_keys.contains(pool)) {
            require(pool != std::numeric_limits<std::uint64_t>::max(),"raw source identity space exhausted");
            ++pool;
        }
    }
    V span(std::uint64_t source, std::size_t offset, std::size_t length) {
        return V::product({number(source),number(offset),number(length)});
    }
    V text(std::string_view bytes) {
        if (const auto found = text_spans.find(bytes); found != text_spans.end())
            return span(pool, found->second.first, found->second.second);
        const auto offset=text_bytes.size();
        (void)machine::checked_size_add(offset,bytes.size());
        for (const unsigned char byte:bytes) text_bytes.push_back(static_cast<std::byte>(byte));
        text_spans.emplace(std::string(bytes), std::make_pair(offset, bytes.size()));
        return span(pool,offset,bytes.size());
    }
    ID shape(std::string key, D definition) {
        const auto found=shapes.find(key);
        if(found!=shapes.end()) return found->second;
        require(next != std::numeric_limits<std::uint64_t>::max(),"native declaration identity space exhausted");
        definition.id={ns,next++,1};
        const auto id=definition.id;defs.push_back(std::move(definition));shapes.emplace(std::move(key),id);
        return id;
    }
    std::uint64_t name(std::string_view label) {
        const auto found=name_keys.find(std::string(label));
        if(found!=name_keys.end())return found->second;
        const auto key=machine::checked_size_add(name_keys.size(),1);
        name_keys.emplace(label,key);labels.push_back(V::product({number(key),text(label)}));return key;
    }
    static std::string signature(ID id) {
        return std::to_string(id.name_space)+":"+std::to_string(id.key)+":"+std::to_string(id.version);
    }
    V pack(NativeRecordValue input, std::size_t depth) {
        auto [type, value] = compile(std::move(input), depth);
        if (!layouts.contains(type))
            layouts.emplace(type, composite::CompiledLayout(
                composite::DeclarationRegistry(defs, native_data_limits().temporal.declarations),
                type, composite::BooleanLayoutProfile::Lsb0OctetsV1,
                native_data_limits().temporal.declarations));
        const auto bits = storage::encode_declared_bits(layouts.at(type), value);
        const auto offset = text_bytes.size();
        (void)machine::checked_size_add(offset, bits.bytes.size());
        text_bytes.insert(text_bytes.end(), bits.bytes.begin(), bits.bytes.end());
        return V::product({number(pool), number(offset), number(bits.bytes.size()),
            number(type.name_space), number(type.key), number(type.version), number(bits.bit_count)});
    }
    std::pair<ID,V> compile(NativeRecordValue input, std::size_t depth=1) {
        require(depth<=native_data_limits().temporal.declarations.max_depth,"native composite nesting profile exceeded");
        return std::visit([&](auto &x)->std::pair<ID,V>{
            using T=std::decay_t<decltype(x)>;
            if constexpr(std::is_same_v<T,std::monostate>) return {null_type,V::product({})};
            else if constexpr(std::is_same_v<T,bool>) return {boolean_type,V::boolean(x)};
            else if constexpr(std::is_same_v<T,numeric::Integer>) return {integer_type,V::integer(std::move(x))};
            else if constexpr(std::is_same_v<T,numeric::Rational>) return {rational_type,V::rational(std::move(x))};
            else if constexpr(std::is_same_v<T,std::string>) return {text_type,text(x)};
            else if constexpr(std::is_same_v<T,NativeRecordValue::Owned>) {
                require(static_cast<bool>(x.payload),"owned declaration requires its complete payload");
                auto [type,value]=compile(*x.payload,depth+1);
                while(value_keys.contains(next_value)) {require(next_value!=std::numeric_limits<std::uint64_t>::max(),"owned value identity exhausted");++next_value;}
                const auto key=next_value;value_keys.insert(key);
                owned_values.push_back({key,type,x.block,0,x.cell,x.occurrence,std::move(value)});
                return {owned_reference_type,V::product({number(key),number(type.name_space),number(type.key),number(type.version)})};
            }
            else if constexpr(std::is_same_v<T,NativeRecordValue::ByteSpan>) {
                require(source_keys.contains(x.source),"native record span has no original byte source");
                return {x.text?text_type:bytes_type,span(x.source,x.offset,x.length)};
            } else if constexpr(std::is_same_v<T,NativeRecordValue::Decimal>) {
                return {decimal_type,V::product({V::rational(x.exact),text(x.literal)})};
            } else if constexpr(std::is_same_v<T,NativeRecordValue::ExactNumber>) {
                return {exact_number_type,V::product({V::rational(x.exact),V::boolean(x.approximate_view)})};
            } else if constexpr(std::is_same_v<T,NativeRecordValue::Object>) {
                std::vector<composite::DeclarationFactor> factors;std::vector<V> values;
                std::string key="product";std::set<std::string> seen;
                values.reserve(x.members.size()); factors.reserve(x.members.size());
                for(auto &[label,item]:x.members) {
                    require(seen.insert(label).second,"duplicate native record member");
                    const auto factor=name(label);auto [type,value]=compile(std::move(item),depth+1);
                    factors.push_back({factor,type,{}});values.push_back(std::move(value));
                    key+="/"+std::to_string(factor)+":"+signature(type);
                }
                return {shape(key,D::product({},std::move(factors))),V::product(std::move(values))};
            } else if constexpr(std::is_same_v<T,NativeRecordValue::Array>) {
                std::vector<std::pair<ID,V>> children;std::map<ID,std::uint64_t> types;
                children.reserve(x.elements.size());
                for(auto &item:x.elements) {children.push_back(compile(std::move(item),depth+1));types.emplace(children.back().first,0);}
                ID element=null_type;
                if(types.size()==1)element=types.begin()->first;
                else if(types.size()>1) {
                    std::vector<composite::DeclarationFactor> alternatives;std::string key="alternative";
                    for(auto &[type,tag]:types) {tag=alternatives.size()+1;alternatives.push_back({tag,type,{}});key+="/"+signature(type);}
                    element=shape(key,D::alternative({},std::move(alternatives)));
                }
                // Cardinality is part of this admitted value's declaration. A host
                // address-space maximum would declare nonexistent fixed-width
                // elements and can overflow native layout bounds before encoding.
                const auto count=x.elements.size();
                const auto type=shape("repeated/"+signature(element)+"/"+std::to_string(count),D::repeated({},element,count,count));
                std::vector<V> values;values.reserve(children.size());
                for(auto &[child,value]:children)values.push_back(types.size()>1?V::alternative(types.at(child),std::move(value)):std::move(value));
                return {type,V::repeated(std::move(values))};
            } else if constexpr(std::is_same_v<T,NativeRecordValue::PackedArray>) {
                require(static_cast<bool>(x.emit), "packed array requires an owning producer");
                if (!packed_values) {
                    defs.push_back(definitions(true, true).back());
                    packed_values = true;
                }
                std::vector<V> values;
                values.reserve(x.count);
                x.emit([&](NativeRecordValue item) {
                    require(values.size() < x.count, "packed array producer exceeded its declared cardinality");
                    values.push_back(pack(std::move(item), depth + 1));
                });
                require(values.size() == x.count, "packed array producer did not supply its declared cardinality");
                const auto type = shape("packed-repeated/" + std::to_string(x.count),
                    D::repeated({}, packed_value_type, x.count, x.count));
                return {type, V::repeated(std::move(values))};
            } else {
                std::vector<composite::DeclarationAxis> axes;std::string key="tensor";
                for(const auto extent:x.domain().extents()) {axes.push_back({axes.size()+1,extent,{}});key+="/"+std::to_string(extent);}
                return {shape(key,D::rational_tensor({},std::move(axes))),V::tensor(std::move(x))};
            }
        },input.data);
    }
};
} // namespace

struct NativeRecord::State {
    NativeStudy study;
    // The study owns this immutable endpoint throughout the State lifetime.
    const storage::TemporalEndpoint &publication;
    std::optional<storage::NativeTemporalSeries> history;
    std::map<std::uint64_t,RawByteSource> sources;
    std::map<std::uint64_t,std::string> labels;
    std::map<std::string,std::uint64_t,std::less<>> names;
    ID root_type;
    const V *root{};
    std::uint64_t ns, version{}, root_value_key{}, parent_value_key{};
    struct SourceOrigin {
        std::uint64_t sequence, key;
        auto operator<=>(const SourceOrigin &) const = default;
    };
    struct SourcePayload {
        const numeric::Integer *integer;
        SourceOrigin origin;
        std::size_t length;
    };
    mutable std::map<std::uint64_t, storage::NativeSeriesEndpoint> source_endpoints;
    std::map<std::uint64_t, SourcePayload> source_payloads;
    SourcePayload source_payload(const V &entry, std::uint64_t sequence) const {
        const V *current = &entry;
        const storage::TemporalEndpoint *endpoint = &publication;
        const auto expected_length = machine::checked_size(unsigned_value(parts(entry).at(1)));
        for (;;) {
            const auto &items = parts(*current);
            const auto key = unsigned_value(items.at(0));
            const auto length = machine::checked_size(unsigned_value(items.at(1)));
            require(length == expected_length, "historical byte source length mismatch");
            if (const auto *integer = std::get_if<numeric::Integer>(&items.at(2).data))
                return {integer, {sequence, key}, length};
            const auto &alternative = std::get<V::Alternative>(items.at(2).data);
            require(alternative.element.size() == 1, "invalid byte source alternative");
            if (alternative.key == 1)
                return {&std::get<numeric::Integer>(alternative.element[0].data),
                        {sequence, key},
                        length};
            require(alternative.key == 2 && history.has_value(),
                    "historical byte source requires its owning history");
            const auto &reference = parts(alternative.element[0]);
            const auto source_namespace = unsigned_value(reference.at(0)),
                       previous = unsigned_value(reference.at(1)),
                       source_key = unsigned_value(reference.at(2));
            require(source_namespace == ns && source_key && previous < sequence &&
                        previous >= history->first_sequence(),
                    "historical byte source is outside its owning prefix");
            if (!source_endpoints.contains(previous))
                source_endpoints.emplace(
                    previous, storage::NativeSeriesPrefix::capture(*history, previous).tail());
            const auto &publication = source_endpoints.at(previous);
            require(publication.source_namespace() == source_namespace,
                    "historical byte source namespace mismatch");
            endpoint = &publication.endpoint();
            sequence = previous;
            require(endpoint->declarations.has_value(),
                    "historical byte source declarations missing");
            const storage::TemporalDeclaredValue *header = nullptr, *sources = nullptr;
            for (const auto &value : endpoint->declarations->values) {
                if (value.key == header_key)
                    header = &value;
                if (value.key == sources_key)
                    sources = &value;
            }
            require(header && sources && header->type == header_type,
                    "historical record control missing");
            const auto &h = parts(header->value);
            const auto revision = unsigned_value(h.at(0));
            require((revision == 2 || revision == 3 || revision == 4) && unsigned_value(h.at(1)) == ns,
                    "historical record identity mismatch");
            for (const auto &definition : definitions(revision >= 3, revision >= 4))
                require(endpoint->declarations->registry.at(definition.id) == definition,
                        "historical record declaration changed");
            require(sources->type == (revision >= 3 ? referenced_sources_type : sources_type),
                    "historical byte source declaration mismatch");
            current = nullptr;
            for (const auto &item : std::get<V::Repeated>(sources->value.data).elements)
                if (unsigned_value(parts(item).at(0)) == source_key) {
                    require(!current, "duplicate historical byte source");
                    current = &item;
                }
            require(current, "historical byte source key unavailable");
        }
    }
    std::map<std::uint64_t, const storage::TemporalDeclaredValue *> values;
    std::map<const V *, const storage::TemporalDeclaredValue *> owners;
    using PackedAddress = std::tuple<ID, std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t>;
    static PackedAddress packed_address(const V &value) {
        const auto &p = parts(value);
        return {{unsigned_value(p.at(3)), unsigned_value(p.at(4)), unsigned_value(p.at(5))},
                unsigned_value(p.at(0)), unsigned_value(p.at(1)), unsigned_value(p.at(2)), unsigned_value(p.at(6))};
    }
    std::pair<ID, std::shared_ptr<const V>> unpack(const V &value) const {
        const auto &[type, source, offset, length, bits] = packed_address(value);
        (void)source; (void)offset; (void)length;
        const auto layout = study.compiled_layout(type);
        const auto bytes = span(value);
        return {type, std::make_shared<const V>(storage::decode_declared_bits(
            *layout, bytes, machine::checked_size(bits)))};
    }
    State(const storage::TemporalEndpoint &endpoint, std::uint64_t source_namespace,
          std::vector<RawByteSource> compiled_sources = {})
        : State(NativeStudy(endpoint, source_namespace), source_namespace,
                std::move(compiled_sources)) {}
    State(const storage::NativeSeriesEndpoint &endpoint,
          std::vector<RawByteSource> compiled_sources, storage::NativeTemporalSeries retained,
          const State *reuse = nullptr)
        : State(NativeStudy::from_publication(endpoint), endpoint.source_namespace(),
                std::move(compiled_sources), std::move(retained), reuse) {}
    State(NativeStudy admitted,std::uint64_t source_namespace,
          std::vector<RawByteSource> compiled_sources,
          std::optional<storage::NativeTemporalSeries> retained={}, const State *reuse=nullptr)
        :study(std::move(admitted)),publication(*study.endpoint()),history(std::move(retained)),ns(source_namespace) {
        const auto &endpoint=publication;
        // Only the private compiler supplies these immutable owners. Every
        // source Integer was built from these bytes, or retained from its
        // admitted parent. Readers still materialize the validated cargo.
        const bool compiled=!compiled_sources.empty();
        for(auto &source:compiled_sources)
            require(sources.emplace(source.key(),std::move(source)).second,"duplicate compiled byte source");
        require(endpoint.declarations.has_value(),"native record has no declaration channel");
        for(const auto &v:endpoint.declarations->values) {values.emplace(v.key,&v);if(v.cell || v.occurrence || v.block)owners.emplace(&v.value,&v);}
        const auto declared=[&](std::uint64_t key) -> const storage::TemporalDeclaredValue & {
            require(values.contains(key),"native record control declaration is missing");
            return *values.at(key);
        };
        const auto &header=declared(header_key);
        require(header.type==header_type,"native record header declaration mismatch");
        const auto &h=parts(header.value);
        version = unsigned_value(h.at(0));
        require((version == 2 || version == 3 || version == 4) && unsigned_value(h.at(1)) == source_namespace,
                "native record version/namespace mismatch");
        for (const auto &definition : definitions(version >= 3, version >= 4))
            require(endpoint.declarations->registry.at(definition.id) == definition,
                    "native record format declaration was changed");
        root_type = {unsigned_value(h.at(2)), unsigned_value(h.at(3)), unsigned_value(h.at(4))};
        root_value_key = unsigned_value(h.at(5));
        parent_value_key = unsigned_value(h.at(6));
        require(root_value_key && root_value_key != header_key && root_value_key != labels_key &&
                    root_value_key != sources_key,
                "invalid native root selector");
        const auto &raw = declared(sources_key);
        require(raw.type == (version >= 3 ? referenced_sources_type : sources_type),
                "native record source declaration mismatch");
        std::size_t raw_materialization{};
        std::map<SourceOrigin, RawByteSource> materialized;
        for (const auto &entry : std::get<V::Repeated>(raw.value.data).elements) {
            const auto &p = parts(entry);
            const auto key = unsigned_value(p.at(0));
            const auto length = machine::checked_size(unsigned_value(p.at(1)));
            const auto payload = source_payload(entry, endpoint.sequence);
            require(source_payloads.emplace(key, payload).second, "duplicate native byte source");
            const auto &integer = *payload.integer;
            require(integer >= 0, "negative opaque native byte payload");
            require(composite::declaration_detail::integer_bits(integer) <=
                        machine::checked_size_mul(length, 8),
                    "opaque native byte payload exceeds length");
            if (!materialized.contains(payload.origin))
                raw_materialization = machine::checked_size_add(raw_materialization, length);
            require(raw_materialization <=
                        endpoint.declarations->resource_limits.max_materialized_bytes,
                    "native raw source materialization exceeds caller resources");
            if (compiled) {
                require(sources.contains(key) && sources.at(key).bytes().size() == length,
                        "compiled byte source identity/length mismatch");
            } else {
                if (materialized.contains(payload.origin)) {
                    require(sources.emplace(key, materialized.at(payload.origin).rekey(key)).second,
                            "duplicate native byte source");
                } else if (reuse && reuse->sources.contains(key) &&
                           reuse->sources.at(key).bytes().size() == length &&
                           *reuse->source_payloads.at(key).integer == integer) {
                    require(sources.emplace(key, reuse->sources.at(key)).second,
                            "duplicate native byte source");
                } else {
                    std::vector<std::byte> bytes(length, std::byte{0});
                    std::size_t count{};
                    if (length)
                        mpz_export(bytes.data(), &count, -1, 1, 0, 0, integer.get_mpz_t());
                    require(sources.emplace(key, RawByteSource::adopt_fresh(key, std::move(bytes)))
                                .second,
                            "duplicate native byte source");
                }
            }
            materialized.emplace(payload.origin, sources.at(key));
        }
        require(!compiled ||
                    sources.size() == std::get<V::Repeated>(raw.value.data).elements.size(),
                "compiled byte source count mismatch");
        const auto &label = declared(labels_key);
        require(label.type == labels_type, "native record label declaration mismatch");
        for (const auto &entry : std::get<V::Repeated>(label.value.data).elements) {
            const auto &p = parts(entry);
            const auto key = unsigned_value(p.at(0));
            const auto s = span(p.at(1));
            const std::string text =
                s.empty() ? std::string{}
                          : std::string(reinterpret_cast<const char *>(s.data()), s.size());
            require(key && labels.emplace(key,text).second && names.emplace(text,key).second,"duplicate native factor name");
        }
        const auto &data=declared(root_value_key);
        require(data.type==root_type,"native record root declaration mismatch");
        root=&data.value;
        // The canonical validator checks all declaration constructors. This
        // traversal additionally checks every physical byte address and label.
        struct Visit {
            ID type; const V *value; std::uint64_t completed_reference{};
            std::shared_ptr<const V> owner;
            std::optional<PackedAddress> completed_packed;
            std::size_t next{};
        };
        std::vector<Visit> pending{{root_type,root,0,{}, {}}};
        // Retained declarations can contain packed values even when a later
        // root no longer selects them. Validate their byte addresses as well.
        std::map<ID, bool> contains_packed;
        const std::function<bool(ID)> has_packed = [&](ID type) {
            if (const auto found = contains_packed.find(type); found != contains_packed.end())
                return found->second;
            const auto &declaration = endpoint.declarations->registry.at(type);
            bool found = type == packed_value_type;
            if (declaration.kind == composite::DeclarationKind::Repeated)
                found = found || has_packed(declaration.element);
            else if (declaration.kind == composite::DeclarationKind::Product ||
                     declaration.kind == composite::DeclarationKind::Alternative)
                for (const auto &factor : declaration.factors)
                    found = found || has_packed(factor.type);
            contains_packed.emplace(type, found);
            return found;
        };
        for (const auto &value : endpoint.declarations->values)
            if (value.key != root_value_key && has_packed(value.type))
                pending.push_back({value.type, &value.value, 0, {}, {}});
        std::set<std::uint64_t> active, completed;
        std::set<PackedAddress> active_packed, completed_packed;
        while(!pending.empty()) {
            auto &visit = pending.back();
            const auto type = visit.type; const auto *value = visit.value;
            const auto exit_key = visit.completed_reference;
            if(exit_key) {
                active.erase(exit_key); completed.insert(exit_key);
                pending.pop_back(); continue;
            }
            if (visit.completed_packed) {
                active_packed.erase(*visit.completed_packed);
                completed_packed.insert(*visit.completed_packed);
                pending.pop_back();
                continue;
            }
            if (type == packed_value_type) {
                require(version >= 4, "packed application values require record format 4");
                const auto address = packed_address(*value);
                require(!active_packed.contains(address), "cyclic packed application value");
                if (!completed_packed.contains(address)) {
                    require(active_packed.size() < endpoint.declarations->resource_limits.max_depth,
                            "packed application reference nesting exceeds declaration profile");
                    const auto [decoded_type, decoded] = unpack(*value);
                    active_packed.insert(address);
                    visit.completed_packed = address;
                    pending.push_back({decoded_type, decoded.get(), 0, decoded, {}});
                } else pending.pop_back();
                continue;
            }
            if(type==owned_reference_type) {
                const auto &p=parts(*value);const auto key=unsigned_value(p.at(0));
                require(values.contains(key),"owned application value absent");const auto &target=*values.at(key);
                require(target.type==ID{unsigned_value(p.at(1)),unsigned_value(p.at(2)),unsigned_value(p.at(3))},"owned application value type differs from reference");
                require(!active.contains(key),"cyclic application value reference");
                if(!completed.contains(key)) {
                    active.insert(key); visit.completed_reference = key;
                    pending.push_back({target.type,&target.value,0,{}, {}});
                } else pending.pop_back();
                continue;
            }
            if(type==text_type || type==bytes_type) {
                (void)span(*value); pending.pop_back(); continue;
            }
            if(type==decimal_type) {
                (void)span(parts(*value).at(1)); pending.pop_back(); continue;
            }
            if(type==exact_number_type) { pending.pop_back(); continue; }
            const auto &d=endpoint.declarations->registry.at(type);
            if(d.kind==composite::DeclarationKind::Product) {
                const auto &xs=parts(*value);
                if(visit.next==d.factors.size()) { pending.pop_back(); continue; }
                const auto i=visit.next++;
                require(labels.contains(d.factors[i].key),"native product has no factor name binding");
                const auto owner=visit.owner;
                pending.push_back({d.factors[i].type,&xs.at(i),0,owner,{}});
            } else if(d.kind==composite::DeclarationKind::Repeated) {
                const auto &xs=std::get<V::Repeated>(value->data).elements;
                if(visit.next==xs.size()) { pending.pop_back(); continue; }
                const auto i=visit.next++;
                const auto owner=visit.owner;
                pending.push_back({d.element,&xs.at(i),0,owner,{}});
            } else if(d.kind==composite::DeclarationKind::Alternative) {
                if(visit.next++) { pending.pop_back(); continue; }
                const auto &a=std::get<V::Alternative>(value->data);
                const auto i=endpoint.declarations->registry.alternative_ordinal(type,a.key);
                const auto owner=visit.owner;
                pending.push_back({d.factors.at(i).type,&a.element.at(0),0,owner,{}});
            } else pending.pop_back();
        }
    }
    std::span<const std::byte> span(const V &value) const {
        const auto &p=parts(value);const auto source=unsigned_value(p.at(0)),offset=unsigned_value(p.at(1)),length=unsigned_value(p.at(2));
        require(sources.contains(source),"native span source unavailable");const auto bytes=sources.at(source).bytes();
        require(offset<=bytes.size() && length<=bytes.size()-offset,"native byte span outside original source");
        return bytes.subspan(machine::checked_size(offset),machine::checked_size(length));
    }
};
NativeRecord::NativeRecord(std::shared_ptr<const State> state,ID type,const V *value,
                          std::shared_ptr<const V> owner)
    :state_(std::move(state)),type_(type),value_(value),value_owner_(std::move(owner)),
      root_view_(value==state_->root) {
    std::set<std::uint64_t> visited;
    std::set<State::PackedAddress> packed;
    while(type_==owned_reference_type || type_==packed_value_type) {
        if (type_ == packed_value_type) {
            require(packed.insert(State::packed_address(*value_)).second,
                    "cyclic packed value reference");
            auto [decoded_type, decoded] = state_->unpack(*value_);
            type_ = decoded_type; value_owner_ = std::move(decoded); value_ = value_owner_.get();
            continue;
        }
        const auto &p=parts(*value_);const auto key=unsigned_value(p.at(0));
        require(visited.insert(key).second && state_->values.contains(key),"cyclic or missing owned value reference");
        const auto &target=*state_->values.at(key);
        require(target.type==ID{unsigned_value(p.at(1)),unsigned_value(p.at(2)),unsigned_value(p.at(3))},"owned value reference type mismatch");
        type_=target.type;value_=&target.value;
    }
}
std::optional<NativeRecord> NativeRecord::parent() const {
    if(!state_->parent_value_key)return {};
    return declared_value(state_->parent_value_key);
}
RawByteSource NativeRecord::byte_source() const {
    const auto constructor=kind();
    require(constructor==NativeRecordKind::Bytes || constructor==NativeRecordKind::Text,"record is not a raw source span");
    const auto &p=parts(*value_);const auto source=state_->sources.at(unsigned_value(p.at(0)));
    require(unsigned_value(p.at(1))==0 && unsigned_value(p.at(2))==source.bytes().size(),"record span is not a complete raw source");
    return source;
}
NativeRecord NativeRecord::declared_value(std::uint64_t key) const {
    const auto found=state_->values.find(key);
    require(found!=state_->values.end(),"unknown native declared value");const auto &value=*found->second;
    return NativeRecord(state_,value.type,&value.value);
}
NativeRecordValue NativeRecord::materialize_value(const std::map<std::uint64_t,std::uint64_t> &source_map) const {return copy_value(source_map,true);}
NativeRecordValue NativeRecord::unbound_value(const std::map<std::uint64_t,std::uint64_t> &source_map) const {return copy_value(source_map,false);}
NativeRecordValue NativeRecord::copy_value(const std::map<std::uint64_t,std::uint64_t> &source_map,bool retain_owner) const {
    using N=NativeRecordValue;
    auto payload=[&]() -> N {
    switch(kind()) {
    case NativeRecordKind::Null:return N{};
    case NativeRecordKind::Boolean:return N(boolean());
    case NativeRecordKind::Integer:return N(integer());
    case NativeRecordKind::Rational:return N(rational());
    case NativeRecordKind::ExactNumber:return N(N::ExactNumber{rational(),approximate_view()});
    case NativeRecordKind::Decimal:return N(N::Decimal{rational(),text()});
    case NativeRecordKind::Text:case NativeRecordKind::Bytes:{
        const auto &p=parts(*value_);const auto key=unsigned_value(p.at(0));
        return N(N::ByteSpan{source_map.empty()?key:source_map.at(key),machine::checked_size(unsigned_value(p.at(1))),machine::checked_size(unsigned_value(p.at(2))),kind()==NativeRecordKind::Text});
    }
    case NativeRecordKind::Object:{
        N::Object out;
        const auto &factors=state_->publication.declarations->registry.at(type_).factors;
        const auto &elements=parts(*value_);
        out.members.reserve(factors.size());
        for(std::size_t i=0;i<factors.size();++i)
            out.members.emplace_back(state_->labels.at(factors[i].key),
                NativeRecord(state_,factors[i].type,&elements.at(i),value_owner_).copy_value(source_map,retain_owner));
        return N(std::move(out));
    }
    case NativeRecordKind::Array:{N::Array out;const auto count=size();out.elements.reserve(count);for(std::size_t i=0;i<count;++i)out.elements.push_back(at(i).copy_value(source_map,retain_owner));return N(std::move(out));}
    case NativeRecordKind::Tensor:return N(tensor());
    }
    throw std::invalid_argument("unknown native record constructor");
    }();
    if(retain_owner)if(const auto found=state_->owners.find(value_);found!=state_->owners.end()) {
        const auto *value=found->second;
        return N(N::Owned{value->block,value->cell,value->occurrence,std::make_shared<const N>(std::move(payload))});
    }
    return payload;
}
std::uint64_t new_native_record_namespace() {
    std::uint64_t key{};
    do {
        auto *p=reinterpret_cast<std::byte *>(&key);std::size_t position{};
        while(position<sizeof key) {
            const auto n=getrandom(p+position,sizeof key-position,0);
            if(n<0 && errno==EINTR)continue;
            if(n<=0)throw std::runtime_error("native namespace allocation failed");
            position+=static_cast<std::size_t>(n);
        }
    }while(!key || key==format_namespace);
    return key;
}
namespace {
struct OperationValue {
    NativeRecordValue value;
    std::set<std::uint64_t> request_sources;
};
OperationValue operation_value(const std::vector<RawByteSource> &sources,
    std::vector<RawByteSource> &raw,const NativeRecord &request,NativeRecordValue result) {
    using N=NativeRecordValue;
    std::set<std::uint64_t> used;for(const auto &s:sources)used.insert(s.key());
    std::set<std::uint64_t> request_sources;
    std::uint64_t key=1;std::map<std::uint64_t,std::uint64_t> bindings;
    const auto next=[&] {
        while(used.contains(key)){require(key!=std::numeric_limits<std::uint64_t>::max(),"operation byte identity exhausted");++key;}
        used.insert(key);request_sources.insert(key);return key;
    };
    for(const auto &s:request.raw_sources()) {
        const auto target=next();bindings.emplace(s.key(),target);raw.push_back(s.rekey(target));
    }
    const auto archive_key=next();auto archive=request.owning_publication().archive_source(archive_key);const auto count=archive.bytes().size();
    raw.push_back(std::move(archive));
    N::Object value;
    value.members.emplace_back("schema", N(std::string("gadget.operation.native.v1")));
    value.members.emplace_back("request", request.unbound_value(bindings));
    value.members.emplace_back("result", std::move(result));
    value.members.emplace_back("request_publication", N(N::ByteSpan{archive_key,0,count,false}));
    return {N(std::move(value)),std::move(request_sources)};
}
NativeRecord checked_operation(NativeRecord publication,
                               const std::set<std::uint64_t> &request_sources) {
    struct Cursor {
        NativeRecord value;
        std::vector<std::string> names;
        std::size_t next{};
    };
    std::vector<Cursor> pending;
    const auto push = [&](NativeRecord value) {
        auto names = value.kind() == NativeRecordKind::Object ? value.names()
                                                            : std::vector<std::string>{};
        pending.push_back({std::move(value),std::move(names),0});
    };
    push(publication.at("result"));
    while (!pending.empty()) {
        auto &cursor = pending.back();
        const auto kind = cursor.value.kind();
        if (kind == NativeRecordKind::Text || kind == NativeRecordKind::Bytes)
            require(!request_sources.contains(unsigned_value(parts(cursor.value.value()).at(0))),
                    "operation result span has no original or added byte source");
        if (kind == NativeRecordKind::Object && cursor.next < cursor.names.size()) {
            auto child = cursor.value.at(cursor.names[cursor.next++]);
            push(std::move(child));
        } else if (kind == NativeRecordKind::Array && cursor.next < cursor.value.size()) {
            auto child = cursor.value.at(cursor.next++);
            push(std::move(child));
        } else pending.pop_back();
    }
    return publication;
}
} // namespace
NativeRecord retain_operation(const NativeRecord &source,const NativeRecord &request,NativeRecordValue result,
                              const NativeRecordObserver &observer) {
    std::vector<RawByteSource> raw;
    auto operation=operation_value(source.raw_sources(),raw,request,std::move(result));
    return checked_operation(source.derive(std::move(operation.value),std::move(raw),observer),
                             operation.request_sources);
}
NativeRecord retain_operation(const NativeRecord &source, std::vector<RawByteSource> added,
    const NativeRecord &request, NativeRecordValue result,
    const NativeRecordObserver &observer) {
    auto sources = source.raw_sources();
    std::set<std::uint64_t> keys;
    for (const auto &bytes : sources) keys.insert(bytes.key());
    for (const auto &bytes : added) {
        require(keys.insert(bytes.key()).second, "operation added byte source identity collides");
        sources.push_back(bytes);
    }
    auto operation = operation_value(sources, added, request, std::move(result));
    return checked_operation(source.derive(std::move(operation.value), std::move(added), observer),
                             operation.request_sources);
}
NativeRecord retain_operation(const storage::TemporalEndpoint &source,std::uint64_t source_namespace,
    std::vector<RawByteSource> sources,const NativeRecord &request,NativeRecordValue result) {
    return retain_operation(storage::TemporalEndpoint(source),source_namespace,
                            std::move(sources),request,std::move(result));
}
NativeRecord retain_operation(storage::TemporalEndpoint &&source,std::uint64_t source_namespace,
    std::vector<RawByteSource> sources,const NativeRecord &request,NativeRecordValue result) {
    std::vector<RawByteSource> added;
    auto operation=operation_value(sources,added,request,std::move(result));
    for(auto &bytes:added)sources.push_back(std::move(bytes));
    return checked_operation(NativeRecord::create_from_endpoint(std::move(operation.value),
        source_namespace,std::move(sources),std::move(source)),operation.request_sources);
}
NativeRecord NativeRecord::create(const NativeRecordValue &value,std::uint64_t ns,
    std::vector<RawByteSource> sources,const storage::TemporalEndpoint *source) {
    auto endpoint=source?*source:storage::TemporalEndpoint([] {
        core::RelationalComplexBuilder builder(core::RelationSignature{});builder.declare_witness_identity(1);return builder.freeze();
    }());
    return create_from_endpoint(value,ns,std::move(sources),std::move(endpoint));
}
NativeRecord NativeRecord::create_from_endpoint(const NativeRecordValue &value,std::uint64_t ns,
    std::vector<RawByteSource> sources,storage::TemporalEndpoint endpoint) {
    return create_owned(value,ns,std::move(sources),std::move(endpoint),nullptr);
}
NativeRecord NativeRecord::create_from_endpoint(NativeRecordValue &&value,std::uint64_t ns,
    std::vector<RawByteSource> sources,storage::TemporalEndpoint endpoint,const NativeRecordObserver &observer) {
    return create_owned(std::move(value),ns,std::move(sources),std::move(endpoint),nullptr,{},observer);
}
NativeRecord NativeRecord::publish_fields(const NativeRecord &source, std::uint64_t ns,
                                          std::vector<storage::TemporalTensor> fields,
                                          NativeRecordValue value,
                                          std::vector<RawByteSource> added_sources,
                                          std::vector<storage::TemporalCarrier> added_carriers) {
    require(ns != source.source_namespace(), "field publication requires a new namespace");
    auto endpoint = *source.endpoint();
    endpoint.carriers.insert(endpoint.carriers.end(), std::make_move_iterator(added_carriers.begin()),
                             std::make_move_iterator(added_carriers.end()));
    endpoint.tensors.insert(endpoint.tensors.end(), std::make_move_iterator(fields.begin()),
                            std::make_move_iterator(fields.end()));
    return create_owned(std::move(value), ns, std::move(added_sources), std::move(endpoint),
                        &source);
}
NativeRecord NativeRecord::create_owned(NativeRecordValue value,std::uint64_t ns,
    std::vector<RawByteSource> sources,storage::TemporalEndpoint endpoint,
    const NativeRecord *retained_parent,std::optional<storage::NativeTemporalSeries> history,
    const NativeRecordObserver &observer) {
    require(ns,"native record requires its actual source namespace");
    std::optional<NativeRecord> previous;
    if(retained_parent)previous=*retained_parent;
    else if(endpoint.declarations)for(const auto &v:endpoint.declarations->values) if(v.key==header_key) {
        require(v.type==header_type,"source value collides with native record header identity");
        previous=history?from_history(*history):from_endpoint(endpoint,ns);
    }
    auto declarations=endpoint.declarations?std::move(*endpoint.declarations):storage::TemporalDeclarationState{};
    declarations.resource_limits=native_data_limits().temporal.declarations;
    if(!declarations.namespace_key)declarations.namespace_key=ns;
    if(!previous)for(const auto &v:declarations.values)
        require(v.key!=root_key && v.key!=labels_key && v.key!=sources_key,
            "source value collides with native record control identity");
    const auto additional_sources=sources;
    if(previous) {
        sources=previous->raw_sources();
        for(const auto &raw:additional_sources) {
            require(std::none_of(sources.begin(),sources.end(),[&](const auto &old){return old.key()==raw.key();}),
                    "derived record requires independent added raw source identities");
            sources.push_back(raw);
        }
    }
    for(auto &bytes:sources)bytes=bytes.freeze();
    Compiler compiler(new_native_record_namespace(),sources);
    for(const auto &v:declarations.values)compiler.value_keys.insert(v.key);
    for(const auto key:{header_key,root_key,labels_key,sources_key})compiler.value_keys.insert(key);
    if(previous) {
        compiler.name_keys=previous->state_->names;
        for(const auto &v:declarations.values)if(v.key==labels_key)
            compiler.labels=std::get<V::Repeated>(v.value.data).elements;
    }
    if(observer)observer("source_declarations");
    auto [type,data]=compiler.compile(std::move(value));
    if(observer)observer("composite_values");
    sources.push_back(RawByteSource::adopt_fresh(compiler.pool,std::move(compiler.text_bytes)));
    auto all_definitions=declarations.registry.definitions();
    std::set<ID> ids;for(const auto &d:all_definitions)ids.insert(d.id);
    for(auto &d:compiler.defs) {
        if(ids.contains(d.id))require(declarations.registry.at(d.id)==d,"immutable application declaration redefined");
        else {
            ids.insert(d.id);
            all_definitions.push_back(std::move(d));
        }
    }
    declarations.registry =
        composite::DeclarationRegistry(std::move(all_definitions), declarations.resource_limits);
    std::vector<V> raw;
    require(!history || history->last_sequence() != std::numeric_limits<std::uint64_t>::max(),
            "native annotation publication sequence exhausted");
    const auto current_sequence = history ? history->last_sequence() + 1 : endpoint.sequence;
    for (const auto &source_bytes : sources) {
        const auto bytes = source_bytes.bytes();
        const RawByteSource *retained = nullptr;
        std::optional<State::SourceOrigin> origin;
        if (previous && history)
            for (const auto &[key, old] : previous->state_->sources)
                if (source_bytes.shares_owner(old)) {
                    retained = &old;
                    origin = previous->state_->source_payloads.at(key).origin;
                    break;
                }
        V payload;
        if (retained) {
            require(origin->sequence < current_sequence,
                    "byte source origin must precede its reference");
            payload = V::alternative(
                2, V::product({number(ns), number(origin->sequence), number(origin->key)}));
        } else {
            numeric::Integer integer;
            if (!bytes.empty())
                mpz_import(integer.get_mpz_t(), bytes.size(), -1, 1, 0, 0, bytes.data());
            payload = V::alternative(1, V::integer(std::move(integer)));
        }
        raw.push_back(
            V::product({number(source_bytes.key()), number(bytes.size()), std::move(payload)}));
    }
    for (auto &v : compiler.owned_values)
        declarations.values.push_back(std::move(v));
    std::uint64_t data_key = root_key;
    if (previous) {
        std::set<std::uint64_t> used;
        for (const auto &v : declarations.values)
            used.insert(v.key);
        while (used.contains(data_key)) {
            require(data_key != std::numeric_limits<std::uint64_t>::max(),
                    "native data identity exhausted");
            ++data_key;
        }
    }
    const auto replace = [&](std::uint64_t key, ID declared_type, V v) {
        auto found = std::find_if(declarations.values.begin(), declarations.values.end(),
                                  [&](const auto &x) { return x.key == key; });
        storage::TemporalDeclaredValue next{key,declared_type,0,0,{},{},std::move(v)};
        if(found==declarations.values.end())
            declarations.values.push_back(std::move(next));
        else
            *found = std::move(next);
    };
    const auto version = compiler.packed_values || (previous && previous->format_version() >= 4) ? 4 : 3;
    replace(header_key, header_type,
            V::product({number(version), number(ns), number(type.name_space), number(type.key),
                        number(type.version), number(data_key),
                        number(previous && previous->source_namespace() == ns
                                   ? previous->state_->root_value_key
                                   : 0)}));
    replace(data_key,type,std::move(data));
    replace(labels_key,labels_type,V::repeated(std::move(compiler.labels)));
    replace(sources_key,referenced_sources_type,V::repeated(std::move(raw)));
    endpoint.declarations=std::move(declarations);
    if(observer)observer("declaration_channels");
    if(history) {
        require(history->last_sequence()!=std::numeric_limits<std::uint64_t>::max(),
                "native annotation publication sequence exhausted");
        endpoint.sequence=history->last_sequence()+1;
        const auto base=history->tail();
        auto maps=base.endpoint().declarations
            ?composite::identity_declaration_maps(base.endpoint().declarations->registry)
            :std::vector<composite::DeclarationMap>{};
        const auto target_hash=storage::temporal_endpoint_hash(endpoint);
        auto update=storage::make_temporal_declaration_update(base.endpoint().declarations,
            std::move(endpoint.declarations),std::move(maps));
        // This operation changes declarations only. Publish that operation
        // directly; a general target difference would compute field transport
        // and innovations for numerical channels that did not change.
        storage::TemporalSegment annotation;
        annotation.base_sequence = base.endpoint().sequence;
        annotation.target_sequence = endpoint.sequence;
        annotation.base_source_hash = storage::encode_rfcx(base.endpoint().source).stats.payload_hash;
        annotation.target_source_hash = annotation.base_source_hash;
        annotation.base_endpoint_hash = storage::temporal_endpoint_hash(base.endpoint());
        annotation.target_endpoint_hash = target_hash;
        annotation.source.target_grade_count =
            static_cast<std::uint32_t>(base.endpoint().source.grade_count());
        annotation.declarations = std::move(update);
        history->append(annotation);
        const auto tail=history->tail();
        if(observer)observer("native_freeze");
        auto state=std::make_shared<const State>(tail,std::move(sources),std::move(*history));
        if(observer)observer("record_bindings");
        return {state,state->root_type,state->root};
    }
    storage::NativeTemporalSeries publication(endpoint,ns,native_data_limits());
    const auto tail=publication.tail();
    if(observer)observer("native_freeze");
    auto state=std::make_shared<const State>(tail,std::move(sources),std::move(publication));
    if(observer)observer("record_bindings");
    return {state,state->root_type,state->root};
}
NativeRecord NativeRecord::annotate(const NativeRecordValue &value,storage::NativeTemporalSeries history,
                                   std::vector<RawByteSource> sources) {
    const auto tail=history.tail();
    const auto ns=history.source_namespace();
    return create_owned(value,ns,std::move(sources),tail.endpoint(),nullptr,std::move(history));
}
NativeRecord NativeRecord::annotate(const NativeRecordValue &value,const NativeRecord &producer,
                                   std::vector<RawByteSource> sources) {
    require(producer.is_root() && producer.state_->history.has_value(),"annotation requires an owning producer history");
    return create_owned(value,producer.source_namespace(),std::move(sources),*producer.endpoint(),
        &producer,*producer.state_->history);
}
NativeRecord NativeRecord::derive(const NativeRecordValue &value,std::vector<RawByteSource> sources) const {
    const auto source=endpoint();
    auto history=state_->history && state_->history->segment_count()!=0?state_->history:std::nullopt;
    return create_owned(value,source_namespace(),std::move(sources),*source,this,std::move(history));
}
NativeRecord NativeRecord::derive(NativeRecordValue &&value,std::vector<RawByteSource> sources,
                                 const NativeRecordObserver &observer) const {
    const auto source=endpoint();
    auto history=state_->history && state_->history->segment_count()!=0?state_->history:std::nullopt;
    return create_owned(std::move(value),source_namespace(),std::move(sources),*source,this,std::move(history),observer);
}
std::vector<std::uint64_t> NativeRecord::source_keys(const storage::TemporalEndpoint &endpoint) {
    std::vector<std::uint64_t> keys;
    if(!endpoint.declarations)return keys;
    for(const auto &record:endpoint.declarations->values)if(record.key==sources_key) {
        require(record.type==sources_type || record.type==referenced_sources_type,"retained byte source declaration mismatch");
        for(const auto &entry:std::get<V::Repeated>(record.value.data).elements)
            keys.push_back(unsigned_value(parts(entry).at(0)));
    }
    return keys;
}
NativeRecord NativeRecord::from_endpoint(const storage::TemporalEndpoint &endpoint,std::uint64_t ns) {
    return from_history(storage::NativeTemporalSeries(endpoint,ns,native_data_limits()));
}
NativeRecord NativeRecord::from_history(storage::NativeTemporalSeries history) {
    const auto tail=history.tail();
    auto state=std::make_shared<const State>(tail,std::vector<RawByteSource>{},std::move(history));
    return {state,state->root_type,state->root};
}
std::optional<storage::NativeTemporalSeries> NativeRecord::publication_history() const {return state_->history;}
std::optional<storage::NativeEncodedTensor> NativeRecord::encoded_tensor(std::uint64_t key) const {
    if(!state_->history)return std::nullopt;
    return storage::NativeEncodedTensor::from_series(*state_->history,key,endpoint()->sequence);
}
storage::PinnedNativeAccession NativeRecord::access(
    const storage::NativeAccessionRequest &request) const {
    return state_->study.access(request);
}
NativeRecord NativeRecord::read(const std::filesystem::path &path) {
    auto limits=native_data_limits();
    if(const auto room=runtime::detect_hardware().advisory_memory_headroom())limits.temporal.declarations.max_materialized_bytes=std::max<std::size_t>(1,*room);
    const storage::ReadOnlyMappedFile file(path);
    return decode(file.bytes(),limits);
}
std::optional<NativeRecord> NativeRecord::read_if_record(const std::filesystem::path &path) {
    auto limits=native_data_limits();
    if(const auto room=runtime::detect_hardware().advisory_memory_headroom())
        limits.temporal.declarations.max_materialized_bytes=std::max<std::size_t>(1,*room);
    const storage::ReadOnlyMappedFile file(path);
    const auto bytes=file.bytes();
    if(detect_native_data(bytes)!=NativeDataFormat::TemporalSeries)return {};
    const auto series=storage::decode_native_series(bytes,limits);
    const auto publication=series.tail();
    if(!publication.endpoint().declarations)return {};
    for(const auto &value:publication.endpoint().declarations->values)
        if(value.key==header_key && value.type==header_type)
            return from_history(series);
    return {};
}
NativeRecord NativeRecord::decode(std::span<const std::byte> bytes,storage::NativeSeriesLimits limits) {
    const auto series=storage::decode_native_series(bytes,limits);
    return from_history(series);
}
std::vector<std::byte> NativeRecord::encode() const {
    require(is_root(),"persist an owning record, not a detached member view");
    return storage::encode_native_series(state_->history?*state_->history:
        storage::NativeTemporalSeries(*endpoint(),state_->ns,native_data_limits()));
}
RawByteSource NativeRecord::archive_source(std::uint64_t source_key) const {
    return RawByteSource::adopt_fresh(source_key,encode());
}
void NativeRecord::write(const std::filesystem::path &path) const {
    require(is_root(),"persist an owning record, not a detached member view");
    const auto series=state_->history?*state_->history:
        storage::NativeTemporalSeries(*endpoint(),state_->ns,native_data_limits());
    storage::publish_file_exact(path,storage::FilePublication::CreateNew,[&](std::ostream &output) {
        storage::stream_native_series(series,[&](auto chunk) { storage::write_stream_exact(output,chunk); });
    });
}
bool NativeRecord::is_root() const noexcept {return root_view_;}
std::uint64_t NativeRecord::format_version() const noexcept {return state_->version;}
NativeRecord NativeRecord::owning_publication() const {return {state_,state_->root_type,state_->root};}
NativeRecord NativeRecord::historical_publication(std::uint64_t cutoff) const {
    require(state_->history.has_value(),"record has no immutable history");
    if(cutoff==state_->history->last_sequence())return owning_publication();
    auto history=storage::NativeSeriesPrefix::capture(*state_->history,cutoff).history();
    const auto tail=history.tail();
    auto state=std::make_shared<const State>(tail,std::vector<RawByteSource>{},std::move(history),state_.get());
    return {state,state->root_type,state->root};
}
NativeRecordKind NativeRecord::kind() const {
    if(type_==null_type)return NativeRecordKind::Null;
    if(type_==integer_type)return NativeRecordKind::Integer;
    if(type_==boolean_type)return NativeRecordKind::Boolean;
    if(type_==rational_type)return NativeRecordKind::Rational;
    if(type_==text_type)return NativeRecordKind::Text;
    if(type_==bytes_type)return NativeRecordKind::Bytes;
    if(type_==decimal_type)return NativeRecordKind::Decimal;
    if(type_==exact_number_type)return NativeRecordKind::ExactNumber;
    const auto &d=state_->publication.declarations->registry.at(type_);
    if(d.kind==composite::DeclarationKind::Product)return NativeRecordKind::Object;
    if(d.kind==composite::DeclarationKind::Repeated)return NativeRecordKind::Array;
    if(d.kind==composite::DeclarationKind::RationalTensor)return NativeRecordKind::Tensor;
    throw std::invalid_argument("unsupported native record view constructor");
}
std::size_t NativeRecord::size() const {
    const auto constructor=kind();
    if(constructor==NativeRecordKind::Object)return parts(*value_).size();
    if(constructor==NativeRecordKind::Array)return std::get<V::Repeated>(value_->data).elements.size();
    throw std::invalid_argument("native record size requires a container");
}
bool NativeRecord::has(std::string_view name) const {
    if(kind()!=NativeRecordKind::Object)return false;
    const auto found=state_->names.find(name);if(found==state_->names.end())return false;
    try {
        (void)state_->publication.declarations->registry.factor_ordinal(type_,found->second);
        return true;
    } catch(const std::out_of_range &) {return false;}
}
NativeRecord NativeRecord::at(std::string_view name) const {
    require(kind()==NativeRecordKind::Object,"native member selection requires a product");
    const auto found=state_->names.find(name);
    if(found==state_->names.end())throw std::out_of_range("unknown native factor name");
    const auto &registry=state_->publication.declarations->registry;
    const auto &d=registry.at(type_);
    const auto index=registry.factor_ordinal(type_,found->second);
    return {state_,d.factors.at(index).type,&parts(*value_).at(index),value_owner_};
}
NativeRecord NativeRecord::at(std::size_t index) const {
    require(kind()==NativeRecordKind::Array,"native index selection requires a repeated declaration");
    const auto &registry=state_->publication.declarations->registry;
    auto type=registry.at(type_).element;
    const auto *value=&std::get<V::Repeated>(value_->data).elements.at(index);
    const auto &d=registry.at(type);
    if(d.kind==composite::DeclarationKind::Alternative) {
        const auto &a=std::get<V::Alternative>(value->data);type=d.factors.at(registry.alternative_ordinal(type,a.key)).type;value=&a.element.at(0);
    }
    return {state_,type,value,value_owner_};
}
std::optional<NativeRecord::PackedEncoding> NativeRecord::packed_encoding(std::size_t index) const {
    require(kind() == NativeRecordKind::Array, "packed encoding requires an array");
    const auto &value = std::get<V::Repeated>(value_->data).elements.at(index);
    if (state_->publication.declarations->registry.at(type_).element != packed_value_type)
        return std::nullopt;
    const auto &p = parts(value);
    return PackedEncoding{{unsigned_value(p.at(3)), unsigned_value(p.at(4)),
                           unsigned_value(p.at(5))},
                          machine::checked_size(unsigned_value(p.at(6))), state_->span(value)};
}
std::vector<std::string> NativeRecord::names() const {
    require(kind()==NativeRecordKind::Object,"native factor names require a product");
    const auto &factors=state_->publication.declarations->registry.at(type_).factors;
    std::vector<std::string> result;result.reserve(factors.size());
    for(const auto &factor:factors)result.push_back(state_->labels.at(factor.key));
    return result;
}
bool NativeRecord::boolean() const {require(kind()==NativeRecordKind::Boolean,"native value is not Boolean");return std::get<bool>(value_->data);}
numeric::Integer NativeRecord::integer() const {require(kind()==NativeRecordKind::Integer,"native value is not Integer");return std::get<numeric::Integer>(value_->data);}
numeric::Rational NativeRecord::rational() const {
    const auto constructor=kind();
    require(constructor==NativeRecordKind::Rational || constructor==NativeRecordKind::Decimal || constructor==NativeRecordKind::ExactNumber,"native value is not Rational");
    return std::get<numeric::Rational>(constructor==NativeRecordKind::Rational?value_->data:parts(*value_).at(0).data);
}
bool NativeRecord::approximate_view() const {require(kind()==NativeRecordKind::ExactNumber,"native value has no approximate presentation flag");return std::get<bool>(parts(*value_).at(1).data);}
std::span<const std::byte> NativeRecord::bytes() const {
    const auto constructor=kind();
    require(constructor==NativeRecordKind::Text || constructor==NativeRecordKind::Decimal || constructor==NativeRecordKind::Bytes,"native value has no byte span");
    return state_->span(constructor==NativeRecordKind::Decimal?parts(*value_).at(1):*value_);
}
std::string NativeRecord::text() const {
    const auto constructor=kind();
    require(constructor==NativeRecordKind::Text || constructor==NativeRecordKind::Decimal,"opaque bytes require bytes(), not a text view");
    const auto b=state_->span(constructor==NativeRecordKind::Decimal?parts(*value_).at(1):*value_);
    if(b.empty())return {};
    return {reinterpret_cast<const char *>(b.data()),b.size()};
}
const numeric::SparseRationalTensor &NativeRecord::tensor() const {require(kind()==NativeRecordKind::Tensor,"native value is not a tensor");return std::get<numeric::SparseRationalTensor>(value_->data);}
std::vector<RawByteSource> NativeRecord::raw_sources() const {std::vector<RawByteSource> result;for(const auto &[key,s]:state_->sources){(void)key;result.push_back(s);}return result;}
std::shared_ptr<const storage::TemporalEndpoint> NativeRecord::endpoint() const {return state_->study.endpoint();}
std::shared_ptr<const storage::NativeAccessionDirectory> NativeRecord::accession_directory() const {return state_->study.accession_directory();}
std::uint64_t NativeRecord::source_namespace() const {return state_->ns;}
} // namespace gadget
