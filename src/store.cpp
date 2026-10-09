#include <charconv>
#include <gadget/store.hpp>
#include <limits>
#include <rfc/storage/temporal_seed.hpp>

namespace gadget {
namespace {
namespace native = rfc::storage;
// Canonical descriptor and sync helpers are confined to this translation unit.
namespace durable_io = native::linux_durable_detail;
constexpr std::uint64_t study_namespace = 0x4741444745540001ULL;
constexpr std::size_t artifact_limit = std::numeric_limits<std::size_t>::max();
const std::array<std::string, 3> artifact_names{"snapshot.json", "metadata.json", "study.rfcxs"};
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::invalid_argument("study store: " + message);
}
void members(const Json &value, std::initializer_list<const char *> names) {
    require(json_object_is_type(value.raw(), json_type_object), "expected object");
    std::set<std::string> expected;
    for (const auto *name : names) {
        expected.insert(name);
        require(value.has(name), "missing member: " + std::string(name));
    }
    json_object_object_foreach(value.raw(), key, item) {
        (void)item;
        require(expected.contains(key), "unknown member: " + std::string(key));
    }
}
std::uint64_t decimal(const std::string &value) {
    require(!value.empty() && value.size() <= 20 && (value.size() == 1 || value.front() != '0'),
            "expected canonical unsigned decimal string");
    std::uint64_t result{};
    auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    require(error == std::errc{} && end == value.data() + value.size(), "invalid unsigned integer");
    return result;
}
std::string hex(std::uint64_t value) {
    constexpr char digits[] = "0123456789abcdef";
    std::string result(16, '0');
    for (std::size_t i = 0; i < 16; ++i)
        result[15 - i] = digits[(value >> (4 * i)) & 15];
    return result;
}
std::uint64_t unhex(const std::string &value) {
    require(value.size() == 16 && value.find_first_not_of("0123456789abcdef") == std::string::npos,
            "invalid lowercase hexadecimal identity");
    std::uint64_t result{};
    auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result, 16);
    require(error == std::errc{} && end == value.data() + value.size(),
            "invalid hexadecimal identity");
    return result;
}
Json root_json(const native::DurableNativeRoot &root) {
    auto value = Json::object();
    value.set("commit", Json::string(hex(root.commit.high) + hex(root.commit.low)));
    value.set("source_namespace", Json::string(std::to_string(root.source_namespace)));
    value.set("generation", Json::string(std::to_string(root.generation)));
    value.set("sequence", Json::string(std::to_string(root.sequence)));
    return value;
}
native::DurableNativeRoot parse_root(const Json &value) {
    members(value, {"commit", "source_namespace", "generation", "sequence"});
    auto id = value.at("commit").text();
    require(id.size() == 32, "invalid commit identity");
    native::DurableNativeRoot root{{unhex(id.substr(0, 16)), unhex(id.substr(16))},
                                   decimal(value.at("source_namespace").text()),
                                   decimal(value.at("generation").text()),
                                   decimal(value.at("sequence").text())};
    require(!root.commit.empty() && root.source_namespace == study_namespace &&
                root.generation == root.sequence,
            "unexpected predecessor profile");
    return root;
}
native::DurableNativeLimits limits() {
    native::DurableNativeLimits result;
    result.series=native_data_limits();
    result.max_object_bytes=std::numeric_limits<std::size_t>::max();
    result.max_commits=std::numeric_limits<std::size_t>::max();
    return result;
}
fs::path store_path(const fs::path &requested) {
    auto descriptor = durable_io::directory(requested); // Reject a symlink at the store root.
    return fs::canonical(requested);
}
std::string string_bytes(const std::vector<std::byte> &value) {
    return {reinterpret_cast<const char *>(value.data()), value.size()};
}
std::span<const std::byte> bytes(const std::string &value) {
    return {reinterpret_cast<const std::byte *>(value.data()), value.size()};
}
struct Frozen {
    Workspace workspace;
    std::map<std::string, std::string> artifacts;
};
Frozen freeze(const fs::path &path) {
    Frozen result{load_snapshot(path), {}};
    auto directory = durable_io::directory(path);
    for (const auto &name : artifact_names)
        if(fs::exists(path/name))result.artifacts.emplace(
            name,
            string_bytes(durable_io::read(directory.get(), name,
                                          artifact_limit)));
    const auto &origin = result.workspace.snapshot_origin;
    if(result.artifacts.contains("snapshot.json"))require("sha256:"+sha256(result.artifacts.at("snapshot.json"))==origin.at("completion_sha256").text(),"snapshot changed during preparation");
    for(const auto &[name,data]:result.artifacts)if(name!="snapshot.json") {
        const auto record=origin.at("artifacts").at(name);
        require("sha256:"+sha256(data)==record.at("sha256").text() && data.size()==static_cast<std::size_t>(record.at("bytes").integer()),"snapshot artifact changed during preparation");
    }
    return result;
}
void private_directory(const fs::path &path) {
    require(fs::create_directory(path), "destination already exists: " + path.string());
    fs::permissions(path, fs::perms::owner_all, fs::perm_options::replace);
}
// Prepare immutable provenance before calling canonical publication. There is no
// application HEAD or journal. Interrupted random stages remain unreachable.
void prepare(const fs::path &root, std::uint64_t transaction, const Frozen &frozen,
             const StudyPreparationHook &hook) {
    auto notify = [&](StudyPreparationStage stage) {
        if (hook)
            hook(stage);
    };
    const auto studies = root / "studies";
    auto parent = durable_io::directory(studies);
    const auto final = studies / hex(transaction);
    if (fs::exists(fs::symlink_status(final))) {
        auto directory = durable_io::directory(final);
        for (const auto &[name,artifact_data] : frozen.artifacts) {
            require(string_bytes(durable_io::read(directory.get(), name, artifact_limit)) ==
                        artifact_data,
                    "transaction already prepared with different snapshot bytes");
            auto file = durable_io::regular(directory.get(), name, O_RDONLY);
            durable_io::synchronize(file.get());
            file.close_checked();
        }
        durable_io::synchronize(directory.get());
    } else {
        const auto random = durable_io::random_name();
        const auto stage = studies / (".stage-" + hex(random[0]) + hex(random[1]));
        private_directory(stage);
        notify(StudyPreparationStage::DirectoryCreated);
        auto directory = durable_io::directory(stage);
        for (const auto &[name,artifact_data] : frozen.artifacts) {
            auto file = durable_io::regular(directory.get(), name, O_WRONLY | O_CREAT | O_EXCL);
            durable_io::write(file.get(), bytes(artifact_data));
            notify(StudyPreparationStage::FileWritten);
            durable_io::synchronize(file.get());
            notify(StudyPreparationStage::FileSynced);
            file.close_checked();
        }
        durable_io::synchronize(directory.get());
        notify(StudyPreparationStage::DirectorySynced);
        fs::rename(stage, final); // Exclusive native coordinator protects cooperating publishers.
        notify(StudyPreparationStage::Renamed);
    }
    durable_io::synchronize(parent.get());
    auto directory = durable_io::directory(root);
    durable_io::synchronize(directory.get());
    notify(StudyPreparationStage::ParentSynced);
}
struct Binding {
    std::uint64_t transaction;
    std::string digest;
    NativeRecord publication;
};
Binding binding(const native::TemporalEndpoint &endpoint) {
    const auto publication=NativeRecord::from_endpoint(endpoint,study_namespace);
    require(publication.at("schema").text()=="gadget.store.publication.native.v1","invalid durable publication declaration");
    const auto transaction=publication.at("transaction").integer();
    require(transaction>=0 && transaction<=rfc::numeric::Integer(std::to_string(std::numeric_limits<std::uint64_t>::max())),"transaction outside native identity domain");
    const auto bytes=publication.at("snapshot").bytes();
    return {decimal(transaction.get_str()),"sha256:"+sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size())),publication};
}
using Keys = std::map<std::string, std::uint64_t>;
Keys tensor_keys(const native::TemporalEndpoint &endpoint) {
    Keys result;
    for (const auto &tensor : endpoint.tensors)
            require(tensor.key > 0 &&
                        result.emplace(tensor.interpretation, tensor.key).second,
                    "invalid study tensor identity");
    return result;
}
native::TemporalEndpoint target(const Workspace &workspace,
                                const native::TemporalEndpoint *previous,
                                std::uint64_t transaction) {
    auto endpoint = *workspace.native_snapshot->endpoint();
    auto keys = previous ? tensor_keys(*previous) : tensor_keys(endpoint);
    std::uint64_t next = 1;
    for (const auto &[name, key] : keys) {
        (void)name;
        require(key!=std::numeric_limits<std::uint64_t>::max(),"native tensor identity space exhausted");
        next = std::max(next, key + 1);
    }
    for (auto &tensor : endpoint.tensors) {
        auto found = keys.find(tensor.interpretation);
        if (found == keys.end()) {
            require(next != std::numeric_limits<std::uint64_t>::max(), "native tensor identity space exhausted");
            found = keys.emplace(tensor.interpretation, next++).first;
        }
        tensor.key = found->second;
    }
    if (previous) {
        require(previous->sequence != std::numeric_limits<std::uint64_t>::max(), "native sequence identity space exhausted");
        endpoint.sequence = previous->sequence + 1;
    }
    const auto source=workspace.native_publication.value();
    auto raw=source.raw_sources();std::map<std::uint64_t,std::uint64_t> source_map;
    if(previous) {
        const auto retained=NativeRecord::from_endpoint(*previous,study_namespace).raw_sources();
        std::set<std::uint64_t> used;for(const auto &r:retained)used.insert(r.key());
        std::uint64_t raw_key=1;
        for(auto &r:raw) {
            while(used.contains(raw_key)){require(raw_key!=std::numeric_limits<std::uint64_t>::max(),"raw source identity exhausted");++raw_key;}
            source_map.emplace(r.key(),raw_key);r=RawByteSource(raw_key,r.bytes());used.insert(raw_key);
        }
        endpoint.declarations=previous->declarations;
    } else endpoint.declarations.reset();
    std::set<std::uint64_t> raw_ids;for(const auto &r:raw)raw_ids.insert(r.key());
    if(previous)for(const auto &r:NativeRecord::from_endpoint(*previous,study_namespace).raw_sources())raw_ids.insert(r.key());
    std::uint64_t archive_key=1;while(raw_ids.contains(archive_key)){require(archive_key!=std::numeric_limits<std::uint64_t>::max(),"snapshot byte identity exhausted");++archive_key;}
    auto archive=source.encode();const auto archive_count=archive.size();raw.push_back(RawByteSource::own(archive_key,std::move(archive)));
    using N=NativeRecordValue;
    const auto publication=NativeRecord::create(N(N::Object{{
        {"schema",N(std::string("gadget.store.publication.native.v1"))},
        {"transaction",N(rfc::numeric::Integer(std::to_string(transaction)))},
        {"snapshot",N(N::ByteSpan{archive_key,0,archive_count,false})},
        {"source",source.materialize_value(source_map)}}}),study_namespace,std::move(raw),&endpoint);
    endpoint=*publication.endpoint();
    native::validate_temporal_endpoint(endpoint);
    return endpoint;
}
struct Admission {
    std::vector<std::byte> source;
    std::vector<std::string> genes, assertions;
    std::map<std::string, std::vector<bool>> contexts;
    std::map<std::string, RetainedRecord> meanings;
};
Admission admission(const Workspace &workspace) {
    Admission result;
    result.source = native::encode_rfcx(*workspace.complex).bytes;
    result.genes = workspace.genes;
    for (const auto &assertion : workspace.assertions)
        result.assertions.push_back(assertion.id);
    for (const auto &[name, context] : workspace.contexts) {
        auto &bits = result.contexts[name];
        for (std::uint16_t grade = 0; grade < workspace.complex->grade_count(); ++grade) {
            const auto &mask = context->grade(rfc::core::Grade{grade});
            for (std::size_t i = 0; i < mask.size(); ++i)
                bits.push_back(mask.get(i));
        }
    }
    for (const auto &[name, meaning] : workspace.field_semantics)
        result.meanings.emplace(name, meaning);
    return result;
}
void admit(const Admission &before, const Admission &after) {
    require(before.source == after.source && before.genes == after.genes &&
                before.assertions == after.assertions && before.contexts == after.contexts,
            "fixed study geometry, biomedical identities or contexts changed");
    for (const auto &[name, meaning] : before.meanings)
        require(after.meanings.contains(name) &&
                    after.meanings.at(name).dump() == meaning.dump(),
                "existing field removed or meaning changed: " + name);
}
void match(const Workspace &workspace, const native::TemporalEndpoint &endpoint,const NativeRecord &durable) {
    const auto portable = workspace.native_snapshot->endpoint();
    require(native::encode_rfcx(endpoint.source).bytes ==
                native::encode_rfcx(portable->source).bytes,
            "durable primary source disagrees with portable study");
    const auto keys = tensor_keys(*portable);
    native::TemporalEndpoint normalized(endpoint.source,0);
    normalized.local_blocks=endpoint.local_blocks;normalized.carriers=endpoint.carriers;
    normalized.tensors=endpoint.tensors;normalized.finite_systems=endpoint.finite_systems;
    normalized.codes=endpoint.codes;
    require(normalized.tensors.size() == portable->tensors.size(), "study tensor count mismatch");
    std::set<std::string> seen;
    for (auto &tensor : normalized.tensors) {
        require(keys.contains(tensor.interpretation) && seen.insert(tensor.interpretation).second,
                "study tensor name mismatch");
        tensor.key = keys.at(tensor.interpretation);
    }
    const auto &source=workspace.native_publication.value();
    require(durable.at("schema").text()=="gadget.store.publication.native.v1" && durable.at("source").at("study").kind()==NativeRecordKind::Object,"invalid native durable source binding");
    require(render_record_view(durable.at("source").at("study")).dump()==render_record_view(source.at("study")).dump(),"durable biomedical declarations differ from admitted source");
    const auto originals=source.at("original_files"),retained=durable.at("source").at("original_files");
    require(originals.size()==retained.size(),"durable original byte map is incomplete");
    for(std::size_t i=0;i<originals.size();++i) {
        const auto a=originals.at(i).at("bytes").bytes(),b=retained.at(i).at("bytes").bytes();
        require(originals.at(i).at("path").text()==retained.at(i).at("path").text() && std::equal(a.begin(),a.end(),b.begin(),b.end()),"durable original source bytes differ");
    }
    normalized.declarations.reset();auto portable_math=*portable;portable_math.declarations.reset();
    require(native::encode_temporal_endpoint(normalized) ==
                native::encode_temporal_endpoint(portable_math),
            "durable endpoint disagrees with portable study");
}
void match(const Workspace &workspace,const native::TemporalEndpoint &endpoint) {
    match(workspace,endpoint,NativeRecord::from_endpoint(endpoint,study_namespace));
}
Workspace endpoint_workspace(const fs::path &root,const native::TemporalEndpoint &endpoint,const Binding &reference) {
    require((endpoint.sequence == 0) == (reference.transaction == 0), "invalid genesis binding");
    auto workspace = load_snapshot(root / "studies" / hex(reference.transaction));
    require("sha256:"+sha256(read_text(root/"studies"/hex(reference.transaction)/"study.rfcxs"))==reference.digest,"reachable native snapshot bytes differ from retained publication");
    match(workspace, endpoint,reference.publication);
    return workspace;
}
Workspace endpoint_workspace(const fs::path &root,const native::TemporalEndpoint &endpoint) {
    return endpoint_workspace(root,endpoint,binding(endpoint));
}
// Canonical recovery validates native history; GADGET additionally validates every
// referenced biomedical snapshot, fixed geometry and retained tensor identities.
Workspace recover_history(
    const fs::path &root, const native::NativeTemporalSeries &series,
    const native::DurableNativeRoot &head,
    const std::function<std::optional<native::DurableNativeReceipt>(std::uint64_t)> &status) {
    require(series.source_namespace() == study_namespace && series.first_sequence() == 0 &&
                head.generation == series.last_sequence() &&
                series.segment_count() == series.last_sequence(),
            "unsupported native study history profile");
    std::optional<Admission> previous;
    Keys previous_keys;
    std::optional<Workspace> recovered;
    series.visit_endpoints([&](const auto &publication) {
        const auto endpoint = publication.pin();
        const auto sequence = endpoint->sequence;
        const auto reference = binding(*endpoint);
        if (sequence) {
            const auto receipt = status(reference.transaction);
            require(receipt && receipt->first_sequence + 1 == sequence &&
                        receipt->frame_count == 1 && receipt->root.sequence == sequence &&
                        receipt->root.generation == sequence,
                    "snapshot transaction disagrees with canonical index");
        }
        auto workspace = endpoint_workspace(root,*endpoint,reference);
        auto current = admission(workspace);
        auto keys = tensor_keys(*endpoint);
        if (previous) {
            admit(*previous, current);
            for (const auto &[name, key] : previous_keys)
                require(keys.contains(name) && keys.at(name) == key,
                        "retained tensor identity changed");
        }
        if (sequence == series.last_sequence()) {
            recovered = std::move(workspace);
            return;
        }
        previous = std::move(current);
        previous_keys = std::move(keys);
    });
    require(recovered.has_value(), "native history has no retained endpoint");
    return std::move(*recovered);
}
Workspace recover(const fs::path &root, native::DurableNativeStore &store) {
    return recover_history(root, store.snapshot(), store.root(),
                           [&](auto transaction) { return store.status(transaction); });
}
Json certificate(const native::DurableNativeRoot &root, const Workspace &workspace,
                 const std::string &profile = "linux_exclusive.v1") {
    auto result = Json::object();
    result.set("schema", Json::string("gadget.store.status.v1"));
    result.set("profile", Json::string(profile));
    result.set("root", root_json(root));
    result.set("source_digest", Json::string(workspace.source_digest));
    result.set("snapshot", workspace.snapshot_origin);
    result.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    result.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
    return result;
}
Json certificate(native::DurableNativeStore &store, const Workspace &workspace) {
    return certificate(store.root(), workspace);
}
void bind_workspace(Workspace &workspace, const native::TemporalEndpoint &endpoint) {
    auto keys = tensor_keys(endpoint);
    std::erase_if(keys, [&](const auto &item) { return !workspace.fields.contains(item.first); });
    workspace.native_snapshot =
        std::make_shared<NativeSnapshot>(endpoint, study_namespace, std::move(keys));
    const auto owner = workspace.native_snapshot->endpoint();
    workspace.complex = std::shared_ptr<const rfc::core::RelationalComplex>(owner, &owner->source);
}
Json receipt_json(const native::DurableNativeReceipt &receipt) {
    auto result = Json::object();
    result.set("root", root_json(receipt.root));
    result.set("transaction", Json::string(std::to_string(receipt.transaction)));
    result.set("first_sequence", Json::string(std::to_string(receipt.first_sequence)));
    result.set("frame_count", Json::string(std::to_string(receipt.frame_count)));
    return result;
}
} // namespace

std::uint64_t store_transaction(const std::string &value) {
    const auto result = decimal(value);
    require(result != 0, "transaction zero is reserved for genesis");
    return result;
}
Json initialize_store(const fs::path &snapshot, const fs::path &destination,
                      native::DurableFaultHook hook, StudyPreparationHook preparation_hook) {
    auto frozen = freeze(snapshot);
    auto endpoint = target(frozen.workspace, nullptr, 0);
    native::validate_file_path(destination);
    private_directory(destination);
    const auto root = fs::canonical(destination);
    private_directory(root / "studies");
    prepare(root, 0, frozen, preparation_hook);
    auto parent = durable_io::directory(root.parent_path());
    durable_io::synchronize(parent.get());
    auto store = native::DurableNativeStore::create(root / "history", endpoint, study_namespace,
                                                    limits(), hook);
    frozen.workspace.snapshot_origin.set("directory",
                                         Json::string((root / "studies" / hex(0)).string()));
    return certificate(store, frozen.workspace);
}
Json store_status(const fs::path &requested, std::optional<std::uint64_t> transaction) {
    if (transaction)
        require(*transaction != 0, "transaction zero is reserved for genesis");
    const auto root = store_path(requested);
    auto store = native::DurableNativeStore::open(root / "history", limits());
    auto workspace = recover(root, store);
    auto result = certificate(store, workspace);
    if (transaction) {
        const auto receipt = store.status(*transaction);
        result.set("transaction", Json::string(std::to_string(*transaction)));
        result.set("status", Json::string(receipt ? "committed" : "absent"));
        if (receipt)
            result.set("receipt", receipt_json(*receipt));
    }
    return result;
}
Json publish_store(const fs::path &requested, const fs::path &snapshot, const Json &request,
                   native::DurableFaultHook hook, StudyPreparationHook preparation_hook) {
    members(request, {"schema", "transaction", "expected"});
    require(request.at("schema").text() == "gadget.store.publish.v1", "unsupported request schema");
    const auto transaction = store_transaction(request.at("transaction").text());
    const auto expected = parse_root(request.at("expected"));
    const auto root = store_path(requested);
    auto store = native::DurableNativeStore::open(root / "history", limits());
    auto current = recover(root, store);
    const auto prior_receipt = store.status(transaction);
    require(prior_receipt || expected == store.root(), "stale expected predecessor root");
    const auto series = store.snapshot();
    require(expected.sequence < series.last_sequence() + (prior_receipt ? 0 : 1),
            "invalid original predecessor sequence");
    const auto previous = NativeStudy::from_history(series, {}, expected.sequence).endpoint();
    const auto predecessor = [&] {
        if (expected.sequence == store.root().sequence)
            return admission(current);
        return admission(endpoint_workspace(root, *previous));
    }();
    if (!prior_receipt)
        current = Workspace{};
    auto frozen = freeze(snapshot);
    admit(predecessor, admission(frozen.workspace));
    const auto endpoint = prior_receipt?*NativeStudy::from_history(series, {}, prior_receipt->root.sequence).endpoint():target(frozen.workspace, previous.get(), transaction);
    if(prior_receipt)match(frozen.workspace,endpoint);
    const std::array frames{
        native::encode_temporal_segment(native::make_temporal_segment(*previous, endpoint,{}, {},native::RelationalCodeEncoding::Automatic,{},native::make_temporal_declaration_update(previous->declarations,endpoint.declarations,rfc::composite::identity_declaration_maps(previous->declarations->registry))))};
    prepare(root, transaction, frozen, preparation_hook);
    const auto receipt = store.append(transaction, expected, frames, hook);
    frozen.workspace.snapshot_origin.set(
        "directory", Json::string((root / "studies" / hex(transaction)).string()));
    auto result = certificate(store, prior_receipt ? current : frozen.workspace);
    result.set("transaction", Json::string(std::to_string(transaction)));
    result.set("status", Json::string("committed"));
    result.set("receipt", receipt_json(receipt));
    return result;
}
Workspace load_store(const fs::path &requested) {
    const auto root = store_path(requested);
    auto store = native::DurableNativeStore::open(root / "history", limits());
    auto workspace = recover(root, store);
    workspace.store_origin = certificate(store, workspace);
    workspace.store_origin.set("directory", Json::string(root.string()));
    const auto endpoint = store.snapshot().tail().pin();
    bind_workspace(workspace, *endpoint);
    return workspace;
}

struct StudySnapshot::State {
    std::shared_ptr<const int> origin;
    native::NativeHistorySnapshot native;
    std::shared_ptr<const Workspace> workspace;
    State(std::shared_ptr<const int> service, native::NativeHistorySnapshot pin,
          std::shared_ptr<const Workspace> study)
        : origin(std::move(service)), native(std::move(pin)), workspace(std::move(study)) {}
};
StudySnapshot::StudySnapshot(std::shared_ptr<const State> state) : state_(std::move(state)) {}
native::DurableNativeRoot StudySnapshot::root() const {
    require(!!state_, "moved-from study snapshot");
    return state_->native.root();
}
std::shared_ptr<const Workspace> StudySnapshot::workspace() const {
    require(!!state_, "moved-from study snapshot");
    return state_->workspace;
}
struct ConcurrentStudyStore::State {
    fs::path root;
    std::unique_ptr<native::ConcurrentNativeHistory> history;
    std::shared_ptr<const int> origin = std::make_shared<const int>(0);
    std::atomic<std::shared_ptr<const StudySnapshot::State>> current;
    mutable std::recursive_mutex gate;
    mutable bool busy{};
    bool uncertain{};
    struct Call {
        State &state;
        std::lock_guard<std::recursive_mutex> lock;
        explicit Call(State &owner) : state(owner), lock(owner.gate) {
            if (state.busy)
                throw std::logic_error("study store writer/status/recovery reentry");
            state.busy = true;
        }
        ~Call() { state.busy = false; }
    };
    void process() const { (void)history->snapshot(); }
    void healthy() const {
        require(!uncertain,
                "publication outcome requires explicit study recovery before writes/status");
    }
    std::shared_ptr<const StudySnapshot::State> validate() {
        const auto pin = history->snapshot();
        auto workspace = std::make_shared<Workspace>(
            recover_history(root, pin.series(), pin.root(),
                            [&](auto transaction) { return history->status(transaction); }));
        bind_workspace(*workspace, *pin.tail().pin());
        workspace->store_origin = certificate(pin.root(), *workspace, "linux_concurrent.v1");
        workspace->store_origin.set("directory", Json::string(root.string()));
        return std::make_shared<const StudySnapshot::State>(origin, pin, std::move(workspace));
    }
};
ConcurrentStudyStore::ConcurrentStudyStore(const fs::path &directory)
    : state_(std::make_unique<State>()) {
    state_->root = store_path(directory);
    state_->history.reset(new native::ConcurrentNativeHistory(
        native::ConcurrentNativeHistory::open(state_->root / "history", limits())));
    state_->current.store(state_->validate(), std::memory_order_release);
}
ConcurrentStudyStore::~ConcurrentStudyStore() = default;
StudySnapshot ConcurrentStudyStore::snapshot() const {
    state_->process();
    return StudySnapshot(state_->current.load(std::memory_order_acquire));
}
Json ConcurrentStudyStore::status(std::optional<std::uint64_t> transaction) const {
    state_->process();
    State::Call call(*state_);
    state_->healthy();
    const auto pin = state_->current.load(std::memory_order_acquire);
    auto result = certificate(pin->native.root(), *pin->workspace, "linux_concurrent.v1");
    if (transaction) {
        const auto receipt = state_->history->status(*transaction);
        result.set("transaction", Json::string(std::to_string(*transaction)));
        result.set("status", Json::string(receipt ? "committed" : "absent"));
        if (receipt)
            result.set("receipt", receipt_json(*receipt));
    }
    return result;
}
Json ConcurrentStudyStore::publish(const StudySnapshot &expected, const fs::path &snapshot_path,
                                   const Json &request, std::stop_token stop,
                                   native::DurableFaultHook hook,
                                   StudyPreparationHook preparation_hook) {
    state_->process();
    State::Call call(*state_);
    state_->healthy();
    require(expected.state_ && expected.state_->origin == state_->origin,
            "foreign study snapshot token");
    members(request, {"schema", "transaction", "expected"});
    require(request.at("schema").text() == "gadget.store.publish.v1", "unsupported request schema");
    const auto transaction = store_transaction(request.at("transaction").text());
    require(parse_root(request.at("expected")) == expected.root(),
            "request disagrees with authentic expected snapshot");
    const auto prior_receipt = state_->history->status(transaction);
    if (!prior_receipt && expected.root() != state_->current.load()->native.root())
        throw native::NativeHistoryConflict();
    if (!prior_receipt && stop.stop_requested())
        throw native::NativeHistoryCancelled();
    auto frozen = freeze(snapshot_path);
    admit(admission(*expected.workspace()), admission(frozen.workspace));
    const auto previous = expected.state_->native.tail().pin();
    const auto endpoint = prior_receipt?*NativeStudy::from_history(state_->history->snapshot().series(), {}, prior_receipt->root.sequence).endpoint():target(frozen.workspace, previous.get(), transaction);
    if(prior_receipt)match(frozen.workspace,endpoint);
    const std::array frames{
        native::encode_temporal_segment(native::make_temporal_segment(*previous, endpoint,{}, {},native::RelationalCodeEncoding::Automatic,{},native::make_temporal_declaration_update(previous->declarations,endpoint.declarations,rfc::composite::identity_declaration_maps(previous->declarations->registry))))};
    prepare(state_->root, transaction, frozen, [&](auto event) {
        if (preparation_hook)
            preparation_hook(event);
        if (!prior_receipt && stop.stop_requested())
            throw native::NativeHistoryCancelled();
    });
    bind_workspace(frozen.workspace, endpoint);
    frozen.workspace.snapshot_origin.set(
        "directory", Json::string((state_->root / "studies" / hex(transaction)).string()));
    auto workspace = std::make_shared<Workspace>(std::move(frozen.workspace));
    auto publication =
        std::make_shared<StudySnapshot::State>(state_->origin, expected.state_->native, workspace);
    native::DurableNativeReceipt receipt;
    try {
        receipt = state_->history->append(
            transaction, expected.state_->native, frames, stop, [&](auto event) {
                if (event.stage == native::DurableStage::LivePublished) {
                    publication->native = state_->history->snapshot();
                    workspace->store_origin =
                        certificate(publication->native.root(), *workspace, "linux_concurrent.v1");
                    workspace->store_origin.set("directory", Json::string(state_->root.string()));
                    state_->current.store(publication, std::memory_order_release);
                }
                if (hook)
                    hook(event);
            });
    } catch (const native::DurableUnknownOutcome &) {
        state_->uncertain = true;
        throw;
    }
    const auto current = state_->current.load(std::memory_order_acquire);
    auto result = certificate(current->native.root(), *current->workspace, "linux_concurrent.v1");
    result.set("transaction", Json::string(std::to_string(transaction)));
    result.set("status", Json::string("committed"));
    result.set("receipt", receipt_json(receipt));
    return result;
}
void ConcurrentStudyStore::recover() {
    state_->process();
    State::Call call(*state_);
    state_->uncertain = true;
    state_->history->recover();
    auto publication = state_->validate();
    state_->current.store(std::move(publication), std::memory_order_release);
    state_->uncertain = false;
}
} // namespace gadget
