#include <gadget/native_pattern.hpp>
#include <gadget/record_view.hpp>
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <rfc/learning/context_transport_pattern_bank_codec.hpp>
#include <limits>
#include <set>

namespace gadget {
namespace {
using namespace rfc;
using N = NativeRecordValue;
using V = numeric::RationalVector;
using M = numeric::ExactMatrix;
using R = scientific_detail::ResultValue;
using Fields = scientific_detail::ResultFields;
using Bank = learning::PreparedTransportPatternBank;
using Guard = numeric::ExactActionGuard;
using scientific_detail::require;
std::size_t number(const NativeRecord &value) {
    const auto n = value.integer();
    require(n >= 0 && n <= numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())),
            "pattern index exceeds machine addressability");
    return static_cast<std::size_t>(std::stoull(n.get_str()));
}
void keys(const NativeRecord &value, std::set<std::string> required,
          const std::set<std::string> &optional = {}) {
    require(value.kind() == NativeRecordKind::Object, "pattern declaration requires an object");
    for (const auto &name : value.names()) {
        if (required.erase(name)) continue;
        require(optional.contains(name), "unknown pattern declaration member: " + name);
    }
    require(required.empty(), "pattern declaration member missing");
}
std::string hash(std::span<const std::byte> bytes) {
    return "sha256:" + sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
}
std::string hash(const NativeRecord &record) { return hash(record.archive_source(1).bytes()); }
std::string hash(const std::array<std::byte,32> &bytes) {
    constexpr char hex[] = "0123456789abcdef";
    std::string out = "sha256:";
    for (auto byte : bytes) {
        const auto value = std::to_integer<unsigned>(byte);
        out.push_back(hex[value >> 4]); out.push_back(hex[value & 15]);
    }
    return out;
}
Guard action_guard(const runtime::ExecutionLease &lease) {
    const auto &limits = lease.policy().limits;
    return Guard({std::numeric_limits<std::size_t>::max(), limits.sparse_entries,
                  limits.scalar_work, limits.integer_bits}, [lease] { lease.check_cancelled(); });
}
M matrix(const NativeRecord &record, Guard &guard) {
    keys(record, {"rows", "columns", "entries"});
    M out(number(record.at("rows")), number(record.at("columns")));
    const auto entries = record.at("entries");
    require(entries.kind() == NativeRecordKind::Array, "pattern map entries require an array");
    guard.entries(entries.size());
    std::set<std::pair<std::size_t,std::size_t>> seen;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        guard.step(); const auto entry = entries.at(i);
        require(entry.kind() == NativeRecordKind::Array && entry.size() == 3,
                "pattern map entry requires row, column and rational");
        const auto row = number(entry.at(0)), column = number(entry.at(1));
        require(row < out.rows() && column < out.cols() && seen.emplace(row,column).second,
                "pattern map coordinate is repeated or outside its declared domain");
        auto value = exact_rational(entry.at(2).text()); guard.value(value);
        out.set(row,column,std::move(value));
    }
    guard.matrix(out); return out;
}
const storage::TemporalCarrier &carrier(const Bank &bank, std::uint64_t key) {
    for (const auto &value : bank.endpoint().carriers) if (value.key == key) return value;
    throw std::invalid_argument("original bank carrier absent");
}
std::uint64_t input_carrier(const Bank &bank, std::size_t pattern) {
    require(pattern < bank.input_images().size(), "unknown original bank pattern");
    return machine::checked_size_add(machine::checked_size_add(3,bank.input_images().size()),pattern);
}
NativeRecord meaning(const NativeRecord &declaration, std::uint64_t frame) {
    const auto meanings = declaration.at("frames");
    for (std::size_t i = 0; i < meanings.size(); ++i)
        if (number(meanings.at(i).at("frame")) == frame) return meanings.at(i);
    throw std::invalid_argument("original pattern frame has no declared quantity or units");
}
void declaration(const Bank &bank, const NativeRecord &value) {
    keys(value, {"schema", "frames", "scope"});
    require(value.at("schema").text() == "gadget.pattern.model.declaration.v1" &&
                !value.at("scope").text().empty(), "pattern model declaration scope");
    const auto frames = value.at("frames");
    require(frames.kind() == NativeRecordKind::Array, "pattern frame meanings require an array");
    std::set<std::uint64_t> seen;
    for (std::size_t i = 0; i < frames.size(); ++i) {
        const auto frame = frames.at(i);
        keys(frame, {"frame", "quantity", "units"});
        require(number(frame.at("frame")) && seen.insert(number(frame.at("frame"))).second &&
                    !frame.at("quantity").text().empty() && !frame.at("units").text().empty(),
                "pattern frame meanings require distinct nonzero frames and physical meanings");
    }
    for (const auto &c : bank.endpoint().carriers)
        if (c.key >= 3 + bank.input_images().size()) (void)meaning(value,c.coordinate_system);
}
struct BoundInput { learning::TransportPatternSparseInput input; M map; NativeRecord request; };
BoundInput bind(const Bank &bank, const NativeRecord &model_declaration,
                const ScientificStudy &study, const NativeRecord &request, Guard &guard) {
    keys(request, {"pattern", "frame", "field", "domain_quantity", "domain_units",
                   "codomain_quantity", "codomain_units", "method", "evidence", "map"});
    const auto pattern = number(request.at("pattern")), frame = number(request.at("frame"));
    const auto target = input_carrier(bank,pattern);
    require(frame == carrier(bank,target).coordinate_system, "pattern input belongs to another bank frame");
    const auto expected = meaning(model_declaration,frame);
    const auto field = request.at("field").text();
    const auto description = study.field_declaration(field);
    require(request.at("domain_quantity").text() == description.at("quantity").text() &&
                request.at("domain_units").text() == description.at("units").text() &&
                request.at("codomain_quantity").text() == expected.at("quantity").text() &&
                request.at("codomain_units").text() == expected.at("units").text() &&
                !request.at("method").text().empty() && !request.at("evidence").text().empty(),
            "pattern input requires an explicit source and model quantity map");
    const auto pin = study.field_accession(field); const auto &tensor = *pin.tensor();
    auto map = matrix(request.at("map"),guard);
    require(map.cols() == tensor.domain().logical_size() && map.rows() == bank.input_images()[pattern].rows(),
            "pattern map does not join the actual source tensor and original bank field");
    for (const auto &[coordinate,value] : map.entries()) {
        (void)value; guard.step();
        require(tensor.support().contains(coordinate.second), "pattern input depends on missing source support");
    }
    guard.entries(machine::checked_size_add(map.nonzeros(),tensor.stored_nonzeros()));
    auto values = map.apply_sparse(tensor.values().view()); guard.sparse_vector(values);
    learning::TransportPatternSparseInput input{pattern,frame,numeric::SparseRationalTensor(
        numeric::IndexDomain::linear(map.rows()), numeric::StructuralSupport::full(map.rows()),std::move(values))};
    return {std::move(input),std::move(map),request};
}
R response(Fields &fields, const Bank &bank, const learning::TransportPatternBankResponse &value) {
    auto out = R::object(); const auto axis = fields.axis(value.origin.size());
    out.set("origin",fields.vector(value.origin,axis));
    out.set("directions",fields.matrix(value.directions,fields.axis(value.directions.cols()),axis));
    out.set("offsets",fields.indices(value.offsets,fields.axis(value.offsets.size())));
    auto coordinates = R::array();
    for (std::size_t i = 0; i < value.coordinates.size(); ++i) {
        const auto &v = value.coordinates[i]; auto coordinate = R::object();
        coordinate.set("pattern",R::count(i)); coordinate.set("carrier",R::count(3+i));
        coordinate.set("value",v.empty() && bank.input_images()[i].cols() ? R{} : fields.vector(v,3+i));
        coordinates.push(std::move(coordinate));
    }
    out.set("coordinates",std::move(coordinates)); return out;
}
}

NativeRecord publish_pattern_model(const Bank &bank, std::uint64_t bank_namespace,
    const NativeRecord &source, const NativeRecord &evidence, const NativeRecord &model_declaration) {
    require(bank_namespace, "pattern bank namespace must be nonzero");
    require(source.is_root() && evidence.is_root() && model_declaration.is_root(),
            "pattern training inputs and declaration require their actual root owners");
    declaration(bank,model_declaration);
    auto source_bytes = source.archive_source(1), evidence_bytes = evidence.archive_source(2);
    require(hash(source_bytes.bytes()) == hash(bank.provenance().source_sha256) &&
                hash(evidence_bytes.bytes()) == hash(bank.provenance().evidence_sha256),
            "bank training provenance differs from actual original source or evidence bytes");
    auto encoded = learning::encode_transport_pattern_bank(bank,bank_namespace,native_data_limits());
    auto bytes = RawByteSource::own(3,std::move(encoded)); auto controls = model_declaration.archive_source(4);
    N value(N::Object{{{"schema",N(std::string("gadget.pattern.model.native.v1"))},
        {"bank_namespace",N(numeric::Integer(std::to_string(bank_namespace)))},
        {"source_publication",N(N::ByteSpan{1,0,source_bytes.bytes().size(),false})},
        {"evidence_publication",N(N::ByteSpan{2,0,evidence_bytes.bytes().size(),false})},
        {"bank_publication",N(N::ByteSpan{3,0,bytes.bytes().size(),false})},
        {"declaration_publication",N(N::ByteSpan{4,0,controls.bytes().size(),false})}}});
    storage::TemporalEndpoint endpoint(bank.endpoint().source,bank.endpoint().sequence);
    endpoint.carriers = bank.endpoint().carriers;
    return NativeRecord::create_from_endpoint(std::move(value),new_native_record_namespace(),
        {std::move(source_bytes),std::move(evidence_bytes),std::move(bytes),std::move(controls)},std::move(endpoint));
}
struct NativePatternModel::State {
    NativeRecord publication, declaration;
    Bank bank;
    std::string digest;
    explicit State(NativeRecord owner, Guard &guard)
        :publication(std::move(owner)),
         declaration(NativeRecord::decode(publication.at("declaration_publication").bytes())),
         bank(learning::decode_transport_pattern_bank(publication.at("bank_publication").bytes(),
            number(publication.at("bank_namespace")),guard,native_data_limits())),digest(hash(publication)) {
        require(publication.at("schema").text() == "gadget.pattern.model.native.v1", "native pattern model schema");
        require(hash(publication.at("source_publication").bytes()) == hash(bank.provenance().source_sha256) &&
                    hash(publication.at("evidence_publication").bytes()) == hash(bank.provenance().evidence_sha256),
                "pattern bank source/evidence provenance mismatch");
        const auto geometry = [](const core::RelationalComplex &source) {
            return storage::encode_native_series(storage::NativeTemporalSeries(
                storage::TemporalEndpoint(source),1,native_data_limits()));
        };
        require(geometry(publication.endpoint()->source) == geometry(bank.endpoint().source),
                "pattern wrapper changed its original bank complex");
        for (const auto &expected : bank.endpoint().carriers) {
            const auto actual = std::find_if(publication.endpoint()->carriers.begin(),
                publication.endpoint()->carriers.end(),[&](const auto &c) { return c.key == expected.key; });
            require(actual != publication.endpoint()->carriers.end() && actual->kind == expected.kind &&
                        actual->grade == expected.grade && actual->cell == expected.cell &&
                        actual->coordinate_system == expected.coordinate_system && actual->domain == expected.domain,
                    "pattern wrapper changed an original bank carrier");
        }
        gadget::declaration(bank,declaration);
    }
};
NativePatternModel::NativePatternModel(NativeRecord publication) {
    require(publication.is_root() && publication.at("schema").text() == "gadget.pattern.model.native.v1",
            "pattern admission requires an actual model root publication");
    const auto lease = scientific_detail::source_execution_lease("gadget-native-pattern-admission");
    auto guard = action_guard(lease); state_ = std::make_shared<const State>(std::move(publication),guard);
}
const NativeRecord &NativePatternModel::native_publication() const { return state_->publication; }
std::string NativePatternModel::digest() const { return state_->digest; }

NativeRecord NativePatternModel::execute_record(std::shared_ptr<const ScientificStudy> source,
    const NativeRecord &request, std::shared_ptr<CpuResources> resources) const {
    require(bool(source),"native pattern request requires an actual scientific source owner");
    keys(request,{"schema","operation","source_digest","model_sha256"},{"inputs","nodes","comparisons"});
    require(request.at("schema").text() == "gadget.pattern.query.v1" &&
                request.at("source_digest").text() == source->digest() &&
                request.at("model_sha256").text() == digest(), "native pattern source/model binding mismatch");
    const auto operation = request.at("operation").text();
    require(operation == "sparse_candidates" || operation == "sparse_program" || operation == "sparse_gluing",
            "unsupported native pattern operation");
    if (!resources) resources = CpuResources::local();
    auto reservation = resources->reserve(1);
    const auto lease = scientific_detail::source_execution_lease("gadget-native-pattern-response");
    auto guard = action_guard(lease); Fields fields(state_->publication); auto view = R::object();
    view.set("operation",R::string(operation)); view.set("model_sha256",R::string(digest()));
    view.set("source_digest",R::string(source->digest()));
    const auto &bank = state_->bank; auto bound = R::array();
    std::map<std::shared_ptr<const void>,std::map<std::uint64_t,std::uint64_t>,
        std::owner_less<std::shared_ptr<const void>>> control_sources;
    const auto controls = [&](const NativeRecord &record) {
        const auto previous = control_sources.find(record.publication_owner());
        if (previous != control_sources.end()) return R(record.unbound_value(previous->second));
        std::map<std::uint64_t,std::uint64_t> keys; std::vector<RawByteSource> sources;
        auto next = fields.next_source_key();
        for (const auto &bytes : record.raw_sources()) {
            keys.emplace(bytes.key(),next); sources.push_back(bytes.rekey(next));
            next = machine::checked_size_add(next,1);
        }
        auto copied = R(record.unbound_value(keys)); fields.retain_sources(std::move(sources));
        control_sources.emplace(record.publication_owner(),std::move(keys)); return copied;
    };

    const auto input = [&](const NativeRecord &record) {
        auto result = bind(bank,state_->declaration,*source,record,guard);
        auto mapping = R::object(); mapping.set("field",R::string(record.at("field").text()));
        const auto pin = source->field_accession(record.at("field").text());
        mapping.set("source_namespace",R::string(std::to_string(pin.source_namespace())));
        mapping.set("cutoff",R::count(pin.sequence()));
        mapping.set("tensor",R::count(pin.tensor_declaration()->key));
        mapping.set("carrier",R::count(pin.carrier()->key));
        fields.bind(mapping);
        mapping.set("declaration",controls(record));
        mapping.set("map",fields.matrix(result.map,fields.axis(result.map.cols()),
            input_carrier(bank,result.input.pattern)));
        bound.push(std::move(mapping)); return result.input;
    };
    if (operation == "sparse_program") {
        require(request.has("nodes") && !request.has("inputs") && !request.has("comparisons"),
                "pattern programme requires only its ordered nodes");
        const auto nodes = request.at("nodes"); require(nodes.kind() == NativeRecordKind::Array,"pattern nodes require an array");
        std::vector<learning::TransportPatternProgramNode> programme; std::map<std::string,std::size_t> indices;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const auto node = nodes.at(i); keys(node,{"id","pattern"},{"input","edges"});
            learning::TransportPatternProgramNode step; step.name = node.at("id").text(); step.pattern = number(node.at("pattern"));
            require(!step.name.empty() && indices.emplace(step.name,i).second,"pattern node identities must be distinct");
            if (node.has("input")) { require(!node.has("edges"),"pattern node cannot mix input and incoming maps"); step.input = input(node.at("input")); }
            else {
                const auto edges = node.at("edges"); require(edges.kind() == NativeRecordKind::Array,"pattern edges require an array");
                for (std::size_t j = 0; j < edges.size(); ++j) {
                    const auto edge = edges.at(j); keys(edge,{"from","occurrence_key","source_frame","target_frame","method","evidence","map"});
                    require(!edge.at("method").text().empty() && !edge.at("evidence").text().empty(),
                            "pattern programme map requires its declared method and evidence");
                    (void)meaning(state_->declaration,number(edge.at("source_frame")));
                    (void)meaning(state_->declaration,number(edge.at("target_frame")));
                    const auto from = edge.at("from").text(); require(indices.contains(from) && indices.at(from) < i,"pattern edge requires an earlier named node");
                    const auto occurrence = storage::resolve_stable_occurrence(bank.endpoint().source,
                        storage::StableOccurrenceKey{core::Grade{1},core::OccurrenceKey{number(edge.at("occurrence_key"))}});
                    require(occurrence.has_value(),"pattern edge original occurrence key absent");
                    step.edges.push_back({indices.at(from),*occurrence,number(edge.at("source_frame")),
                        number(edge.at("target_frame")),matrix(edge.at("map"),guard)});
                }
            }
            programme.push_back(std::move(step));
        }
        const auto values = bank.response_sparse_program(programme,guard); auto stages = R::array();
        for (std::size_t i = 0; i < values.size(); ++i) {
            const auto &value = values[i]; auto stage = R::object(); stage.set("id",R::string(programme[i].name));
            const auto status = [&] {
                switch (value.status) {
                case learning::TransportPatternProgramStatus::Supported: return "supported";
                case learning::TransportPatternProgramStatus::OutsideDomain: return "outside_domain";
                case learning::TransportPatternProgramStatus::OpenInput: return "open_input";
                case learning::TransportPatternProgramStatus::Blocked: return "blocked";
                }
                throw std::logic_error("unknown canonical pattern programme status");
            }();
            stage.set("status",R::string(status));
            stage.set("input",value.input ? Fields::staged(*value.input,{},input_carrier(bank,programme[i].pattern),"vector") : R{});
            stage.set("input_directions",value.status == learning::TransportPatternProgramStatus::Blocked ? R{} :
                fields.matrix(value.input_directions,2,input_carrier(bank,programme[i].pattern)));
            stage.set("parameter_directions",value.response ?
                fields.matrix(value.parameter_directions,2,fields.axis(value.parameter_directions.rows())) : R{});
            if (value.response) stage.set("response",response(fields,bank,*value.response)); else stage.set("response",R{});
            stages.push(std::move(stage));
        }
        view.set("stages",std::move(stages));
    } else {
        require(request.has("inputs") && !request.has("nodes"),"native pattern inputs missing");
        const auto declared = request.at("inputs"); require(declared.kind() == NativeRecordKind::Array,"pattern inputs require an array");
        std::vector<learning::TransportPatternSparseInput> inputs;
        for (std::size_t i = 0; i < declared.size(); ++i) inputs.push_back(input(declared.at(i)));
        if (operation == "sparse_candidates") {
            require(!request.has("comparisons"),"candidate queries do not declare gluing comparisons");
            auto candidates = R::array();
            for (const auto &value : bank.response_sparse_candidates(inputs,guard)) {
                auto candidate = R::object(); candidate.set("pattern",R::count(value.pattern));
                candidate.set("status",R::string(value.response ? "supported" : "outside_domain"));
                candidate.set("response",value.response ? response(fields,bank,*value.response) : R{});
                candidates.push(std::move(candidate));
            }
            view.set("candidates",std::move(candidates));
        } else {
            const auto comparisons = request.at("comparisons"); require(comparisons.kind() == NativeRecordKind::Array,"pattern comparisons require an array");
            std::vector<learning::TransportPatternComparison> joins;
            for (std::size_t i = 0; i < comparisons.size(); ++i) {
                const auto c = comparisons.at(i); keys(c,{"name","left_occurrence_key","right_occurrence_key","left_frame","right_frame","interface_frame","method","evidence","left_map","right_map"});
                require(!c.at("method").text().empty() && !c.at("evidence").text().empty(),
                        "pattern comparison requires its declared method and evidence");
                for (const auto *frame : {"left_frame","right_frame","interface_frame"})
                    (void)meaning(state_->declaration,number(c.at(frame)));
                const auto original = [&](const char *key) {
                    const auto slot = storage::resolve_stable_occurrence(bank.endpoint().source,
                        storage::StableOccurrenceKey{core::Grade{1},core::OccurrenceKey{number(c.at(key))}});
                    require(slot.has_value(),"pattern comparison original occurrence absent"); return *slot;
                };
                joins.push_back({c.at("name").text(),original("left_occurrence_key"),original("right_occurrence_key"),
                    number(c.at("left_frame")),number(c.at("right_frame")),number(c.at("interface_frame")),
                    matrix(c.at("left_map"),guard),matrix(c.at("right_map"),guard)});
            }
            const auto value = bank.glue_sparse(inputs,joins,guard);
            view.set("consistent",R::boolean(value.consistent)); view.set("forced",R::boolean(value.forced));
            const auto defect_axis = fields.axis(value.defect_origin.size());
            view.set("defect_origin",fields.vector(value.defect_origin,defect_axis));
            view.set("defect_directions",fields.matrix(value.defect_directions,2,defect_axis));
            view.set("contradiction",value.consistent ? R{} : fields.vector(value.contradiction,defect_axis));
            view.set("parameter_origin",value.consistent ? fields.vector(value.parameter_origin,2) : R{});
            const auto conditioned = fields.axis(value.parameter_directions.cols());
            view.set("parameter_directions",value.consistent ? fields.matrix(value.parameter_directions,conditioned,2) : R{});
            view.set("offsets",fields.indices(value.offsets,fields.axis(value.offsets.size())));
            view.set("section_cells",fields.indices(value.section_cells,fields.axis(value.section_cells.size())));
            auto names = R::array(); for (const auto &name : value.comparison_names) names.push(R::string(name));
            view.set("comparison_names",std::move(names));
            view.set("boundary_response",value.boundary_response ? response(fields,bank,*value.boundary_response) : R{});
            view.set("section",value.section ? response(fields,bank,*value.section) : R{});
        }
    }
    view.set("input_bindings",std::move(bound));
    view.set("model_declaration",controls(state_->declaration));
    const auto result = fields.finish(std::move(view));
    auto source_bytes = source->native_publication().archive_source(1), model_bytes = state_->publication.archive_source(2);
    std::uint64_t next = 1; for (const auto &bytes : result.raw_sources()) next = std::max(next,machine::checked_size_add(bytes.key(),1));
    source_bytes = source_bytes.rekey(next); const auto source_key = next; next = machine::checked_size_add(next,1);
    model_bytes = model_bytes.rekey(next);
    N body(N::Object{{{"schema",N(std::string("gadget.pattern.result.fields.v1"))},
        {"source_digest",N(source->digest())},{"model_sha256",N(digest())},{"fields",result.unbound_value()},
        {"source_study_publication",N(N::ByteSpan{source_key,0,source_bytes.bytes().size(),false})},
        {"model_publication",N(N::ByteSpan{next,0,model_bytes.bytes().size(),false})}}});
    const auto retained = result.derive(std::move(body),{std::move(source_bytes),std::move(model_bytes)});
    return retain_operation(retained,request,retained.unbound_value());
}
Json pattern_record_view(const NativeRecord &publication) {
    const auto result = publication.at("result");
    require(result.at("schema").text() == "gadget.pattern.result.fields.v1","native pattern view schema");
    auto view = Json::object();
    view.set("schema",Json::string("gadget.pattern.certificate.v1"));
    view.set("source_digest",render_record_view(result.at("source_digest")));
    view.set("model_sha256",render_record_view(result.at("model_sha256")));
    view.set("query",render_record_view(publication.at("request")));
    view.set("result",scientific_detail::result_fields_view(result.at("fields")));
    return scientific_detail::seal(view,"certificate_sha256");
}
NativeRecord replay_pattern_record(const NativeRecord &publication, std::shared_ptr<CpuResources> resources) {
    const auto result = publication.at("result");
    require(result.at("schema").text() == "gadget.pattern.result.fields.v1","native pattern result schema");
    auto source = std::make_shared<const ScientificStudy>(NativeRecord::decode(result.at("source_study_publication").bytes()));
    const NativePatternModel model(NativeRecord::decode(result.at("model_publication").bytes()));
    const auto request = NativeRecord::decode(publication.at("request_publication").bytes());
    const auto replay = model.execute_record(std::move(source),request,std::move(resources));
    require(pattern_record_view(publication).dump() == pattern_record_view(replay).dump(),
            "pattern replay changed a complete result, status or declared binding");
    const auto &ca = publication.endpoint()->carriers, &cb = replay.endpoint()->carriers;
    require(ca.size() == cb.size(),"pattern replay carrier family changed");
    for (std::size_t i = 0; i < ca.size(); ++i)
        require(ca[i].key == cb[i].key && ca[i].kind == cb[i].kind && ca[i].grade == cb[i].grade &&
                    ca[i].cell == cb[i].cell && ca[i].coordinate_system == cb[i].coordinate_system && ca[i].domain == cb[i].domain,
                "pattern replay changed an actual carrier declaration");
    const auto &a = publication.endpoint()->tensors, &b = replay.endpoint()->tensors;
    require(a.size() == b.size(),"pattern replay tensor family changed");
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].role == storage::TemporalTensorRole::OpaqueExact) continue;
        const auto same_tensor = [&] {
            if (a[i].tensor.domain() != b[i].tensor.domain() ||
                a[i].tensor.support().count() != b[i].tensor.support().count() ||
                !std::ranges::equal(a[i].tensor.values().entries(),b[i].tensor.values().entries())) return false;
            bool supported = true;
            if (a[i].tensor.support().count() != a[i].tensor.domain().logical_size())
                a[i].tensor.support().for_each_active([&](auto index) {
                    supported = supported && b[i].tensor.support().contains(index);
                });
            return supported;
        }();
        require(a[i].role == b[i].role && a[i].input == b[i].input && a[i].output == b[i].output &&
                    a[i].interpretation == b[i].interpretation && same_tensor,
                "pattern replay changed an actual result tensor or support");
    }
    return replay;
}
} // namespace gadget
