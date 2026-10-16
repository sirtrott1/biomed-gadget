#include "../examples/scientific/fixture.hpp"
#include "../src/scientific_result_fields.hpp"
#include <gadget/native_pattern.hpp>
#include <gadget/workflow.hpp>
#include <rfc/learning/context_transport_pattern_bank_codec.hpp>
#include <algorithm>
#include <span>
#include <iostream>

namespace {
using namespace gadget;
using namespace rfc;
using N = NativeRecordValue;
using Q = numeric::Rational;
using V = numeric::RationalVector;
using M = numeric::ExactMatrix;
using G = numeric::ExactActionGuard;
using F = scientific_detail::ResultFields;
using R = scientific_detail::ResultValue;
std::size_t checks{};
void check(bool value,const char *message) { ++checks; if (!value) throw std::runtime_error(message); }
template<class Fn> void rejects(Fn fn,const char *message) {
    bool refused = false; try { fn(); } catch (const std::exception &) { refused = true; }
    check(refused,message);
}
std::array<std::byte,32> digest(const NativeRecord &publication) {
    const auto owner = publication.archive_source(1); const auto bytes = owner.bytes();
    const auto hex = sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()),bytes.size()));
    std::array<std::byte,32> out{};
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<std::byte>(std::stoul(hex.substr(2*i,2),nullptr,16));
    return out;
}
struct Training {
    learning::PreparedTransportPatternBank bank;
    NativeRecord source,evidence;
};
Training train() {
    G guard; core::RelationSignature signature;
    const auto context = signature.add_type("context",core::BoundaryPolicy::DeclaredExact,0,1);
    const auto anchor = signature.add_type("anchor",core::BoundaryPolicy::DeclaredExact,0,1);
    core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < 6; ++i)
        builder.add_grade_one(i < 3 ? context : anchor,{{0,true,false,Q{1},{}}});
    const auto source = builder.freeze(); local::LinearStalkBundle stalks(source);
    for (std::size_t i = 0; i < 6; ++i)
        stalks.set_dimension({core::Grade{1},static_cast<std::uint32_t>(i)},2);
    local::LinearRelationalFiberSystem fibers(stalks);
    std::vector<learning::TransportCoordinate> parameters;
    for (std::size_t i = 0; i < 6; ++i) {
        fibers.set_transport(core::Grade{1},i,M(1,2));
        for (std::size_t j = 0; j < 2; ++j) parameters.push_back({i,0,j});
    }
    learning::PreparedTransportObservations observations(local::OwnedLinearFiberSystem(fibers),{},parameters,guard);
    std::vector<learning::TransportObservation> samples;
    for (auto i : {1u,2u}) {
        learning::TransportObservation sample; sample.input = V(12);
        sample.input[0] = sample.input[1] = -1;
        sample.input[2*i] = sample.input[2*i+1] = 1;
        sample.target = {0}; sample.observed = {1}; samples.push_back(std::move(sample));
    }
    for (std::size_t i = 3; i < 6; ++i) {
        learning::TransportObservation sample; sample.input = V(12);
        sample.input[2*i] = sample.input[2*i+1] = Q{1,2};
        sample.target = {1}; sample.observed = {1}; samples.push_back(std::move(sample));
    }
    storage::TemporalEndpoint endpoint(source);
    endpoint.carriers.push_back({99,storage::TemporalCarrierKind::CellStalk,core::Grade{0},
        storage::stable_cell_key(source,{core::Grade{0},0}),1,numeric::IndexDomain::linear(1)});
    for (std::size_t i = 0; i < 6; ++i) {
        endpoint.carriers.push_back({100+i,storage::TemporalCarrierKind::CellStalk,core::Grade{1},
            storage::stable_cell_key(source,{core::Grade{1},static_cast<std::uint32_t>(i)}),2,numeric::IndexDomain::linear(2)});
        storage::TemporalTensor transport; transport.key = i+1;
        transport.role = storage::TemporalTensorRole::PhysicalRealization; transport.input = 100+i; transport.output = 99;
        transport.interpretation = "original declared zero transport; both coordinates variable";
        transport.tensor = storage::temporal_matrix_tensor(M(1,2)); endpoint.tensors.push_back(std::move(transport));
    }
    const std::string note = "Declared synthetic training: six original context stars, two observed input channels; no clinical inference.\n";
    auto bytes = RawByteSource(1,std::as_bytes(std::span(note.data(),note.size())));
    const auto owner = NativeRecord::create_from_endpoint(N(N::Object{{
        {"schema",N(std::string("gadget.test.pattern.source.v1"))},
        {"original",N(N::ByteSpan{1,0,note.size(),false})}}}),new_native_record_namespace(),{bytes},std::move(endpoint));
    learning::PreparedTransportEvidence evidence(observations,guard);
    evidence = evidence.append(samples,guard);
    F fields(owner); auto view = R::object(); auto history = R::array();
    const auto input_axis = fields.axis(12), output_axis = fields.axis(1), parameter_axis = fields.axis(12);
    for (const auto &sample : samples) {
        auto row = R::object(); row.set("input",fields.vector(sample.input,input_axis));
        row.set("target",fields.vector(sample.target,output_axis));
        row.set("observed",fields.vector(V{Q{1}},output_axis)); history.push(std::move(row));
    }
    view.set("observations",std::move(history));
    view.set("design",fields.matrix(evidence.system().design,parameter_axis,fields.axis(evidence.system().design.rows())));
    view.set("target",fields.vector(evidence.system().target,fields.axis(evidence.system().target.size())));
    view.set("augmented",fields.matrix(evidence.system().augmented,fields.axis(evidence.system().augmented.cols()),fields.axis(evidence.system().augmented.rows())));
    view.set("input_scope",R::string("observed_fields"));
    view.set("minimum_peers",R::count(2));
    view.set("training_slots",R(Json::parse("[0,1,2,3,4,5]")));
    view.set("origin_slots",R(Json::parse("[0,1,3]")));
    view.set("lower_atlas_accession",fields.matrix(M::identity(1),output_axis,output_axis));
    view.set("upper_atlas_accession",fields.matrix(M::identity(2),fields.axis(2),fields.axis(2)));
    view.set("particular",fields.vector(evidence.particular(),parameter_axis));
    view.set("full_conditioning_kernel",fields.matrix(evidence.kernel(),fields.axis(evidence.kernel().cols()),parameter_axis));
    auto coordinates = R::array();
    for (const auto &p : parameters) {
        auto entry = R::object(); entry.set("occurrence_slot",R::count(p.occurrence));
        entry.set("row",R::count(p.row)); entry.set("column",R::count(p.column)); coordinates.push(std::move(entry));
    }
    view.set("parameter_coordinates",std::move(coordinates));
    const auto retained = fields.finish(std::move(view));
    std::vector<std::vector<learning::TransportContextView>> atlas(2);
    atlas[0].push_back({1,M::identity(1),M::identity(1)});
    for (std::size_t i = 0; i < 6; ++i) atlas[1].push_back({2,M::identity(2),M::identity(2)});
    learning::TransportPatternProvenance provenance;
    provenance.source_sha256 = digest(owner); provenance.evidence_sha256 = digest(retained);
    const learning::PreparedTransportPatternBankTraining prepared(evidence,{0,1,2,3,4,5},atlas,2,guard,
        learning::TransportContextInputScope::ObservedFields);
    check(prepared.acquisition().proof->evidence_kernel() == evidence.kernel(),"training omitted full conditioning kernel");
    auto bank = prepared.compile({0,1,3},{101,102,{100,103,104}},provenance,guard);
    check(bank.parameter_directions().cols() == 1 && evidence.kernel().cols() > 1,
          "observable quotient silently replaced full evidence uncertainty");
    check(evidence.system().design.multiply(evidence.kernel()).nonzeros() == 0,
          "retained exact conditioning directions do not satisfy training equations");
    return {std::move(bank),owner,retained};
}
Json map(std::initializer_list<std::array<std::string,3>> entries,std::size_t columns = 3,std::size_t rows = 2) {
    auto out = Json::object(), values = Json::array();
    out.set("rows",Json::integer(static_cast<std::int64_t>(rows)));
    out.set("columns",Json::integer(static_cast<std::int64_t>(columns)));
    for (const auto &e : entries) {
        auto entry = Json::array(); entry.push(Json::integer(std::stoll(e[0])));
        entry.push(Json::integer(std::stoll(e[1]))); entry.push(Json::string(e[2])); values.push(entry);
    }
    out.set("entries",values); return out;
}
Json input(std::size_t pattern,Json mapping,const std::string &field = "baseline.G") {
    auto out = Json::parse(R"({"frame":2,"domain_quantity":"synthetic_G","domain_units":"dimensionless","codomain_quantity":"synthetic_signal","codomain_units":"dimensionless","method":"declared synthetic coordinate selection","evidence":"explicit source fixture; no clinical claim"})");
    out.set("pattern",Json::integer(static_cast<std::int64_t>(pattern)));
    out.set("field",Json::string(field)); out.set("map",mapping); return out;
}
Json query(const NativePatternModel &model,const ScientificStudy &study,const std::string &operation) {
    auto out = Json::object(); out.set("schema",Json::string("gadget.pattern.query.v1"));
    out.set("operation",Json::string(operation)); out.set("source_digest",Json::string(study.digest()));
    out.set("model_sha256",Json::string(model.digest())); return out;
}
const storage::TemporalTensor &tensor(const NativeRecord &owner,const NativeRecord &reference) {
    for (const auto &t : owner.endpoint()->tensors)
        if (reference.at("tensor").integer() == t.key) return t;
    throw std::runtime_error("actual pattern result tensor absent");
}
bool mutate(N &value,const std::string &name,const N &replacement) {
    if (auto *o = std::get_if<N::Object>(&value.data)) for (auto &[key,v] : o->members) {
        if (key == name) { v = replacement; return true; }
        if (mutate(v,name,replacement)) return true;
    }
    if (auto *a = std::get_if<N::Array>(&value.data)) for (auto &v : a->elements)
        if (mutate(v,name,replacement)) return true;
    return false;
}
}
int main(int argc, char **argv) try {
    if (argc == 3 && std::string_view(argv[1]) == "--view-workflow") {
        const auto record = NativeRecord::read(argv[2]);
        const auto original = read_text(argv[2]);
        check(std::ranges::equal(record.encode(), std::as_bytes(std::span(original.data(), original.size()))),
              "complete native workflow archive changed on decode");
        std::cout << workflow_record_view(record).dump() << '\n';
        return 0;
    }
    check(argc == 1, "unexpected native pattern test arguments");
    const fs::path root = GADGET_PATTERN_TEST_BUILD;
    fs::remove_all(root); fs::create_directories(root);
    example::scientific_study(root/"source");
    auto study = std::make_shared<const ScientificStudy>(root/"source/study.json");
    study->write_native(root/"study.rfcxs");
    auto training = train();
    const auto declaration = admit_record_view(Json::parse(R"({"schema":"gadget.pattern.model.declaration.v1","scope":"Synthetic exact arithmetic and source ownership qualification","frames":[{"frame":1,"quantity":"synthetic_reading","units":"dimensionless"},{"frame":2,"quantity":"synthetic_signal","units":"dimensionless"}]})"));
    const auto publication = publish_pattern_model(training.bank,479,training.source,training.evidence,declaration);
    publication.write(root/"model.rfcxs"); const NativePatternModel model(NativeRecord::read(root/"model.rfcxs"));
    check(model.native_publication().encode() == publication.encode(),"cold model publication changed bytes");
    check(NativeRecord::decode(publication.at("evidence_publication").bytes()).encode() == training.evidence.encode(),
          "model lost actual original training evidence");
    rejects([&] { (void)publish_pattern_model(training.bank,479,declaration,training.evidence,declaration); },
            "model admitted an unrelated training source");
    auto candidates = query(model,*study,"sparse_candidates"), inputs = Json::array();
    inputs.push(input(0,map({{"0","0","1"},{"1","0","1"}})));
    inputs.push(input(1,map({{"0","0","1"}})));
    inputs.push(input(2,map({{"0","0","1"},{"1","0","1"}})));
    candidates.set("inputs",inputs); write_text(root/"candidates.json",candidates.dump());
    const auto result = model.execute_record(study,admit_record_view(candidates)); result.write(root/"candidates.rfcxs");
    const auto view = pattern_record_view(result).at("result").at("candidates");
    check(view.at(0).at("status").text() == "supported" && view.at(1).at("status").text() == "outside_domain" &&
              view.at(2).at("status").text() == "supported","candidate admissions lost evidenced domains");
    const auto values = result.at("result").at("fields").at("view").at("candidates");
    check(tensor(result,values.at(0).at("response").at("directions")).tensor.values().nonzeros() == 1 &&
              tensor(result,values.at(2).at("response").at("directions")).tensor.domain().logical_size() == 0,
          "candidate unknown and identified response were conflated");
    check(tensor(result,values.at(2).at("response").at("origin")).tensor.state(0) == Q{2},"anchor exact response changed");
    auto zero = query(model,*study,"sparse_candidates"), zero_inputs = Json::array();
    auto zero_input = input(0,map({{"0","0","1"},{"1","0","1"}}),"baseline.E");
    zero_input.set("domain_quantity",Json::string("synthetic_expression"));
    zero_input.set("domain_units",Json::string("synthetic_count")); zero_inputs.push(zero_input); zero.set("inputs",zero_inputs);
    const auto zero_result = model.execute_record(study,admit_record_view(zero));
    const auto zero_response = zero_result.at("result").at("fields").at("view").at("candidates").at(0).at("response");
    check(tensor(zero_result,zero_response.at("origin")).tensor.state(0) == Q{0} &&
              tensor(zero_result,zero_response.at("directions")).tensor.domain().logical_size() == 0,
          "supported zero source became outside domain or unknown");
    for (int which = 0; which < 6; ++which) {
        auto bad = Json::parse(candidates.dump());
        if (which == 0) bad.set("source_digest",Json::string("sha256:foreign"));
        if (which == 1) bad.set("model_sha256",Json::string("sha256:foreign"));
        if (which == 2) bad.at("inputs").at(0).set("codomain_units",Json::string("other"));
        if (which == 3) bad.at("inputs").at(0).set("frame",Json::integer(1));
        if (which == 4) bad.at("inputs").at(0).set("map",map({{"0","0","1"},{"0","0","2"}}));
        if (which == 5) { bad.at("inputs").at(0).set("field",Json::string("baseline.partial"));
            bad.at("inputs").at(0).set("domain_quantity",Json::string("synthetic_observation"));
            bad.at("inputs").at(0).set("map",map({{"0","2","1"}})); }
        rejects([&] { (void)model.execute_record(study,admit_record_view(bad)); },"invalid source/model/frame/units/support binding accepted");
    }
    const auto occurrence = [&](std::size_t slot) {
        return training.bank.endpoint().source.grade(core::Grade{1}).occurrence_identities[slot].value;
    };
    auto program = query(model,*study,"sparse_program"), nodes = Json::array();
    const auto root_node = [&](const char *id,std::size_t pattern,const char *gain) {
        auto node = Json::object(); node.set("id",Json::string(id)); node.set("pattern",Json::integer(static_cast<std::int64_t>(pattern)));
        node.set("input",input(pattern,map({{"0","0",gain},{"1","0",gain}}))); return node;
    };
    const auto edge = [&](const char *from,std::size_t slot) {
        auto out = Json::parse(R"({"source_frame":1,"target_frame":2,"method":"declared synthetic duplication","evidence":"original bank occurrence"})");
        out.set("from",Json::string(from)); out.set("occurrence_key",Json::integer(static_cast<std::int64_t>(occurrence(slot))));
        out.set("map",map({{"0","0","1"},{"1","0","1"}},1)); return out;
    };
    const auto stage = [&](const char *id,Json edges) {
        auto out = Json::object(); out.set("id",Json::string(id)); out.set("pattern",Json::integer(2)); out.set("edges",edges); return out;
    };
    nodes.push(root_node("first",0,"1")); nodes.push(root_node("second",1,"-1")); nodes.push(root_node("anchor",2,"1"));
    auto edges = Json::array(); edges.push(edge("first",0)); edges.push(edge("second",1)); edges.push(edge("anchor",2));
    nodes.push(stage("cancelled",edges)); edges = Json::array(); edges.push(edge("first",0)); nodes.push(stage("open",edges));
    edges = Json::parse(edges.dump()); edges.push(edge("open",2)); nodes.push(stage("blocked",edges)); program.set("nodes",nodes);
    write_text(root/"program.json",program.dump());
    const auto programmed = model.execute_record(study,admit_record_view(program)); programmed.write(root/"program.rfcxs");
    const auto stages = programmed.at("result").at("fields").at("view").at("stages");
    check(stages.at(3).at("status").text() == "supported" && stages.at(4).at("status").text() == "open_input" &&
              stages.at(5).at("status").text() == "blocked","program filled an unresolved input");
    check(tensor(programmed,stages.at(0).at("parameter_directions")).input == 2 &&
              tensor(programmed,stages.at(1).at("parameter_directions")).input == 2,
          "correlated stages lost shared bank parameter carrier");
    check(tensor(programmed,stages.at(3).at("input_directions")).tensor.stored_nonzeros() == 0 &&
              tensor(programmed,stages.at(3).at("response").at("origin")).tensor.state(2) == Q{4},
          "joint cancellation or exact composition changed");
    auto gluing = query(model,*study,"sparse_gluing"); inputs = Json::array(); inputs.push(candidates.at("inputs").at(0)); inputs.push(candidates.at("inputs").at(2));
    gluing.set("inputs",inputs); auto comparisons = Json::array();
    auto comparison = Json::parse(R"({"name":"anchor","left_frame":1,"right_frame":1,"interface_frame":1,"method":"declared synthetic comparison","evidence":"original context and anchor"})");
    comparison.set("left_occurrence_key",Json::integer(static_cast<std::int64_t>(occurrence(0))));
    comparison.set("right_occurrence_key",Json::integer(static_cast<std::int64_t>(occurrence(2))));
    comparison.set("left_map",map({{"0","0","1"}},1,1)); comparison.set("right_map",map({{"0","0","1"}},1,1));
    comparisons.push(comparison); gluing.set("comparisons",comparisons);
    const auto glued = model.execute_record(study,admit_record_view(gluing));
    check(pattern_record_view(glued).at("result").at("consistent").boolean(),"compatible exact patterns did not glue");
    comparison = Json::parse(comparison.dump()); comparison.set("name",Json::string("incompatible")); comparison.set("right_map",map({{"0","0","1/2"}},1,1));
    comparisons.push(comparison); gluing.set("comparisons",comparisons); write_text(root/"gluing.json",gluing.dump());
    const auto conflict = model.execute_record(study,admit_record_view(gluing)); conflict.write(root/"gluing.rfcxs");
    check(!pattern_record_view(conflict).at("result").at("consistent").boolean() &&
              conflict.at("result").at("fields").at("view").at("section").kind() == NativeRecordKind::Null,
          "contradiction fabricated a section");
    const auto cv = conflict.at("result").at("fields").at("view");
    const auto &defect = tensor(conflict,cv.at("defect_origin")).tensor;
    const auto &witness = tensor(conflict,cv.at("contradiction")).tensor;
    Q pairing{}; for (std::size_t i = 0; i < defect.domain().logical_size(); ++i)
        pairing += defect.state(i).value()*witness.state(i).value();
    check(pairing != 0,"contradiction omitted its separating witness");
    for (const auto &owner : {result,programmed,glued,conflict})
        check(pattern_record_view(replay_pattern_record(owner)).dump() == pattern_record_view(owner).dump(),"complete native pattern replay changed");
    auto changed = result.unbound_value(); check(mutate(changed,"status",N(std::string("outside_domain"))),"status tamper not exercised");
    rejects([&] { (void)replay_pattern_record(result.derive(std::move(changed))); },"replay accepted a changed result status");
    WorkflowModels callbacks;
    callbacks.emplace("pattern",WorkflowModel(1,[&](const std::string &,std::shared_ptr<CpuResources> resources) {
        return model.execute_record(study,admit_record_view(candidates),std::move(resources));
    }));
    auto workflow = Json::parse(R"({"schema":"gadget.scientific.workflow.v1","id":"native-pattern","limits":{"depth":8,"nodes":8,"model_calls":1,"parallel_workers":1},"definitions":[],"entry":{"id":"pattern","type":"model","role":"pattern","prompt":"Evaluate the declared original native sparse candidates"}})");
    workflow.set("source_digest",Json::string(study->digest())); write_text(root/"workflow.json",workflow.dump());
    auto steps = Json::array(); steps.push(workflow.at("entry"));
    auto condition = Json::parse(R"({"path":["candidates",2,"response","origin"],"coordinate":0,"value_type":"rational","comparison":"equal","value":"2"})");
    auto branch = Json::parse(R"({"id":"branch","type":"branch","from":"pattern","on_true":{"id":"observed","type":"query","query":{"operation":"coverage","field":"baseline.partial"}},"on_false":{"id":"other","type":"query","query":{"operation":"inspect"}}})");
    branch.set("condition",condition); steps.push(branch);
    auto claims = Json::array(); claims.push(condition);
    auto checking = Json::parse(R"({"id":"checked","type":"check_claims","from":"pattern"})");
    checking.set("claims",claims); steps.push(checking);
    auto sequence = Json::parse(R"({"id":"programme","type":"sequence"})"); sequence.set("steps",steps);
    workflow.set("entry",sequence); write_text(root/"workflow.json",workflow.dump());
    const auto receipt = execute_workflow_record(study,admit_record_view(workflow),callbacks);
    write_text(root/"workflow-view.json",workflow_record_view(receipt).dump());
    check(workflow_record_view(receipt).at("status").text() == "complete","native pattern callback failed workflow ownership admission");
    check(workflow_record_view(receipt).at("result").at("matched").integer() == 1,
          "native pattern claim did not read the original tensor coefficient");
    receipt.write(root/"workflow.rfcxs");
    auto roles = Json::parse(R"({"schema":"gadget.pattern.roles.v1","roles":[{"role":"pattern","model":"model.rfcxs","request":"candidates.json"}]})");
    write_text(root/"roles.json",roles.dump());
    fs::rename(root/"source/design.txt",root/"source/design-moved.txt");
    check(pattern_record_view(replay_pattern_record(NativeRecord::read(root/"program.rfcxs"))).dump() == pattern_record_view(programmed).dump(),
          "cold replay depended on an original input path");
    std::cout << checks << " native pattern training, fields, families, contradictions, workflows and replay checks PASS\n";
    return 0;
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
