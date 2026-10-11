#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/file_io.hpp>
#include <iostream>

using namespace gadget;
namespace {
void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
const rfc::storage::TemporalTensor &field(const NativeRecord &record, const NativeRecord &key) {
    for (const auto &tensor : record.endpoint()->tensors)
        if (key.integer() == rfc::numeric::Integer(std::to_string(tensor.key))) return tensor;
    throw std::runtime_error("native field is absent");
}
void refuse(const std::function<void()> &action) {
    bool rejected = false;
    try {action();} catch (const std::exception &) {rejected = true;}
    check(rejected, "infeasible bounds were admitted as measurements");
}
}
int main(int argc, char **argv) {
    const fs::path root = argc > 1 ? argv[1] : "build/native-sequence-fields";
    fs::remove_all(root); fs::create_directories(root);
    const std::string reads = ">q0\r\nAC\r\n>q1\nAAA\n", reference = ">r0\nAC\n>r1\nAC\n";
    auto resources = Json::array();
    for (const auto &[id, bytes] : std::vector<std::pair<std::string, std::string>>{{"reads", reads}, {"reference", reference}}) {
        write_text(root / (id + ".fa"), bytes);
        auto resource = Json::parse(R"({"id":"","format":"fasta","source":{"path":"","sha256":""},"method":"Authored source fixture"})");
        resource.set("id", Json::string(id)); resource.at("source").set("path", Json::string(id + ".fa"));
        resource.at("source").set("sha256", Json::string("sha256:" + sha256(bytes))); resources.push(resource);
    }
    auto spec = Json::parse(R"({"schema":"gadget.assays.import.v1","id":"sequences","design":{"question":"Preserve ordered source witnesses","population":"Authored fixture","selection":"All source records","limitations":["Verification fixture"]},"assays":[],"resources":[],"links":[]})");
    spec.set("resources", resources); write_text(root / "import.json", spec.dump());
    import_assays(root / "import.json", root / "input");
    auto workspace = std::make_shared<const AssayWorkspace>(root / "input/workspace.rfcxs");
    auto alignment = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"exact_sequence_alignment","query_resource":"reads","query_id":"q0","reference_resource":"reference","reference_id":"r0","mode":"global_linear_gap","ambiguity":"literal","scoring":{"match":"0","mismatch":"2","query_gap":"0","reference_gap":"0","units":"declared cost","evidence":"Authored edit rule"}})");
    const auto aligned = workspace->execute_record(admit_record_view(alignment));
    const auto alignment_fields = aligned.at("result").at("result");
    const auto &complex = aligned.endpoint()->source;
    check(alignment_fields.at("schema").text() == "gadget.sequence.alignment.fields.v1" &&
          complex.grade({0}).cell_count() == 9 && complex.grade({1}).cell_count() == 16,
          "complete edit carrier was replaced");
    const auto &optimal = field(alignment_fields, alignment_fields.at("optimal_support")).tensor;
    check(optimal.support().count() == 16, "supported nonminimizing zero was omitted");
    check(alignment_fields.at("route_witnesses").size() > 0 &&
          alignment_fields.at("route_witnesses").at(0).at("minimum_cost").rational() == 0,
          "ordered minimum route was lost");
    const auto lift_key = alignment_fields.at("ordered_source_lift");
    const auto &lift = field(alignment_fields, lift_key);
    check(lift.role == rfc::storage::TemporalTensorRole::Observation && lift.output == 2,
          "ordered routes have no source lift");
    auto assembly = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"sequence_assembly_graph","resource":"reads","k":2,"ambiguity":"literal_tokens","evidence":"All supplied reads","route":{"start":"A","end":"C","objective":"fewest_steps"}})");
    const auto assembled = workspace->execute_record(admit_record_view(assembly));
    const auto assembly_fields = assembled.at("result").at("result");
    check(assembled.endpoint()->source.grade({1}).cell_count() == 2 &&
          assembled.endpoint()->source.grade({1}).incidence_count() == 4,
          "repeated self-loop boundary occurrences were discarded");
    const auto count_map_key = assembly_fields.at("count_map");
    const auto &count_map = field(assembly_fields, count_map_key);
    check(count_map.tensor.domain().extents() == std::vector<std::size_t>({2, 3}) &&
          count_map.tensor.support().count() == 6 && count_map.tensor.values().nonzeros() == 3,
          "actual assembly occurrence map was not retained");
    check(field(assembly_fields, assembly_fields.at("kmer_occurrence_count")).tensor.state(1).value() == 2 &&
          field(assembly_fields, assembly_fields.at("source_offset")).tensor.state(2).value() == 1 &&
          assembled.endpoint()->source.boundary_apply({1}, Vector{0, 1}) == Vector({0, 0}),
          "repeated read offsets or self-loop cancellation differ");
    auto constraints = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"constrained_allocations","domain":"bounded_integer","variables":[{"id":"x","lower":"0","upper":"2"},{"id":"y","lower":"0","upper":"2"}],"constraints":[{"id":"sum","coefficients":["1","1"],"sense":"eq","rhs":"2"}],"objective":{"sense":"minimize","coefficients":["0","0"]},"units":"counts","evidence":"Authored finite domain","witness_representation":"factorized"})");
    const auto constrained = workspace->execute_record(admit_record_view(constraints));
    const auto family = constrained.at("result").at("result");
    check(family.at("optimal_points_or_vertices").integer() == 3 && family.at("blocks").size() == 1 &&
          field(family, family.at("blocks").at(0).at("optimal_parameter_tuples")).tensor.domain().extents() ==
              std::vector<std::size_t>({1, 3}), "complete optimizing factor family was lost");
    check(field(family, family.at("tensors").at("objective")).tensor.support().count() == 2 &&
          field(family, family.at("tensors").at("objective")).tensor.values().nonzeros() == 0,
          "declared zero objective was confused with absence");
    auto expression = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"expression_allocations","samples":[{"id":"sample","resource":"reads","scale":"2","scale_provenance":"Authored positive scale"}],"evidence_source":{"mode":"transcript_sequences","reference_resource":"reference","strands":["+"],"ambiguity":"literal","maximum_mismatches":0,"evidence":"All declared transcript windows"},"targets":[{"id":"r0","candidates":["r0"]},{"id":"r1","candidates":["r1"]},{"id":"unassigned","candidates":["__unassigned__"]}],"contrasts":[],"map_evidence":"Authored targets","coupling":{"witness_representation":"factorized","evidence":"Authored count constraint","constraints":[{"id":"r0_count","terms":[{"sample":"sample","target":"r0","coefficient":"1"}],"sense":"eq","rhs":"1","evidence":"Authored allocation hypothesis"}]}})");
    const auto expressed = workspace->execute_record(admit_record_view(expression));
    const auto expression_fields = expressed.at("result").at("result"), coupling = expression_fields.at("coupling");
    check(coupling.at("status").text() == "feasible" && coupling.at("solvers").size() > 0,
          "native coupled family was not retained");
    for (std::size_t i = 0; i < coupling.at("solvers").size(); ++i) {
        const auto retained = coupling.at("solvers").at(i);
        const auto solver = NativeRecord::decode(retained.at("constraint_publication").bytes());
        check(solver.source_namespace() == retained.at("variable_binding").at("source_namespace").integer() &&
              solver.at("schema").text() == "gadget.constraints.fields.v1",
              "coupled solver source binding differs");
    }
    const auto target_coordinate = expression_fields.at("allocation_fields").at("target_coordinates").tensor().state(0)->get_num().get_ui();
    check(field(expression_fields, coupling.at("samples").at(0).at("lower")).tensor.state(target_coordinate).value() == 2,
          "coupled expression scale or typed optimum differs");
    std::vector<std::pair<std::string, Json>> recipes{{"alignment", alignment}, {"assembly", assembly},
                                                    {"constraints", constraints}, {"expression", expression}};
    std::vector<std::shared_ptr<const AssayWorkspace>> results;
    for (const auto &[id, recipe] : recipes) {
        auto request = Json::parse(R"({"schema":"gadget.assays.materialize.v1","id":"","source_namespace":"730010","channel":"observation","design":{},"recipe":{}})");
        request.set("id", Json::string(id)); request.set("source_namespace", Json::string(std::to_string(730010 + results.size())));
        request.set("design", spec.at("design")); request.set("recipe", recipe);
        auto result = workspace->materialize_workspace(admit_record_view(request), root / id);
        check(!fs::exists(root / id / "workspace.json"), "native publication created a JSON authority");
        AssayWorkspace reopened(root / id / "workspace.rfcxs");
        check(reopened.native_publication().encode() == result->native_publication().encode(),
              "native publication changed on reopen");
        results.push_back(std::move(result));
    }
    auto impossible = Json::parse(constraints.dump()); impossible.at("constraints").at(0).set("rhs", Json::string("5"));
    const auto failed = workspace->execute_record(admit_record_view(impossible)).at("result").at("result");
    check(failed.at("summary").at("status").text() == "infeasible" &&
          field(failed, failed.at("tensors").at("lower")).tensor.support().count() == 0,
          "infeasible bounds have numerical support");
    auto invalid = Json::parse(R"({"schema":"gadget.assays.materialize.v1","id":"invalid","source_namespace":"730020","channel":"observation","design":{},"recipe":{}})");
    invalid.set("design", spec.at("design")); invalid.set("recipe", impossible);
    refuse([&] {workspace->materialize_workspace(admit_record_view(invalid), root / "invalid");});
    fs::rename(root / "input", root / "input-away"); workspace.reset();
    for (std::size_t i = 0; i < results.size(); ++i)
        check(results[i]->replay_derivation(recipes[i].first, root / (recipes[i].first + "-replay")).at("status").text() == "PASS",
              "native replay depends on the former input directory");
    std::cout << "PASS native ordered routes, occurrence maps, optimizing factors, coupled owners and relocated replay\n";
}
