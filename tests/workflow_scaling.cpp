#include "../examples/scientific/fixture.hpp"
#include <gadget/resources.hpp>
#include <gadget/workflow.hpp>
#include <gadget/record_view.hpp>
#include <iostream>
#include <source_location>

namespace {
using namespace gadget;
std::size_t checks{};
void check(bool value, std::source_location location = std::source_location::current()) {
    if (!value)
        throw std::runtime_error("workflow source/budget contract " + std::to_string(checks) +
                                 " at line " + std::to_string(location.line()));
    ++checks;
}
template<class F> void rejects(F f) {
    bool refused = false;
    try { f(); } catch (const std::invalid_argument &) { refused = true; }
    check(refused);
}
Json leaf(const std::string &id) {
    auto node = Json::parse(
        R"({"id":"","type":"query","query":{"operation":"field_summary","field":"baseline.G","support":"reported_only"}})");
    node.set("id", Json::string(id));
    return node;
}
Json declaration(const ScientificStudy &study, Json entry, std::size_t nodes, std::size_t calls = 0,
                 std::size_t workers = 1) {
    auto d = Json::parse(
        R"({"schema":"gadget.scientific.workflow.v1","id":"source-sized-workflow","limits":{"depth":64,"nodes":1,"model_calls":0,"parallel_workers":1},"definitions":[]})");
    d.set("source_digest", Json::string(study.digest()));
    d.set("entry", entry);
    d.at("limits").set("nodes", Json::integer(static_cast<std::int64_t>(nodes)));
    d.at("limits").set("model_calls", Json::integer(static_cast<std::int64_t>(calls)));
    d.at("limits").set("parallel_workers", Json::integer(static_cast<std::int64_t>(workers)));
    return d;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::invalid_argument("test output directory required");
        const fs::path root = argv[1];
        fs::remove_all(root);
        example::scientific_study(root);
        auto study = std::make_shared<const ScientificStudy>(root / "study.json");
        auto sequence = Json::parse(R"({"id":"sequence","type":"sequence","steps":[]})");
        for (std::size_t i = 0; i < 145; ++i)
            sequence.at("steps").push(leaf(std::string(160, 'x') + std::to_string(i)));
        auto d = declaration(*study, sequence, 146);
        check(execute_workflow(study, d).at("status").text() == "complete");
        d.at("limits").set("nodes", Json::integer(145));
        auto exhausted = execute_workflow(study, d);
        check(exhausted.at("status").text() == "failed");
        check(exhausted.at("error").text().find("node invocation") != std::string::npos);
        std::size_t invoked{};
        WorkflowModels models;
        for (std::size_t i = 0; i < 20; ++i)
            models.emplace(
                "reader" + std::to_string(i), WorkflowModel{[&](const std::string &prompt) {
                    check(prompt.size() > 8192);
                    ++invoked;
                    return Json::parse(
                        R"({"schema":"gadget.agent.response.v1","status":"clarification"})");
                }});
        auto calls = Json::parse(R"({"id":"readers","type":"sequence","steps":[]})");
        for (std::size_t i = 0; i < 130; ++i) {
            auto node = Json::parse(R"({"id":"","type":"model","role":"","prompt":""})");
            node.set("id", Json::string("call" + std::to_string(i)));
            node.set("role", Json::string("reader" + std::to_string(i % 20)));
            node.set("prompt", Json::string(std::string(9000, 'x')));
            calls.at("steps").push(node);
        }
        d = declaration(*study, calls, 131, 130);
        check(d.dump().size() > 1024u * 1024u);
        check(execute_workflow(study, d, models).at("status").text() == "complete");
        check(invoked == 130);
        d.at("limits").set("model_calls", Json::integer(129));
        check(execute_workflow(study, d, models).at("status").text() == "failed");
        check(invoked == 259);
        auto graph = Json::parse(R"({"id":"graph","type":"graph","steps":[]})");
        for (std::size_t i = 0; i < 25; ++i) {
            auto wrapper = Json::parse(R"({"id":"","after":[],"node":{}})");
            wrapper.set("id", Json::string("job" + std::to_string(i)));
            wrapper.set("node", leaf("job" + std::to_string(i)));
            if (i == 24)
                for (std::size_t j = 0; j < 24; ++j)
                    wrapper.at("after").push(Json::string("job" + std::to_string(j)));
            graph.at("steps").push(wrapper);
        }
        auto resources = CpuResources::local();
        const auto workers = std::min(std::size_t{5}, resources->coordinator().capacity());
        d = declaration(*study, graph, 26, 0, workers);
        graph.set("workers", Json::integer(static_cast<std::int64_t>(workers)));
        d.set("entry", graph);
        check(execute_workflow(study, d, {}, resources).at("status").text() == "complete");
        auto call = Json::parse(
            R"({"id":"recur","type":"call","definition":"loop","budget":20,"on_exhausted":{}})");
        call.set("on_exhausted", leaf("base"));
        d = declaration(*study, call, 64);
        auto def = Json::object();
        def.set("name", Json::string("loop"));
        def.set("node", call);
        d.at("definitions").push(def);
        for (std::size_t i = 0; i < 20; ++i) {
            auto unused = Json::object();
            unused.set("name", Json::string("unused" + std::to_string(i)));
            unused.set("node", leaf("unused-leaf"));
            d.at("definitions").push(unused);
        }
        check(execute_workflow(study, d).at("status").text() == "complete");
        auto nested = leaf("bottom");
        for (std::size_t i = 0; i < 18; ++i) {
            auto level = Json::parse(R"({"id":"level","type":"sequence","steps":[]})");
            level.at("steps").push(nested);
            nested = level;
        }
        d = declaration(*study, nested, 19);
        check(execute_workflow(study, d).at("status").text() == "complete");
        d.at("limits").set("depth", Json::integer(17));
        bool refused = false;
        try {
            (void)execute_workflow(study, d);
        } catch (const std::exception &) {
            refused = true;
        }
        check(refused);
        ScientificSession session(study);
        const auto model = session.execute_record(admit_record_view(Json::parse(
            R"({"operation":"linear_model_comparison","responses":["baseline.partial"],"terms":[{"id":"G","factors":["baseline.G"]}],"controls":[],"selection":[],"support":"joint_reported","metric":"identity","intercept":true})")));
        auto native_claims = Json::parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[
          {"path":["strata",0,"models",1,"responses",0,"residual_quadrance"],"value_type":"rational","comparison":"equal","value":"0"},
          {"path":["strata",0,"models",1,"responses",0,"represented"],"coordinate":0,"value_type":"rational","comparison":"equal","value":"0"},
          {"path":["strata",0,"models",1,"responses",0,"represented"],"coordinate":1,"value_type":"rational","comparison":"greater","value":"1"},
          {"path":["strata",0,"models",1,"rank"],"value_type":"count","comparison":"equal","value":3}
        ]})");
        native_claims.set("source_namespace", Json::string(std::to_string(model.source_namespace())));
        native_claims.set("cutoff", Json::integer(static_cast<std::int64_t>(model.publication_history()->last_sequence())));
        const auto request = admit_record_view(native_claims);
        const auto checked = check_certificate_claims_record(model, request);
        const auto result = checked.at("result");
        check(result.at("matched").integer() == 3 && result.at("contradicted").integer() == 1);
        check(result.at("checks").at(1).at("observed").rational() == 0);
        const auto pin = result.at("checks").at(1).at("binding");
        check(pin.at("coordinate").integer() == 0 && pin.at("supported").boolean());
        check(pin.at("source_namespace").text() == std::to_string(model.source_namespace()));
        check(pin.at("output_carrier").integer() == model.at("result").at("fields").at("bindings")
            .at(0).at("source_carrier").integer());
        const auto original = NativeRecord::decode(result.at("source_publication").bytes());
        check(original.encode() == model.encode());
        checked.write(root / "native-claims.rfcxs");
        const auto reopened = NativeRecord::read(root / "native-claims.rfcxs");
        check(reopened.at("result").at("checks").at(1).at("observed").rational() == 0);
        check(NativeRecord::decode(reopened.at("request_publication").bytes()).encode() == request.encode());
        const auto replayed = replay_scientific_record(original);
        check(check_certificate_claims_record(replayed, admit_record_view([&] {
            auto q = Json::parse(native_claims.dump());
            q.set("source_namespace", Json::string(std::to_string(replayed.source_namespace())));
            q.set("cutoff", Json::integer(static_cast<std::int64_t>(replayed.publication_history()->last_sequence())));
            return q;
        }())).at("result").at("matched").integer() == 3);
        auto changed = Json::parse(native_claims.dump());
        changed.set("source_namespace", Json::string(std::to_string(study->native_publication().source_namespace())));
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed = Json::parse(native_claims.dump());
        changed.set("cutoff", Json::integer(changed.at("cutoff").integer() + 1));
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed = Json::parse(native_claims.dump());
        changed.at("items").at(1).set("coordinate", Json::integer(2));
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed.at("items").at(1).set("coordinate", Json::integer(3));
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed = Json::parse(native_claims.dump());
        json_object_object_del(changed.at("items").at(1).raw(), "coordinate");
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed = Json::parse(native_claims.dump());
        changed.at("items").at(0).set("value_type", Json::string("text"));
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        changed = Json::parse(native_claims.dump());
        changed.set("items", Json::array());
        rejects([&] { (void)check_certificate_claims_record(model, admit_record_view(changed)); });
        const auto legacy = session.execute_record(admit_record_view(Json::parse(
            R"({"operation":"field_summary","field":"baseline.G","support":"reported_only"})")));
        changed = Json::parse(native_claims.dump());
        changed.set("source_namespace", Json::string(std::to_string(legacy.source_namespace())));
        changed.set("cutoff", Json::integer(static_cast<std::int64_t>(legacy.publication_history()->last_sequence())));
        rejects([&] { (void)check_certificate_claims_record(legacy, admit_record_view(changed)); });
        auto native_leaf = [&](const std::string &id) {
            auto node = Json::object();
            node.set("id", Json::string(id)); node.set("type", Json::string("query"));
            node.set("query", render_record_view(model.at("request")));
            return node;
        };
        auto native_graph = Json::parse(R"({"id":"native-graph","type":"graph","workers":2,"steps":[
            {"id":"fit1","after":[],"node":{}},
            {"id":"fit2","after":[],"node":{}},
            {"id":"verified","after":["fit1","fit2"],"node":{"id":"verified","type":"check_claims","from":"fit1","claims":[{"path":["strata",0,"models",1,"rank"],"value_type":"count","comparison":"equal","value":99}]}},
            {"id":"repair","after":["verified"],"node":{"id":"repair","type":"branch","from":"verified","condition":{"path":["contradicted"],"value_type":"count","comparison":"greater","value":0},"on_true":{"id":"recur","type":"call","definition":"native-loop","budget":2,"on_exhausted":{}},"on_false":{}}}
        ]})");
        native_graph.at("steps").at(0).set("node", native_leaf("fit1"));
        native_graph.at("steps").at(1).set("node", native_leaf("fit2"));
        auto native_call = native_graph.at("steps").at(3).at("node").at("on_true");
        native_call.set("on_exhausted", native_leaf("native-base"));
        native_graph.at("steps").at(3).at("node").set("on_false", native_leaf("unused-native"));
        auto native_declaration = declaration(*study, native_graph, 32, 0, 2);
        auto native_definition = Json::object();
        native_definition.set("name", Json::string("native-loop"));
        native_definition.set("node", native_call);
        native_declaration.at("definitions").push(native_definition);
        const auto native_request = admit_record_view(native_declaration);
        const auto graph_operation = execute_workflow_record(study, native_request);
        const auto native_receipt = graph_operation.at("result");
        check(native_receipt.at("native_execution").boolean() && native_receipt.at("status").text() == "complete");
        const auto first = NativeRecord::decode(native_receipt.at("result").at(0).at("publication").bytes());
        const auto second = NativeRecord::decode(native_receipt.at("result").at(1).at("publication").bytes());
        check(first.source_namespace() != second.source_namespace());
        check(first.at("result").at("fields").at("view").at("strata").at(0).at("models")
            .at(1).at("coefficient_nullity").integer() == 0);
        check(NativeRecord::decode(first.at("request_publication").bytes()).encode() == native_request.encode());
        bool found_native_base{}, found_native_decision{}, found_work_graph{};
        for (std::size_t i = 0; i < native_receipt.at("trace").size(); ++i) {
            const auto event = native_receipt.at("trace").at(i);
            found_native_base |= event.has("base_case") && event.at("base_case").boolean();
            found_work_graph |= event.has("native_work");
            if (event.at("node_id").text() == "fit1")
                check(event.at("result").at("publication").byte_source().key() ==
                    native_receipt.at("result").at(0).at("publication").byte_source().key());
            if (event.at("node_id").text() == "repair") {
                found_native_decision = event.at("branch_taken").text() == "on_true";
                const auto comparison_owner = NativeRecord::decode(event.at("comparison").at("publication").bytes());
                check(comparison_owner.at("result").at("checks").at(0).at("observed").integer() == 1);
            }
        }
        check(found_native_base && found_native_decision && found_work_graph);
        graph_operation.write(root / "owning-workflow.rfcxs");
        const auto graph_reopened = NativeRecord::read(root / "owning-workflow.rfcxs");
        check(NativeRecord::decode(graph_reopened.at("result").at("result").at(0).at("publication").bytes()).encode() == first.encode());
        check(workflow_record_view(graph_reopened).at("status").text() == "complete");
        (void)replay_scientific_record(first);
        auto native_sequence = Json::parse(R"({"id":"coordinate-branch","type":"sequence","steps":[{},
            {"id":"select-coordinate","type":"branch","from":"coordinate-fit","condition":{"path":["strata",0,"models",1,"responses",0,"represented"],"coordinate":1,"value_type":"rational","comparison":"greater","value":"1"},"on_true":{},"on_false":{}}]})");
        auto native_steps = Json::array();
        native_steps.push(native_leaf("coordinate-fit"));
        auto coordinate_branch = native_sequence.at("steps").at(1);
        coordinate_branch.set("on_true", native_leaf("coordinate-true"));
        coordinate_branch.set("on_false", native_leaf("coordinate-false"));
        native_steps.push(coordinate_branch); native_sequence.set("steps", native_steps);
        auto coordinate_workflow = declaration(*study, native_sequence, 8);
        auto coordinate_operation = execute_workflow_record(study, admit_record_view(coordinate_workflow));
        check(workflow_record_view(coordinate_operation).at("status").text() == "complete");
        coordinate_workflow.at("entry").at("steps").at(1).at("condition").set("coordinate", Json::integer(2));
        const auto missing_operation = execute_workflow_record(study, admit_record_view(coordinate_workflow));
        check(missing_operation.at("result").at("status").text() == "failed" &&
            missing_operation.at("result").at("error").text().find("missing") != std::string::npos);
        coordinate_workflow.at("entry").at("steps").at(1).set("on_false", leaf("unsupported-untaken"));
        rejects([&] { (void)execute_workflow_record(study, admit_record_view(coordinate_workflow)); });
        auto model_node = Json::parse(R"({"id":"native-reader","type":"model","role":"native-reader","prompt":"Read the retained native tool result"})");
        auto native_callback_declaration = declaration(*study, model_node, 1, 1);
        WorkflowModels native_callbacks{{"native-reader", WorkflowModel(1, WorkflowModel::NativeInvoke{
            [&](const std::string &, std::shared_ptr<CpuResources>) { return model; }})}};
        const auto callback_operation = execute_workflow_record(study, admit_record_view(native_callback_declaration), native_callbacks);
        check(callback_operation.at("result").at("status").text() == "complete");
        check(NativeRecord::decode(callback_operation.at("result").at("result").at("publication").bytes()).encode() == model.encode());
        WorkflowModels rendered_callbacks{{"native-reader", [&](const std::string &) {
            auto response = Json::object(); response.set("schema", Json::string("gadget.agent.response.v1"));
            response.set("status", Json::string("executed")); response.set("certificate", Json::object());
            return response;
        }}};
        const auto rendered_refusal = execute_workflow_record(study, admit_record_view(native_callback_declaration), rendered_callbacks);
        check(rendered_refusal.at("result").at("status").text() == "failed");
        check(workflow_record_view(rendered_refusal).at("trace").at(0).at("result").at("status").text() == "executed");
        auto model_claim_graph = Json::parse(R"({"id":"native-claims-graph","type":"graph","workers":2,"steps":[
            {"id":"fit","after":[],"node":{}},
            {"id":"proposal","after":["fit"],"node":{"id":"proposal","type":"model","role":"claimant","prompt":"Propose a literal claim","context":["fit"]}},
            {"id":"checked","after":["fit","proposal"],"node":{"id":"checked","type":"check_claims","from":"fit","claims_from":"proposal"}}
        ]})");
        model_claim_graph.at("steps").at(0).set("node", native_leaf("fit"));
        auto model_claim_workflow = declaration(*study, model_claim_graph, 8, 1, 2);
        bool wrong_native_pin{};
        WorkflowModels native_claim_models{{"claimant", [&](const std::string &prompt) {
            const auto marker = prompt.find("Certified data context");
            const auto context = Json::parse(prompt.substr(prompt.find('[', marker))).at(0);
            auto supplied = Json::parse(R"({"schema":"gadget.claims.native.v1","source_namespace":"","cutoff":0,"items":[{"path":["strata",0,"models",1,"rank"],"value_type":"count","comparison":"equal","value":99}]})");
            supplied.set("source_namespace", wrong_native_pin ? Json::string("wrong-owner") : context.at("native_source_namespace"));
            supplied.set("cutoff", context.at("native_cutoff"));
            auto response = Json::object(), certificate = Json::object();
            response.set("schema", Json::string("gadget.agent.response.v1")); response.set("status", Json::string("prospective"));
            response.set("source_digest", context.at("source_digest")); certificate.set("source_digest", context.at("source_digest"));
            response.set("certificate", certificate); response.set("claims", supplied); return response;
        }}};
        const auto native_model_claims = execute_workflow_record(study, admit_record_view(model_claim_workflow), native_claim_models);
        check(native_model_claims.at("result").at("status").text() == "complete");
        check(NativeRecord::decode(native_model_claims.at("result").at("result").at(2).at("publication").bytes())
            .at("result").at("contradicted").integer() == 1);
        wrong_native_pin = true;
        check(execute_workflow_record(study, admit_record_view(model_claim_workflow), native_claim_models)
            .at("result").at("status").text() == "failed");
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks
                  << ",\"node_leaves\":145,\"model_calls\":130,\"mock_roles\":20,\"dependencies\":"
                     "24,\"recursion_budget\":20,\"workers\":"
                  << workers << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
