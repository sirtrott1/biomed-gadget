#include "../examples/scientific/fixture.hpp"
#include <gadget/resources.hpp>
#include <gadget/workflow.hpp>
#include <iostream>

namespace {
using namespace gadget;
std::size_t checks{};
void check(bool value) {
    if (!value)
        throw std::runtime_error("workflow source/budget contract " + std::to_string(checks));
    ++checks;
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
