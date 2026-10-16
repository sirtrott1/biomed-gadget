#include "../examples/scientific/fixture.hpp"
#include <gadget/model.hpp>
#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <future>
#include <iostream>
#include <thread>

namespace {
using namespace gadget;
void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
template<class E, class F> void rejects(F action, const char *message) {
    bool refused = false;
    try { action(); } catch (const E &) { refused = true; }
    check(refused, message);
}
}
int main(int argc, char **argv) {
    try {
        if (argc != 3) throw std::invalid_argument("model and test destination required");
        fs::create_directories(argv[2]);
        const fs::path root = private_directory(argv[2], "model-resources-");
        example::scientific_study(root / "input");
        const auto study = std::make_shared<const ScientificStudy>(root / "input/study.json");
        const auto hardware = rfc::runtime::detect_hardware();
        rfc::runtime::CpuResourcePolicy policy;
        policy.worker_limit = hardware.allowed_cpus.size();
        policy.declared_buffer_byte_limit = rfc::machine::checked_size(
            hardware.advisory_memory_headroom().value());
        const auto refresh = [] { return rfc::runtime::detect_hardware(); };
        const auto resources = std::make_shared<CpuResources>(hardware, policy, refresh);
        auto &coordinator = resources->coordinator();
        check(coordinator.capacity() >= 2, "real model overlap needs two admitted workers");
        LocalModelOptions options;
        options.model = argv[1]; options.resources = resources;
        options.context_tokens = 8192; options.output_tokens = 384;
        options.threads = static_cast<std::uint32_t>(std::min(std::size_t{8}, coordinator.capacity() - 1));
        options.prompt_cache_mib = 64;
        auto small_policy = policy;
        small_policy.declared_buffer_byte_limit = 1;
        const auto small = std::make_shared<CpuResources>(hardware, small_policy, refresh);
        options.resources = small;
        rejects<rfc::runtime::ExecutionLimitExceeded>([&] { LocalAgent refused(study, options); },
            "model loading ignored a caller byte limit");
        check(small->coordinator().available_workers() == small->coordinator().capacity() &&
                  small->coordinator().available_declared_buffer_bytes() == 1,
              "failed load leaked its workers or bytes");
        options.resources = resources;
        const auto before = coordinator.available_declared_buffer_bytes();
        auto receipt = Json::object();
        {
            LocalAgent agent(study, options);
            const auto idle = coordinator.available_declared_buffer_bytes();
            check(idle < before && before - idle >= fs::file_size(options.model),
                  "loaded weights/context have no retained byte charge");
            check(coordinator.available_workers() == coordinator.capacity(),
                  "idle model retains loading workers");
            receipt.set("model_retained_bytes", Json::string(std::to_string(before - idle)));
            {
                auto occupied = coordinator.reserve(1, idle);
                rejects<rfc::runtime::CpuResourceBusy>([&] { (void)coordinator.reserve(1, 1); },
                    "native work ignored resident model bytes");
            }
            const auto field = study->catalog().at("fields").at(0).at("id").text();
            const auto prompt = "Return action query with operation coverage and field " + field +
                ". Do not add other arguments or compute values yourself.";
            std::stop_source stop;
            auto task = std::async(std::launch::async, [&] {
                try { (void)agent.ask_assigned_record(prompt, {}, {}, stop.get_token()); }
                catch (const rfc::runtime::ExecutionCancelled &) { return true; }
                return false;
            });
            struct StopOnExit {
                std::stop_source &source;
                ~StopOnExit() { source.request_stop(); }
            } stop_on_exit{stop};
            const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
            while (coordinator.available_workers() == coordinator.capacity() &&
                   std::chrono::steady_clock::now() < deadline &&
                   task.wait_for(std::chrono::milliseconds(1)) != std::future_status::ready) {}
            if (coordinator.available_workers() == coordinator.capacity()) {
                stop.request_stop(); task.wait();
                throw std::runtime_error("real model call did not reserve its workers");
            }
            // Run an actual native observation while the neural call owns its workers.
            ScientificSession native(study, resources);
            auto query = Json::parse(R"({"operation":"coverage","field":""})");
            query.set("field", Json::string(field));
            const auto observation = native.execute_record(admit_record_view(query));
            observation.write(root / "overlap.rfcxs");
            check(coordinator.available_workers() == coordinator.capacity() - options.threads &&
                      coordinator.available_declared_buffer_bytes() <= idle,
                  "overlapping native/model ownership disappeared early");
            rejects<rfc::runtime::CpuResourceBusy>([&] {
                (void)coordinator.reserve(coordinator.capacity() - options.threads + 1);
            },
                "overlapping neural workers were overcommitted");
            const auto cancelled_at = std::chrono::steady_clock::now();
            stop.request_stop();
            check(task.get(), "cancelled neural call completed an unrequested proposal");
            receipt.set("cancel_join_seconds", Json::number(std::chrono::duration<double>(
                std::chrono::steady_clock::now() - cancelled_at).count()));
            check(coordinator.available_workers() == coordinator.capacity() &&
                      coordinator.available_declared_buffer_bytes() == idle,
                  "cancelled model leaked workers or released resident bytes");
            rejects<rfc::runtime::ExecutionCancelled>([&] {
                (void)agent.ask_assigned_record(prompt, {}, {}, stop.get_token());
            }, "pre-cancelled neural call ran");
            const auto result = agent.ask_record(prompt);
            result.write(root / "resumed.rfcxs");
            check(result.at("result").at("status").text() == "executed" &&
                      result.at("request").at("query").at("field").text() == field,
                  "neural context did not recover after cancellation");
            const auto direct = native.execute(query);
            check(json_object_equal(render_record_view(result.at("result").at("certificate")).raw(),
                                    direct.raw()),
                  "resumed neural calculation differs from native observation");
            receipt.set("model", render_record_view(result.at("result").at("model")));
            receipt.set("performance", render_record_view(result.at("result").at("performance")));
            const auto memory = result.at("result").at("model").at("memory");
            check(std::stoull(memory.at("llama_buffer_bound_bytes").text()) >=
                      std::stoull(memory.at("llama_reported_buffer_bytes").text()) &&
                      std::stoull(memory.at("declared_bytes").text()) == before - idle,
                  "model allocation report and retained reservation disagree");
            auto workflow = Json::parse(R"({"schema":"gadget.scientific.workflow.v1","id":"resident-model-native-graph","limits":{"depth":2,"nodes":3,"model_calls":0,"parallel_workers":2},"definitions":[],"entry":{"id":"parallel","type":"parallel","steps":[]}})");
            workflow.set("source_digest", Json::string(study->digest()));
            for (const auto *id : {"first", "second"}) {
                auto leaf = Json::object();
                leaf.set("id", Json::string(id)); leaf.set("type", Json::string("query"));
                leaf.set("query", query); workflow.at("entry").at("steps").push(leaf);
            }
            const auto graph = execute_workflow_record(study, admit_record_view(workflow), {}, resources);
            graph.write(root / "resident-graph.rfcxs");
            check(graph.at("result").at("status").text() == "complete",
                  "native graph cannot share remaining bytes with loaded model");
            write_text(root / "receipt.json", receipt.dump() + "\n");
        }
        check(coordinator.available_declared_buffer_bytes() == before &&
                  coordinator.available_workers() == coordinator.capacity(),
              "destroyed neural owner retains model bytes or workers");
        std::cout << "PASS: real model byte admission, native overlap, cancellation, resumed exact "
                     "execution and final owner release\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
