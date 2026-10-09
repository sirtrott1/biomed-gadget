#include <atomic>
#include <gadget/scientific.hpp>
#include <gadget/workflow.hpp>
#include <iostream>
#include <thread>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected admitted study");
        auto owner = std::make_shared<const ScientificStudy>(argv[1]);
        ScientificSession session(owner);
        const auto index = owner->source_index();
        auto request = Json::object();
        request.set("schema", Json::string("gadget.scientific.source.read.v1"));
        request.set("source_digest", index.at("source_digest"));
        request.set("publication", Json::string("self"));
        request.set("source_sequence", index.at("source_sequence"));
        request.set("stream", Json::parse(R"({"kind":"field","id":"x"})"));
        request.set("offset", Json::integer(0));
        request.set("limit", Json::integer(2));
        auto check = [](bool condition, const char *message) {
            if (!condition)
                throw std::logic_error(message);
        };
        auto refuses = [&](auto action) {
            try {
                action();
            } catch (const std::invalid_argument &) {
                return;
            } catch (const rfc::runtime::ExecutionLimitExceeded &) {
                return;
            }
            throw std::logic_error("invalid source read accepted");
        };
        const auto expected = owner->read_source(request).dump();
        auto detached = owner->read_source(request);
        detached.at("records").at(0).set("value", Json::string("999"));
        check(owner->read_source(request).dump() == expected, "page mutated source");
        detached = owner->source_index();
        detached.set("source_digest", Json::string("changed"));
        check(owner->source_index().dump() == index.dump(), "index mutated source");
        const auto maximum = std::numeric_limits<std::size_t>::max();
        rfc::runtime::ExecutionHostPolicy policy;
        policy.name = "source-read-current-policy";
        policy.limits = {maximum, maximum, maximum, maximum, maximum};
        policy.limits.records = 2;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.limits.records = maximum;
        policy.limits.scalar_work = 1;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.limits.scalar_work = maximum;
        policy.limits.rational_slots = 0;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.limits.rational_slots = maximum;
        policy.limits.integer_bits = 0;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.limits.integer_bits = maximum;
        policy.allow_allocation = false;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.allow_allocation = true;
        policy.backend = rfc::runtime::ExecutionBackend::ExactGpu;
        refuses([&] { (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy)); });
        policy.fallback = rfc::runtime::BackendFallback::AllowExactCpuFallback;
        check(owner->read_source(request, rfc::runtime::ExecutionLease(policy))
                  .at("execution")
                  .at("fallback")
                  .boolean(),
              "missing declared fallback");
        auto cancelled = std::make_shared<std::atomic_bool>(true);
        try {
            (void)owner->read_source(request, rfc::runtime::ExecutionLease(policy, cancelled));
            throw std::logic_error("cancelled read accepted");
        } catch (const rfc::runtime::ExecutionCancelled &) {
        }
        auto proposal = Json::object();
        proposal.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
        proposal.set("action", Json::string("retrieve"));
        proposal.set("read", request);
        check(apply_scientific_proposal(session, proposal).at("certificate").dump() == expected,
              "model action changed source read");
        std::vector<std::string> results(4);
        std::vector<std::thread> workers;
        for (std::size_t i = 0; i < results.size(); ++i)
            workers.emplace_back([&, i] { results[i] = owner->read_source(request).dump(); });
        for (auto &worker : workers)
            worker.join();
        for (const auto &result : results)
            check(result == expected, "concurrent read changed source page");
        // A dependency can hand a verified page to a different model role.
        auto workflow = Json::parse(
            R"({"schema":"gadget.scientific.workflow.v1","id":"source-page-transfer","limits":{"depth":3,"nodes":3,"model_calls":2,"parallel_workers":1},"definitions":[],"entry":{"id":"sequence","type":"sequence","steps":[{"id":"reader","type":"model","role":"retriever","prompt":"Read the original coordinates"},{"id":"analyst","type":"model","role":"analyst","prompt":"Read support from the complete source","context":["reader"]}]}})");
        workflow.set("source_digest", Json::string(owner->digest()));
        WorkflowModels models;
        models.emplace("retriever", WorkflowModel{[&](const std::string &) {
                           return apply_scientific_proposal(session, proposal);
                       }});
        models.emplace(
            "analyst", WorkflowModel{[&](const std::string &prompt) {
                check(prompt.find(owner->read_source(request).at("page_sha256").text()) !=
                          std::string::npos,
                      "workflow lost original page witness");
                auto action = Json::parse(
                    R"({"schema":"gadget.scientific.agent.proposal.v1","action":"query","query":{"operation":"coverage","field":"x"}})");
                return apply_scientific_proposal(session, action);
            }});
        check(execute_workflow(owner, workflow, models).at("status").text() == "complete",
              "workflow rejected retrieved source page");
        models.at("retriever") = WorkflowModel{[&](const std::string &) {
            auto forged = apply_scientific_proposal(session, proposal);
            forged.at("certificate").at("records").at(0).set("value", Json::string("999"));
            return forged;
        }};
        check(execute_workflow(owner, workflow, models).at("status").text() == "failed",
              "workflow accepted forged page context");
        owner.reset();
        check(session.read_source(request).dump() == expected, "session lost owning source");
        std::cout << "source read ownership, current-policy, cancellation, concurrency and "
                     "proposal contracts passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
