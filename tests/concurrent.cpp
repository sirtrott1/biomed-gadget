#include <barrier>
#include <gadget/model.hpp>
#include <gadget/store.hpp>
#include <iostream>
#include <sys/wait.h>
#include <thread>

namespace {
using namespace gadget;
namespace native = rfc::storage;
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&f, const char *message) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, message);
}
Json request(const Json &root, std::uint64_t transaction) {
    auto value = Json::object();
    value.set("schema", Json::string("gadget.store.publish.v1"));
    value.set("transaction", Json::string(std::to_string(transaction)));
    value.set("expected", Json::parse(root.dump()));
    return value;
}
Json proposal(const Json &query) {
    auto value = Json::object();
    value.set("schema", Json::string("gadget.agent.proposal.v1"));
    value.set("action", Json::string("query"));
    value.set("query", query);
    return value;
}
} // namespace
int main() {
    try {
        const fs::path root = GADGET_CONCURRENT_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root);
        fs::copy(GADGET_TEST_FIXTURE, root / "input", fs::copy_options::recursive);
        auto first = load_pager(root / "input");
        import_measurements(first, root / "input/measurements.json");
        save_snapshot(first, root / "original", true);
        write_text(root / "input/effects.csv", "gene,log2FC,padj\nA,0.25,1e-12\nB,1,NA\n");
        auto second = load_pager(root / "input");
        import_measurements(second, root / "input/measurements.json");
        save_snapshot(second, root / "updated", true);
        const auto genesis = initialize_store(root / "original", root / "store");
        const auto query = Json::parse(
            R"({"operation":"typed_moment","context":"AD","realization":"support","field":"AD.log2FC","other_field":"AD.G","order":1})");
        std::optional<StudySnapshot> retained;
        {
            ConcurrentStudyStore service(root / "store");
            const auto initial = service.snapshot();
            retained = initial;
            QuerySession old_session(initial.workspace());
            const auto old_result = old_session.execute(query);
            check(old_result.dump() == execute_query(*initial.workspace(), query).dump(),
                  "prepared query changed its certificate");
            old_session.execute(query);
            check(old_session.preparation_count() == 1, "preparation was repeated");
            auto invalid = Json::parse(query.dump());
            invalid.set("other_field", Json::string("HD.G"));
            check(apply_model_proposal(old_session, proposal(invalid)).at("status").text() ==
                      "rejected",
                  "invalid model query executed");
            auto hostile = proposal(
                Json::parse(R"json({"operation":"store-publish","code":"system('false')"})json"));
            check(apply_model_proposal(old_session, hostile).at("status").text() == "rejected",
                  "model mutation/code executed");
            auto clarify = Json::parse(
                R"({"schema":"gadget.agent.proposal.v1","action":"clarify","reason":"Which measured background?"})");
            check(apply_model_proposal(old_session, clarify).at("status").text() == "clarification",
                  "clarification rejected");
            hostile = proposal(query);
            hostile.set("code", Json::string("false"));
            check(apply_model_proposal(old_session, hostile).at("status").text() == "rejected",
                  "unknown envelope member executed");
            std::atomic<unsigned> read_errors{}, reads{}, committed{}, conflicts{};
            std::atomic<bool> done{};
            std::barrier readers_ready(5);
            std::vector<std::jthread> readers;
            for (int worker = 0; worker < 4; ++worker)
                readers.emplace_back([&] {
                    readers_ready.arrive_and_wait();
                    try {
                        do {
                            const auto pin = service.snapshot();
                            const auto study = pin.workspace();
                            const auto sequence = pin.root().sequence;
                            check(sequence <= 1 &&
                                      study->fields == (sequence ? second.fields : first.fields),
                                  "reader saw mixed study state");
                            check(study->store_origin.at("root").at("sequence").text() ==
                                      std::to_string(sequence),
                                  "reader root disagreed with study");
                            QuerySession private_session(study);
                            private_session.execute(query);
                            ++reads;
                        } while (!done.load());
                    } catch (const std::exception &error) {
                        std::cerr << "reader: " << error.what() << '\n';
                        ++read_errors;
                    }
                });
            readers_ready.arrive_and_wait();
            std::barrier start(2);
            {
                std::vector<std::jthread> writers;
                for (std::uint64_t transaction : {101, 102})
                    writers.emplace_back([&, transaction] {
                        start.arrive_and_wait();
                        try {
                            service.publish(initial, root / "updated",
                                            request(genesis.at("root"), transaction));
                            ++committed;
                        } catch (const native::NativeHistoryConflict &) {
                            ++conflicts;
                        } catch (const std::exception &error) {
                            std::cerr << "writer: " << error.what() << '\n';
                            ++read_errors;
                        }
                    });
            }
            done = true;
            readers.clear();
            check(read_errors == 0 && reads >= 4 && committed == 1 && conflicts == 1,
                  "reader/writer race contract failed");
            check(old_session.execute(query).dump() == old_result.dump(),
                  "old query session changed after publication");
            auto latest = service.snapshot();
            QuerySession new_session(latest.workspace());
            check(new_session.execute(query).at("result").dump() != old_result.at("result").dump(),
                  "new query session reused an old result");
            const auto winner = service.status(101).at("status").text() == "committed" ? 101 : 102;
            std::stop_source cancelled;
            cancelled.request_stop();
            const auto retry =
                service.publish(initial, root / "updated", request(genesis.at("root"), winner),
                                cancelled.get_token());
            check(retry.at("receipt").at("root").at("sequence").text() == "1",
                  "cancelled historical retry did not resolve");
            rejects(
                [&] {
                    service.publish(latest, root / "updated",
                                    request(service.status().at("root"), 103),
                                    cancelled.get_token());
                },
                "cancelled new publication committed");
            check(service.status(103).at("status").text() == "absent",
                  "cancelled new transaction indexed");
            initialize_store(root / "original", root / "foreign");
            ConcurrentStudyStore foreign(root / "foreign");
            rejects(
                [&] {
                    service.publish(foreign.snapshot(), root / "updated",
                                    request(genesis.at("root"), 104));
                },
                "foreign token admitted");
            // Reentry rejects before the recursive application gate can deadlock.
            const auto equal_request = request(service.status().at("root"), 105);
            service.publish(latest, root / "updated", equal_request, {}, [&](auto event) {
                if (event.stage == native::DurableStage::BeforeHeadReplace)
                    rejects([&] { service.status(); }, "writer callback reentry admitted");
                if (event.stage == native::DurableStage::LivePublished) {
                    const auto pin = service.snapshot();
                    check(pin.root().sequence == 2 &&
                              pin.workspace()->store_origin.at("root").at("sequence").text() ==
                                  "2" &&
                              pin.workspace()->fields == second.fields,
                          "publication hook did not see the complete new study");
                }
            });
            const auto before_unknown = service.snapshot();
            const auto unknown_request = request(service.status().at("root"), 106);
            bool unknown = false;
            try {
                service.publish(before_unknown, root / "updated", unknown_request, {},
                                [&](auto event) {
                                    if (event.stage == native::DurableStage::HeadReplaced)
                                        throw std::runtime_error("lost acknowledgement");
                                });
            } catch (const native::DurableUnknownOutcome &error) {
                unknown = error.transaction() == 106;
            }
            check(unknown && service.snapshot().root() == before_unknown.root(),
                  "unknown outcome damaged owning readers");
            rejects([&] { service.status(106); }, "unknown status was resolved without recovery");
            const auto manifest_path = root / "store/studies/0000000000000000/snapshot.json";
            const auto manifest_bytes = read_text(manifest_path);
            write_text(manifest_path, manifest_bytes + "corrupt");
            rejects([&] { service.recover(); }, "corrupt recovery succeeded");
            check(service.snapshot().root() == before_unknown.root(),
                  "failed recovery replaced old reader");
            write_text(manifest_path, manifest_bytes);
            service.recover();
            check(service.status(106).at("status").text() == "committed" &&
                      service.snapshot().root().sequence == 3,
                  "explicit recovery did not resolve publication");
            service.publish(before_unknown, root / "updated", unknown_request);
            const auto child = ::fork();
            check(child >= 0, "fork failed");
            if (child == 0) {
                try {
                    service.snapshot();
                } catch (const std::logic_error &) {
                    ::_exit(73);
                }
                ::_exit(74);
            }
            int status{};
            check(::waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                      WEXITSTATUS(status) == 73,
                  "inherited live service use admitted");
        }
        QuerySession survivor(retained->workspace());
        check(survivor.execute(query).at("store").at("root").at("sequence").text() == "0",
              "reader did not survive service destruction");
        std::cout
            << "PASS: concurrent owning studies, writer conflict, query preparation reuse, model "
               "admission, retry/cancellation, reentry/fork, unknown and failed/successful "
               "recovery\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
