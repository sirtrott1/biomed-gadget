#include "native_transition_fields.hpp"
#include <atomic>
#include <chrono>
#include <gadget/assay.hpp>
#include <gadget/resources.hpp>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc < 2 || argc > 4)
            throw std::invalid_argument("expected novelty case, optional repetitions and mode");
        const fs::path root = argv[1];
        if (argc == 3 && std::string(argv[2]) == "native") {
            native_transition_test::run(root, true); return 0;
        }
        auto source = std::make_shared<const ScientificStudy>(root / "study.json");
        ScientificSession session(source);
        auto query = Json::parse(read_text(root / "query.json"));
        const auto stats = [&] { return session.native_novelty_preparation_statistics(); };
        if (argc == 4 && std::string(argv[2]) == "workspace") {
            const fs::path manifest = argv[3];
            AssayWorkspace workspace(manifest);
            const auto request =
                Json::parse(read_text(manifest.parent_path().parent_path() / "query.json"));
            const auto cpu = CpuResources::local();
            const auto first = workspace.execute(request, cpu).dump();
            if (workspace.execute(request, cpu).dump() != first)
                throw std::logic_error("workspace factor reuse changed receipt");
            auto before = workspace.native_novelty_preparation_statistics();
            if (before.size() != 1 ||
                before.at(0).at("preparation").at("preparations").integer() != 1)
                throw std::logic_error("workspace did not share the actual native plan");
            const auto hardware = rfc::runtime::detect_hardware();
            rfc::runtime::CpuResourcePolicy separate_policy;
            separate_policy.worker_limit = hardware.allowed_cpus.size();
            separate_policy.tile_limit = std::numeric_limits<std::size_t>::max();
            separate_policy.declared_buffer_byte_limit =
                rfc::machine::checked_size(hardware.advisory_memory_headroom().value());
            auto separate = std::make_shared<CpuResources>(
                hardware, separate_policy, [] { return rfc::runtime::detect_hardware(); });
            (void)workspace.execute(request, separate);
            if (workspace.native_novelty_preparation_statistics().size() != 2)
                throw std::logic_error("independent CPU owners were identified");
            workspace.clear_native_novelty_preparations();
            if (workspace.native_novelty_preparation_statistics().size() != 0)
                throw std::logic_error("workspace clear retained session factors");
            if (workspace.execute(request, cpu).dump() != first)
                throw std::logic_error("workspace clear changed source or result");
            std::cout << workspace.native_novelty_preparation_statistics().dump() << '\n';
            return 0;
        }
        if (argc >= 3) {
            const auto repeats = std::stoull(argv[2]);
            if (!repeats)
                throw std::invalid_argument("positive repetitions required");
            if (argc == 4)
                query.set("preparation", Json::string(argv[3]));
            const auto start = std::chrono::steady_clock::now();
            const auto first = session.execute(query).dump();
            const auto prepared = std::chrono::steady_clock::now();
            for (std::uint64_t i = 0; i < repeats; ++i)
                if (session.execute(query).dump() != first)
                    throw std::logic_error("reuse changed the complete certificate");
            const auto done = std::chrono::steady_clock::now();
            auto out = stats();
            out.set("first_call_microseconds",
                    Json::integer(
                        std::chrono::duration_cast<std::chrono::microseconds>(prepared - start)
                            .count()));
            out.set(
                "reused_calls_microseconds",
                Json::integer(std::chrono::duration_cast<std::chrono::microseconds>(done - prepared)
                                  .count()));
            out.set("repetitions", Json::string(std::to_string(repeats)));
            out.set("certificate_bytes", Json::string(std::to_string(first.size())));
            std::cout << out.dump() << '\n';
            return 0;
        }
        auto check = [](bool value, const char *message) {
            if (!value)
                throw std::logic_error(message);
        };
        session.validate(query);
        check(stats().at("retained_plans").integer() == 0, "admission retained native factors");
        auto first = session.execute(query);
        const auto bytes = first.dump();
        check(session.execute(query).dump() == bytes, "reuse changed novelty certificate");
        check(stats().at("preparations").integer() == 1, "same action refactored the native plan");
        auto changed = Json::parse(query.dump());
        changed.set("old", Json::parse(query.at("old").dump()));
        changed.at("old").set("field", Json::string("zero"));
        changed.set("causes", Json::parse(R"(["new","zero","zero"] )"));
        auto different = session.execute(changed);
        check(stats().at("preparations").integer() == 1, "input changes rebuilt shared factors");
        check(different.at("result").at("causes").at("word").size() == 3 &&
                  different.at("result").at("causes").at("fields").size() == 2,
              "repeated causes lost ordered axes or duplicated their value storage");
        changed = Json::parse(query.dump());
        changed.at("new").set("publication", Json::string("self"));
        check(session.execute(changed).at("result").at("split").dump() ==
                  first.at("result").at("split").dump(),
              "separately admitted byte-equal source changed its declared action");
        check(stats().at("preparations").integer() == 2 &&
                  stats().at("retained_plans").integer() == 2,
              "independent equal-key sources were pooled into one owner");
        const auto n = std::numeric_limits<std::size_t>::max();
        rfc::runtime::ExecutionHostPolicy policy;
        policy.name = "current-novelty-lease";
        policy.limits = {n, n, n, n, n};
        auto refuse = [&](const rfc::runtime::ExecutionLease &lease) {
            const auto before = stats().dump();
            bool refused = false;
            try {
                (void)session.execute_native_novelty(query, lease);
            } catch (const std::exception &) {
                refused = true;
            }
            check(refused && stats().dump() == before, "cached plan ignored current lease");
        };
        refuse(rfc::runtime::ExecutionLease(policy, std::make_shared<std::atomic_bool>(true)));
        for (auto *limit :
             {&policy.limits.records, &policy.limits.rational_slots, &policy.limits.sparse_entries,
              &policy.limits.scalar_work, &policy.limits.integer_bits}) {
            *limit = 0;
            refuse(rfc::runtime::ExecutionLease(policy));
            *limit = n;
        }
        changed = Json::parse(query.dump());
        changed.set("preparation", Json::string("ephemeral"));
        check(session.execute(changed).at("result").at("split").dump() ==
                  first.at("result").at("split").dump(),
              "ephemeral action changed mathematics");
        const auto before_clear = stats();
        session.clear_native_novelty_preparations();
        check(stats().at("retained_plans").integer() == 0 &&
                  stats().at("numeric_entries").integer() == 0,
              "clear retained factor entries");
        first.at("result").set("split", Json::object());
        check(session.execute(query).dump() == bytes,
              "detached output or cache clear changed source");
        check(stats().at("preparations").integer() == before_clear.at("preparations").integer() + 1,
              "cleared factors were not rebuilt");
        // A plan is an admitted request sequence, not a sixteen step math cap.
        auto plan = Json::object(), steps = Json::array();
        plan.set("schema", Json::string("gadget.scientific.plan.v1"));
        plan.set("id", Json::string("seventeen-complete-actions"));
        plan.set("question", Json::string("Repeat the same complete comparison"));
        plan.set("source_digest", Json::string(source->digest()));
        for (std::size_t i = 0; i < 17; ++i) {
            auto step = Json::object();
            step.set("id", Json::string(std::to_string(i)));
            step.set("query", query);
            steps.push(step);
        }
        plan.set("steps", steps);
        check(session.execute_plan(plan).at("steps").size() == 17,
              "arbitrary plan ceiling remains");
        std::cout << stats().dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
