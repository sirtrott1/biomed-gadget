#include <atomic>
#include <chrono>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2 && argc != 3)
            throw std::invalid_argument("expected observation case");
        const fs::path root = argv[1];
        auto owner = std::make_shared<const ScientificStudy>(root / "study/study.json");
        ScientificSession session(owner);
        const auto query = Json::parse(read_text(root / "query.json"));
        auto stats = [&] { return session.native_observation_preparation_statistics(); };
        if (argc == 3) {
            const auto repeats = std::stoull(argv[2]);
            if (!repeats)
                throw std::invalid_argument("positive benchmark repetitions required");
            const auto start = std::chrono::steady_clock::now();
            const auto reference = session.execute(query).dump();
            const auto prepared = std::chrono::steady_clock::now();
            for (std::uint64_t i = 0; i < repeats; ++i)
                if (session.execute(query).dump() != reference)
                    throw std::logic_error("reused observation changed full certificate");
            const auto done = std::chrono::steady_clock::now();
            auto receipt = stats();
            receipt.set("repetitions", Json::string(std::to_string(repeats)));
            receipt.set("first_call_microseconds",
                        Json::integer(
                            std::chrono::duration_cast<std::chrono::microseconds>(prepared - start)
                                .count()));
            receipt.set(
                "reused_calls_microseconds",
                Json::integer(std::chrono::duration_cast<std::chrono::microseconds>(done - prepared)
                                  .count()));
            std::cout << receipt.dump() << '\n';
            return 0;
        }
        auto check = [](bool value, const char *message) {
            if (!value)
                throw std::logic_error(message);
        };
        session.validate(query);
        check(stats().at("base_solve_preparations").integer() == 0, "admission prepared a solve");
        auto first = session.execute(query);
        check(session.execute(query).dump() == first.dump(), "reuse changed certificate");
        auto changed = Json::parse(query.dump());
        changed.set("workers", Json::integer(1));
        check(session.execute(changed).at("result").at("fields").dump() ==
                  first.at("result").at("fields").dump(),
              "CPU workers changed observations");
        changed.set("fields", Json::parse(R"(["x2","x","x2"])"));
        const auto ordered = session.execute(changed).at("result").at("fields");
        check(ordered.at(0).dump() == ordered.at(2).dump() &&
                  ordered.at(1).at("values").dump() ==
                      first.at("result").at("fields").at(0).at("values").dump(),
              "ordered repeated field slots changed meaning");
        changed = Json::parse(query.dump());
        changed.set("accession", Json::string("full"));
        changed.set("observation", Json::parse(R"({"kind":"linear"})"));
        check(session.execute(changed)
                  .at("result")
                  .at("native_preparation")
                  .at("uses_full_solve")
                  .boolean(),
              "complete reading did not use primal route");
        check(stats().at("base_solve_preparations").integer() == 1 &&
                  stats().at("observation_preparations").integer() == 2,
              "accession change refactored the source operator");
        changed.set("ordering",
                    Json::parse(R"({"kind":"declared","rows":[2,0,3,1],"columns":[1,3,0,2]})"));
        const auto permuted = session.execute(changed);
        check(stats().at("base_solve_preparations").integer() == 2,
              "different declared execution ordering reused wrong owner");
        auto zero_query = Json::parse(query.dump());
        zero_query.set("accession", Json::string("zero_map"));
        zero_query.set("observation", Json::parse(R"({"kind":"linear"})"));
        const auto zero_result = session.execute(zero_query).at("result").at("fields");
        for (std::size_t i = 0; i < zero_result.size(); ++i)
            for (std::size_t j = 0; j < zero_result.at(i).at("values").size(); ++j)
                check(zero_result.at(i).at("values").at(j).text() == "0", "zero row changed value");
        const auto n = std::numeric_limits<std::size_t>::max();
        rfc::runtime::ExecutionHostPolicy policy;
        policy.name = "current-observation-lease";
        policy.limits = {n, n, n, n, n};
        auto refuse = [&](const rfc::runtime::ExecutionLease &lease) {
            const auto before = stats().dump();
            for (const auto &requested : {query, zero_query}) {
                bool refused = false;
                try {
                    (void)session.execute_native_observation(requested, lease);
                } catch (const std::exception &) {
                    refused = true;
                }
                check(refused && stats().dump() == before,
                      "cached observation ignored current lease or changed counters");
            }
        };
        refuse(rfc::runtime::ExecutionLease(policy, std::make_shared<std::atomic_bool>(true)));
        policy.limits.records = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.records = n;
        policy.limits.rational_slots = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.rational_slots = n;
        policy.limits.sparse_entries = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.sparse_entries = n;
        policy.limits.scalar_work = 0;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.scalar_work = n;
        policy.limits.integer_bits = 1;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.limits.integer_bits = n;
        policy.backend = rfc::runtime::ExecutionBackend::ExactGpu;
        refuse(rfc::runtime::ExecutionLease(policy));
        policy.fallback = rfc::runtime::BackendFallback::AllowExactCpuFallback;
        check(session.execute_native_observation(query, rfc::runtime::ExecutionLease(policy))
                  .at("result")
                  .at("fields")
                  .at(0)
                  .at("native_execution")
                  .at("fallback")
                  .boolean(),
              "explicit CPU fallback omitted execution receipt");
        first.at("result")
            .at("publication")
            .at("declaration")
            .set("id", Json::string("changed returned metadata"));
        check(session.execute(query)
                      .at("result")
                      .at("publication")
                      .at("declaration")
                      .at("id")
                      .text() != "changed returned metadata",
              "returned JSON aliases immutable metadata");
        const auto bytes = read_text(root / "study/study.rfcxs");
        const auto retained = session.execute(query);
        write_text(root / "study/study.rfcxs", "changed after admission");
        const auto after = session.execute(query);
        write_text(root / "study/study.rfcxs", bytes);
        check(after.dump() == retained.dump(), "plan borrowed mutable source files");
        owner.reset();
        std::cout << stats().dump() << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
