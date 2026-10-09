#include <chrono>
#include <gadget/execution.hpp>
#include <gadget/store.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("usage: query_bench STORE QUERY_JSON REPORT_JSON");
        using namespace gadget;
        using Clock = std::chrono::steady_clock;
        const auto start = Clock::now();
        ConcurrentStudyStore store(argv[1]);
        const auto recovered = Clock::now();
        const auto study = store.snapshot().workspace();
        QuerySession session(study);
        const auto query = Json::parse(read_text(argv[2]));
        const auto first_started = Clock::now();
        const auto first = session.execute(query);
        const auto first_completed = Clock::now();
        auto warm_seconds = Json::array();
        for (int i = 0; i < 5; ++i) {
            const auto before = Clock::now();
            const auto value = session.execute(query);
            if (value.dump() != first.dump())
                throw std::runtime_error("warm preparation changed result/certificate");
            warm_seconds.push(
                Json::number(std::chrono::duration<double>(Clock::now() - before).count()));
        }
        const auto fresh_started = Clock::now();
        const auto fresh = execute_query(*study, query);
        const auto fresh_completed = Clock::now();
        if (fresh.dump() != first.dump())
            throw std::runtime_error("fresh execution disagrees with prepared session");
        auto report = Json::object();
        report.set("status", Json::string("PASS"));
        report.set("root", study->store_origin.at("root"));
        report.set("source_digest", Json::string(study->source_digest));
        report.set("query", query);
        report.set("certificate_sha256", first.at("certificate_sha256"));
        report.set("cold_recovery_seconds",
                   Json::number(std::chrono::duration<double>(recovered - start).count()));
        report.set(
            "first_query_seconds",
            Json::number(std::chrono::duration<double>(first_completed - first_started).count()));
        report.set("warm_query_seconds", warm_seconds);
        report.set(
            "fresh_query_seconds",
            Json::number(std::chrono::duration<double>(fresh_completed - fresh_started).count()));
        report.set("prepared_contexts", Json::integer(session.preparation_count()));
        report.set("rexgraph_headers_sha256", Json::string(GADGET_REXGRAPH_SHA256));
        report.set("gadget_sources_sha256", Json::string(GADGET_SOURCE_SHA256));
        write_text(argv[3], report.dump() + "\n");
        std::cout << report.dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
