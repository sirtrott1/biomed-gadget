#include <charconv>
#include <chrono>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        using Clock = std::chrono::steady_clock;
        if (argc != 4)
            throw std::invalid_argument("usage: gadget_temporal_bench STUDY QUERY REPEATS");
        const std::string count(argv[3]);
        std::size_t repeats{};
        const auto [end, error] = std::from_chars(count.data(), count.data() + count.size(), repeats);
        if (error != std::errc{} || end != count.data() + count.size() || !repeats)
            throw std::invalid_argument("REPEATS must be a positive addressable integer");
        const auto start = Clock::now();
        auto study = std::make_shared<const ScientificStudy>(argv[1]);
        ScientificSession session(study);
        const auto query = Json::parse(read_text(argv[2]));
        const auto admitted = Clock::now();
        const auto first = session.execute(query);
        const auto first_done = Clock::now();
        const auto counts = session.temporal_preparation_statistics();
        const auto expected = first.dump_compact();
        auto timings = Json::array();
        for (std::size_t i = 0; i < repeats; ++i) {
            const auto before = Clock::now();
            const auto actual = session.execute(query);
            const auto done = Clock::now();
            if (actual.dump_compact() != expected ||
                session.temporal_preparation_statistics().dump_compact() != counts.dump_compact())
                throw std::runtime_error("reused temporal execution changed certificate or preparations");
            timings.push(Json::number(std::chrono::duration<double>(done - before).count()));
        }
        auto report = Json::object();
        report.set("status", Json::string("PASS"));
        report.set("admission_seconds", Json::number(std::chrono::duration<double>(admitted - start).count()));
        report.set("first_query_seconds", Json::number(std::chrono::duration<double>(first_done - admitted).count()));
        report.set("repeated_query_seconds", timings);
        report.set("preparations", counts);
        report.set("certificate_sha256", first.at("certificate_sha256"));
        report.set("gadget_sources_sha256", first.at("gadget_sources_sha256"));
        report.set("rexgraph_headers_sha256", first.at("rexgraph_headers_sha256"));
        std::cout << report.dump() << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
