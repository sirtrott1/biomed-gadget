#include "sequence_detail.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>

namespace {
std::size_t checks = 0;
void check(bool value, const char *message) {
    ++checks;
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&f) {
    bool failed = false;
    try { f(); } catch (const std::exception &) { failed = true; }
    check(failed, "invalid representation/source accepted");
}
}
int main() {
    try {
        using namespace gadget::assay_detail;
        using namespace gadget::assay_detail::sequence;
        using gadget::Rational;
        const auto start = std::chrono::steady_clock::now();
        check(token_bases == 32, "intrinsic uint64 seed domain");
        check(token(std::string(32, 'A')) == 0 &&
                  token(std::string(32, 'T')) == std::numeric_limits<std::uint64_t>::max() &&
                  token(std::string(32, 'C')) == 0x5555555555555555ULL,
              "injective full-word two-bit encoding");
        rejects([] { (void)token(std::string(33, 'A')); });
        rejects([] { (void)token("AN"); });
        SequenceRecord reference;
        reference.id = "large";
        reference.bases = std::string(1000032, 'A') + 'N' + std::string(32, 'T');
        SeedIndex index;
        index.request(std::string(32, 'A'));
        index.request(std::string(32, 'T'));
        index.build({reference});
        const auto &hits = index.find(std::string(32, 'A'));
        check(hits.size() == 1000001, "complete index beyond former million-hit cap");
        for (std::size_t i = 0; i < hits.size(); ++i)
            check(hits[i].reference == 0 && hits[i].offset == i, "every seed source offset retained");
        const auto &last = index.find(std::string(32, 'T'));
        check(last.size() == 1 && last[0].offset == 1000033, "word mask and ambiguity reset");
        const auto receipt = index.receipt();
        check(receipt.at("retained_seed_hits").integer() == 1000002 &&
                  receipt.at("reference_base_visits").integer() == 1000065,
              "complete seed work counters");
        Routes routes(3, "qualification.parallel.steps");
        for (std::size_t i = 0; i < 100001; ++i)
            routes.add(i % 2, i % 2 + 1, Rational{static_cast<unsigned long>(i % 3)});
        rejects([&] { routes.add(3, 0, 1); });
        rejects([&] { routes.add(0, 1, -1); });
        const auto source = routes.freeze();
        const auto &grade = source.grade({1});
        check(grade.cell_count() == 100001 && grade.incidence_count() == 200002,
              "complete native carrier beyond former edge cap");
        for (std::size_t i = 0; i < grade.cell_count(); ++i) {
            check(grade.relation_identities.at(i).value == 1000000 + i,
                  "source relation identity order");
            check(source.occurrence_identity({1}, 2 * i).value == 2 * i + 1 &&
                      source.occurrence_identity({1}, 2 * i + 1).value == 2 * i + 2,
                  "retained occurrence identity order");
            check(source.incidence_coefficient({1}, 2 * i) == -1 &&
                      source.incidence_coefficient({1}, 2 * i + 1) == 1,
                  "native oriented boundary coefficients");
        }
        const auto graph = rfc::view::project_pairwise_directed(source, routes.costs());
        const auto route = rfc::path::dijkstra_exact(graph, graph.object_index(1), graph.object_index(3));
        check(route.reachable && route.distance == 0 && route.relations.size() == 2,
              "native shortest route across complete parallel source");
        const auto realization = realized_receipt(source, routes, route.relations, 0, 2, route.distance);
        check(realization.at("boundary_verified").boolean() &&
                  realization.at("cost_verified").boolean(), "complete native route realization");
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout << "{\"status\":\"PASS\",\"independent_checks\":" << checks
                  << ",\"runtime_seconds\":" << seconds << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
