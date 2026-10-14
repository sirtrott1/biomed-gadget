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
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, "invalid representation/source accepted");
}
} // namespace
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
        {
            const std::vector<std::string> expected{"ACGTACGTNAAAA", std::string(70, 'T'), "AC",
                                                    "GT"};
            ResourceImage image{gadget::Json::parse(R"({"format":"fasta"})"),
                                ">first original header\r\nac\r\ngtac\r\ngtnAA\r\naa\r\n>repeat\n" +
                                    std::string(35, 't') + "\n" + std::string(35, 't') +
                                    "\n>left\nac\n>right\ngt",
                                {}};
            const auto original = image.bytes;
            const SequenceSource source(image);
            check(source.size() == expected.size(), "every original sequence record");
            for (std::size_t r = 0; r < source.size(); ++r) {
                check(source.record(r).bases == expected[r], "original wrapped/lowercase bases");
                for (std::size_t begin = 0; begin <= expected[r].size(); ++begin)
                    for (std::size_t length = 0; length <= expected[r].size() - begin; ++length)
                        check(source.bases(r, begin, length) == expected[r].substr(begin, length),
                              "every original coordinate window across source lines");
            }
            const std::vector<std::string> seeds{"A",    "AC",   "ACGT",
                                                 "AAAA", "TTTT", std::string(32, 'T')};
            for (const auto width : {1u, 2u, 3u, 31u, 32u, 33u, 100u}) {
                SeedIndex streamed;
                for (const auto &seed : seeds)
                    streamed.request(seed);
                streamed.build(source, width);
                for (const auto &seed : seeds) {
                    std::vector<std::pair<std::size_t, std::size_t>> oracle, actual;
                    for (std::size_t r = 0; r < expected.size(); ++r)
                        for (std::size_t offset = 0; offset + seed.size() <= expected[r].size();
                             ++offset)
                            if (expected[r].substr(offset, seed.size()) == seed)
                                oracle.emplace_back(r, offset);
                    for (const auto &hit : streamed.find(seed))
                        actual.emplace_back(hit.reference, hit.offset);
                    check(actual == oracle,
                          "independent complete seed positions across every seam");
                }
                check(streamed.receipt().at("peak_normalized_reference_block_bases").integer() <=
                          width,
                      "normalized reference blocks respect caller width");
            }
            rejects([&] { (void)source.bases(0, expected[0].size(), 1); });
            rejects([&] { source.blocks(0, 0, [](auto, auto) {}); });
            check(image.bytes == original, "all original bytes unchanged");
            ResourceImage fastq{gadget::Json::parse(R"({"format":"fastq"})"),
                                "@first desc\r\nacgt\r\n+first desc\r\nACGT\r\n@second\nTT\n+\n!!",
                                {}};
            const SequenceSource reads(fastq);
            check(reads.size() == 2 && reads.record(0).bases == "ACGT" &&
                      reads.record(1).bases == "TT",
                  "FASTQ bases exclude header and qualities");
            SeedIndex quality_control;
            quality_control.request("ACGT");
            quality_control.build(reads, 1);
            check(quality_control.find("ACGT").size() == 1,
                  "FASTQ quality tokens never create seed occurrences");
            const std::string repeated_bytes =
                "@same first lane\r\nac\r\n+same first lane\r\n!!\r\n"
                "@same second lane\nAC\n+\n55\n@same\nCG\n+\n~~";
            ResourceImage repeated{gadget::Json::parse(R"({"format":"fastq"})"), repeated_bytes, {}};
            const SequenceSource occurrences(repeated);
            check(occurrences.size() == 3, "repeated FASTQ labels discarded original records");
            for (std::size_t i = 0; i < 3; ++i) {
                check(occurrences.id(i) == "same" && occurrences.length(i) == 2 &&
                          occurrences.evidence(i).at("record_index").integer() == static_cast<std::int64_t>(i),
                      "original FASTQ label and occurrence index differ");
                check(occurrences.record(i).bases == (i < 2 ? "AC" : "CG"),
                      "repeated label replaced distinct original bases");
                const auto &span = occurrences.original_span(i);
                check(repeated_bytes.substr(span.start, span.end - span.start).starts_with("@same"),
                      "repeated FASTQ occurrence lost its original byte span");
                for (const auto width : {1u, 2u, 7u}) {
                    std::string joined;
                    occurrences.blocks(i, width, [&](auto offset, auto block) {
                        check(offset == joined.size(), "repeated read block coordinate differs");
                        joined += block;
                    });
                    check(joined == occurrences.record(i).bases,
                          "repeated read block partition differs");
                }
            }
            SeedIndex repeated_index;
            repeated_index.request("AC"); repeated_index.build(occurrences, 1);
            check(repeated_index.find("AC").size() == 2 &&
                      repeated_index.find("AC")[0].reference == 0 &&
                      repeated_index.find("AC")[1].reference == 1,
                  "seed index merged repeated source labels");
            check(repeated.bytes == repeated_bytes, "repeated original FASTQ bytes changed");
            ResourceImage duplicate_fasta{gadget::Json::parse(R"({"format":"fasta"})"),
                                         ">same\nAC\n>same\nCG\n", {}};
            rejects([&] { (void)SequenceSource(duplicate_fasta); });
            ResourceImage empty_identifier{gadget::Json::parse(R"({"format":"fastq"})"),
                                          "@\nAC\n+\n!!\n", {}};
            rejects([&] { (void)SequenceSource(empty_identifier); });
        }
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
            check(hits[i].reference == 0 && hits[i].offset == i,
                  "every seed source offset retained");
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
        const auto route =
            rfc::path::dijkstra_exact(graph, graph.object_index(1), graph.object_index(3));
        check(route.reachable && route.distance == 0 && route.relations.size() == 2,
              "native shortest route across complete parallel source");
        auto retained = route_source(source, routes);
        SourceGradeFields fields(retained.endpoint);
        member(std::get<N::Object>(retained.declaration.data), "witness",
               retain_route(fields, source, routes, route.relations, 0, 2, route.distance));
        const auto realization = retained_route_view(source_record(std::move(retained)).at("witness"));
        check(realization.at("boundary_verified").boolean() &&
                  realization.at("cost_verified").boolean(),
              "complete native route realization");
        const auto seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cout << "{\"status\":\"PASS\",\"independent_checks\":" << checks
                  << ",\"runtime_seconds\":" << seconds << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
