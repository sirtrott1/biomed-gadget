#include "genomic_detail.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

namespace {
std::size_t checks = 0;
void check(bool condition, const char *message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
}
int main() {
    try {
        using namespace gadget;
        using namespace gadget::assay_detail;
        {
            const auto request = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"sequence_motifs","resource":"sequence","motifs":[],"strands":"forward"})");
            for (std::size_t i = 0; i < 129; ++i) {
                auto pattern = Json::object();
                pattern.set("id", Json::string(i == 0 ? std::string(5001, 'x') : std::to_string(i)));
                pattern.set("pattern", Json::string(i == 0 ? "A" : "C"));
                request.at("motifs").push(pattern);
            }
            const auto declaration = Json::parse(R"({"format":"fasta"})");
            const auto result = motif_analysis(">source\n" + std::string(10003, 'A') + "N\n", declaration, request);
            check(result.at("summaries").size() == 129, "every declared pattern retained");
            check(result.at("hits").size() == 10132, "every matching window retained");
            for (std::size_t i = 0; i < 10004; ++i) {
                const auto hit = result.at("hits").at(i);
                check(hit.at("start").integer() == static_cast<std::int64_t>(i) &&
                          hit.at("end").integer() == static_cast<std::int64_t>(i + 1),
                      "complete source window order");
                check(hit.at("definite").boolean() == (i < 10003), "IUPAC definite versus possible");
                check(hit.at("motif").text() == std::string(5001, 'x'), "full pattern identity");
            }
            for (std::size_t i = 1; i < 129; ++i) {
                const auto summary = result.at("summaries").at(i);
                check(summary.at("definite_occurrences").integer() == 0 &&
                          summary.at("possible_occurrences").integer() == 1,
                      "all marginal ambiguity counts");
            }
        }
        {
            auto request = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"sequence_motifs","resource":"sequence","motifs":[{"id":"long","pattern":"A"}],"strands":"forward"})");
            request.at("motifs").at(0).set("pattern", Json::string(std::string(513, 'A')));
            const auto result = motif_analysis(">source\n" + std::string(200003, 'A') + "\n",
                                               Json::parse(R"({"format":"fasta"})"), request);
            check(result.at("hits").size() == 199491, "complete long-pattern scan above 100 million comparisons");
            for (std::size_t i = 0; i < result.at("hits").size(); ++i) {
                const auto hit = result.at("hits").at(i);
                check(hit.at("start").integer() == static_cast<std::int64_t>(i) &&
                          hit.at("end").integer() == static_cast<std::int64_t>(i + 513) &&
                          hit.at("definite").boolean(), "every long-pattern source window");
            }
        }
        {
            std::string sequence, qualities, bytes;
            for (std::size_t i = 0; i < 100001; ++i) {
                sequence += "CSWN"[i % 4];
                qualities += "!5I~"[i % 4];
            }
            for (std::size_t i = 0; i < 101; ++i)
                bytes += "@read" + std::to_string(i) + "\n" + sequence + "\n+\n" + qualities + "\n";
            const auto result = quality_analysis(
                {Json::parse(R"({"format":"fastq"})"), bytes, {}},
                Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_quality","resource":"reads","quality_thresholds":[0,20,40,93]})"));
            check(result.at("total_bases").integer() == 10100101, "all bases beyond former ten-million ceiling");
            check(result.at("reads").size() == 101 && result.at("cycle_histograms").size() == 100001,
                  "all source reads and long cycles retained");
            const std::array<std::int64_t, 4> scores{0,20,40,93}, below{0,25001,50001,75001};
            for (std::size_t i = 0; i < 101; ++i) {
                const auto read = result.at("reads").at(i);
                check(read.at("id").text() == "read" + std::to_string(i), "each original read identity");
                check(read.at("gc_lower_bases").text() == "50001" &&
                          read.at("gc_upper_bases").text() == "75001" &&
                          read.at("quality_score_sum").text() == "3825000", "independent complete GC and stored-score sums");
                for (std::size_t j = 0; j < below.size(); ++j)
                    check(read.at("quality_thresholds").at(j).at("bases_below").text() == std::to_string(below[j]),
                          "independent stored-score threshold counts");
            }
            for (std::size_t i = 0; i < 100001; ++i) {
                const auto histogram = result.at("cycle_histograms").at(i);
                check(histogram.at("cycle").integer() == static_cast<std::int64_t>(i + 1), "complete physical cycle axis");
                check(histogram.at("phred33_counts").size() == 94, "intrinsic printable Phred+33 alphabet");
                for (std::size_t q = 0; q < 94; ++q)
                    check(histogram.at("phred33_counts").at(q).integer() ==
                              (static_cast<std::int64_t>(q) == scores[i % 4] ? 101 : 0),
                          "every complete cycle histogram coordinate");
            }
        }
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
