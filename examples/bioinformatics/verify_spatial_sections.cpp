#include <fstream>
#include <gadget/assay.hpp>
#include <iostream>
#include <map>
#include <rfc/storage/rfcx.hpp>

using namespace gadget;
using Q = rfc::numeric::Rational;
namespace {
void check(bool value, const char *message) { if (!value) throw std::runtime_error(message); }
Q number(const NativeRecord &record) {
    if (record.kind() == NativeRecordKind::Integer) return Q(record.integer());
    if (record.kind() == NativeRecordKind::Rational || record.kind() == NativeRecordKind::ExactNumber ||
        record.kind() == NativeRecordKind::Decimal) return record.rational();
    return Q(record.kind() == NativeRecordKind::Object ? record.at("rational").text() : record.text());
}
struct Summary {
    std::string context, field;
    std::uint64_t anchor{};
    std::size_t cells{};
    Q input{}, mean{}, quadrance{}, residual{}, hodge{}, green{};
    bool has_hodge{}, has_green{};
};
} // namespace
int main(int argc, char **argv) try {
    check(argc == 5, "usage: verify_spatial_sections RESULT.rfcxs WORKSPACE.rfcxs ASSAY SUMMARY.tsv");
    const auto publication = NativeRecord::read(argv[1]);
    const auto result = publication.at("result");
    check(result.at("schema").text() == "gadget.scientific.sequence.result.v1", "sequence result schema");
    AssayWorkspace workspace(argv[2]);
    const auto study = workspace.study(argv[3]);
    check(result.at("source_digest").text() == study->digest(), "sequence source witness differs");
    const auto original = study->endpoint();
    check(rfc::storage::encode_rfcx(original->source).bytes ==
          rfc::storage::encode_rfcx(publication.endpoint()->source).bytes, "sequence changed its source topology");
    const rfc::storage::StableSourceIndex index(original->source);
    const auto producer = study->producer_publication().at("result").at("result");
    check(producer.at("schema").text() == "gadget.spatial.fields.v1", "spatial producer required");
    std::map<std::string, std::pair<std::uint64_t, std::uint64_t>> charts;
    for (std::size_t i = 0; i < producer.at("charts").size(); ++i) {
        const auto chart = producer.at("charts").at(i);
        check(charts.emplace(chart.at("id").text(), std::pair{std::stoull(chart.at("block").text()),
                                                            std::stoull(chart.at("anchor").text())}).second,
              "repeated native spatial chart");
    }
    std::map<std::pair<std::string, std::string>, Summary> summaries;
    std::size_t checked = 0, incompatible = 0, green = 0, hodge = 0;
    const auto answers = result.at("answers");
    for (std::size_t a = 0; a < answers.size(); ++a) {
        const auto certificate = answers.at(a).at("certificate"), query = certificate.at("query");
        check(certificate.at("source_digest").text() == study->digest(), "query source witness differs");
        const auto context = query.at("context").text(), operation = query.at("operation").text();
        check(charts.contains(context), "query chart is not source declared");
        const auto [block_key, anchor] = charts.at(context);
        const auto block = std::find_if(original->local_blocks.begin(), original->local_blocks.end(),
                                       [&](const auto &b) { return b.block_key == block_key; });
        check(block != original->local_blocks.end(), "source chart missing");
        std::vector<std::uint32_t> cells, edges;
        for (const auto &cell : block->cells) {
            check(cell.grade.value < 2, "spatial oracle requires grade-one charts");
            (cell.grade.value == 0 ? cells : edges).push_back(index.cell(cell).value().index);
        }
        std::sort(cells.begin(), cells.end()); std::sort(edges.begin(), edges.end());
        std::map<std::uint32_t, std::size_t> local;
        for (std::size_t i = 0; i < cells.size(); ++i) local.emplace(cells[i], i);
        const auto &arena = original->source.grade({1});
        std::vector<std::pair<std::size_t, std::size_t>> pairs;
        for (auto edge : edges) {
            const auto [begin, end] = arena.incidence_range(edge);
            check(end == begin + 2 && original->source.incidence_coefficient({1}, begin) == -1 &&
                  original->source.incidence_coefficient({1}, begin + 1) == 1, "oracle requires oriented physical pairs");
            pairs.emplace_back(local.at(arena.lower_index[begin]), local.at(arena.lower_index[begin + 1]));
        }
        const auto anchor_index = local.at(index.cell({{0}, anchor}).value().index);
        const auto values = [&](const std::string &field) {
            const auto pin = study->field_accession(field);
            check(pin.carrier()->grade.value == 0, "section source field grade differs");
            const auto tensor = pin.tensor();
            std::vector<Q> out;
            for (auto cell : cells) {
                const auto value = tensor->state(cell); check(value.has_value(), "unobserved source field coordinate");
                out.push_back(*value);
            }
            return out;
        };
        const auto laplace = [&](const std::vector<Q> &input) {
            std::vector<Q> out(input.size());
            for (const auto &[i, j] : pairs) { const Q delta = input[j] - input[i]; out[i] -= delta; out[j] += delta; }
            return out;
        };
        if (operation == "local_sections") {
            check(query.at("construction").text() == "cellular_sheaf" && query.at("metric").text() == "identity",
                  "section declaration differs");
            const auto maps = query.at("maps");
            check(maps.size() == block->occurrences.size(), "all original flags need restriction maps");
            for (std::size_t i = 0; i < maps.size(); ++i)
                check(maps.at(i).at("matrix").size() == 1 && maps.at(i).at("matrix").at(0).size() == 1 &&
                      number(maps.at(i).at("matrix").at(0).at(0)) == 1, "restriction is not scalar identity");
            const auto output = certificate.at("result"), sections = output.at("sections");
            check(output.at("cohomology").at(0).at("dimension").integer() == 1,
                  "physical chart constant-section kernel differs");
            for (std::size_t s = 0; s < sections.size(); ++s) {
                const auto section = sections.at(s);
                const auto id = section.at("fields").at(0).text();
                const auto input = values(id); Q sum{0}; for (const auto &q : input) sum += q;
                const Q mean = sum / Q(std::to_string(input.size())); Q quadrance{0};
                check(section.at("values").size() == cells.size(), "complete restricted input differs");
                for (std::size_t i = 0; i < input.size(); ++i) {
                    const Q residual = input[i] - mean;
                    check(number(section.at("values").at(i)) == input[i] && number(section.at("glued").at(i)) == mean &&
                          number(section.at("residual").at(i)) == residual, "exact local projection differs");
                    quadrance += residual * residual;
                }
                const auto defects = section.at("coboundary_defect");
                check(defects.size() == pairs.size(), "complete sheaf defect carrier differs");
                for (std::size_t i = 0; i < pairs.size(); ++i)
                    check(number(defects.at(i)) == input[pairs[i].second] - input[pairs[i].first],
                          "original physical sheaf defect differs");
                check(number(section.at("residual_quadrance")) == quadrance &&
                      section.at("compatible").boolean() == (quadrance == 0), "compatibility or quadrance differs");
                auto &summary = summaries[{context, id}];
                summary.context = context; summary.field = id; summary.anchor = anchor; summary.cells = cells.size();
                summary.input = input[anchor_index]; summary.mean = mean; summary.quadrance = quadrance;
                summary.residual = input[anchor_index] - mean;
                ++checked; incompatible += quadrance != 0;
            }
        } else {
            check(operation == "green_action" || operation == "hodge", "unsupported spatial oracle query");
            const auto id = query.at("field").text(); const auto input = values(id);
            const auto output = certificate.at("result");
            check(output.size() == input.size(), "complete operator carrier differs");
            std::vector<Q> actual;
            for (std::size_t i = 0; i < output.size(); ++i) {
                check(output.at(i).at("native_key").text() == std::to_string(original->source.grade({0}).identities[cells[i]]),
                      "operator output source identity differs");
                actual.push_back(number(output.at(i).at("value")));
            }
            auto &summary = summaries[{context, id}];
            if (operation == "hodge") {
                check(actual == laplace(input), "source-bound Hodge action differs");
                summary.hodge = actual[anchor_index]; summary.has_hodge = true; ++hodge;
            } else {
                const Q lambda(query.at("lambda").text()); check(lambda > 0, "Green lambda must be positive");
                const auto action = laplace(actual); Q before{0}, after{0};
                for (std::size_t i = 0; i < input.size(); ++i) {
                    check(actual[i] + lambda * action[i] == input[i], "independent Green identity differs");
                    before += input[i]; after += actual[i];
                }
                check(before == after, "Green constant section conservation differs");
                summary.green = actual[anchor_index]; summary.has_green = true; ++green;
            }
        }
    }
    std::ofstream table(argv[4]); check(bool(table), "cannot write presentation table");
    table << "context\tfield\tanchor\tlocal_cells\tcount\tlocal_mean\tresidual_quadrance\tanchor_residual\tanchor_hodge\tanchor_green\n";
    for (const auto &[key, summary] : summaries) {
        (void)key;
        check(summary.cells > 0 && summary.has_green && summary.has_hodge, "incomplete local method family");
        table << summary.context << '\t' << summary.field << '\t' << summary.anchor << '\t' << summary.cells << '\t'
              << summary.input.get_str() << '\t' << summary.mean.get_str() << '\t' << summary.quadrance.get_str() << '\t'
              << summary.residual.get_str() << '\t' << summary.hodge.get_str() << '\t' << summary.green.get_str() << '\n';
    }
    check(bool(table), "presentation table write failed");
    std::cout << "PASS\nsource_digest\t" << study->digest() << "\nsections\t" << checked << "\nincompatible\t" << incompatible
              << "\nGreen_identities\t" << green << "\nHodge_identities\t" << hodge << '\n';
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
