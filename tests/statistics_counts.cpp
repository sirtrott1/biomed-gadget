#include "../src/statistics/count_answers.hpp"
#include <cmath>
#include <gadget/assay.hpp>
#include <gadget/io.hpp>
#include <gadget/record_view.hpp>
#include <gadget/statistics/count_analysis.hpp>
#include <gadget/statistics/counts.hpp>
#include <iostream>
#include <rfc/runtime/execution_policy.hpp>
#include <sstream>

using namespace gadget;
using namespace gadget::statistics;
namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void near(double a, double b, double tolerance = 1e-8) {
    if (!std::isfinite(a) || std::abs(a - b) > tolerance * (1 + std::abs(b)))
        throw std::runtime_error("numerical reference differs: " + std::to_string(a) + " / " +
                                 std::to_string(b));
}
Json real(double value) { return std::isfinite(value) ? Json::number(value) : Json{}; }
std::vector<double> values(const Json &array) {
    std::vector<double> out;
    for (std::size_t i = 0; i < array.size(); ++i)
        out.push_back(std::stod(array.at(i).text()));
    return out;
}
std::vector<rfc::numeric::RationalVector> rows(const Json &array) {
    std::vector<rfc::numeric::RationalVector> out;
    for (std::size_t i = 0; i < array.size(); ++i) {
        rfc::numeric::RationalVector row;
        for (std::size_t j = 0; j < array.at(i).size(); ++j)
            row.emplace_back(array.at(i).at(j).text());
        out.push_back(std::move(row));
    }
    return out;
}
Json list(std::span<const double> values) {
    auto out = Json::array();
    for (const auto value : values)
        out.push(real(value));
    return out;
}
Json evaluate(const Json &q) {
    const auto op = q.at("operation").text();
    auto out = Json::object();
    CountWork work(1000000000);
    if (op == "mass")
        out.set("value",
                real(count_log_mass(std::stod(q.at("y").text()), std::stod(q.at("mu").text()),
                                    std::stod(q.at("alpha").text()))));
    else if (op == "tails") {
        out.set("normal", real(normal_log_two_sided(std::stod(q.at("z").text()))));
        out.set("chi_square", real(chi_square_log_survival(std::stod(q.at("statistic").text()),
                                                           q.at("degrees").integer())));
    } else if (op == "support") {
        const PreparedCountDesign design(rows(q.at("design")));
        out.set("status", Json::string(design.support_status(values(q.at("y")), work)));
    } else if (op == "adjust")
        out.set("values",
                list(adjust_log_probabilities(values(q.at("values")), q.at("policy").text())));
    else {
        PreparedCountDesign design(rows(q.at("design")));
        const auto y = values(q.at("y")), offset = values(q.at("offset")),
                   contrast = values(q.at("contrast"));
        double alpha = std::stod(q.at("alpha").text());
        if (op == "profile") {
            DispersionControl control;
            control.minimum = 1e-5;
            control.maximum = 10;
            control.grid_points = 20;
            control.log_tolerance = 1e-6;
            auto fit = profile_count_dispersion(
                design, 1, [&](std::size_t) { return y; }, offset, control, {}, work);
            alpha = fit.dispersion;
            out.set("profile_converged", Json::boolean(fit.converged));
            out.set("profile_objective", real(fit.objective));
        }
        const auto fit = design.fit(y, offset, alpha, contrast, {}, work);
        out.set("converged", Json::boolean(fit.converged));
        out.set("status", Json::string(fit.status));
        out.set("alpha", real(alpha));
        out.set("coefficients", list(fit.coefficients));
        out.set("means", list(fit.means));
        out.set("score", list(fit.score));
        out.set("log_likelihood", real(fit.log_likelihood));
        out.set("effect", fit.effect ? real(*fit.effect) : Json{});
        out.set("variance", fit.variance ? real(*fit.variance) : Json{});
        out.set("work", Json::integer(work.used()));
    }
    return out;
}
const rfc::storage::TemporalTensor &tensor(const NativeRecord &owner, const NativeRecord &view) {
    const auto key = view.at("tensor").integer();
    for (const auto &field : owner.endpoint()->tensors)
        if (key == rfc::numeric::Integer(std::to_string(field.key)))
            return field;
    throw std::runtime_error("native field absent");
}
void native(const fs::path &workspace, const fs::path &request, const fs::path &destination) {
    AssayWorkspace assays(workspace);
    const auto study = assays.study("measurements");
    auto q = Json::parse(read_text(request));
    q.set("source_digest", Json::string(study->digest()));
    const auto admitted = admit_record_view(q);
    const auto before = study->native_publication().encode();
    auto prepared = PreparedCountAnalysis::prepare(*study, admitted);
    auto result = prepared.execute();
    const auto restored = NativeRecord::decode(result.encode());
    const auto view = restored.at("result").at("view");
    check(restored.publication_history() && view.at("features").size() == q.at("features").size(),
          "native result lost features");
    check(view.at("source_digest").text() == study->digest(), "native source digest changed");
    for (std::size_t i = 0; i < q.at("samples").size(); ++i) {
        const auto coordinate = q.at("samples").at(i).at("coordinate").integer();
        check(tensor(restored, view.at("source_coordinates")).tensor.state(i) ==
                  rfc::numeric::Rational(std::to_string(coordinate)),
              "sample coordinate map changed");
        check(view.at("samples").at(i).at("native_id").text() ==
                  study->cell_identifier({0}, static_cast<std::size_t>(coordinate)),
              "original source identity changed");
        check(tensor(restored, view.at("source_selection_map"))
                      .tensor.state(coordinate * q.at("samples").size() + i) ==
                  rfc::numeric::Rational(1),
              "native source lift changed");
    }
    const auto original = NativeRecord::decode(view.at("source_study_publication").bytes());
    const ScientificStudy cold(original);
    const auto replayed = count_model_analysis(cold, admitted);
    check(count_result_view(replayed).at("features").dump() ==
              count_result_view(restored).at("features").dump(),
          "cold native replay changed model results");
    check(study->native_publication().encode() == before, "statistics mutated source publication");
    for (std::size_t i = 0; i < q.at("features").size(); ++i) {
        const auto pin = study->field_accession(q.at("features").at(i).at("field").text());
        const auto &original = *pin.tensor(),
                   &retained = tensor(restored, view.at("features").at(i).at("input")).tensor;
        check(retained.domain().extents() == original.domain().extents(),
              "original count shape changed");
        for (std::size_t j = 0; j < original.domain().logical_size(); ++j)
            check(retained.state(j) == original.state(j), "original count support/value changed");
    }
    result.write(destination);
    std::cout << count_result_view(restored).dump() << '\n';
}
void contracts() {
    std::vector<rfc::numeric::RationalVector> rows{{1, 0}, {1, 0}, {1, 0}, {1, 1}, {1, 1}, {1, 1}};
    PreparedCountDesign design(rows);
    std::vector<double> y{1, 2, 6, 4, 9, 14}, offset(6), contrast{0, 1};
    CountWork work(10000000);
    for (double alpha : {0.0, 0.1, 2.0}) {
        const auto fit = design.fit(y, offset, alpha, contrast, {}, work);
        check(fit.converged, "two group model did not converge");
        near(fit.coefficients[0], std::log(3.0));
        near(*fit.effect, std::log(3.0));
        near(*fit.variance, (1.0 / 3 + alpha) / 3 + (1.0 / 9 + alpha) / 3);
        near(fit.means[0], 3);
        near(fit.means[5], 9);
    }
    near(count_log_mass(0, 2, 0), -2);
    near(count_log_mass(3, 2, 1), std::log(8.0 / 81));
    near(normal_critical(0.95), 1.959963984540054);
    near(chi_square_log_survival(6, 2), -3);
    check(normal_log_two_sided(40) < -800, "normal underflow lost log tail");
    const auto posterior = shrink_normal_effect(2, 3, 0, 1, .95);
    near(posterior.mean, .5);
    near(posterior.variance, .75);
    const std::vector<rfc::numeric::Rational> p{rfc::numeric::Rational(1, 100),
                                                rfc::numeric::Rational(1, 20),
                                                rfc::numeric::Rational(1, 5)};
    const auto adjusted = adjust_probabilities(p, "benjamini_hochberg");
    check(adjusted[0] == rfc::numeric::Rational(3, 100) &&
              adjusted[1] == rfc::numeric::Rational(3, 40),
          "exact BH differs");
    const auto logs =
        adjust_log_probabilities(std::vector<double>{-750, -750, -1}, "benjamini_hochberg");
    near(logs[0], -750 + std::log(1.5));
    PreparedCountDesign singular({{1, 1}, {1, 1}});
    const auto invalid =
        singular.fit(std::vector<double>{1, 2}, std::vector<double>(2), 0, {}, {}, work);
    check(!invalid.converged && invalid.status == "rank_deficient", "singular model was accepted");
    PreparedCountDesign separated({{1, 0}, {1, 0}, {1, 1}, {1, 1}});
    const auto boundary = separated.fit(std::vector<double>{0, 0, 5, 6}, std::vector<double>(4), 0,
                                        std::vector<double>{0, 1}, {}, work);
    check(!boundary.converged && boundary.status == "boundary_likelihood",
          "boundary coefficients acquired finite inference");
    PreparedCountDesign balanced({{1}, {-1}});
    const auto finite_zero = balanced.fit(std::vector<double>(2), std::vector<double>(2), 0,
                                          std::vector<double>{1}, {}, work);
    check(finite_zero.converged, "finite all zero model refused");
    near(*finite_zero.effect, 0);
    near(*finite_zero.variance, .5);
    PreparedCountDesign finite_support({{1, 0}, {-1, 0}, {0, 1}});
    check(finite_support.support_status(std::vector<double>{0, 0, 5}, work) ==
              "finite_coefficients",
          "zero observations lost a finite fit");
    const auto zero = design.fit(std::vector<double>(6), offset, 0, contrast, {}, work);
    check(!zero.converged && zero.status == "all_zero_counts",
          "zero model acquired finite coefficients");
    bool refused{};
    try {
        CountWork tiny(1);
        (void)design.fit(y, offset, 0, contrast, {}, tiny);
    } catch (const CountWorkLimit &) {
        refused = true;
    }
    check(refused, "work admission ignored");
    std::stop_source cancellation;
    cancellation.request_stop();
    refused = false;
    try {
        work.charge(0, cancellation.get_token());
    } catch (const rfc::runtime::ExecutionCancelled &) {
        refused = true;
    }
    check(refused, "cancellation ignored");
    std::cout << "count model contracts passed\n";
}
} // namespace
int main(int argc, char **argv) try {
    if (argc == 3 && std::string(argv[1]) == "--answer") {
        const auto value = Json::parse(read_text(argv[2]));
        std::cout << gadget::statistics::count_answer_view(value.at("directory").text(),
                                                           value.at("reference"))
                         .dump()
                  << '\n';
    } else if (argc == 2 && std::string(argv[1]) == "--evaluate") {
        for (std::string line; std::getline(std::cin, line);)
            std::cout << evaluate(Json::parse(line)).dump_compact() << '\n';
    } else if (argc == 4)
        native(argv[1], argv[2], argv[3]);
    else
        contracts();
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
