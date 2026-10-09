// Standalone acceptance consumer for the sealed I-SPY2 relational recipe.
// Its input frames and lifted fields were produced natively and are separately
// checked against the XLSX source by the independent oracle.
#include <gadget/io.hpp>
#include <iostream>
#include <rfc/field/exact_hodge.hpp>

int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::invalid_argument("usage: relational_winding RESULTS_JSON");
        using namespace gadget;
        using namespace rfc;
        const auto bytes = read_text(argv[1]);
        const auto result = Json::parse(bytes);
        auto open = result.at("topology").at("unfilled_design");
        auto design = result.at("topology").at("design");
        if (open.at("to_ambient").at(1).dump() != design.at("to_ambient").at(1).dump())
            throw std::invalid_argument("different primary relation carriers");
        auto frame = open.at("harmonic_frame");
        const auto rows = frame.size(), columns = frame.at(0).size();
        if (columns != 6)
            throw std::invalid_argument("expected six open-design harmonic directions");
        numeric::ExactMatrix harmonic(rows, columns);
        Vector weights;
        auto map = open.at("to_ambient").at(1);
        for (std::size_t i = 0; i < rows; ++i) {
            const auto index = static_cast<std::size_t>(map.at(i).integer());
            weights.emplace_back(
                result.at("declaration").at("relations").at(index).at("occurrences").size());
            for (std::size_t j = 0; j < columns; ++j)
                harmonic.set(i, j, Rational{frame.at(i).at(j).text()});
        }
        const auto metric = field::ExactMetric::diagonal(weights);
        auto readings = Json::array();
        auto fields = result.at("fields");
        for (std::size_t i = 0; i < fields.size(); ++i) {
            auto f = fields.at(i);
            auto v = f.at("lifted_field");
            Vector lifted;
            for (std::size_t j = 0; j < v.size(); ++j)
                lifted.emplace_back(v.at(j).text());
            auto pairing = field::metric_frame_dual_apply(harmonic, metric, lifted);
            if (!numeric::is_zero(pairing))
                throw std::logic_error("coboundary has nonzero harmonic pairing");
            auto record = Json::object(), values = Json::array();
            record.set("field", f.at("id"));
            for (const auto &q : pairing)
                values.push(Json::string(numeric::str(q)));
            record.set("pairings", values);
            readings.push(record);
        }
        // Positive control: a primary relation direction genuinely excites
        // this frame. Empty harmonic output cannot pass the test.
        bool excited = false;
        for (std::size_t i = 0; i < rows && !excited; ++i) {
            Vector seed(rows, Rational{0});
            seed[i] = 1;
            excited = !numeric::is_zero(field::metric_frame_dual_apply(harmonic, metric, seed));
        }
        if (!excited)
            throw std::logic_error("harmonic positive control failed");
        auto out = Json::object();
        out.set("status", Json::string("PASS"));
        out.set("result_sha256", Json::string(sha256(bytes)));
        out.set("harmonic_directions", Json::integer(static_cast<std::int64_t>(columns)));
        out.set("exact_pairings_checked",
                Json::integer(static_cast<std::int64_t>(columns * fields.size())));
        out.set("primary_seed_positive_control", Json::boolean(true));
        out.set("fields", readings);
        std::cout << out.dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
