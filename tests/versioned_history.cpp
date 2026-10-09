#include <gadget/scientific.hpp>
#include <gadget/study.hpp>
#include <iostream>
#include <rfc/temporal/transfer.hpp>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected qualification directory");
        const fs::path root = argv[1];
        auto study = std::make_shared<const ScientificStudy>(root / "study.json");
        ScientificSession session(study);
        const auto query = Json::parse(read_text(root / "history-query.json"));
        const auto before = session.execute(query).dump();
        const auto preparation = session.temporal_preparation_statistics().dump();
        if (session.temporal_preparation_statistics().at("green_preparations").integer() != 3 ||
            session.temporal_preparation_statistics().at("sector_preparations").integer() != 3)
            throw std::runtime_error("temporal endpoint preparation was repeated within the history");
        if (session.execute(query).dump() != before ||
            session.temporal_preparation_statistics().dump() != preparation)
            throw std::runtime_error("temporal preparation or certificate changed on session reuse");
        auto equivalent = Json::parse(query.dump());
        equivalent.set("lambda", Json::string("1/1"));
        for (std::size_t i = 0; i < equivalent.at("word").size(); ++i)
            if (equivalent.at("word").at(i).at("kind").text() == "green")
                equivalent.at("word").at(i).set("lambda", Json::string("1/1"));
        (void)session.execute(equivalent);
        if (session.temporal_preparation_statistics().dump() != preparation)
            throw std::runtime_error("equivalent rational parameter duplicated temporal factors");
        if (session.preparation_count() != 3)
            throw std::runtime_error("publication preparations omitted from session accounting");
        // Compare the composed operator-first path with the SDK's independent
        // materialized Green-delta reference on both original source grades.
        for (const auto &file_name : {"certificate.json", "relations.json"}) {
            const auto result = Json::parse(read_text(root / file_name)).at("result");
            for (std::size_t t = 0; t < 2; ++t) {
                ScientificStudy a(root / ("p" + std::to_string(t)) / "study.json");
                ScientificStudy b(root / ("p" + std::to_string(t + 1)) / "study.json");
                const auto old_endpoint = a.endpoint(), new_endpoint = b.endpoint();
                std::vector<rfc::field::ExactMetric> old_metrics, new_metrics;
                for (std::size_t g = 0; g < old_endpoint->source.grade_count(); ++g) {
                    const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(g)};
                    old_metrics.push_back(rfc::field::ExactMetric::identity(
                        old_endpoint->source.grade(grade).cell_count()));
                    new_metrics.push_back(rfc::field::ExactMetric::identity(
                        new_endpoint->source.grade(grade).cell_count()));
                }
                NativeStudy old_study(old_endpoint->source, old_metrics),
                    new_study(new_endpoint->source, new_metrics);
                const auto old = result.at("endpoints").at(t),
                           current = result.at("endpoints").at(t + 1);
                const auto grade =
                    rfc::core::Grade{static_cast<std::uint16_t>(old.at("grade").integer())};
                const auto decode = [](const Json &v) {
                    Vector x;
                    for (std::size_t i = 0; i < v.size(); ++i)
                        x.push_back(exact_rational(v.at(i).text()));
                    return x;
                };
                const auto x = decode(old.at("field")), y = decode(current.at("field"));
                const auto step = result.at("steps").at(t),
                           encoded = step.at("source_maps").at(grade.value);
                rfc::numeric::ExactMatrix j(
                    static_cast<std::size_t>(encoded.at("rows").integer()),
                    static_cast<std::size_t>(encoded.at("columns").integer()));
                for (std::size_t i = 0; i < encoded.at("entries").size(); ++i) {
                    const auto e = encoded.at("entries").at(i);
                    j.set(static_cast<std::size_t>(e.at("row").integer()),
                          static_cast<std::size_t>(e.at("column").integer()),
                          exact_rational(e.at("value").text()));
                }
                const auto delta = rfc::temporal::green_response_delta(
                    old_study.hodge(grade).materialize(), new_study.hodge(grade).materialize(), j,
                    rfc::numeric::ExactMatrix::identity(x.size()),
                    rfc::numeric::ExactMatrix::identity(y.size()), j, x, y, Rational{1});
                const auto structure = step.at("structure");
                if (decode(structure.at("realization_term")) != delta.realization_term ||
                    decode(structure.at("innovation_term")) != delta.innovation_term ||
                    decode(structure.at("operator_term")) != delta.operator_term ||
                    decode(structure.at("total")) != delta.aggregate)
                    throw std::runtime_error("operator-first Green history differs from native "
                                             "GreenResponseDelta reference");
            }
        }
        // The admitted children own all bytes. Later file edits cannot substitute
        // a new publication under the old certificate or invalidate its owners.
        const auto file = root / "p1/study.json";
        const auto bytes = read_text(file);
        struct Restore {
            fs::path file;
            std::string bytes;
            ~Restore() { write_text(file, bytes); }
        } restore{file, bytes};
        write_text(file, "{}\n");
        const auto after = session.execute(query).dump();
        if (before != after)
            throw std::runtime_error("retained publication changed after file edit");
        const auto original_catalog = study->catalog().dump();
        auto catalog = study->catalog();
        catalog.at("publications")
            .at(1)
            .at("catalog")
            .set("source_digest", Json::string("changed"));
        if (study->catalog().dump() != original_catalog)
            throw std::runtime_error("detached publication catalog mutated source owner");
        if (session.execute(query).dump() != before)
            throw std::runtime_error("source certificate changed after detached catalog edit");
        bool refused = false;
        try {
            ScientificStudy fresh(root / "study.json");
        } catch (const std::exception &) {
            refused = true;
        }
        if (!refused)
            throw std::runtime_error("fresh admission accepted an altered publication");
        std::cout << "PASS: owning versioned publications, detached catalogs, preparation reuse "
                     "and fresh hash refusal\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
