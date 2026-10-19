#include <gadget/scientific.hpp>
#include <gadget/record_view.hpp>
#include <iostream>

using namespace gadget;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
const rfc::storage::TemporalTensor *tensor(const NativeRecord &owner, NativeRecord value) {
    const auto key = value.at("tensor").integer();
    for (const auto &item : owner.endpoint()->tensors) if (key == rfc::numeric::Integer(std::to_string(item.key))) return &item;
    throw std::runtime_error("native tensor absent");
}
}
int main(int argc, char **argv) try {
    if (argc != 4) throw std::invalid_argument("publication, queries and output required");
    const auto previous = NativeRecord::read(argv[1]);
    const auto original = NativeRecord::decode(previous.at("result").at("source_study_publication").bytes());
    const auto study = std::make_shared<const ScientificStudy>(original);
    const auto queries = admit_record_view(Json::parse(read_text(argv[2])));
    const fs::path output(argv[3]);
    check(!fs::exists(output), "output must be new"); fs::create_directories(output);
    ScientificSession session(study);
    std::vector<NativeRecord> saved;
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name);
        const auto pin = study->field_accession(query.at("field").text());
        const auto publication = session.execute_record(query);
        const auto owner = publication.at("result").at("fields");
        const auto view = owner.at("view"), binding = owner.at("bindings").at(0);
        const auto &indices = *tensor(owner,binding.at("source_coordinates"));
        const auto &input = *tensor(owner,binding.at("source_field"));
        check(binding.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
              binding.at("cutoff").integer() == rfc::numeric::Integer(std::to_string(pin.sequence())), "actual source pin changed");
        check(binding.at("input_carrier").integer() == rfc::numeric::Integer(std::to_string(*input.output)) &&
              binding.at("window_carrier").integer() == rfc::numeric::Integer(std::to_string(*indices.output)) &&
              input.output != indices.output, "window coordinates lost their independent carrier");
        check(view.at("observations").size() == indices.tensor.domain().logical_size(), "window extent lost");
        Rational sum{}; std::size_t reported{};
        const auto lo = query.at("lower").at(0).integer().get_ui();
        const auto hi = query.at("upper").at(0).integer().get_ui();
        check(view.at("observations").size() == hi-lo, "original interval changed");
        for (std::size_t i = 0; i < hi-lo; ++i) {
            const auto coordinate = lo+i;
            const auto row = view.at("observations").at(i);
            check(indices.tensor.state(i) == Rational(coordinate) &&
                  row.at("source_coordinate").integer() == coordinate, "original coordinate lift changed");
            const auto expected = pin.tensor()->state(coordinate);
            check(row.at("reported").boolean() == expected.has_value() &&
                  input.tensor.state(coordinate) == expected, "missing and supported zero merged");
            check(row.at("value").at("tensor").integer() == rfc::numeric::Integer(std::to_string(input.key)) &&
                  row.at("value").at("index").integer() == coordinate, "result copied into an unrelated field");
            if (expected) {sum += *expected; ++reported;}
        }
        check(view.at("reported").integer() == reported && view.at("missing").integer() == hi-lo-reported,
              "complete source support counts changed");
        const auto &total = *tensor(owner,view.at("sum"));
        check(total.tensor.support().contains(0) == (reported != 0), "undefined sum became supported zero");
        if (reported) check(total.tensor.state(0) == sum, "independent exact sum differs");
        check(NativeRecord::decode(publication.at("result").at("source_study_publication").bytes()).encode() ==
              original.encode(), "original archive bytes changed");
        publication.write(output/(name+".rfcxs"));
        check(NativeRecord::read(output/(name+".rfcxs")).encode() == publication.encode(), "archive reopen changed bytes");
        saved.push_back(publication);
    }
    session.clear_tensor_window_preparations();
    check(session.tensor_window_preparation_statistics().at("retained_preparations").integer() == 0,
          "window cache not released");
    for (const auto &publication : saved) {
        const auto replayed = replay_scientific_record(publication);
        check(publication.endpoint()->tensors.size() == replayed.endpoint()->tensors.size(), "replay lost a field");
        for (std::size_t i = 0; i < publication.endpoint()->tensors.size(); ++i)
            check(rfc::storage::temporal_endpoint_detail::tensor_equal(publication.endpoint()->tensors[i],
                  replayed.endpoint()->tensors[i]), "portable replay altered native support or values");
        check(std::ranges::equal(publication.at("request_publication").bytes(),replayed.at("request_publication").bytes()),
              "replay replaced the original request owner");
    }
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    return 0;
} catch (const std::exception &error) {std::cerr << error.what() << '\n'; return 1;}
