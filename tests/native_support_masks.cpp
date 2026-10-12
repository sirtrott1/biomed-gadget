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
const rfc::storage::TemporalTensor *tensor(const NativeRecord &owner, NativeRecord descriptor) {
    for (const auto &value : owner.endpoint()->tensors)
        if (descriptor.at("tensor").integer() == value.key) return &value;
    throw std::runtime_error("native mask tensor absent");
}
}
int main(int argc, char **argv) try {
    if (argc != 4) throw std::invalid_argument("publication, queries and output required");
    const auto previous = NativeRecord::read(argv[1]);
    const auto original = NativeRecord::decode(previous.at("result").at("source_study_publication").bytes());
    const auto study = std::make_shared<const ScientificStudy>(original);
    ScientificSession session(study);
    const auto queries = admit_record_view(Json::parse(read_text(argv[2])));
    const fs::path output(argv[3]);
    check(!fs::exists(output), "output must be new"); fs::create_directories(output);
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name), ids = query.at("fields");
        std::vector<rfc::storage::PinnedNativeAccession> pins;
        for (std::size_t i = 0; i < ids.size(); ++i) pins.push_back(study->field_accession(ids.at(i).text()));
        const auto publication = session.execute_record(query);
        const auto owner = publication.at("result").at("fields"), view = owner.at("view");
        const auto binding = owner.at("bindings").at(0), inputs = binding.at("inputs");
        const auto &mask = *tensor(owner,binding.at("output_masks").at(0));
        const auto n = mask.tensor.domain().logical_size();
        check(inputs.size() == pins.size() && mask.tensor.support().count() == n,
              "complete eligibility domain or original planes lost");
        std::size_t selected{};
        for (std::size_t i = 0; i < n; ++i) {
            bool all = true, any = false;
            for (std::size_t f = 0; f < pins.size(); ++f) {
                const auto &pin = pins[f]; const auto input = inputs.at(f);
                const auto &plane = *tensor(owner,input.at("support_plane"));
                const auto present = pin.tensor()->support().contains(i);
                check(input.at("id").text() == ids.at(f).text() &&
                      input.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
                      input.at("cutoff").integer() == pin.sequence(), "original ordered field owner changed");
                check(plane.output == mask.output && plane.tensor.support().contains(i) &&
                      plane.tensor.state(i) == Rational(present ? 1 : 0), "original support plane changed");
                const auto word = input.at("packed_words").at(i/64).integer().get_ui();
                check(((word >> (i%64)) & 1u) == present, "actual packed input word differs");
                const auto &source = *tensor(owner,input.at("source_field"));
                check(source.tensor.state(i) == pin.tensor()->state(i), "original supported zero or missing value changed");
                all = all && present; any = any || present;
            }
            const auto mode = query.at("mode").text();
            const bool expected = mode == "all_reported" ? all : mode == "any_reported" ? any : !any;
            check(mask.tensor.state(i) == Rational(expected ? 1 : 0), "independent complete support oracle differs");
            const auto word = binding.at("output_words").at(0).at(i/64).integer().get_ui();
            check(((word >> (i%64)) & 1u) == expected, "actual packed output word differs");
            if (expected) {
                check(view.at("selection").at(selected).at("coordinate").integer() == i,
                      "selected source coordinate order changed");
                ++selected;
            }
        }
        check(view.at("coordinates").integer() == n && view.at("selected").integer() == selected &&
              view.at("selection").size() == selected, "complete selected family changed");
        const auto instructions = binding.at("instructions");
        const auto count = 1+2*(ids.size()-1)+(query.at("mode").text() == "none_reported" ? 1 : 0);
        check(instructions.size() == count && binding.at("outputs").size() == 1 &&
              binding.at("outputs").at(0).integer() == 0, "ordered instruction or output occurrence lost");
        check(instructions.at(0).at("opcode").integer() == 0 && instructions.at(0).at("aspect").integer() == 0,
              "initial original support plane load changed");
        for (std::size_t i = 1; i < ids.size(); ++i) {
            const auto load = instructions.at(2*i-1), combine = instructions.at(2*i);
            check(load.at("opcode").integer() == 0 && load.at("dst").integer() == 1 &&
                  load.at("aspect").integer() == i, "ordered support aspect changed");
            check(combine.at("opcode").integer() == (query.at("mode").text() == "all_reported" ? 4 : 5) &&
                  combine.at("a").integer() == 0 && combine.at("b").integer() == 1 &&
                  combine.at("dst").integer() == 0, "actual combine instruction changed");
        }
        if (query.at("mode").text() == "none_reported")
            check(instructions.at(count-1).at("opcode").integer() == 3, "final support complement lost");
        check(NativeRecord::decode(publication.at("result").at("source_study_publication").bytes()).encode() ==
              original.encode(), "original source archive changed");
        publication.write(output/(name+".rfcxs"));
        check(NativeRecord::read(output/(name+".rfcxs")).encode() == publication.encode(), "archive reopen changed bytes");
        const auto repeated = session.execute_record(query);
        check(repeated.at("result").at("fields").at("execution")
                  .at("prepared_operation_reused").boolean(), "actual prepared owner not reused");
        for (auto replay = replay_scientific_record(publication);;) {
            check(replay.endpoint()->tensors.size() == publication.endpoint()->tensors.size(), "portable replay lost complete masks");
            for (std::size_t i = 0; i < publication.endpoint()->tensors.size(); ++i)
                check(rfc::storage::temporal_endpoint_detail::tensor_equal(replay.endpoint()->tensors[i],
                      publication.endpoint()->tensors[i]), "portable replay changed original support or exact output");
            check(std::ranges::equal(replay.at("request_publication").bytes(),publication.at("request_publication").bytes()),
                  "portable replay changed the request owner");
            break;
        }
    }
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    return 0;
} catch (const std::exception &error) {std::cerr << error.what() << '\n'; return 1;}
