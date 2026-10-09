#include "fixture.hpp"
#include <iostream>
int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::invalid_argument("usage: gadget_scientific_example NEW_DIRECTORY");
        gadget::example::scientific_study(argv[1]);
        auto study = std::make_shared<const gadget::ScientificStudy>(gadget::fs::path(argv[1]) /
                                                                     "study.json");
        auto plan = gadget::Json::object();
        plan.set("schema", gadget::Json::string("gadget.scientific.plan.v1"));
        plan.set("id", gadget::Json::string("synthetic-get-change"));
        plan.set("question", study->declaration().at("design").at("question"));
        plan.set("source_digest", gadget::Json::string(study->digest()));
        plan.set("steps", gadget::Json::parse(R"([
          {"id":"baseline-GET","query":{"operation":"typed_moment","context":"baseline","field":"baseline.G","other_field":"baseline.E","order":1}},
          {"id":"longitudinal-change","query":{"operation":"temporal_moment_change","context":"baseline","other_context":"followup","field":"baseline.G","other_field":"baseline.E","new_field":"followup.G","new_other_field":"followup.E","order":1}},
          {"id":"missingness","query":{"operation":"coverage","field":"baseline.partial"}},
          {"id":"declared-batch","query":{"operation":"coverage","field":"baseline.batch"}},
          {"id":"native-CPU","query":{"operation":"native_boundary","field_set":"synthetic.assays","backend":"cpu_tiles","workers":4,"tile_rows":1}}
        ])"));
        gadget::write_text(gadget::fs::path(argv[1]) / "plan.json", plan.dump() + "\n");
        std::cout << study->digest() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
