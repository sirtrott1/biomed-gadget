#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 5)
            throw std::invalid_argument("usage: gadget_biomedical_fields EVIDENCE STUDY REQUEST NEW_DIRECTORY");
        const gadget::fs::path destination(argv[4]);
        if (gadget::fs::exists(destination))
            throw std::invalid_argument("biomedical field destination exists");
        const gadget::EvidenceStore evidence(argv[1]);
        const gadget::ScientificStudy study(argv[2]);
        const auto request = gadget::read_record(argv[3]).native();
        const auto binding = evidence.bind_reported_field(study, request);
        const auto result = binding.at("result");
        const auto field = gadget::NativeRecord::decode(result.at("field_publication").bytes());
        const auto selector = binding.derive(result.at("selector").unbound_value());
        gadget::fs::create_directories(destination);
        binding.write(destination / "binding.rfcxs");
        field.write(destination / "field-source.rfcxs");
        selector.write(destination / "selector.rfcxs");
        std::cout << "Retained biomedical binding, field source and native selector\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
