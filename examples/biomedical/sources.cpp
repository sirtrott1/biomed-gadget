#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("usage: gadget_biomedical_sources STUDY REQUEST NEW_RECORD");
        const gadget::fs::path destination(argv[3]);
        if (gadget::fs::exists(destination))
            throw std::invalid_argument("source record destination exists");
        const gadget::ScientificStudy study(argv[1]);
        const auto request = gadget::read_record(argv[2]).native();
        const auto publication = gadget::retain_study_sources(study, request);
        const gadget::ScientificStudy admitted(publication);
        if (admitted.digest() != study.digest())
            throw std::runtime_error("original provenance changed admitted measurements");
        publication.write(destination);
        std::cout << "Retained original source bytes and admitted study\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
