#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4 && argc != 5)
            throw std::invalid_argument("usage: proposals evidence-directory context-request output [proposal]");
        const gadget::EvidenceStore evidence(argv[1]);
        const auto request = gadget::read_record(argv[2]).native();
        auto result = argc == 4
            ? evidence.source_context(request)
            : evidence.admit_proposal(request, gadget::read_record(argv[4]).native(),
                                      request.at("question").text());
        result.write(argv[3]);
        std::cout << gadget::render_record_view(result.at("result")).dump() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
