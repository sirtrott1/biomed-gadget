#include <charconv>
#include <gadget/relational.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4 && argc != 5)
            throw std::invalid_argument(
                "usage: gadget_relational INPUT_JSON NEW_DIRECTORY cohort_design_v1 [WORKERS]");
        std::size_t workers = 1;
        if (argc == 5) {
            const std::string text = argv[4];
            auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), workers);
            if (error != std::errc{} || end != text.data() + text.size())
                throw std::invalid_argument("invalid worker count");
        }
        // argv[3] is deliberately reserved for an explicit profile token.
        if (std::string(argv[3]) != "cohort_design_v1")
            throw std::invalid_argument("unknown relational audit profile");
        gadget::run_relational_audit(argv[1], argv[2], workers);
        std::cout << "completed retained-complex audit\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
