#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected native cutoff case directory");
        const fs::path root = argv[1];
        ScientificStudy first(root / "p0/study.json"), middle(root / "p1/study.json");
        const auto prefix = middle.native_history();
        if (first.native_history()->segment_count() != 0 || prefix->segment_count() != 1 ||
            prefix->first_sequence() != 0 || prefix->cutoff() != 1 ||
            prefix->tail().sequence() != middle.endpoint()->sequence)
            throw std::runtime_error("native cutoff ownership mismatch");
        bool refused = false;
        try {
            (void)prefix->at(2);
        } catch (const std::invalid_argument &) {
            refused = true;
        }
        if (!refused)
            throw std::runtime_error("historical prefix exposed a later endpoint");
        const auto interval = prefix->difference(0, 1);
        if (interval.segment_count() != 1 || prefix->segment_bytes(0).empty())
            throw std::runtime_error("historical native frame was lost");
        const auto file = root / "p1/history.rfcxs";
        const auto bytes = read_text(file);
        struct Restore {
            fs::path file;
            std::string bytes;
            ~Restore() { write_text(file, bytes); }
        } restore{file, bytes};
        write_text(file, "changed after admission");
        if (prefix->at(0).sequence() != 0 || prefix->tail().sequence() != 1)
            throw std::runtime_error("historical owner depends on mutable archive file");
        refused = false;
        try {
            ScientificStudy changed(root / "p1/study.json");
        } catch (const std::exception &) {
            refused = true;
        }
        if (!refused)
            throw std::runtime_error("fresh admission accepted an altered native archive");
        std::cout << "PASS: exact cutoff, actual frames, exclusion of later endpoints and owning "
                     "replay\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
