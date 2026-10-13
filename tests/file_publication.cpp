#include <rfc/storage/file_io.hpp>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;
using namespace rfc::storage;
void check(bool condition, const char *message) {
    if (!condition) throw std::logic_error(message);
}
template<class Action> void refuses(Action action) {
    bool refused = false;
    try { action(); } catch (const std::exception &) { refused = true; }
    check(refused, "file publication admitted a failed producer");
}
void clean(const fs::path &root) {
    for (const auto &entry : fs::directory_iterator(root))
        check(!entry.path().filename().string().starts_with(".rfc-stage-"),
              "ordinary failure retained a staging directory");
}
int main(int argc, char **argv) {
    try {
        check(argc == 2, "new owned fixture directory required");
        const fs::path root(argv[1]);
        check(!fs::exists(root), "owned fixture directory already exists");
        fs::create_directories(root);
        const auto data = std::as_bytes(std::span("original\0complete", 17));
        const auto path = root / "publication";
        refuses([&] {
            publish_file_exact(path, FilePublication::CreateNew, [&](std::ostream &out) {
                write_stream_exact(out, data.first(4));
                throw std::runtime_error("controlled producer failure");
            });
        });
        check(!fs::exists(path), "producer exception published a partial destination"); clean(root);
        write_file_exact(path, data, FilePublication::CreateNew);
        check(read_file_exact(path, data.size()) == std::vector(data.begin(), data.end()),
              "complete source bytes changed after retry"); clean(root);
        bool ran = false;
        refuses([&] { publish_file_exact(path, FilePublication::CreateNew,
                                        [&](std::ostream &) { ran = true; }); });
        check(!ran, "existing destination reached the producer");
        refuses([&] { publish_file_exact(path, FilePublication::ReplaceRegular,
                                        [&](std::ostream &out) { out.setstate(std::ios::badbit); }); });
        check(read_file_exact(path, data.size()) == std::vector(data.begin(), data.end()),
              "failed replacement altered an existing publication"); clean(root);
        const auto winner = root / "winner";
        refuses([&] {
            publish_file_exact(winner, FilePublication::CreateNew, [&](std::ostream &out) {
                write_stream_exact(out, data);
                write_file_exact(winner, data.first(3), FilePublication::CreateNew);
            });
        });
        check(read_file_exact(winner, 3) == std::vector(data.begin(), data.begin() + 3),
              "late exclusive creation replaced another publication"); clean(root);
        for (const auto &invalid : {fs::path{}, fs::path(std::string("prefix\0suffix", 13))})
            refuses([&] { write_file_exact(invalid, data, FilePublication::CreateNew); });
        refuses([&] { (void)read_file_exact(path, data.size() - 1); });
        std::ostringstream failed;
        failed.setstate(std::ios::badbit);
        refuses([&] { write_stream_exact(failed, data); });
        std::cout << "PASS: canonical producer/write/flush failure, complete retry, existing and late destination preservation, path and actual byte representation controls\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
