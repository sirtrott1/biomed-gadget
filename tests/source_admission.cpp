#include <gadget/scientific.hpp>
#include <iostream>
#include <thread>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected version-two study");
        const fs::path path = argv[1];
        auto owner = std::make_shared<const ScientificStudy>(path);
        const auto original = read_text(path.parent_path() / "table.txt");
        if (sha256_file(path.parent_path() / "table.txt") != sha256(original))
            throw std::logic_error("incremental identity differs from byte identity");
        const auto index = owner->source_index();
        auto request = Json::object();
        request.set("schema", Json::string("gadget.scientific.source.read.v1"));
        request.set("source_digest", index.at("source_digest"));
        request.set("publication", Json::string("self"));
        request.set("source_sequence", index.at("source_sequence"));
        request.set("stream", Json::parse(R"({"kind":"document","id":"table"})"));
        request.set("offset", Json::integer(0));
        request.set("limit", Json::integer(original.size()));
        const auto expected = owner->read_source(request).dump();
        write_text(path.parent_path() / "table.txt", "replacement after admission\n");
        std::exception_ptr failure[2];
        std::thread first([&] {
            try {
                if (owner->read_source(request).dump() != expected)
                    throw std::logic_error("admitted source followed external file mutation");
            } catch (...) {
                failure[0] = std::current_exception();
            }
        });
        std::thread second([&] {
            try {
                ScientificSession session(owner);
                if (session.read_source(request).dump() != expected)
                    throw std::logic_error("concurrent immutable source changed");
            } catch (...) {
                failure[1] = std::current_exception();
            }
        });
        first.join();
        second.join();
        write_text(path.parent_path() / "table.txt", original);
        for (const auto &error : failure)
            if (error)
                std::rethrow_exception(error);
        auto detached = owner->declaration();
        detached.at("cells").at(0).set("id", Json::string("changed"));
        if (owner->declaration().at("cells").at(0).at("id").text() == "changed")
            throw std::logic_error("detached declaration changed owner");
        std::cout << "PASS: streaming file identity, immutable admitted bytes, concurrent readers "
                     "and detached metadata\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
