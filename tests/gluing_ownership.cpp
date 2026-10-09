#include <gadget/assay.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 3)
            throw std::invalid_argument("expected workspace and gluing request");
        AssayWorkspace workspace(argv[1]);
        auto request = Json::parse(read_text(argv[2]));
        const auto original = request.dump();
        auto result = workspace.execute(request);
        const auto expected = result.dump();
        auto check = [](bool condition, const char *message) {
            if (!condition) throw std::logic_error(message);
        };
        request.at("declaration").set("evidence", Json::string("changed caller"));
        request.at("sections").at(0).set("method", Json::string("changed caller"));
        request.at("embeddings").at("common_left").at("cells").at(0)
            .set("to", Json::string("changed caller"));
        check(result.dump() == expected, "caller mutation changed the sealed gluing result");
        auto binding = result.at("result").at("source_bindings").at("left");
        auto owner = workspace.study("left", 0);
        const auto index = owner->source_index().dump();
        auto read = Json::object(), stream = Json::object();
        read.set("schema", Json::string("gadget.scientific.source.read.v1"));
        read.set("source_digest", Json::string(owner->digest()));
        read.set("publication", Json::string("self"));
        read.set("source_sequence", binding.at("source_sequence"));
        stream.set("kind", Json::string("catalog")); stream.set("section", Json::string("cells"));
        read.set("stream", stream); read.set("offset", Json::integer(0));
        read.set("limit", Json::integer(1));
        const auto page = owner->read_source(read).dump();
        binding.at("native").set("sha256", Json::string("changed result"));
        binding.set("declaration_sha256", Json::string("changed result"));
        result.at("result").at("sections").at(0).at("left_meaning")
            .set("units", Json::string("changed result"));
        result.at("result").at("declaration").set("evidence", Json::string("changed result"));
        result.at("query").at("sources").at("left").set("key", Json::string("999"));
        check(owner->source_index().dump() == index, "gluing result mutated the admitted index");
        check(owner->read_source(read).dump() == page, "gluing result mutated original source records");
        request = Json::parse(original);
        check(workspace.execute(request).dump() == expected,
              "gluing result mutation changed subsequent native execution");
        check(request.dump() == original, "result mutation changed the caller request");
        std::cout << "gluing query, result and publication ownership contracts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
