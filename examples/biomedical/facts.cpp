#include <gadget/facts.hpp>
#include <gadget/record_view.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--inspect") {
            const auto operation = gadget::NativeRecord::read(argv[2]);
            const auto result = operation.at("result");
            const auto metadata = result.at("schema").text() == "gadget.scientific.result.fields.v1"
                ? result.at("certificate") : result;
            auto info = gadget::Json::object();
            info.set("source_namespace", gadget::Json::string(std::to_string(operation.source_namespace())));
            info.set("cutoff", gadget::Json::integer(static_cast<std::int64_t>(operation.publication_history()->last_sequence())));
            info.set("source_digest", gadget::Json::string(metadata.at("source_digest").text()));
            std::cout << info.dump() << '\n';
            return 0;
        }
        if (argc != 4 && argc != 5)
            throw std::invalid_argument("usage: gadget_native_facts OPERATION REQUEST NEW_DIRECTORY [CLAIMS]");
        const gadget::fs::path destination(argv[3]);
        if (gadget::fs::exists(destination))
            throw std::invalid_argument("native fact destination exists");
        const auto operation = gadget::NativeRecord::read(argv[1]);
        const auto request = gadget::read_record(argv[2]).native();
        const gadget::ScientificFacts facts(operation, request);
        gadget::fs::create_directories(destination);
        facts.native_publication().write(destination / "facts.rfcxs");
        gadget::write_text(destination / "catalog.json", facts.catalog().dump() + "\n");
        gadget::write_text(destination / "reading.json", facts.model_view().dump() + "\n");
        if (argc == 5)
            gadget::write_text(destination / "checked.json",
                facts.check(gadget::Json::parse(gadget::read_text(argv[4]))).dump() + "\n");
        std::cout << "Retained native fact selections and reading view\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
