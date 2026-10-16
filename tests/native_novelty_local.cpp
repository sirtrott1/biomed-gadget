#include <gadget/model.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 4)
            throw std::invalid_argument("expected case, local GGUF model, output");
        const fs::path root = argv[1];
        auto source = std::make_shared<const ScientificStudy>(root / "study.json");
        auto query = Json::parse(read_text(root / "query.json"));
        LocalModelOptions options;
        options.model = argv[2];
        options.context_tokens = 8192;
        options.output_tokens = 1536;
        options.threads = 8;
        options.paged_catalog = true;
        auto bindings = Json::array();
        for (const auto *name : {"old", "new"}) {
            const auto endpoint = query.at(name);
            auto binding = Json::object();
            binding.set("publication", endpoint.at("publication"));
            binding.set("source_sequence", endpoint.at("source_sequence"));
            for (const auto *section : {"fields", "operators"}) {
                auto read = Json::object();
                read.set("schema", Json::string("gadget.scientific.source.read.v1"));
                read.set("source_digest", Json::string(source->digest()));
                read.set("publication", endpoint.at("publication"));
                read.set("source_sequence", endpoint.at("source_sequence"));
                auto stream = Json::object();
                stream.set("kind", Json::string("catalog"));
                stream.set("section", Json::string(section));
                read.set("stream", stream);
                read.set("offset", Json::integer(0));
                read.set("limit", Json::integer(64));
                auto descriptors = Json::array();
                while (true) {
                    const auto page = source->read_source(read);
                    for (std::size_t j = 0; j < page.at("records").size(); ++j) {
                        const auto record = page.at("records").at(j);
                        auto item = Json::object();
                        for (const auto *key :
                             {"id", "key", "carrier", "input_carrier", "output_carrier", "extents",
                              "value_type", "quantity", "units", "grade", "context", "method"})
                            if (record.has(key))
                                item.set(key, record.at(key));
                        descriptors.push(item);
                    }
                    if (!page.has("next_request") ||
                        json_object_get_type(page.at("next_request").raw()) == json_type_null)
                        break;
                    read = page.at("next_request");
                }
                binding.set(section, descriptors);
            }
            bindings.push(binding);
        }
        LocalAgent agent(source, options);
        ScientificSession exact(source);
        auto receipts = Json::array();
        for (std::size_t i = 0; i < 3; ++i) {
            auto requested = Json::parse(query.dump());
            if (i == 0)
                json_object_object_del(requested.raw(), "causes");
            else if (i == 2)
                requested.set("preparation", Json::string("ephemeral"));
            auto proposal = Json::object();
            proposal.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
            proposal.set("action", Json::string("query"));
            proposal.set("query", requested);
            auto read = Json::object();
            read.set("schema", Json::string("gadget.scientific.source.read.v1"));
            read.set("source_digest", Json::string(source->digest()));
            read.set("publication", query.at("new").at("publication"));
            read.set("source_sequence", query.at("new").at("source_sequence"));
            const auto evidence = query.at("comparison").at("evidence").at(0);
            auto stream = Json::object();
            stream.set("kind", Json::string("document"));
            stream.set("id", evidence.at("document"));
            read.set("stream", stream);
            read.set("offset", evidence.at("start"));
            read.set("limit",
                     Json::integer(evidence.at("end").integer() - evidence.at("start").integer()));
            const auto page = source->read_source(read);
            const auto prompt =
                "The caller has already loaded the actual field and operator catalogue IDs and the "
                "complete declared comparison evidence below. "
                "All supplied IDs are admitted; the evidence is original source data. No further "
                "source read is needed for this already declared query. "
                "Return action query using exactly the supplied proposal. Preserve every ordered "
                "cause including repetitions, every map entry and all evidence offsets. "
                "Publication IDs are p0 and p1. declarations is the document ID, not a "
                "publication. "
                "Actual source catalogue descriptors=" +
                bindings.dump() + "\nComplete original comparison evidence=" + page.dump() +
                "\nExact requested query proposal=" + proposal.dump();
            const auto native = agent.ask_record(prompt);
            const auto receipt = render_record_view(
                native.at("schema").text() == "gadget.operation.native.v1" ? native.at("result")
                                                                           : native);
            if (receipt.at("status").text() != "executed" ||
                !json_object_equal(receipt.at("proposal").at("query").raw(), requested.raw()) ||
                receipt.at("certificate").dump() != exact.execute(requested).dump())
                throw std::logic_error("local decoder changed the complete admitted native query");
            receipts.push(receipt);
        }
        auto out = Json::object();
        out.set("status", Json::string("PASS"));
        out.set("purpose", Json::string("local structured proposal/execution integration"));
        out.set("receipts", receipts);
        write_text(argv[3], out.dump() + "\n");
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
