#include <gadget/model.hpp>
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
        options.threads = 2;
        options.paged_catalog = true;
        LocalAgent agent(source, options);
        ScientificSession exact(source);
        auto receipts = Json::array();
        for (std::size_t i = 0; i < 3; ++i) {
            auto requested = Json::parse(query.dump());
            if (i == 0) {
                json_object_object_del(requested.raw(), "causes");
            } else if (i == 2)
                requested.set("preparation", Json::string("ephemeral"));
            auto proposal = Json::object();
            proposal.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
            proposal.set("action", Json::string("query"));
            proposal.set("query", requested);
            const auto receipt = agent.ask(
                "Execute this already declared source-bound comparison. "
                "Return precisely this proposal, retaining all comparison entries, field ids, "
                "cause slots and evidence bytes: " +
                proposal.dump());
            if (receipt.at("status").text() != "executed" ||
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
