#include <gadget/model.hpp>
#include <gadget/scientific.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 4)
            throw std::invalid_argument("expected study, model and receipt output");
        auto study = std::make_shared<const ScientificStudy>(argv[1]);
        LocalModelOptions options;
        options.model = argv[2];
        options.context_tokens = 8192;
        options.output_tokens = 768;
        options.threads = 2;
        options.paged_catalog = true;
        LocalAgent agent(study, options);
        const auto index = study->source_index();
        auto request = Json::object();
        request.set("schema", Json::string("gadget.scientific.source.read.v1"));
        request.set("source_digest", index.at("source_digest"));
        request.set("publication", Json::string("self"));
        request.set("source_sequence", index.at("source_sequence"));
        request.set("stream", Json::parse(R"({"kind":"field","id":"x"})"));
        request.set("offset", Json::integer(0));
        request.set("limit", Json::integer(2));
        auto proposal = Json::object();
        proposal.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
        proposal.set("action", Json::string("retrieve"));
        proposal.set("read", request);
        auto check = [](bool condition, const char *message) {
            if (!condition)
                throw std::logic_error(message);
        };
        auto receipts = Json::array();
        const auto first = agent.ask("Read the first two original x coordinates. Return this "
                                     "exact admitted retrieval proposal: " +
                                     proposal.dump());
        check(first.at("status").text() == "retrieved", "model did not retrieve first page");
        check(first.at("certificate").dump() == study->read_source(request).dump(),
              "model changed requested first page");
        check(!first.at("model_context").at("catalog").has("retrieved_page"),
              "first context included future receipt");
        receipts.push(first);
        const auto second = agent.ask("Continue the same coordinate read. Return action retrieve "
                                      "using the exact next_request in retrieved_page.");
        check(second.at("status").text() == "retrieved", "model did not continue source read");
        const auto next = first.at("certificate").at("next_request");
        check(second.at("certificate").dump() == study->read_source(next).dump(),
              "model changed original continuation");
        check(second.at("model_context").at("catalog").at("retrieved_page").dump() ==
                  first.at("certificate").dump(),
              "second model context lost actual receipt");
        check(second.at("certificate").at("complete").boolean(), "continuation did not finish");
        receipts.push(second);
        const auto third =
            agent.ask("Now measure the complete reported support of field x. "
                      "Return action query with query {\"operation\":\"coverage\","
                      "\"field\":\"x\"}. This operates on the complete native source.");
        ScientificSession exact(study);
        check(third.at("status").text() == "executed", "model did not request native coverage");
        check(third.at("certificate").dump() ==
                  exact.execute(Json::parse(R"({"operation":"coverage","field":"x"})")).dump(),
              "retrieved pages replaced native source");
        receipts.push(third);
        try {
            (void)agent.ask(std::string(100000, 'x'));
            throw std::logic_error("caller context overflow was accepted");
        } catch (const std::invalid_argument &) {
        }
        agent.clear_source_context();
        const auto cleared =
            agent.ask("Return action query with query {\"operation\":\"inspect\"}.");
        check(!cleared.at("model_context").at("catalog").has("retrieved_page"),
              "context clear retained old receipt");
        receipts.push(cleared);
        agent.bind_scientific(std::make_shared<const ScientificStudy>(argv[1]));
        const auto rebound =
            agent.ask("Return action query with query {\"operation\":\"inspect\"}.");
        check(!rebound.at("model_context").at("catalog").has("retrieved_page"),
              "publication bind retained preceding receipt");
        receipts.push(rebound);
        auto result = Json::object();
        result.set("status", Json::string("PASS"));
        result.set("receipts", receipts);
        write_text(argv[3], result.dump() + "\n");
        std::cout << "local source retrieval, continuation, native query and bind reset passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
