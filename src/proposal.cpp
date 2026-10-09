#include <gadget/model.hpp>
#include <set>

namespace gadget {
Json apply_model_proposal(QuerySession &session, const Json &proposal) {
    auto result = Json::object();
    result.set("schema", Json::string("gadget.agent.response.v1"));
    result.set("proposal", proposal);
    try {
        if (!json_object_is_type(proposal.raw(), json_type_object))
            throw std::invalid_argument("proposal must be an object");
        if (proposal.at("schema").text() != "gadget.agent.proposal.v1")
            throw std::invalid_argument("unsupported model proposal schema");
        const auto action = proposal.at("action").text();
        const std::set<std::string> allowed =
            action == "query" ? std::set<std::string>{"schema", "action", "query"}
                              : std::set<std::string>{"schema", "action", "reason"};
        json_object_object_foreach(proposal.raw(), key, value) {
            (void)value;
            if (!allowed.contains(key))
                throw std::invalid_argument("unknown model proposal member: " + std::string(key));
        }
        if (action == "clarify") {
            const auto reason = proposal.at("reason").text();
            if (reason.empty() || reason.size() > 4096)
                throw std::invalid_argument("clarification must be 1..4096 bytes");
            result.set("status", Json::string("clarification"));
            result.set("reason", Json::string(reason));
        } else if (action == "query") {
            result.set("certificate", session.execute(proposal.at("query")));
            result.set("status", Json::string("executed"));
        } else
            throw std::invalid_argument("unsupported model action");
    } catch (const std::invalid_argument &error) {
        result.set("status", Json::string("rejected"));
        result.set("error", Json::string(error.what()));
    }
    return result;
}
} // namespace gadget
