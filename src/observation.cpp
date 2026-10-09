#include "scientific_detail.hpp"
#include <gadget/biomedical.hpp>

namespace gadget {
Json admit_observation(const Json &input, const std::string &digest, const Json &media) {
    using namespace scientific_detail;
    const auto proposal = Json::parse(input.dump());
    members(proposal, {"schema", "observations", "limitations"});
    require(proposal.at("schema").text() == "gadget.model.observation.v1",
            "unsupported observation schema");
    const auto observations = proposal.at("observations");
    require(json_object_is_type(observations.raw(), json_type_array) && observations.size() > 0 &&
                observations.size() <= 16,
            "observation count budget");
    (void)text(proposal.at("limitations"), 1024);
    auto claims = Json::array();
    for (std::size_t i = 0; i < observations.size(); ++i) {
        auto claim = Json::parse(observations.at(i).dump());
        members(claim, {"statement", "image_index"});
        require(!text(claim.at("statement"), 1024).empty(), "empty observation statement");
        if (claim.at("image_index").raw()) {
            require(media.size() > 0, "visual claim requires an attached image");
            const auto index = count(claim.at("image_index"), media.size() - 1, "image index");
            claim.set("image_sha256", media.at(index).at("sha256"));
        }
        claim.set("status", Json::string("prospective_unverified_model_statement"));
        claims.push(claim);
    }
    auto cert = Json::object();
    cert.set("schema", Json::string("gadget.model.observation.receipt.v1"));
    cert.set("source_digest", Json::string(digest));
    cert.set("proposal", proposal);
    cert.set("claims", claims);
    cert.set("media", Json::parse(media.dump()));
    cert.set(
        "scope",
        Json::string(
            "Input identities and output schema admitted; claim content is not verified and never "
            "becomes a native measurement, relation, identity map, time or clinical conclusion."));
    cert = seal(cert, "certificate_sha256");
    auto response = Json::object();
    response.set("schema", Json::string("gadget.agent.response.v1"));
    response.set("status", Json::string("prospective"));
    response.set("proposal", proposal);
    response.set("certificate", cert);
    return response;
}
} // namespace gadget
