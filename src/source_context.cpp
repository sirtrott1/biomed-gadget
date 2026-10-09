#include "scientific_detail.hpp"
#include <gadget/assay.hpp>

namespace gadget {
using namespace scientific_detail;

Json ScientificStudy::source_target(const Json &input, const Json &evidence) const {
    const auto target = Json::parse(input.dump());
    const auto kind = target.at("kind").text();
    members(target, {"kind", "source_digest"},
            kind == "occurrence" ? std::set<std::string>{"grade", "key"}
                                 : std::set<std::string>{"id"});
    require(target.at("source_digest").text() == state_->digest,
            "source target belongs to a different native study");
    auto descriptor = Json::object();
    auto declaration = [&](const Json &value) {
        json_object_object_foreach(value.raw(), name, member) {
            (void)member;
            if (std::string(name) != "evidence")
                descriptor.set(name, value.at(name));
        }
    };
    if (kind == "cell") {
        const auto id = text(target.at("id"), source_byte_limit());
        require(state_->has_cell(id), "source cell is undeclared");
        declaration(state_->cell_record(state_->resolve_cell(id)));
    } else if (kind == "occurrence") {
        const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(
            count(target.at("grade"), state_->grade_pins.size() - 1, "source occurrence grade"))};
        const rfc::storage::StableOccurrenceKey stable{
            grade, rfc::core::OccurrenceKey{key(target.at("key"))}};
        if(!state_->compact)require(state_->occurrences.contains(stable), "source occurrence is undeclared");
        declaration(state_->occurrence_record(stable));
        // Canonical pinned access retains zero coefficients and inactive flags.
        const auto occurrence = state_->directory->occurrence(stable);
        descriptor.set("upper", state_->cell_record(occurrence.parent()).at("id"));
        descriptor.set("lower", state_->cell_record(occurrence.participant()).at("id"));
        descriptor.set("coefficient", Json::string(occurrence.coefficient().get_str()));
        auto aspects = Json::object();
        const auto endpoint = state_->directory->pin_snapshot();
        const auto &schema = endpoint->source.signature().aspects();
        for (std::size_t a = 0; a < schema.size(); ++a) {
            const auto &aspect = schema.at(static_cast<rfc::composite::AspectId>(a));
            aspects.set(aspect.name, Json::boolean(occurrence.aspect(aspect.name)));
        }
        descriptor.set("aspects", aspects);
    } else if (kind == "field") {
        const auto id = text(target.at("id"), source_byte_limit());
        require(state_->fields.contains(id), "source field is undeclared");
        const auto &field = state_->fields.at(id);
        declaration(field.declaration);
        const auto &pin = state_->field_pins.at(field.tensor);
        const auto tensor = pin.tensor();
        const auto carrier = pin.carrier();
        descriptor.set("grade", cardinality(carrier->grade.value));
        descriptor.set("carrier", Json::string(std::to_string(carrier->key)));
        descriptor.set("coordinate_system",
                       Json::string(std::to_string(carrier->coordinate_system)));
        auto extents = Json::array();
        for (const auto extent : tensor->domain().extents())
            extents.push(cardinality(extent));
        descriptor.set("extents", extents);
        descriptor.set("coordinates", cardinality(tensor->domain().logical_size()));
        descriptor.set("reported", cardinality(tensor->support().count()));
    } else {
        require(kind == "context", "unknown source target kind");
        const auto id = text(target.at("id"), source_byte_limit());
        require(state_->contexts.contains(id), "source context is undeclared");
        declaration(state_->contexts.at(id).declaration);
    }
    auto out = Json::object();
    out.set("target", target);
    out.set("descriptor", descriptor);
    out.set("source_digest", Json::string(state_->digest));
    out.set("native_sha256", Json::string(state_->native_digest));
    out.set("source_namespace", state_->declaration.at("native").at("source_namespace"));
    out.set("evidence", source_evidence(evidence, state_->documents, true));
    return seal(out, "target_sha256");
}

Json admit_source_annotations(const AssayWorkspace &workspace, const Json &request,
                              const Json &input, const std::string &question) {
    const auto context = workspace.source_context(request);
    require(question == context.at("question").text(), "annotation question binding mismatch");
    const auto proposal = Json::parse(input.dump());
    members(proposal, {"schema", "action", "source_digest", "context_sha256", "question_sha256"},
            {"mentions", "reason"});
    require(proposal.at("schema").text() == "gadget.assays.annotation.proposal.v1" &&
                proposal.at("source_digest").text() == workspace.digest() &&
                proposal.at("context_sha256").text() == context.at("context_sha256").text() &&
                proposal.at("question_sha256").text() == context.at("question_sha256").text(),
            "annotation proposal source/context/question/schema mismatch");
    const auto action = proposal.at("action").text();
    auto mentions = Json::array();
    if (action == "clarify") {
        require(proposal.has("reason") && !proposal.has("mentions"),
                "annotation clarification requires only a reason");
        (void)text(proposal.at("reason"), source_byte_limit());
    } else {
        require(action == "annotate" && proposal.has("mentions") && !proposal.has("reason"),
                "source annotation requires mentions");
        const auto candidates = proposal.at("mentions");
        source_array(candidates, "annotation mentions");
        require(candidates.size() > 0, "source annotation requires mentions");
        std::map<std::string, Json> targets;
        for (std::size_t i = 0; i < context.at("targets").size(); ++i) {
            const auto entry = context.at("targets").at(i);
            targets.emplace(entry.at("id").text(), entry);
        }
        std::set<std::string> ids;
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            auto mention = Json::parse(candidates.at(i).dump());
            members(mention, {"id", "target", "interpretation", "polarity"});
            require(ids.insert(text(mention.at("id"), source_byte_limit())).second,
                    "duplicate annotation mention identity");
            const auto target = text(mention.at("target"), source_byte_limit());
            require(targets.contains(target), "annotation target is undeclared");
            (void)text(mention.at("interpretation"), source_byte_limit());
            const auto polarity = mention.at("polarity").text();
            require(polarity == "asserted" || polarity == "negated" || polarity == "uncertain",
                    "annotation polarity must be explicit");
            mention.set("binding_sha256", targets.at(target).at("binding_sha256"));
            mention.set("interpretation_verified", Json::boolean(false));
            mentions.push(mention);
        }
    }
    auto receipt = Json::object();
    receipt.set("schema", Json::string("gadget.assays.annotation.receipt.v1"));
    receipt.set("status", Json::string(action == "clarify" ? "clarification" : "prospective"));
    receipt.set("source_digest", Json::string(workspace.digest()));
    receipt.set("context_sha256", context.at("context_sha256"));
    receipt.set("question_sha256", context.at("question_sha256"));
    receipt.set("context", context);
    receipt.set("proposal", proposal);
    receipt.set("mentions", mentions);
    receipt = seal(receipt, "receipt_sha256");
    auto out = Json::object();
    out.set("schema", Json::string("gadget.agent.response.v1"));
    out.set("status", receipt.at("status"));
    out.set("proposal", proposal);
    out.set("certificate", receipt);
    return out;
}
} // namespace gadget
