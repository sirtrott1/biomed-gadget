#include "research_calibration.hpp"
#include "research_fields.hpp"
#include "research_linear.hpp"
#include "scientific_reference.hpp"
#include <algorithm>
#include <gadget/resources.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/statistics/ensemble.hpp>
#include <tuple>

namespace gadget::research_detail {
namespace {
Json rational(Rational value) {
    value.canonicalize();
    return Json::string(value.get_str());
}
struct State {
    std::string id;
    Rational mass;
    std::vector<std::string> responses;
};
struct Hypothesis {
    std::string id, stratum, response, reduced;
    std::size_t response_index;
};
} // namespace

std::set<std::string> finite_reference_field_ids(const Json &query) {
    auto ids = linear_model_field_ids(query.at("model"));
    const auto states = query.at("law").at("states");
    source_array(states, "finite reference states");
    for (std::size_t s = 0; s < states.size(); ++s)
        for (const auto &id :
             names(states.at(s).at("responses"), states.at(s).at("responses").size()))
            ids.insert(id);
    return ids;
}

Json finite_reference_calibration(const ScientificStudy &study, const Json &query, bool validate,
                                  const std::shared_ptr<CpuResources> &resources) {
    members(query, {"operation", "model", "law", "hypotheses", "workers"}, {"alpha"});
    auto model_query = Json::parse(query.at("model").dump_compact());
    require(model_query.at("operation").text() == "linear_model_comparison",
            "finite calibration needs a declared linear model comparison");
    if (model_query.has("workers"))
        require(model_query.at("workers").integer() == query.at("workers").integer(),
                "finite calibration and model must share one worker budget");
    model_query.set("workers", query.at("workers"));
    const auto observed = names(model_query.at("responses"), model_query.at("responses").size());
    const auto law = query.at("law");
    members(law, {"kind", "observed_state", "states", "null_hypothesis", "population",
                  "assumptions", "derivation", "evidence"});
    require(law.at("kind").text() == "declared_joint_finite", "finite reference law kind");
    for (const auto &key_ : {"null_hypothesis", "population", "derivation"})
        (void)text(law.at(key_), source_byte_limit());
    const auto assumptions = names(law.at("assumptions"), law.at("assumptions").size());
    (void)assumptions;
    source_array(law.at("evidence"), "finite law source evidence");
    require(law.at("evidence").size() > 0, "finite reference law requires source evidence");
    const auto declarations = law.at("states");
    source_array(declarations, "finite reference states");
    require(declarations.size() > 0, "finite reference law must have positive support");
    std::vector<State> states;
    std::set<std::string> state_ids, all_responses(observed.begin(), observed.end());
    Vector weights;
    Rational total{0};
    std::optional<std::size_t> observed_state;
    const auto observed_id = text(law.at("observed_state"), source_byte_limit());
    for (std::size_t s = 0; s < declarations.size(); ++s) {
        const auto declaration = declarations.at(s);
        members(declaration, {"id", "mass", "responses"});
        State state{text(declaration.at("id"), source_byte_limit()),
                    exact_rational(text(declaration.at("mass"), source_byte_limit())),
                    names(declaration.at("responses"), observed.size())};
        require(state_ids.insert(state.id).second && state.mass > 0 &&
                    state.responses.size() == observed.size(),
                "finite reference state identity, positive mass and response axes");
        if (state.id == observed_id)
            observed_state = s;
        all_responses.insert(state.responses.begin(), state.responses.end());
        total += state.mass;
        weights.push_back(state.mass);
        states.push_back(std::move(state));
    }
    require(total == 1 && observed_state.has_value(),
            "finite reference law must sum exactly to one and contain its observed state");
    const Inputs inputs(study, finite_reference_field_ids(query), false, std::nullopt);
    for (const auto &state : states)
        for (std::size_t h = 0; h < observed.size(); ++h) {
            const auto &origin = observed[h], &reference = state.responses[h];
            for (const auto &key_ : {"units", "quantity", "value_type"})
                require(inputs.meaning.at(origin).at(key_).text() ==
                            inputs.meaning.at(reference).at(key_).text(),
                        "finite reference response meaning differs from observed axis");
            for (std::size_t i = 0; i < inputs.entities.size(); ++i) {
                const bool reported = inputs.fields.at(origin)->tensor.support().contains(i);
                require(reported == inputs.fields.at(reference)->tensor.support().contains(i),
                        "finite reference response must retain observed support");
                if (state.id == observed_id && reported)
                    require(inputs.value(origin, i) == inputs.value(reference, i),
                            "finite law observed state differs from original observations");
            }
        }
    std::set<std::string> stratum_ids{"all"}, reduced_ids{"controls"};
    if (model_query.has("strata"))
        for (std::size_t s = 0; s < model_query.at("strata").size(); ++s)
            stratum_ids.insert(text(model_query.at("strata").at(s).at("id"), source_byte_limit()));
    if (model_query.has("ablations"))
        for (std::size_t s = 0; s < model_query.at("ablations").size(); ++s)
            reduced_ids.insert(
                text(model_query.at("ablations").at(s).at("id"), source_byte_limit()));
    const auto family = query.at("hypotheses");
    source_array(family, "finite calibration hypothesis family");
    require(family.size() > 0, "finite calibration requires a declared hypothesis family");
    std::vector<Hypothesis> hypotheses;
    std::set<std::string> hypothesis_ids;
    std::set<std::tuple<std::string, std::string, std::string>> hypothesis_axes;
    for (std::size_t h = 0; h < family.size(); ++h) {
        const auto declared = family.at(h);
        members(declared, {"id", "stratum", "response", "reduced", "statistic"});
        Hypothesis hypothesis{text(declared.at("id"), source_byte_limit()),
                              text(declared.at("stratum"), source_byte_limit()),
                              text(declared.at("response"), source_byte_limit()),
                              text(declared.at("reduced"), source_byte_limit()), 0};
        const auto found = std::find(observed.begin(), observed.end(), hypothesis.response);
        require(
            found != observed.end() && stratum_ids.contains(hypothesis.stratum) &&
                reduced_ids.contains(hypothesis.reduced) &&
                declared.at("statistic").text() == "fraction_of_reduced_residual" &&
                hypothesis_ids.insert(hypothesis.id).second &&
                hypothesis_axes.emplace(hypothesis.stratum, hypothesis.response, hypothesis.reduced)
                    .second,
            "finite hypothesis must identify unique admitted model axes and statistic");
        hypothesis.response_index = static_cast<std::size_t>(found - observed.begin());
        hypotheses.push_back(std::move(hypothesis));
    }
    std::optional<Rational> alpha;
    if (query.has("alpha")) {
        alpha = exact_rational(text(query.at("alpha"), source_byte_limit()));
        require(*alpha > 0 && *alpha < 1,
                "finite calibration alpha must lie strictly between zero and one");
    }
    auto joint_responses = Json::array();
    for (const auto &id : all_responses)
        joint_responses.push(Json::string(id));
    model_query.set("responses", joint_responses);
    auto cpu = resources ? resources : CpuResources::local();
    auto models = linear_model_comparison(study, model_query, validate, cpu);
    if (validate)
        return Json{};
    auto contraction_reservation = cpu->reserve(1);
    std::map<std::tuple<std::string, std::string, std::string>, Json> readings;
    for (std::size_t s = 0; s < models.at("strata").size(); ++s) {
        const auto stratum = models.at("strata").at(s);
        const auto id = stratum.at("id").text();
        if (!stratum.has("comparisons"))
            continue;
        for (std::size_t c = 0; c < stratum.at("comparisons").size(); ++c) {
            const auto comparison = stratum.at("comparisons").at(c);
            for (std::size_t r = 0; r < comparison.at("responses").size(); ++r) {
                const auto reading = comparison.at("responses").at(r);
                readings.emplace(std::make_tuple(id, comparison.at("reduced").text(),
                                                 reading.at("field").text()),
                                 reading);
            }
        }
    }
    std::vector<Vector> scores(states.size(), Vector(hypotheses.size(), Rational{0}));
    std::vector<std::vector<bool>> defined(states.size(),
                                           std::vector<bool>(hypotheses.size(), true));
    auto state_records = Json::array();
    for (std::size_t s = 0; s < states.size(); ++s) {
        auto record = Json::object(), values = Json::array();
        record.set("id", Json::string(states[s].id));
        record.set("mass", rational(states[s].mass));
        record.set("responses", declarations.at(s).at("responses"));
        for (std::size_t h = 0; h < hypotheses.size(); ++h) {
            const auto &hypothesis = hypotheses[h];
            const auto found =
                readings.find(std::make_tuple(hypothesis.stratum, hypothesis.reduced,
                                              states[s].responses[hypothesis.response_index]));
            if (found == readings.end() ||
                !found->second.at("fraction_of_reduced_residual").raw()) {
                values.push(Json{});
                defined[s][h] = false;
            } else {
                scores[s][h] =
                    exact_rational(found->second.at("fraction_of_reduced_residual").text());
                require(scores[s][h] >= 0 && scores[s][h] <= 1,
                        "nested dimensionless reference statistic");
                values.push(rational(scores[s][h]));
            }
        }
        record.set("statistics", values);
        state_records.push(record);
    }
    auto calibration =
        calibrate_finite_scores(scores, defined, weights, *observed_state, family, alpha);
    auto result = Json::object();
    result.set("schema", Json::string("gadget.finite.reference.calibration.v1"));
    result.set("status", calibration.at("status"));
    result.set("arithmetic", Json::string("exact_rational"));
    result.set("law", law);
    result.set("states", state_records);
    result.set("hypotheses", calibration.at("hypotheses"));
    result.set("expanded_model_query", model_query);
    result.set("model_comparison", models);
    result.set("total_mass", rational(total));
    if (alpha)
        result.set("alpha", rational(*alpha));
    result.set("joint_family", calibration.at("joint_family"));
    result.set("tail_rule", Json::string("Sum declared masses where statistic >= observed; include "
                                         "all ties and the observed state"));
    result.set("calibration_scope",
               Json::string(
                   "Conditional on the complete declared joint finite law. Marginal tails and the "
                   "family maximum are superuniform under that law. Source checks do not establish "
                   "exchangeability, population validity or strong control under partial nulls."));
    return result;
}
} // namespace gadget::research_detail
