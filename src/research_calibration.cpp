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

std::optional<NativeRecord> evaluate_finite_reference(const ScientificStudy &study, const Json &query, bool validate,
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
    if (validate) {
        (void)linear_model_comparison(study, model_query, true, cpu);
        return std::nullopt;
    }
    auto models = linear_model_fields(study, model_query, cpu);
    auto contraction_reservation = cpu->reserve(1);
    scientific_detail::ResultFields publication(models.publication);
    using Value = scientific_detail::ResultValue;
    const auto state_carrier = publication.axis(states.size()),
               hypothesis_carrier = publication.axis(hypotheses.size());
    const auto rational = [&](const Rational &value) {return publication.scalar(value);};
    auto binding = Value::object();
    binding.set("state_carrier", Value::count(state_carrier));
    binding.set("hypothesis_carrier", Value::count(hypothesis_carrier));
    binding.set("states", Value(declarations));
    binding.set("hypotheses", Value(family));
    binding.set("masses", publication.vector(weights, state_carrier));
    binding.set("model_bindings", Value(models.publication.at("bindings").unbound_value()));
    std::vector<Vector> scores(states.size(), Vector(hypotheses.size(), Rational{0}));
    std::vector<std::vector<bool>> defined(states.size(),
                                           std::vector<bool>(hypotheses.size(), true));
    for (std::size_t s = 0; s < states.size(); ++s) {
        for (std::size_t h = 0; h < hypotheses.size(); ++h) {
            const auto &hypothesis = hypotheses[h];
            const auto found = models.comparisons.find(LinearModelAxis{
                hypothesis.stratum, hypothesis.reduced, states[s].responses[hypothesis.response_index]});
            if (found == models.comparisons.end() || !found->second) defined[s][h] = false;
            else {
                scores[s][h] = *found->second;
                require(scores[s][h] >= 0 && scores[s][h] <= 1,
                        "nested dimensionless reference statistic");
            }
        }
    }
    const auto size = rfc::machine::checked_size_mul(states.size(), hypotheses.size());
    std::vector<std::size_t> support;
    std::vector<rfc::numeric::SparseRationalEntry> entries;
    for (std::size_t s = 0; s < states.size(); ++s)
        for (std::size_t h = 0; h < hypotheses.size(); ++h) if (defined[s][h]) {
            const auto index = s * hypotheses.size() + h;
            support.push_back(index);
            if (scores[s][h] != 0) entries.push_back({index, scores[s][h]});
        }
    const auto score_key = publication.observation(rfc::numeric::SparseRationalTensor(
        rfc::numeric::IndexDomain(std::vector<std::size_t>{states.size(), hypotheses.size()}),
        rfc::numeric::StructuralSupport::sparse(size, std::move(support)),
        rfc::numeric::SparseRationalVector(size, std::move(entries))), hypothesis_carrier, state_carrier);
    binding.set("statistics", Value(score_key));
    publication.bind(std::move(binding));
    auto state_records = Value::array();
    for (std::size_t s = 0; s < states.size(); ++s) {
        auto record = Value::object(), values = Value::array();
        record.set("id", Value::string(states[s].id));
        record.set("mass", rational(states[s].mass));
        record.set("responses", declarations.at(s).at("responses"));
        for (std::size_t h = 0; h < hypotheses.size(); ++h)
            values.push(publication.element(score_key, s * hypotheses.size() + h));
        record.set("statistics", values);
        state_records.push(record);
    }
    auto calibration =
        calibrate_finite_scores_fields(scores, defined, weights, *observed_state, family, alpha, publication, state_carrier, hypothesis_carrier);
    auto result = Value::object();
    result.set("schema", Value::string("gadget.finite.reference.calibration.v1"));
    result.set("status", calibration.at("status"));
    result.set("arithmetic", Value::string("exact_rational"));
    result.set("law", law);
    result.set("states", state_records);
    result.set("hypotheses", calibration.at("hypotheses"));
    result.set("expanded_model_query", model_query);
    result.set("model_comparison", Value(models.publication.at("view").unbound_value()));
    result.set("total_mass", rational(total));
    if (alpha)
        result.set("alpha", rational(*alpha));
    result.set("joint_family", calibration.at("joint_family"));
    result.set("tail_rule", Value::string("Sum declared masses where statistic >= observed; include "
                                         "all ties and the observed state"));
    result.set("calibration_scope",
               Value::string(
                   "Conditional on the complete declared joint finite law. Marginal tails and the "
                   "family maximum are superuniform under that law. Source checks do not establish "
                   "exchangeability, population validity or strong control under partial nulls."));
    return publication.finish(std::move(result));
}
NativeRecord finite_reference_fields(const ScientificStudy &study, const Json &query,
                                      const std::shared_ptr<CpuResources> &resources) {
    return evaluate_finite_reference(study, query, false, resources).value();
}
Json finite_reference_calibration(const ScientificStudy &study, const Json &query, bool validate,
                                 const std::shared_ptr<CpuResources> &resources) {
    auto result = evaluate_finite_reference(study, query, validate, resources);
    return result ? render_record_view(*result) : Json{};
}
} // namespace gadget::research_detail
