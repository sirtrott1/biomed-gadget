#include "research_joint_endpoints.hpp"
#include "scientific_ensemble.hpp"
#include "scientific_relational_field.hpp"
#include <rfc/local/constraint_system.hpp>
#include <rfc/temporal/multitime.hpp>

namespace gadget::research_detail {
using namespace rfc;
using namespace scientific_detail;
namespace {
using Matrix = numeric::ExactMatrix;
struct State {
    std::string id;
    Rational weight;
    std::vector<std::optional<Rational>> events;
    local::FiniteConstraintSystem<std::string>::Assignment assignment;
};
} // namespace
scientific_detail::ResultValue joint_endpoint_family(const NativeStudy &study, const Json &declaration, const Vector &lower,
                           const Vector &upper, const Vector &lower_closed,
                           const Vector &upper_closed, const Vector &upper_finite,
                           const Vector &selected, const std::vector<Rational> &times,
                           const Rational &origin, const Matrix &readout, const Matrix &contrasts,
                           const Json &coordinates, bool validate, CensoredResultFields *native) {
    members(declaration, {"meaning", "evidence", "states"}, {"constraints"});
    (void)text(declaration.at("meaning"));
    const auto n = lower.size();
    std::map<std::string, std::size_t> entities;
    for (std::size_t i = 0; i < n; ++i)
        if (selected[i] != 0)
            entities.emplace(coordinates.at(i).at("id").text(), i);
    auto event_value = [&](const Json &value, std::size_t i) -> std::optional<Rational> {
        if (json_object_is_type(value.raw(), json_type_null)) {
            require(upper_finite[i] == 0, "no-event state conflicts with finite event interval");
            return {};
        }
        auto t = exact_rational(value.text());
        require(
            (t > lower[i] || (t == lower[i] && lower_closed[i] == 1)) &&
                (upper_finite[i] == 0 || t < upper[i] || (t == upper[i] && upper_closed[i] == 1)),
            "joint state lies outside the recorded event interval");
        return t;
    };
    const auto token = [](const std::optional<Rational> &value) {
        return value ? value->get_str() : std::string("no_event");
    };
    source_array(declaration.at("states"), "joint event states");
    require(declaration.at("states").size(), "joint event family is empty");
    std::vector<State> states;
    Rational mass{0};
    std::set<std::string> ids;
    std::vector<std::vector<std::string>> domains(n);
    for (std::size_t s = 0; s < declaration.at("states").size(); ++s) {
        const auto d = declaration.at("states").at(s);
        members(d, {"id", "weight", "values"});
        State state{text(d.at("id")),
                    exact_rational(d.at("weight").text()),
                    std::vector<std::optional<Rational>>(n),
                    {}};
        require(ids.insert(state.id).second && state.weight > 0,
                "joint state identity/positive law mass required");
        mass += state.weight;
        source_array(d.at("values"), "joint event assignment");
        require(d.at("values").size() == entities.size(),
                "joint state must retain every selected original coordinate");
        std::set<std::size_t> seen;
        for (std::size_t v = 0; v < d.at("values").size(); ++v) {
            const auto entry = d.at("values").at(v);
            members(entry, {"entity", "event_time"});
            require(entities.contains(entry.at("entity").text()), "unknown joint event entity");
            const auto i = entities.at(entry.at("entity").text());
            require(seen.insert(i).second, "duplicate joint event entity");
            state.events[i] = event_value(entry.at("event_time"), i);
            const auto value = token(state.events[i]);
            domains[i].push_back(value);
            state.assignment.emplace(core::CellId{core::Grade{0}, static_cast<std::uint32_t>(i)},
                                     value);
        }
        states.push_back(std::move(state));
    }
    require(mass == 1, "joint observation law must sum exactly to one");
    local::FiniteConstraintSystem<std::string> constraints(study.source());
    for (const auto &[entity, i] : entities) {
        (void)entity;
        constraints.set_domain({core::Grade{0}, static_cast<std::uint32_t>(i)}, domains[i]);
    }
    std::set<std::string> chart_ids;
    if (declaration.has("constraints")) {
        source_array(declaration.at("constraints"), "joint event constraints");
        for (std::size_t c = 0; c < declaration.at("constraints").size(); ++c) {
            const auto d = declaration.at("constraints").at(c);
            members(d, {"id", "entities", "allowed", "meaning", "evidence"});
            require(chart_ids.insert(text(d.at("id"))).second, "duplicate joint constraint");
            (void)text(d.at("meaning"));
            source_array(d.at("entities"), "joint constraint entities");
            source_array(d.at("allowed"), "joint constraint tuples");
            algebra::SubcomplexMask patch(study.source());
            std::vector<core::CellId> support;
            for (std::size_t j = 0; j < d.at("entities").size(); ++j) {
                const auto entity = d.at("entities").at(j).text();
                require(entities.contains(entity), "joint constraint entity outside selection");
                core::CellId cell{core::Grade{0}, static_cast<std::uint32_t>(entities.at(entity))};
                support.push_back(cell);
                patch.select(cell);
            }
            std::vector<std::vector<std::string>> allowed;
            for (std::size_t a = 0; a < d.at("allowed").size(); ++a) {
                const auto tuple = d.at("allowed").at(a);
                source_array(tuple, "joint constraint tuple");
                require(tuple.size() == support.size(), "joint constraint tuple arity");
                std::vector<std::string> values;
                for (std::size_t j = 0; j < tuple.size(); ++j)
                    values.push_back(token(event_value(tuple.at(j), support[j].index)));
                allowed.push_back(std::move(values));
            }
            constraints.add_chart(d.at("id").text(), std::move(patch), std::move(support),
                                  std::move(allowed));
        }
    }
    for (const auto &state : states)
        require(constraints.globally_valid(state.assignment),
                "joint observation law includes a state violating its native constraints");
    if (validate)
        return ResultValue::object();
    require(native, "joint endpoint native result owner required");
    const auto source_axis = native->source, group_axis = native->groups, contrast_axis = native->contrasts;
    const auto sample_axis = native->publication.axis(states.size());
    const auto rat = [&](const Rational &value) { return native->publication.stage_scalar(value); };
    const auto vec = [&](const Vector &value, std::uint64_t axis) { return ResultFields::stage_vector(value, axis); };
    const auto mat = [&](const Matrix &value, std::uint64_t input, std::uint64_t output) {
        return ResultFields::stage_matrix(value, input, output);
    };
    const PreparedRelationalField hodge(study, core::Grade{0});
    const auto &metric = study.tower().metrics()[0];
    std::vector<temporal::MultiTimeChannel> event_means, restricted_means;
    const auto channel_count = machine::checked_size_mul(times.size(), std::size_t{4});
    Matrix event_moment(channel_count, channel_count), restricted_moment = event_moment;
    constexpr const char *names[] = {"full", "gradient", "curl", "harmonic"};
    for (std::size_t t = 0; t < times.size(); ++t)
        for (const auto *name : names) {
            const auto label = std::to_string(t) + ":" + name;
            event_means.push_back({label, Vector(n, Rational{0})});
            restricted_means.push_back({label, Vector(n, Rational{0})});
        }
    const auto channel_axis = native->publication.axis(channel_count);
    auto law_binding = ResultValue::object(), state_ids = ResultValue::array(), channel_labels = ResultValue::array();
    law_binding.set("kind", ResultValue::string("finite_joint_endpoint_law"));
    law_binding.set("source_carrier", ResultValue::count(source_axis));
    law_binding.set("state_carrier", ResultValue::count(sample_axis));
    law_binding.set("channel_carrier", ResultValue::count(channel_axis));
    for (const auto &state : states) state_ids.push(ResultValue::string(state.id));
    for (const auto &channel : event_means) channel_labels.push(ResultValue::string(channel.label));
    law_binding.set("states", std::move(state_ids));
    law_binding.set("channels", std::move(channel_labels));
    auto outputs = ResultValue::array();
    Vector weights;
    std::vector<std::vector<Vector>> event_samples(channel_count), duration_samples(channel_count);
    for (const auto &state : states) {
        weights.push_back(state.weight);
        auto record = ResultValue::object(), readings = ResultValue::array();
        record.set("id", ResultValue::string(state.id));
        record.set("weight", rat(state.weight));
        std::vector<std::size_t> selected_rows, event_rows;
        for (std::size_t i = 0; i < n; ++i) {
            if (selected[i] == 0) continue;
            selected_rows.push_back(i);
            if (state.events[i]) event_rows.push_back(i);
        }
        record.set("event_times", ResultFields::staged(assay_detail::SourceGradeFields::tensor(n,
            [&](std::size_t i) {return state.events[i].value_or(Rational{0});},
            numeric::StructuralSupport::sparse(n, event_rows)), {}, source_axis, "vector"));
        record.set("event_presence", ResultFields::staged(assay_detail::SourceGradeFields::tensor(n,
            [&](std::size_t i) {return Rational{state.events[i] ? 1 : 0};},
            numeric::StructuralSupport::sparse(n, selected_rows)), {}, source_axis, "vector"));
        std::vector<temporal::MultiTimeChannel> events, durations;
        for (std::size_t t = 0; t < times.size(); ++t) {
            Vector event(n, Rational{0}), duration = event;
            for (const auto &[entity, i] : entities) {
                (void)entity;
                event[i] = state.events[i] && *state.events[i] <= times[t] ? 1 : 0;
                duration[i] =
                    (state.events[i] ? std::min(*state.events[i], times[t]) : times[t]) - origin;
            }
            auto reading = ResultValue::object(), event_sectors = ResultValue::object(),
                 duration_sectors = ResultValue::object();
            const auto e = hodge.apply(event), r = hodge.apply(duration);
            for (std::size_t k = 0; k < e.size(); ++k) {
                const auto index = events.size();
                event_sectors.set(names[k], vec(e[k], source_axis));
                duration_sectors.set(names[k], vec(r[k], source_axis));
                events.push_back({event_means[index].label, e[k]});
                durations.push_back({restricted_means[index].label, r[k]});
                event_samples[index].push_back(e[k]);
                duration_samples[index].push_back(r[k]);
                event_means[index].target_field = numeric::add(event_means[index].target_field,
                                                               numeric::scale(state.weight, e[k]));
                restricted_means[index].target_field = numeric::add(
                    restricted_means[index].target_field, numeric::scale(state.weight, r[k]));
            }
            reading.set("time", rat(times[t]));
            reading.set("event_sectors", event_sectors);
            reading.set("restricted_time_sectors", duration_sectors);
            reading.set("event_fraction", vec(readout.apply(event), group_axis));
            reading.set("restricted_time", vec(readout.apply(duration), group_axis));
            reading.set("event_fraction_contrast", vec(contrasts.apply(event), contrast_axis));
            reading.set("restricted_time_contrast", vec(contrasts.apply(duration), contrast_axis));
            readings.push(reading);
        }
        event_moment =
            event_moment.add(temporal::multi_time_channel_gram(events, metric).scale(state.weight));
        restricted_moment = restricted_moment.add(
            temporal::multi_time_channel_gram(durations, metric).scale(state.weight));
        record.set("readings", readings);
        outputs.push(record);
    }
    const auto moment = [&](const std::vector<temporal::MultiTimeChannel> &means,
                            const Matrix &second, const std::vector<std::vector<Vector>> &samples) {
        auto out = ResultValue::object(), channels = ResultValue::array();
        for (const auto &mean : means) {
            auto channel = ResultValue::object();
            channel.set("label", ResultValue::string(mean.label));
            channel.set("mean_field", vec(mean.target_field, source_axis));
            channels.push(channel);
        }
        const auto mean_gram = temporal::multi_time_channel_gram(means, metric);
        out.set("channels", channels);
        out.set("expected_gram", mat(second, channel_axis, channel_axis));
        out.set("mean_gram", mat(mean_gram, channel_axis, channel_axis));
        out.set("centered_gram", mat(second.add(mean_gram.scale(Rational{-1})), channel_axis, channel_axis));
        std::vector<WeightedFieldEnsemble> ensembles;
        auto families = ResultValue::array();
        for (const auto &fields : samples) {
            ensembles.emplace_back(fields, weights, metric);
            const auto &ensemble = ensembles.back();
            auto family = ResultValue::object();
            family.set("centered_fields", mat(ensemble.centered, source_axis, sample_axis));
            family.set("rank", ResultValue::count(ensemble.rank));
            family.set("variance",
                       rat(ensemble.moments ? ensemble.moments->variance : Rational{0}));
            families.push(family);
        }
        auto coupling = ResultValue::array();
        for (std::size_t a = 0; a < ensembles.size(); ++a)
            for (std::size_t b = a; b < ensembles.size(); ++b) {
                auto pair = ResultValue::object();
                const bool defined = ensembles[a].moments && ensembles[b].moments;
                pair.set("left", ResultValue::count(a));
                pair.set("right", ResultValue::count(b));
                pair.set("defined", ResultValue::boolean(defined));
                pair.set("value", defined ? rat(weighted_field_coupling(ensembles[a], ensembles[b],
                                                                        weights, metric))
                                          : ResultValue{});
                coupling.push(pair);
            }
        out.set("field_families", families);
        out.set("ensemble_coupling", coupling);
        return out;
    };
    law_binding.set("weights", vec(weights, sample_axis));
    native->publication.bind(std::move(law_binding));
    auto out = ResultValue::object();
    out.set("declaration", ResultValue(declaration));
    out.set("scope", ResultValue::string("The supplied finite joint observation law; recorded event "
                                  "intervals remain unchanged"));
    out.set("native_constraints_checked", ResultValue::count(chart_ids.size()));
    out.set("states", outputs);
    out.set("event_moments", moment(event_means, event_moment, event_samples));
    out.set("restricted_time_moments",
            moment(restricted_means, restricted_moment, duration_samples));
    return out;
}
} // namespace gadget::research_detail
