#include "research_temporal.hpp"
#include "research_factored.hpp"
#include "research_joint_endpoints.hpp"
#include "research_prepared.hpp"
#include "scientific_relational_field.hpp"
#include <algorithm>
#include <gadget/research.hpp>
#include <rfc/temporal/moment_change.hpp>
#include <rfc/temporal/multitime.hpp>
#include <rfc/temporal/operation_word.hpp>
#include <rfc/temporal/source_delta.hpp>
#include <rfc/temporal/transfer.hpp>

namespace gadget::research_detail {
using namespace rfc;
namespace {
using namespace scientific_detail;
using Matrix = rfc::numeric::ExactMatrix;
Json rat(Rational x) {
    x.canonicalize();
    return Json::string(x.get_str());
}
Json vec(const Vector &x) {
    auto out = Json::array();
    for (const auto &v : x)
        out.push(rat(v));
    return out;
}
Json mat(const Matrix &x) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(x.rows()));
    out.set("columns", cardinality(x.cols()));
    for (const auto &[p, v] : x.entries()) {
        auto e = Json::object();
        e.set("row", cardinality(p.first));
        e.set("column", cardinality(p.second));
        e.set("value", rat(v));
        entries.push(e);
    }
    out.set("entries", entries);
    return out;
}
struct Word {
    rfc::temporal::OperationWord operations;
    std::vector<std::uint16_t> grades;
    std::vector<std::function<Vector(const Vector &)>> actions;
    std::vector<Json> factors;
    Vector apply_range(Vector x, std::size_t begin, std::size_t end) const {
        for (std::size_t i = begin; i < end; ++i)
            x = actions.at(i)(x);
        return x;
    }
    Vector apply(const Vector &x) const {
        return actions.empty() ? operations.apply(x) : apply_range(x, 0, actions.size());
    }
};
Word word(const TemporalReading &reading, const Json &declarations, bool validate,
          FactoredTemporalMaps *registry, std::size_t endpoint) {
    const auto &study = *reading.study;
    const auto initial = reading.grade;
    std::vector<Matrix> ops;
    std::vector<std::uint16_t> grades{initial.value};
    auto grade = initial;
    Word result;
    source_array(declarations, "temporal operation word");
    require(declarations.size(), "operation word requires declared factors");
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        auto d = declarations.at(i);
        members(d, {"kind"}, {"lambda", "fields"});
        const auto kind = d.at("kind").text();
        Matrix op;
        const auto input_grade = grade;
        Json factor;
        std::function<Vector(const Vector &)> action;
        if (kind == "boundary") {
            require(grade.value > 0, "word boundary at grade zero");
            if (!validate) {
                if (registry) {
                    factor = registry->boundary(endpoint, grade);
                    action = [&study, input_grade](const Vector &x) {
                        return study.source().boundary_apply(input_grade, x);
                    };
                } else
                    op = field::materialize_boundary_matrix(study.source(), grade);
            }
            --grade.value;
        } else if (kind == "adjoint") {
            require(grade.value < study.source().max_grade().value,
                    "word adjoint has no upper grade");
            auto next = core::Grade{static_cast<std::uint16_t>(grade.value + 1)};
            if (!validate) {
                if (registry) {
                    factor = registry->adjoint(endpoint, next);
                    action = [&study, input_grade, next](const Vector &x) {
                        return field::BoundaryOperator(study.source(), next)
                            .adjoint(x, study.tower().metrics()[input_grade.value],
                                     study.tower().metrics()[next.value]);
                    };
                } else
                    op = field::metric_adjoint_matrix(
                        field::materialize_boundary_matrix(study.source(), next),
                        study.tower().metrics()[next.value], study.tower().metrics()[grade.value]);
            }
            grade = next;
        } else if (kind == "diagonal_field") {
            require(d.has("fields") && reading.diagonal_fields.contains(i),
                    "word diagonal field is not admitted");
            if (!validate) {
                op = Matrix::diagonal(reading.diagonal_fields.at(i));
                if (registry) {
                    factor = registry->matrix(op);
                    action = [op](const Vector &x) { return op.apply(x); };
                }
            }
        } else if (kind == "hodge") {
            if (!validate) {
                if (registry) {
                    factor = registry->hodge(endpoint, grade);
                    action = [&study, input_grade](const Vector &x) {
                        return study.hodge(input_grade).apply(x);
                    };
                } else
                    op = reading.prepared->hodge(grade);
            }
        } else if (kind == "green") {
            require(d.has("lambda"), "word Green parameter missing");
            auto lambda = exact_rational(d.at("lambda").text());
            require(lambda > 0, "word Green parameter must be positive");
            if (!validate) {
                if (registry) {
                    factor = registry->green(endpoint, grade, lambda);
                    auto *prepared = reading.prepared;
                    action = [prepared, input_grade, lambda](const Vector &x) {
                        return prepared->green(input_grade, lambda).apply(x);
                    };
                } else
                    op = reading.prepared->green(grade, lambda).matrix();
            }
        } else
            require(false, "unsupported temporal word factor");
        require(kind == "green" || !d.has("lambda"), "lambda belongs to Green factors");
        require(kind == "diagonal_field" || !d.has("fields"), "fields belong to diagonal factors");
        if (!validate && registry) {
            result.actions.push_back(std::move(action));
            result.factors.push_back(factor);
        } else
            ops.push_back(std::move(op));
        grades.push_back(grade.value);
    }
    if (!validate && !registry)
        result.operations = temporal::OperationWord(std::move(ops));
    result.grades = std::move(grades);
    return result;
}
} // namespace
Json temporal_field_history(const Json &query, const std::vector<TemporalReading> &readings,
                            const std::vector<TemporalMaps> &maps, bool validate) {
    members(query, {"operation", "endpoints", "correspondences", "word", "lambda"},
            {"lineages", "map_output"});
    require(readings.size() >= 2 && maps.size() + 1 == readings.size(),
            "history requires consecutive endpoints and correspondences");
    const auto lambda = exact_rational(query.at("lambda").text());
    require(lambda > 0, "history Green parameter must be positive");
    const bool factored = FactoredTemporalMaps::requested(query);
    std::vector<PreparedTemporalStudy *> owners;
    if (!validate)
        for (const auto &r : readings)
            owners.push_back(r.prepared);
    FactoredTemporalMaps registry(std::move(owners));
    std::vector<Word> words;
    for (const auto &r : readings) {
        require(validate || (r.prepared && &r.prepared->study() == r.study),
                "temporal preparation must belong to the admitted endpoint");
        require(core::bianchi_through(r.study->source(), r.study->source().max_grade()).satisfied(),
                "temporal Hodge sectors require source closure");
        if (!words.empty() && r.shared_identity)
            (void)temporal::RelationalSourceDelta(readings[words.size() - 1].study->source(),
                                                  r.study->source());
        words.push_back(word(r, query.at("word"), validate,
                             factored && !validate ? &registry : nullptr, words.size()));
    }
    if (validate)
        return Json::object();
    auto out = Json::object(), endpoints = Json::array(), steps = Json::array();
    std::vector<Matrix> js;
    std::vector<Vector> innovations;
    for (std::size_t t = 0; t < readings.size(); ++t) {
        const auto &r = readings[t];
        auto e = Json::object();
        e.set("declaration", r.declaration);
        e.set("coordinates", r.coordinates);
        e.set("field", vec(r.values));
        e.set("word_field", vec(words[t].apply(r.values)));
        e.set("grade", cardinality(r.grade.value));
        endpoints.push(e);
        if (!t)
            continue;
        const auto &old = readings[t - 1];
        const auto &j = maps[t - 1].at(r.grade.value);
        js.push_back(j);
        innovations.push_back(temporal::field_delta(j, old.values, r.values).innovation);
        std::vector<Matrix> word_maps;
        for (auto k : words[t].grades)
            word_maps.push_back(maps[t - 1].at(k));
        Json direct_map;
        std::vector<Json> position_maps, contribution_maps;
        std::vector<Vector> contribution_fields;
        if (factored) {
            std::vector<Json> declared;
            for (const auto &map : word_maps)
                declared.push_back(registry.matrix(map));
            const auto &before = words[t - 1], &after = words[t];
            direct_map =
                registry.defect(registry.compose(before.factors), registry.compose(after.factors),
                                declared.front(), declared.back());
            Vector prefix_field = old.values;
            for (std::size_t i = 0; i < before.factors.size(); ++i) {
                const auto position = registry.defect(before.factors[i], after.factors[i],
                                                      declared[i], declared[i + 1]);
                std::vector<Json> order(before.factors.begin(), before.factors.begin() + i);
                order.push_back(position);
                order.insert(order.end(), after.factors.begin() + i + 1, after.factors.end());
                position_maps.push_back(position);
                contribution_maps.push_back(registry.compose(order));
                const auto old_next = before.actions[i](prefix_field);
                const auto position_field =
                    numeric::sub(after.actions[i](word_maps[i].apply(prefix_field)),
                                 word_maps[i + 1].apply(old_next));
                contribution_fields.push_back(
                    after.apply_range(position_field, i + 1, after.factors.size()));
                prefix_field = old_next;
            }
        } else {
            const auto defect = temporal::finite_word_defect(words[t - 1].operations,
                                                             words[t].operations, word_maps);
            require(defect.exact_identity(), "native ordered-word identity");
            direct_map = mat(defect.direct);
            for (std::size_t i = 0; i < defect.contributions.size(); ++i) {
                position_maps.push_back(mat(defect.position_defects[i]));
                contribution_maps.push_back(mat(defect.contributions[i]));
                contribution_fields.push_back(defect.contributions[i].apply(old.values));
            }
        }
        auto step = Json::object(), factors = Json::array();
        step.set("from", cardinality(t - 1));
        step.set("to", cardinality(t));
        step.set("source_maps", Json::array());
        for (const auto &map : maps[t - 1])
            step.at("source_maps").push(mat(map));
        step.set("structure", compare_prepared_structures(*old.prepared, *r.prepared, r.grade, {},
                                                          {}, old.values, r.values, lambda,
                                                          &maps[t - 1], true, factored));
        step.set("source_identity_comparable", Json::boolean(r.shared_identity));
        if (r.declaration.has("lineage"))
            step.set("lineage", r.declaration.at("lineage"));
        auto source_changes = Json::array();
        if (r.shared_identity) {
            const temporal::RelationalSourceDelta source_delta(old.study->source(),
                                                               r.study->source());
            for (const auto &g : source_delta.grades()) {
                if (g.grade == r.grade) {
                    auto observations = Json::array();
                    const auto before = old.declaration.at("observations");
                    const auto after = r.declaration.at("observations");
                    for (std::size_t i = 0; i < g.old_to_new.size(); ++i) {
                        const auto target = g.old_to_new[i];
                        if (target == temporal::GradeSourceDelta::absent_index)
                            continue;
                        auto observation = Json::object();
                        observation.set("old_coordinate", cardinality(i));
                        observation.set("new_coordinate", cardinality(target));
                        observation.set("old", before.at(i));
                        observation.set("new", after.at(target));
                        observations.push(observation);
                    }
                    step.set("retained_observation_alignment", observations);
                }

                auto record = Json::object(), born = Json::array(), removed = Json::array(),
                     born_occ = Json::array(), removed_occ = Json::array();
                record.set("grade", cardinality(g.grade.value));
                for (const auto &c : g.born_cells)
                    born.push(Json::string(std::to_string(c.value)));
                for (const auto &c : g.removed_cells)
                    removed.push(Json::string(std::to_string(c.value)));
                for (const auto &c : g.born_occurrences)
                    born_occ.push(Json::string(std::to_string(c.value)));
                for (const auto &c : g.removed_occurrences)
                    removed_occ.push(Json::string(std::to_string(c.value)));
                record.set("born_cell_keys", born);
                record.set("removed_cell_keys", removed);
                record.set("born_occurrence_keys", born_occ);
                record.set("removed_occurrence_keys", removed_occ);
                record.set("stable_identity_correspondence",
                           mat(source_delta.correspondence(g.grade)));
                record.set("structure_unchanged", Json::boolean(g.source_structure_unchanged()));
                auto occurrence_changes = Json::array();
                for (const auto &change : g.occurrence_changes) {
                    auto c = Json::object();
                    c.set("key", Json::string(std::to_string(change.occurrence.value)));
                    c.set("old_lower_key", Json::string(std::to_string(change.old_lower.value)));
                    c.set("new_lower_key", Json::string(std::to_string(change.new_lower.value)));
                    c.set("old_parent_key", Json::string(std::to_string(change.old_parent.value)));
                    c.set("new_parent_key", Json::string(std::to_string(change.new_parent.value)));
                    c.set("old_coefficient", rat(change.old_coefficient));
                    c.set("new_coefficient", rat(change.new_coefficient));
                    c.set("attachment_changed", Json::boolean(change.attachment_changed()));
                    c.set("existence_term", rat(change.composite.existence));
                    c.set("orientation_term", rat(change.composite.orientation));
                    c.set("magnitude_term", rat(change.composite.magnitude));
                    c.set("coefficient_change", rat(change.composite.aggregate));
                    auto aspects = Json::array();
                    for (auto aspect : change.changed_aspects)
                        aspects.push(cardinality(aspect));
                    c.set("changed_aspects", aspects);
                    occurrence_changes.push(c);
                }
                auto type_changes = Json::array(), order_changes = Json::array();
                for (const auto &change : g.relation_type_changes) {
                    auto c = Json::object();
                    c.set("key", Json::string(std::to_string(change.relation.value)));
                    c.set("old_type", cardinality(change.old_type));
                    c.set("new_type", cardinality(change.new_type));
                    type_changes.push(c);
                }
                for (const auto &change : g.attaching_order_changes) {
                    auto c = Json::object(), old_order = Json::array(), new_order = Json::array();
                    c.set("key", Json::string(std::to_string(change.relation.value)));
                    for (auto k : change.old_order)
                        old_order.push(Json::string(std::to_string(k.value)));
                    for (auto k : change.new_order)
                        new_order.push(Json::string(std::to_string(k.value)));
                    c.set("old_order", old_order);
                    c.set("new_order", new_order);
                    order_changes.push(c);
                }
                record.set("occurrence_changes", occurrence_changes);
                record.set("relation_type_changes", type_changes);
                record.set("attaching_order_changes", order_changes);
                source_changes.push(record);
            }
        }
        step.set("source_changes", source_changes);
        step.set("word_defect", direct_map);
        std::vector<temporal::MultiTimeChannel> delta_channels;
        for (std::size_t i = 0; i < contribution_fields.size(); ++i) {
            auto f = Json::object();
            const auto &v = contribution_fields[i];
            f.set("position", cardinality(i));
            f.set("defect", position_maps[i]);
            f.set("contribution", contribution_maps[i]);
            f.set("field", vec(v));
            factors.push(f);
            delta_channels.push_back({"factor_" + std::to_string(i), v});
        }
        delta_channels.push_back({"field_innovation", words[t].apply(innovations.back())});
        const auto &metric = r.study->tower().metrics()[words[t].grades.back()];
        const auto gram = temporal::multi_time_channel_gram(delta_channels, metric);
        Vector delta(metric.size(), Rational{0});
        for (const auto &c : delta_channels)
            delta = numeric::add(delta, c.target_field);
        const auto direct = numeric::sub(words[t].apply(r.values),
                                         word_maps.back().apply(words[t - 1].apply(old.values)));
        require(delta == direct, "word field-delta reconstruction");
        Rational gram_sum{0};
        for (const auto &[p, v] : gram.entries()) {
            (void)p;
            gram_sum += v;
        }
        require(gram_sum == metric.inner(delta, delta), "complete field-delta Gram contraction");
        step.set("word_factors", factors);
        step.set("field_delta", vec(delta));
        step.set("field_delta_channel_gram", mat(gram));
        auto retained_channels = Json::array();
        for (const auto &channel : delta_channels) {
            auto e = Json::object();
            e.set("label", Json::string(channel.label));
            e.set("field", vec(channel.target_field));
            retained_channels.push(e);
        }
        step.set("field_delta_channels", retained_channels);
        step.set("field_delta_quadrance", rat(metric.inner(delta, delta)));
        const auto old_word_field = words[t - 1].apply(old.values),
                   new_word_field = words[t].apply(r.values);
        // Contract the two complete cause fields, rather than constructing
        // W* M W on every input coordinate. These are coefficient channels
        // (transported field and innovation), not a replacement data carrier.
        const auto transported_word_field = words[t].apply(j.apply(old.values));
        const auto innovation_word_field = words[t].apply(innovations.back());
        Matrix old_form(2, 2);
        old_form.set(0, 0,
                     old.study->tower().metrics()[words[t - 1].grades.back()].inner(
                         old_word_field, old_word_field));
        const auto new_form =
            temporal::multi_time_channel_gram({{"transported_field", transported_word_field},
                                               {"field_innovation", innovation_word_field}},
                                              metric);
        const auto channel_identity = Matrix::identity(2);
        const auto change = temporal::changing_dimension_cross_moment(
            {Rational{1}, Rational{0}}, {Rational{1}, Rational{0}}, old_form,
            {Rational{1}, Rational{1}}, {Rational{1}, Rational{1}}, new_form, channel_identity,
            channel_identity);
        require(change.total == metric.inner(new_word_field, new_word_field) -
                                    old.study->tower().metrics()[words[t - 1].grades.back()].inner(
                                        old_word_field, old_word_field),
                "native endpoint quadrance-change identity");
        auto terms = Json::object();
        terms.set("pullback_form", rat(change.pullback_form));
        terms.set("left_innovation", rat(change.left_innovation));
        terms.set("right_innovation", rat(change.right_innovation));
        terms.set("innovation_interaction", rat(change.innovation_interaction));
        terms.set("total", rat(change.total));
        const Rational interval = exact_decimal(r.declaration.at("time").at("position").text()) -
                                  exact_decimal(old.declaration.at("time").at("position").text());
        terms.set("finite_rate", rat(temporal::finite_rate(change.total, interval)));
        terms.set("interval", rat(interval));
        terms.set("time_unit", r.declaration.at("time").at("unit"));
        step.set("endpoint_quadrance_change", terms);
        steps.push(step);
    }
    const auto accumulation = temporal::accumulate_innovations(readings[0].values, js, innovations);
    require(accumulation.terminal == readings.back().values,
            "native accumulated innovation identity");
    std::vector<temporal::MultiTimeChannel> channels{{"initial", accumulation.transported_initial}};
    for (std::size_t t = 0; t < accumulation.propagated_innovations.size(); ++t)
        channels.push_back(
            {"innovation_" + std::to_string(t), accumulation.propagated_innovations[t]});
    auto fields = Json::array();
    for (const auto &c : channels) {
        auto e = Json::object();
        e.set("label", Json::string(c.label));
        e.set("field", vec(c.target_field));
        fields.push(e);
    }
    const auto &terminal = readings.back();
    const auto &metric = terminal.study->tower().metrics()[terminal.grade.value];
    const auto gram = temporal::multi_time_channel_gram(channels, metric);
    Rational sum{0};
    for (const auto &[p, v] : gram.entries()) {
        (void)p;
        sum += v;
    }
    require(sum == metric.inner(terminal.values, terminal.values),
            "complete terminal Gram contraction");
    out.set("endpoints", endpoints);
    out.set("steps", steps);
    out.set("terminal_channels", fields);
    out.set("terminal_channel_gram", mat(gram));
    out.set("terminal_quadrance", rat(sum));
    if (factored)
        out.set("factored_maps", registry.archive());
    out.set("identity_verified", Json::boolean(true));
    return out;
}
Json censored_endpoint_bounds(const NativeStudy &study, const Json &query,
                              const std::map<std::string, Vector> &fields,
                              const std::vector<bool> &upper_support,
                              const std::map<std::string, Json> &meanings, const Json &coordinates,
                              bool validate) {
    members(query,
            {"operation", "context", "lower", "upper", "lower_closed", "upper_closed",
             "upper_finite", "selection", "strata", "times", "origin", "time_unit", "population",
             "event", "evidence"},
            {"contrasts", "joint_family"});
    const auto lower_id = query.at("lower").text(), upper_id = query.at("upper").text();
    const auto &l = fields.at(lower_id), &u = fields.at(upper_id);
    const auto n = l.size();
    const auto origin = exact_rational(query.at("origin").text());
    const auto unit = query.at("time_unit").text();
    require(!unit.empty() && !query.at("event").text().empty() &&
                !query.at("population").text().empty(),
            "endpoint time, event and population meanings required");
    for (const auto *key : {"channel", "quantity", "method"})
        require(meanings.at(lower_id).at(key).text() == meanings.at(upper_id).at(key).text(),
                "lower and upper event-time meanings differ");
    for (const auto &id : {lower_id, upper_id})
        require(meanings.at(id).at("value_type").text() == "quantitative" &&
                    meanings.at(id).at("units").text() == unit,
                "censored time field units/type mismatch");
    std::vector<std::string> indicator_ids{query.at("lower_closed").text(),
                                           query.at("upper_closed").text(),
                                           query.at("upper_finite").text()};
    for (const auto *key : {"selection", "strata"}) {
        source_array(query.at(key), key);
        std::set<std::string> seen;
        for (std::size_t i = 0; i < query.at(key).size(); ++i) {
            auto id = query.at(key).at(i).text();
            require(seen.insert(id).second, "duplicate endpoint indicator");
            indicator_ids.push_back(id);
        }
    }
    for (const auto &id : indicator_ids)
        require(meanings.at(id).at("value_type").text() == "indicator",
                "censoring and grouping require explicit indicator fields");
    const auto &lc = fields.at(query.at("lower_closed").text()),
               &uc = fields.at(query.at("upper_closed").text()),
               &uf = fields.at(query.at("upper_finite").text());
    std::vector<Rational> times;
    source_array(query.at("times"), "endpoint evaluation times");
    require(query.at("times").size(), "endpoint evaluation times empty");
    for (std::size_t t = 0; t < query.at("times").size(); ++t) {
        auto s = exact_rational(query.at("times").at(t).text());
        require(s >= origin && (times.empty() || times.back() < s),
                "evaluation times must increase from origin");
        times.push_back(s);
    }
    Vector selected(n, Rational{1});
    std::size_t population = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t f = 0; f < query.at("selection").size(); ++f)
            selected[i] *= fields.at(query.at("selection").at(f).text())[i];
        if (selected[i] == 0)
            continue;
        ++population;
        require(l[i] >= origin, "event lower bound precedes origin");
        require(uf[i] == 0 || upper_support[i], "finite event upper bound is unobserved");
        require(uf[i] != 0 || uc[i] == 0, "infinite upper endpoint cannot be closed");
        if (uf[i] != 0)
            require(u[i] > l[i] || (u[i] == l[i] && lc[i] == 1 && uc[i] == 1),
                    "empty event-time admissible interval");
    }
    require(population > 0, "endpoint selection empty");

    Matrix readout(1 + query.at("strata").size(), n);
    auto groups = Json::array();
    for (std::size_t g = 0; g < readout.rows(); ++g) {
        Vector membership = selected;
        if (g)
            for (std::size_t i = 0; i < n; ++i)
                membership[i] *= fields.at(query.at("strata").at(g - 1).text())[i];
        Rational mass{0};
        for (const auto &v : membership)
            mass += v;
        auto group = Json::object();
        group.set("id", g ? query.at("strata").at(g - 1) : Json::string("selected_population"));
        group.set("size", rat(mass));
        group.set("defined", Json::boolean(mass != 0));
        groups.push(group);
        if (mass != 0)
            for (std::size_t i = 0; i < n; ++i)
                if (membership[i] != 0)
                    readout.set(g, i, membership[i] / mass);
    }
    auto contrasts = Json::array();
    const std::size_t contrasts_count = query.has("contrasts") ? query.at("contrasts").size() : 0;
    Matrix contrast_readout(contrasts_count, n);
    std::set<std::string> contrast_ids;
    if (query.has("contrasts"))
        source_array(query.at("contrasts"), "endpoint contrasts");
    for (std::size_t c = 0; c < contrasts_count; ++c) {
        auto d = query.at("contrasts").at(c);
        members(d, {"id", "left", "right", "meaning"});
        require(contrast_ids.insert(d.at("id").text()).second && !d.at("meaning").text().empty(),
                "endpoint contrast identity/meaning");
        std::optional<std::size_t> left, right;
        for (std::size_t g = 0; g < groups.size(); ++g) {
            if (groups.at(g).at("id").text() == d.at("left").text())
                left = g;
            if (groups.at(g).at("id").text() == d.at("right").text())
                right = g;
        }
        require(left && right, "endpoint contrast group unknown");
        auto record = Json::parse(d.dump());
        const bool defined =
            groups.at(*left).at("defined").boolean() && groups.at(*right).at("defined").boolean();
        record.set("defined", Json::boolean(defined));
        contrasts.push(record);
        if (defined)
            for (std::size_t i = 0; i < n; ++i)
                contrast_readout.set(c, i, readout.get(*left, i) - readout.get(*right, i));
    }
    Matrix positive(contrasts_count, n), negative(contrasts_count, n);
    for (const auto &[p, v] : contrast_readout.entries())
        (v > 0 ? positive : negative).set(p.first, p.second, v);
    require(core::bianchi_through(study.source(), study.source().max_grade()).satisfied(),
            "censored relational field requires source closure");
    Json joint;
    if (query.has("joint_family"))
        joint =
            joint_endpoint_family(study, query.at("joint_family"), l, u, lc, uc, uf, selected,
                                  times, origin, readout, contrast_readout, coordinates, validate);
    if (validate)
        return Json::object();
    auto observations = Json::array();
    for (std::size_t i = 0; i < n; ++i)
        if (selected[i] != 0) {
            auto r = Json::object();
            r.set("coordinate", coordinates.at(i));
            r.set("lower", rat(l[i]));
            r.set("upper", uf[i] != 0 ? rat(u[i]) : Json{});
            r.set("lower_closed", Json::boolean(lc[i] == 1));
            r.set("upper_closed", Json::boolean(uc[i] == 1));
            r.set("first_possible_grid_index", cardinality(times.size()));
            r.set("first_definite_grid_index", cardinality(times.size()));
            for (std::size_t t = 0; t < times.size(); ++t) {
                if (l[i] < times[t] || (l[i] == times[t] && lc[i] == 1)) {
                    r.set("first_possible_grid_index", cardinality(t));
                    break;
                }
            }
            for (std::size_t t = 0; t < times.size(); ++t)
                if (uf[i] != 0 && u[i] <= times[t]) {
                    r.set("first_definite_grid_index", cardinality(t));
                    break;
                }
            observations.push(r);
        }
    const PreparedRelationalField hodge(study, core::Grade{0});
    auto readings = Json::array();
    for (const auto &s : times) {
        Vector emin(n, Rational{0}), emax(n, Rational{0}), rmin(n, Rational{0}),
            rmax(n, Rational{0});
        auto lower_attained = Json::array(), upper_attained = Json::array();
        for (std::size_t i = 0; i < n; ++i) {
            if (selected[i] != 0) {
                emin[i] = uf[i] != 0 && u[i] <= s ? 1 : 0;
                emax[i] = l[i] < s || (l[i] == s && lc[i] == 1) ? 1 : 0;
                rmin[i] = std::min(l[i], s) - origin;
                rmax[i] = (uf[i] != 0 ? std::min(u[i], s) : s) - origin;
            }
            lower_attained.push(Json::boolean(selected[i] == 0 || lc[i] == 1 || l[i] >= s));
            upper_attained.push(
                Json::boolean(selected[i] == 0 || uf[i] == 0 || uc[i] == 1 || u[i] > s));
        }
        auto r = Json::object();
        r.set("time", rat(s));
        r.set("event_fraction_contrast_lower",
              vec(numeric::add(positive.apply(emin), negative.apply(emax))));
        r.set("event_fraction_contrast_upper",
              vec(numeric::add(positive.apply(emax), negative.apply(emin))));
        r.set("restricted_time_contrast_infimum",
              vec(numeric::add(positive.apply(rmin), negative.apply(rmax))));
        r.set("restricted_time_contrast_supremum",
              vec(numeric::add(positive.apply(rmax), negative.apply(rmin))));
        r.set("event_field_lower", vec(emin));
        r.set("event_field_upper", vec(emax));
        r.set("cumulative_event_fraction_lower", vec(readout.apply(emin)));
        r.set("cumulative_event_fraction_upper", vec(readout.apply(emax)));
        r.set("restricted_event_free_time_infimum", vec(readout.apply(rmin)));
        r.set("restricted_event_free_time_supremum", vec(readout.apply(rmax)));
        r.set("restricted_time_lower_field", vec(rmin));
        r.set("restricted_time_upper_field", vec(rmax));
        r.set("lower_coordinate_bound_attained", lower_attained);
        r.set("upper_coordinate_bound_attained", upper_attained);
        std::vector<std::size_t> ambiguous;
        for (std::size_t i = 0; i < n; ++i)
            if (emin[i] != emax[i])
                ambiguous.push_back(i);
        Matrix frame(n, ambiguous.size());
        for (std::size_t j = 0; j < ambiguous.size(); ++j)
            frame.set(ambiguous[j], j, Rational{1});
        auto family = Json::object();
        family.set("offset", vec(emin));
        family.set("frame", mat(frame));
        family.set(
            "parameter_domain",
            Json::string("independent zero-or-one event indicators at this evaluation time"));
        auto sectors = Json::object();
        const auto offsets = hodge.apply(emin);
        constexpr const char *names[] = {"full", "gradient", "curl", "harmonic"};
        std::array<Matrix, 4> sector_frames{
            Matrix(n, ambiguous.size()), Matrix(n, ambiguous.size()), Matrix(n, ambiguous.size()),
            Matrix(n, ambiguous.size())};
        for (std::size_t j = 0; j < ambiguous.size(); ++j) {
            Vector e(n, Rational{0});
            e[ambiguous[j]] = 1;
            const auto parts = hodge.apply(e);
            for (std::size_t k = 0; k < 4; ++k)
                for (std::size_t i = 0; i < n; ++i)
                    if (parts[k][i] != 0)
                        sector_frames[k].set(i, j, parts[k][i]);
        }
        for (std::size_t k = 0; k < 4; ++k) {
            auto a = Json::object();
            a.set("offset", vec(offsets[k]));
            a.set("frame", mat(sector_frames[k]));
            Vector lower = offsets[k], upper = offsets[k];
            for (const auto &[p, v] : sector_frames[k].entries()) {
                if (v < 0)
                    lower[p.first] += v;
                else
                    upper[p.first] += v;
            }
            a.set("coordinate_lower", vec(lower));
            a.set("coordinate_upper", vec(upper));
            sectors.set(names[k], a);
        }
        family.set("sectors", sectors);
        r.set("admissible_event_family", family);

        readings.push(r);
    }
    auto out = Json::object();
    out.set("population", query.at("population"));
    out.set("event", query.at("event"));
    out.set("origin", rat(origin));
    out.set("time_unit", query.at("time_unit"));
    out.set("groups", groups);
    out.set("contrasts", contrasts);
    out.set("contrast_readout", mat(contrast_readout));
    out.set("readout", mat(readout));
    out.set("observations", observations);
    out.set("readings", readings);
    out.set("coordinate_carrier", coordinates);
    auto reader = Json::object(), aggregate_readings = Json::array();
    reader.set(
        "scope",
        Json::string("Finite-population marginal event and restricted-time bounds; original "
                     "coordinates and coupled field families remain in the source certificate"));
    reader.set("groups", groups);
    reader.set("contrasts", contrasts);
    reader.set("time_unit", query.at("time_unit"));
    for (std::size_t i = 0; i < readings.size(); ++i) {
        auto e = Json::object();
        e.set("full_result_index", cardinality(i));
        for (const auto *key :
             {"time", "cumulative_event_fraction_lower", "cumulative_event_fraction_upper",
              "restricted_event_free_time_infimum", "restricted_event_free_time_supremum",
              "event_fraction_contrast_lower", "event_fraction_contrast_upper"})
            e.set(key, readings.at(i).at(key));
        aggregate_readings.push(e);
    }
    reader.set("readings", aggregate_readings);
    out.set("reader", reader);

    out.set("metric", mat(study.tower().metrics()[0].materialize()));
    out.set("sharp_for",
            Json::string("Each displayed marginal over the independent declared event-time "
                         "intervals; bounds at different times or for different groups need not "
                         "share one attaining event-time assignment."));
    std::size_t unidentified = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (selected[i] != 0 && (uf[i] == 0 || u[i] != l[i]))
            ++unidentified;
    out.set("unidentified_event_count", cardinality(unidentified));
    auto facts = Json::object();
    facts.set("unidentified_event_count", cardinality(unidentified));
    out.set("facts", facts);
    if (query.has("joint_family"))
        out.set("joint_family", joint);
    return out;
}
} // namespace gadget::research_detail
