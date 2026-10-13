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
template<class Value> struct BasicWord {
    rfc::temporal::OperationWord operations;
    std::vector<std::uint16_t> grades;
    std::vector<std::function<Vector(const Vector &)>> actions;
    std::vector<Value> factors;
    std::optional<rfc::numeric::ExactActionLimits> limits;
    std::function<void()> checkpoint;
    Vector apply_range(Vector x, std::size_t begin, std::size_t end) const {
        for (std::size_t i = begin; i < end; ++i)
            x = actions.at(i)(x);
        return x;
    }
    Vector apply(const Vector &x) const {
        if(!actions.empty())return apply_range(x,0,actions.size());
        if(!limits)return operations.apply(x);
        rfc::numeric::ExactActionGuard guard(*limits,checkpoint);auto value=x;
        for(std::size_t i=0;i<operations.length();++i)value=rfc::numeric::checked_apply(operations.op(i),value,guard);
        return value;
    }
};
template<class Value> BasicWord<Value> word(const TemporalReading &reading, const Json &declarations, bool validate,
          BasicFactoredTemporalMaps<Value> *registry, std::size_t endpoint,
          TemporalResultFields *native) {
    const auto &study = *reading.study;
    const auto initial = reading.grade;
    std::vector<Matrix> ops;
    std::vector<std::uint16_t> grades{initial.value};
    auto grade = initial;
    BasicWord<Value> result;
    if(native) {
        std::size_t dimension=0;
        for(std::size_t k=0;k<study.source().grade_count();++k)dimension=std::max(dimension,study.source().grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count());
        const auto &policy=native->lease.policy().limits;
        result.limits=numeric::ExactActionLimits{dimension,policy.sparse_entries,policy.scalar_work,policy.integer_bits};
        result.checkpoint=[lease=native->lease] {lease.check_cancelled();};
    }
    source_array(declarations, "temporal operation word");
    require(declarations.size(), "operation word requires declared factors");
    for (std::size_t i = 0; i < declarations.size(); ++i) {
        auto d = declarations.at(i);
        members(d, {"kind"}, {"lambda", "fields"});
        const auto kind = d.at("kind").text();
        Matrix op;
        const auto input_grade = grade;
        Value factor;
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
                    factor = registry->matrix(op, native ? std::optional<std::uint64_t>(native->axes.at(endpoint).at(input_grade.value)) : std::nullopt,
                        native ? std::optional<std::uint64_t>(native->axes.at(endpoint).at(input_grade.value)) : std::nullopt);
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
            if (native) {
                auto original_action=std::move(action);const auto lease=native->lease;
                action=[original_action=std::move(original_action),lease](const Vector &x) {
                    lease.check_cancelled();auto value=original_action(x);lease.check_cancelled();return value;
                };
            }
            result.actions.push_back(std::move(action));
            result.factors.push_back(factor);
        } else {
            if (native && !validate) admit_structure_matrix(op,native->lease);
            ops.push_back(std::move(op));
        }
        grades.push_back(grade.value);
    }
    if (!validate && !registry)
        result.operations = temporal::OperationWord(std::move(ops));
    result.grades = std::move(grades);
    return result;
}
} // namespace
template<class Value> Value collect_temporal_history(const Json &query, const std::vector<TemporalReading> &readings,
                            const std::vector<TemporalMaps> &maps, bool validate, TemporalResultFields *native) {
    const auto axis=[&](std::size_t endpoint,std::uint16_t grade) {return native ? native->axes.at(endpoint).at(grade) : 0;};
    const auto scalar=[&](const Rational &value) -> Value {
        if constexpr (std::is_same_v<Value,ResultValue>) {
            StructuralResultFields admission{native->publication,{},{},native->lease};admission.admit(value);
            return native->publication.stage_scalar(value);
        } else return rat(value);
    };
    const auto vector=[&](const Vector &value,std::uint64_t carrier) -> Value {
        if constexpr (std::is_same_v<Value,ResultValue>) {
            StructuralResultFields admission{native->publication,{},{},native->lease};admission.admit(value);
            return ResultFields::stage_vector(value,carrier);
        } else { (void)carrier;return vec(value); }
    };
    const auto matrix=[&](const Matrix &value,std::uint64_t input,std::uint64_t output) -> Value {
        if constexpr (std::is_same_v<Value,ResultValue>) {
            admit_structure_matrix(value,native->lease);auto field=ResultFields::stage_matrix(value,input,output);
            field.set("display",ResultValue::string("sparse_matrix"));return field;
        } else { (void)input;(void)output;return mat(value); }
    };
    const auto count=[](std::size_t value) -> Value {
        if constexpr (std::is_same_v<Value,ResultValue>) return ResultValue::count(value);
        else return cardinality(value);
    };
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
    BasicFactoredTemporalMaps<Value> registry(std::move(owners),native ? &native->publication : nullptr,
        native ? native->axes : std::vector<std::vector<std::uint64_t>>{},
        [native](const Matrix &m) {if(native)admit_structure_matrix(m,native->lease);},
        [native](const Rational &v) {if(native) {StructuralResultFields admission{native->publication,{},{},native->lease};admission.admit(v);}});
    std::vector<BasicWord<Value>> words;
    for (const auto &r : readings) {
        require(validate || (r.prepared && &r.prepared->study() == r.study),
                "temporal preparation must belong to the admitted endpoint");
        require(core::bianchi_through(r.study->source(), r.study->source().max_grade()).satisfied(),
                "temporal Hodge sectors require source closure");
        if (!words.empty() && r.shared_identity)
            (void)temporal::RelationalSourceDelta(readings[words.size() - 1].study->source(),
                                                  r.study->source());
        words.push_back(word<Value>(r, query.at("word"), validate,
                             factored && !validate ? &registry : nullptr, words.size(),native));
    }
    if (validate)
        return Value::object();
    auto out = Value::object(), endpoints = Value::array(), steps = Value::array();
    std::vector<Matrix> js;
    std::vector<Vector> innovations;
    for (std::size_t t = 0; t < readings.size(); ++t) {
        const auto &r = readings[t];
        auto e = Value::object();
        if constexpr (std::is_same_v<Value,ResultValue>) {
            e.set("declaration",native->declarations.at(t));
            e.set("coordinates",native->coordinates.at(t));
        } else {
            e.set("declaration",r.declaration);e.set("coordinates",r.coordinates);
        }
        e.set("field", vector(r.values,axis(t,r.grade.value)));
        e.set("word_field", vector(words[t].apply(r.values),axis(t,words[t].grades.back())));
        e.set("grade", count(r.grade.value));
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
        Value direct_map;
        std::vector<Value> position_maps, contribution_maps;
        std::vector<Vector> contribution_fields;
        if (factored) {
            std::vector<Value> declared;
            for (std::size_t position=0;position<word_maps.size();++position)
                declared.push_back(registry.matrix(word_maps[position],axis(t-1,words[t].grades[position]),axis(t,words[t].grades[position])));
            const auto &before = words[t - 1], &after = words[t];
            direct_map =
                registry.defect(registry.compose(before.factors), registry.compose(after.factors),
                                declared.front(), declared.back());
            Vector prefix_field = old.values;
            for (std::size_t i = 0; i < before.factors.size(); ++i) {
                const auto position = registry.defect(before.factors[i], after.factors[i],
                                                      declared[i], declared[i + 1]);
                std::vector<Value> order(before.factors.begin(), before.factors.begin() + i);
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
            direct_map = matrix(defect.direct,axis(t-1,words[t-1].grades.front()),axis(t,words[t].grades.back()));
            for (std::size_t i = 0; i < defect.contributions.size(); ++i) {
                position_maps.push_back(matrix(defect.position_defects[i],axis(t-1,words[t-1].grades[i]),axis(t,words[t].grades[i+1])));
                contribution_maps.push_back(matrix(defect.contributions[i],axis(t-1,words[t-1].grades.front()),axis(t,words[t].grades.back())));
                contribution_fields.push_back(defect.contributions[i].apply(old.values));
            }
        }
        auto step = Value::object(), factors = Value::array();
        step.set("from", count(t - 1));
        step.set("to", count(t));
        step.set("source_maps", Value::array());
        for (std::size_t k=0;k<maps[t-1].size();++k) {
            auto source_maps=step.at("source_maps");source_maps.push(matrix(maps[t-1][k],axis(t-1,k),axis(t,k)));step.set("source_maps",std::move(source_maps));
        }
        if constexpr (std::is_same_v<Value,ResultValue>) {
            StructuralResultFields fields{native->publication,native->axes.at(t-1),native->axes.at(t),native->lease};
            step.set("structure",compare_prepared_structure_fields(*old.prepared,*r.prepared,r.grade,{}, {},
                old.values,r.values,lambda,&maps[t-1],factored,fields,true));
        } else step.set("structure",compare_prepared_structures(*old.prepared,*r.prepared,r.grade,{}, {},
            old.values,r.values,lambda,&maps[t-1],true,factored));
        step.set("source_identity_comparable", Value::boolean(r.shared_identity));
        if (r.declaration.has("lineage"))
            step.set("lineage", Value(r.declaration.at("lineage")));
        auto source_changes = Value::array();
        if (r.shared_identity) {
            const temporal::RelationalSourceDelta source_delta(old.study->source(),
                                                               r.study->source());
            for (const auto &g : source_delta.grades()) {
                if (g.grade == r.grade) {
                    auto observations = Value::array();
                    const auto before = old.declaration.at("observations");
                    const auto after = r.declaration.at("observations");
                    for (std::size_t i = 0; i < g.old_to_new.size(); ++i) {
                        const auto target = g.old_to_new[i];
                        if (target == temporal::GradeSourceDelta::absent_index)
                            continue;
                        auto observation = Value::object();
                        observation.set("old_coordinate", count(i));
                        observation.set("new_coordinate", count(target));
                        if constexpr (std::is_same_v<Value,ResultValue>) {
                            observation.set("old",native->declarations.at(t-1).at("observations").at(i));
                            observation.set("new",native->declarations.at(t).at("observations").at(target));
                        } else {
                            observation.set("old",before.at(i));observation.set("new",after.at(target));
                        }
                        observations.push(observation);
                    }
                    step.set("retained_observation_alignment", observations);
                }

                auto record = Value::object(), born = Value::array(), removed = Value::array(),
                     born_occ = Value::array(), removed_occ = Value::array();
                record.set("grade", count(g.grade.value));
                for (const auto &c : g.born_cells)
                    born.push(Value::string(std::to_string(c.value)));
                for (const auto &c : g.removed_cells)
                    removed.push(Value::string(std::to_string(c.value)));
                for (const auto &c : g.born_occurrences)
                    born_occ.push(Value::string(std::to_string(c.value)));
                for (const auto &c : g.removed_occurrences)
                    removed_occ.push(Value::string(std::to_string(c.value)));
                record.set("born_cell_keys", born);
                record.set("removed_cell_keys", removed);
                record.set("born_occurrence_keys", born_occ);
                record.set("removed_occurrence_keys", removed_occ);
                record.set("stable_identity_correspondence",
                           matrix(source_delta.correspondence(g.grade),axis(t-1,g.grade.value),axis(t,g.grade.value)));
                record.set("structure_unchanged", Value::boolean(g.source_structure_unchanged()));
                auto occurrence_changes = Value::array();
                for (const auto &change : g.occurrence_changes) {
                    auto c = Value::object();
                    c.set("key", Value::string(std::to_string(change.occurrence.value)));
                    c.set("old_lower_key", Value::string(std::to_string(change.old_lower.value)));
                    c.set("new_lower_key", Value::string(std::to_string(change.new_lower.value)));
                    c.set("old_parent_key", Value::string(std::to_string(change.old_parent.value)));
                    c.set("new_parent_key", Value::string(std::to_string(change.new_parent.value)));
                    c.set("old_coefficient", scalar(change.old_coefficient));
                    c.set("new_coefficient", scalar(change.new_coefficient));
                    c.set("attachment_changed", Value::boolean(change.attachment_changed()));
                    c.set("existence_term", scalar(change.composite.existence));
                    c.set("orientation_term", scalar(change.composite.orientation));
                    c.set("magnitude_term", scalar(change.composite.magnitude));
                    c.set("coefficient_change", scalar(change.composite.aggregate));
                    auto aspects = Value::array();
                    for (auto aspect : change.changed_aspects)
                        aspects.push(count(aspect));
                    c.set("changed_aspects", aspects);
                    occurrence_changes.push(c);
                }
                auto type_changes = Value::array(), order_changes = Value::array();
                for (const auto &change : g.relation_type_changes) {
                    auto c = Value::object();
                    c.set("key", Value::string(std::to_string(change.relation.value)));
                    c.set("old_type", count(change.old_type));
                    c.set("new_type", count(change.new_type));
                    type_changes.push(c);
                }
                for (const auto &change : g.attaching_order_changes) {
                    auto c = Value::object(), old_order = Value::array(), new_order = Value::array();
                    c.set("key", Value::string(std::to_string(change.relation.value)));
                    for (auto k : change.old_order)
                        old_order.push(Value::string(std::to_string(k.value)));
                    for (auto k : change.new_order)
                        new_order.push(Value::string(std::to_string(k.value)));
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
            auto f = Value::object();
            const auto &v = contribution_fields[i];
            f.set("position", count(i));
            f.set("defect", position_maps[i]);
            f.set("contribution", contribution_maps[i]);
            f.set("field", vector(v,axis(t,words[t].grades.back())));
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
        step.set("field_delta", vector(delta,axis(t,words[t].grades.back())));
        const auto channels_axis=native ? native->publication.axis(delta_channels.size()) : 0;
        step.set("field_delta_channel_gram",matrix(gram,channels_axis,channels_axis));
        auto retained_channels = Value::array();
        for (const auto &channel : delta_channels) {
            auto e = Value::object();
            e.set("label", Value::string(channel.label));
            e.set("field", vector(channel.target_field,axis(t,words[t].grades.back())));
            retained_channels.push(e);
        }
        step.set("field_delta_channels", retained_channels);
        step.set("field_delta_quadrance", scalar(metric.inner(delta, delta)));
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
        auto terms = Value::object();
        if constexpr (std::is_same_v<Value,ResultValue>) {
            const auto causes=native->publication.axis(2);
            terms.set("old_form",matrix(old_form,causes,causes));
            terms.set("new_form",matrix(new_form,causes,causes));
            terms.set("correspondence",matrix(channel_identity,causes,causes));
            terms.set("old_left",vector(Vector{Rational{1},Rational{0}},causes));
            terms.set("old_right",vector(Vector{Rational{1},Rational{0}},causes));
            terms.set("new_left",vector(Vector{Rational{1},Rational{1}},causes));
            terms.set("new_right",vector(Vector{Rational{1},Rational{1}},causes));
        }
        terms.set("pullback_form", scalar(change.pullback_form));
        terms.set("left_innovation", scalar(change.left_innovation));
        terms.set("right_innovation", scalar(change.right_innovation));
        terms.set("innovation_interaction", scalar(change.innovation_interaction));
        terms.set("total", scalar(change.total));
        const Rational interval = exact_decimal(r.declaration.at("time").at("position").text()) -
                                  exact_decimal(old.declaration.at("time").at("position").text());
        terms.set("finite_rate", scalar(temporal::finite_rate(change.total, interval)));
        terms.set("interval", scalar(interval));
        terms.set("time_unit",Value(r.declaration.at("time").at("unit")));
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
    auto fields = Value::array();
    for (const auto &c : channels) {
        auto e = Value::object();
        e.set("label", Value::string(c.label));
        e.set("field",vector(c.target_field,axis(readings.size()-1,readings.back().grade.value)));
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
    const auto terminal_channels_axis=native ? native->publication.axis(channels.size()) : 0;
    out.set("terminal_channel_gram",matrix(gram,terminal_channels_axis,terminal_channels_axis));
    out.set("terminal_quadrance", scalar(sum));
    if (factored)
        out.set("factored_maps", registry.archive());
    out.set("identity_verified", Value::boolean(true));
    return out;
}
Json temporal_field_history(const Json &query,const std::vector<TemporalReading> &readings,
    const std::vector<TemporalMaps> &maps,bool validate) {
    return collect_temporal_history<Json>(query,readings,maps,validate,nullptr);
}
ResultValue temporal_history_fields(const Json &query,const std::vector<TemporalReading> &readings,
    const std::vector<TemporalMaps> &maps,TemporalResultFields &fields) {
    return collect_temporal_history<ResultValue>(query,readings,maps,false,&fields);
}
scientific_detail::ResultValue censored_endpoint_bounds(const NativeStudy &study, const Json &query,
                              const std::map<std::string, Vector> &fields,
                              const std::vector<bool> &upper_support,
                              const std::map<std::string, Json> &meanings, const Json &coordinates,
                              bool validate, CensoredResultFields *native) {
    using scientific_detail::ResultValue;
    using scientific_detail::ResultFields;
    require(validate || native, "censored native result owner required");
    const auto rat = [&](const Rational &value) {
        return native ? native->publication.stage_scalar(value) : ResultValue{};
    };
    const auto vec = [&](const Vector &value, std::uint64_t axis) {
        return ResultFields::stage_vector(value, axis);
    };
    const auto mat = [&](const Matrix &value, std::uint64_t input, std::uint64_t output) {
        return ResultFields::stage_matrix(value, input, output);
    };
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
    auto groups = ResultValue::array();
    std::vector<bool> group_defined;
    for (std::size_t g = 0; g < readout.rows(); ++g) {
        Vector membership = selected;
        if (g)
            for (std::size_t i = 0; i < n; ++i)
                membership[i] *= fields.at(query.at("strata").at(g - 1).text())[i];
        Rational mass{0};
        for (const auto &v : membership)
            mass += v;
        auto group = ResultValue::object();
        group.set("id", g ? ResultValue(query.at("strata").at(g - 1)) : ResultValue::string("selected_population"));
        group.set("size", rat(mass));
        group_defined.push_back(mass != 0);
        group.set("defined", ResultValue::boolean(mass != 0));
        groups.push(group);
        if (mass != 0)
            for (std::size_t i = 0; i < n; ++i)
                if (membership[i] != 0)
                    readout.set(g, i, membership[i] / mass);
    }
    auto contrasts = ResultValue::array();
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
        auto record = ResultValue(d);
        const bool defined =
            group_defined.at(*left) && group_defined.at(*right);
        record.set("defined", ResultValue::boolean(defined));
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
    ResultValue joint;
    if (query.has("joint_family"))
        joint =
            joint_endpoint_family(study, query.at("joint_family"), l, u, lc, uc, uf, selected,
                                  times, origin, readout, contrast_readout, coordinates, validate, native);
    if (validate)
        return ResultValue::object();
    auto observations = ResultValue::array();
    for (std::size_t i = 0; i < n; ++i)
        if (selected[i] != 0) {
            auto r = ResultValue::object();
            r.set("coordinate", coordinates.at(i));
            r.set("lower", rat(l[i]));
            r.set("upper", uf[i] != 0 ? rat(u[i]) : ResultValue{});
            r.set("lower_closed", ResultValue::boolean(lc[i] == 1));
            r.set("upper_closed", ResultValue::boolean(uc[i] == 1));
            r.set("first_possible_grid_index", ResultValue::count(times.size()));
            r.set("first_definite_grid_index", ResultValue::count(times.size()));
            for (std::size_t t = 0; t < times.size(); ++t) {
                if (l[i] < times[t] || (l[i] == times[t] && lc[i] == 1)) {
                    r.set("first_possible_grid_index", ResultValue::count(t));
                    break;
                }
            }
            for (std::size_t t = 0; t < times.size(); ++t)
                if (uf[i] != 0 && u[i] <= times[t]) {
                    r.set("first_definite_grid_index", ResultValue::count(t));
                    break;
                }
            observations.push(r);
        }
    const auto source_axis = native->source, group_axis = native->groups, contrast_axis = native->contrasts;
    const PreparedRelationalField hodge(study, core::Grade{0});
    auto readings = ResultValue::array();
    for (const auto &s : times) {
        Vector emin(n, Rational{0}), emax(n, Rational{0}), rmin(n, Rational{0}),
            rmax(n, Rational{0});
        auto lower_attained = ResultValue::array(), upper_attained = ResultValue::array();
        for (std::size_t i = 0; i < n; ++i) {
            if (selected[i] != 0) {
                emin[i] = uf[i] != 0 && u[i] <= s ? 1 : 0;
                emax[i] = l[i] < s || (l[i] == s && lc[i] == 1) ? 1 : 0;
                rmin[i] = std::min(l[i], s) - origin;
                rmax[i] = (uf[i] != 0 ? std::min(u[i], s) : s) - origin;
            }
            lower_attained.push(ResultValue::boolean(selected[i] == 0 || lc[i] == 1 || l[i] >= s));
            upper_attained.push(
                ResultValue::boolean(selected[i] == 0 || uf[i] == 0 || uc[i] == 1 || u[i] > s));
        }
        auto r = ResultValue::object();
        r.set("time", rat(s));
        r.set("event_fraction_contrast_lower",
              vec(numeric::add(positive.apply(emin), negative.apply(emax)), contrast_axis));
        r.set("event_fraction_contrast_upper",
              vec(numeric::add(positive.apply(emax), negative.apply(emin)), contrast_axis));
        r.set("restricted_time_contrast_infimum",
              vec(numeric::add(positive.apply(rmin), negative.apply(rmax)), contrast_axis));
        r.set("restricted_time_contrast_supremum",
              vec(numeric::add(positive.apply(rmax), negative.apply(rmin)), contrast_axis));
        r.set("event_field_lower", vec(emin, source_axis));
        r.set("event_field_upper", vec(emax, source_axis));
        r.set("cumulative_event_fraction_lower", vec(readout.apply(emin), group_axis));
        r.set("cumulative_event_fraction_upper", vec(readout.apply(emax), group_axis));
        r.set("restricted_event_free_time_infimum", vec(readout.apply(rmin), group_axis));
        r.set("restricted_event_free_time_supremum", vec(readout.apply(rmax), group_axis));
        r.set("restricted_time_lower_field", vec(rmin, source_axis));
        r.set("restricted_time_upper_field", vec(rmax, source_axis));
        r.set("lower_coordinate_bound_attained", lower_attained);
        r.set("upper_coordinate_bound_attained", upper_attained);
        std::vector<std::size_t> ambiguous;
        for (std::size_t i = 0; i < n; ++i)
            if (emin[i] != emax[i])
                ambiguous.push_back(i);
        const auto parameter_axis = native->publication.axis(ambiguous.size());
        auto parameter_binding = ResultValue::object();
        parameter_binding.set("kind", ResultValue::string("independent_event_indicator_parameters"));
        parameter_binding.set("time", rat(s));
        parameter_binding.set("source_carrier", ResultValue::count(source_axis));
        parameter_binding.set("parameter_carrier", ResultValue::count(parameter_axis));
        parameter_binding.set("source_coordinates", native->publication.indices(ambiguous, parameter_axis));
        native->publication.bind(std::move(parameter_binding));
        Matrix frame(n, ambiguous.size());
        for (std::size_t j = 0; j < ambiguous.size(); ++j)
            frame.set(ambiguous[j], j, Rational{1});
        auto family = ResultValue::object();
        family.set("offset", vec(emin, source_axis));
        family.set("frame", mat(frame, parameter_axis, source_axis));
        family.set(
            "parameter_domain",
            ResultValue::string("independent zero-or-one event indicators at this evaluation time"));
        auto sectors = ResultValue::object();
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
            auto a = ResultValue::object();
            a.set("offset", vec(offsets[k], source_axis));
            a.set("frame", mat(sector_frames[k], parameter_axis, source_axis));
            Vector lower = offsets[k], upper = offsets[k];
            for (const auto &[p, v] : sector_frames[k].entries()) {
                if (v < 0)
                    lower[p.first] += v;
                else
                    upper[p.first] += v;
            }
            a.set("coordinate_lower", vec(lower, source_axis));
            a.set("coordinate_upper", vec(upper, source_axis));
            sectors.set(names[k], a);
        }
        family.set("sectors", sectors);
        r.set("admissible_event_family", family);

        readings.push(r);
    }
    auto out = ResultValue::object();
    out.set("population", query.at("population"));
    out.set("event", query.at("event"));
    out.set("origin", rat(origin));
    out.set("time_unit", query.at("time_unit"));
    out.set("groups", groups);
    out.set("contrasts", contrasts);
    out.set("contrast_readout", mat(contrast_readout, source_axis, contrast_axis));
    out.set("readout", mat(readout, source_axis, group_axis));
    out.set("observations", observations);
    out.set("readings", readings);
    out.set("coordinate_carrier", coordinates);
    auto reader = ResultValue::object(), aggregate_readings = ResultValue::array();
    reader.set(
        "scope",
        ResultValue::string("Finite-population marginal event and restricted-time bounds; original "
                     "coordinates and coupled field families remain in the source certificate"));
    reader.set("groups", groups);
    reader.set("contrasts", contrasts);
    reader.set("time_unit", query.at("time_unit"));
    for (std::size_t i = 0; i < readings.size(); ++i) {
        auto e = ResultValue::object();
        e.set("full_result_index", ResultValue::count(i));
        for (const auto *key :
             {"time", "cumulative_event_fraction_lower", "cumulative_event_fraction_upper",
              "restricted_event_free_time_infimum", "restricted_event_free_time_supremum",
              "event_fraction_contrast_lower", "event_fraction_contrast_upper"})
            e.set(key, readings.at(i).at(key));
        aggregate_readings.push(e);
    }
    reader.set("readings", aggregate_readings);
    out.set("reader", reader);

    out.set("metric", mat(study.tower().metrics()[0].materialize(), source_axis, source_axis));
    out.set("sharp_for",
            ResultValue::string("Each displayed marginal over the independent declared event-time "
                         "intervals; bounds at different times or for different groups need not "
                         "share one attaining event-time assignment."));
    std::size_t unidentified = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (selected[i] != 0 && (uf[i] == 0 || u[i] != l[i]))
            ++unidentified;
    out.set("unidentified_event_count", ResultValue::count(unidentified));
    auto facts = ResultValue::object();
    facts.set("unidentified_event_count", ResultValue::count(unidentified));
    out.set("facts", facts);
    if (query.has("joint_family"))
        out.set("joint_family", joint);
    return out;
}
} // namespace gadget::research_detail
