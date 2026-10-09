#include "research_response.hpp"
#include "research_fields.hpp"
#include "scientific_model.hpp"
#include "scientific_reference.hpp"
#include "scientific_relational_field.hpp"
#include "scientific_shape.hpp"
#include <algorithm>
#include <array>
#include <gadget/resources.hpp>
#include <memory>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::research_detail {
namespace {
using Matrix = numeric::ExactMatrix;
using Model = scientific_detail::PreparedObservedModel;
constexpr std::array<const char *, 4> sectors{"full", "gradient", "curl", "harmonic"};
Json rat(Rational x) {
    x.canonicalize();
    return Json::string(x.get_str());
}
Json vec(const Vector &x) {
    auto a = Json::array();
    for (const auto &v : x)
        a.push(rat(v));
    return a;
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
struct Response {
    std::string id, kind, field, denominator;
    Json declaration;
    Vector codes;
};
struct Axis {
    std::size_t response, cut;
    std::string id;
};
struct Variant {
    std::string id;
    std::set<std::string> remove;
};
struct State {
    std::string id;
    std::vector<std::string> fields;
    Rational mass;
};
struct Prepared {
    Matrix readout;
    std::vector<std::unique_ptr<Model>> variants;
};
Json reading(const scientific_detail::PreparedRelationalField &hodge, const Vector &x,
             const std::optional<Rational> &lambda, const runtime::CpuTileStop &stop) {
    auto out = Json::object(), fs = Json::object(), qs = Json::object();
    const auto parts = hodge.apply(x, [&] { stop.check(); });
    for (std::size_t k = 0; k < parts.size(); ++k) {
        fs.set(sectors[k], vec(parts[k]));
        qs.set(sectors[k], rat(hodge.metric.inner(parts[k], parts[k])));
    }
    out.set("sectors", fs);
    out.set("quadrances", qs);
    out.set("hodge_action", vec(hodge.study.hodge(hodge.grade).apply(x)));
    if (hodge.grade.value)
        out.set("lower_boundary", vec(hodge.study.source().boundary_apply(hodge.grade, x)));
    if (hodge.grade.value < hodge.study.source().max_grade().value) {
        const auto upper = core::Grade{static_cast<std::uint16_t>(hodge.grade.value + 1u)};
        out.set("upper_adjoint",
                vec(field::BoundaryOperator(hodge.study.source(), upper)
                        .adjoint(x, hodge.metric, hodge.study.tower().metrics()[upper.value])));
    }
    if (lambda) {
        const auto y = hodge.study.green_action(hodge.grade, *lambda, x).value;
        require(numeric::add(y, numeric::scale(*lambda, hodge.study.hodge(hodge.grade).apply(y))) ==
                    x,
                "response Green action residual");
        out.set("green", vec(y));
    }
    return out;
}
} // namespace

std::set<std::string> response_model_field_ids(const Json &query) {
    std::set<std::string> ids;
    const auto responses = query.at("responses");
    source_array(responses, "response declarations");
    for (std::size_t r = 0; r < responses.size(); ++r) {
        ids.insert(text(responses.at(r).at("field"), source_byte_limit()));
        if (responses.at(r).has("denominator"))
            ids.insert(text(responses.at(r).at("denominator"), source_byte_limit()));
    }
    const auto terms = query.at("model").at("terms");
    source_array(terms, "response model terms");
    for (std::size_t t = 0; t < terms.size(); ++t)
        for (std::size_t f = 0; f < terms.at(t).at("factors").size(); ++f)
            ids.insert(text(terms.at(t).at("factors").at(f), source_byte_limit()));
    if (query.has("dose_comparisons")) {
        const auto doses = query.at("dose_comparisons");
        source_array(doses, "dose comparisons");
        for (std::size_t d = 0; d < doses.size(); ++d) {
            ids.insert(text(doses.at(d).at("field"), source_byte_limit()));
            const auto groups = doses.at(d).at("groups");
            source_array(groups, "dose groups");
            for (std::size_t g = 0; g < groups.size(); ++g)
                ids.insert(text(groups.at(g).at("indicator"), source_byte_limit()));
        }
    }
    if (query.has("reference_law")) {
        const auto states = query.at("reference_law").at("states");
        source_array(states, "reference states");
        for (std::size_t s = 0; s < states.size(); ++s) {
            const auto bindings = states.at(s).at("response_fields");
            require(json_object_is_type(bindings.raw(), json_type_object),
                    "joint response bindings must be an object");
            json_object_object_foreach(bindings.raw(), key, value) {
                (void)key;
                ids.insert(text(Json(json_object_get(value)), source_byte_limit()));
            }
        }
    }
    return ids;
}

Json observed_response_model(const NativeStudy &study, const Json &query, core::Grade grade,
                             const std::map<std::string, Vector> &fields,
                             const std::map<std::string, Json> &meanings, const Json &coordinates,
                             bool validate, const std::shared_ptr<CpuResources> &resources) {
    const bool constrained = query.at("operation").text() == "constrained_response_model";
    std::set<std::string> required{"operation", "context", "responses",
                                   "model",     "workers", "interpretation"};
    if (constrained)
        required.insert("physical_shape");
    members(query, required, {"dose_comparisons", "reference_law", "green_parameter"});
    const auto interpretation = query.at("interpretation");
    members(interpretation, {"population", "meaning", "evidence"});
    (void)text(interpretation.at("population"), source_byte_limit());
    (void)text(interpretation.at("meaning"), source_byte_limit());
    const auto n = study.source().grade(grade).cell_count();
    require(n && coordinates.size() == n, "response source coordinate family");
    const auto &metric = study.tower().metrics().at(grade.value);
    for (const auto &[id, x] : fields) {
        (void)id;
        require(x.size() == n, "complete response source field");
    }
    std::optional<Rational> lambda;
    if (query.has("green_parameter")) {
        lambda = exact_rational(text(query.at("green_parameter"), source_byte_limit()));
        require(*lambda > 0, "response Green parameter must be positive");
    }
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        count(query.at("workers"), cpu->coordinator().capacity(), "response workers");
    require(workers > 0, "response workers must be positive");
    const auto declarations = query.at("responses");
    require(declarations.size(), "response declarations required");
    std::vector<Response> responses;
    std::vector<Axis> axes;
    std::set<std::string> response_ids, response_fields;
    for (std::size_t r = 0; r < declarations.size(); ++r) {
        const auto d = declarations.at(r);
        members(d, {"id", "kind", "field", "meaning", "evidence"},
                {"denominator", "rate_units", "levels"});
        Response response{text(d.at("id"), source_byte_limit()),
                          d.at("kind").text(),
                          d.at("field").text(),
                          {},
                          d,
                          {}};
        require(response_ids.insert(response.id).second, "duplicate response identity");
        response_fields.insert(response.field);
        (void)text(d.at("meaning"), source_byte_limit());
        if (response.kind == "count" || response.kind == "binomial_count") {
            require(d.has("denominator") && d.has("rate_units") && !d.has("levels"),
                    "count requires an explicit denominator and rate units");
            response.denominator = d.at("denominator").text();
            require(response.denominator != response.field,
                    "count cannot be its own exposure or trials");
            (void)text(d.at("rate_units"), source_byte_limit());
            require(meanings.at(response.field).at("value_type").text() == "quantitative" &&
                        meanings.at(response.denominator).at("value_type").text() == "quantitative",
                    "count and denominator must be quantitative source fields");
            if (response.kind == "binomial_count")
                require(d.at("rate_units").text() == "dimensionless" &&
                            meanings.at(response.field).at("units").text() ==
                                meanings.at(response.denominator).at("units").text(),
                        "successes and trials must have common counting units");
            for (const auto &e : fields.at(response.denominator))
                require(e > 0 && (response.kind != "binomial_count" || e.get_den() == 1),
                        "exposure must be positive; trials must be integral");
            axes.push_back({r, 0, response.id});
        } else if (response.kind == "ordinal") {
            require(d.has("levels") && !d.has("denominator") && !d.has("rate_units"),
                    "ordinal requires explicit ordered levels only");
            require(meanings.at(response.field).at("value_type").text() == "categorical",
                    "ordinal response must be declared categorical");
            const auto levels = d.at("levels");
            source_array(levels, "ordered levels");
            require(levels.size() > 1, "ordinal requires at least two levels");
            std::set<Rational> codes;
            std::set<std::string> labels;
            for (std::size_t k = 0; k < levels.size(); ++k) {
                const auto level = levels.at(k);
                members(level, {"code", "label"});
                auto code = exact_rational(text(level.at("code"), source_byte_limit()));
                require(codes.insert(code).second &&
                            labels.insert(text(level.at("label"), source_byte_limit())).second,
                        "duplicate ordinal code or label");
                response.codes.push_back(code);
                if (k)
                    axes.push_back({r, k, response.id + ":" + std::to_string(k)});
            }
        } else {
            require(response.kind == "continuous" && !d.has("levels") && !d.has("denominator") &&
                        !d.has("rate_units"),
                    "unknown or incompatible response kind");
            require(meanings.at(response.field).at("value_type").text() == "quantitative",
                    "continuous response must be quantitative");
            axes.push_back({r, 0, response.id});
        }
        responses.push_back(std::move(response));
    }
    const auto model_query = query.at("model");
    members(model_query, {"intercept", "terms"}, {"ablations"});
    require(json_object_is_type(model_query.at("intercept").raw(), json_type_boolean),
            "model intercept must be Boolean");
    const auto terms = model_query.at("terms");
    std::vector<std::string> column_ids;
    std::vector<std::vector<std::string>> factors;
    if (model_query.at("intercept").boolean()) {
        column_ids.push_back("intercept");
        factors.push_back({});
    }
    std::set<std::string> term_ids, all_ids(column_ids.begin(), column_ids.end());
    for (std::size_t t = 0; t < terms.size(); ++t) {
        const auto term = terms.at(t);
        members(term, {"id", "factors"});
        const auto id = text(term.at("id"), source_byte_limit());
        require(all_ids.insert(id).second, "duplicate model column");
        term_ids.insert(id);
        column_ids.push_back(id);
        // Repeated factors declare powers; identity repetition is meaningful.
        const auto fs = term.at("factors");
        source_array(fs, "term factors");
        require(fs.size(), "term requires factors");
        std::vector<std::string> declared;
        for (std::size_t f = 0; f < fs.size(); ++f) {
            const auto field = text(fs.at(f), source_byte_limit());
            require(!response_fields.contains(field) &&
                        meanings.at(field).at("value_type").text() != "categorical",
                    "response or categorical code cannot be a model factor");
            declared.push_back(field);
        }
        factors.push_back(std::move(declared));
    }
    std::vector<Variant> variants{{"baseline", term_ids}, {"full", {}}};
    std::set<std::string> variant_ids{"baseline", "full"};
    if (model_query.has("ablations")) {
        const auto ablations = model_query.at("ablations");
        source_array(ablations, "model ablations");
        for (std::size_t a = 0; a < ablations.size(); ++a) {
            const auto d = ablations.at(a);
            members(d, {"id", "remove"});
            Variant v{text(d.at("id"), source_byte_limit()), {}};
            require(variant_ids.insert(v.id).second, "duplicate model variant");
            for (const auto &id : names(d.at("remove"), term_ids.size())) {
                require(term_ids.contains(id), "ablation must name a term");
                v.remove.insert(id);
            }
            variants.push_back(std::move(v));
        }
    }
    Matrix all(n, column_ids.size());
    auto columns = Json::array();
    for (std::size_t j = 0; j < column_ids.size(); ++j) {
        auto c = Json::object(), units = Json::array();
        c.set("id", Json::string(column_ids[j]));
        for (const auto &f : factors[j])
            units.push(meanings.at(f).at("units"));
        c.set("factor_units", units);
        columns.push(c);
        for (std::size_t i = 0; i < n; ++i) {
            Rational value{1};
            for (const auto &f : factors[j])
                value *= fields.at(f)[i];
            all.set(i, j, value);
        }
    }
    std::vector<State> states{{"observed", {}, Rational{1}}};
    for (const auto &r : responses)
        states[0].fields.push_back(r.field);
    std::size_t observed = 0;
    std::optional<Rational> alpha;
    if (query.has("reference_law")) {
        const auto law = query.at("reference_law");
        members(law, {"states", "observed_state", "meaning", "evidence"}, {"alpha"});
        (void)text(law.at("meaning"), source_byte_limit());
        const auto ds = law.at("states");
        require(ds.size(), "finite response law needs states");
        states.clear();
        std::set<std::string> ids;
        Rational total{0};
        bool found = false;
        for (std::size_t s = 0; s < ds.size(); ++s) {
            const auto d = ds.at(s);
            members(d, {"id", "mass", "response_fields"});
            State state{text(d.at("id"), source_byte_limit()),
                        {},
                        exact_rational(text(d.at("mass"), source_byte_limit()))};
            require(ids.insert(state.id).second && state.mass > 0,
                    "finite response law unique states and positive masses");
            const auto bs = d.at("response_fields");
            members(bs, response_ids);
            for (const auto &r : responses) {
                const auto f = bs.at(r.id).text();
                for (const auto *m : {"channel", "quantity", "units", "method", "value_type"})
                    require(meanings.at(f).at(m).dump() == meanings.at(r.field).at(m).dump(),
                            "reference response meaning mismatch");
                state.fields.push_back(f);
            }
            if (state.id == law.at("observed_state").text()) {
                observed = s;
                found = true;
            }
            total += state.mass;
            states.push_back(std::move(state));
        }
        require(found && total == 1, "observed state and complete unit mass required");
        for (std::size_t r = 0; r < responses.size(); ++r)
            require(fields.at(states[observed].fields[r]) == fields.at(responses[r].field),
                    "reference observed response must equal the original observation");
        if (law.has("alpha")) {
            alpha = exact_rational(text(law.at("alpha"), source_byte_limit()));
            require(*alpha > 0 && *alpha <= 1, "reference alpha domain");
        }
    }
    std::vector<std::vector<Vector>> observations(states.size());
    for (std::size_t s = 0; s < states.size(); ++s) {
        for (const auto &axis : axes) {
            const auto &r = responses[axis.response];
            const auto &x = fields.at(states[s].fields[axis.response]);
            Vector y = x;
            for (std::size_t i = 0; i < n; ++i) {
                if (r.kind == "ordinal") {
                    const auto found = std::find(r.codes.begin(), r.codes.end(), x[i]);
                    require(found != r.codes.end(),
                            "observed ordinal code outside declared levels");
                    y[i] = static_cast<std::size_t>(found - r.codes.begin()) >= axis.cut
                               ? Rational{1}
                               : Rational{0};
                } else if (!r.denominator.empty()) {
                    require(x[i] >= 0 && x[i].get_den() == 1,
                            "observed counts must be nonnegative integers");
                    if (r.kind == "binomial_count")
                        require(x[i] <= fields.at(r.denominator)[i],
                                "successes cannot exceed trials");
                }
            }
            observations[s].push_back(std::move(y));
        }
    }
    // Admit dose groups even during validate-only, before preparation/execution.
    struct Dose {
        Json declaration;
        Vector doses;
        Matrix means, lift, differences;
    };
    std::vector<Dose> doses;
    std::set<std::string> dose_ids;
    if (query.has("dose_comparisons")) {
        const auto ds = query.at("dose_comparisons");
        source_array(ds, "dose comparisons");
        for (std::size_t d = 0; d < ds.size(); ++d) {
            const auto decl = ds.at(d);
            members(decl, {"id", "field", "groups", "direction", "meaning", "evidence"});
            require(dose_ids.insert(text(decl.at("id"), source_byte_limit())).second,
                    "duplicate dose comparison");
            (void)text(decl.at("meaning"), source_byte_limit());
            const auto field = decl.at("field").text();
            require(meanings.at(field).at("value_type").text() == "quantitative",
                    "dose must be quantitative");
            const auto direction = decl.at("direction").text();
            require(direction == "nondecreasing" || direction == "nonincreasing" ||
                        direction == "unrestricted",
                    "dose direction must be declared");
            const auto groups = decl.at("groups");
            source_array(groups, "dose groups");
            require(groups.size() > 1, "dose comparison needs multiple dose groups");
            Dose dose{decl,
                      {},
                      Matrix(groups.size(), n),
                      Matrix(n, groups.size()),
                      Matrix(groups.size() - 1, groups.size())};
            std::set<std::string> group_ids;
            std::set<std::size_t> covered;
            for (std::size_t g = 0; g < groups.size(); ++g) {
                const auto group = groups.at(g);
                members(group, {"id", "indicator"});
                require(group_ids.insert(text(group.at("id"), source_byte_limit())).second,
                        "duplicate dose group");
                const auto selector = group.at("indicator").text();
                require(meanings.at(selector).at("value_type").text() == "indicator",
                        "dose group requires declared indicator");
                std::vector<std::size_t> rows;
                for (std::size_t i = 0; i < n; ++i)
                    if (fields.at(selector)[i] == 1) {
                        require(covered.insert(i).second, "dose groups must be disjoint");
                        rows.push_back(i);
                    }
                require(!rows.empty(), "empty dose group");
                const auto level = fields.at(field)[rows[0]];
                require(level >= 0, "dose must be nonnegative");
                if (g)
                    require(level > dose.doses.back(),
                            "dose groups must be ordered with distinct increasing dose levels");
                for (auto i : rows) {
                    require(fields.at(field)[i] == level, "dose group mixes measured dose levels");
                    dose.means.set(g, i, Rational{1} / Rational{numeric::Integer{rows.size()}});
                    dose.lift.set(i, g, Rational{1});
                }
                dose.doses.push_back(level);
                if (g) {
                    const Rational spacing = level - dose.doses[g - 1];
                    dose.differences.set(g - 1, g, Rational{1} / spacing);
                    dose.differences.set(g - 1, g - 1, -Rational{1} / spacing);
                }
            }
            require(covered.size() == n, "dose groups must cover every context coordinate");
            require(dose.means.multiply(dose.lift) == Matrix::identity(groups.size()),
                    "dose mean/lift reconstruction");
            doses.push_back(std::move(dose));
        }
    }
    std::vector<std::pair<std::size_t, std::size_t>> shape_doses;
    if (constrained) {
        const auto shape = query.at("physical_shape");
        members(shape, {"meaning", "evidence", "dose_constraints"});
        (void)text(shape.at("meaning"), source_byte_limit());
        source_array(shape.at("dose_constraints"), "physical dose constraints");
        std::set<std::pair<std::string, std::string>> seen;
        for (std::size_t i = 0; i < shape.at("dose_constraints").size(); ++i) {
            const auto d = shape.at("dose_constraints").at(i);
            members(d, {"response", "comparison"});
            const auto r = d.at("response").text(), c = d.at("comparison").text();
            require(seen.emplace(r, c).second && response_ids.contains(r) && dose_ids.contains(c),
                    "physical dose constraint requires unique declared response/comparison");
            const auto ri = std::find_if(responses.begin(), responses.end(),
                                         [&](const auto &v) { return v.id == r; });
            const auto di = std::find_if(doses.begin(), doses.end(), [&](const auto &v) {
                return v.declaration.at("id").text() == c;
            });
            require(di->declaration.at("direction").text() != "unrestricted",
                    "physical dose constraint needs a declared direction");
            shape_doses.emplace_back(static_cast<std::size_t>(ri - responses.begin()),
                                     static_cast<std::size_t>(di - doses.begin()));
        }
    }
    if (validate) {
        auto out = Json::object();
        out.set("status", Json::string("validated"));
        return out;
    }
    auto reservation = cpu->reserve(workers);
    const scientific_detail::PreparedRelationalField hodge(study, grade);
    std::map<std::string, std::unique_ptr<Prepared>> prepared;
    auto preparation_records = Json::array();
    for (const auto &r : responses) {
        if (prepared.contains(r.denominator))
            continue;
        Vector exposure(n, Rational{1}), inverse(n, Rational{1});
        if (!r.denominator.empty())
            exposure = fields.at(r.denominator);
        for (std::size_t i = 0; i < n; ++i)
            inverse[i] = Rational{1} / exposure[i];
        const auto e = Matrix::diagonal(exposure);
        auto p = std::make_unique<Prepared>();
        p->readout = Matrix::diagonal(inverse);
        auto record = Json::object(), models = Json::array();
        record.set("denominator_field",
                   r.denominator.empty() ? Json{} : Json::string(r.denominator));
        record.set("exposure", vec(exposure));
        record.set("readout", mat(p->readout));
        record.set("readout_metric_pullback_map", mat(e));
        for (const auto &variant : variants) {
            std::vector<std::size_t> selected;
            auto axis = Json::array();
            for (std::size_t j = 0; j < column_ids.size(); ++j)
                if (!variant.remove.contains(column_ids[j])) {
                    selected.push_back(j);
                    axis.push(Json::string(column_ids[j]));
                }
            Matrix a(n, selected.size());
            for (std::size_t j = 0; j < selected.size(); ++j)
                for (std::size_t i = 0; i < n; ++i)
                    a.set(i, j, all.get(i, selected[j]));
            p->variants.push_back(std::make_unique<Model>(e.multiply(a), metric));
            const auto &m = *p->variants.back();
            auto model = Json::object();
            model.set("id", Json::string(variant.id));
            model.set("columns", axis);
            model.set("design", mat(m.design));
            model.set("independent_frame", mat(m.projection.frame));
            model.set("coefficient_kernel", mat(m.kernel));
            model.set("rank", cardinality(m.pivots.size()));
            model.set("coefficient_nullity", cardinality(m.kernel.cols()));
            auto identified = Json::array();
            for (std::size_t j = 0; j < m.design.cols(); ++j) {
                bool unique = true;
                for (std::size_t k = 0; k < m.kernel.cols(); ++k)
                    unique = unique && m.kernel.get(j, k) == 0;
                identified.push(Json::boolean(unique));
            }
            model.set("coefficient_identified", identified);
            models.push(model);
        }
        record.set("models", models);
        preparation_records.push(record);
        prepared.emplace(r.denominator, std::move(p));
    }
    struct Job {
        Model::Fit fit;
        Vector readout;
        Json record;
    };
    const auto jobs = machine::checked_size_mul(
        machine::checked_size_mul(states.size(), axes.size()), variants.size());
    std::vector<Job> completed(jobs);
    auto index = [&](std::size_t s, std::size_t a, std::size_t v) {
        return (s * axes.size() + a) * variants.size() + v;
    };
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), jobs});
    const auto report = executor.run_tiles(jobs, [&](std::size_t j, std::size_t,
                                                     const runtime::CpuTileStop &stop) {
        const auto v = j % variants.size(), a = (j / variants.size()) % axes.size(),
                   s = j / (variants.size() * axes.size());
        const auto &axis = axes[a];
        const auto &r = responses[axis.response];
        const auto &p = *prepared.at(r.denominator);
        const auto &model = *p.variants[v];
        auto &job = completed[j];
        job.fit = model.apply(observations[s][a], [&] { stop.check(); });
        const auto x = p.readout.apply(observations[s][a]);
        job.readout = p.readout.apply(job.fit.represented);
        const auto residual = p.readout.apply(job.fit.residual);
        auto out = Json::object();
        out.set("axis", Json::string(axis.id));
        out.set("response", Json::string(r.id));
        out.set("model", Json::string(variants[v].id));
        out.set("readout_units", r.kind == "ordinal"
                                     ? Json::string("dimensionless")
                                     : (!r.denominator.empty() ? r.declaration.at("rate_units")
                                                               : meanings.at(r.field).at("units")));
        out.set("objective_input", vec(observations[s][a]));
        out.set("objective_represented", vec(job.fit.represented));
        out.set("objective_residual", vec(job.fit.residual));
        out.set("particular_coefficients", vec(job.fit.coefficients));
        out.set("residual_quadrance", rat(job.fit.quadrance));
        out.set("normal_equations",
                vec(model.design.transpose().apply(metric.apply(job.fit.residual))));
        out.set("observed", reading(hodge, x, lambda, stop));
        out.set("represented", reading(hodge, job.readout, lambda, stop));
        out.set("residual", reading(hodge, residual, lambda, stop));
        auto violations = Json::array();
        if (!r.denominator.empty() || r.kind == "ordinal")
            for (std::size_t i = 0; i < n; ++i)
                if (job.readout[i] < 0 ||
                    ((r.kind == "ordinal" || r.kind == "binomial_count") && job.readout[i] > 1)) {
                    auto e = Json::object();
                    e.set("coordinate", cardinality(i));
                    e.set("value", rat(job.readout[i]));
                    violations.push(e);
                }
        out.set("range_violations", violations);
        out.set("range_admissible", Json::boolean(violations.size() == 0));
        job.record = out;
    });
    auto state_records = Json::array(), family = Json::array();
    std::vector<Vector> scores(states.size(), Vector(axes.size(), Rational{0}));
    std::vector<std::vector<bool>> defined(states.size(), std::vector<bool>(axes.size(), false));
    for (const auto &axis : axes) {
        auto f = Json::object();
        f.set("id", Json::string(axis.id));
        f.set("statistic", Json::string("baseline_to_full_objective_fraction"));
        family.push(f);
    }
    for (std::size_t s = 0; s < states.size(); ++s) {
        auto state = Json::object(), fits = Json::array(), comparisons = Json::array(),
             categories = Json::array(), dose_records = Json::array();
        state.set("id", Json::string(states[s].id));
        state.set("mass", query.has("reference_law") ? rat(states[s].mass) : Json{});
        auto bindings = Json::object();
        for (std::size_t r = 0; r < responses.size(); ++r)
            bindings.set(responses[r].id, Json::string(states[s].fields[r]));
        state.set("response_fields", bindings);
        for (std::size_t a = 0; a < axes.size(); ++a) {
            for (std::size_t v = 0; v < variants.size(); ++v)
                fits.push(completed[index(s, a, v)].record);
            const auto &full = completed[index(s, a, 1)].fit;
            for (std::size_t v = 0; v < variants.size(); ++v)
                if (v != 1) {
                    const auto &reduced = completed[index(s, a, v)].fit;
                    const auto d = numeric::sub(full.represented, reduced.represented);
                    const Rational gain = reduced.quadrance - full.quadrance;
                    require(gain >= 0 && metric.inner(d, d) == gain &&
                                metric.inner(d, full.residual) == 0,
                            "nested response projection reduction");
                    auto c = Json::object();
                    c.set("axis", Json::string(axes[a].id));
                    c.set("reduced", Json::string(variants[v].id));
                    c.set("gain", rat(gain));
                    c.set("represented_change", vec(d));
                    c.set("fraction",
                          reduced.quadrance > 0 ? rat(gain / reduced.quadrance) : Json{});
                    comparisons.push(c);
                    if (v == 0 && reduced.quadrance > 0) {
                        defined[s][a] = true;
                        scores[s][a] = gain / reduced.quadrance;
                    }
                }
        }
        for (std::size_t r = 0; r < responses.size(); ++r)
            if (responses[r].kind == "ordinal") {
                std::vector<std::size_t> cut_axes;
                for (std::size_t a = 0; a < axes.size(); ++a)
                    if (axes[a].response == r)
                        cut_axes.push_back(a);
                for (std::size_t v = 0; v < variants.size(); ++v) {
                    auto c = Json::object(), fields_out = Json::array(), violations = Json::array();
                    c.set("response", Json::string(responses[r].id));
                    c.set("model", Json::string(variants[v].id));
                    for (std::size_t k = 0; k < responses[r].codes.size(); ++k) {
                        Vector values(n, Rational{0});
                        for (std::size_t i = 0; i < n; ++i) {
                            const Rational lower =
                                k ? completed[index(s, cut_axes[k - 1], v)].readout[i]
                                  : Rational{1};
                            const Rational upper =
                                k < cut_axes.size() ? completed[index(s, cut_axes[k], v)].readout[i]
                                                    : Rational{0};
                            values[i] = lower - upper;
                            if (values[i] < 0 || values[i] > 1) {
                                auto e = Json::object();
                                e.set("level", cardinality(k));
                                e.set("coordinate", cardinality(i));
                                e.set("value", rat(values[i]));
                                violations.push(e);
                            }
                        }
                        fields_out.push(vec(values));
                    }
                    c.set("level_fields", fields_out);
                    c.set("violations", violations);
                    c.set("admissible", Json::boolean(violations.size() == 0));
                    categories.push(c);
                }
            }
        for (const auto &dose : doses) {
            auto dr = Json::object(), values = Json::array();
            dr.set("id", dose.declaration.at("id"));
            dr.set("dose_levels", vec(dose.doses));
            dr.set("dose_units", meanings.at(dose.declaration.at("field").text()).at("units"));
            dr.set("mean_map", mat(dose.means));
            dr.set("lift_map", mat(dose.lift));
            dr.set("divided_difference_map", mat(dose.differences));
            for (std::size_t a = 0; a < axes.size(); ++a) {
                const auto &r = responses[axes[a].response];
                const auto &p = *prepared.at(r.denominator);
                const auto x = p.readout.apply(observations[s][a]);
                const auto means = dose.means.apply(x), slopes = dose.differences.apply(means);
                const auto within = numeric::sub(x, dose.lift.apply(means));
                require(numeric::is_zero(dose.means.apply(within)),
                        "dose within-group mean residual");
                auto vr = Json::object();
                vr.set("axis", Json::string(axes[a].id));
                vr.set("means", vec(means));
                vr.set("divided_differences", vec(slopes));
                vr.set("within_group_field", vec(within));
                bool valid = true;
                const auto direction = dose.declaration.at("direction").text();
                for (const auto &v : slopes)
                    valid = valid && (direction == "unrestricted" ||
                                      (direction == "nondecreasing" && v >= 0) ||
                                      (direction == "nonincreasing" && v <= 0));
                vr.set("declared_direction_satisfied", Json::boolean(valid));
                values.push(vr);
            }
            dr.set("responses", values);
            dose_records.push(dr);
        }
        state.set("fits", fits);
        state.set("nested_comparisons", comparisons);
        state.set("ordinal_level_readings", categories);
        state.set("dose_comparisons", dose_records);
        state_records.push(state);
    }
    auto out = Json::object();
    out.set("status", Json::string("completed"));
    out.set("grade", cardinality(grade.value));
    auto dimensions = Json::object();
    dimensions.set("gradient", cardinality(hodge.projections[0].frame.cols()));
    dimensions.set("curl", cardinality(hodge.projections[1].frame.cols()));
    dimensions.set("harmonic", cardinality(n - hodge.projections[0].frame.cols() -
                                           hodge.projections[1].frame.cols()));
    out.set("hodge_dimensions", dimensions);
    auto boundaries = Json::array();
    for (std::size_t k = 1; k < study.source().grade_count(); ++k)
        boundaries.push(mat(field::materialize_boundary_matrix(
            study.source(), core::Grade{static_cast<std::uint16_t>(k)})));
    out.set("source_boundaries", boundaries);
    out.set("coordinates", coordinates);
    out.set("responses", declarations);
    out.set("columns", columns);
    out.set("preparations", preparation_records);
    out.set("states", state_records);
    out.set("observed_state_index", cardinality(observed));
    if (constrained) {
        struct ShapeCarrier {
            Matrix frame, source_normals, normals;
            Vector rhs;
            Json rows;
            field::ExactMetric gram;
            std::size_t family;
            std::vector<std::size_t> representatives;
        };
        struct Shape {
            std::size_t response, variant;
            std::vector<std::size_t> axes;
            std::shared_ptr<const ShapeCarrier> maps;
        };
        std::vector<Shape> shapes;
        std::map<std::string, std::shared_ptr<const ShapeCarrier>> families;
        auto shape_preparations = Json::array(), shape_families = Json::array();
        for (std::size_t r = 0; r < responses.size(); ++r) {
            const auto &response = responses[r];
            const auto &p = *prepared.at(response.denominator);
            std::vector<std::size_t> selected_axes;
            for (std::size_t a = 0; a < axes.size(); ++a)
                if (axes[a].response == r)
                    selected_axes.push_back(a);
            const auto dimension = machine::checked_size_mul(n, selected_axes.size());
            std::vector<std::map<std::size_t, Rational>> source_rows;
            Vector rhs;
            auto row_records = Json::array();
            const auto row = [&](std::map<std::size_t, Rational> values, Rational bound,
                                 const std::string &kind, std::size_t coordinate, std::size_t axis,
                                 Json evidence = Json{}) {
                source_rows.push_back(std::move(values));
                rhs.push_back(bound);
                auto record = Json::object();
                record.set("kind", Json::string(kind));
                record.set("coordinate", cardinality(coordinate));
                record.set("axis_or_level", cardinality(axis));
                record.set("declaration", evidence);
                row_records.push(record);
            };
            for (std::size_t i = 0; i < n; ++i) {
                if (response.kind == "count" || response.kind == "binomial_count") {
                    row({{i, -1}}, 0, "nonnegative_count", i, 0);
                    if (response.kind == "binomial_count")
                        row({{i, 1}}, fields.at(response.denominator)[i], "trial_upper_bound", i,
                            0);
                } else if (response.kind == "ordinal") {
                    row({{i, 1}}, 1, "ordinal_category_nonnegative", i, 0);
                    for (std::size_t k = 1; k < selected_axes.size(); ++k)
                        row({{(k - 1) * n + i, -1}, {k * n + i, 1}}, 0,
                            "ordinal_category_nonnegative", i, k);
                    row({{(selected_axes.size() - 1) * n + i, -1}}, 0,
                        "ordinal_category_nonnegative", i, selected_axes.size());
                }
            }
            for (const auto &[ri, di] : shape_doses) {
                if (ri != r)
                    continue;
                const auto &dose = doses[di];
                const auto map = dose.differences.multiply(dose.means).multiply(p.readout);
                const Rational sign = dose.declaration.at("direction").text() == "nondecreasing"
                                          ? Rational{-1}
                                          : Rational{1};
                for (std::size_t a = 0; a < selected_axes.size(); ++a)
                    for (std::size_t k = 0; k < map.rows(); ++k) {
                        std::map<std::size_t, Rational> entries;
                        for (std::size_t i = 0; i < n; ++i)
                            if (map.get(k, i) != 0)
                                entries.emplace(a * n + i, sign * map.get(k, i));
                        row(std::move(entries), 0, "measured_dose_direction", k, a,
                            dose.declaration);
                    }
            }
            Matrix c(source_rows.size(), dimension);
            for (std::size_t i = 0; i < source_rows.size(); ++i)
                for (const auto &[j, value] : source_rows[i])
                    c.set(i, j, value);
            for (std::size_t v = 0; v < variants.size(); ++v) {
                const auto &model = *p.variants[v];
                auto signature = Json::array();
                for (const auto &name : {response.kind, response.denominator, variants[v].id})
                    signature.push(Json::string(name));
                signature.push(mat(c));
                signature.push(vec(rhs));
                signature.push(row_records);
                const auto key = signature.dump();
                if (!families.contains(key)) {
                    const auto rank = model.projection.frame.cols();
                    Matrix frame(dimension, machine::checked_size_mul(rank, selected_axes.size()));
                    std::vector<field::ExactMetric> grams;
                    const auto base_gram =
                        field::metric_gram(model.projection.frame, metric, model.projection.frame);
                    for (std::size_t a = 0; a < selected_axes.size(); ++a) {
                        for (const auto &[coordinate, value] : model.projection.frame.entries())
                            frame.set(a * n + coordinate.first, a * rank + coordinate.second,
                                      value);
                        grams.emplace_back(base_gram);
                    }
                    const auto normals = c.multiply(frame);
                    auto maps = std::make_shared<const ShapeCarrier>(ShapeCarrier{
                        frame, c, normals, rhs, row_records,
                        field::ExactMetric::block_diagonal(std::move(grams)), shape_families.size(),
                        scientific_detail::physical_shape_representatives(normals, rhs)});
                    auto owner = Json::object();
                    owner.set("frame", mat(maps->frame));
                    owner.set("source_normals", mat(maps->source_normals));
                    owner.set("normals", mat(maps->normals));
                    owner.set("rhs", vec(maps->rhs));
                    owner.set("constraint_rows", maps->rows);
                    owner.set("gram", mat(maps->gram.materialize()));
                    auto representatives = Json::array();
                    for (auto i : maps->representatives)
                        representatives.push(cardinality(i));
                    owner.set("constraint_representatives", representatives);
                    shape_families.push(owner);
                    families.emplace(key, std::move(maps));
                }
                Shape shape{r, v, selected_axes, families.at(key)};
                auto record = Json::object(), ids = Json::array();
                record.set("response", Json::string(response.id));
                record.set("model", Json::string(variants[v].id));
                for (auto a : selected_axes)
                    ids.push(Json::string(axes[a].id));
                record.set("axes", ids);
                record.set("family", cardinality(shape.maps->family));
                shape_preparations.push(record);
                shapes.push_back(std::move(shape));
            }
        }
        const auto shape_jobs = machine::checked_size_mul(states.size(), shapes.size());
        std::vector<Json> shape_results(shape_jobs);
        std::vector<Rational> shape_objectives(shape_jobs);
        std::vector<Vector> shape_represented(shape_jobs);
        const auto shape_report = executor.run_tiles(
            shape_jobs, [&](std::size_t j, std::size_t, const runtime::CpuTileStop &stop) {
                const auto s = j / shapes.size(), owner = j % shapes.size();
                const auto &shape = shapes[owner];
                const auto &maps = *shape.maps;
                const auto &response = responses[shape.response];
                const auto &p = *prepared.at(response.denominator);
                const auto &model = *p.variants[shape.variant];
                const auto rank = model.pivots.size();
                Vector v0, input;
                Rational unbounded_objective{0};
                for (auto a : shape.axes) {
                    const auto &fit = completed[index(s, a, shape.variant)].fit;
                    for (auto pivot : model.pivots)
                        v0.push_back(fit.coefficients[pivot]);
                    input.insert(input.end(), observations[s][a].begin(), observations[s][a].end());
                    unbounded_objective += fit.quadrance;
                }
                const auto fit = scientific_detail::physical_shape_fit(
                    maps.normals, maps.rhs, maps.gram, v0, maps.representatives,
                    [&] { stop.check(); });
                const auto represented = maps.frame.apply(fit.coordinates);
                require(maps.normals.apply(fit.coordinates) ==
                            maps.source_normals.apply(represented),
                        "physical source/model constraint commutation");
                auto result = Json::object(), fits = Json::array(), active = Json::array();
                Rational objective{0};
                for (std::size_t a = 0; a < shape.axes.size(); ++a) {
                    Vector original(input.begin() + a * n, input.begin() + (a + 1) * n);
                    Vector predicted(represented.begin() + a * n,
                                     represented.begin() + (a + 1) * n);
                    const auto residual = numeric::sub(original, predicted);
                    Vector coefficients(model.design.cols(), Rational{0});
                    for (std::size_t k = 0; k < rank; ++k)
                        coefficients[model.pivots[k]] = fit.coordinates[a * rank + k];
                    require(model.design.apply(coefficients) == predicted,
                            "physical optimizer coefficient fibre");
                    const auto quadrance = metric.inner(residual, residual);
                    objective += quadrance;
                    auto f = Json::object();
                    f.set("axis", Json::string(axes[shape.axes[a]].id));
                    f.set("objective_represented", vec(predicted));
                    f.set("objective_residual", vec(residual));
                    f.set("particular_coefficients", vec(coefficients));
                    f.set("coefficient_kernel", mat(model.kernel));
                    f.set("residual_quadrance", rat(quadrance));
                    f.set("represented", reading(hodge, p.readout.apply(predicted), lambda, stop));
                    f.set("residual", reading(hodge, p.readout.apply(residual), lambda, stop));
                    fits.push(f);
                }
                const auto correction = numeric::sub(fit.coordinates, v0);
                require(objective == unbounded_objective + maps.gram.inner(correction, correction),
                        "physical original objective pullback identity");
                for (auto i : fit.active)
                    active.push(cardinality(i));
                result.set("state", Json::string(states[s].id));
                result.set("response", Json::string(response.id));
                result.set("model", Json::string(variants[shape.variant].id));
                result.set("preparation", cardinality(owner));
                result.set("fits", fits);
                if (response.kind == "ordinal") {
                    auto categories = Json::array();
                    for (std::size_t k = 0; k <= shape.axes.size(); ++k) {
                        Vector values(n, Rational{0});
                        for (std::size_t i = 0; i < n; ++i) {
                            const Rational lower = k ? represented[(k - 1) * n + i] : Rational{1};
                            const Rational upper =
                                k < shape.axes.size() ? represented[k * n + i] : Rational{0};
                            values[i] = lower - upper;
                            require(values[i] >= 0 && values[i] <= 1,
                                    "joint physical ordinal category range");
                        }
                        categories.push(vec(values));
                    }
                    result.set("ordinal_level_fields", categories);
                }
                result.set("joint_residual_quadrance", rat(objective));
                result.set("unconstrained_joint_quadrance", rat(unbounded_objective));
                result.set("constraint_cost", rat(objective - unbounded_objective));
                result.set("coordinates", vec(fit.coordinates));
                result.set("unconstrained_coordinates", vec(v0));
                result.set("active_rows", active);
                result.set("active_rank", cardinality(fit.active.size()));
                result.set("multipliers", vec(fit.multipliers));
                result.set("slack", vec(fit.slack));
                result.set("half_gradient_stationarity", vec(fit.stationarity));
                result.set("tested_independent_faces", exact_count(fit.tested));
                result.set("dependent_faces", exact_count(fit.dependent));
                result.set(
                    "optimality",
                    Json::string(
                        "primal feasibility, nonnegative active duals and exact stationarity"));
                shape_results[j] = result;
                shape_objectives[j] = objective;
                shape_represented[j] = represented;
            });
        auto results = Json::array(), comparisons = Json::array();
        for (const auto &result : shape_results)
            results.push(result);
        for (std::size_t s = 0; s < states.size(); ++s)
            for (std::size_t r = 0; r < responses.size(); ++r)
                for (std::size_t v = 0; v < variants.size(); ++v) {
                    if (v == 1)
                        continue;
                    const auto reduced = s * shapes.size() + r * variants.size() + v;
                    const auto full = s * shapes.size() + r * variants.size() + 1;
                    Rational displacement{0}, cross{0};
                    const auto delta =
                        numeric::sub(shape_represented[full], shape_represented[reduced]);
                    for (std::size_t a = 0; a < shapes[r * variants.size()].axes.size(); ++a) {
                        Vector d(delta.begin() + a * n, delta.begin() + (a + 1) * n);
                        const auto residual =
                            numeric::sub(observations[s][shapes[r * variants.size()].axes[a]],
                                         Vector(shape_represented[full].begin() + a * n,
                                                shape_represented[full].begin() + (a + 1) * n));
                        displacement += metric.inner(d, d);
                        cross += 2 * metric.inner(residual, d);
                    }
                    const Rational gain = shape_objectives[reduced] - shape_objectives[full];
                    require(gain >= 0 && gain == displacement + cross,
                            "nested physical objective finite cross-term identity");
                    auto comparison = Json::object();
                    comparison.set("state", Json::string(states[s].id));
                    comparison.set("response", Json::string(responses[r].id));
                    comparison.set("reduced", Json::string(variants[v].id));
                    comparison.set("joint_gain", rat(gain));
                    comparison.set("displacement_quadrance", rat(displacement));
                    comparison.set("residual_cross_term", rat(cross));
                    comparisons.push(comparison);
                }
        auto physical = Json::object();
        physical.set("declaration", query.at("physical_shape"));
        physical.set("preparations", shape_preparations);
        physical.set("families", shape_families);
        physical.set("prepared_families", cardinality(shape_families.size()));
        physical.set("results", results);
        physical.set("nested_comparisons", comparisons);
        physical.set("completed_tiles", cardinality(shape_report.completed_tiles));
        physical.set("requested_workers", cardinality(shape_report.requested_workers));
        physical.set(
            "owner",
            Json::string(
                "rfc::path::least_quadrance_feasible, exact affine algebra and CpuTileExecutor"));
        physical.set("reference_calibration",
                     Json::string("none; any top-level finite calibration describes the "
                                  "unconstrained reference objective"));
        out.set("physical_shape", physical);
    }
    if (query.has("reference_law")) {
        Vector masses;
        for (const auto &state : states)
            masses.push_back(state.mass);
        out.set("reference_calibration", scientific_detail::calibrate_finite_scores(
                                             scores, defined, masses, observed, family, alpha));
    }
    auto execution = Json::object();
    execution.set("model_carriers", cardinality(prepared.size()));
    execution.set("model_preparations",
                  cardinality(machine::checked_size_mul(prepared.size(), variants.size())));
    execution.set("hodge_preparations", cardinality(2));
    execution.set("requested_workers", cardinality(report.requested_workers));
    execution.set("worker_slots", cardinality(report.worker_slots));
    execution.set("completed_tiles", cardinality(report.completed_tiles));
    execution.set(
        "owner",
        Json::string("rfc::numeric::PreparedExactSolve and rfc::runtime::CpuTileExecutor"));
    out.set("execution", execution);
    return out;
}
} // namespace gadget::research_detail
