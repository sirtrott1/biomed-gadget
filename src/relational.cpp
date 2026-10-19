#include <gadget/build.hpp>
#include <algorithm>
#include <chrono>
#include <gadget/relational.hpp>
#include <gadget/resources.hpp>
#include <rfc/algebra/materialize.hpp>
#include <rfc/field/resolvent.hpp>
#include <rfc/runtime/work_graph.hpp>
#include <rfc/storage/rfcx.hpp>
#include <rfc/topology/integral_homology.hpp>
#include <set>

namespace gadget {
namespace {
using namespace rfc;
using Matrix = numeric::ExactMatrix;
void require(bool ok, const std::string &reason) {
    if (!ok)
        throw std::invalid_argument("relational audit: " + reason);
}
Json count(std::size_t n) { return Json::integer(static_cast<std::int64_t>(n)); }
Json qjson(Rational q) {
    q.canonicalize();
    return Json::string(numeric::str(q));
}
Json vjson(const Vector &v) {
    auto out = Json::array();
    for (const auto &q : v)
        out.push(qjson(q));
    return out;
}
template <class M> Json mjson(const M &m) {
    auto out = Json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        auto row = Json::array();
        for (std::size_t j = 0; j < m.cols(); ++j)
            row.push(qjson(Rational{m.get(i, j)}));
        out.push(row);
    }
    return out;
}
std::uint64_t key(const Json &j) {
    const auto k = j.integer();
    require(k > 0, "keys must be positive signed-64 integers");
    return static_cast<std::uint64_t>(k);
}
std::size_t index(const Json &j, std::size_t size) {
    auto n = j.integer();
    require(n >= 0 && static_cast<std::uint64_t>(n) < size, "lower index outside carrier");
    return static_cast<std::size_t>(n);
}
Rational coefficient(const Json &j) {
    auto s = j.text();
    require(s.size() <= 128 && s.find('/') == std::string::npos,
            "integer incidence coefficient required");
    auto q = exact_decimal(s);
    require(numeric::denominator(q) == 1 && q != 0, "nonzero integral incidence required");
    return q;
}
std::vector<field::ExactMetric> metrics(const core::RelationalComplex &c) {
    std::vector<field::ExactMetric> out;
    for (std::size_t k = 0; k < c.grade_count(); ++k) {
        const auto &g = c.grade(core::Grade{static_cast<std::uint16_t>(k)});
        Vector diagonal(g.cell_count(), Rational{1});
        if (k == 1)
            for (std::size_t j = 0; j < diagonal.size(); ++j) {
                auto [b, e] = g.incidence_range(j);
                diagonal[j] = Rational{e - b};
            }
        out.push_back(field::ExactMetric::diagonal(std::move(diagonal)));
    }
    return out;
}
struct Context {
    std::string id;
    algebra::MaterializedSubcomplex source;
    std::vector<field::ExactMetric> metric;
    Context(std::string name, const core::RelationalComplex &c, const algebra::SubcomplexMask &mask)
        : id(std::move(name)), source(algebra::materialize_subcomplex(c, mask)),
          metric(metrics(source.complex)) {}
    field::ExactHodgeOperator hodge(unsigned k) const {
        return {source.complex, core::Grade{static_cast<std::uint16_t>(k)}, metric};
    }
};
Json topology_reading(const Context &c) {
    auto out = Json::object(), grades = Json::array();
    for (unsigned k = 1; k < c.source.complex.grade_count(); ++k) {
        auto h = topology::integral_homology(c.source.complex,
                                             core::Grade{static_cast<std::uint16_t>(k)});
        auto g = Json::object();
        g.set("grade", count(k));
        g.set("chain_rank", count(h.chain_rank));
        g.set("cycle_rank", count(h.cycle_rank));
        g.set("boundary_rank", count(h.boundary_rank));
        g.set("free_rank", count(h.free_rank));
        g.set("cycle_basis", mjson(h.cycle_basis));
        auto torsion = Json::array();
        for (auto &t : h.torsion_invariants)
            torsion.push(qjson(Rational{t}));
        g.set("torsion", torsion);
        if (k == 1)
            out.set("patient_harmonic_dimension",
                    count(c.source.complex.grade(core::Grade{0}).cell_count() -
                          (h.chain_rank - h.cycle_rank)));
        grades.push(g);
    }
    out.set("homology", grades);
    out.set("boundary_1",
            mjson(field::materialize_boundary_matrix(c.source.complex, core::Grade{1})));
    out.set("boundary_2",
            mjson(c.source.complex.grade_count() > 2
                      ? field::materialize_boundary_matrix(c.source.complex, core::Grade{2})
                      : Matrix(c.hodge(1).size(), 0)));
    auto map = Json::array();
    for (const auto &grade : c.source.to_ambient) {
        auto a = Json::array();
        for (auto i : grade)
            a.push(count(i));
        map.push(a);
    }
    out.set("to_ambient", map);
    auto op = c.hodge(1);
    auto p = field::exact_hodge_projectors(op);
    out.set("gradient_dimension", count(p.gradient.cols()));
    out.set("curl_dimension", count(p.curl.cols()));
    out.set("harmonic_dimension", count(p.harmonic.cols()));
    out.set("harmonic_frame", mjson(p.harmonic));
    // Read every primary relation direction, including directions a lifted
    // patient field cannot excite. No biological circulation is implied.
    auto seeds = Json::array();
    for (std::size_t i = 0; i < op.size(); ++i) {
        Vector e(op.size(), Rational{0});
        e[i] = 1;
        auto a = p.project_gradient(e), b = p.project_curl(e), h = p.project_harmonic(e);
        require(numeric::add(numeric::add(a, b), h) == e, "Hodge sector reconstruction");
        require(op.metric().inner(a, b) == 0 && op.metric().inner(a, h) == 0 &&
                    op.metric().inner(b, h) == 0,
                "Hodge sector orthogonality");
        Vector upper;
        if (c.source.complex.grade_count() > 2)
            upper = field::BoundaryOperator(c.source.complex, core::Grade{2})
                        .adjoint(e, c.metric[1], c.metric[2]);
        auto coordinates = field::CompleteFieldCoordinates{
            field::BoundaryOperator(c.source.complex, core::Grade{1}).apply(e), upper, h};
        // One complete inverse certificate per context; all directions still
        // retain their three exact sector vectors below.
        if (i == 0)
            require(field::reconstruct_complete_field(c.source.complex, core::Grade{1}, c.metric,
                                                      coordinates) == e,
                    "complete field inverse");
        auto green = field::exact_green_apply(op, Rational{1}, e);
        require(green.exact && numeric::add(green.value, op.apply(green.value)) == e,
                "grade-one Green residual");
        auto s = Json::object();
        s.set("ambient_relation", count(c.source.to_ambient[1][i]));
        s.set("gradient", vjson(a));
        s.set("curl", vjson(b));
        s.set("harmonic", vjson(h));
        s.set("green", vjson(green.value));
        seeds.push(s);
    }
    out.set("relation_seeds", seeds);
    return out;
}
Json field_reading(const Context &design, const Context &outcome, const Context &control,
                   const Vector &x, const Vector &selector, const Json &declaration) {
    auto a = design.hodge(0), b = outcome.hodge(0);
    auto ga = field::exact_green_apply(a, Rational{1}, x);
    auto gb = field::exact_green_apply(b, Rational{1}, x);
    require(ga.exact && gb.exact, "patient Green solve incomplete");
    require(numeric::add(ga.value, a.apply(ga.value)) == x &&
                numeric::add(gb.value, b.apply(gb.value)) == x,
            "patient Green residual");
    const auto delta = numeric::sub(gb.value, ga.value);
    const auto forcing = numeric::sub(b.apply(ga.value), a.apply(ga.value));
    auto driven = field::exact_green_apply(b, Rational{1}, forcing);
    require(driven.exact && numeric::scale(Rational{-1}, driven.value) == delta,
            "cross-context Green identity");
    auto co = control.hodge(0);
    auto gc = field::exact_green_apply(co, Rational{1}, x);
    require(gc.exact && numeric::add(gc.value, co.apply(gc.value)) == x,
            "population-control Green residual");
    const auto split = numeric::sub(gb.value, gc.value);
    const auto split_forcing = numeric::sub(b.apply(gc.value), co.apply(gc.value));
    auto split_driven = field::exact_green_apply(b, Rational{1}, split_forcing);
    require(split_driven.exact && numeric::scale(Rational{-1}, split_driven.value) == split,
            "split Green identity");
    auto restriction = Matrix::diagonal(selector);
    const auto selected = restriction.apply(x);
    auto selected_response = field::exact_green_apply(a, Rational{1}, selected);
    require(selected_response.exact, "selected Green solve incomplete");
    const auto commutator = numeric::sub(restriction.apply(ga.value), selected_response.value);
    const auto B = field::BoundaryOperator(design.source.complex, core::Grade{1});
    const auto lifted = B.adjoint(x, design.metric[0], design.metric[1]);
    const auto patient_frame =
        field::materialize_boundary_matrix(design.source.complex, core::Grade{1}).image_frame();
    const auto represented = field::projector_apply(patient_frame, design.metric[0], x);
    const auto residual = numeric::sub(x, represented);
    require(numeric::is_zero(a.apply(residual)) &&
                design.metric[0].inner(represented, residual) == 0,
            "patient harmonic decomposition");
    const auto winding =
        field::metric_frame_dual_apply(design.hodge(1).harmonic_frame(), design.metric[1], lifted);
    require(numeric::is_zero(winding), "coboundary must annihilate harmonic cycles");
    auto out = Json::object();
    out.set("id", declaration.at("id"));
    out.set("unit", declaration.at("unit"));
    out.set("visit", declaration.at("visit"));
    out.set("input", vjson(x));
    out.set("design_response", vjson(ga.value));
    out.set("outcome_response", vjson(gb.value));
    out.set("population_control_response", vjson(gc.value));
    out.set("outcome_split_difference", vjson(split));
    out.set("outcome_split_quadrance", qjson(numeric::dot(split, split)));
    out.set("context_difference", vjson(delta));
    out.set("context_difference_quadrance", qjson(numeric::dot(delta, delta)));
    out.set("context_difference_nonzero",
            count(static_cast<std::size_t>(
                std::count_if(delta.begin(), delta.end(), [](const auto &q) { return q != 0; }))));
    out.set("selection_commutator", vjson(commutator));
    out.set("selection_commutator_quadrance", qjson(numeric::dot(commutator, commutator)));
    out.set("lifted_field", vjson(lifted));
    out.set("patient_context_component", vjson(represented));
    out.set("patient_harmonic_component", vjson(residual));
    out.set("patient_harmonic_quadrance", qjson(numeric::dot(residual, residual)));
    out.set("lifted_harmonic_pairings", vjson(winding));
    out.set("exact_residuals", Json::boolean(true));
    out.set("green_iterations_design", count(ga.iterations));
    out.set("green_iterations_outcome", count(gb.iterations));
    return out;
}
} // namespace

void run_relational_audit(const fs::path &input, const fs::path &output, std::size_t workers) {
    require(!fs::exists(output), "output already exists");
    require(workers >= 1 && workers <= 4, "workers must be 1..4");
    const auto started = std::chrono::steady_clock::now();
    auto bytes = read_text(input, 16 * 1024 * 1024);
    auto packet = Json::parse(bytes);
    require(packet.at("schema").text() == "gadget.relational.cohort_design.v1", "schema");
    require(packet.at("metric").text() == "patient_identity_relation_cardinality_face_identity",
            "metric declaration");
    require(packet.at("lambda").text() == "1", "profile requires lambda=1");
    auto artifacts = packet.at("artifacts");
    require(artifacts.size() >= 1 && artifacts.size() <= 16, "artifact count");
    for (std::size_t i = 0; i < artifacts.size(); ++i) {
        auto a = artifacts.at(i);
        fs::path path(a.at("path").text());
        require(!path.is_absolute() && !path.empty(), "artifact relative path");
        for (const auto &part : path)
            require(part != "..", "artifact path escape");
        auto resolved = fs::canonical(input.parent_path() / path),
             root = fs::canonical(input.parent_path());
        auto relative = resolved.lexically_relative(root);
        require(!relative.empty() && *relative.begin() != "..", "artifact symlink escape");
        require(sha256(read_text(resolved)) == a.at("sha256").text(), "artifact hash");
    }
    auto vertices = packet.at("vertices"), relations = packet.at("relations"),
         faces = packet.at("faces"), fields = packet.at("fields");
    require(vertices.size() > 0 && vertices.size() <= 512 && relations.size() > 0 &&
                relations.size() <= 64 && faces.size() <= 32 && fields.size() > 0 &&
                fields.size() <= 32,
            "profile size limit");
    core::RelationSignature signature;
    auto membership = signature.add_type("declared_category_membership",
                                         core::BoundaryPolicy::DeclaredExact, 1, 0);
    auto proof = signature.add_type("authored_partition_identity",
                                    core::BoundaryPolicy::DeclaredExact, 1, 0);
    core::RelationalComplexBuilder builder(signature);
    std::set<std::uint64_t> vertex_keys, relation_keys, face_keys, occurrence_keys;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        auto v = vertices.at(i);
        require(vertex_keys.insert(key(v.at("key"))).second, "duplicate vertex");
        builder.declare_witness_identity(core::IdentityKey{key(v.at("key"))});
    }
    for (std::size_t i = 0; i < relations.size(); ++i) {
        auto r = relations.at(i);
        require(relation_keys.insert(key(r.at("key"))).second, "duplicate relation");
        require(r.at("stage").text() == "design" || r.at("stage").text() == "outcome" ||
                    r.at("stage").text() == "control",
                "relation stage");
        auto terms = r.at("occurrences");
        require(terms.size() > 0 && terms.size() <= 512, "membership arity");
        std::vector<core::BoundaryOccurrenceInput> occurrences;
        for (std::size_t j = 0; j < terms.size(); ++j) {
            auto t = terms.at(j);
            auto lower = index(t.at("lower"), vertices.size());
            auto q = coefficient(t.at("coefficient"));
            require(q == 1, "membership coefficient must be +1");
            auto id = key(t.at("key"));
            require(occurrence_keys.insert(id).second, "duplicate occurrence key");
            occurrences.push_back({core::IdentityKey{key(vertices.at(lower).at("key"))},
                                   true,
                                   false,
                                   q,
                                   {},
                                   core::OccurrenceKey{id}});
        }
        builder.add_grade_one(core::RelationKey{key(r.at("key"))}, membership, occurrences);
    }
    const auto lower_source = builder.freeze();
    require(faces.size() > 0, "at least one authored identity required");
    for (std::size_t i = 0; i < faces.size(); ++i) {
        auto f = faces.at(i);
        require(face_keys.insert(key(f.at("key"))).second, "duplicate face");
        require(f.at("stage").text() == "design" || f.at("stage").text() == "outcome" ||
                    f.at("stage").text() == "control",
                "face stage");
        auto terms = f.at("occurrences");
        require(terms.size() > 0 && terms.size() <= 64, "identity arity");
        std::vector<core::HigherIncidenceInput> occurrences;
        for (std::size_t j = 0; j < terms.size(); ++j) {
            auto t = terms.at(j);
            auto lower = index(t.at("lower"), relations.size());
            require(relations.at(lower).at("stage").text() == "design" ||
                        relations.at(lower).at("stage").text() == f.at("stage").text(),
                    "identity crosses unavailable stages");
            auto q = coefficient(t.at("coefficient"));
            auto id = key(t.at("key"));
            require(occurrence_keys.insert(id).second, "duplicate occurrence key");
            occurrences.push_back({static_cast<std::uint32_t>(lower),
                                   true,
                                   q < 0,
                                   q < 0 ? -q : q,
                                   {},
                                   core::OccurrenceKey{id}});
        }
        builder.add_higher(core::Grade{2}, core::RelationKey{key(f.at("key"))}, proof, occurrences);
    }
    const auto ambient = builder.freeze(); // Validates both exact and integral chain families.
    algebra::SubcomplexMask base(lower_source), filled(ambient), all(ambient),
        control_mask(ambient);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        core::CellId cell{core::Grade{0}, static_cast<std::uint32_t>(i)};
        base.select(cell);
        filled.select(cell);
        all.select(cell);
        control_mask.select(cell);
    }
    for (std::size_t i = 0; i < relations.size(); ++i) {
        core::CellId cell{core::Grade{1}, static_cast<std::uint32_t>(i)};
        if (relations.at(i).at("stage").text() != "control")
            all.select(cell);
        if (relations.at(i).at("stage").text() != "outcome")
            control_mask.select(cell);
        if (relations.at(i).at("stage").text() == "design") {
            base.select(cell);
            filled.select(cell);
        }
    }
    for (std::size_t i = 0; i < faces.size(); ++i) {
        core::CellId cell{core::Grade{2}, static_cast<std::uint32_t>(i)};
        if (faces.at(i).at("stage").text() != "control")
            all.select(cell);
        if (faces.at(i).at("stage").text() != "outcome")
            control_mask.select(cell);
        if (faces.at(i).at("stage").text() == "design")
            filled.select(cell);
    }
    require(base.boundary_closed(lower_source) && filled.boundary_closed(ambient) &&
                all.boundary_closed(ambient) && control_mask.boundary_closed(ambient),
            "context closure");
    Context open("unfilled_design", lower_source, base), design("design", ambient, filled),
        outcome("with_outcome", ambient, all),
        control("with_population_control", ambient, control_mask);
    require(design.source.to_ambient[0] == outcome.source.to_ambient[0], "common patient carrier");
    auto glued = algebra::pushout_inside(ambient, all, control_mask);
    auto lattice = Json::object(), overlap = Json::array(), union_counts = Json::array(),
         closures = Json::array();
    for (std::size_t k = 0; k < ambient.grade_count(); ++k) {
        core::Grade grade{static_cast<std::uint16_t>(k)};
        overlap.push(count(glued.overlap.grade(grade).count()));
        union_counts.push(count(glued.union_mask.grade(grade).count()));
        for (std::size_t i = 0; i < ambient.grade(grade).cell_count(); ++i)
            require(glued.overlap.contains({grade, static_cast<std::uint32_t>(i)}) ==
                        filled.contains({grade, static_cast<std::uint32_t>(i)}),
                    "source-context intersection differs from shared design");
    }
    for (std::size_t i = 0; i < faces.size(); ++i) {
        algebra::SubcomplexMask closure(ambient);
        closure.select({core::Grade{2}, static_cast<std::uint32_t>(i)});
        closure.close_boundaries(ambient);
        auto row = Json::array();
        for (std::size_t k = 0; k < ambient.grade_count(); ++k)
            row.push(count(closure.grade(core::Grade{static_cast<std::uint16_t>(k)}).count()));
        closures.push(row);
    }
    lattice.set("overlap_counts", overlap);
    lattice.set("union_counts", union_counts);
    lattice.set("face_boundary_closure_counts", closures);
    lattice.set("glued_patient_identities",
                count(glued.glued.complex.grade(core::Grade{0}).cell_count()));
    auto selector_input = packet.at("selector");
    auto selector = Vector(vertices.size(), Rational{0});
    require(selector_input.size() == vertices.size(), "selector carrier");
    for (std::size_t i = 0; i < selector.size(); ++i) {
        auto q = exact_decimal(selector_input.at(i).text());
        require(q == 0 || q == 1, "selector expects 0/1");
        selector[i] = q;
    }
    auto resources = CpuResources::local();
    require(workers <= resources->coordinator().capacity(), "worker admission");
    auto top = Json::object();
    top.set(open.id, topology_reading(open));
    top.set(design.id, topology_reading(design));
    top.set(outcome.id, topology_reading(outcome));
    top.set(control.id, topology_reading(control));
    std::vector<runtime::WorkJob> jobs;
    std::set<std::string> field_ids;
    for (std::size_t i = 0; i < fields.size(); ++i) {
        auto f = fields.at(i);
        require(field_ids.insert(f.at("id").text()).second, "duplicate field ID");
        auto v = f.at("values");
        require(v.size() == vertices.size(), "full supported field required");
        Vector x;
        x.reserve(v.size());
        for (std::size_t j = 0; j < v.size(); ++j) {
            require(v.at(j).text().size() <= 128, "decimal size");
            x.push_back(exact_decimal(v.at(j).text()));
        }
        jobs.push_back(runtime::WorkJob::make<Json>(
            {i + 1, 1}, {1, 8 * 1024 * 1024, 8 * 1024 * 1024}, {},
            [&, x = std::move(x), f](const auto &, auto &, std::stop_token stop) {
                if (stop.stop_requested())
                    throw runtime::ExecutionCancelled{};
                return field_reading(design, outcome, control, x, selector, f);
            }));
    }
    runtime::WorkGraphLimits limits;
    limits.in_flight = workers;
    auto graph = runtime::WorkGraph(std::move(jobs), limits);
    auto completed = graph.submit(resources->coordinator());
    completed.wait();
    auto observations = Json::array();
    for (std::size_t i = 0; i < fields.size(); ++i)
        observations.push(*completed.completion({i + 1, 1}).result<Json>());
    auto result = Json::object();
    result.set("schema", Json::string("gadget.relational.audit.v1"));
    result.set("input_sha256", Json::string(sha256(bytes)));
    result.set("declaration", packet);
    result.set("rexgraph_headers_sha256", Json::string(gadget::build::rexgraph_sha256()));
    result.set("gadget_sources_sha256", Json::string(gadget::build::source_sha256()));
    result.set("topology", top);
    result.set("native_lattice", lattice);
    result.set("fields", observations);
    auto jobs_receipt = Json::array();
    for (std::size_t i = 0; i < fields.size(); ++i) {
        const auto receipt = completed.completion({i + 1, 1}).receipt();
        auto entry = Json::object();
        entry.set("field", fields.at(i).at("id"));
        entry.set("workers", count(receipt.demand.workers));
        entry.set("execution_ns", count(receipt.execution_ns));
        entry.set("queue_ns", count(receipt.queue_ns));
        jobs_receipt.push(entry);
    }
    auto runtime = Json::object();
    runtime.set("workers", count(workers));
    runtime.set("jobs", jobs_receipt);
    runtime.set(
        "seconds",
        Json::number(
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()));
    result.set("runtime", runtime);
    fs::create_directories(output);
    storage::write_rfcx(ambient, (output / "retained.rfcx").string());
    result.set("native_sha256", Json::string(sha256(read_text(output / "retained.rfcx"))));
    write_text(output / "results.json", result.dump() + "\n");
}
} // namespace gadget
