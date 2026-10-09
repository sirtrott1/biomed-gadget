#include "scientific_detail.hpp"
#include <algorithm>
#include <gadget/structure.hpp>
#include <rfc/path/realization.hpp>
#include <rfc/topology/integral_homology.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc;
using numeric::ExactMatrix;
Json vector_json(const Vector &v) {
    auto out = Json::array();
    for (const auto &x : v) {
        out.push(Json::string(x.get_str()));
    }
    return out;
}
Json integers(const std::vector<numeric::Integer> &v) {
    auto out = Json::array();
    for (const auto &x : v)
        out.push(Json::string(x.get_str()));
    return out;
}
template <class Matrix> Json matrix_json(const Matrix &m) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(m.rows()));
    out.set("columns", cardinality(m.cols()));
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (m.get(i, j) != 0) {
                auto e = Json::array();
                e.push(cardinality(i));
                e.push(cardinality(j));
                e.push(Json::string(m.get(i, j).get_str()));
                entries.push(e);
            }
    out.set("entries", entries);
    return out;
}
bool zero(const Vector &v) {
    return std::all_of(v.begin(), v.end(), [](const auto &x) { return x == 0; });
}
ExactMatrix boundary(const core::RelationalComplex &c, std::size_t k) {
    if (k == 0)
        return ExactMatrix(0, c.grade(core::Grade{0}).cell_count());
    if (k == c.grade_count())
        return ExactMatrix(c.grade(c.max_grade()).cell_count(), 0);
    return field::materialize_boundary_matrix(c, core::Grade{static_cast<std::uint16_t>(k)});
}
void admit_source(const core::RelationalComplex &c) {
    std::size_t cells = 0, occurrences = 0;
    for (std::size_t k = 0; k < c.grade_count(); ++k) {
        cells = machine::checked_size_add(cells, c.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count());
        occurrences = machine::checked_size_add(occurrences, c.grade(core::Grade{static_cast<std::uint16_t>(k)}).incidence_count());
    }
    (void)numeric::IndexDomain({cells});
    (void)numeric::IndexDomain({occurrences});
}
void check_chain_lattice(const core::RelationalComplex &c) {
    for (std::size_t k = 1; k < c.grade_count(); ++k) {
        const auto g = core::Grade{static_cast<std::uint16_t>(k)};
        const auto B = boundary(c, k);
        const auto I = topology::integral_boundary_matrix(c, g);
        ExactMatrix upper(B.cols(), B.cols()), lower(B.rows(), B.rows()),
            lattice(B.rows(), B.cols());
        for (std::size_t i = 0; i < B.cols(); ++i)
            upper.set(i, i, c.grade(g).cell_scales[i]);
        for (std::size_t i = 0; i < B.rows(); ++i) {
            lower.set(i, i, c.grade(core::Grade{static_cast<std::uint16_t>(k - 1)}).cell_scales[i]);
            for (std::size_t j = 0; j < B.cols(); ++j)
                if (I.get(i, j) != 0)
                    lattice.set(i, j, Rational{I.get(i, j)});
        }
        require(B.multiply(upper) == lower.multiply(lattice),
                "retained coordinated lattice identity");
        if (k >= 2) {
            require(boundary(c, k - 1).multiply(B).nonzeros() == 0, "rational chain law");
            const auto prior = topology::integral_boundary_matrix(
                c, core::Grade{static_cast<std::uint16_t>(k - 1)});
            require(prior * I == numeric::IntegerMatrix(prior.rows(), I.cols()),
                    "integral chain law");
        }
    }
}
Vector read_values(const Json &values, std::size_t n) {
    require(values.size() == n, "complete chain coordinates required");
    Vector out;
    for (std::size_t i = 0; i < n; ++i) {
        out.push_back(exact_rational(text(values.at(i), 4096)));
        require(storage::native_execution_detail::integer_bits(out.back()) <= 65536,
                "chain coefficient integer-bit budget");
    }
    return out;
}
Json coordinates(const core::RelationalComplex &c, core::Grade grade) {
    auto out = Json::array();
    const auto &arena = c.grade(grade);
    for (std::size_t i = 0; i < arena.cell_count(); ++i) {
        auto item = Json::object();
        item.set("index", cardinality(i));
        item.set("key",
                 Json::string(std::to_string(
                     grade.value == 0 ? arena.identities[i] : arena.relation_identities[i].value)));
        item.set("scale", Json::string(arena.cell_scales[i].get_str()));
        out.push(item);
    }
    return out;
}
} // namespace

Json analyze_chain(const NativeStudy &study, const Json &query, bool validate_only) {
    members(query, {"operation", "context", "grade", "coefficient_domain", "values"});
    require(text(query.at("operation")) == "chain_certificate", "chain certificate operation");
    const auto &c = study.source();
    admit_source(c);
    const auto grade = core::Grade{
        static_cast<std::uint16_t>(count(query.at("grade"), c.grade_count() - 1, "chain grade"))};
    const auto domain = text(query.at("coefficient_domain"));
    require(domain == "rational" || domain == "retained_lattice", "chain coefficient domain");
    const auto &arena = c.grade(grade);
    const auto input = read_values(query.at("values"), arena.cell_count());
    Vector x = input;
    if (domain == "retained_lattice")
        for (std::size_t i = 0; i < x.size(); ++i) {
            require(numeric::denominator(x[i]) == 1, "lattice chains require integral coordinates");
            x[i] *= arena.cell_scales[i];
        }
    if (validate_only)
        return Json::object();
    check_chain_lattice(c);
    const auto B = boundary(c, grade.value), U = boundary(c, grade.value + 1);
    const auto bx = B.apply(x);
    const auto hodge = study.hodge(grade);
    const auto projectors = field::exact_hodge_projectors(hodge);
    const auto gradient = projectors.project_gradient(x), curl = projectors.project_curl(x),
               harmonic = projectors.project_harmonic(x);
    require(numeric::add(numeric::add(gradient, curl), harmonic) == x, "Hodge reconstruction");
    const auto &metric = hodge.metric();
    require(metric.inner(gradient, curl) == 0 && metric.inner(gradient, harmonic) == 0 &&
                metric.inner(curl, harmonic) == 0,
            "Hodge sector orthogonality");
    const auto complete = field::complete_field_coordinates(c, grade, study.tower().metrics(), x);
    require(field::reconstruct_complete_field(c, grade, study.tower().metrics(), complete) == x,
            "complete field reconstruction");
    const bool cycle = zero(bx), rational_boundary = cycle && zero(harmonic);
    Vector preimage(U.cols(), Rational{0});
    if (rational_boundary && !zero(x)) {
        const auto response = field::generalized_hodge_inverse_apply(hodge, x);
        preimage = field::metric_adjoint_apply(U, study.tower().metrics().at(grade.value + 1),
                                               metric, response);
    }
    if (rational_boundary)
        require(U.apply(preimage) == x, "rational boundary preimage");
    auto out = Json::object(), facts = Json::object(), sectors = Json::object(),
         energy = Json::object();
    out.set("grade", cardinality(grade.value));
    out.set("coefficient_domain", Json::string(domain));
    out.set("coordinates", coordinates(c, grade));
    out.set("input_coordinates", vector_json(input));
    out.set("rational_coordinates", vector_json(x));
    out.set("boundary", vector_json(bx));
    out.set("boundary_matrix", matrix_json(B));
    out.set("upper_boundary_matrix", matrix_json(U));
    out.set("boundary_preimage", rational_boundary ? vector_json(preimage) : Json{});
    sectors.set("gradient", vector_json(gradient));
    sectors.set("curl", vector_json(curl));
    sectors.set("harmonic", vector_json(harmonic));
    for (const auto &[name, v] : std::vector<std::pair<std::string, Vector>>{
             {"gradient", gradient}, {"curl", curl}, {"harmonic", harmonic}, {"total", x}})
        energy.set(name, Json::string(metric.inner(v, v).get_str()));
    out.set("sectors", sectors);
    out.set("quadrance", energy);
    auto reconstruction = Json::object();
    reconstruction.set("lower_boundary", vector_json(complete.lower_boundary));
    reconstruction.set("upper_adjoint", vector_json(complete.upper_adjoint));
    reconstruction.set("harmonic", vector_json(complete.harmonic));
    reconstruction.set("reconstructed", vector_json(field::reconstruct_complete_field(
                                            c, grade, study.tower().metrics(), complete)));
    out.set("complete_field_coordinates", reconstruction);
    facts.set("zero_chain", Json::boolean(zero(x)));
    facts.set("cycle", Json::boolean(cycle));
    facts.set("rational_boundary", Json::boolean(rational_boundary));
    facts.set("nonzero_rational_homology_class", Json::boolean(cycle && !rational_boundary));
    facts.set("hodge_reconstruction_verified", Json::boolean(true));
    facts.set("coordinated_lattice_verified", Json::boolean(true));
    auto occurrences = Json::array(), sums = Json::array();
    std::size_t cancelled = 0, active = 0, nonzero_contributions = 0;
    if (grade.value > 0) {
        ExactMatrix glue(B.rows(), arena.incidence_count());
        Vector contributions(arena.incidence_count(), Rational{0});
        std::vector<std::size_t> nonzero_at_lower(B.rows());
        for (std::size_t p = 0; p < arena.incidence_count(); ++p) {
            const bool exists = arena.aspects.get(c.signature().existence_aspect(), p);
            active += exists;
            const auto coefficient = c.incidence_coefficient(grade, p);
            contributions[p] = exists ? coefficient * x[arena.upper_index[p]] : Rational{0};
            if (contributions[p] != 0) {
                ++nonzero_contributions;
                ++nonzero_at_lower[arena.lower_index[p]];
            }
            glue.set(arena.lower_index[p], p, Rational{1});
            auto o = Json::object(), aspects = Json::object();
            o.set("key", Json::string(std::to_string(arena.occurrence_identities[p].value)));
            o.set("upper", cardinality(arena.upper_index[p]));
            o.set("lower", cardinality(arena.lower_index[p]));
            o.set("exists", Json::boolean(exists));
            o.set("lattice_coefficient", Json::string(arena.lattice_coefficient[p].get_str()));
            o.set("rational_coefficient", Json::string(coefficient.get_str()));
            o.set("contribution", Json::string(contributions[p].get_str()));
            for (std::size_t a = 0; a < c.signature().aspects().size(); ++a) {
                const auto id = static_cast<composite::AspectId>(a);
                aspects.set(c.signature().aspects().at(id).name,
                            Json::boolean(arena.aspects.get(id, p)));
            }
            o.set("aspects", aspects);
            occurrences.push(o);
        }
        require(glue.apply(contributions) == bx,
                "retained occurrence gluing reproduces native boundary");
        const auto lower_coordinates =
            coordinates(c, core::Grade{static_cast<std::uint16_t>(grade.value - 1)});
        for (std::size_t i = 0; i < B.rows(); ++i) {
            auto s = lower_coordinates.at(i);
            const bool cancels = nonzero_at_lower[i] > 0 && bx[i] == 0;
            cancelled += cancels;
            s.set("sum", Json::string(bx[i].get_str()));
            s.set("nonzero_contributions", cardinality(nonzero_at_lower[i]));
            s.set("cancellation", Json::boolean(cancels));
            sums.push(s);
        }
    }
    out.set("occurrences", occurrences);
    out.set("lower_stalk_sums", sums);
    facts.set("retained_occurrences", cardinality(arena.incidence_count()));
    facts.set("active_occurrences", cardinality(active));
    facts.set("nonzero_occurrence_contributions", cardinality(nonzero_contributions));
    facts.set("cancelled_lower_coordinates", cardinality(cancelled));
    Vector lattice_values(x.size());
    std::vector<numeric::Integer> lattice_chain;
    bool belongs = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
        lattice_values[i] = x[i] / arena.cell_scales[i];
        belongs = belongs && numeric::denominator(lattice_values[i]) == 1;
        lattice_chain.push_back(numeric::numerator(lattice_values[i]));
    }
    out.set("lattice_coordinates", belongs ? vector_json(lattice_values) : Json{});
    const auto homology = topology::integral_homology(c, grade);
    auto integral = Json::object();
    integral.set("free_rank", cardinality(homology.free_rank));
    integral.set("torsion_invariants", integers(homology.torsion_invariants));
    integral.set("cycle_basis", matrix_json(homology.cycle_basis));
    integral.set("boundaries_in_cycle_coordinates",
                 matrix_json(homology.boundaries_in_cycle_coordinates));
    integral.set("smith_diagonal", matrix_json(homology.boundary_smith.diagonal));
    integral.set("smith_left", matrix_json(homology.boundary_smith.left));
    integral.set("smith_right", matrix_json(homology.boundary_smith.right));
    facts.set("belongs_to_retained_lattice", Json::boolean(belongs));
    if (belongs) {
        const auto certificate = topology::integral_chain_certificate(c, grade, lattice_chain);
        require(certificate.cycle == cycle, "coordinated cycle membership");
        facts.set("integral_boundary", Json::boolean(certificate.boundary));
        facts.set("nonzero_integral_homology_class",
                  Json::boolean(certificate.represents_nonzero_homology_class()));
        integral.set("boundary_preimage",
                     certificate.boundary ? integers(certificate.boundary_preimage) : Json{});
        integral.set("boundary",
                     integers(topology::integral_boundary_matrix(c, grade).apply(lattice_chain)));
    } else {
        facts.set("integral_boundary", Json{});
        facts.set("nonzero_integral_homology_class", Json{});
        integral.set("boundary_preimage", Json{});
    }
    out.set("integral_homology", integral);
    out.set("facts", facts);
    auto reader = Json::object();
    reader.set("facts", facts);
    reader.set("rational_coordinates", vector_json(x));
    reader.set("boundary", vector_json(bx));
    reader.set("sectors", sectors);
    reader.set("quadrance", energy);
    std::string explanation = "The declared grade-" + std::to_string(grade.value) + " chain has " +
                              std::to_string(nonzero_contributions) +
                              " nonzero retained incidence contributions. ";
    if (cancelled)
        explanation += "Oppositely signed contributions cancel at " + std::to_string(cancelled) +
                       " lower coordinates; their occurrence identities remain in the source. ";
    explanation += cycle ? "Its boundary is the zero chain. " : "Its boundary is nonzero. ";
    if (cycle)
        explanation += rational_boundary
                           ? "It is a boundary over Q; a higher-chain preimage is recorded. "
                           : "It represents a nonzero homology class over Q. ";
    if (belongs && cycle && rational_boundary && !facts.at("integral_boundary").boolean())
        explanation += "It is not a boundary in the retained integral lattice. ";
    explanation +=
        "The gradient, curl and harmonic fields reconstruct the input in the declared metric.";
    reader.set("explanation", Json::string(explanation));
    out.set("reader", reader);
    return out;
}

Json compare_realized_paths(const NativeStudy &study, const Json &query, bool validate_only) {
    members(query, {"operation", "context", "grade", "paths"});
    require(text(query.at("operation")) == "realized_path_comparison", "path comparison operation");
    const auto &c = study.source();
    admit_source(c);
    const auto grade = core::Grade{
        static_cast<std::uint16_t>(count(query.at("grade"), c.grade_count() - 1, "path grade"))};
    require(grade.value > 0, "path realizations need a positive grade");
    const auto entries = query.at("paths");
    require(entries.size() > 0, "path family requires paths");
    std::vector<path::RealizedPath> paths;
    std::set<std::string> ids;
    auto records = Json::array();
    std::size_t total_steps = 0;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto p = entries.at(i);
        members(p, {"id", "steps"});
        require(ids.insert(text(p.at("id"), 128)).second, "duplicate path identity");
        std::vector<path::TraversalStep> steps;
        total_steps = machine::checked_size_add(total_steps, p.at("steps").size());
        for (std::size_t j = 0; j < p.at("steps").size(); ++j) {
            const auto s = p.at("steps").at(j);
            members(s, {"index", "coefficient"});
            require(c.grade(grade).cell_count() > 0, "empty path carrier");
            steps.push_back({static_cast<std::uint32_t>(count(
                                 s.at("index"), c.grade(grade).cell_count() - 1, "path relation")),
                             exact_rational(text(s.at("coefficient")))});
        }
        paths.push_back(path::realize(c, grade, std::move(steps)));
        auto r = Json::parse(p.dump());
        r.set("field", vector_json(paths.back().field));
        r.set("boundary", vector_json(paths.back().boundary));
        records.push(r);
    }
    if (validate_only)
        return Json::object();
    check_chain_lattice(c);
    const auto P = path::realization_matrix(paths),
               endpoints = boundary(c, grade.value).multiply(P);
    const auto gram = path::inherited_gram(paths, study.tower().metrics().at(grade.value));
    const auto endpoint_gram =
        field::metric_gram(endpoints, study.tower().metrics().at(grade.value - 1), endpoints);
    const auto hodge = study.hodge(grade);
    const auto projectors = field::exact_hodge_projectors(hodge);
    auto pairs = Json::array();
    for (std::size_t i = 0; i < paths.size(); ++i)
        for (std::size_t j = i + 1; j < paths.size(); ++j) {
            const auto difference = numeric::sub(paths[i].field, paths[j].field);
            const auto g = projectors.project_gradient(difference),
                       u = projectors.project_curl(difference),
                       h = projectors.project_harmonic(difference);
            require(numeric::add(numeric::add(g, u), h) == difference,
                    "path difference reconstruction");
            auto p = Json::object(), sectors = Json::object();
            p.set("left", entries.at(i).at("id"));
            p.set("right", entries.at(j).at("id"));
            p.set("equal_fields", Json::boolean(zero(difference)));
            p.set("equal_boundaries", Json::boolean(paths[i].boundary == paths[j].boundary));
            p.set("difference", vector_json(difference));
            sectors.set("gradient", vector_json(g));
            sectors.set("curl", vector_json(u));
            sectors.set("harmonic", vector_json(h));
            p.set("sectors", sectors);
            pairs.push(p);
        }
    auto out = Json::object();
    out.set("paths", records);
    out.set("coordinates", coordinates(c, grade));
    out.set("realization", matrix_json(P));
    out.set("endpoint_realization", matrix_json(endpoints));
    out.set("shared_field_form", matrix_json(gram));
    out.set("endpoint_form", matrix_json(endpoint_gram));
    out.set("realization_rank", cardinality(P.rank()));
    out.set("endpoint_rank", cardinality(endpoints.rank()));
    out.set("directions_lost_by_endpoint_view", cardinality(P.rank() - endpoints.rank()));
    out.set("pairs", pairs);
    out.set("scope",
            Json::string("Ordered declared relational realizations; algebraic fields and boundary "
                         "views do not certify simple, directed or biologically identified walks"));
    return out;
}
} // namespace gadget
