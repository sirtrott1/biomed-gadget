#include "scientific_structure.hpp"
#include <algorithm>
#include <gadget/structure.hpp>
#include <rfc/path/realization.hpp>
#include <rfc/topology/integral_homology.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc;
using numeric::ExactMatrix;
ResultValue structural_matrix(ResultFields &fields, const ExactMatrix &matrix,
    std::uint64_t input, std::uint64_t output) {
    auto value = fields.matrix(matrix, input, output);
    value.set("display", ResultValue::string("matrix_triplets"));
    return value;
}
ResultValue integer_vector(const std::vector<numeric::Integer> &values,
                           ResultFields &fields, std::uint64_t axis) {
    Vector rational;
    for (const auto &value : values) rational.emplace_back(value);
    return fields.vector(rational, axis);
}
template <class Matrix> ResultValue integer_matrix(const Matrix &matrix,
    ResultFields &fields, std::uint64_t input, std::uint64_t output) {
    ExactMatrix rational(matrix.rows(), matrix.cols());
    for (std::size_t i = 0; i < matrix.rows(); ++i)
        for (std::size_t j = 0; j < matrix.cols(); ++j)
            if (matrix.get(i, j) != 0) rational.set(i, j, Rational{matrix.get(i, j)});
    return structural_matrix(fields, rational, input, output);
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
        out.push_back(exact_rational(text(values.at(i))));
    }
    return out;
}
ResultValue coordinates(const core::RelationalComplex &c, core::Grade grade, ResultFields &fields) {
    auto out = ResultValue::array();
    const auto &arena = c.grade(grade);
    for (std::size_t i = 0; i < arena.cell_count(); ++i) {
        auto item = ResultValue::object();
        item.set("index", ResultValue::count(i));
        item.set("key",
                 ResultValue::string(std::to_string(
                     grade.value == 0 ? arena.identities[i] : arena.relation_identities[i].value)));
        item.set("scale", fields.scalar(arena.cell_scales[i]));
        out.push(item);
    }
    return out;
}
} // namespace

ResultValue scientific_detail::chain_fields(const NativeStudy &study, const Json &query,
    const Vector &input, ResultFields &fields, const std::vector<std::uint64_t> &axes,
    const runtime::ExecutionLease &lease, bool validate_only) {
    members(query, {"operation", "context", "grade", "coefficient_domain"}, {"field", "values"});
    require(query.has("field") != query.has("values"), "chain requires values or a source field");
    require(text(query.at("operation")) == "chain_certificate", "chain certificate operation");
    const auto &c = study.source();
    admit_source(c);
    const auto grade = core::Grade{
        static_cast<std::uint16_t>(count(query.at("grade"), c.grade_count() - 1, "chain grade"))};
    const auto domain = text(query.at("coefficient_domain"));
    require(domain == "rational" || domain == "retained_lattice", "chain coefficient domain");
    const auto &arena = c.grade(grade);
    require(input.size() == arena.cell_count(), "complete chain coordinates required");
    lease.check_cancelled();
    std::size_t bits = 0;
    for (const auto &value : input) bits = machine::checked_size_add(bits,
        storage::native_execution_detail::integer_bits(value));
    runtime::admit_execution({1, input.size(), 0, 0, bits}, lease.policy());
    Vector x = input;
    if (domain == "retained_lattice")
        for (std::size_t i = 0; i < x.size(); ++i) {
            require(numeric::denominator(x[i]) == 1, "lattice chains require integral coordinates");
            x[i] *= arena.cell_scales[i];
        }
    if (validate_only)
        return ResultValue::object();
    check_chain_lattice(c);
    lease.check_cancelled();
    const auto axis = axes.at(grade.value);
    const auto lower = grade.value ? axes.at(grade.value - 1) : fields.axis(0);
    const auto upper = grade.value + 1u < axes.size() ? axes.at(grade.value + 1u) : fields.axis(0);
    const auto lattice_axis = fields.axis(arena.cell_count());
    const auto lower_lattice = fields.axis(grade.value ? c.grade(core::Grade{
        static_cast<std::uint16_t>(grade.value - 1)}).cell_count() : 0);
    const auto upper_lattice = fields.axis(grade.value + 1u < c.grade_count() ? c.grade(core::Grade{
        static_cast<std::uint16_t>(grade.value + 1u)}).cell_count() : 0);
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
    auto out = ResultValue::object(), facts = ResultValue::object(), sectors = ResultValue::object(),
         energy = ResultValue::object();
    out.set("grade", ResultValue::count(grade.value));
    out.set("coefficient_domain", ResultValue::string(domain));
    out.set("coordinates", coordinates(c, grade, fields));
    out.set("input_coordinates", fields.vector(input, (domain == "retained_lattice" ? lattice_axis : axis)));
    out.set("rational_coordinates", fields.vector(x, axis));
    out.set("boundary", fields.vector(bx, lower));
    out.set("boundary_matrix", structural_matrix(fields, B, axis, lower));
    out.set("upper_boundary_matrix", structural_matrix(fields, U, upper, axis));
    out.set("boundary_preimage", rational_boundary ? fields.vector(preimage, upper) : ResultValue{});
    sectors.set("gradient", fields.vector(gradient, axis));
    sectors.set("curl", fields.vector(curl, axis));
    sectors.set("harmonic", fields.vector(harmonic, axis));
    for (const auto &[name, v] : std::vector<std::pair<std::string, Vector>>{
             {"gradient", gradient}, {"curl", curl}, {"harmonic", harmonic}, {"total", x}})
        energy.set(name, fields.scalar(metric.inner(v, v)));
    out.set("sectors", sectors);
    out.set("quadrance", energy);
    auto reconstruction = ResultValue::object();
    reconstruction.set("lower_boundary", fields.vector(complete.lower_boundary, lower));
    reconstruction.set("upper_adjoint", fields.vector(complete.upper_adjoint, upper));
    reconstruction.set("harmonic", fields.vector(complete.harmonic, axis));
    reconstruction.set("reconstructed", fields.vector(field::reconstruct_complete_field(
                                            c, grade, study.tower().metrics(), complete), axis));
    out.set("complete_field_coordinates", reconstruction);
    facts.set("zero_chain", ResultValue::boolean(zero(x)));
    facts.set("cycle", ResultValue::boolean(cycle));
    facts.set("rational_boundary", ResultValue::boolean(rational_boundary));
    facts.set("nonzero_rational_homology_class", ResultValue::boolean(cycle && !rational_boundary));
    facts.set("hodge_reconstruction_verified", ResultValue::boolean(true));
    facts.set("coordinated_lattice_verified", ResultValue::boolean(true));
    auto occurrences = ResultValue::array(), sums = ResultValue::array();
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
            auto o = ResultValue::object(), aspects = ResultValue::object();
            o.set("key", ResultValue::string(std::to_string(arena.occurrence_identities[p].value)));
            o.set("upper", ResultValue::count(arena.upper_index[p]));
            o.set("lower", ResultValue::count(arena.lower_index[p]));
            o.set("exists", ResultValue::boolean(exists));
            o.set("lattice_coefficient", fields.scalar(arena.lattice_coefficient[p]));
            o.set("rational_coefficient", fields.scalar(coefficient));
            o.set("contribution", fields.scalar(contributions[p]));
            for (std::size_t a = 0; a < c.signature().aspects().size(); ++a) {
                const auto id = static_cast<composite::AspectId>(a);
                aspects.set(c.signature().aspects().at(id).name,
                            ResultValue::boolean(arena.aspects.get(id, p)));
            }
            o.set("aspects", aspects);
            occurrences.push(o);
        }
        const auto occurrence_axis = fields.axis(arena.incidence_count());
        ExactMatrix contribution_map(arena.incidence_count(), x.size());
        for (std::size_t p = 0; p < arena.incidence_count(); ++p)
            if (arena.aspects.get(c.signature().existence_aspect(), p))
                contribution_map.set(p, arena.upper_index[p], c.incidence_coefficient(grade, p));
        out.set("occurrence_contribution_map", structural_matrix(fields, contribution_map, axis, occurrence_axis));
        out.set("occurrence_contributions", fields.vector(contributions, occurrence_axis));
        out.set("occurrence_gluing", structural_matrix(fields, glue, occurrence_axis, lower));
        require(glue.apply(contributions) == bx,
                "retained occurrence gluing reproduces native boundary");
        const auto lower_coordinates =
            coordinates(c, core::Grade{static_cast<std::uint16_t>(grade.value - 1)}, fields);
        for (std::size_t i = 0; i < B.rows(); ++i) {
            auto s = lower_coordinates.at(i);
            const bool cancels = nonzero_at_lower[i] > 0 && bx[i] == 0;
            cancelled += cancels;
            s.set("sum", fields.scalar(bx[i]));
            s.set("nonzero_contributions", ResultValue::count(nonzero_at_lower[i]));
            s.set("cancellation", ResultValue::boolean(cancels));
            sums.push(s);
        }
    }
    out.set("occurrences", occurrences);
    out.set("lower_stalk_sums", sums);
    facts.set("retained_occurrences", ResultValue::count(arena.incidence_count()));
    facts.set("active_occurrences", ResultValue::count(active));
    facts.set("nonzero_occurrence_contributions", ResultValue::count(nonzero_contributions));
    facts.set("cancelled_lower_coordinates", ResultValue::count(cancelled));
    Vector lattice_values(x.size());
    std::vector<numeric::Integer> lattice_chain;
    bool belongs = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
        lattice_values[i] = x[i] / arena.cell_scales[i];
        belongs = belongs && numeric::denominator(lattice_values[i]) == 1;
        lattice_chain.push_back(numeric::numerator(lattice_values[i]));
    }
    out.set("lattice_coordinates", belongs ? fields.vector(lattice_values, lattice_axis) : ResultValue{});
    const auto homology = topology::integral_homology(c, grade);
    auto integral = ResultValue::object();
    const auto cycle_axis = fields.axis(homology.cycle_rank);
    const auto smith_cycles = fields.axis(homology.cycle_rank);
    const auto smith_upper = fields.axis(homology.boundaries_in_cycle_coordinates.cols());
    integral.set("coefficient_domain", ResultValue::string("integer"));
    integral.set("chain_rank", ResultValue::count(homology.chain_rank));
    integral.set("cycle_rank", ResultValue::count(homology.cycle_rank));
    integral.set("boundary_rank", ResultValue::count(homology.boundary_rank));
    integral.set("free_rank", ResultValue::count(homology.free_rank));
    integral.set("torsion_invariants", integer_vector(homology.torsion_invariants, fields, fields.axis(homology.torsion_invariants.size())));
    integral.set("cycle_basis", integer_matrix(homology.cycle_basis, fields, cycle_axis, lattice_axis));
    integral.set("boundaries_in_cycle_coordinates",
                 integer_matrix(homology.boundaries_in_cycle_coordinates, fields, upper_lattice, cycle_axis));
    integral.set("smith_diagonal", integer_matrix(homology.boundary_smith.diagonal, fields, smith_upper, smith_cycles));
    integral.set("smith_left", integer_matrix(homology.boundary_smith.left, fields, cycle_axis, smith_cycles));
    integral.set("smith_right", integer_matrix(homology.boundary_smith.right, fields, smith_upper, upper_lattice));
    facts.set("belongs_to_retained_lattice", ResultValue::boolean(belongs));
    bool integral_boundary = false;
    if (belongs) {
        const auto certificate = topology::integral_chain_certificate(c, grade, lattice_chain);
        integral_boundary = certificate.boundary;
        require(certificate.cycle == cycle, "coordinated cycle membership");
        facts.set("integral_boundary", ResultValue::boolean(certificate.boundary));
        facts.set("nonzero_integral_homology_class",
                  ResultValue::boolean(certificate.represents_nonzero_homology_class()));
        integral.set("boundary_preimage",
                     certificate.boundary ? integer_vector(certificate.boundary_preimage, fields, upper_lattice) : ResultValue{});
        integral.set("boundary",
                     integer_vector(topology::integral_boundary_matrix(c, grade).apply(lattice_chain), fields, lower_lattice));
    } else {
        facts.set("integral_boundary", ResultValue{});
        facts.set("nonzero_integral_homology_class", ResultValue{});
        integral.set("boundary_preimage", ResultValue{});
    }
    const auto upper_kernel = U.kernel_frame();
    const auto parameter_axis = fields.axis(upper_kernel.cols());
    out.set("boundary_preimage_directions", structural_matrix(fields, upper_kernel, parameter_axis, upper));
    if (grade.value + 1u < c.grade_count()) {
        const auto upper_homology = topology::integral_homology(c,
            core::Grade{static_cast<std::uint16_t>(grade.value + 1u)});
        integral.set("boundary_preimage_directions", integer_matrix(upper_homology.cycle_basis,
            fields, fields.axis(upper_homology.cycle_rank), upper_lattice));
    } else {
        integral.set("boundary_preimage_directions", structural_matrix(fields, ExactMatrix(0, 0),
            fields.axis(0), upper_lattice));
    }
    const auto scales = [&](core::Grade g, std::uint64_t from, std::uint64_t to) {
        const auto &a = c.grade(g);
        ExactMatrix map(a.cell_count(), a.cell_count());
        for (std::size_t i = 0; i < a.cell_count(); ++i) map.set(i, i, a.cell_scales[i]);
        return structural_matrix(fields, map, from, to);
    };
    integral.set("lattice_realization", scales(grade, lattice_axis, axis));
    integral.set("lower_boundary_matrix", integer_matrix(topology::integral_boundary_matrix(c, grade),
        fields, lattice_axis, lower_lattice));
    if (grade.value) integral.set("lower_lattice_realization", scales(
        core::Grade{static_cast<std::uint16_t>(grade.value - 1)}, lower_lattice, lower));
    if (grade.value + 1u < c.grade_count()) {
        const core::Grade g{static_cast<std::uint16_t>(grade.value + 1u)};
        integral.set("upper_boundary_matrix", integer_matrix(topology::integral_boundary_matrix(c, g),
            fields, upper_lattice, lattice_axis));
        integral.set("upper_lattice_realization", scales(g, upper_lattice, upper));
    } else integral.set("upper_boundary_matrix", structural_matrix(fields, ExactMatrix(x.size(), 0),
        upper_lattice, lattice_axis));
    out.set("integral_homology", integral);
    out.set("facts", facts);
    auto reader = ResultValue::object();
    reader.set("facts", facts);
    reader.set("rational_coordinates", out.at("rational_coordinates"));
    reader.set("boundary", out.at("boundary"));
    reader.set("sectors", sectors);
    reader.set("quadrance", energy);
    std::string explanation = "The declared grade " + std::to_string(grade.value) + " chain has " +
                              std::to_string(nonzero_contributions) +
                              " nonzero retained incidence contributions. ";
    if (cancelled)
        explanation += "Oppositely signed contributions cancel at " + std::to_string(cancelled) +
                       " lower coordinates; their occurrence identities remain in the source. ";
    explanation += cycle ? "Its boundary is the zero chain. " : "Its boundary is nonzero. ";
    if (cycle)
        explanation += rational_boundary
                           ? "It is a boundary over Q; a preimage in the higher grade is recorded. "
                           : "It represents a nonzero homology class over Q. ";
    if (belongs && cycle && rational_boundary && !integral_boundary)
        explanation += "It is not a boundary in the retained integral lattice. ";
    explanation +=
        "The gradient, curl and harmonic fields reconstruct the input in the declared metric.";
    reader.set("explanation", ResultValue::string(explanation));
    lease.check_cancelled();
    out.set("reader", reader);
    return out;
}

ResultValue scientific_detail::path_fields(const NativeStudy &study, const Json &query,
    ResultFields &fields, const std::vector<std::uint64_t> &axes,
    const runtime::ExecutionLease &lease, bool validate_only) {
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
    auto records = ResultValue::array();
    const auto axis = axes.at(grade.value), lower = axes.at(grade.value - 1);
    const auto family_axis = fields.axis(entries.size());
    lease.check_cancelled();
    std::size_t total_steps = 0;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto p = entries.at(i);
        members(p, {"id", "steps"});
        require(ids.insert(text(p.at("id"))).second, "duplicate path identity");
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
        std::size_t bits = 0;
        for (const auto &step : steps) bits = machine::checked_size_add(bits,
            storage::native_execution_detail::integer_bits(step.coefficient));
        runtime::admit_execution({1, steps.size(), 0, 0, bits}, lease.policy());
        lease.check_cancelled();
        paths.push_back(path::realize(c, grade, std::move(steps)));
        auto r = ResultValue::object(), ordered = ResultValue::array();
        const auto &realized = paths.back();
        const auto step_axis = fields.axis(realized.ordered_steps.size());
        Vector coefficients;
        ExactMatrix step_realization(realized.field.size(), realized.ordered_steps.size());
        for (std::size_t j = 0; j < realized.ordered_steps.size(); ++j) {
            const auto &step = realized.ordered_steps[j];
            coefficients.push_back(step.coefficient);
            step_realization.set(step.relation, j, Rational{1});
            auto declaration = ResultValue::object();
            declaration.set("index", ResultValue::count(step.relation));
            declaration.set("coefficient", fields.scalar(step.coefficient));
            ordered.push(std::move(declaration));
        }
        r.set("id", ResultValue::string(text(p.at("id"))));
        r.set("steps", std::move(ordered));
        r.set("ordered_coefficients", fields.vector(coefficients, step_axis));
        r.set("step_realization", structural_matrix(fields, step_realization, step_axis, axis));
        r.set("field", fields.vector(realized.field, axis));
        r.set("boundary", fields.vector(realized.boundary, lower));
        records.push(r);
    }
    if (validate_only)
        return ResultValue::object();
    check_chain_lattice(c);
    const auto P = path::realization_matrix(paths),
               endpoints = boundary(c, grade.value).multiply(P);
    const auto gram = path::inherited_gram(paths, study.tower().metrics().at(grade.value));
    const auto endpoint_gram =
        field::metric_gram(endpoints, study.tower().metrics().at(grade.value - 1), endpoints);
    const auto hodge = study.hodge(grade);
    const auto projectors = field::exact_hodge_projectors(hodge);
    auto pairs = ResultValue::array();
    for (std::size_t i = 0; i < paths.size(); ++i)
        for (std::size_t j = i + 1; j < paths.size(); ++j) {
            lease.check_cancelled();
            const auto difference = numeric::sub(paths[i].field, paths[j].field);
            const auto g = projectors.project_gradient(difference),
                       u = projectors.project_curl(difference),
                       h = projectors.project_harmonic(difference);
            require(numeric::add(numeric::add(g, u), h) == difference,
                    "path difference reconstruction");
            auto p = ResultValue::object(), sectors = ResultValue::object();
            p.set("left", ResultValue::string(text(entries.at(i).at("id"))));
            p.set("right", ResultValue::string(text(entries.at(j).at("id"))));
            p.set("equal_fields", ResultValue::boolean(zero(difference)));
            p.set("equal_boundaries", ResultValue::boolean(paths[i].boundary == paths[j].boundary));
            p.set("difference", fields.vector(difference, axis));
            sectors.set("gradient", fields.vector(g, axis));
            sectors.set("curl", fields.vector(u, axis));
            sectors.set("harmonic", fields.vector(h, axis));
            p.set("sectors", sectors);
            pairs.push(p);
        }
    auto out = ResultValue::object();
    out.set("paths", records);
    out.set("coordinates", coordinates(c, grade, fields));
    out.set("realization", structural_matrix(fields, P, family_axis, axis));
    out.set("endpoint_realization", structural_matrix(fields, endpoints, family_axis, lower));
    out.set("shared_field_form", structural_matrix(fields, gram, family_axis, family_axis));
    out.set("endpoint_form", structural_matrix(fields, endpoint_gram, family_axis, family_axis));
    out.set("realization_rank", ResultValue::count(P.rank()));
    out.set("endpoint_rank", ResultValue::count(endpoints.rank()));
    out.set("directions_lost_by_endpoint_view", ResultValue::count(P.rank() - endpoints.rank()));
    out.set("pairs", pairs);
    const auto realization_kernel = P.kernel_frame(), endpoint_kernel = endpoints.kernel_frame();
    out.set("realization_kernel", structural_matrix(fields, realization_kernel,
        fields.axis(realization_kernel.cols()), family_axis));
    out.set("endpoint_kernel", structural_matrix(fields, endpoint_kernel,
        fields.axis(endpoint_kernel.cols()), family_axis));
    out.set("boundary_matrix", structural_matrix(fields, boundary(c, grade.value), axis, lower));
    lease.check_cancelled();
    out.set("scope",
            ResultValue::string("Ordered declared relational realizations; algebraic fields and boundary "
                         "views do not certify simple, directed or biologically identified walks"));
    return out;
}
namespace scientific_detail {
Vector declared_chain_values(const Json &query, std::size_t size) {
    return read_values(query.at("values"), size);
}
} // namespace scientific_detail
namespace {
struct StandaloneFields {
    NativeRecord source;
    ResultFields fields;
    std::vector<std::uint64_t> axes;
    explicit StandaloneFields(const NativeStudy &study)
        : source(NativeRecord::create_from_endpoint(NativeRecordValue(NativeRecordValue::Object{}),
            new_native_record_namespace(), {}, storage::TemporalEndpoint(study.source()))), fields(source) {
        for (std::size_t k = 0; k < study.source().grade_count(); ++k) {
            const auto axis = fields.source_grade(k);
            axes.push_back(axis);
            auto metric = ResultValue::object();
            metric.set("kind", ResultValue::string("declared_metric"));
            metric.set("grade", ResultValue::count(k));
            metric.set("metric", ResultFields::stage_matrix(study.tower().metrics().at(k).materialize(),
                axis, axis, storage::TemporalTensorRole::Metric));
            fields.bind(std::move(metric));
        }
    }
};
} // namespace
void scientific_detail::validate_chain_values(const NativeStudy &study, const Json &query,
    const Vector &input) {
    StandaloneFields publication(study);
    (void)chain_fields(study, query, input, publication.fields, publication.axes,
        source_execution_lease("scientific-chain-admission"), true);
}
Json analyze_chain(const NativeStudy &study, const Json &query, bool validate_only) {
    StandaloneFields publication(study);
    const auto grade = count(query.at("grade"), study.source().grade_count() - 1, "chain grade");
    const auto input = scientific_detail::declared_chain_values(query,
        study.source().grade(core::Grade{static_cast<std::uint16_t>(grade)}).cell_count());
    auto out = scientific_detail::chain_fields(study, query, input, publication.fields,
        publication.axes, source_execution_lease("scientific-chain"), validate_only);
    return validate_only ? Json::object() : result_fields_view(publication.fields.finish(std::move(out)));
}
Json compare_realized_paths(const NativeStudy &study, const Json &query, bool validate_only) {
    StandaloneFields publication(study);
    auto out = scientific_detail::path_fields(study, query, publication.fields,
        publication.axes, source_execution_lease("scientific-paths"), validate_only);
    return validate_only ? Json::object() : result_fields_view(publication.fields.finish(std::move(out)));
}
} // namespace gadget
