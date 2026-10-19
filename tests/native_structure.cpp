#include "../examples/scientific/fixture.hpp"
#include "../src/scientific_structure.hpp"
#include <gadget/assay.hpp>
#include <iostream>
#include <rfc/storage/rfcx.hpp>

using namespace gadget;
using namespace rfc::storage;
using Vector = rfc::numeric::RationalVector;
using Matrix = rfc::numeric::ExactMatrix;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
const TemporalTensor *tensor(const NativeRecord &owner, NativeRecord reference) {
    const auto key = reference.at("tensor").integer();
    for (const auto &field : owner.endpoint()->tensors) if (rfc::numeric::Integer(std::to_string(field.key)) == key) return &field;
    throw std::runtime_error("native structural tensor absent");
}
Vector values(const NativeRecord &owner, NativeRecord reference) {
    const auto &field = (*tensor(owner, reference)).tensor;
    check(field.support().count() == field.domain().logical_size(), "complete structural support lost");
    Vector out;
    for (std::size_t i = 0; i < field.domain().logical_size(); ++i) out.push_back(*field.state(i));
    return out;
}
Matrix matrix(const NativeRecord &owner, NativeRecord reference) {
    const auto &field = (*tensor(owner, reference));
    check(field.tensor.support().count() == field.tensor.domain().logical_size(), "complete native matrix support lost");
    return temporal_tensor_matrix(field);
}
Matrix original_boundary(const rfc::core::RelationalComplex &source, std::size_t grade) {
    const auto size = source.grade({static_cast<std::uint16_t>(grade)}).cell_count();
    const auto lower = grade ? source.grade({static_cast<std::uint16_t>(grade - 1)}).cell_count() : 0;
    Matrix out(lower, size);
    if (!grade) return out;
    const auto &arena = source.grade({static_cast<std::uint16_t>(grade)});
    for (std::size_t i = 0; i < arena.incidence_count(); ++i) {
        if (!arena.aspects.get(source.signature().existence_aspect(), i)) continue;
        const auto coefficient = source.incidence_coefficient({static_cast<std::uint16_t>(grade)}, i);
        out.set(arena.lower_index[i], arena.upper_index[i],
            out.get(arena.lower_index[i], arena.upper_index[i]) + coefficient);
    }
    return out;
}
void chain(const NativeRecord &owner) {
    const auto out = owner.at("view");
    const auto grade = out.at("grade").integer().get_ui();
    const auto &source = owner.endpoint()->source;
    const auto &arena = source.grade({static_cast<std::uint16_t>(grade)});
    const auto x = values(owner, out.at("rational_coordinates"));
    const auto b = matrix(owner, out.at("boundary_matrix"));
    const auto u = matrix(owner, out.at("upper_boundary_matrix"));
    check(b == original_boundary(source, grade), "boundary differs from all original occurrences");
    check(b.apply(x) == values(owner, out.at("boundary")), "native boundary witness differs");
    auto expected_upper = grade + 1 < source.grade_count() ? original_boundary(source, grade + 1)
                                                         : Matrix(x.size(), 0);
    check(u == expected_upper, "upper boundary source differs");
    const auto &chain_field = (*tensor(owner, out.at("rational_coordinates")));
    const auto &boundary_field = (*tensor(owner, out.at("boundary")));
    check(chain_field.output != boundary_field.output && (*tensor(owner, out.at("boundary_matrix"))).input == chain_field.output &&
          (*tensor(owner, out.at("boundary_matrix"))).output == boundary_field.output,
          "independent grades identified by dimension");
    const auto sectors = out.at("sectors");
    check(rfc::numeric::add(rfc::numeric::add(values(owner, sectors.at("gradient")),
        values(owner, sectors.at("curl"))), values(owner, sectors.at("harmonic"))) == x,
        "native sectors do not reconstruct");
    check(values(owner, out.at("complete_field_coordinates").at("reconstructed")) == x,
          "complete native coordinates do not reconstruct");
    const auto k = matrix(owner, out.at("boundary_preimage_directions"));
    check(k.rows() == u.cols() && k.cols() == u.cols() - u.rank() && u.multiply(k).nonzeros() == 0,
          "complete rational preimage directions missing");
    if (out.at("facts").at("rational_boundary").boolean())
        check(u.apply(values(owner, out.at("boundary_preimage"))) == x, "rational preimage differs");
    const auto integral = out.at("integral_homology");
    const auto ib = matrix(owner, integral.at("lower_boundary_matrix"));
    const auto iu = matrix(owner, integral.at("upper_boundary_matrix"));
    const auto scale = matrix(owner, integral.at("lattice_realization"));
    check(ib.multiply(iu).nonzeros() == 0, "integral chain law differs");
    check((*tensor(owner, integral.at("lattice_realization"))).input != chain_field.output &&
          (*tensor(owner, integral.at("lattice_realization"))).output == chain_field.output,
          "integer and rational carriers identified");
    if (grade) check(b.multiply(scale) == matrix(owner, integral.at("lower_lattice_realization")).multiply(ib),
                     "lower coordinated lattice identity differs");
    if (grade + 1 < source.grade_count()) check(u.multiply(matrix(owner, integral.at("upper_lattice_realization"))) == scale.multiply(iu),
                                               "upper coordinated lattice identity differs");
    const auto cycles = matrix(owner, integral.at("cycle_basis"));
    const auto presentation = matrix(owner, integral.at("boundaries_in_cycle_coordinates"));
    check(ib.multiply(cycles).nonzeros() == 0 && cycles.cols() == ib.cols() - ib.rank() &&
          cycles.multiply(presentation) == iu, "complete retained integral presentation differs");
    const auto &left = (*tensor(owner, integral.at("smith_left")));
    const auto &right = (*tensor(owner, integral.at("smith_right")));
    const auto &diagonal = (*tensor(owner, integral.at("smith_diagonal")));
    check(left.input != left.output && right.input != right.output && diagonal.input == right.input &&
          diagonal.output == left.output, "Smith coordinate carriers identified by dimension");
    check(temporal_tensor_matrix(left).multiply(presentation).multiply(temporal_tensor_matrix(right)) ==
          temporal_tensor_matrix(diagonal), "complete Smith factorization differs");
    for (const auto *transform : {&left, &right}) {
        const auto inverse = temporal_tensor_matrix(*transform).inverse();
        for (std::size_t i = 0; i < inverse.rows(); ++i) for (std::size_t j = 0; j < inverse.cols(); ++j)
            check(rfc::numeric::denominator(inverse.get(i, j)) == 1, "Smith transform lost unimodularity");
    }
    const auto ik = matrix(owner, integral.at("boundary_preimage_directions"));
    check(ik.rows() == iu.cols() && ik.cols() == iu.cols() - iu.rank() && iu.multiply(ik).nonzeros() == 0,
          "complete integral preimage directions missing");
    for (std::size_t i = 0; i < ik.rows(); ++i) for (std::size_t j = 0; j < ik.cols(); ++j)
        check(rfc::numeric::denominator(ik.get(i, j)) == 1, "integral direction coefficient changed domain");
    if (out.at("facts").at("belongs_to_retained_lattice").boolean()) {
        check(scale.apply(values(owner, out.at("lattice_coordinates"))) == x, "retained lattice membership differs");
        if (out.at("facts").at("integral_boundary").boolean()) check(iu.apply(values(owner,
            integral.at("boundary_preimage"))) == values(owner, out.at("lattice_coordinates")),
            "integral preimage differs");
    }
    const auto occurrences = out.at("occurrences");
    check(occurrences.size() == arena.incidence_count(), "ordered occurrences missing");
    Vector contributions;
    for (std::size_t i = 0; i < occurrences.size(); ++i) {
        const auto original = occurrences.at(i);
        const bool active = arena.aspects.get(source.signature().existence_aspect(), i);
        check(original.at("key").text() == std::to_string(arena.occurrence_identities[i].value) &&
              original.at("exists").boolean() == active && original.at("upper").integer() == arena.upper_index[i] &&
              original.at("lower").integer() == arena.lower_index[i], "original ordered occurrence declaration differs");
        const auto expected = active ? source.incidence_coefficient({static_cast<std::uint16_t>(grade)}, i) * x[arena.upper_index[i]] : Rational{0};
        check(values(owner, original.at("contribution")) == Vector{expected}, "original signed contribution differs");
        contributions.push_back(expected);
    }
    if (grade) {
        check(matrix(owner, out.at("occurrence_contribution_map")).apply(x) == contributions &&
              contributions == values(owner, out.at("occurrence_contributions")), "occurrence realization differs");
        check(matrix(owner, out.at("occurrence_gluing")).apply(contributions) == b.apply(x), "original gluing differs");
    }
}
void paths(const NativeRecord &owner) {
    const auto out = owner.at("view"), family = out.at("paths");
    const auto p = matrix(owner, out.at("realization"));
    const auto endpoints = matrix(owner, out.at("endpoint_realization"));
    const auto boundary = matrix(owner, out.at("boundary_matrix"));
    check(boundary.multiply(p) == endpoints, "path endpoints differ from boundary realization");
    const auto &map = (*tensor(owner, out.at("realization")));
    check(map.input != map.output && (*tensor(owner, out.at("endpoint_realization"))).input == map.input &&
          (*tensor(owner, out.at("endpoint_realization"))).output != map.output,
          "independent path family and source carriers identified");
    for (std::size_t i = 0; i < family.size(); ++i) {
        const auto path = family.at(i), steps = path.at("steps");
        const auto coefficients = values(owner, path.at("ordered_coefficients"));
        Vector expected(p.rows(), Rational{0});
        check(coefficients.size() == steps.size(), "ordered path step support missing");
        for (std::size_t j = 0; j < steps.size(); ++j) {
            check(values(owner, steps.at(j).at("coefficient")) == Vector{coefficients[j]}, "ordered step coefficient differs");
            expected.at(steps.at(j).at("index").integer().get_ui()) += coefficients[j];
        }
        check(expected == [&] { Vector column(p.rows()); for (std::size_t row = 0; row < p.rows(); ++row) column[row] = p.get(row, i); return column; }() && expected == values(owner, path.at("field")) &&
              matrix(owner, path.at("step_realization")).apply(coefficients) == expected,
              "ordered repeated steps differ from native realization");
        check(boundary.apply(expected) == values(owner, path.at("boundary")), "path boundary differs");
    }
    const auto metric = [&](std::uint64_t axis) {
        for (const auto &field : owner.endpoint()->tensors)
            if (field.role == TemporalTensorRole::Metric && field.input == axis && field.output == axis)
                return temporal_tensor_matrix(field);
        throw std::runtime_error("declared path metric absent");
    };
    check(p.transpose().multiply(metric(*map.output)).multiply(p) == matrix(owner, out.at("shared_field_form")),
          "shared field form lost its original metric");
    const auto &endpoint_map = (*tensor(owner, out.at("endpoint_realization")));
    check(endpoints.transpose().multiply(metric(*endpoint_map.output)).multiply(endpoints) == matrix(owner, out.at("endpoint_form")),
          "endpoint form lost its independent declared metric");
    const auto pairs = out.at("pairs");
    std::size_t pair = 0;
    for (std::size_t i = 0; i < family.size(); ++i) for (std::size_t j = i + 1; j < family.size(); ++j) {
        const auto item = pairs.at(pair++);
        const auto difference = rfc::numeric::sub(values(owner, family.at(i).at("field")), values(owner, family.at(j).at("field")));
        check(item.at("left").text() == family.at(i).at("id").text() && item.at("right").text() == family.at(j).at("id").text() &&
              values(owner, item.at("difference")) == difference, "ordered path pair difference changed");
        const auto sectors = item.at("sectors");
        check(rfc::numeric::add(rfc::numeric::add(values(owner, sectors.at("gradient")), values(owner, sectors.at("curl"))),
            values(owner, sectors.at("harmonic"))) == difference, "native path sectors do not reconstruct");
        check(item.at("equal_fields").boolean() == rfc::numeric::is_zero(difference) &&
              item.at("equal_boundaries").boolean() == rfc::numeric::is_zero(boundary.apply(difference)), "path pair facts changed");
    }
    check(pair == pairs.size(), "complete ordered path pairs missing");
    for (const auto *name : {"realization_kernel", "endpoint_kernel"}) {
        const auto base = std::string(name) == "realization_kernel" ? p : endpoints;
        const auto kernel = matrix(owner, out.at(name));
        check(base.multiply(kernel).nonzeros() == 0 && kernel.cols() == base.cols() - base.rank(),
              "complete path family directions missing");
        check((*tensor(owner, out.at(name))).input != map.input && (*tensor(owner, out.at(name))).output == map.input,
              "path ambiguity parameter carrier identified");
    }
}
void verify(const NativeRecord &publication, const NativeRecord &original, const NativeRecord &query) {
    check(has_native_scientific_publisher(query.at("operation").text()), "structural producer not native");
    const auto owner = publication.at("result").at("fields");
    check(owner.at("schema").text() == "gadget.scientific.fields.v1", "structural native owner absent");
    check(NativeRecord::decode(publication.at("result").at("source_study_publication").bytes()).encode() == original.encode(),
          "whole original study publication changed");
    if (query.at("operation").text() == "chain_certificate") chain(owner); else paths(owner);
    const auto binding = owner.at("bindings").at(0);
    check(NativeRecord::decode(binding.at("original_study_publication").bytes()).encode() == original.encode(),
          "context lost original byte owner");
    if (query.has("field")) {
        const auto &field = (*tensor(owner, binding.at("source_field")));
        const auto key = binding.at("source_tensor").integer();
        const TemporalTensor *actual = nullptr;
        for (const auto &item : original.endpoint()->tensors) if (rfc::numeric::Integer(std::to_string(item.key)) == key) actual = &item;
        check(actual && temporal_endpoint_detail::tensor_equal(actual->tensor, field.tensor) && actual->output && rfc::numeric::Integer(std::to_string(*actual->output)) == binding.at("source_carrier").integer(),
              "authentic original source pin differs");
        const ScientificStudy admitted(original);
        const auto input_pin = admitted.field_accession(query.at("field").text());
        check(binding.at("source_namespace").text() == std::to_string(input_pin.source_namespace()) &&
              binding.at("source_cutoff").integer() == rfc::numeric::Integer(std::to_string(input_pin.sequence())) &&
              temporal_endpoint_detail::tensor_equal(*actual, *input_pin.tensor_declaration()),
              "authentic source namespace, cutoff or tensor differs");
        check(matrix(owner, binding.at("context_realization")).apply(values(owner, binding.at("source_field"))) ==
              values(owner, owner.at("view").at("input_coordinates")), "native source context realization differs");
    }
}
void portable(const NativeRecord &publication, const fs::path &path) {
    publication.write(path);
    const auto restored = NativeRecord::read(path);
    check(restored.encode() == publication.encode(), "portable native bytes changed");
    const auto replayed = replay_scientific_record(restored);
    check(replayed.endpoint()->tensors.size() == publication.endpoint()->tensors.size(), "portable replay field count changed");
    for (std::size_t i = 0; i < publication.endpoint()->tensors.size(); ++i)
        check(temporal_endpoint_detail::tensor_equal(publication.endpoint()->tensors[i], replayed.endpoint()->tensors[i]),
              "portable replay changed native carriers, coefficients or support");
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 5) {
            const auto study = std::string(argv[4]) == "-" ? std::make_shared<const ScientificStudy>(fs::path(argv[1])) :
                AssayWorkspace(argv[1]).study(argv[4]);
            const auto original = study->native_publication();
            ScientificSession session(study);
            const auto queries = admit_record_view(Json::parse(read_text(argv[2])));
            fs::create_directories(argv[3]);
            for (const auto &name : queries.names()) {
                const auto query = queries.at(name), publication = session.execute_record(query);
                verify(publication, original, query);
                portable(publication, fs::path(argv[3]) / (name + ".rfcxs"));
            }
        } else {
            if (argc != 2) throw std::invalid_argument("usage: native_structure OWNED_PARENT");
            const fs::path root = fs::path(argv[1]) / ("run-" + std::to_string(new_native_record_namespace()));
            fs::create_directories(root);
            example::scientific_study(root / "source");
            const auto study = std::make_shared<const ScientificStudy>(root / "source/study.json");
            const auto original = study->native_publication();
            ScientificSession session(study);
            std::vector<NativeRecord> publications;
            for (const auto *context : {"baseline", "followup"}) {
                for (const auto *field : {"G", "E"}) {
                    auto q = Json::parse(R"({"operation":"chain_certificate","context":"baseline","grade":0,"coefficient_domain":"rational","field":"baseline.G"})");
                    q.set("context", Json::string(context)); q.set("field", Json::string(std::string(context) + "." + field));
                    const auto query = admit_record_view(q), publication = session.execute_record(query);
                    verify(publication, original, query); publications.push_back(publication);
                }
                for (const auto *input : {"[\"0\",\"0\",\"0\"]", "[\"1\",\"0\",\"2\"]"}) {
                    auto q = Json::parse(R"({"operation":"chain_certificate","context":"baseline","grade":1,"coefficient_domain":"rational","values":[]})");
                    q.set("context", Json::string(context));
                    q.set("values", Json::parse(std::string(context) == "baseline" ?
                        (std::string(input) == "[\"0\",\"0\",\"0\"]" ? "[\"0\",\"0\"]" : "[\"1\",\"2\"]") : input));
                    for (const auto *domain : {"rational", "retained_lattice"}) {
                        q.set("coefficient_domain", Json::string(domain));
                        const auto query = admit_record_view(q), publication = session.execute_record(query);
                        verify(publication, original, query); publications.push_back(publication);
                    }
                }
                auto path_query = Json::parse(R"({"operation":"realized_path_comparison","context":"baseline","grade":1,"paths":[{"id":"cancel","steps":[{"index":2,"coefficient":"1"},{"index":2,"coefficient":"-1"}]},{"id":"zero","steps":[]},{"id":"closed","steps":[{"index":2,"coefficient":"3/2"}]}]})");
                path_query.set("context", Json::string(context));
                if (std::string(context) == "baseline") {
                    path_query.at("paths").at(0).at("steps").at(0).set("index", Json::integer(1));
                    path_query.at("paths").at(0).at("steps").at(1).set("index", Json::integer(1));
                    path_query.at("paths").at(2).at("steps").at(0).set("index", Json::integer(1));
                }
                const auto query = admit_record_view(path_query);
                const auto publication = session.execute_record(query);
                verify(publication, original, query); publications.push_back(publication);
            }
            {
                using namespace scientific_detail;
                std::vector<rfc::field::ExactMetric> metrics;
                for (std::size_t i = 0; i < original.endpoint()->source.grade_count(); ++i)
                    metrics.push_back(rfc::field::ExactMetric::identity(original.endpoint()->source.grade(
                        {static_cast<std::uint16_t>(i)}).cell_count()));
                const NativeStudy native(original.endpoint()->source, std::move(metrics));
                ResultFields fields(original);
                std::vector<std::uint64_t> axes;
                for (std::size_t i = 0; i < native.source().grade_count(); ++i) axes.push_back(fields.source_grade(i));
                const auto request = Json::parse(R"({"operation":"chain_certificate","context":"global","grade":1,"coefficient_domain":"rational","values":["1","0","2"]})");
                const auto baseline_lease = source_execution_lease("structural-test-policy");
                auto policy = baseline_lease.policy();
                policy.limits.rational_slots = 2;
                bool resource_refused = false;
                try { (void)chain_fields(native, request, Vector{1, 0, 2}, fields, axes, rfc::runtime::ExecutionLease(policy)); }
                catch (const rfc::runtime::ExecutionLimitExceeded &) { resource_refused = true; }
                check(resource_refused, "declared structural coefficient policy was ignored");
                auto cancelled = std::make_shared<std::atomic_bool>(true);
                bool cancellation_refused = false;
                try { (void)chain_fields(native, request, Vector{1, 0, 2}, fields, axes,
                    rfc::runtime::ExecutionLease(baseline_lease.policy(), cancelled)); }
                catch (const rfc::runtime::ExecutionCancelled &) { cancellation_refused = true; }
                check(cancellation_refused, "structural cancellation was ignored");
            }
            bool refused = false;
            try { session.execute_record(admit_record_view(Json::parse(R"({"operation":"chain_certificate","context":"baseline","grade":0,"coefficient_domain":"rational","field":"baseline.partial"})"))); }
            catch (const std::invalid_argument &) { refused = true; }
            check(refused, "absent chain input was synthesized as zero");
            rfc::numeric::Integer large = rfc::numeric::Integer{1} << 65537;
            auto large_query = Json::parse(R"({"operation":"chain_certificate","context":"followup","grade":1,"coefficient_domain":"rational","values":[]})");
            auto coordinates = Json::array();
            coordinates.push(Json::string(large.get_str())); coordinates.push(Json::string("0")); coordinates.push(Json::string("0"));
            large_query.set("values", coordinates);
            const auto large_request = admit_record_view(large_query), large_publication = session.execute_record(large_request);
            verify(large_publication, original, large_request); publications.push_back(large_publication);
            fs::rename(root / "source", root / "moved-source");
            for (std::size_t i = 0; i < publications.size(); ++i) portable(publications[i], root / (std::to_string(i) + ".rfcxs"));
            check(study->native_publication().encode() == original.encode(), "execution mutated original study");
        }
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
