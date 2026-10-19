#include "../src/scientific_result_fields.hpp"
#include <gadget/statistics/execution.hpp>
#include <gadget/snapshot.hpp>
#include <iostream>
#include <rfc/storage/rfcx.hpp>

namespace {
using namespace gadget;
using gadget::statistics::Execution;
using namespace rfc;
using numeric::ExactMatrix;
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
const storage::TemporalTensor &field(const NativeRecord &owner, const NativeRecord &reference) {
    return assay_detail::retained_tensor(owner, reference.at("tensor"));
}
Vector column(const ExactMatrix &value, std::size_t index) {
    Vector out(value.rows());
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = value.get(i, index);
    return out;
}
Vector vector(const NativeRecord &owner, const NativeRecord &reference) {
    const auto &t = field(owner, reference);
    check(t.output && t.tensor.support().count() == t.tensor.domain().logical_size(), "defined field lost its carrier or supported zero");
    Vector out(t.tensor.domain().logical_size());
    for (std::size_t i = 0; i < out.size(); ++i) out[i] = *t.tensor.state(i);
    return out;
}
ExactMatrix matrix(const NativeRecord &owner, const NativeRecord &reference) {
    const auto &t = field(owner, reference);
    check(t.input && t.output && t.tensor.support().count() == t.tensor.domain().logical_size(), "declared matrix lost full support or carrier axes");
    return storage::temporal_tensor_matrix(t);
}
NativeRecord publication(const NativeRecord &reference) {
    check(reference.at("path").size() == 0, "analysis archive unexpectedly selects another owner");
    return NativeRecord::decode(reference.at("publication").bytes());
}
void context(const NativeRecord &owner, const Workspace &original) {
    const auto view = owner.at("view");
    Execution actual(original, view.at("context").text(), view.at("realization").text());
    check(storage::encode_rfcx(owner.endpoint()->source).bytes == storage::encode_rfcx(actual.complex()).bytes, "context source or original ordered occurrences changed");
    const auto b = matrix(owner, view.at("boundary"));
    check(b == field::boundary_matrix(actual.complex(), {1}), "actual retained boundary differs");
    const auto m0 = matrix(owner, view.at("gene_metric")), m1 = matrix(owner, view.at("relation_metric"));
    check(m0 == actual.metrics()[0].materialize() && m1 == actual.metrics()[1].materialize(), "actual metrics were omitted or coerced");
    const auto gram = matrix(owner, view.at("relation_gram"));
    const auto exact_gram = b.transpose().multiply(b);
    for (std::size_t i = 0; i < gram.rows(); ++i) for (std::size_t j = 0; j < gram.cols(); ++j)
        check(gram.get(i,j) == Rational(exact_gram.get(i,j).get_d()), "Gram is not the actual lifted binary64 coefficient");
    const auto gene_indices = vector(owner, view.at("source_genes")), relation_indices = vector(owner, view.at("source_relations"));
    for (std::size_t i = 0; i < gene_indices.size(); ++i) check(gene_indices[i] == i, "ambient gene identity map is incomplete");
    for (std::size_t i = 0; i < relation_indices.size(); ++i) check(relation_indices[i] == actual.source_relations.at(i), "original relation identity map is incomplete");
    check(view.at("channel_carrier").integer() != view.at("gene_carrier").integer() &&
              view.at("gene_carrier").integer() != view.at("relation_carrier").integer(), "equal dimensions collapsed independent channel/grade carriers");
    const auto pulls = matrix(owner, view.at("gene_pulls")), coefficients = matrix(owner, view.at("gene_adjoint"));
    check(pulls == b.transpose(), "gene sweep did not retain all actual pulled columns");
    for (std::size_t gene = 0; gene < pulls.cols(); ++gene)
        check(column(coefficients, gene) == actual.metrics()[1].solve(column(pulls, gene)), "gene adjoint column does not reconstruct");
    for (std::size_t p = 0; p < view.at("moments").size(); ++p) {
        const auto moment = view.at("moments").at(p);
        const auto values = matrix(owner, moment.at("matrix")), cosine = matrix(owner, moment.at("cosine"));
        for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) {
            const auto expected = p == 0 ? actual.metrics()[0].inner(vector(owner, view.at("source_fields").at(i)), vector(owner, view.at("source_fields").at(j))) :
                numeric::dot(vector(owner, view.at("pulled_fields").at(i)), vector(owner, moment.at("relation_coefficients").at(j)));
            check(values.get(i,j) == expected, "moment lost complete ordered channel pairings");
            const auto denominator = std::sqrt(values.get(i,i).get_d() * values.get(j,j).get_d());
            check(cosine.get(i,j) == Rational(denominator > 0 ? expected.get_d()/denominator : 0), "computed normalization differs from its source field");
        }
    }
    const field::ExactHodgeOperator hodge(owner.endpoint()->source, {0}, actual.metrics());
    std::optional<std::uint64_t> shared_system;
    for (std::size_t i = 0; i < view.at("green_solves").size(); ++i) {
        const auto solve = view.at("green_solves").at(i);
        const auto source = vector(owner, solve.at("source")), pulled = vector(owner, solve.at("pulled")), solved = vector(owner, solve.at("solved"));
        const auto response = vector(owner, solve.at("response")), residual = vector(owner, solve.at("residual"));
        check(pulled == b.transpose().apply(source), "Green pull differs from actual carrier action");
        check(response == numeric::sub(source, b.apply(solved)), "Green response does not reconstruct from retained solve");
        check(residual == numeric::sub(numeric::add(response, hodge.apply(response)), source), "Green residual was reconstructed from a display");
        const auto system = matrix(owner, solve.at("system"));
        const auto key = solve.at("system").at("tensor").integer().get_ui();
        if (shared_system) check(*shared_system == key, "same declared context/lambda system stored repeatedly");
        shared_system = key;
        for (std::size_t r = 0; r < system.rows(); ++r) for (std::size_t c = 0; c < system.cols(); ++c)
            check(system.get(r,c) == Rational(m1.get(r,c).get_d() + gram.get(r,c).get_d()), "Green assembled system differs");
        const auto diagnostic = solve.at("diagnostic");
        check(diagnostic.at("absolute_residual").rational() == Rational(gadget::statistics::norm(residual)) &&
                  diagnostic.at("relative_residual").rational() <= diagnostic.at("tolerance").rational(), "Green residual quantity or bound is not retained exactly");
    }
    for (std::size_t i = 0; i < view.at("exact_checks").size(); ++i) {
        const auto row = view.at("exact_checks").at(i);
        const auto gene = row.at("gene").integer().get_ui();
        Vector basis(original.genes.size()); basis.at(gene) = 1;
        check(vector(owner, row.at("hodge_action")) == hodge.apply(basis), "selected exact Hodge field lost full support");
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc < 3 || argc > 4) throw std::invalid_argument("analysis archive, original input and optional measurements required");
        auto original = load_pager(argv[2]);
        if (argc == 4) import_measurements(original, argv[3]);
        const auto record = NativeRecord::read(argv[1]), fields = record.at("summary");
        check(record.at("schema").text() == "gadget.analysis.native.v1" && fields.at("schema").text() == "gadget.scientific.fields.v1", "analysis is a numerical report import");
        check(storage::encode_rfcx(record.endpoint()->source).bytes == storage::encode_rfcx(*original.complex).bytes, "analysis lost the complete original source");
        const auto original_publication = capture_workspace(original);
        for (const auto &bytes : original_publication.raw_sources()) {
            const auto sources = record.raw_sources();
            const auto found = std::find_if(sources.begin(), sources.end(), [&](const auto &x) { return x.key() == bytes.key(); });
            check(found != sources.end() && std::ranges::equal(found->bytes(), bytes.bytes()), "analysis lost original immutable source bytes");
        }
        const auto results = fields.at("bindings").at(2);
        check(results.at("contexts").size() == 4 && results.at("gene_sweeps").size() == 2 && results.at("cross_context").size() == 3, "analysis result family was truncated");
        std::set<std::uint64_t> namespaces;
        for (std::size_t i = 0; i < results.at("contexts").size(); ++i) {
            const auto child = publication(results.at("contexts").at(i));
            check(namespaces.insert(child.source_namespace()).second, "analysis contexts share an invented mathematical owner");
            context(child, original);
        }
        for (std::size_t i = 0; i < results.at("gene_sweeps").size(); ++i) {
            const auto sweep = results.at("gene_sweeps").at(i);
            for (const auto &name : {"mu1_AD", "mu2_AD", "mu1_HD", "mu2_HD", "cross_AD_HD", "delta_L_squared", "unclamped_delta_L_squared", "delta_L_norm"})
                check(vector(fields, sweep.at(name)).size() == original.genes.size(), "complete original gene sweep omitted rows");
        }
        for (std::size_t i = 0; i < results.at("cross_context").size(); ++i) {
            const auto cross = results.at("cross_context").at(i);
            const auto delta = vector(fields, cross.at("delta")), predicted = vector(fields, cross.at("predicted"));
            check(delta == numeric::sub(vector(fields, cross.at("response_HD")), vector(fields, cross.at("response_AD"))) &&
                      vector(fields, cross.at("residual")) == numeric::sub(delta, predicted), "finite cross context identity lost its actual fields");
        }
        const auto parent = fs::path(argv[1]).parent_path();
        for (std::size_t i = 0; i < record.at("outputs").size(); ++i) {
            const auto item = record.at("outputs").at(i);
            const auto path = parent / item.at("path").text();
            const auto bytes = rfc::storage::read_file_exact(path, source_byte_limit());
            check(std::ranges::equal(item.at("bytes").bytes(), bytes), "generated compatibility export lost its actual bytes");
        }
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
