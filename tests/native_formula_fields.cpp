#include "../src/source_grade_fields.hpp"
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/rfcx.hpp>
#include <iostream>
#include <set>

namespace {
using namespace gadget;
using namespace rfc;
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
const storage::TemporalTensor *tensor(const NativeRecord &owner, NativeRecord key) {
    return &assay_detail::retained_tensor(owner, key);
}
void vector(const storage::TemporalTensor &field, const numeric::RationalVector &expected) {
    check(field.tensor.domain().logical_size() == expected.size(), "field dimension changed");
    for (std::size_t i = 0; i < expected.size(); ++i)
        check(field.tensor.state(i) == expected[i], "native field value or support changed");
}
void run(const fs::path &root, const fs::path &expectation) {
    const auto expected = read_record(expectation).native();
    const AssayWorkspace workspace(root / "workspace.rfcxs");
    const auto publication = NativeRecord::read(root / "result.rfcxs");
    check(NativeRecord::decode(publication.encode()).encode() == publication.encode(),
          "native result changed on its own archive roundtrip");
    const auto certificate = publication.at("result"), fields = certificate.at("result");
    check(certificate.at("schema").text() == "gadget.assay.result.native.v1" &&
              fields.at("schema").text() == "gadget.sections.fields.v1", "result imports a numerical display");
    check(NativeRecord::decode(certificate.at("workspace_publication").bytes()).encode() ==
              workspace.native_publication().encode(), "complete original workspace changed");
    const auto n = expected.at("axes").size();
    std::vector<numeric::RationalVector> original(2), compatible(2), residual(2);
    for (std::size_t v = 0; v < 2; ++v)
        for (std::size_t i = 0; i < n; ++i)
            original[v].emplace_back(expected.at(v ? "right" : "left").at(i).text());
    for (std::size_t i = 0; i < n; ++i) {
        const numeric::Rational mean = (original[0][i] + original[1][i]) / 2;
        for (std::size_t v = 0; v < 2; ++v) {
            compatible[v].push_back(mean);
            residual[v].push_back(original[v][i] - mean);
        }
    }
    std::set<std::uint64_t> namespaces, carriers;
    std::vector<std::uint64_t> view_carriers;
    const auto section = fields.at("sections").at(0);
    check(fields.at("sections").size() == 1 && section.at("fields").size() == 2 &&
              section.at("complete").boolean() && section.at("reconstruction_verified").boolean(),
          "section omitted a source or supported coordinate");
    for (std::size_t v = 0; v < 2; ++v) {
        const std::string name = v ? "pubchem" : "chembl";
        const auto study = workspace.study(name);
        const auto pin = study->field_accession("all.reported_formula");
        check(namespaces.insert(pin.source_namespace()).second, "independent source identities were collapsed");
        vector(*pin.tensor_declaration(), original[v]);
        const auto retained = study->native_publication().at("source_provenance").at("sources");
        check(retained.size() == expected.at("sources").size(), "original byte declarations omitted");
        for (std::size_t s = 0; s < retained.size(); ++s) {
            const auto control = expected.at("sources").at(s);
            const auto raw = retained.at(s).at("bytes").bytes();
            const std::string_view bytes(reinterpret_cast<const char *>(raw.data()),raw.size());
            check(retained.at(s).at("id").text() == control.at("id").text() &&
                      "sha256:" + sha256(bytes) == control.at("sha256").text() &&
                      retained.at(s).at("sha256").text() == control.at("sha256").text(),
                  "retained original source bytes or identities changed");
        }
        const auto row = section.at("fields").at(v);
        check(row.at("view").text() == name && row.at("source_digest").text() == study->digest() &&
                  row.at("source_grade").integer() == pin.carrier()->grade.value,
              "section source binding changed");
        for (std::size_t i = 0; i < n; ++i)
            check(row.at("ids").at(i).text() == expected.at("axes").at(i).text(),
                  "element coordinate order changed");
        const auto input_key = row.at("input");
        const auto &input = *tensor(fields,input_key);
        check(input.output && carriers.insert(*input.output).second,
              "equal dimensions conflated independent view carriers");
        view_carriers.push_back(*input.output);
        vector(input,original[v]);
        vector(*tensor(fields,row.at("compatible")),compatible[v]);
        vector(*tensor(fields,row.at("residual")),residual[v]);
        for (const auto *component : {"compatible","residual"})
            check(tensor(fields,row.at(component))->output == input.output, "component left its source stalk");
    }
    numeric::RationalVector defect;
    numeric::Rational quadrance{};
    for (std::size_t i = 0; i < n; ++i) {
        defect.push_back(original[0][i]-original[1][i]);
        quadrance += residual[0][i]*residual[0][i] + residual[1][i]*residual[1][i];
    }
    vector(*tensor(fields,section.at("constraints").at(0).at("field")),defect);
    check(section.at("residual_quadrance").rational() == quadrance &&
              section.at("orthogonality").rational() == 0, "exact metric witnesses changed");
    vector(*tensor(fields,fields.at("metric_diagonal")),numeric::RationalVector(2*n,1));
    const auto delta_key = fields.at("coboundary");
    const auto &delta_tensor = *tensor(fields,delta_key);
    const auto delta = storage::temporal_tensor_matrix(delta_tensor);
    check(delta.rows() == n && delta.cols() == 2*n, "complete interface map shape changed");
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < 2*n; ++j)
            check(delta.get(i,j) == (j == i ? 1 : j == n+i ? -1 : 0), "signed element restriction changed");
    check(fields.at("restrictions").size() == 2, "original two role occurrences omitted");
    for (std::size_t v = 0; v < 2; ++v) {
        const auto restriction_key = fields.at("restrictions").at(v);
        const auto &restriction = *tensor(fields,restriction_key);
        const auto map = storage::temporal_tensor_matrix(restriction);
        check(restriction.input == view_carriers[v] && map.rows() == n && map.cols() == n,
              "restriction lost its actual stalk carrier");
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                check(map.get(i,j) == (i == j ? 1 : 0), "complete identity restriction changed");
    }
    check(fields.at("prepared_groups").size() == 1, "full support unexpectedly split");
    const auto group = fields.at("prepared_groups").at(0);
    const auto action = storage::temporal_tensor_matrix(*tensor(fields,group.at("action")));
    const auto kernel = storage::temporal_tensor_matrix(*tensor(fields,group.at("kernel_frame")));
    const auto dependencies = storage::temporal_tensor_matrix(*tensor(fields,group.at("row_dependencies")));
    check(action == delta && kernel.rows() == 2*n && kernel.cols() == n && kernel.rank() == n &&
              action.multiply(kernel).nonzeros() == 0 && dependencies.rows() == n && dependencies.cols() == 0,
          "complete compatible family or row dependence space omitted");
    const auto &source = fields.endpoint()->source.grade({1});
    check(source.cell_count() == 1 && source.incidence_count() == 2,
          "interface relation lost full ordered source arity");
    const auto replayed = workspace.execute_record(certificate.at("query"));
    check(render_record_view(replayed.at("result")).dump() == render_record_view(publication.at("result")).dump(),
          "cold native workspace replay changed complete result");
}
}
int main(int argc, char **argv) try {
    if (argc != 3) throw std::invalid_argument("usage: native_formula_fields DIRECTORY EXPECTED");
    run(argv[1],argv[2]);
    std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
    return 0;
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
