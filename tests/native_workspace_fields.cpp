#include "../src/scientific_result_fields.hpp"
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <iostream>
#include <rfc/algebra/materialize.hpp>
#include <rfc/field/exact_hodge.hpp>
#include <rfc/storage/rfcx.hpp>
#include <set>

namespace {
using namespace gadget;
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
ExactMatrix matrix(const NativeRecord &owner, const NativeRecord &reference) {
    const auto &tensor = field(owner, reference);
    check(tensor.input && tensor.output && tensor.tensor.support().count() == tensor.tensor.domain().logical_size(),
          "matrix lost its typed axes or supported zeros");
    return storage::temporal_tensor_matrix(tensor);
}
void same_source(const NativeRecord &a, const NativeRecord &b) {
    check(storage::encode_rfcx(a.endpoint()->source).bytes == storage::encode_rfcx(b.endpoint()->source).bytes,
          "original complex, occurrences or attachments differ");
    check(a.endpoint()->tensors.size() <= b.endpoint()->tensors.size(), "original fields omitted");
    for (std::size_t i = 0; i < a.endpoint()->tensors.size(); ++i)
        check(storage::temporal_endpoint_detail::tensor_equal(a.endpoint()->tensors[i], b.endpoint()->tensors[i]),
              "original field coefficients, support, owner or role differ");
    for (const auto &original : a.raw_sources()) {
        const auto sources = b.raw_sources();
        const auto found = std::find_if(sources.begin(), sources.end(),
            [&](const auto &bytes) { return bytes.key() == original.key(); });
        check(found != sources.end() && std::ranges::equal(found->bytes(), original.bytes()),
              "original immutable byte owner omitted or changed");
    }
}
void persistence(const AssayWorkspace &workspace, const NativeRecord &request,
                 const NativeRecord &publication) {
    const auto retained = publication.at("result").at("result");
    check(retained.at("schema").text() == "gadget.persistence.fields.v1", "persistence is a display import");
    const auto study = workspace.study(request.at("assay").text());
    const auto original = NativeRecord::decode(retained.at("study_publication").bytes());
    check(original.encode() == study->native_publication().encode(), "persistence original study differs");
    const auto fields = retained.at("fields"), view = fields.at("view");
    same_source(original, fields);
    std::vector<std::vector<std::size_t>> births;
    for (std::size_t k = 0; k < original.endpoint()->source.grade_count(); ++k)
        births.emplace_back(original.endpoint()->source.grade({static_cast<std::uint16_t>(k)}).cell_count());
    for (std::size_t i = 0; i < request.at("births").size(); ++i) {
        const auto b = request.at("births").at(i);
        births.at(b.at("grade").integer().get_ui()).at(b.at("index").integer().get_ui()) =
            b.at("level").integer().get_ui();
    }
    std::vector<NativeRecord> stages;
    std::vector<algebra::MaterializedSubcomplex> materialized;
    std::set<std::uint64_t> namespaces;
    for (std::size_t t = 0; t < view.at("stages").size(); ++t) {
        const auto stage = NativeRecord::decode(view.at("stages").at(t).at("publication").bytes());
        stages.push_back(stage);
        check(namespaces.insert(stage.source_namespace()).second,
              "independent stage namespace collapsed by equal dimensions");
        const auto binding = fields.at("bindings").at(t);
        check(binding.at("source_namespace").text() == std::to_string(stage.source_namespace()) &&
                  binding.at("cutoff").integer() == rfc::numeric::Integer(std::to_string(stage.publication_history()->last_sequence())),
              "stage publication binding lost namespace or cutoff");
        algebra::SubcomplexMask mask(original.endpoint()->source);
        for (std::size_t k = 0; k < births.size(); ++k)
            for (std::size_t i = 0; i < births[k].size(); ++i)
                if (births[k][i] <= t) mask.select({{static_cast<std::uint16_t>(k)}, static_cast<std::uint32_t>(i)});
        materialized.push_back(algebra::materialize_subcomplex(original.endpoint()->source, mask));
        check(storage::encode_rfcx(stage.endpoint()->source).bytes ==
                  storage::encode_rfcx(materialized.back().complex).bytes,
              "actual stage lost an original selected cell or ordered occurrence");
        for (std::size_t k = 0; k < births.size(); ++k) {
            const auto grade = stage.at("view").at("grades").at(k), witness = stage.at("bindings").at(k);
            const auto b = matrix(stage, grade.at("boundary"));
            const auto c = matrix(stage, grade.at("cycle_frame"));
            const auto d = matrix(stage, grade.at("boundary_frame"));
            const auto upper = matrix(stage, witness.at("upper_boundary"));
            check(b.multiply(c).nonzeros() == 0 && b.multiply(d).nonzeros() == 0 &&
                      upper.rank() == d.cols() && c.rank() == c.cols(),
                  "stage cycle or image witnesses incomplete");
            const auto axes = binding.at("axes").at(k);
            check(field(stage, grade.at("cycle_frame")).output == witness.at("source_carrier").integer().get_ui() &&
                      field(stage, grade.at("cycle_frame")).input == witness.at("cycle_carrier").integer().get_ui(),
                  "stage cycle frame lost authentic source carrier");
            const auto source_carrier = std::find_if(stage.endpoint()->carriers.begin(), stage.endpoint()->carriers.end(),
                [&](const auto &axis) { return rfc::numeric::Integer(std::to_string(axis.key)) == axes.at("source_carrier").integer(); });
            check(source_carrier != stage.endpoint()->carriers.end() &&
                      source_carrier->kind == storage::TemporalCarrierKind::SourceGrade && source_carrier->grade.value == k &&
                      source_carrier->domain.logical_size() == births[k].size() -
                          static_cast<std::size_t>(std::count_if(births[k].begin(), births[k].end(), [&](auto level) { return level > t; })),
                  "stage source grade carrier is fabricated");
            const auto lift = matrix(fields, axes.at("ambient_lift"));
            const auto &to_ambient = materialized.back().to_ambient[k];
            check(lift.rows() == births[k].size() && lift.cols() == to_ambient.size(), "ambient lift shape differs");
            for (std::size_t i = 0; i < lift.rows(); ++i)
                for (std::size_t j = 0; j < lift.cols(); ++j)
                    check(lift.get(i, j) == (i == to_ambient[j] ? 1 : 0), "ambient source lift differs");
            if (witness.has("integral_homology")) {
                const auto integral = witness.at("integral_homology"), smith = integral.at("boundary_smith");
                const auto cycles = matrix(stage, integral.at("cycle_basis"));
                const auto coordinates = matrix(stage, integral.at("boundaries_in_cycle_coordinates"));
                const auto lower = matrix(stage, integral.at("lower_boundary"));
                const auto higher = matrix(stage, integral.at("upper_boundary"));
                const auto left = matrix(stage, smith.at("left")), right = matrix(stage, smith.at("right"));
                check(lower.multiply(cycles).nonzeros() == 0 && cycles.multiply(coordinates) == higher,
                      "integer cycle lattice witnesses fail their original boundary equations");
                check(left.multiply(coordinates).multiply(right) == matrix(stage, smith.at("diagonal")),
                      "Smith transformations do not reconstruct the retained integer presentation");
                check(left.determinant() == 1 || left.determinant() == -1, "Smith left map is not unimodular");
                check(right.determinant() == 1 || right.determinant() == -1, "Smith right map is not unimodular");
            }
        }
    }
    for (std::size_t at = 0; at < view.at("inclusions").size(); ++at) {
        const auto map = view.at("inclusions").at(at);
        const auto i = map.at("from").integer().get_ui(), j = map.at("to").integer().get_ui(), k = map.at("grade").integer().get_ui();
        const auto inclusion = matrix(fields, map.at("inclusion"));
        const auto &from = materialized[i].to_ambient[k], &to = materialized[j].to_ambient[k];
        for (std::size_t r = 0; r < to.size(); ++r)
            for (std::size_t c = 0; c < from.size(); ++c)
                check(inclusion.get(r, c) == (to[r] == from[c] ? 1 : 0), "inclusion identifies independent coordinate by dimension");
        check(inclusion.multiply(matrix(stages[i], stages[i].at("view").at("grades").at(k).at("cycle_frame"))) ==
                  matrix(fields, map.at("mapped_cycles")), "mapped cycle family incomplete");
        if (k != 0) {
            const auto lower_inclusion = matrix(fields, view.at("inclusions").at(at - 1).at("inclusion"));
            check(matrix(stages[j], stages[j].at("view").at("grades").at(k).at("boundary")).multiply(inclusion) ==
                      lower_inclusion.multiply(matrix(stages[i], stages[i].at("view").at("grades").at(k).at("boundary"))),
                  "retained inclusion is not an actual chain map");
        }
        bool found = false;
        for (std::size_t q = 0; q < fields.at("bindings").size(); ++q) {
            const auto binding = fields.at("bindings").at(q);
            if (!binding.has("rank_invariant") || binding.at("grade").integer() != k) continue;
            const auto &rank = field(fields, binding.at("rank_invariant")).tensor;
            const auto levels = stages.size();
            for (std::size_t a = 0; a < levels; ++a)
                for (std::size_t b = 0; b < levels; ++b)
                    check(rank.support().contains(a * levels + b) == (a <= b), "rank tensor confused an unrequested direction with supported zero");
            check(rank.state(i * levels + j).value() == Rational(map.at("homology_image_rank").integer()),
                  "native homology rank invariant differs from its actual inclusion");
            found = true;
            break;
        }
        check(found, "complete rank invariant tensor absent");
    }
}
void native_query(const AssayWorkspace &workspace, const NativeRecord &request,
                  const NativeRecord &publication) {
    const auto collection = publication.at("result").at("result");
    check(collection.at("schema").text() == "gadget.scientific.fields.v1", "partition answer is a display import");
    const auto answers = collection.at("view").at("answers");
    for (std::size_t i = 0; i < answers.size(); ++i) {
        const auto answer = answers.at(i);
        const auto p = answer.at("partition").integer().get_ui();
        const auto child = NativeRecord::decode(answer.at("certificate").at("publication").bytes());
        const auto study = workspace.study(request.at("assay").text(), p);
        const auto original = study->native_publication();
        const auto index = study->source_index_record();
        same_source(original, index);
        std::size_t cells{}, occurrences{};
        for (std::size_t k = 0; k < original.endpoint()->source.grade_count(); ++k) {
            const auto &grade = original.endpoint()->source.grade({static_cast<std::uint16_t>(k)});
            cells += grade.cell_count();
            occurrences += grade.incidence_count();
        }
        check(index.at("catalog_sections").at("cells").integer() == cells &&
                  index.at("catalog_sections").at("occurrences").integer() == occurrences &&
                  render_record_view(index).dump() == study->source_index().dump(),
              "native catalog count or compatibility view differs from its actual source");
        check(child.at("result").at("schema").text() == "gadget.scientific.result.fields.v1", "child lost native producer");
        check(NativeRecord::decode(child.at("request_publication").bytes()).encode() == request.owning_publication().encode(),
              "child lost the original whole outer request");
        same_source(original, child.at("result").at("fields"));
        check(NativeRecord::decode(child.at("result").at("source_study_publication").bytes()).encode() == original.encode(),
              "selected partition source publication differs");
        (void)replay_scientific_record(child);
        if (answer.has("source_index")) {
            const auto index = NativeRecord::decode(answer.at("source_index").at("publication").bytes());
            same_source(original, index);
            check(index.at("catalog_sections").at("carriers").integer() == original.endpoint()->carriers.size(),
                  "source index counts were reconstructed from a numerical display");
        }
        if (request.has("partitions"))
            check(p == request.at("partitions").at(i).integer(), "selected original partition order changed");
    }
}
void native_run(const AssayWorkspace &workspace, const NativeRecord &request,
                const NativeRecord &publication, const fs::path &out) {
    const auto fields = publication.at("result");
    check(fields.at("schema").text() == "gadget.scientific.fields.v1", "batch run imported numerical metadata");
    same_source(workspace.native_publication(), fields);
    const auto results = fields.at("view").at("results");
    const auto view = render_record_view(fields);
    check(view.at("schema").text() == "gadget.assays.run.v1" &&
              Json::parse(read_text(out / "results.json")).dump() == view.dump(), "run native and compatibility receipts differ");
    check(NativeRecord::read(out / "results.rfcxs").encode() == publication.encode(), "run root archive roundtrip differs");
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto receipt = results.at(i);
        const auto p = receipt.at("partition").integer().get_ui();
        check(p == request.at("partitions").at(i).integer(), "run reordered selected partitions");
        const auto path = out / receipt.at("native_path").text();
        check(receipt.at("bytes").integer() == fs::file_size(path) &&
                  receipt.at("sha256").text() == "sha256:" + sha256_file(path), "run file hash or exact bytes differs");
        const auto child = NativeRecord::read(path);
        check(child.at("result").at("schema").text() == "gadget.scientific.result.fields.v1", "run child lost actual native mathematics");
        check(NativeRecord::decode(child.at("request_publication").bytes()).encode() == request.owning_publication().encode(), "run recreated the original worker request");
        same_source(workspace.study(request.at("assay").text(), p)->native_publication(), child.at("result").at("fields"));
        check(child.at("request").at("operation").text() == request.at("query").at("operation").text(), "run selected the wrong original request child");
        (void)replay_scientific_record(child);
        if (receipt.has("source_index")) {
            const auto index = NativeRecord::decode(receipt.at("source_index").at("publication").bytes());
            check(index.at("catalog_sections").at("cells").kind() == NativeRecordKind::Integer, "run catalog was a numeric JSON import");
        }
    }
}

} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 4) throw std::invalid_argument("workspace, queries, output directory required");
        const AssayWorkspace workspace(argv[1]);
        const auto queries = read_record(argv[2]).native();
        const std::filesystem::path out(argv[3]);
        std::filesystem::create_directories(out);
        for (const auto &name : queries.names()) {
            std::cerr << "checking " << name << "\n";
            const auto request = queries.at(name);
            const bool is_run = name == "partition_run";
            const auto run_out = out / (name + "-files");
            const auto answer = is_run ? workspace.run_record(request, run_out) : workspace.execute_record(request);
            const auto operation = request.at("operation").text();
            check(AssayWorkspace::supports_native_publication(operation), "native dispatcher predicate excludes actual producer");
            check(NativeRecord::decode(answer.at("request_publication").bytes()).encode() == queries.encode(),
                  "whole original request publication was recreated");
            if (is_run) native_run(workspace, request, answer, run_out);
            else if (operation == "native_query") native_query(workspace, request, answer);
            else if (operation == "relational_persistence") persistence(workspace, request, answer);
            else if (operation == "inspect") same_source(workspace.native_publication(), answer);
            else if (operation == "query_batch") {
                const auto children = answer.at("result").at("result").at("results");
                for (std::size_t i = 0; i < children.size(); ++i) {
                    const auto child = NativeRecord::decode(children.at(i).at("child_publication").bytes());
                    check(child.at("result").at("schema").text() == "gadget.assay.result.native.v1", "batch child used a numerical display");
                }
            } else throw std::invalid_argument("unsupported test query");
            const auto path = out / (name + ".rfcxs");
            answer.write(path);
            check(NativeRecord::read(path).encode() == answer.encode(), "complete native archive roundtrip differs");
            write_text(out / (name + ".json"), render_record_view(answer.at("result")).dump());
        }
        std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
