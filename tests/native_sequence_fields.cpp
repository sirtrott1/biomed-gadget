#include <gadget/assay.hpp>
#include <gadget/pipeline.hpp>
#include "../src/genomic_detail.hpp"
#include <gadget/record_view.hpp>
#include <atomic>
#include <rfc/storage/file_io.hpp>
#include <rfc/field/boundary.hpp>
#include <iostream>
#include <tuple>

using namespace gadget;
namespace {
void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
const rfc::storage::TemporalTensor &field(const NativeRecord &record, const NativeRecord &key) {
    for (const auto &tensor : record.endpoint()->tensors)
        if (key.integer() == rfc::numeric::Integer(std::to_string(tensor.key))) return tensor;
    throw std::runtime_error("native field is absent");
}
void refuse(const std::function<void()> &action) {
    bool rejected = false;
    try {action();} catch (const std::exception &) {rejected = true;}
    check(rejected, "invalid native request was admitted");
}
bool source_bytes_equal(const NativeRecord &record, std::string_view bytes) {
    const auto expected = std::as_bytes(std::span(bytes));
    return std::ranges::equal(record.bytes(), expected);
}
void verify_routes(const NativeRecord &record, std::string_view name = "routes") {
    const auto &costs = field(record, record.at("costs")).tensor;
    check(costs.support().count() == record.endpoint()->source.grade({1}).cell_count(),
          "route cost support omitted declared zero");
    const auto boundary = rfc::field::boundary_matrix(record.endpoint()->source, {1});
    for (std::size_t i = 0; i < record.at(name).size(); ++i) {
        const auto route = record.at(name).at(i);
        const auto &values = field(record, route.at("realized_field")).tensor;
        const auto &demand = field(record, route.at("boundary")).tensor;
        const auto dense = values.values().dense();
        Vector ordered(dense.size(), Rational{0});
        for (std::size_t j = 0; j < route.at("ordered_relations").size(); ++j)
            ordered.at(route.at("ordered_relations").at(j).integer().get_ui()) += 1;
        check(dense == ordered && boundary.apply(dense) == demand.values().dense() &&
                  rfc::numeric::dot(costs.values().dense(), dense) == route.at("minimum_cost").rational(),
              "native ordered route field, boundary or cost changed");
    }
}
}
int main(int argc, char **argv) {
    const fs::path root = argc > 1 ? argv[1] : "build/native-sequence-fields";
    fs::remove_all(root); fs::create_directories(root);
    const std::string reads = ">q0\r\nAC\r\n>q1\nAAA\n", reference = ">r0\nAC\n>r1\nAC\n";
    auto resources = Json::array();
    const std::string comparisons = "@q0\r\nAC\r\n+\r\nII\r\n@q0\nGT\n+\nII\n";
    const std::string edited_queries = "@same\r\nAACC\r\n+\r\nIIII\r\n@same\nAACC\n+\nIIII\n@inserted\nAATCC\n+\nIIIII\n";
    const std::string edited_reference = ">r0\r\nAAGT\r\nAGCC\r\n>r1\nAAGTAGCC\n";
    const std::string consensus_reads = "@same\r\nAC\r\n+\r\nII\r\n@same\nAG\n+\nII\n";
    for (const auto &[id, bytes] : std::vector<std::pair<std::string, std::string>>{{"reads", reads}, {"reference", reference}, {"ambiguous", ">a\nAN\n>b\nTA\n"}, {"quality", "@q0\nAC\n+\n!I\n"}, {"comparisons", comparisons}, {"history-query", "@same\nAAA\n+\nIII\n@same\nAAA\n+\nIII\n"}, {"history-reference", ">r0\nAAA\n>r1\nAAA\n"}, {"cover", ">a\nACA\n>b\nCAC\n>c\nACA\n"}, {"classes", "2\n1\nA\nB\n2 0 1 3\n"}}) {
        const auto filename = id + (id == "quality" ? ".fastq" : ".fa");
        write_text(root / filename, bytes);
        auto resource = Json::parse(R"({"id":"","format":"fasta","source":{"path":"","sha256":""},"method":"Authored source fixture"})");
        resource.set("id", Json::string(id)); resource.at("source").set("path", Json::string(filename));
        if (id == "quality" || id == "comparisons" || id == "history-query") resource.set("format", Json::string("fastq"));
        if (id == "classes") resource.set("format", Json::string("salmon_eq"));
        resource.at("source").set("sha256", Json::string("sha256:" + sha256(bytes))); resources.push(resource);
    }
    for (const auto &[id, bytes, format] : std::vector<std::tuple<std::string, std::string, std::string>>{
            {"edited-query", edited_queries, "fastq"}, {"edited-reference", edited_reference, "fasta"},
            {"consensus", consensus_reads, "fastq"}}) {
        const auto filename = id + ".source"; write_text(root / filename, bytes);
        auto resource = Json::parse(R"({"id":"","format":"","source":{"path":"","sha256":""},"method":"Original edited splice fixture"})");
        resource.set("id", Json::string(id)); resource.set("format", Json::string(format));
        resource.at("source").set("path", Json::string(filename));
        resource.at("source").set("sha256", Json::string("sha256:" + sha256(bytes))); resources.push(resource);
    }
    auto spec = Json::parse(R"({"schema":"gadget.assays.import.v1","id":"sequences","design":{"question":"Preserve ordered source witnesses","population":"Authored fixture","selection":"All source records","limitations":["Verification fixture"]},"assays":[],"resources":[],"links":[]})");
    spec.set("resources", resources); write_text(root / "import.json", spec.dump());
    import_assays(root / "import.json", root / "input");
    auto workspace = std::make_shared<const AssayWorkspace>(root / "input/workspace.rfcxs");
    auto alignment = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"exact_sequence_alignment","query_resource":"reads","query_id":"q0","reference_resource":"reference","reference_id":"r0","mode":"global_linear_gap","ambiguity":"literal","scoring":{"match":"0","mismatch":"2","query_gap":"0","reference_gap":"0","units":"declared cost","evidence":"Authored edit rule"}})");
    const auto aligned = workspace->execute_record(admit_record_view(alignment));
    const auto alignment_fields = aligned.at("result").at("result");
    const auto &complex = aligned.endpoint()->source;
    check(alignment_fields.at("schema").text() == "gadget.sequence.alignment.fields.v1" &&
          complex.grade({0}).cell_count() == 9 && complex.grade({1}).cell_count() == 16,
          "complete edit carrier was replaced");
    const auto &optimal = field(alignment_fields, alignment_fields.at("optimal_support")).tensor;
    check(optimal.support().count() == 16, "supported nonminimizing zero was omitted");
    check(alignment_fields.at("route_witnesses").size() > 0 &&
          alignment_fields.at("route_witnesses").at(0).at("minimum_cost").rational() == 0,
          "ordered minimum route was lost");
    const auto lift_key = alignment_fields.at("ordered_source_lift");
    const auto &lift = field(alignment_fields, lift_key);
    check(lift.role == rfc::storage::TemporalTensorRole::Observation && lift.output == 2,
          "ordered routes have no source lift");
    auto assembly = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"sequence_assembly_graph","resource":"reads","k":2,"ambiguity":"literal_tokens","evidence":"All supplied reads","route":{"start":"A","end":"C","objective":"fewest_steps"}})");
    const auto assembled = workspace->execute_record(admit_record_view(assembly));
    const auto assembly_fields = assembled.at("result").at("result");
    const auto motif_request = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"sequence_motifs","resource":"reads","motifs":[{"id":"AC","pattern":"AC"},{"id":"GT","pattern":"GT"},{"id":"long","pattern":"CCCC"},{"id":"any","pattern":"N"}],"strands":"both"})");
    const auto motif_publication = workspace->execute_record(admit_record_view(motif_request));
    const auto motif_fields = motif_publication.at("result").at("result");
    check(motif_fields.at("schema").text() == "gadget.sequence.motifs.fields.v1" &&
              motif_publication.endpoint()->source.grade({0}).cell_count() == 6 &&
              motif_publication.endpoint()->source.grade({1}).cell_count() == 12 &&
              motif_publication.endpoint()->source.grade({1}).incidence_count() == 24,
          "motif original sequence/pattern/hit occurrences changed");
    const auto possible_key = motif_fields.at("possible_counts");
    const auto &possible = field(motif_fields, possible_key);
    check(possible.output == motif_fields.at("summary_carrier").integer().get_ui() &&
              possible.tensor.values().dense() == Vector({1, 1, 0, 4, 0, 0, 0, 6}) &&
              possible.tensor.support().count() == 8,
          "motif complete summary carrier or supported zero changed");
    const auto motif_view = render_record_view(motif_publication.at("result"));
    check(motif_view.at("result").at("hits").at(1).at("strand").text() == "-" &&
              motif_view.at("result").at("hits").at(1).at("motif").text() == "GT" &&
              motif_view.at("result").at("summaries").at(2).at("possible_occurrences").integer() == 0,
          "motif hit order/orientation/zero view changed");
    check(workspace->execute(motif_request).at("result").dump() ==
              motif_view.at("result").dump(), "motif display/native producers differ");
    motif_publication.write(root / "motifs.rfcxs");
    check(NativeRecord::read(root / "motifs.rfcxs").encode() == motif_publication.encode(),
          "motif native roundtrip changed");
    auto ambiguous_request = Json::parse(motif_request.dump());
    ambiguous_request.set("resource", Json::string("ambiguous"));
    const auto ambiguous_publication = workspace->execute_record(admit_record_view(ambiguous_request));
    const auto ambiguous_fields = ambiguous_publication.at("result").at("result");
    check(field(ambiguous_fields, ambiguous_fields.at("possible_counts")).tensor.values().dense() ==
              Vector({1, 1, 0, 4, 0, 0, 0, 4}) &&
              field(ambiguous_fields, ambiguous_fields.at("definite_counts")).tensor.values().dense() ==
              Vector({0, 0, 0, 4, 0, 0, 0, 4}),
          "ambiguous motif choices were replaced by a literal base or joint extrema");
    auto invalid_motif = Json::parse(motif_request.dump());
    invalid_motif.at("motifs").at(1).set("id", Json::string("AC"));
    refuse([&] { workspace->execute_record(admit_record_view(invalid_motif)); });
    invalid_motif = Json::parse(motif_request.dump());
    invalid_motif.at("motifs").at(0).set("pattern", Json::string("?"));
    refuse([&] { workspace->execute_record(admit_record_view(invalid_motif)); });
    check(assembled.endpoint()->source.grade({1}).cell_count() == 2 &&
          assembled.endpoint()->source.grade({1}).incidence_count() == 4,
          "repeated self-loop boundary occurrences were discarded");
    const auto count_map_key = assembly_fields.at("count_map");
    const auto &count_map = field(assembly_fields, count_map_key);
    check(count_map.tensor.domain().extents() == std::vector<std::size_t>({2, 3}) &&
          count_map.tensor.support().count() == 6 && count_map.tensor.values().nonzeros() == 3,
          "actual assembly occurrence map was not retained");
    check(field(assembly_fields, assembly_fields.at("kmer_occurrence_count")).tensor.state(1).value() == 2 &&
          field(assembly_fields, assembly_fields.at("source_offset")).tensor.state(2).value() == 1 &&
          assembled.endpoint()->source.boundary_apply({1}, Vector{0, 1}) == Vector({0, 0}),
          "repeated read offsets or self-loop cancellation differ");
    auto constraints = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"constrained_allocations","domain":"bounded_integer","variables":[{"id":"x","lower":"0","upper":"2"},{"id":"y","lower":"0","upper":"2"}],"constraints":[{"id":"sum","coefficients":["1","1"],"sense":"eq","rhs":"2"}],"objective":{"sense":"minimize","coefficients":["0","0"]},"units":"counts","evidence":"Authored finite domain","witness_representation":"factorized"})");
    const auto constrained = workspace->execute_record(admit_record_view(constraints));
    const auto family = constrained.at("result").at("result");
    check(family.at("optimal_points_or_vertices").integer() == 3 && family.at("blocks").size() == 1 &&
          field(family, family.at("blocks").at(0).at("optimal_parameter_tuples")).tensor.domain().extents() ==
              std::vector<std::size_t>({1, 3}), "complete optimizing factor family was lost");
    check(field(family, family.at("tensors").at("objective")).tensor.support().count() == 2 &&
          field(family, family.at("tensors").at("objective")).tensor.values().nonzeros() == 0,
          "declared zero objective was confused with absence");
    auto expression = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"expression_allocations","samples":[{"id":"sample","resource":"reads","scale":"2","scale_provenance":"Authored positive scale"}],"evidence_source":{"mode":"transcript_sequences","reference_resource":"reference","strands":["+"],"ambiguity":"literal","maximum_mismatches":0,"evidence":"All declared transcript windows"},"targets":[{"id":"r0","candidates":["r0"]},{"id":"r1","candidates":["r1"]},{"id":"unassigned","candidates":["__unassigned__"]}],"contrasts":[],"map_evidence":"Authored targets","coupling":{"witness_representation":"factorized","evidence":"Authored count constraint","constraints":[{"id":"r0_count","terms":[{"sample":"sample","target":"r0","coefficient":"1"}],"sense":"eq","rhs":"1","evidence":"Authored allocation hypothesis"}]}})");
    const auto expressed = workspace->execute_record(admit_record_view(expression));
    const auto expression_fields = expressed.at("result").at("result"), coupling = expression_fields.at("coupling");
    check(coupling.at("status").text() == "feasible" && coupling.at("solvers").size() > 0,
          "native coupled family was not retained");
    for (std::size_t i = 0; i < coupling.at("solvers").size(); ++i) {
        const auto retained = coupling.at("solvers").at(i);
        const auto solver = NativeRecord::decode(retained.at("constraint_publication").bytes());
        check(rfc::numeric::Integer(std::to_string(solver.source_namespace())) == retained.at("variable_binding").at("source_namespace").integer() &&
              solver.at("schema").text() == "gadget.constraints.fields.v1",
              "coupled solver source binding differs");
    }
    const auto target_coordinate = expression_fields.at("allocation_fields").at("target_coordinates").tensor().state(0)->get_num().get_ui();
    check(field(expression_fields, coupling.at("samples").at(0).at("lower")).tensor.state(target_coordinate).value() == 2,
          "coupled expression scale or typed optimum differs");
    std::vector<std::pair<std::string, Json>> recipes{{"alignment", alignment}, {"assembly", assembly},
                                                    {"constraints", constraints}, {"expression", expression}};
    std::vector<std::shared_ptr<const AssayWorkspace>> results;
    for (const auto &[id, recipe] : recipes) {
        auto request = Json::parse(R"({"schema":"gadget.assays.materialize.v1","id":"","source_namespace":"730010","channel":"observation","design":{},"recipe":{}})");
        request.set("id", Json::string(id)); request.set("source_namespace", Json::string(std::to_string(730010 + results.size())));
        request.set("design", spec.at("design")); request.set("recipe", recipe);
        auto result = workspace->materialize_workspace(admit_record_view(request), root / id);
        check(!fs::exists(root / id / "workspace.json"), "native publication created a JSON authority");
        AssayWorkspace reopened(root / id / "workspace.rfcxs");
        check(reopened.native_publication().encode() == result->native_publication().encode(),
              "native publication changed on reopen");
        results.push_back(std::move(result));
    }
    auto impossible = Json::parse(constraints.dump()); impossible.at("constraints").at(0).set("rhs", Json::string("5"));
    const auto failed = workspace->execute_record(admit_record_view(impossible)).at("result").at("result");
    check(failed.at("summary").at("status").text() == "infeasible" &&
          field(failed, failed.at("tensors").at("lower")).tensor.support().count() == 0,
          "infeasible bounds have numerical support");
    auto batch = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"query_batch","concurrent_jobs":2,"queries":[]})");
    const std::vector<std::pair<std::string, Json>> batch_jobs{
        {"motifs", motif_request}, {"missing", impossible}, {"alignment", alignment},
        {"motifs_again", motif_request}, {"assembly", assembly},
        {"quality", Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_quality","resource":"quality","quality_thresholds":[0,20]})")}};
    for (const auto &[id, query] : batch_jobs) {
        auto job = Json::object(); job.set("id", Json::string(id)); job.set("query", query);
        batch.at("queries").push(job);
    }
    const auto hardware = rfc::runtime::detect_hardware();
    rfc::runtime::CpuResourcePolicy policy;
    policy.worker_limit = 2; policy.declared_buffer_byte_limit = 128u * 1024u * 1024u;
    const auto cpu = std::make_shared<CpuResources>(hardware, policy, [hardware] { return hardware; });
    const auto batch_request = admit_record_view(batch);
    const auto batched = workspace->execute_record(batch_request, cpu);
    const auto batch_fields = batched.at("result").at("result");
    const auto entries = batch_fields.at("results");
    const auto batch_view = render_record_view(batch_fields);
    check(workspace->execute(batch, cpu).at("result").dump() == batch_view.dump(),
          "public batch display differs from the owning native producer");
    check(AssayWorkspace::supports_native_publication("query_batch") &&
              batch_fields.at("schema").text() == "gadget.assay.batch.fields.v1" &&
              entries.size() == batch_jobs.size() &&
              batched.source_namespace() == workspace->native_publication().source_namespace() &&
              batched.endpoint()->tensors.size() == workspace->native_publication().endpoint()->tensors.size(),
          "batch substituted a common mathematical carrier for independent children");
    check(NativeRecord::decode(batched.at("request_publication").bytes()).encode() == batch_request.encode(),
          "batch lost the complete original request publication");
    std::vector<NativeRecord> children;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto entry = entries.at(i);
        const auto child = NativeRecord::decode(entry.at("child_publication").bytes());
        check(entry.at("id").text() == batch_jobs[i].first &&
                  entry.at("admitted_workers").integer() == 1 &&
                  render_record_view(child.at("request")).dump() ==
                      render_record_view(batch_request.at("queries").at(i).at("query")).dump(),
              "batch changed the admitted child selection or declaration order");
        check(NativeRecord::decode(child.at("request_publication").bytes()).encode() == batch_request.encode(),
              "selected child lost its original outer publication");
        check(render_record_view(child.at("result")).dump() ==
                  batch_view.at("results").at(i).at("answer").dump(),
              "batch view did not read the retained child owner");
        check(render_record_view(child.at("result")).at("result").dump() ==
                  workspace->execute(batch_jobs[i].second, cpu).at("result").dump(),
              "batch result differs from the complete direct result");
        child.write(root / ("batch-child-" + std::to_string(i) + ".rfcxs"));
        children.push_back(child);
    }
    const auto first_motifs = children.at(0).at("result").at("result");
    const auto repeated_motifs = children.at(3).at("result").at("result");
    const auto missing_fields = children.at(1).at("result").at("result");
    check(first_motifs.source_namespace() != repeated_motifs.source_namespace() &&
              first_motifs.at("possible_counts").integer() == repeated_motifs.at("possible_counts").integer() &&
              field(first_motifs, first_motifs.at("possible_counts")).tensor.state(2).value() == 0 &&
              field(missing_fields, missing_fields.at("tensors").at("lower")).tensor.support().count() == 0,
          "independent equal keys, supported zero and missing support were collapsed");
    rfc::storage::NativeAccessionRequest selected;
    selected.subject = first_motifs.at("possible_counts").integer().get_ui();
    const auto pin = children.at(0).access(selected);
    check(pin.source_namespace() == children.at(0).source_namespace() &&
              pin.sequence() == children.at(0).endpoint()->sequence &&
              pin.tensor().get() == &field(first_motifs, first_motifs.at("possible_counts")).tensor &&
              children.at(3).access(selected).source_namespace() != pin.source_namespace(),
          "child field accession used another source or prefix");
    refuse([&] { batched.access(selected); });
    refuse([&] { children.at(0).historical_publication(children.at(0).endpoint()->sequence + 1); });
    auto reversed = Json::parse(batch.dump());
    auto reverse_jobs = Json::array();
    for (std::size_t i = batch_jobs.size(); i > 0; --i) reverse_jobs.push(batch.at("queries").at(i - 1));
    reversed.set("queries", reverse_jobs);
    const auto reversed_view = assay_detail::query_batch(*workspace, reversed, cpu);
    for (std::size_t i = 0; i < batch_jobs.size(); ++i)
        check(reversed_view.at("results").at(i).at("id").text() == batch_jobs[batch_jobs.size() - i - 1].first &&
                  reversed_view.at("results").at(i).at("answer").at("result").dump() ==
                      batch_view.at("results").at(batch_jobs.size() - i - 1).at("answer").at("result").dump(),
              "reordered batch changed source results or publication order");
    auto invalid_batch = Json::parse(batch.dump());
    invalid_batch.at("queries").at(1).set("id", Json::string("motifs"));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch = Json::parse(batch.dump());
    invalid_batch.at("queries").at(0).set("query", batch);
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch = Json::parse(batch.dump()); invalid_batch.set("concurrent_jobs", Json::integer(0));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch.set("concurrent_jobs", Json::integer(3));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch = Json::parse(batch.dump()); invalid_batch.set("queries", Json::array());
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch = Json::parse(batch.dump());
    invalid_batch.at("queries").at(0).at("query").set("workers", Json::integer(0));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch.at("queries").at(0).at("query").set("workers", Json::integer(3));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    invalid_batch = Json::parse(batch.dump());
    invalid_batch.at("queries").at(2).at("query").set("operation", Json::string("unknown"));
    refuse([&] { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); });
    const auto unconverted = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"allocation_bounds","resource":"classes","targets":[{"id":"A","candidates":["A"]}],"candidate_search":{"method":"Authored candidate set","completeness":"All declared candidates","reference_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000"},"scale":"1","scale_provenance":"Original raw count"})");
    check(workspace->execute(unconverted).at("result").at("samples").size() == 1,
          "unconverted-child refusal fixture is not an otherwise valid query");
    invalid_batch = Json::parse(batch.dump());
    invalid_batch.at("queries").at(batch_jobs.size() - 1).set("query", unconverted);
    if (!AssayWorkspace::supports_native_publication("allocation_bounds")) {
        bool unconverted_refused = false;
        try { assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu); }
        catch (const std::invalid_argument &error) {
            unconverted_refused = std::string(error.what()) ==
                                  "scientific study: source query child has no native publisher";
        }
        check(unconverted_refused, "mixed batch admitted an unconverted mathematical child");
    } else {
        const auto admitted = assay_detail::query_batch_record(*workspace, admit_record_view(invalid_batch), cpu);
        const auto child = NativeRecord::decode(admitted.at("result").at("result").at("results")
            .at(batch_jobs.size() - 1).at("child_publication").bytes());
        check(render_record_view(child.at("result")).at("result").dump() == workspace->execute(unconverted).at("result").dump(),
              "new native publisher batch child differs from its complete direct result");
    }
    auto small_policy = policy; small_policy.declared_buffer_byte_limit = sizeof(std::size_t) - 1;
    const auto small_cpu = std::make_shared<CpuResources>(hardware, small_policy, [hardware] { return hardware; });
    refuse([&] { assay_detail::query_batch_record(*workspace, batch_request, small_cpu); });
    std::atomic_size_t observations{};
    const auto cancelled_cpu = std::make_shared<CpuResources>(hardware, policy, [&] {
        if (++observations > 1) throw rfc::runtime::ExecutionCancelled{};
        return hardware;
    });
    bool cancelled = false;
    try { assay_detail::query_batch_record(*workspace, batch_request, cancelled_cpu); }
    catch (const rfc::runtime::ExecutionCancelled &) { cancelled = true; }
    check(cancelled, "native child cancellation was converted to a published partial batch");
    { auto recovered = cancelled_cpu->coordinator().reserve(2); }
    { auto recovered = cpu->coordinator().reserve(2); }
    const std::vector<Json> searches{
        Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"indexed_sequence_alignment","query_resource":"reads","reference_resource":"reference","strands":["+","-"],"maximum_edits":1,"scoring":"unit_edit_cost","evidence":"Every original placement"})"),
        Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"reference_sequence_matches","query_resource":"comparisons","reference_resource":"reference","strands":["+","-"],"ambiguity":"literal","maximum_mismatches":1,"junctions":[],"minimum_overhang":1,"evidence":"Every original comparison"})"),
        Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"sequence_splice_discovery","query_resource":"reads","reference_resource":"reference","strands":["+"],"minimum_exon_bases":2,"minimum_intron_bases":4,"maximum_intron_bases":10,"maximum_junctions":1,"motifs":["GT-AG"],"evidence":"Every original exon chain"})"),
        Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_overlap_reconstruction","resource":"cover","minimum_overlap":2,"objective":"minimum_suffix_prefix_read_cover","evidence":"Every original read"})")};
    std::vector<NativeRecord> search_publications;
    auto history_request = Json::parse(searches[0].dump());
    history_request.set("query_resource", Json::string("history-query"));
    history_request.set("reference_resource", Json::string("history-reference"));
    history_request.set("strands", Json::parse(R"(["+"])"));
    const auto one_history = workspace->execute(history_request, cpu).at("result");
    history_request.set("edit_histories", Json::string("one_minimum"));
    check(workspace->execute(history_request, cpu).at("result").dump() == one_history.dump(),
          "explicit one-history mode changed the compatibility result");
    history_request.set("edit_histories", Json::string("all_minimum"));
    const auto all_history_publication = workspace->execute_record(admit_record_view(history_request), cpu);
    const auto all_history = all_history_publication.at("result").at("result");
    check(all_history.at("view").at("matches").size() == 12 &&
              all_history.at("view").at("work").at("minimum_edit_histories").integer() == 28,
          "complete tied edit history family or repeated original records changed");
    bool absent_suffix = false;
    for (std::size_t c = 0; c < all_history.at("children").size(); ++c)
        verify_routes(NativeRecord::decode(all_history.at("children").at(c).bytes()));
    for (std::size_t p = 0; p < all_history.at("placements").size(); ++p) {
        const auto placement = all_history.at("placements").at(p);
        const auto child = NativeRecord::decode(placement.at("child_publication").bytes());
        const auto evidence = placement.at("minimum_history_evidence");
        const auto &suffix = field(child, evidence.at("suffix")).tensor;
        absent_suffix |= suffix.support().count() != suffix.domain().logical_size();
        const auto operations = placement.at("minimum_history_operations");
        const auto count = all_history.at("view").at("matches").at(p).at("minimum_edit_history_count").integer();
        check(count == operations.size() && (count == 1 || count == 3), "native family selectors lost minimum histories");
        const auto primary = all_history.at("view").at("matches").at(p).at("route_index").integer();
        check(primary == placement.at("operation_index").integer(), "primary route context became a family ordinal");
    }
    check(absent_suffix, "unreachable edit suffix became supported zero");
    all_history_publication.write(root / "all-history.rfcxs");
    check(NativeRecord::read(root / "all-history.rfcxs").encode() == all_history_publication.encode(),
          "all-history native archive roundtrip changed");
    auto invalid_history = Json::parse(history_request.dump());
    invalid_history.set("edit_histories", Json::string("sampled"));
    refuse([&] { workspace->execute_record(admit_record_view(invalid_history), cpu); });
    auto edited_request = Json::parse(searches[2].dump());
    edited_request.set("query_resource", Json::string("edited-query"));
    edited_request.set("reference_resource", Json::string("edited-reference"));
    edited_request.set("maximum_edits", Json::integer(1));
    edited_request.set("scoring", Json::string("unit_edit_cost"));
    const auto edited_publication = workspace->execute_record(admit_record_view(edited_request), cpu);
    const auto edited_fields = edited_publication.at("result").at("result");
    check(source_bytes_equal(edited_fields.at("query_bytes"), edited_queries) &&
              source_bytes_equal(edited_fields.at("reference_bytes"), edited_reference) &&
              edited_fields.at("queries").size() == 3 && edited_fields.at("references").size() == 2 &&
              edited_fields.at("queries").at(0).at("id").text() == edited_fields.at("queries").at(1).at("id").text(),
          "edited splice changed original repeated records, line framing or quality bytes");
    bool edited_intron = false, edited_missing = false, edited_error = false;
    for (std::size_t c = 0; c < edited_fields.at("children").size(); ++c) {
        const auto child = NativeRecord::decode(edited_fields.at("children").at(c).bytes());
        verify_routes(child);
        check(child.source_namespace() != edited_fields.source_namespace() &&
                  field(child, child.at("operation_codes")).tensor.support().count() ==
                      child.endpoint()->source.grade({1}).cell_count(),
              "edited splice operation fields lost native ownership or supported zero");
        const auto &prefix = field(child, child.at("query_prefix")).tensor;
        check(prefix.support().count() + child.at("view").at("terminal_vertices").size() ==
                  child.endpoint()->source.grade({0}).cell_count(),
              "terminal observation falsely acquired a biological prefix coordinate");
    }
    rfc::numeric::Integer edited_count{0};
    for (std::size_t g = 0; g < edited_fields.at("view").at("minimum_history_groups").size(); ++g)
        edited_count += edited_fields.at("view").at("minimum_history_groups").at(g).at("histories").integer();
    check(edited_count == edited_fields.at("view").at("work").at("minimum_edit_histories").integer() &&
              edited_count == rfc::numeric::Integer(std::to_string(edited_fields.at("placements").size())),
          "edited splice exact complete family counts changed");
    for (std::size_t p = 0; p < edited_fields.at("placements").size(); ++p) {
        const auto placement = edited_fields.at("placements").at(p), match = edited_fields.at("view").at("matches").at(p);
        const auto child = NativeRecord::decode(placement.at("child_publication").bytes());
        check(placement.at("operation_index").integer() == match.at("route_index").integer() &&
                  placement.at("reference_introns").size() == match.at("junctions").size() &&
                  placement.at("minimum_history_group").integer() == match.at("minimum_history_group").integer() &&
                  !placement.has("minimum_history_evidence") && !match.has("minimum_history_evidence"),
              "edited splice parent does not select its actual route or intron spans");
        const auto group = edited_fields.at("view").at("minimum_history_groups")
            .at(placement.at("minimum_history_group").integer().get_ui());
        check(group.at("carrier").integer() == placement.at("child_context").integer(),
              "edited splice certificate group refers to a different native child owner");
        const auto evidence = group.at("minimum_history_evidence");
        const auto &suffix = field(child, evidence.at("suffix")).tensor;
        edited_missing |= suffix.support().count() < child.endpoint()->source.grade({0}).cell_count();
        edited_error |= match.at("minimum_edits").rational() != 0;
        for (std::size_t j = 0; j < placement.at("reference_introns").size(); ++j) {
            const auto intron = placement.at("reference_introns").at(j);
            std::string bytes;
            for (std::size_t s = 0; s < intron.at("bytes").size(); ++s) bytes += intron.at("bytes").at(s).text();
            check(bytes == "GTAG", "edited splice intron lost original multi-line byte spans"); edited_intron = true;
        }
    }
    check(edited_intron && edited_missing && edited_error, "edited native splice witness controls were not exercised");
    check(workspace->execute(edited_request, cpu).at("result").dump() == render_record_view(edited_fields).dump(),
          "edited splice display recomputed or changed retained native mathematics");
    edited_publication.write(root / "edited-splice.rfcxs");
    check(NativeRecord::read(root / "edited-splice.rfcxs").encode() == edited_publication.encode(),
          "edited splice complete archive I/O changed");
    const auto consensus_request = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_edit_consensus","resource":"consensus","strands":["+"],"maximum_edits":1,"minimum_paired_bases":2,"objective":"minimum_length_then_total_unit_edits","ambiguity":"literal_ACGT","evidence":"Every joint source alignment"})");
    const auto consensus_publication = workspace->execute_record(admit_record_view(consensus_request), cpu);
    const auto consensus_fields = consensus_publication.at("result").at("result");
    check(source_bytes_equal(consensus_fields.at("source_bytes"), consensus_reads) &&
              consensus_fields.at("original_reads").size() == 2 && consensus_fields.at("original_tokens").size() == 4 &&
              consensus_fields.at("view").at("family_count").integer() == 2 &&
              consensus_fields.at("view").at("minimum_edits").integer() == 1,
          "consensus original repeated read/quality bytes or tied minimizing families changed");
    const auto consensus_optimization = NativeRecord::decode(consensus_fields.at("optimization_publication").bytes());
    verify_routes(consensus_optimization);
    check(consensus_optimization.source_namespace() != consensus_fields.source_namespace(),
          "consensus binary optimization replaced its independent biological source");
    const auto association_key = consensus_fields.at("paired_occurrence_map");
    const auto &association = field(consensus_fields, association_key);
    const auto &biological_costs = field(consensus_fields, consensus_fields.at("edit_costs")).tensor;
    const auto &biological = consensus_fields.endpoint()->source.grade({1});
    check(association.input == 1 && association.output == 1 &&
              biological_costs.support().count() == biological.cell_count(),
          "consensus source maps or supported zero edit costs lost their carriers");
    for (std::size_t f = 0; f < consensus_fields.at("view").at("families").size(); ++f) {
        const auto family = consensus_fields.at("view").at("families").at(f);
        check(family.at("read_alignments").size() == 2 && family.at("tokens").size() == 2,
              "consensus discarded complete original read contributions");
        for (std::size_t p = 0; p < family.at("tokens").size(); ++p) {
            const auto relation = family.at("tokens").at(p).at("relation").integer().get_ui() - 1;
            const auto range = biological.incidence_range(relation);
            check(range.second - range.first == 6, "joint consensus evidence was expanded into binary pairs");
            const auto output = biological.lower_index[range.first + 1];
            const auto width = consensus_fields.endpoint()->source.grade({0}).cell_count();
            Rational total{0};
            for (std::size_t input = 0; input < width; ++input) total += association.tensor.state(output * width + input).value();
            check(total == 2, "consensus paired occurrence association changed or erased a contradictory token");
        }
    }
    check(workspace->execute(consensus_request, cpu).at("result").dump() == render_record_view(consensus_fields).dump(),
          "consensus compatibility view replaced retained native results");
    consensus_publication.write(root / "consensus.rfcxs");
    check(NativeRecord::read(root / "consensus.rfcxs").encode() == consensus_publication.encode(),
          "complete consensus archive changed during I/O");
    auto consensus_defaults = Json::parse(R"({"schema":"gadget.pipeline.default.v1","profile":"sequence","question":"Every joint source alignment","goals":[{"id":"consensus","method":"read_edit_consensus","parameters":{}}]})");
    auto consensus_parameters = Json::object();
    json_object_object_foreach(consensus_request.raw(), key, value) {
        (void)value;
        if (std::string_view(key) != "schema" && std::string_view(key) != "operation")
            consensus_parameters.set(key, consensus_request.at(key));
    }
    consensus_defaults.at("goals").at(0).set("parameters", consensus_parameters);
    const auto consensus_plan = admit_record_view(default_pipeline(*workspace, consensus_defaults));
    const auto consensus_pipeline = execute_pipeline_record(workspace, consensus_plan,
        root / "consensus-pipeline", {}, cpu);
    check(consensus_pipeline.at("result").at("status").text() == "complete",
          "default consensus pipeline did not execute");
    const auto consensus_artifact = NativeRecord::read(root / "consensus-pipeline/consensus/result.rfcxs");
    check(source_bytes_equal(consensus_artifact.at("result").at("result").at("source_bytes"), consensus_reads) &&
              render_record_view(consensus_artifact.at("result").at("result")).dump() ==
                  render_record_view(consensus_fields).dump(),
          "pipeline changed original consensus bytes or minimizing families");
    consensus_artifact.write(root / "consensus-pipeline-roundtrip.rfcxs");
    check(NativeRecord::read(root / "consensus-pipeline-roundtrip.rfcxs").encode() == consensus_artifact.encode(),
          "pipeline consensus changed complete native owners during archive I/O");
    for (std::size_t i = 0; i < searches.size(); ++i) {
        const auto publication = workspace->execute_record(admit_record_view(searches[i]), cpu);
        const auto result = publication.at("result").at("result");
        check(AssayWorkspace::supports_native_publication(searches[i].at("operation").text()) &&
                  workspace->execute(searches[i], cpu).at("result").dump() == render_record_view(result).dump(),
              "sequence native admission or same-owner display failed");
        if (i == 3) {
            verify_routes(result);
            const auto overlap = NativeRecord::decode(result.at("overlap_publication").bytes());
            verify_routes(overlap, "reachability_routes");
            check(source_bytes_equal(result.at("source_bytes"), ">a\nACA\n>b\nCAC\n>c\nACA\n") &&
                      result.at("view").at("layouts").size() == 4,
                  "cover discarded repeated original reads or minimum layouts");
            for (std::size_t g = 0; g < result.at("gluing_publications").size(); ++g) {
                const auto glued = NativeRecord::decode(result.at("gluing_publications").at(g).bytes());
                check(glued.at("tokens").size() == glued.endpoint()->source.grade({0}).cell_count(),
                      "categorical gluing tokens lost their source cells");
            }
        } else {
            check(source_bytes_equal(result.at("query_bytes"), i == 1 ? comparisons : reads) && source_bytes_equal(result.at("reference_bytes"), reference) &&
                      result.at("placements").size() == result.at("view").at("matches").size(),
                  "sequence original bytes or placement couplings changed");
            for (std::size_t c = 0; c < result.at("children").size(); ++c) {
                const auto child = NativeRecord::decode(result.at("children").at(c).bytes());
                check(child.source_namespace() != result.source_namespace(), "child lost its independent source namespace");
                if (i != 1) verify_routes(child);
                else {
                    check(child.endpoint()->source.grade({1}).cell_count() == 0 &&
                              child.endpoint()->carriers.at(2).domain.logical_size() == 0,
                          "comparison lost its explicitly empty source grade");
                    const auto map_key = child.at("count_map");
                    const auto &map = field(child, map_key);
                    check(map.input == 1 && map.output == child.at("scalar_carrier").integer().get_ui() &&
                              map.tensor.support().count() == child.endpoint()->source.grade({0}).cell_count(),
                          "count map is detached from its declared source and scalar carriers");
                    for (std::size_t p = 0; p < child.at("comparisons").size(); ++p) {
                        const auto comparison = child.at("comparisons").at(p);
                        const auto &disagreement = field(child, comparison.at("disagreement")).tensor;
                        Rational sum{0};
                        for (std::size_t b = 0; b < disagreement.domain().logical_size(); ++b) {
                            sum += disagreement.state(b).value();
                            check(map.tensor.state(b).value() == 1 &&
                                      field(child, comparison.at("oriented_positions")).tensor.state(b).value() ==
                                          Rational(std::to_string(comparison.at("reverse").boolean() ? disagreement.domain().logical_size() - 1 - b : b)),
                                  "comparison map or reverse original-coordinate permutation changed");
                        }
                        check(sum == field(child, comparison.at("count")).tensor.state(0).value(),
                              "retained comparison operands disagree with native scalar output");
                    }
                }
            }
        }
        publication.write(root / ("search-" + std::to_string(i) + ".rfcxs"));
        check(NativeRecord::read(root / ("search-" + std::to_string(i) + ".rfcxs")).encode() == publication.encode(),
              "native sequence I/O changed complete archive bytes");
        search_publications.push_back(publication);
    }
    batched.write(root / "batch.rfcxs");
    auto invalid = Json::parse(R"({"schema":"gadget.assays.materialize.v1","id":"invalid","source_namespace":"730020","channel":"observation","design":{},"recipe":{}})");
    invalid.set("design", spec.at("design")); invalid.set("recipe", impossible);
    refuse([&] {workspace->materialize_workspace(admit_record_view(invalid), root / "invalid");});
    fs::rename(root / "input", root / "input-away"); workspace.reset();
    NativeRecord::decode(consensus_publication.at("result").at("workspace_publication").bytes())
        .write(root / "consensus-workspace.rfcxs");
    const AssayWorkspace consensus_workspace(root / "consensus-workspace.rfcxs");
    const auto consensus_replay = consensus_workspace.execute_record(
        NativeRecord::decode(consensus_publication.at("request_publication").bytes()), cpu);
    check(render_record_view(consensus_replay.at("result").at("result")).dump() == render_record_view(consensus_fields).dump(),
          "consensus native replay depends on original intake paths");
    NativeRecord::decode(consensus_artifact.at("result").at("workspace_publication").bytes())
        .write(root / "pipeline-consensus-workspace.rfcxs");
    const AssayWorkspace pipeline_consensus_workspace(root / "pipeline-consensus-workspace.rfcxs");
    const auto pipeline_consensus_replay = pipeline_consensus_workspace.execute_record(
        NativeRecord::decode(consensus_artifact.at("request_publication").bytes()), cpu);
    const auto before_consensus = consensus_artifact.endpoint(), after_consensus = pipeline_consensus_replay.endpoint();
    check(before_consensus->tensors.size() == after_consensus->tensors.size(),
          "pipeline consensus replay changed its native field family");
    for (std::size_t t = 0; t < before_consensus->tensors.size(); ++t)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(before_consensus->tensors.at(t), after_consensus->tensors.at(t)),
              "pipeline consensus replay changed values, support, roles or carriers");
    check(source_bytes_equal(pipeline_consensus_replay.at("result").at("result").at("source_bytes"), consensus_reads) &&
              render_record_view(pipeline_consensus_replay.at("result").at("result")).dump() == render_record_view(consensus_fields).dump(),
          "pipeline consensus source/request replay depends on former intake paths");
    NativeRecord::decode(edited_publication.at("result").at("workspace_publication").bytes())
        .write(root / "edited-splice-workspace.rfcxs");
    const AssayWorkspace edited_workspace(root / "edited-splice-workspace.rfcxs");
    const auto edited_replay = edited_workspace.execute_record(
        NativeRecord::decode(edited_publication.at("request_publication").bytes()), cpu);
    check(render_record_view(edited_replay.at("result").at("result")).dump() == render_record_view(edited_fields).dump(),
          "edited splice native replay depends on former intake paths");
    NativeRecord::decode(all_history_publication.at("result").at("workspace_publication").bytes())
        .write(root / "all-history-workspace.rfcxs");
    const AssayWorkspace history_workspace(root / "all-history-workspace.rfcxs");
    const auto history_replay = history_workspace.execute_record(
        NativeRecord::decode(all_history_publication.at("request_publication").bytes()), cpu);
    check(render_record_view(history_replay.at("result").at("result")).dump() == render_record_view(all_history).dump(),
          "all-history embedded replay depends on original intake paths");
    for (std::size_t i = 0; i < searches.size(); ++i) {
        const auto publication = NativeRecord::read(root / ("search-" + std::to_string(i) + ".rfcxs"));
        NativeRecord::decode(publication.at("result").at("workspace_publication").bytes())
            .write(root / ("search-workspace-" + std::to_string(i) + ".rfcxs"));
        const AssayWorkspace replay_workspace(root / ("search-workspace-" + std::to_string(i) + ".rfcxs"));
        const auto replay = replay_workspace.execute_record(NativeRecord::decode(publication.at("request_publication").bytes()), cpu);
        check(render_record_view(replay.at("result").at("result")).dump() ==
                  render_record_view(search_publications[i].at("result").at("result")).dump(),
              "sequence embedded replay depends on the former intake paths");
    }
    const auto relocated_batch = NativeRecord::read(root / "batch.rfcxs");
    check(relocated_batch.encode() == batched.encode() &&
              render_record_view(relocated_batch.at("result").at("result")).dump() == batch_view.dump(),
          "batch native roundtrip or relocated retained view changed");
    NativeRecord::decode(relocated_batch.at("result").at("workspace_publication").bytes())
        .write(root / "batch-workspace.rfcxs");
    const AssayWorkspace batch_workspace(root / "batch-workspace.rfcxs");
    const auto replayed = batch_workspace.execute_record(
        NativeRecord::decode(relocated_batch.at("request_publication").bytes()), cpu);
    for (std::size_t i = 0; i < batch_jobs.size(); ++i) {
        const auto child = NativeRecord::decode(
            replayed.at("result").at("result").at("results").at(i).at("child_publication").bytes());
        const auto before = children.at(i).endpoint(), after = child.endpoint();
        check(before->tensors.size() == after->tensors.size(), "batch replay changed the field family");
        for (std::size_t t = 0; t < before->tensors.size(); ++t)
            check(rfc::storage::temporal_endpoint_detail::tensor_equal(before->tensors.at(t), after->tensors.at(t)),
                  "batch replay changed field values, support, roles or carriers");
        check(render_record_view(child.at("result")).at("result").dump() ==
                  batch_view.at("results").at(i).at("answer").at("result").dump(),
              "batch replay depends on the former intake paths");
    }
    const auto relocated_motifs = NativeRecord::read(root / "motifs.rfcxs");
    check(render_record_view(relocated_motifs.at("result")).at("result").dump() ==
              motif_view.at("result").dump(), "motif view depends on original intake paths");
    for (std::size_t i = 0; i < results.size(); ++i)
        check(results[i]->replay_derivation(recipes[i].first, root / (recipes[i].first + "-replay")).at("status").text() == "PASS",
              "native replay depends on the former input directory");
    std::cout << "PASS native ordered routes, occurrence maps, optimizing factors, coupled owners and relocated replay\n";
}
