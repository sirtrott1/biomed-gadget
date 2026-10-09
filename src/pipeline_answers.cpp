#include "pipeline_answers.hpp"
#include <gadget/record_view.hpp>
#include <algorithm>
#include <limits>
#include <optional>
#include <rfc/core/builder.hpp>
#include <rfc/local/constraint_system.hpp>
#include <set>

namespace gadget::pipeline_detail {
namespace {
void require(bool value, const std::string &message) {
    if (!value)
        throw std::invalid_argument(message);
}
void keys(const Json &value, const std::set<std::string> &required) {
    require(json_object_is_type(value.raw(), json_type_object), "answer declaration object");
    for (const auto &key : required)
        require(value.has(key), "missing answer declaration member: " + key);
    json_object_object_foreach(value.raw(), key, member) {
        (void)member;
        require(required.contains(key), "unknown answer declaration member: " + std::string(key));
    }
}
std::string text(const Json &value) {
    const auto out = value.text();
    require(!out.empty() && out.find('\0') == std::string::npos,
            "answer declaration requires nonempty text without NUL");
    return out;
}
struct Observable {
    std::vector<std::string> methods;
    std::string meaning;
};
// Output meanings of reviewed, shipped methods. Purposes and model text cannot
// amend this registry. The complete bound goal owns the scientific scope.
const std::map<std::string, Observable> observables{
    {"read_quality",
     {{"quality"}, "Retained Phred+33 scores, cycle histograms and marginal IUPAC GC counts"}},
    {"alignment_records",
     {{"alignment_evidence"},
      "Selected source alignment records and retained CIGAR junction occurrences"}},
    {"annotation_features",
     {{"annotation_features"},
      "Declared annotation groups, source records and physical interval unions"}},
    {"feature_allocations",
     {{"feature_counts"},
      "Annotation-compatible segment candidates and sharp class-allocation bounds"}},
    {"allocation_ranges",
     {{"allocation_bounds"},
      "Sharp marginal target and contrast ranges over declared nonnegative class allocations"}},
    {"expression_allocation_ranges",
     {{"expression_allocations"},
      "Read-supported target and contrast ranges under declared scales and optional coupling "
      "constraints"}},
    {"field_support", {{"coverage"}, "Reported and missing coordinates of the declared fields"}},
    {"paired_field_contrast",
     {{"expression_contrast"},
      "Full matched-field differences under declared pairings and scales"}},
    {"controlled_field_relations",
     {{"controlled_relations"}, "Declared control-span decomposition and exact field pairings"}},
    {"linear_model_comparison",
     {{"linear_model_comparison"},
      "Joint field projections, complete coefficient fibres, nested term comparisons and estimable "
      "contrasts"}},
    {"finite_reference_calibration",
     {{"finite_reference_calibration"},
      "Exact inclusive tails, joint maximum-statistic comparisons and weighted variation under a "
      "declared finite reference law"}},
    {"structural_reference_calibration",
     {{"structural_reference_calibration"},
      "Source-checked exchange laws, conditional structural tails and joint measurement scenario "
      "sensitivity"}},
    {"observed_response_model",
     {{"observed_response_model"},
      "Observed count/exposure and ordinal fields, exact model coefficient fibres, full residual "
      "sectors, measured dose contrasts and declared joint finite-law calibration"}},
    {"constrained_response_model",
     {{"constrained_response_model"},
      "Physically admissible represented counts and joint ordinal cuts, authored dose directions, "
      "complete coefficient fibres and exact primal/dual optimality certificates"}},
    {"source_structure_ablation",
     {{"source_structure_ablation"},
      "Source-bound role ablations, strict inclusions, relative homology, persistent "
      "representative "
      "images, complete sector and Green changes and declared sheaf restrictions"}},
    {"relational_field_ensemble",
     {{"relational_field_ensemble"},
      "Complete source field sectors, Dirac and Green actions, weighted spread, sector coupling "
      "and declared group contrasts"}},
    {"hodge_action",
     {{"hodge"}, "Full declared metric Hodge operator action on the supplied field"}},
    {"green_action", {{"green_action"}, "Exact shifted Green field on the declared carrier"}},
    {"local_sections", {{"local_sections"}, "Compatibility under declared local restrictions"}},
    {"relational_gluing",
     {{"relational_gluing"}, "Declared source pushout, complete occurrence maps and compatible scalar sections"}},
    {"assay_sections",
     {{"assay_sections"}, "Compatible and residual assay fields with assessed-support witnesses"}},
    {"linear_aggregation",
     {{"aggregation", "design_aggregation"}, "Full declared feature and observation linear maps"}},
    {"interval_coverage",
     {{"interval_coverage"},
      "Complete physical coverage atoms, length integral and declared threshold components"}},
    {"chromatin_contrast",
     {{"chromatin_contrast"},
      "Matched sample and input coverage differences on a shared physical reference"}},
    {"controlled_coverage",
     {{"controlled_coverage"},
      "Full physical contrast, control-span projection, residual and original interval carrier"}},
    {"spatial_relations",
     {{"spatial_relations"},
      "All physical neighbor relations within the declared exact metric quadrance"}},
    {"field_transfer",
     {{"field_transfer"},
      "Full mapped field, residual, exact pairing and ambiguity of the authored inter-assay map"}},
    {"genotype_observations",
     {{"genotypes"},
      "Reported allele-copy counts and incomplete-call support from retained genotype records"}},
    {"motif_matches",
     {{"motifs"},
      "Every definite and possible source window under declared IUPAC patterns and strands"}},
    {"global_edit_routes",
     {{"exact_sequence_alignment"},
      "Exact declared global edit objective, minimizing-step union and complete prefix-state "
      "carrier"}},
    {"reference_placements",
     {{"reference_sequence_matches"},
      "Every admitted ungapped or authored-junction placement and source mismatch witness"}},
    {"indexed_edit_placements",
     {{"indexed_sequence_alignment"},
      "Every placement within the declared unit-edit bound, with a native minimum route for each "
      "placement"}},
    {"splice_chains",
     {{"sequence_splice_discovery"},
      "Every exact exon chain permitted by the declared motifs, exon lengths and junction bounds"}},
    {"kmer_structure",
     {{"sequence_assembly_graph"},
      "All observed literal k-mer occurrences, source-read walks and retained branching carrier"}},
    {"read_cover_reconstructions",
     {{"read_overlap_reconstruction"},
      "All minimum suffix-prefix read-cover layouts and spellings with complete source-read "
      "gluing"}},
    {"constrained_allocations",
     {{"constrained_allocations"},
      "Feasibility, sharp marginal and optimal-face ranges in the declared bounded integer or "
      "rational constraint domain"}},
    {"chain_structure",
     {{"chain_certificate"}, "Retained boundary tower and exact integral chain structure"}},
    {"realized_path_comparison",
     {{"path_comparison"},
      "Complete declared path fields, endpoint boundaries and shared-field metric form"}},
    {"typed_moment",
     {{"typed_moment"}, "Exact pairing in the declared carrier, channels and metric"}},
    {"finite_moment_change",
     {{"temporal_moment_change"}, "Declared finite operator, observation and interaction terms"}},
    {"temporal_field_history",
     {{"temporal_field_history"},
      "Publication-bound ordered operation defects, transported innovations and common-target "
      "interactions"}},
    {"censored_endpoint_bounds",
     {{"censored_endpoint_bounds"},
      "Sharp marginal finite-population event fractions and restricted-time bounds over declared "
      "observation intervals"}},
    {"structure_transition",
     {{"structure_transition"}, "Declared transport, sector and Green response changes"}},
    {"declared_response_transition",
     {{"declared_response_transition"},
      "Source-bound realization, innovation, resolving and observation changes"}},
    {"hodge_frame_decomposition",
     {{"hodge_frame_decomposition"},
      "Complete source-grade gradient, curl and harmonic fields with independent frame "
      "coordinates"}},
    {"linear_response_observation",
     {{"linear_response_observation"},
      "Declared linear observations of source-bound operator solves"}},
    {"metric_relative_novelty",
     {{"metric_relative_novelty"},
      "Metric-orthogonal change within the declared old image, its complement and complete "
      "ordered change-field pairings"}},
    {"green_response", {{"green_response"}, "Declared multi-order Green responses"}},
    {"relational_persistence",
     {{"relational_persistence"}, "Classes and induced maps through declared retained contexts"}},
    {"population_p_value",
     {{}, "Requires an authored population, sampling/null family and an implemented calibration"}}};
Json strings(const std::vector<std::string> &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(Json::string(value));
    return out;
}
Json hashed(Json value, const std::string &key) {
    value.set(key, Json::string("sha256:" + sha256(value.dump())));
    return value;
}
Json reference(const std::string &path, const std::string &digest, const Json &certificate) {
    require(certificate.has("result") && certificate.has("query") &&
                certificate.has("source_digest") && certificate.has("certificate_sha256"),
            "answer requires a complete mathematical certificate");
    auto out = Json::object();
    out.set("path", Json::string(path));
    out.set("sha256", Json::string(digest));
    out.set("pointer", Json::string("/result"));
    out.set("source_digest", certificate.at("source_digest"));
    out.set("query", certificate.at("query"));
    out.set("certificate_sha256", certificate.at("certificate_sha256"));
    return out;
}
} // namespace
fs::path checked_answer_path(const fs::path &directory, const Json &record) {
    const fs::path path(record.at("path").text());
    require(!path.empty() && !path.is_absolute(), "answer artifact must be relative");
    for (const auto &part : path)
        require(part != ".." && part != ".", "answer artifact path traversal");
    const auto resolved = fs::canonical(directory / path), root = fs::canonical(directory);
    const auto relative = resolved.lexically_relative(root);
    require(!relative.empty() && *relative.begin() != "..", "answer artifact escapes owner");
    return resolved;
}
std::string checked_answer_bytes(const fs::path &directory, const Json &record) {
    const auto resolved = checked_answer_path(directory, record);
    const auto bytes = read_text(resolved);
    require(record.at("sha256").text() == "sha256:" + sha256(bytes),
            "answer artifact hash mismatch");
    return bytes;
}
Json checked_answer_view(const fs::path &directory,const Json &record) {
    if (record.has("native_selection")) {
        const auto selection = record.at("native_selection");
        keys(selection, {"kind", "assay"});
        require(selection.at("kind").text() == "generating_result",
                "unsupported native answer selection");
        const auto resolved = checked_answer_path(directory, record);
        const AssayWorkspace workspace(resolved);
        const auto history = *workspace.native_publication().publication_history();
        const auto digest = "sha256:" + sha256_chunks([&](const auto &sink) {
                                rfc::storage::stream_native_series(history, sink);
                            });
        require(record.at("sha256").text() == digest, "native answer publication hash mismatch");
        return workspace.derivation_certificate(selection.at("assay").text());
    }
    const auto bytes = checked_answer_bytes(directory, record);
    const auto retained = decode_record(std::span<const std::byte>{
        reinterpret_cast<const std::byte *>(bytes.data()), bytes.size()});
    if (retained.has("schema") && retained.at("schema").text() == "gadget.operation.native.v1")
        return retained.at("result");
    return static_cast<Json>(retained);
}
Json observable_catalog() {
    auto out = Json::array();
    for (const auto &[id, observable] : observables) {
        auto entry = Json::object();
        entry.set("id", Json::string(id));
        entry.set("methods", strings(observable.methods));
        entry.set("meaning", Json::string(observable.meaning));
        entry.set("available", Json::boolean(!observable.methods.empty()));
        out.push(entry);
    }
    return out;
}
Json answer_admission(const Json &intent, const Json &catalog, const Json *included,
                      const std::string *question) {
    auto out = Json::object();
    out.set("schema", Json::string("gadget.pipeline.answer.admission.v1"));
    out.set("source_digest", catalog.at("source_digest"));
    out.set("intent_sha256", catalog.at("intent_sha256"));
    if (!intent.has("answer_contract")) {
        out.set("status", Json::string("untyped"));
        out.set("authority",
                Json::string("Recipe provenance only; question coverage is undeclared"));
        return out;
    }
    const auto contract = intent.at("answer_contract");
    keys(contract, {"schema", "question", "requirements"});
    require(text(contract.at("schema")) == "gadget.pipeline.answer.contract.v1",
            "answer contract schema");
    const auto declared_question = text(contract.at("question"));
    if (question)
        require(*question == declared_question,
                "question differs from the researcher answer contract");
    const auto requests = contract.at("requirements");
    require(json_object_is_type(requests.raw(), json_type_array) && requests.size() > 0,
            "answer requirements must be a nonempty array");
    require(requests.size() <= std::numeric_limits<std::uint32_t>::max(),
            "answer requirements exceed the native cell representation");
    std::map<std::string, Json> bindings;
    for (std::size_t i = 0; i < catalog.at("recipes").size(); ++i) {
        const auto binding = catalog.at("recipes").at(i);
        bindings.emplace(binding.at("id").text(), binding);
    }
    std::set<std::string> selection;
    if (included)
        for (std::size_t i = 0; i < included->size(); ++i)
            selection.insert(included->at(i).text());
    using namespace rfc;
    core::RelationSignature signature;
    core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < requests.size(); ++i)
        builder.declare_witness_identity(i + 1);
    const auto carrier = builder.freeze();
    local::FiniteConstraintSystem<std::uint8_t> constraints(carrier);
    local::FiniteConstraintSystem<std::uint8_t>::Assignment mandatory;
    std::optional<core::CellId> missing;
    algebra::SubcomplexMask full(carrier);
    auto requirements = Json::array();
    std::set<std::string> identities;
    for (std::size_t i = 0; i < requests.size(); ++i) {
        const auto request = requests.at(i);
        keys(request, {"id", "observable", "bindings"});
        const auto id = text(request.at("id")), observable = text(request.at("observable"));
        require(identities.insert(id).second, "duplicate answer requirement");
        require(observables.contains(observable), "undeclared answer observable: " + observable);
        const auto &definition = observables.at(observable);
        const auto candidates = request.at("bindings");
        require(json_object_is_type(candidates.raw(), json_type_array), "answer bindings array");
        auto eligible = Json::array(), supplied = Json::array(), witnesses = Json::array();
        std::set<std::string> unique;
        for (std::size_t b = 0; b < candidates.size(); ++b) {
            const auto name = text(candidates.at(b));
            require(unique.insert(name).second && bindings.contains(name),
                    "answer requires distinct, declared recipe bindings");
            const auto binding = bindings.at(name);
            const auto method = binding.at("method").text();
            if (std::find(definition.methods.begin(), definition.methods.end(), method) ==
                definition.methods.end())
                continue;
            eligible.push(Json::string(name));
            auto witness = Json::object();
            witness.set("binding", Json::string(name));
            witness.set("goal_sha256", binding.at("goal_sha256"));
            witness.set("method", binding.at("method"));
            witnesses.push(witness);
            if (!included || selection.contains(name))
                supplied.push(Json::string(name));
        }
        const core::CellId cell{core::Grade{0}, static_cast<std::uint32_t>(i)};
        constraints.set_boolean_domain(cell);
        mandatory.emplace(cell, std::uint8_t{1});
        if (!supplied.size() && !missing)
            missing = cell;
        algebra::SubcomplexMask patch(carrier);
        patch.select(cell);
        full.select(cell);
        constraints.add_chart("require/" + id, patch, {cell}, {{1}});
        constraints.add_chart("supply/" + id, patch, {cell},
                              {{static_cast<std::uint8_t>(supplied.size() != 0)}});
        auto row = Json::parse(request.dump());
        row.set("eligible_bindings", eligible);
        row.set("supplied_bindings", supplied);
        row.set("scope_witnesses", witnesses);
        row.set("meaning", Json::string(definition.meaning));
        row.set("fulfilled", Json::boolean(supplied.size() != 0));
        requirements.push(row);
    }
    // Every authored require chart fixes its variable to one. Validate that
    // unique possible section directly, avoiding a redundant recursive search.
    // A local unsatisfiable chart pair certifies failure of the global family.
    const bool valid = constraints.globally_valid(mandatory);
    auto native = Json::object(), charts = Json::array(), assignment = Json::array();
    native.set("owner", Json::string("rfc::local::FiniteConstraintSystem<uint8_t>"));
    native.set("carrier",
               Json::string("Separate grade-zero answer obligations; no biological incidence"));
    native.set("cells", Json::integer(requests.size()));
    native.set("constraint_complete",
               Json::boolean(constraints.cover_is_constraint_complete({full})));
    native.set("global_section_exists", Json::boolean(valid));
    for (const auto &chart : constraints.charts()) {
        auto row = Json::object();
        row.set("name", Json::string(chart.name));
        row.set("cell", Json::integer(chart.support.at(0).index));
        row.set("allowed", Json::integer(chart.allowed.at(0).at(0)));
        charts.push(row);
    }
    native.set("charts", charts);
    native.set("validation", Json::string("mandatory section and exact local conflict core"));
    if (valid) {
        for (std::size_t i = 0; i < requests.size(); ++i)
            assignment.push(
                Json::integer(mandatory.at({core::Grade{0}, static_cast<std::uint32_t>(i)})));
        native.set("assignment", assignment);
    } else {
        require(missing.has_value(), "native answer failure requires an unsupplied obligation");
        algebra::SubcomplexMask patch(carrier);
        patch.select(*missing);
        const auto core = constraints.exact_patch_unsat_core(patch);
        require(core.has_value() && core->subset_minimal, "native answer conflict certificate");
        auto conflict = Json::object(), names = Json::array(), cells = Json::array();
        for (const auto index : core->chart_indices)
            names.push(Json::string(constraints.charts().at(index).name));
        for (const auto cell : core->cells)
            cells.push(Json::integer(cell.index));
        conflict.set("charts", names);
        conflict.set("cells", cells);
        conflict.set("subset_minimal", Json::boolean(core->subset_minimal));
        native.set("unsat_core", conflict);
    }
    out.set("status", Json::string(valid ? "satisfied" : "unsatisfied"));
    out.set("question_sha256", Json::string("sha256:" + sha256(declared_question)));
    out.set("question_bytes", Json::integer(declared_question.size()));
    out.set("contract_sha256", Json::string("sha256:" + sha256(contract.dump())));
    out.set("authority", Json::string("Researcher-declared observables and bound goal parameters"));
    out.set("requirements", requirements);
    out.set("native", native);
    return hashed(out, "admission_sha256");
}
Json answer_witnesses(const Json &planning, const Json &pipeline, const Json &results,
                      const fs::path &directory) {
    const auto admission = planning.at("scientific_admission");
    require(admission.at("status").text() == "satisfied", "execution answer admission");
    std::map<std::string, Json> receipts;
    for (std::size_t i = 0; i < results.size(); ++i)
        receipts.emplace(results.at(i).at("id").text(), results.at(i));
    auto answers = Json::array();
    const auto requirements = admission.at("requirements"), nodes = pipeline.at("steps");
    for (std::size_t i = 0; i < requirements.size(); ++i) {
        const auto requirement = requirements.at(i);
        auto answer = Json::object(), artifacts = Json::array();
        answer.set("id", requirement.at("id"));
        answer.set("observable", requirement.at("observable"));
        answer.set("scope_witnesses", requirement.at("scope_witnesses"));
        std::set<std::string> supplied;
        for (std::size_t b = 0; b < requirement.at("supplied_bindings").size(); ++b)
            supplied.insert(requirement.at("supplied_bindings").at(b).text());
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            const auto node = nodes.at(n);
            if (!node.has("recipe_binding") || !supplied.contains(node.at("recipe_binding").text()))
                continue;
            const auto receipt = receipts.at(node.at("id").text());
            auto record = Json::object();
            record.set("path", receipt.at("result"));
            record.set("sha256", receipt.at("result_sha256"));
            const auto result = checked_answer_view(directory,record);
            auto effective_record = Json::object();
            effective_record.set("path", receipt.at("request"));
            effective_record.set("sha256", receipt.at("request_sha256"));
            const auto effective = checked_answer_view(directory,effective_record);
            const auto job = directory / node.at("id").text();
            auto references = Json::array();
            if (node.at("type").text() == "derive") {
                const AssayWorkspace workspace(job / "workspace/workspace.rfcxs");
                require(workspace.digest() == receipt.at("derived_digest").text(),
                        "answer derived workspace owner mismatch");
                if (workspace.native_publication().at("schema").text() ==
                    "gadget.assay.workspace.native.v2") {
                    const auto assay = effective.at("id").text();
                    const auto certificate = workspace.derivation_certificate(assay);
                    require(certificate.at("source_digest").text() ==
                                    receipt.at("source_digest").text() &&
                                certificate.at("query").dump() == effective.at("recipe").dump(),
                            "native generating certificate owner/request mismatch");
                    const auto path = job / "workspace/workspace.rfcxs";
                    auto ref = reference(node.at("id").text() + "/workspace/workspace.rfcxs",
                                         "sha256:" + sha256_file(path), certificate);
                    auto selection = Json::object();
                    selection.set("kind", Json::string("generating_result"));
                    selection.set("assay", Json::string(assay));
                    ref.set("native_selection", selection);
                    references.push(ref);
                } else {
                    const auto original =
                        workspace.catalog().at("assays").at(0).at("originals").at("provenance-0");
                    const auto certificate = checked_answer_view(job / "workspace", original);
                    require(certificate.at("schema").text() == "gadget.assay.result.v1" &&
                                certificate.at("source_digest").text() ==
                                    receipt.at("source_digest").text() &&
                                certificate.at("query").dump() == effective.at("recipe").dump(),
                            "answer generating certificate owner/request mismatch");
                    references.push(
                        reference(node.at("id").text() + "/workspace/" + original.at("path").text(),
                                  original.at("sha256").text(), certificate));
                }
            } else if (effective.at("operation").text() == "native_query") {
                require(result.at("schema").text() == "gadget.assays.run.v1" &&
                            result.at("source_digest").text() ==
                                receipt.at("source_digest").text() &&
                            result.at("query").dump() == effective.dump(),
                        "answer native result owner/request mismatch");
                for (std::size_t p = 0; p < result.at("results").size(); ++p) {
                    const auto native = result.at("results").at(p);
                    const auto certificate =
                        checked_answer_view(job/"results",native);
                    require(certificate.at("source_digest").text() ==
                                    native.at("source_digest").text() &&
                                certificate.at("query").dump() ==
                                    result.at("query").at("query").dump(),
                            "answer native certificate owner/request mismatch");
                    references.push(
                        reference(node.at("id").text() + "/results/" + native.at("path").text(),
                                  native.at("sha256").text(), certificate));
                }
            } else {
                require(result.at("schema").text() == "gadget.assay.result.v1" &&
                            result.at("source_digest").text() ==
                                receipt.at("source_digest").text() &&
                            result.at("query").dump() == effective.dump(),
                        "answer certificate owner/request mismatch");
                references.push(
                    reference(record.at("path").text(), record.at("sha256").text(), result));
            }
            require(references.size() > 0, "answer method emitted no result artifacts");
            auto artifact = Json::object();
            artifact.set("binding", node.at("recipe_binding"));
            artifact.set("step", node.at("id"));
            artifact.set("source_digest", receipt.at("source_digest"));
            artifact.set("effective_request_sha256", receipt.at("request_sha256"));
            artifact.set("results", references);
            artifacts.push(artifact);
        }
        require(artifacts.size() > 0, "answer obligation has no completed method output");
        answer.set("artifacts", artifacts);
        if (planning.has("source_context")) {
            const auto mappings = planning.at("source_context").at("requirements");
            for (std::size_t m = 0; m < mappings.size(); ++m)
                if (mappings.at(m).at("id").text() == requirement.at("id").text())
                    answer.set("source_mentions", mappings.at(m).at("mentions"));
        }
        answers.push(answer);
    }
    auto out = Json::object();
    out.set("schema", Json::string("gadget.pipeline.answers.v1"));
    out.set("source_digest", planning.at("source_digest"));
    out.set("intent_sha256", planning.at("intent_sha256"));
    out.set("question_sha256", admission.at("question_sha256"));
    out.set("admission_sha256", admission.at("admission_sha256"));
    out.set("answers", answers);
    if (planning.has("source_context"))
        out.set("source_context", planning.at("source_context"));
    return hashed(out, "answers_sha256");
}
} // namespace gadget::pipeline_detail
