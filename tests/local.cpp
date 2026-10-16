#include "../examples/scientific/fixture.hpp"
#include <algorithm>
#include <gadget/assay.hpp>
#include <gadget/model.hpp>
#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <gadget/snapshot.hpp>
#include <gadget/workflow.hpp>
#include "../src/scientific_result_fields.hpp"
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

namespace {
using namespace gadget;
void check(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void agrees(const Json &response, const Workspace &study, const Json &query) {
    check(response.at("status").text() == "executed", "model did not execute requested query");
    const auto proposed = response.at("proposal").at("query");
    auto effective = Json::parse(proposed.dump());
    // The shipped typed-moment query defaults to first order. Preserve every
    // other requested member and compare objects without depending on key order.
    if (query.at("operation").text() == "typed_moment" && !effective.has("order"))
        effective.set("order", Json::integer(1));
    check(json_object_equal(effective.raw(), query.raw()),
          "model substituted a different observable");
    check(response.at("certificate").dump() == execute_query(study, proposed).dump(),
          "model result disagrees with canonical query certificate");
    auto explicit_certificate = execute_query(study, query);
    explicit_certificate.set("query", proposed);
    json_object_object_del(explicit_certificate.raw(), "certificate_sha256");
    explicit_certificate.set("certificate_sha256",
                             Json::string("sha256:" + sha256(explicit_certificate.dump())));
    check(json_object_equal(response.at("certificate").raw(), explicit_certificate.raw()),
          "model result differs from the explicitly requested order");
}
void evidence_role(const fs::path &model, const fs::path &evidence_path,
                   const fs::path &request_path) {
    const fs::path root = fs::path(GADGET_LOCAL_TEST_BUILD) / "evidence-role";
    fs::create_directories(root);
    const auto evidence = std::make_shared<const EvidenceStore>(evidence_path);
    const auto request = read_record(request_path).native();
    const auto context = evidence->source_context(request).at("result");
    LocalModelOptions options;
    options.model = model;
    options.threads = 2;
    options.context_tokens = 8192;
    options.output_tokens = 768;
    LocalAgent agent(evidence, request, options);
    check(!agent.prospective_publication(), "model fabricated an initial proposal");
    const auto response = agent.ask(context.at("question").text());
    write_text(root / "response.json", response.dump() + "\n");
    check(response.at("status").text() == "prospective" &&
              response.at("proposal").at("action").text() == "query",
          "evidence model did not admit the requested prospective query");
    const auto expected = Json::parse(
        R"({"operation":"affinity_evidence","endpoint":"Ki","units":"nM","maximum":"0.5"})");
    check(json_object_equal(response.at("proposal").at("query").raw(), expected.raw()),
          "model substituted a different threshold question");
    const auto publication = agent.prospective_publication();
    check(publication.has_value(), "model lost its native prospective publication");
    publication->write(root / "proposal.rfcxs");
    const auto proposed = publication->at("request").at("query");
    const auto executed = evidence->execute_record(proposed);
    executed.write(root / "result.rfcxs");
    const auto replay = NativeRecord::read(root / "proposal.rfcxs");
    const auto disabled = evidence->execute_record(replay.at("request").at("query"));
    check(json_object_equal(render_record_view(executed.at("result").at("result")).raw(),
                            render_record_view(disabled.at("result").at("result")).raw()),
          "model-disabled native query differs from retained proposal execution");
    const auto original = evidence->native_publication().raw_sources();
    const auto retained = executed.raw_sources();
    for (const auto &source : original) {
        const auto found = std::find_if(retained.begin(), retained.end(),
                                       [&](const auto &row) { return row.key() == source.key(); });
        check(found != retained.end() && std::ranges::equal(found->bytes(), source.bytes()),
              "model query changed original provider byte owners");
    }
    bool rejected = false;
    try { agent.ask("changed question"); }
    catch (const std::invalid_argument &) { rejected = true; }
    check(rejected, "model admitted an unbound question");
    check(agent.prospective_publication()->encode() == publication->encode(),
          "failed model call replaced the previously admitted publication");
    std::cout << "PASS: real local evidence planning, native prospective owner, separate native "
                 "execution, complete model-disabled equality, original bytes and replay\n";
}
void native_workflow_role(const fs::path &model, const fs::path &destination) {
    fs::create_directories(destination);
    const auto root = private_directory(destination, "native-");
    example::scientific_study(root / "input");
    auto study = std::make_shared<const ScientificStudy>(root / "input/study.json");
    LocalModelOptions options;
    options.model = model;
    options.context_tokens = 8192;
    options.output_tokens = 384;
    options.threads = 2;
    options.prompt_cache_mib = 64;
    LocalAgent agent(study, options);
    const auto field = study->catalog().at("fields").at(0).at("id").text();
    auto workflow = Json::object();
    workflow.set("schema", Json::string("gadget.scientific.workflow.v1"));
    workflow.set("id", Json::string("local-native-query"));
    workflow.set("source_digest", Json::string(study->digest()));
    workflow.set("limits", Json::parse(R"({"depth":2,"nodes":4,"model_calls":1,"parallel_workers":2})"));
    workflow.set("definitions", Json::array());
    auto entry = Json::object();
    entry.set("id", Json::string("reader")); entry.set("type", Json::string("model"));
    entry.set("role", Json::string("reader"));
    entry.set("prompt", Json::string("Return action query with operation coverage and field " + field +
        ". Do not add other arguments or compute values yourself."));
    workflow.set("entry", entry);
    std::optional<NativeRecord> response;
    WorkflowModels models{{"reader", WorkflowModel(2, WorkflowModel::NativeInvoke{
        [&](const std::string &prompt, std::shared_ptr<CpuResources> assigned) {
            response = agent.ask_assigned_record(prompt, {}, std::move(assigned));
            return *response;
        }})}};
    const auto result = execute_workflow_record(study, admit_record_view(workflow), models);
    result.write(root / "workflow.rfcxs");
    check(result.at("result").at("status").text() == "complete", "real local native workflow failed");
    check(response && response->at("result").at("status").text() == "executed" &&
              response->at("request").at("query").at("operation").text() == "coverage" &&
              response->at("request").at("query").at("field").text() == field,
          "real native local model substituted the requested source observation");
    response->write(root / "response.rfcxs");
    check(response->at("result").has("model") && response->at("result").has("performance"),
          "native local response lost model identity or execution receipts");
    const auto leaf = NativeRecord::decode(result.at("result").at("result").at("publication").bytes());
    const auto view = leaf.at("result").at("fields").at("view");
    const auto pin = study->field_accession(field);
    check(view.at("reported").integer() == pin.tensor_declaration()->tensor.support().count(),
          "native model coverage disagrees with original supported coordinates");
    leaf.write(root / "leaf.rfcxs");
    fs::rename(root / "input", root / "moved-input");
    const auto replay = replay_scientific_record(NativeRecord::read(root / "leaf.rfcxs"));
    check(replay.at("result").at("schema").text() == "gadget.scientific.result.fields.v1",
          "model-disabled relocated replay did not retain native mathematical fields");
    const auto trace = NativeRecord::decode(result.at("result").at("trace").at(0).at("result").at("publication").bytes());
    check(trace.encode() == response->encode(), "native workflow lost the exact original model response");
    std::cout << "PASS: real local model query, native workflow ownership, original support and model-disabled relocated replay\n";
}
void native_image_role(const fs::path &model, const fs::path &projector, const fs::path &image,
                       const fs::path &workspace_path, const std::string &assay,
                       const fs::path &destination) {
    fs::create_directories(destination);
    const auto root = private_directory(destination, "native-image-");
    const AssayWorkspace workspace(workspace_path);
    const auto study = workspace.study(assay);
    const auto original = read_text(image, 16u * 1024u * 1024u);
    const auto image_digest = "sha256:" + sha256(original);
    const auto catalog = workspace.catalog(); bool declared_image = false, retained_image = false;
    for (std::size_t i = 0; i < catalog.at("assays").size(); ++i) {
        const auto declaration = catalog.at("assays").at(i).at("declaration");
        if (declaration.at("id").text() != assay || !declaration.has("provenance")) continue;
        const auto provenance = declaration.at("provenance");
        for (std::size_t j = 0; j < provenance.size(); ++j)
            declared_image = declared_image || provenance.at(j).at("sha256").text() == image_digest;
    }
    for (const auto &source : workspace.native_publication().raw_sources()) {
        const auto bytes = source.bytes();
        retained_image = retained_image || (bytes.size() == original.size() &&
            std::equal(bytes.begin(),bytes.end(),reinterpret_cast<const std::byte *>(original.data())));
    }
    check(declared_image && retained_image,"image is not an original retained source of the selected assay");
    write_text(root / "original.png", original);
    LocalModelOptions options;
    options.model = model; options.projector = projector; options.task = "observer";
    options.context_tokens = 8192; options.output_tokens = 512;
    options.threads = 2; options.prompt_cache_mib = 64;
    LocalAgent observer(study, options);
    const auto response = observer.ask_record(
        "Describe one visible staining region in image 0. Include region "
        "[left,top,right,bottom] in original pixel edge coordinates. "
        "Do not infer diagnosis, cell type or gene expression.", {root / "original.png"});
    response.write(root / "response.rfcxs");
    check(response.at("status").text() == "prospective", "source image observer did not admit its proposal");
    check(response.at("media_sources").size() == 1 &&
              response.at("media_sources").at(0).at("image_index").integer() == 0,
          "source image occurrence was lost");
    const auto bytes = response.at("media_sources").at(0).at("bytes").bytes();
    check(bytes.size() == original.size() &&
              std::equal(bytes.begin(), bytes.end(), reinterpret_cast<const std::byte *>(original.data())),
          "native image input bytes changed or were replaced by a digest");
    const auto claims = response.at("certificate").at("claims");
    bool region = false;
    for (std::size_t i = 0; i < claims.size(); ++i) {
        check(claims.at(i).at("status").text() == "prospective_unverified_model_statement",
              "visual model statement became a measurement");
        if (claims.at(i).has("region") && claims.at(i).at("region").kind() != NativeRecordKind::Null) {
            check(claims.at(i).at("region_space").text() == "original_image_pixel_edges" &&
                      claims.at(i).at("image_sha256").text() == "sha256:" + sha256(original),
                  "visual region lost its original source image");
            region = true;
        }
    }
    check(region, "source image question did not retain a requested region");
    fs::rename(root / "original.png", root / "moved-original.png");
    const auto reopened = NativeRecord::read(root / "response.rfcxs");
    check(reopened.encode() == response.encode() &&
              std::ranges::equal(reopened.at("media_sources").at(0).at("bytes").bytes(), bytes),
          "relocated image response lost its original byte owner");
    std::cout << "PASS: real original source image, native input byte owner, original region and relocated replay; visual content remains unverified\n";
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 8 && std::string(argv[1]) == "--native-image") {
            native_image_role(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
            return 0;
        }
        if (argc == 4 && std::string(argv[1]) == "--native-workflow") {
            native_workflow_role(argv[2], argv[3]);
            return 0;
        }
        if (argc == 4) {
            evidence_role(argv[1], argv[2], argv[3]);
            return 0;
        }
        if (argc != 2)
            throw std::invalid_argument("local test requires one existing GGUF path");
        const fs::path root = GADGET_LOCAL_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root);
        fs::copy(GADGET_TEST_FIXTURE, root / "input", fs::copy_options::recursive);
        auto imported = load_pager(root / "input");
        import_measurements(imported, root / "input/measurements.json");
        save_snapshot(imported, root / "original", true);
        const auto original = std::make_shared<const Workspace>(load_snapshot(root / "original"));
        write_text(root / "input/effects.csv", "gene,log2FC,padj\nA,0.25,1e-12\nB,1,NA\n");
        auto updated = load_pager(root / "input");
        import_measurements(updated, root / "input/measurements.json");
        save_snapshot(updated, root / "updated", true);
        const auto next = std::make_shared<const Workspace>(load_snapshot(root / "updated"));
        LocalModelOptions options;
        options.model = argv[1];
        const auto hardware = rfc::runtime::detect_hardware();
        rfc::runtime::CpuResourcePolicy policy;
        policy.worker_limit = hardware.allowed_cpus.size();
        policy.declared_buffer_byte_limit = rfc::machine::checked_size(
            hardware.advisory_memory_headroom().value());
        options.resources = std::make_shared<CpuResources>(hardware, policy,
            [] { return rfc::runtime::detect_hardware(); });
        bool invalid_profile = false;
        try {
            LocalAgent rejected_planner(
                Json::parse(R"({"source_digest":"sha256:unused","selected_profile":"unknown"})"),
                options);
        } catch (const std::invalid_argument &) {
            invalid_profile = true;
        }
        check(invalid_profile, "pipeline planner admitted an unknown selected profile");
        LocalAgent agent(original, options);
        const std::string prompt = "Compute the first-order typed moment coupling AD.G and AD.E "
                                   "in AD using the support realization.";
        const auto query = Json::parse(
            R"({"operation":"typed_moment","context":"AD","realization":"support","field":"AD.G","other_field":"AD.E","order":1})");
        const auto first = agent.ask(prompt);
        agrees(first, *original, query);
        const auto repeated = agent.ask(prompt);
        agrees(repeated, *original, query);
        check(repeated.at("performance").at("prepared_contexts").integer() == 1,
              "repeated model query rebuilt context");
        check(repeated.at("performance").at("prefix_tokens_reused").integer() > 0,
              "repeated model query did not restore canonical prompt state");
        {
            const auto requested = Json::parse(
                R"({"operation":"green","context":"AD","realization":"support","grade":0,"field":"AD.G","lambda":"1","mode":"numerical"})");
            const std::string constrained_prompt =
                "Return action query using exactly this query without changing its parameters: " +
                requested.dump();
            const auto permitted = agent.ask(constrained_prompt);
            write_text(root / "assigned-positive.json", permitted.dump() + "\n");
            check(permitted.at("status").text() == "executed",
                  "legacy model did not execute the declared numerical Green query");
            const auto proposed = permitted.at("proposal").at("query");
            auto effective = Json::parse(proposed.dump());
            auto expected = Json::parse(requested.dump());
            // Compare exact decimal values and the existing reporting default.
            // Keep the original proposal and complete source calculation unchanged.
            effective.set("lambda", Json::string(exact_decimal(proposed.at("lambda").text()).get_str()));
            if (!effective.has("top_k")) effective.set("top_k", Json::integer(10));
            expected.set("top_k", Json::integer(10));
            check(json_object_equal(effective.raw(), expected.raw()),
                  "legacy model changed the declared numerical Green controls");
            check(permitted.at("certificate").dump() == execute_query(*original, proposed).dump(),
                  "legacy model numerical Green differs from the original source calculation");
            const auto bytes = options.resources->coordinator().available_declared_buffer_bytes();
            auto delegated = rfc::runtime::CpuResourceCoordinator::delegate(
                options.resources->coordinator().reserve(options.threads, 1));
            auto assigned = options.resources->assigned(delegated);
            bool limited = false;
            try {
                (void)agent.ask_assigned(constrained_prompt, {}, assigned);
            } catch (const rfc::runtime::ExecutionLimitExceeded &error) {
                const Execution reference(*original, "AD", "support", options.resources);
                const auto n = reference.boundary().domain_size();
                const auto required = rfc::machine::checked_size_mul(
                    rfc::machine::checked_size_mul(n, n), sizeof(double));
                check(error.resource() == "declared buffer bytes" &&
                          error.required() == required && error.limit() == 1,
                      "legacy model reported the wrong assigned resource refusal");
                limited = true;
                auto refusal = Json::object();
                refusal.set("resource", Json::string(error.resource()));
                refusal.set("required", Json::string(std::to_string(error.required())));
                refusal.set("limit", Json::string(std::to_string(error.limit())));
                write_text(root / "assigned-refusal.json", refusal.dump() + "\n");
            }
            check(limited, "legacy model bypassed its actual assigned dense byte limit");
            check(assigned->coordinator().available_workers() == options.threads &&
                      assigned->coordinator().available_declared_buffer_bytes() == 1 &&
                      options.resources->coordinator().available_declared_buffer_bytes() == bytes - 1,
                  "legacy assigned query retained or bypassed its resources");
        }
        bool rejected = false;
        try {
            agent.bind(nullptr);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        check(rejected, "invalid study switch succeeded");
        agrees(agent.ask(prompt), *original, query);
        agent.bind(next);
        const auto switched = agent.ask(prompt);
        agrees(switched, *next, query);
        check(switched.at("model").dump() == first.at("model").dump(),
              "publication switch changed/reloaded model identity");
        check(switched.at("source_digest").text() != first.at("source_digest").text(),
              "publication switch retained old source");
        check(switched.at("performance").at("prepared_contexts").integer() == 1,
              "publication switch reused preparation from a different source");
        check(switched.at("performance").at("prefix_tokens_reused").integer() == 0,
              "publication switch reused old model prompt state");
        options.cache_prompt = false;
        LocalAgent uncached(original, options);
        agrees(uncached.ask(prompt), *original, query);
        check(uncached.ask(prompt).at("performance").at("prefix_tokens_reused").integer() == 0,
              "disabled prompt cache restored state");
        options.output_tokens = 32;
        LocalAgent bounded(original, options);
        const auto incomplete = bounded.ask(prompt);
        check(incomplete.at("status").text() == "rejected" && !incomplete.has("certificate"),
              "incomplete model proposal executed");
        const auto child = ::fork();
        check(child >= 0, "model service fork test failed");
        if (child == 0) {
            try {
                agent.ask(prompt);
            } catch (const std::logic_error &) {
                ::_exit(73);
            }
            ::_exit(74);
        }
        int status{};
        check(::waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                  WEXITSTATUS(status) == 73,
              "inherited model service use admitted");
        example::scientific_study(root / "scientific");
        auto scientific = std::make_shared<const ScientificStudy>(root / "scientific/study.json");
        LocalModelOptions scientific_options;
        scientific_options.model = argv[1];
        scientific_options.context_tokens = 8192;
        LocalAgent scientific_agent(scientific, scientific_options);
        const auto inspect = scientific_agent.ask(
            "Return a query proposal with operation inspect to inspect this study catalog.");
        ScientificSession scientific_session(scientific);
        check(inspect.at("status").text() == "executed" &&
                  inspect.at("proposal").at("query").dump() ==
                      Json::parse(R"({"operation":"inspect"})").dump() &&
                  inspect.at("certificate").dump() ==
                      scientific_session.execute(Json::parse(R"({"operation":"inspect"})")).dump(),
              "general-profile model inference differs from direct scientific execution");
        bool profile_rejected = false;
        try {
            scientific_agent.bind(original);
        } catch (const std::invalid_argument &) {
            profile_rejected = true;
        }
        check(profile_rejected, "cross-profile model binding admitted");
        scientific_agent.bind_scientific(scientific);
        const auto rebound = scientific_agent.ask(
            "Return a query proposal with operation inspect to inspect this study catalog.");
        check(rebound.at("certificate").dump() == inspect.at("certificate").dump() &&
                  rebound.at("model").dump() == inspect.at("model").dump() &&
                  rebound.at("performance").at("prefix_tokens_reused").integer() == 0,
              "scientific study binding failed to retain weights or clear preparation");
        std::cout << "PASS: real local GGUF inference, requested observable and full certificate "
                     "parity, warm reuse, explicit study switch, failed switch, output bound "
                     "and general scientific inference/binding\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
