#include <charconv>
#include <gadget/assay.hpp>
#include <gadget/biomedical.hpp>
#include <gadget/execution.hpp>
#include <gadget/facts.hpp>
#include <gadget/model.hpp>
#include <gadget/pipeline.hpp>
#include <gadget/resources.hpp>
#include <gadget/store.hpp>
#include <gadget/table.hpp>
#include <gadget/workflow.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
[[maybe_unused]] void object_keys(const gadget::Json &object, const std::set<std::string> &required,
                                  const std::set<std::string> &optional = {}) {
    if (!json_object_is_type(object.raw(), json_type_object))
        throw std::invalid_argument("expected a configuration object");
    for (const auto &key : required)
        if (!object.has(key))
            throw std::invalid_argument("missing configuration member: " + key);
    json_object_object_foreach(object.raw(), key, value) {
        (void)value;
        if (!required.contains(key) && !optional.contains(key))
            throw std::invalid_argument("unknown configuration member: " + std::string(key));
    }
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout
                << "usage: gadget audit|analyze|query (--input DATA_ROOT | --snapshot STUDY | "
                   "--store DIRECTORY)\n"
                   "       snapshot (--input DATA_ROOT | --snapshot STUDY) --out NEW_DIRECTORY\n"
                   "       native-record --input FILE --out FILE.rfcxs\n"
                   "       record-view --input FILE\n"
                   "       replay-scientific --record FILE.rfcxs --out-record FILE.rfcxs\n"
                   "       native-study --study FILE --out FILE.rfcxs\n"
                   "       native-workspace --workspace FILE --out FILE.rfcxs\n"
                   "       [--measurements DECLARATION] (raw input only)\n"
                   "       analyze --out NEW_DIRECTORY; query --query FILE\n"
                   "       audit|query --study DECLARATION [--query FILE]\n"
                   "       source-index --study DECLARATION; source-read --study DECLARATION "
                   "--request FILE\n"
                   "       assay-source-read --workspace MANIFEST --request FILE\n"
                   "       plan|annotate|proposal --study DECLARATION --request FILE\n"
                   "       import-table --spec DECLARATION --out NEW_DIRECTORY; resources|devices\n"
                   "       bioinformatics; inspect-file --spec DECLARATION [--out RECORD.rfcxs]\n"
                   "       import-assays --spec DECLARATION --out NEW_DIRECTORY\n"
                   "       assay-query --workspace MANIFEST --request FILE [--out NEW_DIRECTORY | --out-record FILE.rfcxs]\n"
                   "       derive-assays --workspace MANIFEST --request FILE --out NEW_DIRECTORY\n"
                   "       replay-derivation --workspace MANIFEST --assay ID --out NEW_DIRECTORY\n"
                   "       allocation-bounds --request FILE\n"
                   "       pipelines; pipeline-prepare --workspace MANIFEST --request DESIGN\n"
                   "       pipeline --workspace MANIFEST --request PLAN --out NEW_DIRECTORY\n"
                   "       plan-pipeline --workspace MANIFEST --model LOCAL_GGUF --prompt TEXT "
                   "[--out NEW_DIRECTORY]\n"
                   "       recipe-catalog --workspace MANIFEST --intent FILE\n"
                   "       fact-catalog --workspace MANIFEST --request FILE "
                   "[--artifacts ROOT --view full|model]\n"
                   "       claim-proposal --workspace MANIFEST --facts FILE --request FILE "
                   "[--artifacts ROOT]\n"
                   "       read-facts --workspace MANIFEST --facts FILE --model LOCAL_GGUF "
                   "--prompt QUESTION [--artifacts ROOT --threads N --context-tokens N "
                   "--output-tokens N --gpu-layers N --catalog-mode compact|full]\n"
                   "       source-context --workspace MANIFEST --request FILE\n"
                   "       annotation-proposal --workspace MANIFEST --context FILE --request FILE\n"
                   "       annotate-source --workspace MANIFEST --context FILE --model LOCAL_GGUF "
                   "--prompt QUESTION [--threads N --context-tokens N --output-tokens N "
                   "--gpu-layers N]\n"
                   "       recipe-proposal --workspace MANIFEST --intent FILE --request FILE "
                   "--prompt TEXT [--out NEW_DIRECTORY]\n"
                   "       plan-recipes --workspace MANIFEST --intent FILE --model LOCAL_GGUF "
                   "--prompt TEXT [--out NEW_DIRECTORY]\n"
                   "       workflow --workspace MANIFEST --assay ID [--partition 0] --request "
                   "FILE\n"
                   "       workflow --study DECLARATION --request FILE [--models LOCAL_ROLES]\n"
                   "       databases; database-request --request FILE\n"
                   "       fetch-databases --request FILE --out NEW_DIRECTORY\n"
                   "       import-evidence --snapshot DIRECTORY --out NEW_DIRECTORY\n"
                   "       evidence-query --evidence DIRECTORY --query FILE [--out-record FILE]\n"
                   "       workflow [--evidence DIRECTORY]; ask|session [--task "
                   "planner|observer|claim_reader]\n"
                   "       compute options: [--compute-device 0] [--allow-software-device off]\n"
                   "       store-init --snapshot STUDY --out NEW_DIRECTORY\n"
                   "       store-status --store DIRECTORY [--transaction DECIMAL_ID]\n"
                   "       store-publish --store DIRECTORY --snapshot STUDY --request FILE\n";
            std::cout << "       ask|session (--snapshot STUDY | --store DIRECTORY | --study "
                         "DECLARATION) --model LOCAL_GGUF\n"
                         "       ask --prompt TEXT; session reads JSON lines {\"prompt\":\"...\"}\n"
                         "       scientific session also accepts {\"study\":\"DECLARATION\"}\n"
                         "       [--threads 8] [--context-tokens 4096] [--output-tokens 256] "
                         "[--gpu-layers 0]\n"
                         "       [--prompt-cache on|off] [--prompt-cache-mib 256] "
                         "[--catalog-mode compact|full|paged]\n";
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--version") {
            std::cout << "GADGET 0.1.0 (canonical RexGraph consumer)\n";
            return 0;
        }
        if (argc < 2)
            throw std::invalid_argument("use gadget --help for commands");
        const std::string command = argv[1];
        std::map<std::string, std::set<std::string>> allowed{
            {"databases", {}},
            {"bioinformatics", {}},
            {"pipelines", {}},
            {"pipeline-prepare", {"--workspace", "--request"}},
            {"pipeline", {"--workspace", "--request", "--out"}},
            {"pipeline-proposal", {"--workspace", "--request", "--out"}},
            {"recipe-catalog", {"--workspace", "--intent", "--profile"}},
            {"fact-catalog", {"--workspace", "--request", "--artifacts", "--view", "--out"}},
            {"claim-proposal", {"--workspace", "--facts", "--request", "--artifacts"}},
            {"read-facts",
             {"--workspace", "--facts", "--model", "--prompt", "--artifacts", "--threads",
              "--context-tokens", "--output-tokens", "--gpu-layers", "--catalog-mode"}},
            {"source-index", {"--study"}},
            {"source-read", {"--study", "--request"}},
            {"assay-source-read", {"--workspace", "--request"}},
            {"source-context", {"--workspace", "--request"}},
            {"annotation-proposal", {"--workspace", "--context", "--request"}},
            {"annotate-source",
             {"--workspace", "--context", "--model", "--prompt", "--threads", "--context-tokens",
              "--output-tokens", "--gpu-layers"}},
            {"recipe-proposal",
             {"--workspace", "--intent", "--request", "--prompt", "--out", "--profile"}},
            {"plan-recipes",
             {"--workspace", "--intent", "--model", "--prompt", "--out", "--threads",
              "--context-tokens", "--output-tokens", "--gpu-layers", "--profile"}},
            {"plan-pipeline",
             {"--workspace", "--model", "--prompt", "--out", "--threads", "--context-tokens",
              "--output-tokens", "--gpu-layers", "--profile"}},
            {"inspect-file", {"--spec", "--out"}},
            {"import-assays", {"--spec", "--out"}},
            {"assay-query", {"--workspace", "--request", "--out", "--out-record"}},
            {"derive-assays", {"--workspace", "--request", "--out"}},
            {"replay-derivation", {"--workspace", "--assay", "--out"}},
            {"allocation-bounds", {"--request"}},
            {"fetch-databases", {"--request", "--out"}},
            {"database-request", {"--request"}},
            {"import-evidence", {"--snapshot", "--out", "--adapter-version"}},
            {"evidence-query", {"--evidence", "--query", "--out-record"}},
            {"import-table", {"--spec", "--out"}},
            {"native-record", {"--input","--out"}},
            {"record-view", {"--input","--member","--format","--assay"}},
            {"replay-scientific", {"--record", "--out-record"}},
            {"native-study", {"--study","--out"}},
            {"native-workspace", {"--workspace","--out"}},
            {"resources", {}},
            {"devices", {}},
            {"workflow",
             {"--study", "--workspace", "--assay", "--partition", "--request", "--models",
              "--evidence", "--compute-device", "--allow-software-device","--out"}},
            {"audit", {"--input", "--snapshot", "--store", "--study", "--measurements"}},
            {"analyze", {"--input", "--snapshot", "--store", "--measurements", "--out"}},
            {"query", {"--input", "--snapshot", "--store", "--study", "--measurements", "--query","--out"}},
            {"plan", {"--study", "--request", "--out"}},
            {"annotate", {"--study", "--request", "--out"}},
            {"proposal", {"--study", "--request", "--out"}},
            {"snapshot", {"--input", "--snapshot", "--measurements", "--out"}},
            {"store-init", {"--snapshot", "--out"}},
            {"store-status", {"--store", "--transaction"}},
            {"store-publish", {"--store", "--snapshot", "--request"}},
            {"ask",
             {"--snapshot", "--store", "--study", "--model", "--prompt", "--threads",
              "--context-tokens", "--output-tokens", "--gpu-layers", "--prompt-cache",
              "--prompt-cache-mib", "--catalog-mode"}},
            {"session",
             {"--snapshot", "--store", "--study", "--model", "--threads", "--context-tokens",
              "--output-tokens", "--gpu-layers", "--prompt-cache", "--prompt-cache-mib",
              "--catalog-mode"}}};
        for (const auto &name : {"query", "plan", "proposal", "ask", "session"})
            for (const auto &option : {"--compute-device", "--allow-software-device"})
                allowed.at(name).insert(option);
        allowed.at("ask").insert("--image");
        for (const auto &name : {"ask", "session"})
            allowed.at(name).insert("--task");
        for (const auto &name : {"ask", "session"})
            allowed.at(name).insert("--mmproj");
        if (!allowed.contains(command))
            throw std::invalid_argument("unknown command: " + command);
        std::map<std::string, std::string> options;
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (!allowed.at(command).contains(option))
                throw std::invalid_argument("unsupported option for " + command + ": " + option);
            if (options.contains(option))
                throw std::invalid_argument("duplicate option: " + option);
            if (i + 1 >= argc)
                throw std::invalid_argument("missing argument: " + option);
            const std::string value = argv[++i];
            if (value.empty())
                throw std::invalid_argument("empty argument: " + option);
            options.emplace(option, value);
        }
        auto required = [&](const std::string &option) -> const std::string & {
            if (!options.contains(option))
                throw std::invalid_argument(command + " requires " + option);
            return options.at(option);
        };
        gadget::ComputeOptions compute;
        if (options.contains("--compute-device")) {
            const auto &s = options.at("--compute-device");
            auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), compute.device_index);
            if (ec != std::errc{} || end != s.data() + s.size() || s.empty())
                throw std::invalid_argument("invalid compute device index");
        }
        if (options.contains("--allow-software-device")) {
            const auto value = options.at("--allow-software-device");
            if (value != "on" && value != "off")
                throw std::invalid_argument("software device option requires on|off");
            compute.allow_software_device = value == "on";
        }
        if (command == "devices") {
            std::cout << gadget::compute_devices().dump() << '\n';
            return 0;
        }
        if (command == "bioinformatics") {
            std::cout << gadget::bioinformatics_capabilities().dump() << '\n';
            return 0;
        }
        if (command == "inspect-file") {
            if(options.contains("--out")) {
                const auto publication=gadget::inspect_bioinformatics_record(required("--spec"));
                publication.write(required("--out"));return 0;
            }
            std::cout << gadget::inspect_bioinformatics_file(required("--spec")).dump() << '\n';
            return 0;
        }
        if (command == "import-assays") {
            const auto workspace = gadget::import_assay_workspace(required("--spec"), required("--out"));
            std::cout << workspace->catalog().dump() << '\n';
            return 0;
        }
        if (command == "pipelines") {
            std::cout << gadget::pipeline_profiles().dump() << '\n';
            return 0;
        }
        if (command == "pipeline-prepare" || command == "pipeline" || command == "plan-pipeline" ||
            command == "pipeline-proposal" || command == "recipe-catalog" ||
            command == "recipe-proposal" || command == "plan-recipes" ||
            command == "source-context" || command == "annotation-proposal" ||
            command == "annotate-source" || command == "fact-catalog" ||
            command == "claim-proposal" || command == "read-facts") {
            auto workspace =
                std::make_shared<const gadget::AssayWorkspace>(required("--workspace"));
            gadget::Json intent;
            std::shared_ptr<const gadget::ScientificFacts> facts;
            if (command == "fact-catalog" || command == "claim-proposal" ||
                command == "read-facts") {
                const gadget::fs::path path =
                    required(command == "fact-catalog" ? "--request" : "--facts");
                const gadget::fs::path root = options.contains("--artifacts")
                                                  ? gadget::fs::path(options.at("--artifacts"))
                                                  : gadget::fs::absolute(path).parent_path();
                facts = std::make_shared<const gadget::ScientificFacts>(
                    *workspace, gadget::read_record(path).native(), root);
                if (command == "fact-catalog") {
                    if(options.contains("--out"))facts->native_publication().write(required("--out"));
                    const auto view = options.contains("--view") ? options.at("--view") : "full";
                    if (view != "full" && view != "model")
                        throw std::invalid_argument("fact catalog view must be full or model");
                    std::cout << (view == "model" ? facts->model_view() : facts->catalog()).dump()
                              << '\n';
                    return 0;
                }
                if (command == "claim-proposal") {
                    std::cout << facts
                                     ->check(gadget::Json::parse(
                                         gadget::read_text(required("--request"))))
                                     .dump()
                              << '\n';
                    return 0;
                }
            }
            if (command == "source-context") {
                std::cout << workspace
                                 ->source_context(
                                     static_cast<gadget::Json>(gadget::read_record(required("--request"))))
                                 .dump()
                          << '\n';
                return 0;
            }
            if (command == "annotation-proposal" || command == "annotate-source") {
                intent = static_cast<gadget::Json>(gadget::read_record(required("--context")));
                if (command == "annotation-proposal") {
                    std::cout << gadget::admit_source_annotations(
                                     *workspace, intent,
                                     static_cast<gadget::Json>(gadget::read_record(required("--request"))),
                                     intent.at("question").text())
                                     .dump()
                              << '\n';
                    return 0;
                }
            }
            if (command == "recipe-catalog" || command == "recipe-proposal" ||
                command == "plan-recipes") {
                intent = static_cast<gadget::Json>(gadget::read_record(required("--intent")));
                if (options.contains("--profile") &&
                    options.at("--profile") != intent.at("profile").text())
                    throw std::invalid_argument("selected profile differs from recipe intent");
                if (command == "recipe-catalog") {
                    std::cout << gadget::recipe_model_catalog(*workspace, intent).dump() << '\n';
                    return 0;
                }
            }
            if (command == "plan-pipeline" || command == "plan-recipes" ||
                command == "annotate-source" || command == "read-facts") {
#if GADGET_HAS_LOCAL_MODEL
                gadget::LocalModelOptions settings;
                settings.model = required("--model");
                auto integer = [&](const std::string &name, std::uint32_t fallback) {
                    if (!options.contains(name))
                        return fallback;
                    const auto &s = options.at(name);
                    std::uint32_t n{};
                    const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), n);
                    if (error != std::errc{} || end != s.data() + s.size())
                        throw std::invalid_argument("invalid pipeline model option: " + name);
                    return n;
                };
                settings.threads = integer("--threads", 8);
                settings.context_tokens = integer("--context-tokens", 8192);
                settings.output_tokens = integer("--output-tokens", 1024);
                const auto gpu = integer("--gpu-layers", 0);
                if (gpu > static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max()))
                    throw std::invalid_argument(
                        "GPU-layer value exceeds llama.cpp integer representation");
                settings.gpu_layers = gpu;
                if (options.contains("--catalog-mode")) {
                    const auto mode = options.at("--catalog-mode");
                    if (mode != "compact" && mode != "full" && mode != "paged")
                        throw std::invalid_argument("catalog mode must be compact or full");
                    settings.compact_catalog = mode != "full";
                    settings.paged_catalog = mode == "paged";
                }
                std::unique_ptr<gadget::LocalAgent> planner;
                if (command == "read-facts")
                    planner = std::make_unique<gadget::LocalAgent>(facts, settings);
                else if (command == "plan-recipes" || command == "annotate-source")
                    planner = std::make_unique<gadget::LocalAgent>(workspace, intent, settings);
                else {
                    auto catalog = gadget::pipeline_model_catalog(*workspace);
                    if (options.contains("--profile"))
                        catalog.set("selected_profile",
                                    gadget::Json::string(options.at("--profile")));
                    planner = std::make_unique<gadget::LocalAgent>(catalog, settings);
                }
                auto answer = planner->ask(required("--prompt"));
                const bool typed =
                    command != "plan-recipes" ||
                    (answer.has("scientific_admission") &&
                     answer.at("scientific_admission").at("status").text() == "satisfied");
                if (answer.at("status").text() == "admitted" && options.contains("--out") && typed)
                    answer.set("execution",
                               gadget::execute_pipeline(
                                   workspace, answer.at("pipeline"), required("--out"), {}, {},
                                   nullptr, command == "plan-recipes" ? &answer : nullptr));
                if (answer.at("status").text() == "admitted" && options.contains("--out") && !typed)
                    answer.set("execution_withheld",
                               gadget::Json::string("Automatic recipe execution requires a "
                                                    "researcher-declared answer contract"));
                std::cout << answer.dump() << '\n';
                return answer.at("status").text() == "rejected" ? 2 : 0;
#else
                throw std::invalid_argument("local llama.cpp support is disabled in this build");
#endif
            }
            const auto retained_declaration=gadget::read_record(required("--request"));
            const auto declaration = static_cast<gadget::Json>(retained_declaration);
            if (command == "recipe-proposal") {
                auto answer = gadget::admit_recipe_proposal(*workspace, intent, declaration,
                                                            required("--prompt"));
                if (answer.at("status").text() == "admitted" && options.contains("--out"))
                    answer.set("execution", gadget::execute_pipeline(
                                                workspace, answer.at("pipeline"), required("--out"),
                                                {}, {}, nullptr, &answer));
                std::cout << answer.dump() << '\n';
                return 0;
            }
            if (command == "pipeline-proposal") {
                auto answer = gadget::admit_pipeline_proposal(declaration, workspace->digest());
                if (answer.at("status").text() == "admitted" && options.contains("--out"))
                    answer.set("execution",
                               gadget::execute_pipeline(workspace, answer.at("pipeline"),
                                                        required("--out")));
                std::cout << answer.dump() << '\n';
                return 0;
            }
            const auto schema = declaration.at("schema").text();
            const auto plan = schema == "gadget.pipeline.default.v1"
                                  ? gadget::default_pipeline(*workspace, declaration)
                              : schema == "gadget.pipeline.design.v1"
                                  ? gadget::prepare_pipeline(*workspace, declaration)
                                  : declaration;
            if (command == "pipeline-prepare") {
                if (schema != "gadget.pipeline.design.v1" && schema != "gadget.pipeline.default.v1")
                    throw std::invalid_argument(
                        "pipeline preparation requires a design or selected defaults");
                std::cout << plan.dump() << '\n';
            } else {
                const auto native_plan=(schema=="gadget.pipeline.default.v1" || schema=="gadget.pipeline.design.v1")
                    ? retained_declaration.native().derive(gadget::import_record_view(plan)) : retained_declaration.native();
                const auto publication=gadget::execute_pipeline_record(workspace,native_plan,required("--out"));
                std::cout << gadget::render_record_view(publication.at("result")).dump() << '\n';
            }
            return 0;
        }
        if (command == "replay-derivation") {
            gadget::AssayWorkspace workspace(required("--workspace"));
            std::cout << workspace.replay_derivation(required("--assay"), required("--out")).dump()
                      << '\n';
            return 0;
        }
        if (command == "derive-assays") {
            gadget::AssayWorkspace workspace(required("--workspace"));
            std::cout << workspace
                             .materialize_record(
                                 gadget::read_record(required("--request")).native(),
                                 required("--out"))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "assay-query") {
            gadget::AssayWorkspace workspace(required("--workspace"));
            const auto request = gadget::read_record(required("--request"));
            if(options.contains("--out-record")) {
                if(options.contains("--out"))throw std::invalid_argument("choose --out or --out-record");
                workspace.execute_record(request.native()).write(required("--out-record"));
                return 0;
            }
            std::cout << (options.contains("--out")
                ? gadget::render_record_view(workspace.run_record(request.native(),required("--out")).at("result"))
                : workspace.execute(static_cast<gadget::Json>(request)))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "allocation-bounds") {
            std::cout << gadget::analyze_allocations(
                             static_cast<gadget::Json>(gadget::read_record(required("--request"))))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "databases") {
            std::cout << gadget::database_capabilities().dump() << '\n';
            return 0;
        }
        if (command == "database-request") {
            std::cout << gadget::plan_database_fetch(
                             static_cast<gadget::Json>(gadget::read_record(required("--request"))))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "fetch-databases") {
            std::cout << gadget::fetch_databases(
                             static_cast<gadget::Json>(gadget::read_record(required("--request"))),
                             required("--out"))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "import-evidence") {
            const auto version =
                options.contains("--adapter-version") ? options.at("--adapter-version") : "2";
            if (version != "1" && version != "2")
                throw std::invalid_argument("biomedical adapter version requires 1 or 2");
            gadget::import_evidence(required("--snapshot"), required("--out"),
                                    version == "1" ? 1 : 2);
            gadget::EvidenceStore store(required("--out"));
            std::cout << store.execute(gadget::Json::parse("{\"operation\":\"inspect\"}")).dump()
                      << '\n';
            return 0;
        }
        if (command == "evidence-query") {
            gadget::EvidenceStore store(required("--evidence"));
            if (options.contains("--out-record")) {
                store.execute_record(gadget::read_record(required("--query")).native())
                    .write(options.at("--out-record"));
                return 0;
            }
            std::cout
                << store.execute(static_cast<gadget::Json>(gadget::read_record(required("--query")))).dump()
                << '\n';
            return 0;
        }
        if (command == "workflow") {
            std::shared_ptr<const gadget::ScientificStudy> study;
            if (options.contains("--workspace")) {
                if (options.contains("--study"))
                    throw std::invalid_argument("workflow source must be unambiguous");
                std::size_t partition{};
                if (options.contains("--partition")) {
                    const auto &s = options.at("--partition");
                    const auto [end, ec] =
                        std::from_chars(s.data(), s.data() + s.size(), partition);
                    if (ec != std::errc{} || end != s.data() + s.size())
                        throw std::invalid_argument("invalid assay partition");
                }
                gadget::AssayWorkspace workspace(required("--workspace"));
                study = workspace.study(required("--assay"), partition);
            } else {
                if (options.contains("--assay") || options.contains("--partition"))
                    throw std::invalid_argument("assay selection requires a workspace");
                study = std::make_shared<const gadget::ScientificStudy>(required("--study"));
            }
            gadget::WorkflowModels callbacks;
#if GADGET_HAS_LOCAL_MODEL
            std::map<std::string, std::unique_ptr<gadget::LocalAgent>> agents;
            if (options.contains("--models")) {
                const auto config = static_cast<gadget::Json>(gadget::read_record(options.at("--models")));
                object_keys(config, {"schema", "models"});
                if (config.at("schema").text() != "gadget.local.roles.v1" ||
                    !json_object_is_type(config.at("models").raw(), json_type_array))
                    throw std::invalid_argument("local role configuration requires "
                                                "gadget.local.roles.v1 and a model array");
                for (std::size_t i = 0; i < config.at("models").size(); ++i) {
                    const auto item = config.at("models").at(i);
                    object_keys(item, {"role", "model"},
                                {"threads", "context_tokens", "output_tokens", "prompt_cache_mib",
                                 "projector", "images", "task"});
                    const auto role = item.at("role").text();
                    if (role.empty() || role.size() > 128 || agents.contains(role))
                        throw std::invalid_argument("duplicate/empty local model role");
                    gadget::LocalModelOptions settings;
                    if (item.has("task"))
                        settings.task = item.at("task").text();
                    settings.model = item.at("model").text();
                    if (settings.model.is_relative())
                        settings.model =
                            gadget::fs::path(options.at("--models")).parent_path() / settings.model;
                    for (const auto &key :
                         {"threads", "context_tokens", "output_tokens", "prompt_cache_mib"})
                        if (item.has(key)) {
                            const auto n = item.at(key).integer();
                            if (n < 1 || static_cast<std::uint64_t>(n) >
                                             std::numeric_limits<std::uint32_t>::max())
                                throw std::invalid_argument(
                                    "local role numeric setting does not fit its API type");
                            if (std::string(key) == "threads")
                                settings.threads = static_cast<std::uint32_t>(n);
                            else if (std::string(key) == "context_tokens")
                                settings.context_tokens = static_cast<std::uint32_t>(n);
                            else if (std::string(key) == "output_tokens")
                                settings.output_tokens = static_cast<std::uint32_t>(n);
                            else
                                settings.prompt_cache_mib = static_cast<std::uint32_t>(n);
                        }
                    settings.compute = compute;
                    if (item.has("projector")) {
                        settings.projector = item.at("projector").text();
                        if (settings.projector.is_relative())
                            settings.projector =
                                gadget::fs::path(options.at("--models")).parent_path() /
                                settings.projector;
                    }
                    std::vector<gadget::fs::path> images;
                    if (item.has("images")) {
                        if (!json_object_is_type(item.at("images").raw(), json_type_array))
                            throw std::invalid_argument("local role images must be an array");
                        for (std::size_t j = 0; j < item.at("images").size(); ++j) {
                            gadget::fs::path image = item.at("images").at(j).text();
                            if (image.is_relative())
                                image =
                                    gadget::fs::path(options.at("--models")).parent_path() / image;
                            images.push_back(image);
                        }
                    }
                    agents.emplace(role, std::make_unique<gadget::LocalAgent>(study, settings));
                    callbacks.emplace(role,
                                      gadget::WorkflowModel(
                                          settings.threads,
                                          [agent = agents.at(role).get(),
                                           images](const std::string &prompt,
                                                   std::shared_ptr<gadget::CpuResources> assigned) {
                                              return agent->ask_assigned_record(prompt, images, assigned);
                                          }));
                }
            }
#else
            if (options.contains("--models"))
                throw std::invalid_argument("local llama.cpp support is disabled in this build");
#endif
            const auto request = gadget::read_record(required("--request"));
            std::unique_ptr<gadget::EvidenceStore> evidence;
            if (options.contains("--evidence"))
                evidence = std::make_unique<gadget::EvidenceStore>(options.at("--evidence"));
            gadget::Json result;
            if (options.contains("--out")) {
                const auto publication = gadget::execute_workflow_record(study, request.native(),
                    callbacks, gadget::CpuResources::local(), compute, evidence.get());
                publication.write(required("--out"));
                result = gadget::workflow_record_view(publication);
            } else {
                result = gadget::execute_workflow(study, request.native(), callbacks,
                    gadget::CpuResources::local(), compute, evidence.get());
            }
            std::cout << result.dump() << '\n';
            return result.at("status").text() == "complete" ? 0 : 2;
        }
        if (command == "assay-source-read") {
            gadget::AssayWorkspace workspace(required("--workspace"));
            std::cout << workspace
                             .read_source(
                                 static_cast<gadget::Json>(gadget::read_record(required("--request"))))
                             .dump()
                      << '\n';
            return 0;
        }
        if (command == "import-table") {
            gadget::import_table(required("--spec"), required("--out"));
            gadget::ScientificStudy study(gadget::fs::path(required("--out")) / "study.json");
            std::cout << study.catalog().dump() << '\n';
            return 0;
        }
        if(command=="native-record") {gadget::read_record(required("--input")).native().write(required("--out"));return 0;}
        if (command == "replay-scientific") {
            const auto replayed = gadget::replay_scientific_record(
                gadget::NativeRecord::read(required("--record")), gadget::CpuResources::local());
            replayed.write(required("--out-record"));
            return 0;
        }
        if(command=="record-view") {
            const auto retained=options.contains("--assay")
                ? gadget::RetainedRecord(gadget::AssayWorkspace(required("--input")).derivation_certificate(options.at("--assay")))
                : gadget::read_record(required("--input"));
            const auto view=options.contains("--member")?static_cast<gadget::Json>(retained.at(options.at("--member"))):static_cast<gadget::Json>(retained);
            const auto format=options.contains("--format")?options.at("--format"):"pretty";
            if(format!="pretty" && format!="compact")throw std::invalid_argument("record view format must be pretty or compact");
            std::cout<<(format=="compact"?view.dump_compact():view.dump())<<'\n';return 0;
        }
        if(command=="native-study") {gadget::ScientificStudy(required("--study")).write_native(required("--out"));return 0;}
        if(command=="native-workspace") {gadget::AssayWorkspace(required("--workspace")).write_native(required("--out"));return 0;}
        if (command == "resources") {
            std::cout << gadget::CpuResources::local()->report().dump() << '\n';
            return 0;
        }
        if (command == "ask" || command == "session") {
#if GADGET_HAS_LOCAL_MODEL
            if (options.contains("--snapshot") + options.contains("--store") +
                    options.contains("--study") !=
                1)
                throw std::invalid_argument(
                    "local model commands require exactly one of --snapshot, --store or --study");
            gadget::LocalModelOptions settings;
            if (options.contains("--task"))
                settings.task = options.at("--task");
            settings.model = required("--model");
            settings.compute = compute;
            if (options.contains("--mmproj"))
                settings.projector = options.at("--mmproj");
            auto integer = [&](const std::string &key, std::uint32_t fallback) {
                if (!options.contains(key))
                    return fallback;
                const auto &text = options.at(key);
                std::uint32_t value{};
                auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
                if (error != std::errc{} || end != text.data() + text.size() || text.empty() ||
                    (text.size() > 1 && text[0] == '0'))
                    throw std::invalid_argument("expected unsigned decimal option: " + key);
                return value;
            };
            settings.threads = integer("--threads", 8);
            settings.context_tokens = integer("--context-tokens", 4096);
            settings.output_tokens = integer("--output-tokens", 256);
            settings.prompt_cache_mib = integer("--prompt-cache-mib", 256);
            if (options.contains("--catalog-mode")) {
                const auto &mode = options.at("--catalog-mode");
                if (mode != "compact" && mode != "full" && mode != "paged")
                    throw std::invalid_argument("--catalog-mode expects compact, full or paged");
                settings.compact_catalog = mode != "full";
                settings.paged_catalog = mode == "paged";
            }
            if (options.contains("--prompt-cache")) {
                const auto &value = options.at("--prompt-cache");
                if (value != "on" && value != "off")
                    throw std::invalid_argument("--prompt-cache expects on or off");
                settings.cache_prompt = value == "on";
            }
            const auto gpu_layers = integer("--gpu-layers", 0);
            if (gpu_layers > 999)
                throw std::invalid_argument("GPU-layer bound exceeded");
            settings.gpu_layers = gpu_layers;
            const auto prompt = command == "ask" ? required("--prompt") : std::string{};
#if GADGET_HAS_DURABLE_STORE
            std::unique_ptr<gadget::ConcurrentStudyStore> store;
#endif
            std::shared_ptr<const gadget::Workspace> source;
            std::unique_ptr<gadget::LocalAgent> agent;
            if (options.contains("--study"))
                agent = std::make_unique<gadget::LocalAgent>(
                    std::make_shared<const gadget::ScientificStudy>(options.at("--study")),
                    settings);
            else if (options.contains("--snapshot"))
                source = std::make_shared<const gadget::Workspace>(
                    gadget::load_snapshot(options.at("--snapshot")));
            else {
#if GADGET_HAS_DURABLE_STORE
                store = std::make_unique<gadget::ConcurrentStudyStore>(options.at("--store"));
                source = store->snapshot().workspace();
#else
                throw std::invalid_argument("durable store support is disabled in this build");
#endif
            }
            if (!agent)
                agent =
                    std::make_unique<gadget::LocalAgent>(std::move(source), std::move(settings));
            if (command == "ask") {
                auto result = options.contains("--image")
                                  ? agent->ask(prompt, {options.at("--image")})
                                  : agent->ask(prompt);
                std::cout << result.dump() << '\n';
                return result.at("status").text() == "rejected" ? 2 : 0;
            }
            while (true) {
                std::string line;
                bool consumed = false, oversized = false;
                char character;
                while (std::cin.get(character)) {
                    consumed = true;
                    if (character == '\n')
                        break;
                    if (line.size() < gadget::source_byte_limit())
                        line += character;
                    else
                        oversized = true;
                }
                if (!consumed)
                    break;
                gadget::Json result;
                try {
                    if (oversized)
                        throw std::invalid_argument("session input byte budget exceeded");
                    const auto request = gadget::Json::parse(line);
                    if (request.has("action") &&
                        request.at("action").text() == "clear_source_context") {
                        object_keys(request, {"action"});
                        agent->clear_source_context();
                        result = gadget::Json::object();
                        result.set("schema",
                                   gadget::Json::string("gadget.source.context.cleared.v1"));
                        result.set("status", gadget::Json::string("cleared"));
                        std::cout << result.dump_compact() << '\n' << std::flush;
                        continue;
                    }
                    if (!json_object_is_type(request.raw(), json_type_object) ||
                        (!request.has("prompt") && !request.has("study")))
                        throw std::invalid_argument(
                            "session requires a prompt or scientific study path per line");
                    object_keys(request,
                                request.has("study") ? std::set<std::string>{"study"}
                                                     : std::set<std::string>{"prompt"},
                                request.has("study") ? std::set<std::string>{}
                                                     : std::set<std::string>{"images"});
                    if (request.has("study")) {
                        if (!options.contains("--study"))
                            throw std::invalid_argument(
                                "study switching requires a scientific --study session");
                        auto next = std::make_shared<const gadget::ScientificStudy>(
                            request.at("study").text());
                        agent->bind_scientific(next);
                        result = gadget::Json::object();
                        result.set("schema", gadget::Json::string("gadget.study.bound.v1"));
                        result.set("status", gadget::Json::string("bound"));
                        result.set("source_digest", gadget::Json::string(next->digest()));
                        result.set("study_id", next->source_index().at("study_id"));
                    } else {
                        std::vector<gadget::fs::path> images;
                        if (request.has("images")) {
                            const auto list = request.at("images");
                            if (!json_object_is_type(list.raw(), json_type_array))
                                throw std::invalid_argument("session images must be an array");
                            for (std::size_t i = 0; i < list.size(); ++i)
                                images.emplace_back(list.at(i).text());
                        }
                        result = agent->ask(request.at("prompt").text(), images);
                    }
                } catch (const std::exception &error) {
                    result = gadget::Json::object();
                    result.set("schema", gadget::Json::string("gadget.agent.response.v1"));
                    result.set("status", gadget::Json::string("rejected"));
                    result.set("error", gadget::Json::string(error.what()));
                }
                std::cout << json_object_to_json_string_ext(result.raw(), JSON_C_TO_STRING_PLAIN)
                          << '\n'
                          << std::flush;
            }
            if (std::cin.bad())
                throw std::runtime_error("session input read failed");
            return 0;
#else
            throw std::invalid_argument("local llama.cpp support is disabled in this build");
#endif
        }
        if (options.contains("--study")) {
            if (options.contains("--input") || options.contains("--snapshot") ||
                options.contains("--store") || options.contains("--measurements"))
                throw std::invalid_argument(
                    "scientific declaration must be the sole admitted source");
            auto study = std::make_shared<const gadget::ScientificStudy>(required("--study"));
            gadget::ScientificSession session(study, {}, compute);
            gadget::Json result;
            if (command == "source-index")
                result = study->source_index();
            else if (command == "source-read")
                result = study->read_source(
                    static_cast<gadget::Json>(gadget::read_record(required("--request"))));
            else if (command == "audit")
                result = study->catalog();
            else if (command == "query") {
                const auto publication=session.execute_record(gadget::read_record(required("--query")).native());
                if(options.contains("--out"))publication.write(required("--out"));
                result=gadget::render_record_view(publication.at("result"));
            }
            else {
                const auto request=gadget::read_record(required("--request"));
                if (command == "plan") {
                    const auto publication=session.execute_plan_record(request.native());
                    if(options.contains("--out"))publication.write(required("--out"));
                    result=gadget::render_record_view(publication.at("result"));
                } else if (command == "proposal" && request.native().has("action") &&
                           request.native().at("action").kind() == gadget::NativeRecordKind::Text &&
                           (request.native().at("action").text() == "query" ||
                            request.native().at("action").text() == "plan")) {
                    const auto publication = session.apply_proposal_record(request.native());
                    if (options.contains("--out")) publication.write(required("--out"));
                    result = gadget::render_record_view(publication.at("result"));
                } else {
                    result=command=="annotate"
                        ?session.propose_annotations(static_cast<gadget::Json>(request))
                        :gadget::apply_scientific_proposal(session,static_cast<gadget::Json>(request));
                    if(options.contains("--out"))
                        gadget::retain_operation(study->native_publication(),request.native(),
                            gadget::import_record_view(result)).write(required("--out"));
                }
            }
            std::cout << result.dump() << '\n';
            return result.has("status") && result.at("status").text() == "rejected" ? 2 : 0;
        }
        if (command.starts_with("store-")) {
#if GADGET_HAS_DURABLE_STORE
            gadget::Json result;
            if (command == "store-init")
                result = gadget::initialize_store(required("--snapshot"), required("--out"));
            else if (command == "store-publish")
                result = gadget::publish_store(
                    required("--store"), required("--snapshot"),
                    static_cast<gadget::Json>(gadget::read_record(required("--request"))));
            else {
                std::optional<std::uint64_t> transaction;
                if (options.contains("--transaction"))
                    transaction = gadget::store_transaction(options.at("--transaction"));
                result = gadget::store_status(required("--store"), transaction);
            }
            std::cout << result.dump() << '\n';
            return 0;
#else
            throw std::invalid_argument("durable store support is disabled in this build");
#endif
        }
        if (options.contains("--input") + options.contains("--snapshot") +
                options.contains("--store") !=
            1)
            throw std::invalid_argument("choose exactly one of --input, --snapshot or --store");
        if (options.contains("--measurements") && !options.contains("--input"))
            throw std::invalid_argument("import measurements with --input before snapshot capture");
        if (command == "analyze" || command == "snapshot")
            required("--out");
        if (command == "query")
            required("--query");
        auto source = [&] {
            if (options.contains("--input"))
                return gadget::load_pager(options.at("--input"));
            if (options.contains("--snapshot"))
                return gadget::load_snapshot(options.at("--snapshot"));
#if GADGET_HAS_DURABLE_STORE
            return gadget::load_store(options.at("--store"));
#else
            throw std::invalid_argument("durable store support is disabled in this build");
#endif
        }();
        if (options.contains("--measurements"))
            gadget::import_measurements(source, options.at("--measurements"));
        if (command == "analyze")
            gadget::analyze(source, options.at("--out"));
        else if (command == "snapshot")
            std::cout << gadget::save_snapshot(source, options.at("--out")).dump() << '\n';
        else if (command == "audit") {
            auto result = source.manifest();
            if (source.store_origin.has("root"))
                result.set("store", source.store_origin);
            std::cout << result.dump() << '\n';
        } else
            std::cout << gadget::execute_query(
                             source, static_cast<gadget::Json>(gadget::read_record(options.at("--query"))))
                             .dump()
                      << '\n';
        return 0;
#if GADGET_HAS_DURABLE_STORE
    } catch (const rfc::storage::DurableUnknownOutcome &e) {
        auto result = gadget::Json::object();
        result.set("schema", gadget::Json::string("gadget.store.unknown.v1"));
        result.set("transaction", gadget::Json::string(std::to_string(e.transaction())));
        result.set("status", gadget::Json::string("unknown"));
        result.set(
            "action",
            gadget::Json::string(
                e.transaction()
                    ? "Reopen with store-status for this transaction; retain the identical request "
                      "and snapshot for retry."
                    : "Reopen the existing directory with store-status; never overwrite genesis."));
        std::cerr << result.dump() << '\n';
        return 3;
#endif
    } catch (const std::exception &e) {
        std::cerr << "gadget: " << e.what() << '\n';
        return 2;
    }
}
