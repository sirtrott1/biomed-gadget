#include "pipeline_answers.hpp"
#include <algorithm>
#include <chrono>
#include <gadget/pipeline.hpp>
#include <gadget/record_view.hpp>
#include <mutex>
#include <optional>
#include <rfc/runtime/work_graph.hpp>
#include <set>

namespace gadget {
namespace {
void require(bool value, const std::string &message) {
    if (!value)
        throw std::invalid_argument(message);
}
void keys(const Json &v, const std::set<std::string> &required,
          const std::set<std::string> &optional = {}) {
    require(json_object_is_type(v.raw(), json_type_object), "pipeline member must be an object");
    for (const auto &key : required)
        require(v.has(key), "missing pipeline member: " + key);
    json_object_object_foreach(v.raw(), key, value) {
        (void)value;
        require(required.contains(key) || optional.contains(key),
                "unknown pipeline member: " + std::string(key));
    }
}
std::string identifier(const Json &v) {
    const auto s = v.text();
    require(!s.empty() && s.size() <= 128 && s != "source" &&
                std::all_of(s.begin(), s.end(),
                            [](unsigned char c) {
                                return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                                       (c >= '0' && c <= '9') || c == '_' || c == '-';
                            }),
            "invalid pipeline identity");
    return s;
}
Json sealed(Json value, const std::string &key) {
    value.set(key, Json::string("sha256:" + sha256(value.dump())));
    return value;
}
const std::map<std::string, std::vector<std::string>> profiles{
    {"rna_seq",
     {"bind_assays", "design_aggregation", "read_quality", "alignment_evidence",
      "annotation_features", "feature_counts", "allocation_bounds", "expression_allocations",
      "constrained_allocations", "relational_persistence", "native_query"}},
    {"single_cell",
     {"bind_assays", "design_aggregation", "assay_sections", "relational_gluing", "aggregate_assay",
      "relational_persistence", "native_query"}},
    {"atac_seq",
     {"bind_assays", "design_aggregation", "read_quality", "alignment_evidence",
      "interval_coverage", "chromatin_contrast", "controlled_coverage", "sequence_motifs",
      "relational_persistence", "native_query"}},
    {"chip_seq",
     {"bind_assays", "design_aggregation", "read_quality", "alignment_evidence",
      "interval_coverage", "chromatin_contrast", "controlled_coverage", "sequence_motifs",
      "relational_persistence", "native_query"}},
    {"spatial",
     {"bind_assays", "design_aggregation", "spatial_relations", "relational_gluing", "aggregate_assay", "field_transfer",
      "relational_persistence", "native_query"}},
    {"multi_omics",
     {"bind_assays", "design_aggregation", "assay_sections", "relational_gluing", "field_transfer", "aggregate_assay",
      "relational_persistence", "native_query"}},
    {"variants",
     {"bind_assays", "design_aggregation", "variant_genotypes", "relational_persistence",
      "native_query"}},
    {"sequence",
     {"bind_assays", "read_quality", "sequence_motifs", "exact_sequence_alignment",
      "reference_sequence_matches", "indexed_sequence_alignment", "sequence_splice_discovery",
      "sequence_assembly_graph", "read_overlap_reconstruction", "relational_persistence",
      "native_query"}},
    {"custom", {}}};
struct Recipe {
    std::string operation;
    bool native{}, materializable{}, binding{};
};
// Recipes assemble existing contracts. Scientific operators and source
// materialization remain in their existing RexGraph-backed services.
const std::map<std::string, Recipe> recipes{
    {"quality", {"read_quality"}},
    {"coverage", {"coverage", true}},
    {"alignment_evidence", {"alignment_evidence"}},
    {"annotation_features", {"annotation_features"}},
    {"feature_counts", {"feature_counts", false, true}},
    {"allocation_bounds", {"allocation_bounds"}},
    {"expression_allocations", {"expression_allocations", false, true}},
    {"expression_contrast", {"paired_field_contrasts", true}},
    {"controlled_relations", {"controlled_pairings", true}},
    {"linear_model_comparison", {"linear_model_comparison", true}},
    {"finite_reference_calibration", {"finite_reference_calibration", true}},
    {"structural_reference_calibration", {"structural_reference_calibration", true}},
    {"observed_response_model", {"observed_response_model", true}},
    {"constrained_response_model", {"constrained_response_model", true}},
    {"source_structure_ablation", {"source_structure_ablation", true}},
    {"relational_field_ensemble", {"relational_field_ensemble", true}},
    {"aggregation", {"aggregate_assay", false, true}},
    {"design_aggregation", {"design_aggregation", false, true}},
    {"bind_assays", {"bind_assays", false, false, true}},
    {"interval_coverage", {"interval_coverage"}},
    {"chromatin_contrast", {"chromatin_contrast", false, true}},
    {"controlled_coverage", {"controlled_coverage", false, true}},
    {"spatial_relations", {"spatial_relations", false, true}},
    {"field_transfer", {"field_transfer"}},
    {"assay_sections", {"assay_sections", false, true}},
    {"relational_gluing", {"relational_gluing"}},
    {"genotypes", {"variant_genotypes", false, true}},
    {"motifs", {"sequence_motifs"}},
    {"hodge", {"hodge", true}},
    {"green_action", {"green_action", true}},
    {"local_sections", {"local_sections", true}},
    {"exact_sequence_alignment", {"exact_sequence_alignment", false, true}},
    {"reference_sequence_matches", {"reference_sequence_matches"}},
    {"indexed_sequence_alignment", {"indexed_sequence_alignment"}},
    {"sequence_splice_discovery", {"sequence_splice_discovery"}},
    {"sequence_assembly_graph", {"sequence_assembly_graph", false, true}},
    {"read_overlap_reconstruction", {"read_overlap_reconstruction"}},
    {"constrained_allocations", {"constrained_allocations", false, true}},
    {"relational_persistence", {"relational_persistence"}},
    {"chain_certificate", {"chain_certificate", true}},
    {"path_comparison", {"realized_path_comparison", true}},
    {"typed_moment", {"typed_moment", true}},
    {"temporal_moment_change", {"temporal_moment_change", true}},
    {"temporal_field_history", {"temporal_field_history", true}},
    {"censored_endpoint_bounds", {"censored_endpoint_bounds", true}},
    {"declared_response_transition", {"declared_response_transition", true}},
    {"hodge_frame_decomposition", {"hodge_frame_decomposition", true}},
    {"linear_response_observation", {"linear_response_observation", true}},
    {"metric_relative_novelty", {"metric_relative_novelty", true}},
    {"structure_transition", {"structure_transition", true}},
    {"green_response", {"green_response", true}}};
bool allowed_recipe(const std::string &profile, const Recipe &recipe) {
    if (recipe.operation == "paired_field_contrasts" && profile != "custom" &&
        profile != "rna_seq" && profile != "single_cell" && profile != "spatial" &&
        profile != "multi_omics")
        return false;
    const auto &allowed = profiles.at(profile);
    const auto operation = recipe.native ? "native_query" : recipe.operation;
    return allowed.empty() || std::find(allowed.begin(), allowed.end(), operation) != allowed.end();
}
void evidence_binding(const Json &node) {
    if (!node.has("normalization_evidence_field"))
        return;
    const auto request = node.at("request");
    require(
        node.at("type").text() == "query" && request.at("operation").text() == "native_query" &&
            request.at("query").at("operation").text() == "paired_field_contrasts" &&
            !request.at("query").at("normalization").has("evidence") &&
            !node.at("normalization_evidence_field").text().empty(),
        "normalization evidence binding requires an expression contrast without supplied evidence");
}
std::set<std::string> binding_inputs(const Json &request) {
    keys(request, {"schema", "id", "design", "assays", "resources", "links", "method", "evidence"});
    require(request.at("schema").text() == "gadget.assays.bindings.v1", "bindings request schema");
    std::set<std::string> inputs;
    for (const auto *axis : {"assays", "resources"}) {
        const auto entries = request.at(axis);
        require(json_object_is_type(entries.raw(), json_type_array), "bindings must be arrays");
        for (std::size_t i = 0; i < entries.size(); ++i) {
            keys(entries.at(i),
                 {"input", axis == std::string("assays") ? "assay" : "resource", "id"});
            inputs.insert(entries.at(i).at("input").text());
        }
    }
    require(!inputs.empty(), "bindings require selected inputs");
    return inputs;
}
void validate(const Json &plan, const std::string &digest) {
    keys(plan, {"schema", "source_digest", "profile", "question", "steps"},
         {"concurrent_jobs", "plan_sha256"});
    require(plan.at("schema").text() == "gadget.pipeline.v1", "pipeline schema");
    require(plan.at("source_digest").text() == digest, "pipeline source digest mismatch");
    require(profiles.contains(plan.at("profile").text()), "unknown pipeline profile");
    require(!plan.at("question").text().empty() &&
                plan.at("question").text().size() <= source_byte_limit(),
            "pipeline question");
    const auto steps = plan.at("steps");
    require(json_object_is_type(steps.raw(), json_type_array) && steps.size() > 0,
            "pipeline requires selected steps");
    std::map<std::string, std::string> declared{{"source", "workspace"}};
    for (std::size_t i = 0; i < steps.size(); ++i) {
        const auto node = steps.at(i);
        keys(node, {"id", "type", "input", "request"},
             {"depends_on", "assay", "partition", "workers", "normalization_evidence_field",
              "recipe_binding"});
        if (node.has("recipe_binding"))
            (void)identifier(node.at("recipe_binding"));
        if (node.has("workers"))
            require(node.at("workers").integer() > 0, "pipeline workers must be positive");
        const auto id = identifier(node.at("id")), type = node.at("type").text(),
                   input = node.at("input").text();
        require(!declared.contains(id), "duplicate pipeline step");
        require(declared.contains(input) && declared.at(input) == "workspace",
                "pipeline input must name source or an earlier derived workspace");
        require(type == "query" || type == "derive" || type == "workflow" || type == "bind",
                "unsupported pipeline step type");
        if (node.has("depends_on")) {
            require(json_object_is_type(node.at("depends_on").raw(), json_type_array),
                    "pipeline dependencies must be an array");
            std::set<std::string> dependencies;
            for (std::size_t j = 0; j < node.at("depends_on").size(); ++j) {
                const auto dependency = node.at("depends_on").at(j).text();
                require(dependency != "source" && declared.contains(dependency) &&
                            dependencies.insert(dependency).second,
                        "pipeline dependency must name a distinct earlier step");
            }
        }
        const auto request = node.at("request");
        require(json_object_is_type(request.raw(), json_type_object),
                "pipeline request must be an object");
        evidence_binding(node);
        if (type == "bind") {
            require(input == "source" && !node.has("assay") && !node.has("partition") &&
                        !node.has("normalization_evidence_field"),
                    "bind step selects inputs inside its binding request");
            for (const auto &selected : binding_inputs(request))
                require(declared.contains(selected) && declared.at(selected) == "workspace",
                        "binding must select source or an earlier workspace-producing step");
            const auto &allowed = profiles.at(plan.at("profile").text());
            require(allowed.empty() ||
                        std::find(allowed.begin(), allowed.end(), "bind_assays") != allowed.end(),
                    "bindings do not belong to selected assay profile");
        } else if (type == "workflow") {
            require(node.has("assay"), "workflow step requires an assay");
            require(request.at("schema").text() == "gadget.scientific.workflow.v1",
                    "pipeline workflow schema");
        } else {
            require(!node.has("assay") && !node.has("partition"),
                    "assay/partition belong inside query requests");
            require(
                request.at("schema").text() ==
                    (type == "derive" ? "gadget.assays.materialize.v1" : "gadget.assays.query.v1"),
                "pipeline request schema");
            const auto op = type == "derive" ? request.at("recipe").at("operation").text()
                                             : request.at("operation").text();
            require(op != "query_batch",
                    "use pipeline dependencies and concurrency for query batches");
            const auto &allowed = profiles.at(plan.at("profile").text());
            require(allowed.empty() ||
                        std::find(allowed.begin(), allowed.end(), op) != allowed.end(),
                    "operation does not belong to selected assay profile");
        }
        declared.emplace(id, type == "derive" || type == "bind" ? "workspace" : "certificate");
    }
    if (plan.has("concurrent_jobs"))
        require(plan.at("concurrent_jobs").integer() > 0, "pipeline concurrency must be positive");
    if (plan.has("plan_sha256")) {
        auto unsigned_plan = Json::object();
        json_object_object_foreach(plan.raw(), key, value) {
            (void)value;
            if (std::string(key) != "plan_sha256")
                unsigned_plan.set(key, plan.at(key));
        }
        require(plan.at("plan_sha256").text() == "sha256:" + sha256(unsigned_plan.dump()),
                "pipeline plan seal mismatch");
    }
}
struct WorkspaceLease {
    std::mutex mutex;
    std::shared_ptr<const AssayWorkspace> workspace;
    std::size_t consumers;
    WorkspaceLease(std::shared_ptr<const AssayWorkspace> source, std::size_t users)
        : workspace(std::move(source)), consumers(users) {}
    std::shared_ptr<const AssayWorkspace> acquire() {
        std::lock_guard lock(mutex);
        require(workspace && consumers, "pipeline workspace has no remaining consumer");
        auto source = workspace;
        if (!--consumers) workspace.reset();
        return source;
    }
};
struct Artifact {
    std::shared_ptr<WorkspaceLease> workspace;
    RetainedRecord receipt;
};
} // namespace

Json pipeline_profiles() {
    auto out = Json::object(), entries = Json::array(), definitions = Json::array();
    out.set("schema", Json::string("gadget.pipeline.profiles.v1"));
    for (const auto &[name, methods] : profiles) {
        auto entry = Json::object(), operations = Json::array(), defaults = Json::array();
        entry.set("profile", Json::string(name));
        for (const auto &method : methods)
            operations.push(Json::string(method));
        entry.set("operations", operations);
        for (const auto &[method, recipe] : recipes)
            if (allowed_recipe(name, recipe))
                defaults.push(Json::string(method));
        entry.set("recipes", defaults);
        entry.set("selection",
                  Json::string("Only the researcher's selected steps execute; assay-specific "
                               "parameters and mappings must be declared"));
        entries.push(entry);
    }
    for (const auto &[method, recipe] : recipes) {
        auto choice = Json::object();
        choice.set("method", Json::string(method));
        choice.set("operation", Json::string(recipe.operation));
        choice.set("native_query", Json::boolean(recipe.native));
        choice.set("materializable", Json::boolean(recipe.materializable));
        choice.set("native_materialization", Json::boolean(
            AssayWorkspace::supports_native_materialization(recipe.operation)));
        choice.set("workspace_binding", Json::boolean(recipe.binding));
        choice.set("parameters_schema",
                   Json::string(recipe.binding  ? "assays-bindings.schema.json"
                                : recipe.native ? "scientific-query.schema.json"
                                                : "assays-query.schema.json"));
        definitions.push(choice);
    }
    out.set("observables", pipeline_detail::observable_catalog());
    out.set("profiles", entries);
    out.set("recipes", definitions);
    return out;
}
Json prepare_pipeline(const AssayWorkspace &workspace, const Json &design) {
    keys(design, {"schema", "profile", "question", "steps"}, {"source_digest", "concurrent_jobs"});
    require(design.at("schema").text() == "gadget.pipeline.design.v1", "pipeline design schema");
    auto plan = Json::parse(design.dump());
    plan.set("schema", Json::string("gadget.pipeline.v1"));
    if (design.has("source_digest"))
        require(design.at("source_digest").text() == workspace.digest(),
                "pipeline design source mismatch");
    plan.set("source_digest", Json::string(workspace.digest()));
    validate(plan, workspace.digest());
    return sealed(plan, "plan_sha256");
}
Json default_pipeline(const AssayWorkspace &workspace, const Json &settings) {
    keys(settings, {"schema", "profile", "question", "goals"},
         {"assay", "fields", "reads", "quality_thresholds", "concurrent_jobs"});
    require(settings.at("schema").text() == "gadget.pipeline.default.v1",
            "default pipeline schema");
    const auto profile = settings.at("profile").text();
    require(profiles.contains(profile), "unknown pipeline profile");
    auto design = Json::object(), steps = Json::array();
    design.set("schema", Json::string("gadget.pipeline.design.v1"));
    design.set("profile", settings.at("profile"));
    design.set("question", settings.at("question"));
    if (settings.has("concurrent_jobs"))
        design.set("concurrent_jobs", settings.at("concurrent_jobs"));
    const auto goals = settings.at("goals");
    require(json_object_is_type(goals.raw(), json_type_array), "default goals must be an array");
    std::set<std::string> selected;
    std::map<std::string, std::string> derived;
    for (std::size_t i = 0; i < goals.size(); ++i) {
        auto goal = goals.at(i);
        if (json_object_is_type(goal.raw(), json_type_string)) {
            const auto method = goal.text();
            goal = Json::object();
            goal.set("method", Json::string(method));
            auto parameters = Json::object();
            if (method == "quality") {
                require(settings.has("reads") && settings.has("quality_thresholds"),
                        "quality goal requires reads and explicit quality thresholds");
                parameters.set("resource", settings.at("reads"));
                parameters.set("quality_thresholds", settings.at("quality_thresholds"));
            } else {
                require(method == "coverage", "parameterized recipe requires a goal object");
                if (settings.has("fields"))
                    parameters.set("fields", settings.at("fields"));
            }
            goal.set("parameters", parameters);
        }
        keys(goal, {"method", "parameters"},
             {"id", "input", "assay", "output", "depends_on", "workers", "partitions",
              "partition_jobs", "normalization_evidence_field"});
        const auto method = goal.at("method").text();
        require(recipes.contains(method), "unknown default recipe: " + method);
        const auto &recipe = recipes.at(method);
        require(allowed_recipe(profile, recipe),
                "recipe does not belong to selected assay profile");
        const auto id = identifier(goal.has("id") ? goal.at("id") : goal.at("method"));
        require(selected.insert(id).second, "duplicate default pipeline goal identity");
        const auto input = goal.has("input") ? goal.at("input").text() : "source";
        require(input == "source" || derived.contains(input),
                "recipe input must name source or an earlier materialized goal");
        auto parameters = Json::parse(goal.at("parameters").dump());
        require(json_object_is_type(parameters.raw(), json_type_object) &&
                    !parameters.has("schema") && !parameters.has("operation"),
                "recipe parameters contain method arguments, without schema or operation");
        auto node = Json::object(), request = Json::object();
        node.set("id", Json::string(id));
        node.set("recipe_binding", Json::string(id));
        node.set("type", Json::string("query"));
        node.set("input", Json::string(input));
        for (const auto *member : {"depends_on", "workers", "normalization_evidence_field"})
            if (goal.has(member))
                node.set(member, goal.at(member));
        if (recipe.binding) {
            require(input == "source" && !goal.has("output") && !goal.has("assay") &&
                        !goal.has("partitions") && !goal.has("partition_jobs") &&
                        !goal.has("normalization_evidence_field"),
                    "binding recipe declares its selected inputs and workspace in parameters");
            request = parameters;
            request.set("schema", Json::string("gadget.assays.bindings.v1"));
            for (const auto &selected_input : binding_inputs(request))
                require(selected_input == "source" || derived.contains(selected_input),
                        "bindings require earlier workspace-producing goals");
            node.set("type", Json::string("bind"));
            derived.emplace(id, std::string{});
        } else if (recipe.native) {
            std::string assay;
            if (goal.has("assay"))
                assay = goal.at("assay").text();
            else if (input != "source")
                assay = derived.at(input);
            else if (settings.has("assay"))
                assay = settings.at("assay").text();
            else {
                const auto assays = workspace.catalog().at("assays");
                require(assays.size() == 1,
                        "recipe requires an assay when the source contains multiple assays");
                assay = assays.at(0).at("declaration").at("id").text();
            }
            require(!assay.empty(),
                    "native query on a bound collection requires an explicit assay");
            request.set("schema", Json::string("gadget.assays.query.v1"));
            request.set("operation", Json::string("native_query"));
            request.set("assay", Json::string(assay));
            if (goal.has("partitions"))
                request.set("partitions", goal.at("partitions"));
            if (goal.has("partition_jobs"))
                request.set("concurrent_jobs", goal.at("partition_jobs"));
            if (method == "coverage") {
                keys(parameters, {}, {"fields"});
                auto fields = parameters.has("fields") ? parameters.at("fields") : Json::array();
                if (!parameters.has("fields")) {
                    require(input == "source",
                            "derived coverage requires explicit field identities");
                    const auto declared = workspace.study(assay)->declaration().at("fields");
                    for (std::size_t f = 0; f < declared.size(); ++f)
                        fields.push(declared.at(f).at("id"));
                }
                require(json_object_is_type(fields.raw(), json_type_array) && fields.size() > 0,
                        "coverage goal requires measured fields");
                std::set<std::string> unique;
                for (std::size_t f = 0; f < fields.size(); ++f) {
                    require(unique.insert(fields.at(f).text()).second, "duplicate coverage field");
                    auto query = Json::object();
                    query.set("operation", Json::string("coverage"));
                    query.set("field", fields.at(f));
                    auto selected_request = Json::parse(request.dump());
                    selected_request.set("query", query);
                    auto selected_node = Json::parse(node.dump());
                    selected_node.set("request", selected_request);
                    if (fields.size() > 1)
                        selected_node.set("id", Json::string(id + "-" + std::to_string(f)));
                    steps.push(selected_node);
                }
                continue;
            }
            parameters.set("operation", Json::string(recipe.operation));
            request.set("query", parameters);
        } else {
            require(!goal.has("assay") && !goal.has("partitions") && !goal.has("partition_jobs"),
                    "workspace recipe assay belongs in parameters; partition options require a "
                    "native recipe");
            request = parameters;
            request.set("schema", Json::string("gadget.assays.query.v1"));
            request.set("operation", Json::string(recipe.operation));
        }
        if (goal.has("output")) {
            require(recipe.materializable, "recipe cannot publish a derived assay");
            const auto output = goal.at("output");
            keys(output, {"id", "source_namespace", "design", "channel"}, {"view", "component"});
            auto producing = Json::parse(output.dump());
            producing.set("schema", Json::string("gadget.assays.materialize.v1"));
            producing.set("recipe", request);
            request = producing;
            node.set("type", Json::string("derive"));
            derived.emplace(id, output.at("id").text());
        }
        node.set("request", request);
        steps.push(node);
    }
    design.set("steps", steps);
    return prepare_pipeline(workspace, design);
}
Json admit_pipeline_proposal(const Json &proposal, const std::string &source_digest) {
    keys(proposal, {"schema", "action", "source_digest"},
         {"profile", "question", "steps", "reason"});
    require(proposal.at("schema").text() == "gadget.pipeline.proposal.v1" &&
                proposal.at("source_digest").text() == source_digest,
            "pipeline proposal source/schema mismatch");
    auto out = Json::object();
    if (proposal.at("action").text() == "clarify") {
        require(proposal.has("reason") && !proposal.has("steps") && !proposal.has("profile") &&
                    !proposal.has("question"),
                "clarification requires only a reason");
        require(!proposal.at("reason").text().empty(), "empty pipeline clarification");
        out.set("status", Json::string("clarification"));
        out.set("reason", proposal.at("reason"));
    } else {
        require(proposal.at("action").text() == "pipeline" && !proposal.has("reason"),
                "pipeline proposal action");
        auto plan = Json::object(), steps = Json::array();
        plan.set("schema", Json::string("gadget.pipeline.v1"));
        plan.set("source_digest", Json::string(source_digest));
        plan.set("profile", proposal.at("profile"));
        plan.set("question", proposal.at("question"));
        for (std::size_t i = 0; i < proposal.at("steps").size(); ++i) {
            const auto entry = proposal.at("steps").at(i);
            keys(entry, {"id", "type", "input"},
                 {"request", "request_json", "assay", "partition", "depends_on", "workers"});
            require(entry.has("request") != entry.has("request_json"),
                    "pipeline proposal requires one structured or encoded request");
            auto node = Json::object();
            for (const auto &key :
                 {"id", "type", "input", "assay", "partition", "depends_on", "workers"})
                if (entry.has(key))
                    node.set(key, entry.at(key));
            node.set("request", entry.has("request")
                                    ? Json::parse(entry.at("request").dump())
                                    : Json::parse(entry.at("request_json").text()));
            steps.push(node);
        }
        plan.set("steps", steps);
        validate(plan, source_digest);
        out.set("status", Json::string("admitted"));
        out.set("pipeline", sealed(plan, "plan_sha256"));
    }
    out.set("proposal", proposal);
    return out;
}
Json pipeline_model_catalog(const AssayWorkspace &workspace) {
    auto out = Json::object(), assays = Json::array(), resources = Json::array();
    const auto catalog = workspace.catalog();
    for (std::size_t i = 0; i < catalog.at("assays").size(); ++i) {
        const auto a = catalog.at("assays").at(i);
        auto item = Json::object();
        item.set("declaration", a.at("declaration"));
        const auto field_catalog =
            workspace.study(a.at("declaration").at("id").text())->catalog().at("fields");
        auto fields = Json::array();
        for (std::size_t j = 0; j < std::min<std::size_t>(field_catalog.size(), 32); ++j)
            fields.push(field_catalog.at(j));
        item.set("fields", fields);
        item.set("field_catalog_complete", Json::boolean(field_catalog.size() <= 32));
        item.set("partitions", Json::integer(a.at("partitions").size()));
        assays.push(item);
    }
    for (std::size_t i = 0; i < catalog.at("resources").size(); ++i)
        resources.push(catalog.at("resources").at(i).at("declaration"));
    out.set("source_digest", Json::string(workspace.digest()));
    // The current model decoder proposes query steps, not default-recipe goal
    // objects. Present that executable interface; recipe configurations have
    // their own schema and must not be confused with workspace operations.
    const auto declared_profiles = pipeline_profiles();
    auto query_profiles = Json::object(), choices = Json::array();
    query_profiles.set("schema", declared_profiles.at("schema"));
    for (std::size_t i = 0; i < declared_profiles.at("profiles").size(); ++i) {
        auto choice = Json::object();
        const auto profile = declared_profiles.at("profiles").at(i);
        for (const auto *member : {"profile", "selection"})
            choice.set(member, profile.at(member));
        auto operations = Json::array();
        for (std::size_t j = 0; j < profile.at("operations").size(); ++j)
            if (profile.at("operations").at(j).text() != "bind_assays")
                operations.push(profile.at("operations").at(j));
        choice.set("operations", operations);
        choices.push(choice);
    }
    query_profiles.set("profiles", choices);
    out.set("profiles", query_profiles);
    out.set("assays", assays);
    out.set("resources", resources);
    return out;
}
namespace {
struct RecipeIntent {
    RetainedRecord declaration, defaults, catalog;
    std::vector<std::string> names;
    std::map<std::string, std::size_t> indices;
    std::vector<std::vector<std::size_t>> dependencies;
};
RecipeIntent prepare_intent(const AssayWorkspace &workspace, const Json &declaration) {
    RecipeIntent prepared;
    prepared.declaration = Json::parse(declaration.dump());
    const auto &bound = prepared.declaration;
    keys(bound, {"schema", "profile", "recipes"},
         {"assay", "concurrent_jobs", "answer_contract", "source_context"});
    require(bound.at("schema").text() == "gadget.pipeline.intent.v1" &&
                profiles.contains(bound.at("profile").text()),
            "recipe intent schema/profile");
    const auto entries = bound.at("recipes");
    require(entries.native().kind()==NativeRecordKind::Array, "intent recipes must be an array");
    prepared.defaults = Json::object();
    prepared.defaults.set("schema", Json::string("gadget.pipeline.default.v1"));
    prepared.defaults.set("profile", bound.at("profile"));
    prepared.defaults.set("question", Json::string("Validate declared recipe bindings"));
    for (const auto *key : {"assay", "concurrent_jobs"})
        if (bound.has(key))
            prepared.defaults.set(key, bound.at(key));
    auto goals = Json::array(), choices = Json::array();
    std::map<std::string, std::size_t> earlier;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto entry = entries.at(i);
        keys(entry, {"id", "purpose", "goal"});
        const auto id = identifier(entry.at("id"));
        require(!earlier.contains(id) && !entry.at("purpose").text().empty(),
                "recipe binding requires a unique identity and purpose");
        auto goal = Json::parse(entry.at("goal").dump());
        require(!goal.has("id"), "recipe goal identity is supplied by its binding");
        goal.set("id", Json::string(id));
        std::vector<std::size_t> dependencies;
        auto declared_dependencies = Json::array();
        auto dependency = [&](const std::string &name) {
            require(earlier.contains(name), "recipe dependency must name an earlier binding");
            if (std::find(dependencies.begin(), dependencies.end(), earlier.at(name)) ==
                dependencies.end()) {
                dependencies.push_back(earlier.at(name));
                declared_dependencies.push(Json::string(name));
            }
        };
        if (goal.has("input") && goal.at("input").text() != "source")
            dependency(goal.at("input").text());
        if (goal.at("method").text() == "bind_assays") {
            auto request = Json::parse(goal.at("parameters").dump());
            request.set("schema", Json::string("gadget.assays.bindings.v1"));
            for (const auto &input : binding_inputs(request))
                if (input != "source")
                    dependency(input);
        }
        if (goal.has("depends_on")) {
            require(json_object_is_type(goal.at("depends_on").raw(), json_type_array),
                    "recipe dependencies must be an array");
            for (std::size_t j = 0; j < goal.at("depends_on").size(); ++j)
                dependency(goal.at("depends_on").at(j).text());
        }
        auto choice = Json::object();
        choice.set("id", Json::string(id));
        choice.set("purpose", entry.at("purpose"));
        choice.set("method", goal.at("method"));
        choice.set("dependencies", declared_dependencies);
        choice.set("goal_sha256", Json::string("sha256:" + sha256(entry.at("goal").dump())));
        goals.push(goal);
        choices.push(choice);
        prepared.names.push_back(id);
        prepared.dependencies.push_back(std::move(dependencies));
        earlier.emplace(id, i);
    }
    prepared.defaults.set("goals", goals);
    prepared.indices = std::move(earlier);
    // The same compiler checks all declared bindings, even currently unselected
    // ones. Operation-level numerical admission remains with its existing owner.
    if (!prepared.names.empty())
        (void)default_pipeline(workspace, prepared.defaults);
    prepared.catalog = Json::object();
    prepared.catalog.set("schema", Json::string("gadget.pipeline.intent.catalog.v1"));
    prepared.catalog.set("source_digest", Json::string(workspace.digest()));
    auto binding = Json::object();
    binding.set("source_digest", Json::string(workspace.digest()));
    binding.set("intent", prepared.declaration);
    prepared.catalog.set("intent_sha256", Json::string("sha256:" + sha256(binding.dump())));
    prepared.catalog.set("profile", bound.at("profile"));
    prepared.catalog.set("recipes", choices);
    prepared.catalog.set("observables", pipeline_detail::observable_catalog());
    prepared.catalog.set("scientific_admission",
                         pipeline_detail::answer_admission(prepared.declaration, prepared.catalog));
    if (bound.has("source_context")) {
        require(bound.has("answer_contract"),
                "pipeline source context requires an answer contract");
        const auto context = bound.at("source_context");
        keys(context, {"schema", "context", "annotations", "requirements"});
        require(context.at("schema").text() == "gadget.pipeline.source.context.v1",
                "pipeline source context schema");
        const auto annotated =
            admit_source_annotations(workspace, context.at("context"), context.at("annotations"),
                                     bound.at("answer_contract").at("question").text());
        require(annotated.at("status").text() == "prospective",
                "pipeline source mapping requires annotation mentions");
        std::set<std::string> mentions, requirements, mapped;
        const auto admitted_mentions = annotated.at("certificate").at("mentions");
        for (std::size_t i = 0; i < admitted_mentions.size(); ++i)
            mentions.insert(admitted_mentions.at(i).at("id").text());
        const auto declared_requirements = bound.at("answer_contract").at("requirements");
        for (std::size_t i = 0; i < declared_requirements.size(); ++i)
            requirements.insert(declared_requirements.at(i).at("id").text());
        const auto mappings = context.at("requirements");
        require(mappings.native().kind()==NativeRecordKind::Array && mappings.size() > 0,
                "pipeline source context requires request mappings");
        for (std::size_t i = 0; i < mappings.size(); ++i) {
            const auto mapping = mappings.at(i);
            keys(mapping, {"id", "mentions"});
            const auto id = mapping.at("id").text();
            require(requirements.contains(id) && mapped.insert(id).second,
                    "source mapping requires a unique declared answer requirement");
            const auto references = mapping.at("mentions");
            require(references.native().kind()==NativeRecordKind::Array && references.size() > 0,
                    "source mapping requires mention references");
            std::set<std::string> selected;
            for (std::size_t m = 0; m < references.size(); ++m) {
                const auto mention = references.at(m).text();
                require(mentions.contains(mention) && selected.insert(mention).second,
                        "source mapping has an unknown or duplicate mention");
            }
        }
        auto receipt = Json::object();
        receipt.set("schema", Json::string("gadget.pipeline.source.context.receipt.v1"));
        receipt.set("annotations", annotated.at("certificate"));
        receipt.set("requirements", mappings);
        receipt.set("mapping_sha256", Json::string("sha256:" + sha256(receipt.dump())));
        prepared.catalog.set("source_context", receipt);
    }
    return prepared;
}
} // namespace
Json recipe_model_catalog(const AssayWorkspace &workspace, const Json &intent) {
    return prepare_intent(workspace, intent).catalog;
}
Json admit_recipe_proposal(const AssayWorkspace &workspace, const Json &intent,
                           const Json &proposal, const std::string &question) {
    const auto prepared = prepare_intent(workspace, intent);
    keys(proposal, {"schema", "action", "source_digest", "intent_sha256"}, {"goals", "reason"});
    require(proposal.at("schema").text() == "gadget.pipeline.recipe.proposal.v1" &&
                proposal.at("source_digest").text() == workspace.digest() &&
                proposal.at("intent_sha256").text() == prepared.catalog.at("intent_sha256").text(),
            "recipe proposal source/intent/schema mismatch");
    require(!question.empty() && question.size() <= source_byte_limit() &&
                question.find('\0') == std::string::npos,
            "recipe question must follow the pipeline question contract");
    auto out = Json::object(), included = Json::array();
    const auto available = pipeline_detail::answer_admission(prepared.declaration, prepared.catalog,
                                                             nullptr, &question);
    Json scientific = available;
    if (proposal.at("action").text() == "clarify") {
        require(proposal.has("reason") && !proposal.has("goals") &&
                    !proposal.at("reason").text().empty(),
                "recipe clarification requires only a reason");
        out.set("status",
                Json::string(available.at("status").text() == "satisfied"     ? "selection_required"
                             : available.at("status").text() == "unsatisfied" ? "unsupported"
                                                                              : "clarification"));
        out.set("reason", proposal.at("reason"));
    } else {
        require(proposal.at("action").text() == "recipes" && proposal.has("goals") &&
                    !proposal.has("reason"),
                "recipe selection requires goals");
        const auto requested = proposal.at("goals");
        require(json_object_is_type(requested.raw(), json_type_array) && requested.size() > 0,
                "recipe selection must contain named goals");
        std::vector<bool> selected(prepared.names.size());
        for (std::size_t i = 0; i < requested.size(); ++i) {
            const auto found = prepared.indices.find(requested.at(i).text());
            require(found != prepared.indices.end(), "undeclared recipe binding");
            const auto index = found->second;
            require(!selected[index], "duplicate requested recipe binding");
            selected[index] = true;
        }
        // Every edge names an earlier binding. A reverse pass includes the full
        // prerequisite closure without recursion or a fixed graph-depth limit.
        for (std::size_t i = selected.size(); i > 0; --i)
            if (selected[i - 1])
                for (const auto dependency : prepared.dependencies[i - 1])
                    selected[dependency] = true;
        auto settings = Json::parse(prepared.defaults.dump()), goals = Json::array();
        settings.set("question", Json::string(question));
        for (std::size_t i = 0; i < selected.size(); ++i)
            if (selected[i]) {
                goals.push(prepared.defaults.at("goals").at(i));
                included.push(Json::string(prepared.names[i]));
            }
        settings.set("goals", goals);
        scientific = pipeline_detail::answer_admission(prepared.declaration, prepared.catalog,
                                                       &included, &question);
        if (scientific.at("status").text() == "unsatisfied")
            out.set("status",
                    Json::string(available.at("status").text() == "unsatisfied" ? "unsupported"
                                                                                : "incompatible"));
        else {
            out.set("pipeline", default_pipeline(workspace, settings));
            out.set("status", Json::string("admitted"));
        }
    }
    out.set("scientific_admission", scientific);
    out.set("available_answers", available);
    out.set("source_digest", prepared.catalog.at("source_digest"));
    out.set("intent_sha256", prepared.catalog.at("intent_sha256"));
    out.set("intent", prepared.declaration);
    out.set("question", Json::string(question));
    out.set("proposal", Json::parse(proposal.dump()));
    out.set("included_goals", included);
    if (prepared.catalog.has("source_context"))
        out.set("source_context", prepared.catalog.at("source_context"));
    return out;
}
static Json execute_pipeline_impl(std::shared_ptr<const AssayWorkspace>,const Json &,
    const fs::path &,const WorkflowModels &,std::shared_ptr<CpuResources>,
    const EvidenceStore *,const Json *,const NativeRecord *,const NativeRecord *);
NativeRecord execute_pipeline_record(std::shared_ptr<const AssayWorkspace> workspace,const NativeRecord &pipeline,
    const fs::path &destination,const WorkflowModels &models,std::shared_ptr<CpuResources> resources,
    const EvidenceStore *evidence,const NativeRecord *planning) {
    require(static_cast<bool>(workspace),"pipeline requires its actual native workspace owner");
    const auto plan=render_record_view(pipeline);std::optional<Json> planning_view;
    if(planning)planning_view=render_record_view(*planning);
    (void)execute_pipeline_impl(workspace,plan,destination,models,std::move(resources),evidence,
        planning_view?&*planning_view:nullptr,&pipeline,planning);
    return NativeRecord::read(destination/"pipeline.rfcxs");
}
static Json execute_pipeline_impl(std::shared_ptr<const AssayWorkspace> workspace, const Json &pipeline,
                      const fs::path &destination, const WorkflowModels &models,
                      std::shared_ptr<CpuResources> resources, const EvidenceStore *evidence,
                      const Json *planning,const NativeRecord *original_request,const NativeRecord *original_planning) {
    require(static_cast<bool>(workspace), "pipeline requires an owning workspace");
    const RetainedRecord native_plan=original_request?RetainedRecord(*original_request):RetainedRecord(pipeline);
    const auto admitted = static_cast<Json>(native_plan);
    validate(admitted, workspace->digest());
    std::optional<RetainedRecord> planning_record;
    if (planning) {
        planning_record = original_planning?RetainedRecord(*original_planning):RetainedRecord(Json::parse(planning->dump()));
        require(planning_record->at("status").text() == "admitted" &&
                    planning_record->at("source_digest").text() == workspace->digest() &&
                    planning_record->at("pipeline").dump() == admitted.dump() &&
                    !planning_record->has("execution"),
                "planning receipt does not bind this pipeline");
        const auto checked = admit_recipe_proposal(*workspace, planning_record->at("intent"),
                                                   planning_record->at("proposal"),
                                                   planning_record->at("question").text());
        require(checked.at("pipeline").dump() == admitted.dump() &&
                    checked.at("intent_sha256").text() ==
                        planning_record->at("intent_sha256").text() &&
                    checked.at("included_goals").dump() ==
                        planning_record->at("included_goals").dump() &&
                    checked.at("scientific_admission").dump() ==
                        planning_record->at("scientific_admission").dump() &&
                    checked.at("available_answers").dump() ==
                        planning_record->at("available_answers").dump(),
                "planning declaration mismatch");
        require(checked.has("source_context") == planning_record->has("source_context") &&
                    (!checked.has("source_context") ||
                     checked.at("source_context").dump() ==
                         planning_record->at("source_context").dump()),
                "planning source context mismatch");
    }
    require(!fs::exists(destination), "pipeline destination already exists");
    if (!resources)
        resources = CpuResources::local();
    const auto nodes = admitted.at("steps");
    const auto concurrency =
        admitted.has("concurrent_jobs")
            ? static_cast<std::size_t>(admitted.at("concurrent_jobs").integer())
            : std::min(nodes.size(), resources->coordinator().capacity());
    require(concurrency <= resources->coordinator().capacity(),
            "pipeline concurrency exceeds native CPU capacity");
    fs::create_directories(fs::absolute(destination).parent_path());
    const auto stage =
        fs::absolute(destination).parent_path() /
        (".gadget-pipeline-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    require(fs::create_directory(stage), "pipeline staging collision");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code ec;
            fs::remove_all(path, ec);
        }
    } cleanup{stage};
    std::vector<rfc::runtime::WorkJob> jobs;
    std::map<std::string, rfc::runtime::WorkKey> identities;
    std::map<std::string, std::size_t> workspace_consumers;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const auto node = nodes.at(i);
        std::set<std::string> selected{node.at("input").text()};
        if (node.at("type").text() == "bind") {
            const auto inputs = binding_inputs(node.at("request"));
            selected.insert(inputs.begin(), inputs.end());
        }
        selected.erase("source");
        for (const auto &input : selected)
            workspace_consumers[input] = rfc::machine::checked_size_add(workspace_consumers[input], 1);
    }
    const auto buffer = resources->coordinator().policy().declared_buffer_byte_limit /
                        std::min(concurrency, nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const auto node = native_plan.member("steps").element(i);
        const auto id = node.at("id").text(), input = node.at("input").text();
        const rfc::runtime::WorkKey key{i + 1, 1};
        std::vector<rfc::runtime::WorkKey> dependencies;
        if (input != "source")
            dependencies.push_back(identities.at(input));
        std::map<std::string, rfc::runtime::WorkKey> selected_inputs;
        if (node.at("type").text() == "bind")
            for (const auto &name : binding_inputs(node.at("request"))) {
                const auto selected =
                    name == "source" ? rfc::runtime::WorkKey{} : identities.at(name);
                selected_inputs.emplace(name, selected);
                if (selected.id && std::find(dependencies.begin(), dependencies.end(), selected) ==
                                       dependencies.end())
                    dependencies.push_back(selected);
            }
        if (node.has("depends_on"))
            for (std::size_t j = 0; j < node.at("depends_on").size(); ++j) {
                const auto dependency = identities.at(node.at("depends_on").at(j).text());
                if (std::find(dependencies.begin(), dependencies.end(), dependency) ==
                    dependencies.end())
                    dependencies.push_back(dependency);
            }
        const auto input_key = input == "source" ? rfc::runtime::WorkKey{} : identities.at(input);
        const auto request = node.at("request");
        const auto mathematical_request =
            node.at("type").text() == "derive" ? request.at("recipe") : request;
        const auto worker_request =
            mathematical_request.has("operation") &&
                    mathematical_request.at("operation").text() == "native_query"
                ? mathematical_request.at("query")
                : mathematical_request;
        const auto workers =
            node.has("workers")
                ? node.at("workers").integer()
                : (worker_request.has("workers") ? worker_request.at("workers").integer() : 1);
        require(workers > 0 &&
                    static_cast<std::size_t>(workers) <= resources->coordinator().capacity(),
                "pipeline step workers exceed native CPU capacity");
        // Each pipeline job delegates its assigned native workers. Inner workflows
        // use that same coordinator rather than a competing application pool.
        jobs.push_back(rfc::runtime::WorkJob::make<Artifact>(
            key, {static_cast<std::size_t>(workers), buffer, sizeof(Artifact)}, dependencies,
            [workspace, resources, stage, node, input_key, selected_inputs,
             consumers = workspace_consumers[id], &models,
             evidence](const auto &inputs, auto &assigned, std::stop_token stop) {
                if (stop.stop_requested())
                    throw rfc::runtime::ExecutionCancelled{};
                const auto parent =
                    input_key.id ? inputs.at(input_key).template as<Artifact>() : nullptr;
                std::map<rfc::runtime::WorkKey, std::shared_ptr<const AssayWorkspace>> owners;
                const auto acquire = [&](rfc::runtime::WorkKey key) {
                    if (!owners.contains(key)) {
                        const auto artifact = inputs.at(key).template as<Artifact>();
                        require(static_cast<bool>(artifact->workspace),
                                "pipeline input has no admitted workspace");
                        auto owner = artifact->workspace->acquire();
                        require(owner->digest() == artifact->receipt.at("derived_digest").text(),
                                "pipeline input digest mismatch");
                        owners.emplace(key, std::move(owner));
                    }
                    return owners.at(key);
                };
                const auto source = parent ? acquire(input_key) : workspace;
                require(static_cast<bool>(source), "pipeline input has no admitted workspace");
                if (parent)
                    require(source->digest() == parent->receipt.at("derived_digest").text(),
                            "derived pipeline input digest mismatch");
                const auto job = stage / node.at("id").text();
                fs::create_directory(job);
                auto request = Json::parse(node.at("request").dump());
                if (node.has("normalization_evidence_field")) {
                    const auto partition = request.has("partitions") && request.at("partitions").size()
                        ? request.at("partitions").at(0).integer() : 0;
                    require(partition >= 0, "negative normalization evidence partition");
                    const auto fields =
                        source
                            ->study(request.at("assay").text(), static_cast<std::size_t>(partition))
                            ->declaration()
                            .at("fields");
                    bool found = false;
                    for (std::size_t f = 0; f < fields.size(); ++f)
                        if (fields.at(f).at("id").text() ==
                            node.at("normalization_evidence_field").text()) {
                            request.at("query")
                                .at("normalization")
                                .set("evidence", fields.at(f).at("evidence"));
                            found = true;
                            break;
                        }
                    require(found, "normalization evidence field is absent from the actual input");
                }
                auto assigned_resources = resources->assigned(assigned);
                Json certificate;
                std::optional<NativeRecord> direct_result;
                std::shared_ptr<const AssayWorkspace> produced;
                if (node.at("type").text() == "bind") {
                    std::map<std::string, std::shared_ptr<const AssayWorkspace>> bindings;
                    for (const auto &[name, key] : selected_inputs) {
                        if (!key.id)
                            bindings.emplace(name, workspace);
                        else {
                            bindings.emplace(name, acquire(key));
                        }
                    }
                    const auto native_request = admit_record_view(request);
                    native_request.write(job / "request.rfcxs");
                    produced = AssayWorkspace::bind_workspaces(
                        bindings, native_request, job / "workspace");
                    certificate = produced->catalog();
                    direct_result = produced->native_publication();
                } else if (node.at("type").text() == "derive") {
                    if (AssayWorkspace::supports_native_materialization(
                            admit_record_view(request.at("recipe")))) {
                        const auto native_request = admit_record_view(request);
                        native_request.write(job / "request.rfcxs");
                        produced = source->materialize_workspace(
                            native_request, job / "workspace", assigned_resources);
                        certificate = produced->catalog();
                        direct_result = produced->study(request.at("id").text())->producer_publication();
                    } else {
                        certificate = source->materialize(request, job / "workspace", assigned_resources);
                        produced = std::make_shared<const AssayWorkspace>(job / "workspace" /
                                                                          "workspace.rfcxs");
                    }
                } else if (node.at("type").text() == "workflow") {
                    const auto partition =
                        node.has("partition") ? node.at("partition").integer() : 0;
                    require(partition >= 0, "negative pipeline partition");
                    const auto study =
                        source->study(node.at("assay").text(), static_cast<std::size_t>(partition));
                    auto bound = Json::parse(request.dump());
                    if (bound.has("source_digest"))
                        require(bound.at("source_digest").text() == study->digest(),
                                "pipeline workflow source mismatch");
                    bound.set("source_digest", Json::string(study->digest()));
                    request = bound;
                    certificate =
                        execute_workflow(study, request, models, assigned_resources, {}, evidence);
                    require(certificate.at("status").text() == "complete",
                            "pipeline workflow did not complete");
                } else if (request.at("operation").text() == "native_query") {
                    certificate = source->run(request, job / "results", assigned_resources);
                } else {
                    write_record(job / "request.rfcxs", request);
                    const auto typed_request = read_record(job / "request.rfcxs");
                    direct_result = source->execute_record(typed_request.native(), assigned_resources);
                    certificate = render_record_view(direct_result->at("result"));
                }
                if (!direct_result) write_record(job / "request.rfcxs",request);
                const auto native_request=read_record(job/"request.rfcxs");
                const auto native_result=direct_result ? std::move(*direct_result) :
                    retain_operation(source->native_publication(),native_request.native(),import_record_view(certificate));
                native_result.write(job/"result.rfcxs");
                write_text(job / "request.json", request.dump());
                write_text(job / "result.json", certificate.dump());
                auto receipt = Json::object();
                receipt.set("id", node.at("id"));
                receipt.set("type", node.at("type"));
                receipt.set("input", node.at("input"));
                receipt.set("source_digest", Json::string(source->digest()));
                receipt.set("request", Json::string(node.at("id").text() + "/request.rfcxs"));
                receipt.set("request_sha256", Json::string("sha256:" + sha256(read_text(job/"request.rfcxs"))));
                receipt.set("result", Json::string(node.at("id").text() + "/result.rfcxs"));
                receipt.set("result_sha256", Json::string("sha256:" + sha256(read_text(job/"result.rfcxs"))));
                if (produced) {
                    receipt.set("workspace",
                                Json::string(node.at("id").text() + "/workspace/workspace.rfcxs"));
                    receipt.set("derived_digest", Json::string(produced->digest()));
                }
                return Artifact{produced && consumers ?
                    std::make_shared<WorkspaceLease>(std::move(produced), consumers) : nullptr, receipt};
            }));
        identities.emplace(id, key);
    }
    rfc::runtime::WorkGraphLimits limits;
    limits.jobs = nodes.size();
    limits.in_flight = concurrency;
    limits.depth = nodes.size();
    limits.edges = 0;
    for (const auto &job : jobs)
        limits.edges = rfc::machine::checked_size_add(limits.edges, job.dependencies().size());
    limits.retained_result_bytes = rfc::machine::checked_size_mul(nodes.size(), sizeof(Artifact));
    const rfc::runtime::WorkGraph graph(std::move(jobs), limits);
    auto work = graph.submit(resources->coordinator());
    work.wait();
    auto results = Json::array();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const auto completion = work.completion(rfc::runtime::WorkKey{i + 1, 1});
        results.push(completion.result<Artifact>()->receipt);
    }
    auto index = Json::object();
    index.set("schema", Json::string("gadget.pipeline.result.v1"));
    index.set("status", Json::string("complete"));
    index.set("source_digest", Json::string(workspace->digest()));
    index.set("pipeline", admitted);
    index.set("results", results);
    index.set("concurrent_jobs", Json::integer(concurrency));
    index.set("dispatcher", Json::string("rfc::runtime::WorkGraph"));
    if (planning_record) {
        if (planning_record->at("scientific_admission").at("status").text() == "satisfied") {
            const auto answers =
                pipeline_detail::answer_witnesses(*planning_record, admitted, results, stage);
            const auto answer_bytes = answers.dump();
            write_text(stage / "answers.json", answer_bytes);
            retain_operation(workspace->native_publication(),native_plan.native(),import_record_view(answers)).write(stage/"answers.rfcxs");
            index.set("answers", Json::string("answers.rfcxs"));
            index.set("answers_sha256", Json::string("sha256:" + sha256(read_text(stage/"answers.rfcxs"))));
        }
        const auto bytes = planning_record->dump();
        write_text(stage / "planning.json", bytes);
        retain_operation(workspace->native_publication(),native_plan.native(),import_record_view(static_cast<Json>(*planning_record))).write(stage/"planning.rfcxs");
        index.set("planning", Json::string("planning.rfcxs"));
        index.set("planning_sha256", Json::string("sha256:" + sha256(read_text(stage/"planning.rfcxs"))));
    }
    index = sealed(index, "result_sha256");
    write_text(stage / "pipeline.json", index.dump());
    retain_operation(workspace->native_publication(),native_plan.native(),import_record_view(index)).write(stage/"pipeline.rfcxs");
    fs::rename(stage, destination);
    return index;
}
Json execute_pipeline(std::shared_ptr<const AssayWorkspace> workspace,const Json &pipeline,
    const fs::path &destination,const WorkflowModels &models,std::shared_ptr<CpuResources> resources,
    const EvidenceStore *evidence,const Json *planning) {
    return execute_pipeline_impl(std::move(workspace),pipeline,destination,models,std::move(resources),evidence,planning,nullptr,nullptr);
}
} // namespace gadget
