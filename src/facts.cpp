#include "claim_detail.hpp"
#include "pipeline_answers.hpp"
#include "scientific_detail.hpp"
#include <gadget/facts.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace claim_detail;
using Names = std::map<std::string, Json>;
Names index_names(const Json &entries) {
    source_array(entries, "named fact source entries");
    Names out;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto entry = entries.at(i);
        require(out.emplace(text(entry.at("id"), source_byte_limit()), entry).second,
                "duplicate fact source identity");
    }
    return out;
}
Json named(const Names &entries, const std::string &id) {
    require(entries.contains(id), "undeclared fact source: " + id);
    return entries.at(id);
}
Json reading_view(const Json &catalog) {
    auto view = Json::object(), task = Json::object();
    view.set("schema", Json::string("gadget.scientific.facts.model.view.v2"));
    task.set("action", Json::string("report_values"));
    task.set("comparison", Json::string("equal"));
    task.set("scope", Json::string("all_selected_facts"));
    view.set("task", task);
    for (const auto *name : {"question", "source_digest", "question_sha256", "facts_sha256"})
        view.set(name, catalog.at(name));
    auto selected = Json::array(), owners = Json::array();
    std::map<std::string, std::string> retained_owners;
    const auto facts = catalog.at("facts");
    view.set("fact_count", cardinality(facts.size()));
    for (std::size_t i = 0; i < facts.size(); ++i) {
        const auto fact = facts.at(i), owner = fact.at("owner");
        const auto bytes = owner.dump(), digest = "sha256:" + sha256(bytes);
        const auto [entry, inserted] = retained_owners.emplace(digest, bytes);
        require(inserted || entry->second == bytes, "scientific fact owner digest collision");
        if (inserted) {
            auto retained = Json::object(), shown = Json::parse(bytes);
            retained.set("owner_sha256", Json::string(digest));
            retained.set("kind", fact.at("from").at("kind"));
            if (shown.has("reference")) {
                auto reference = shown.at("reference");
                reference.set("query_sha256",
                              Json::string("sha256:" + sha256(reference.at("query").dump())));
                json_object_object_del(reference.raw(), "query");
            }
            retained.set("owner", shown);
            owners.push(retained);
        }
        auto shown = Json::object();
        for (const auto *name :
             {"id", "value_type", "value", "meaning", "from", "path", "fact_sha256"})
            shown.set(name, fact.at(name));
        if (fact.has("coordinate")) shown.set("coordinate", fact.at("coordinate"));
        shown.set("owner_sha256", Json::string(digest));
        selected.push(shown);
    }
    view.set("facts", selected);
    view.set("owners", owners);
    return view;
}
} // namespace
struct ScientificFacts::State {
    std::optional<NativeRecord> publication;
    RetainedRecord catalog, model_view;
    std::map<std::string, RetainedRecord> facts;
    bool native_literals{};
};
ScientificFacts::ScientificFacts(const AssayWorkspace &workspace, const Json &input,
                                 const fs::path &artifact_root)
    : state_(std::make_unique<State>()) {
    const auto request = Json::parse(input.dump());
    members(request, {"schema", "source_digest", "question", "facts"}, {"context", "answers"});
    require(request.at("schema").text() == "gadget.scientific.facts.v1" &&
                request.at("source_digest").text() == workspace.digest(),
            "scientific facts schema/workspace mismatch");
    const auto question = text(request.at("question"), source_byte_limit());
    auto catalog = Json::object(), selected = Json::array(), certificates = Json::array();
    catalog.set("schema", Json::string("gadget.scientific.facts.catalog.v1"));
    catalog.set("source_digest", Json::string(workspace.digest()));
    catalog.set("question", Json::string(question));
    catalog.set("question_sha256", Json::string("sha256:" + sha256(question)));
    catalog.set("request", request);
    Json context, answers;
    Names targets, answer_names;
    using ArtifactKey = std::pair<std::string, std::string>;
    std::map<std::string, std::map<ArtifactKey, Json>> answer_artifacts;
    fs::path answer_root;
    using N=NativeRecordValue;
    std::vector<RawByteSource> retained_bytes;
    N::Array originals;
    std::map<fs::path, std::size_t> retained_paths;
    std::uint64_t next_source=1;
    for(const auto &raw:workspace.native_publication().raw_sources()) {
        require(raw.key()!=std::numeric_limits<std::uint64_t>::max(),"fact byte identity space exhausted");
        next_source = std::max(next_source, raw.key() + 1);
    }
    const auto retain_artifact = [&](const fs::path &root, const Json &record) {
        const auto path = pipeline_detail::checked_answer_path(root, record);
        if (!retained_paths.contains(path)) {
            require(next_source != std::numeric_limits<std::uint64_t>::max(),
                    "fact byte identity space exhausted");
            const auto key = next_source++;
            const auto source = RawByteSource::read(key, path);
            retained_paths.emplace(path, retained_bytes.size());
            retained_bytes.push_back(source);
            originals.elements.emplace_back(
                N::Object{{{"reference", import_record_view(record)},
                           {"bytes", N(N::ByteSpan{key, 0, source.bytes().size(), false})}}});
        }
        const auto source = retained_bytes.at(retained_paths.at(path));
        const auto bytes = source.bytes();
        require(record.at("sha256").text() ==
                    "sha256:" + sha256(std::string_view(
                                    reinterpret_cast<const char *>(bytes.data()), bytes.size())),
                "fact artifact byte hash mismatch");
        return source;
    };
    if (request.has("context")) {
        context = workspace.source_context(request.at("context"));
        require(context.at("question").text() == question, "fact source question mismatch");
        catalog.set("context", context);
        targets = index_names(context.at("targets"));
    }
    if (request.has("answers")) {
        require(!artifact_root.empty(), "answer facts require an artifact root");
        const auto record = request.at("answers");
        members(record, {"path", "sha256"});
        retain_artifact(artifact_root,record);
        answers = pipeline_detail::checked_answer_view(artifact_root,record);
        verify_seal(answers, "answers_sha256");
        require(answers.at("schema").text() == "gadget.pipeline.answers.v1" &&
                    answers.at("source_digest").text() == workspace.digest() &&
                    answers.at("question_sha256").text() == catalog.at("question_sha256").text(),
                "completed answer owner/question mismatch");
        if (context.raw() && answers.has("source_context"))
            require(answers.at("source_context").at("annotations").at("context").dump() ==
                        context.dump(),
                    "completed answer source context mismatch");
        catalog.set("answers", answers);
        answer_names = index_names(answers.at("answers"));
        for (const auto &[id, answer] : answer_names) {
            const auto artifacts = answer.at("artifacts");
            source_array(artifacts, "completed answer artifacts");
            for (std::size_t a = 0; a < artifacts.size(); ++a) {
                const auto artifact = artifacts.at(a);
                const ArtifactKey key{text(artifact.at("binding"), source_byte_limit()),
                                      text(artifact.at("step"), source_byte_limit())};
                require(answer_artifacts[id].emplace(key, artifact).second,
                        "duplicate completed answer binding/step");
            }
        }
        answer_root =
            fs::canonical(artifact_root / fs::path(record.at("path").text())).parent_path();
    }
    const auto selections = request.at("facts");
    source_array(selections, "scientific fact selections");
    require(selections.size(), "scientific facts require selected values");
    // Each referenced certificate is captured once; the admitted object owns it
    // independently of later filesystem edits or model execution.
    using ResultKey = std::pair<std::string, std::string>;
    std::map<ResultKey, Json> loaded;
    std::map<ResultKey, std::string> loaded_hashes;
    std::set<std::string> selected_ids;
    for (std::size_t i = 0; i < selections.size(); ++i) {
        const auto selection = selections.at(i), from = selection.at("from");
        members(selection, {"id", "meaning", "from", "path", "value_type"});
        const auto id = text(selection.at("id"), source_byte_limit());
        require(selected_ids.insert(id).second, "duplicate selected fact identity");
        (void)text(selection.at("meaning"), source_byte_limit());
        const auto kind = from.at("kind").text();
        Json value, owner;
        if (kind == "source") {
            members(from, {"kind", "target"});
            require(context.raw(), "source facts require a declared source context");
            const auto target = named(targets, text(from.at("target"), source_byte_limit()));
            value = member_path(target.at("resolved").at("descriptor"), selection.at("path"));
            owner = Json::object();
            owner.set("context_sha256", context.at("context_sha256"));
            owner.set("binding_sha256", target.at("binding_sha256"));
            owner.set("target", target.at("target"));
        } else {
            require(kind == "answer", "unknown scientific fact source kind");
            members(from, {"kind", "answer", "binding", "step", "result"});
            require(answers.raw(), "result facts require completed answers");
            const auto id = text(from.at("answer"), source_byte_limit());
            const auto answer = named(answer_names, id);
            const ArtifactKey key{text(from.at("binding"), source_byte_limit()),
                                  text(from.at("step"), source_byte_limit())};
            require(answer_artifacts.contains(id) && answer_artifacts.at(id).contains(key),
                    "undeclared completed answer binding/step");
            const auto artifact = answer_artifacts.at(id).at(key);
            const auto results = artifact.at("results");
            source_array(results, "completed answer results");
            require(results.size(), "completed answer has no result");
            const auto reference = results.at(
                count(from.at("result"), results.size() - 1, "answer result coordinate"));
            require(reference.at("pointer").text() == "/result",
                    "answer fact requires complete result pointer");
            const auto path = reference.at("path").text();
            const ResultKey result_key{path, reference.has("native_selection")
                                                 ? reference.at("native_selection").dump()
                                                 : std::string{}};
            if (!loaded.contains(result_key)) {
                const auto retained_source = retain_artifact(answer_root, reference);
                const auto certificate =
                    reference.has("native_selection")
                        ? pipeline_detail::checked_answer_view(answer_root, reference)
                        : [&] {
                              const auto value = decode_record(retained_source.bytes());
                              if (value.has("schema") &&
                                  value.at("schema").text() == "gadget.operation.native.v1")
                                  return static_cast<Json>(value.at("result"));
                              return static_cast<Json>(value);
                          }();
                verify_seal(certificate, "certificate_sha256");
                const auto schema = certificate.at("schema").text();
                require(schema == "gadget.assay.result.v1" ||
                            schema == "gadget.scientific.certificate.v1" ||
                            schema == "gadget.biomedical.certificate.v1",
                        "unsupported answer certificate");
                loaded.emplace(result_key, certificate);
                loaded_hashes.emplace(result_key, reference.at("sha256").text());
                auto retained = Json::object();
                retained.set("reference", reference);
                retained.set("certificate", certificate);
                certificates.push(retained);
            }
            const auto certificate = loaded.at(result_key);
            require(reference.at("sha256").text() == loaded_hashes.at(result_key) &&
                        certificate.at("source_digest").dump() ==
                            reference.at("source_digest").dump() &&
                        certificate.at("query").dump() == reference.at("query").dump() &&
                        certificate.at("certificate_sha256").dump() ==
                            reference.at("certificate_sha256").dump(),
                    "answer result byte/source/query/certificate mismatch");
            value = member_path(certificate.at("result"), selection.at("path"));
            owner = Json::object();
            owner.set("answers_sha256", answers.at("answers_sha256"));
            owner.set("observable", answer.at("observable"));
            owner.set("reference", reference);
        }
        auto fact = Json::parse(selection.dump());
        fact.set("observed", value);
        fact.set("value", literal_value(value, selection.at("value_type").text()));
        fact.set("owner", owner);
        fact = seal(fact, "fact_sha256");
        selected.push(fact);
    }
    catalog.set("certificates", certificates);
    catalog.set("facts", selected);
    state_->publication=workspace.native_publication().derive(N(N::Object{{
        {"schema",N(std::string("gadget.scientific.facts.native.v1"))},
        {"catalog",import_record_view(seal(catalog,"facts_sha256"))},
        {"originals",N(std::move(originals))}}}),std::move(retained_bytes));
    state_->catalog=state_->publication->at("catalog");
    state_->facts.clear();
    for(std::size_t i=0;i<selected.size();++i)state_->facts.emplace(selected.at(i).at("id").text(),state_->catalog.member("facts").element(i));
    state_->model_view = reading_view(state_->catalog);
}
ScientificFacts::ScientificFacts(const AssayWorkspace &workspace,const NativeRecord &request,
    const fs::path &artifact_root)
    :ScientificFacts(workspace,render_record_view(request),artifact_root) {
    const auto result=state_->publication->materialize_value();
    state_->publication=retain_operation(*state_->publication,request,result);
    state_->catalog=state_->publication->at("result").at("catalog");
    state_->facts.clear();
    const auto selected=state_->catalog.member("facts");
    for(std::size_t i=0;i<selected.size();++i)
        state_->facts.emplace(selected.element(i).at("id").text(),selected.element(i));
}
ScientificFacts::ScientificFacts(const NativeRecord &operation, const NativeRecord &request)
    : state_(std::make_unique<State>()) {
    using N = NativeRecordValue;
    require(operation.is_root() && operation.at("schema").text() == "gadget.operation.native.v1",
            "native facts require an owning native operation");
    const auto input = render_record_view(request);
    members(input, {"schema", "source_digest", "source_namespace", "cutoff", "question", "facts"});
    require(input.at("schema").text() == "gadget.scientific.facts.request.v2" &&
                input.at("source_namespace").text() == std::to_string(operation.source_namespace()),
            "native facts request schema/namespace mismatch");
    const auto history = operation.publication_history();
    require(history.has_value() && request.at("cutoff").kind() == NativeRecordKind::Integer &&
                request.at("cutoff").integer() == rfc::numeric::Integer(std::to_string(history->last_sequence())),
            "native facts request cutoff mismatch");
    const auto result = operation.at("result");
    const auto schema = result.at("schema").text();
    require(schema == "gadget.scientific.result.fields.v1" || schema == "gadget.workflow.native.receipt.v1",
            "native facts refuse an unmigrated operation");
    const bool workflow = schema == "gadget.workflow.native.receipt.v1";
    if (workflow)
        require(result.at("native_execution").boolean() && result.at("status").text() == "complete",
                "native facts require a completed native workflow");
    const auto digest = (workflow ? result : result.at("certificate")).at("source_digest").text();
    require(input.at("source_digest").text() == digest, "native facts source digest mismatch");
    const auto question = text(input.at("question"), source_byte_limit());
    const auto selections = input.at("facts");
    source_array(selections, "native fact selections");
    require(selections.size(), "native facts require selected values");
    std::map<std::string, NativeRecord> owners;
    std::set<std::string> identities;
    N::Array facts;
    for (std::size_t i = 0; i < selections.size(); ++i) {
        const auto selection = selections.at(i);
        members(selection, {"id", "meaning", "from", "path", "value_type"}, {"coordinate"});
        const auto id = text(selection.at("id"), source_byte_limit());
        require(identities.insert(id).second, "duplicate native fact identity");
        (void)text(selection.at("meaning"), source_byte_limit());
        const auto from = selection.at("from");
        members(from, {"kind", "path"});
        require(from.at("kind").text() == "native", "native fact source kind");
        const auto path = from.at("path");
        source_array(path, "native operation path");
        const auto key = path.dump();
        if (!owners.contains(key)) {
            auto selected = operation;
            for (std::size_t j = 0; j < path.size(); ++j) {
                const auto part = path.at(j);
                if (json_object_is_type(part.raw(), json_type_string))
                    selected = selected.at(text(part, source_byte_limit()));
                else
                    selected = selected.at(count(part, selected.size() ? selected.size()-1 : 0,
                                                  "native operation path coordinate"));
            }
            if (selected.has("schema") && selected.at("schema").text() == "gadget.workflow.publication.v1") {
                const auto owner = NativeRecord::decode(selected.at("publication").bytes());
                require(selected.at("source_namespace").text() == std::to_string(owner.source_namespace()) &&
                            selected.at("cutoff").integer() == rfc::numeric::Integer(std::to_string(owner.publication_history()->last_sequence())),
                        "native workflow reference owner mismatch");
                selected = owner;
            }
            require(selected.is_root() && selected.at("schema").text() == "gadget.operation.native.v1" &&
                        selected.at("result").at("schema").text() == "gadget.scientific.result.fields.v1" &&
                        selected.at("result").at("fields").at("schema").text() == "gadget.scientific.fields.v1" &&
                        selected.at("result").at("certificate").at("source_digest").text() == digest,
                    "native fact requires the completed source-bound scientific owner");
            owners.emplace(key, std::move(selected));
        }
        const auto &owner = owners.at(key);
        const auto view = owner.at("result").at("fields").at("view");
        const auto selected = select_scientific_result(view, request.at("facts").at(i).at("path"));
        std::optional<std::size_t> coordinate;
        if (selection.has("coordinate"))
            coordinate = count(selection.at("coordinate"), std::numeric_limits<std::size_t>::max(),
                               "native fact original flat coordinate");
        auto observed = read_native_literal(selected, selection.at("value_type").text(), coordinate);
        auto fact = import_record_view(selection);
        auto &members = std::get<N::Object>(fact.data).members;
        members.emplace_back("observed", observed.value);
        members.emplace_back("value", std::move(observed.value));
        members.emplace_back("owner", N(N::Object{{
            {"enclosing_namespace", N(std::to_string(operation.source_namespace()))},
            {"enclosing_cutoff", N(rfc::numeric::Integer(std::to_string(history->last_sequence())))},
            {"operation_namespace", N(std::to_string(owner.source_namespace()))},
            {"operation_cutoff", N(rfc::numeric::Integer(std::to_string(owner.publication_history()->last_sequence())))},
            {"source_digest", N(digest)}, {"binding", std::move(observed.binding)}}}));
        const auto publication = NativeRecord::create(fact, new_native_record_namespace());
        const auto shown = render_record_view(publication);
        members.emplace_back("fact_sha256", N(std::string("sha256:" + sha256(shown.dump()))));
        facts.elements.push_back(std::move(fact));
    }
    auto catalog = import_record_view(input);
    auto &members = std::get<N::Object>(catalog.data).members;
    for (auto &[name, value] : members) {
        if (name == "schema") value = N(std::string("gadget.scientific.facts.catalog.v1"));
        if (name == "facts") value = N(std::move(facts));
    }
    members.emplace_back("request", import_record_view(input));
    members.emplace_back("certificates", N(N::Array{}));
    members.emplace_back("question_sha256", N(std::string("sha256:" + sha256(question))));
    const auto shown = render_record_view(NativeRecord::create(catalog, new_native_record_namespace()));
    members.emplace_back("facts_sha256", N(std::string("sha256:" + sha256(shown.dump()))));
    auto source = NativeRecord::annotate(N(N::Object{{
        {"schema", N(std::string("gadget.scientific.facts.native.v2"))},
        {"source_operation", operation.materialize_value()}, {"catalog", std::move(catalog)}}}), operation);
    state_->publication = retain_operation(source, request, source.unbound_value());
    state_->catalog = state_->publication->at("result").at("catalog");
    const auto retained = state_->catalog.member("facts");
    for (std::size_t i = 0; i < retained.size(); ++i)
        state_->facts.emplace(retained.element(i).at("id").text(), retained.element(i));
    state_->model_view = reading_view(state_->catalog);
    state_->native_literals = true;
}
ScientificFacts::~ScientificFacts() = default;
const NativeRecord &ScientificFacts::native_publication() const {return state_->publication.value();}
Json ScientificFacts::catalog() const { return Json::parse(state_->catalog.dump()); }
Json ScientificFacts::model_view() const { return Json::parse(state_->model_view.dump()); }
Json ScientificFacts::check(const Json &input) const {
    const auto proposal = Json::parse(input.dump());
    members(proposal, {"schema", "action", "source_digest", "facts_sha256", "question_sha256"},
            {"claims", "reason"});
    require(proposal.at("schema").text() == "gadget.scientific.claims.proposal.v1",
            "scientific claims schema");
    for (const auto *name : {"source_digest", "facts_sha256", "question_sha256"})
        require(proposal.at(name).dump() == state_->catalog.at(name).dump(),
                "scientific claim catalog binding mismatch");
    auto checks = Json::array(), missing = Json::array();
    std::set<std::string> ids, supplied;
    std::size_t matched = 0;
    const auto action = proposal.at("action").text();
    if (action == "clarify") {
        require(proposal.has("reason") && !proposal.has("claims"),
                "scientific clarification requires only reason");
        (void)text(proposal.at("reason"), source_byte_limit());
    } else {
        require(action == "claims" && proposal.has("claims") && !proposal.has("reason"),
                "scientific claims require comparisons");
        const auto claims = proposal.at("claims");
        source_array(claims, "scientific literal claims");
        require(claims.size(), "scientific claims require at least one comparison");
        for (std::size_t i = 0; i < claims.size(); ++i) {
            const auto claim = claims.at(i);
            members(claim, {"id", "fact", "comparison", "value"});
            require(ids.insert(text(claim.at("id"), source_byte_limit())).second,
                    "duplicate literal claim identity");
            const auto id = text(claim.at("fact"), source_byte_limit());
            require(state_->facts.contains(id), "literal claim names an undeclared fact");
            supplied.insert(id);
            const auto fact = state_->facts.at(id);
            if (fact.at("value_type").text() == "rational")
                require(json_object_is_type(claim.at("value").raw(), json_type_string),
                        "rational claim values must be literal strings");
            auto literal = Json::object(), path = Json::array();
            path.push(Json::string("value"));
            literal.set("path", path);
            literal.set("value_type", fact.at("value_type"));
            literal.set("comparison", claim.at("comparison"));
            literal.set("value", claim.at("value"));
            auto check = Json::object();
            auto comparison = Json::object();
            if (state_->native_literals) {
                comparison.set("claim", literal);
                comparison.set("observed", fact.at("value"));
                comparison.set("binding", fact.at("owner").at("binding"));
                comparison.set("agrees", Json::boolean(compare_native_literal(fact.native().at("observed"),
                    claim.at("value"), fact.at("value_type").text(), claim.at("comparison").text())));
                comparison = seal(comparison, "check_sha256");
            } else comparison = check_literal_claim(fact, literal);
            if (comparison.at("agrees").boolean())
                ++matched;
            check.set("claim", claim);
            check.set("fact_sha256", fact.at("fact_sha256"));
            check.set("comparison", comparison);
            checks.push(seal(check, "check_sha256"));
        }
    }
    for (std::size_t i = 0; i < state_->catalog.at("facts").size(); ++i) {
        const auto id = state_->catalog.at("facts").at(i).at("id").text();
        if (!supplied.contains(id))
            missing.push(Json::string(id));
    }
    auto receipt = Json::object();
    receipt.set("schema", Json::string("gadget.scientific.claims.receipt.v1"));
    receipt.set("status", Json::string(action == "clarify"        ? "selection_required"
                                       : checks.size() != matched ? "contradicted"
                                       : missing.size()           ? "incomplete"
                                                                  : "matched"));
    receipt.set("source_digest", state_->catalog.at("source_digest"));
    receipt.set("facts_sha256", state_->catalog.at("facts_sha256"));
    receipt.set("question_sha256", state_->catalog.at("question_sha256"));
    receipt.set("catalog", state_->catalog);
    receipt.set("proposal", proposal);
    receipt.set("checks", checks);
    receipt.set("missing", missing);
    receipt.set("matched", cardinality(matched));
    receipt.set("contradicted", cardinality(checks.size() - matched));
    receipt.set("interpretation_verified", Json::boolean(false));
    receipt = seal(receipt, "receipt_sha256");
    auto out = Json::object();
    out.set("schema", Json::string("gadget.agent.response.v1"));
    out.set("status", receipt.at("status"));
    out.set("proposal", proposal);
    out.set("certificate", receipt);
    return out;
}
} // namespace gadget
