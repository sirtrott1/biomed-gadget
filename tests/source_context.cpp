#include "../examples/scientific/fixture.hpp"
#include <gadget/pipeline.hpp>
#include <iostream>

namespace {
using namespace gadget;
std::size_t checks{};
void check(bool condition, const std::string &message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&action) {
    bool failed = false;
    try {
        action();
    } catch (const std::invalid_argument &) {
        failed = true;
    }
    check(failed, "invalid source binding was accepted");
}
Json read(const fs::path &path) { return Json::parse(read_text(path)); }
void native(const fs::path &root) {
    fs::remove_all(root);
    example::scientific_study(root);
    const ScientificStudy study(root / "study.json");
    const auto declaration = study.declaration();
    const auto before = declaration.dump();
    auto target = Json::object();
    target.set("source_digest", Json::string(study.digest()));
    auto evidence = declaration.at("cells").at(0).at("evidence");
    for (const auto *kind : {"cell", "field", "context"}) {
        auto selector = Json::parse(target.dump());
        selector.set("kind", Json::string(kind));
        selector.set("id", Json::string(std::string(kind) == "cell"    ? "gene:A"
                                        : std::string(kind) == "field" ? "baseline.E"
                                                                       : "baseline"));
        auto result = study.source_target(selector, evidence);
        check(result.at("source_namespace").dump() ==
                  declaration.at("native").at("source_namespace").dump(),
              "native owner mismatch");
        check(result.at("evidence").at(0).has("document_sha256"), "quote hash absent");
        result.at("descriptor").set("id", Json::string("mutated consumer JSON"));
        check(study.source_target(selector, evidence).at("descriptor").at("id").text() !=
                  "mutated consumer JSON",
              "mutable source descriptor");
        selector.set("source_digest", Json::string("sha256:" + std::string(64, '0')));
        rejects([&] { study.source_target(selector, evidence); });
    }
    std::size_t zero{}, inactive{};
    const auto occurrences = declaration.at("occurrences");
    for (std::size_t i = 0; i < occurrences.size(); ++i) {
        const auto occurrence = occurrences.at(i);
        auto selector = Json::parse(target.dump());
        selector.set("kind", Json::string("occurrence"));
        selector.set("grade", occurrence.at("grade"));
        selector.set("key", occurrence.at("key"));
        const auto result = study.source_target(selector, occurrence.at("evidence"));
        const auto descriptor = result.at("descriptor");
        if (descriptor.at("coefficient").text() == "0")
            ++zero;
        if (descriptor.at("aspects").has("exists") &&
            !descriptor.at("aspects").at("exists").boolean())
            ++inactive;
        check(descriptor.has("upper") && descriptor.has("lower"), "lost incidence endpoints");
        check(descriptor.at("key").text() == occurrence.at("key").text(),
              "lost occurrence identity");
    }
    check(zero > 0, "fixture zero occurrence not retained");
    check(inactive > 0, "fixture inactive occurrence not retained");
    target.set("kind", Json::string("cell"));
    target.set("id", Json::string("gene:A"));
    auto invalid = Json::parse(evidence.dump());
    invalid.at(0).set("text", Json::string("invented evidence"));
    rejects([&] { study.source_target(target, invalid); });
    const auto document = read_text(root / "design.txt");
    const auto position = document.find("β");
    invalid = Json::array();
    auto split = Json::object();
    split.set("document", Json::string("design"));
    split.set("start", Json::integer(position));
    split.set("end", Json::integer(position + 1));
    split.set("text", Json::string(document.substr(position, 1)));
    invalid.push(split);
    rejects([&] { study.source_target(target, invalid); });
    check(study.declaration().dump() == before, "annotation changed native declaration");
    // The earlier cell-only interface also accepts source-sized mention counts.
    auto proposal = Json::object(), mentions = Json::array();
    proposal.set("schema", Json::string("gadget.semantic.proposal.v1"));
    proposal.set("source_digest", Json::string(study.digest()));
    for (std::size_t i = 0; i < 257; ++i) {
        auto mention = Json::object();
        mention.set("entity", Json::string("gene:A"));
        mention.set("kind", Json::string("gene"));
        mention.set("polarity", Json::string("uncertain"));
        mention.set("evidence", evidence);
        mentions.push(mention);
    }
    proposal.set("mentions", mentions);
    ScientificSession session(std::make_shared<const ScientificStudy>(root / "study.json"));
    check(session.propose_annotations(proposal).at("mentions").size() == 257,
          "legacy annotation ceiling retained");
}
void pipeline(const fs::path &manifest, const fs::path &intent_path, const fs::path &proposal_path,
              const fs::path &root) {
    auto workspace = std::make_shared<const AssayWorkspace>(manifest);
    const auto intent = read(intent_path), proposal = read(proposal_path);
    const auto question = intent.at("answer_contract").at("question").text();
    const auto planning = admit_recipe_proposal(*workspace, intent, proposal, question);
    check(planning.at("status").text() == "admitted", "test plan was not admitted");
    const auto before = workspace->catalog().dump();
    auto view = workspace->source_context(intent.at("source_context").at("context"));
    view.at("targets").at(0).at("resolved").at("descriptor").set("id", Json::string("mutated"));
    check(workspace->catalog().dump() == before, "context view aliases workspace");
    for (const auto *member : {"mapping_sha256", "requirements", "annotations"}) {
        auto forged = Json::parse(planning.dump());
        forged.at("source_context").set(member, Json::string("forged"));
        rejects([&] {
            execute_pipeline(workspace, planning.at("pipeline"), root, {}, {}, nullptr, &forged);
        });
        check(!fs::exists(root), "invalid context published execution");
    }
    auto omitted = Json::parse(planning.dump());
    json_object_object_del(omitted.raw(), "source_context");
    rejects([&] {
        execute_pipeline(workspace, planning.at("pipeline"), root, {}, {}, nullptr, &omitted);
    });
    check(!fs::exists(root), "omitted context published execution");
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "native")
            native(argv[2]);
        else if (argc == 6 && std::string(argv[1]) == "pipeline")
            pipeline(argv[2], argv[3], argv[4], argv[5]);
        else
            throw std::invalid_argument("source context test arguments");
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
