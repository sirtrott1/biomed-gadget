#include <gadget/biomedical.hpp>
#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <algorithm>
#include <iostream>

using namespace gadget;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, "invalid evidence admitted");
}
Json parse(const char *s) { return Json::parse(s); }
int main() {
    try {
        const fs::path root = GADGET_BIOMEDICAL_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root / "snapshot");
        auto manifest = Json::object(), responses = Json::array();
        manifest.set("schema", Json::string("gadget.database.snapshot.v1"));
        manifest.set("id", Json::string("authored-binding-fixture"));
        const auto add = [&](const char *id, const char *provider, const char *resource,
                             const Json &body) {
            const auto bytes = body.dump();
            write_text(root / "snapshot" / (std::string(id) + ".json"), bytes);
            write_text(root / "snapshot" / (std::string(id) + ".headers.txt"),
                       "authored offline fixture\n");
            auto r = Json::object(), q = Json::object();
            q.set("id", Json::string(id));
            q.set("provider", Json::string(provider));
            q.set("resource", Json::string(resource));
            q.set("identifier", Json::string("P00519"));
            r.set("request", q);
            r.set("url", Json::string("offline:authored-fixture"));
            r.set("path", Json::string(std::string(id) + ".json"));
            r.set("sha256", Json::string("sha256:" + sha256(bytes)));
            r.set("headers_path", Json::string(std::string(id) + ".headers.txt"));
            r.set("headers_sha256", Json::string("sha256:" + sha256("authored offline fixture\n")));
            responses.push(r);
        };
        auto activities = Json::object(), rows = Json::array();
        for (std::size_t i = 0; i < 6; ++i) {
            auto row = parse(
                R"({"molecule_chembl_id":"CHEMBL941","target_chembl_id":"CHEMBL1862","assay_chembl_id":"CHEMBLTEST","document_chembl_id":"CHEMBLDOC","standard_type":"Kd","standard_relation":"=","standard_units":"nM"})");
            row.set("activity_id", Json::integer(static_cast<std::int64_t>(100 + i)));
            row.set("standard_value", i == 2 ? Json{}
                                             : Json::string(i == 0   ? "1.25"
                                                            : i == 1 ? "2.5"
                                                            : i == 3 ? "0"
                                                            : i == 4 ? "5"
                                                                     : "1.5"));
            if (i == 1)
                row.set("standard_relation", Json::string("<"));
            if (i == 2)
                row.set("standard_relation", Json{});
            if (i == 4)
                row.set("standard_units", Json::string("uM"));
            rows.push(row);
        }
        activities.set("activities", rows);
        add("activity", "chembl", "activity", activities);
        add("binding", "bindingdb", "uniprot",
            parse(
                R"({"getLindsByUniprotsResponse":{"affinities":[{"monomerid":"10","affinity":"<1.000","affinity_type":"Ki","doi":"10.test/example"},{"monomerid":"10","affinity":"1.000","affinity_type":"Ki","doi":"10.test/example"}]}})"));
        add("target", "chembl", "target",
            parse(
                R"({"target_chembl_id":"CHEMBL1862","target_components":[{"accession":"P00519"}]})"));
        add("molecule", "chembl", "molecule",
            parse(
                R"({"molecule_chembl_id":"CHEMBL941","cross_references":[{"xref_src":"PubChem","xref_id":"100"}]})"));
        add("compound", "pubchem", "compound",
            parse(
                R"({"PropertyTable":{"Properties":[{"CID":100,"InChIKey":"AUTHORED-TEST-KEY"}]}})"));
        manifest.set("records", responses);
        manifest.set("snapshot_sha256", Json::string("sha256:" + sha256(manifest.dump())));
        write_text(root / "snapshot/snapshot.json", manifest.dump());
        const auto relative_snapshot =
            capture_database_snapshot(fs::relative(root / "snapshot", fs::current_path()));
        check(relative_snapshot.at("artifacts").size() == responses.size() * 2,
              "relative database capture lost complete response/header bytes");
        const auto retained_response = relative_snapshot.at("artifacts").at(0).at("bytes").bytes();
        const auto original_response = read_text(root / "snapshot/activity.json");
        check(std::string(reinterpret_cast<const char *>(retained_response.data()),
                          retained_response.size()) == original_response,
              "relative database capture changed original provider bytes");
        import_evidence(root / "snapshot", root / "store");
        import_evidence(root / "snapshot", root / "legacy", 1);
        auto legacy_receipt = Json::parse(read_text(root / "legacy/store.json"));
        json_object_object_del(legacy_receipt.raw(), "adapter_version");
        json_object_object_del(legacy_receipt.raw(), "source_digest");
        legacy_receipt.set("source_digest",
                           Json::string("sha256:" + sha256(legacy_receipt.dump())));
        write_text(root / "legacy/store.json", legacy_receipt.dump());
        EvidenceStore legacy(root / "legacy");
        check(legacy.execute(parse(R"({"operation":"coverage"})"))
                      .at("result")
                      .at("binding_records")
                      .integer() == 8,
              "older frozen interpretation changed on load");
        rejects([&] { import_evidence(root / "snapshot", root / "unsupported", 3); });
        EvidenceStore store(root / "store");
        auto coverage = store.execute(parse(R"({"operation":"coverage"})")).at("result");
        check(coverage.at("binding_records").integer() == 8 &&
                  coverage.at("reported").integer() == 7 && coverage.at("missing").integer() == 1,
              "missing/censoring/duplicates lost");
        check(coverage.at("equality_relation_records").integer() == 5,
              "relation categories coerced");
        const auto selected =
            store
                .execute(parse(
                    R"({"operation":"select","field":"relation","comparison":"equal","value":"<"})"))
                .at("result");
        check(selected.at("selected").integer() == 2, "native exact category selection");
        const auto bounds =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Ki","units":"nM","maximum":"0.5"})"))
                .at("result");
        check(bounds.at("selected").integer() == 0 && bounds.at("contradicted").integer() == 1 &&
                  bounds.at("indeterminate").integer() == 1,
              "censored threshold treated as a point");
        const auto exact =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Ki","units":"nM","maximum":"1"})"))
                .at("result");
        check(exact.at("selected").integer() == 2, "strict bound entailment at cutoff");
        const auto kd =
            store
                .execute(parse(
                    R"({"operation":"affinity_evidence","endpoint":"Kd","units":"nM","maximum":"1"})"))
                .at("result");
        check(kd.at("selected").integer() == 1 && kd.at("contradicted").integer() == 2 &&
                  kd.at("indeterminate").integer() == 2,
              "missing/unit/zero threshold classification");
        auto inspect = store.execute(parse(R"({"operation":"inspect"})"));
        std::string cid, sid, protein;
        const auto entities = inspect.at("result").at("entities");
        for (std::size_t i = 0; i < entities.size(); ++i) {
            const auto j = entities.at(i);
            const auto id = j.at("id").text();
            if (id == "pubchem:cid:100")
                cid = j.at("native_subject_key").text();
            if (id == "pubchem:sid:100")
                sid = j.at("native_subject_key").text();
            if (id == "uniprot:P00519")
                protein = j.at("native_subject_key").text();
        }
        check(!cid.empty() && !sid.empty() && cid != sid, "CID/SID identity collapsed");
        auto neighbors = parse(R"({"operation":"neighbors","subject_key":"1"})");
        neighbors.set("subject_key", Json::string(protein));
        const auto bindings = store.execute(neighbors).at("result").at("records");
        check(bindings.size() == 3, "native role traversal/cross-reference lost");
        auto points =
            store
                .execute(parse(
                    R"({"operation":"select","field":"provider","comparison":"equal","value":"chembl"})"))
                .at("result")
                .at("records");
        auto compare = parse(R"({"operation":"compare","left":"1","right":"2"})");
        compare.set("left", points.at(0).at("native_assertion_key"));
        compare.set("right", points.at(5).at("native_assertion_key"));
        const auto difference = store.execute(compare).at("result");
        check(difference.at("reported_difference").at("rational").text() == "-1/4" &&
                  difference.at("reported_ratio").at("rational").text() == "5/6",
              "reported exact arithmetic");
        compare.set("right", points.at(1).at("native_assertion_key"));
        rejects([&] { store.execute(compare); });
        compare.set("right", points.at(4).at("native_assertion_key"));
        rejects([&] { store.execute(compare); });
        auto study = std::make_shared<const ScientificStudy>(root / "store/study/study.json");
        ScientificSession session(study);
        const auto summary =
            session
                .execute(parse(
                    R"({"operation":"field_summary","field":"all.Kd_nM","support":"reported_only"})"))
                .at("result");
        check(summary.at("reported").integer() == 3 && summary.at("missing").integer() == 3 &&
                  summary.at("sum").at("rational").text() == "11/4",
              "censored threshold or other unit became point/zero lost");
        auto workflow = parse(
            R"({"schema":"gadget.scientific.workflow.v1","id":"evidence-to-math","source_digest":"","evidence_digest":"","limits":{"depth":8,"nodes":16,"model_calls":0,"parallel_workers":2},"definitions":[],"entry":{"id":"root","type":"parallel","workers":2,"steps":[{"id":"evidence","type":"evidence","query":{"operation":"coverage"}},{"id":"math","type":"query","query":{"operation":"support_mask","fields":["all.Kd_nM"],"mode":"all_reported","backend":"cpu_tiles","workers":2,"tile_words":1}}]}})");
        workflow.set("source_digest", Json::string(study->digest()));
        workflow.set("evidence_digest", Json::string(store.digest()));
        const auto receipt =
            execute_workflow(study, workflow, {}, CpuResources::local(), {}, &store);
        check(receipt.at("status").text() == "complete", "native delegated evidence/math workflow");
        const auto operation = execute_workflow_record(
            study, admit_record_view(workflow), {}, CpuResources::local(), {}, &store);
        operation.write(root / "complete-operation.rfcxs");
        fs::rename(root / "store", root / "unavailable-store");
        fs::rename(root / "snapshot", root / "unavailable-snapshot");
        const auto restored_operation = NativeRecord::read(root / "complete-operation.rfcxs");
        const auto sources = restored_operation.parent();
        check(sources && sources->at("evidence_digest").text() == store.digest(),
              "workflow lost its actual evidence source binding");
        const auto restored_evidence =
            NativeRecord::decode(sources->at("evidence_publication").bytes());
        check(restored_evidence.encode() == store.native_publication().encode(),
              "workflow changed or omitted its native evidence publication");
        const auto expected_raw = store.native_publication().raw_sources();
        const auto restored_raw = restored_evidence.raw_sources();
        check(expected_raw.size() == restored_raw.size(), "workflow lost original byte owners");
        for (std::size_t i = 0; i < expected_raw.size(); ++i)
            check(expected_raw[i].key() == restored_raw[i].key() &&
                      std::ranges::equal(expected_raw[i].bytes(), restored_raw[i].bytes()),
                  "workflow changed original database response/header bytes");
        fs::rename(root / "unavailable-store", root / "store");
        fs::rename(root / "unavailable-snapshot", root / "snapshot");
        workflow.set("evidence_digest", Json::string(study->digest()));
        rejects([&] { execute_workflow(study, workflow, {}, CpuResources::local(), {}, &store); });
        auto
            media = parse(R"([{"sha256":"sha256:authored-image"}])"),
            observation = parse(
                R"({"schema":"gadget.model.observation.v1","observations":[{"statement":"a white cell","image_index":0}],"limitations":"authored contract fixture"})");
        check(admit_observation(observation, study->digest(), media).at("status").text() ==
                  "prospective",
              "observation promoted or rejected");
        observation.at("observations").at(0).set("image_index", Json::integer(1));
        rejects([&] { admit_observation(observation, study->digest(), media); });
        const auto before = store.execute(parse(R"({"operation":"coverage"})")).dump();
        coverage.set("missing", Json::integer(0));
        check(before == store.execute(parse(R"({"operation":"coverage"})")).dump(),
              "public JSON mutated source");
        auto corrupt = read_text(root / "store/evidence.rfcxi");
        corrupt.back() ^= 1;
        write_text(root / "store/evidence.rfcxi", corrupt);
        rejects([&] { EvidenceStore damaged(root / "store"); });
        write_text(root / "snapshot/activity.json", "{}");
        rejects([&] { import_evidence(root / "snapshot", root / "corrupt"); });
        check(!fs::exists(root / "corrupt"), "failed import published output");
        std::cout << "PASS: retained binding evidence, native queries, censoring, missing/zero, "
                     "CID/SID, exact comparison, source binding, delegated workflow and "
                     "prospective observations\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
