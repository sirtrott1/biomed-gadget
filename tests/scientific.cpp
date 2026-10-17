#include "../examples/scientific/fixture.hpp"
#include <algorithm>
#include <future>
#include <gadget/model.hpp>
#include <gadget/resources.hpp>
#include <gadget/structure.hpp>
#include <gadget/workflow.hpp>
#include <iostream>

namespace {
using namespace gadget;
void check(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&action, const char *message) {
    bool failed = false;
    try {
        action();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, message);
}
Json query(const std::string &source) { return Json::parse(source); }
} // namespace
int main() {
    try {
        const fs::path root = GADGET_SCIENTIFIC_TEST_BUILD;
        fs::remove_all(root);
        example::scientific_study(root);
        auto study = std::make_shared<const ScientificStudy>(root / "study.json");
        ScientificSession session(study);
        const auto old_moment = query(
            R"({"operation":"typed_moment","context":"baseline","field":"baseline.G","other_field":"baseline.E","order":1})");
        const auto new_moment = query(
            R"({"operation":"typed_moment","context":"followup","field":"followup.G","other_field":"followup.E","order":1})");
        const auto temporal = query(
            R"({"operation":"temporal_moment_change","context":"baseline","other_context":"followup","field":"baseline.G","other_field":"baseline.E","new_field":"followup.G","new_other_field":"followup.E","order":1})");
        const auto old = session.execute(old_moment);
        const auto next = session.execute(new_moment);
        const auto change = session.execute(temporal);
        check(old.at("result").at("rational").text() == "-1", "first declared moment");
        check(next.at("result").at("rational").text() == "0", "second declared moment");
        auto result = change.at("result");
        check(result.at("pullback_form").at("rational").text() == "-2/3" &&
                  result.at("left_innovation").at("rational").text() == "5/3" &&
                  result.at("right_innovation").at("rational").text() == "1/3" &&
                  result.at("innovation_interaction").at("rational").text() == "-1/3" &&
                  result.at("total").at("rational").text() == "1" &&
                  result.at("finite_rate").at("rational").text() == "1/2",
              "complete exact longitudinal identity");
        auto zero = Json::parse(temporal.dump());
        zero.set("order", Json::integer(0));
        check(session.execute(zero).at("result").at("total").at("rational").text() == "6",
              "zeroth longitudinal moment");
        const auto structure = query(
            R"({"operation":"structure_transition","context":"baseline","other_context":"followup","field":"baseline.G","new_field":"followup.G","lambda":"1"})");
        const auto structural = session.execute(structure);
        check(structural.at("result").at("identity_exact").boolean() &&
                  structural.at("result").at("sector_transport").size() == 3,
              "native structural transport and Green change");
        const auto structural_preparations = session.temporal_preparation_statistics().dump();
        check(session.temporal_preparation_statistics().at("green_preparations").integer() == 2 &&
                  session.execute(structure).dump() == structural.dump() &&
                  session.temporal_preparation_statistics().dump() == structural_preparations,
              "two-context transition did not reuse original endpoint preparations");
        ScientificSession validation_session(study);
        validation_session.validate(structure);
        check(validation_session.temporal_preparation_statistics()
                      .at("green_preparations")
                      .integer() == 0,
              "validation prepared numerical Green operations");
        auto factored_structure = Json::parse(structure.dump());
        factored_structure.set("map_output", Json::string("factored"));
        ScientificSession factored_session(study);
        factored_session.validate(factored_structure);
        check(
            factored_session.temporal_preparation_statistics().at("green_preparations").integer() ==
                0,
            "factored validation prepared numerical actions");
        const auto factored_result = factored_session.execute(factored_structure);
        const auto factored_counts = factored_session.temporal_preparation_statistics().dump();
        check(
            factored_result.at("result").at("factored_maps").at("format").text() ==
                    "gadget-factored-maps-v1" &&
                factored_result.at("result").at("sector_maps").size() == 3 &&
                factored_session.temporal_preparation_statistics().at("green_matrices").integer() ==
                    0 &&
                factored_result.at("result").at("total").dump() ==
                    structural.at("result").at("total").dump() &&
                factored_session.execute(factored_structure).dump() == factored_result.dump() &&
                factored_session.temporal_preparation_statistics().dump() == factored_counts,
            "factored transition lost complete fields or prepared full Green matrices");
        auto bad_structure = Json::parse(structure.dump());
        bad_structure.set("lambda", Json::string("0"));
        rejects([&] { session.execute(bad_structure); },
                "zero structural Green parameter admitted");
        const auto controlled = query(
            R"({"operation":"controlled_pairings","fields":["baseline.G"],"endpoint":"baseline.T","controls":["baseline.E"],"selection":[],"support":"joint_reported","metric":"identity","intercept":true})");
        const auto conditional = session.execute(controlled);
        check(conditional.at("result").at("n").integer() == 3,
              "endpoint query must work with the generic native fixture");
        auto changed_conditional = Json::parse(conditional.dump());
        changed_conditional.at("field_semantics")
            .at("baseline.G")
            .set("units", Json::string("edited"));
        check(session.execute(controlled).dump() == conditional.dump(),
              "endpoint certificate exposed mutable admitted metadata");
        const auto coverage =
            session.execute(query(R"({"operation":"coverage","field":"baseline.partial"})"));
        check(coverage.at("result").at("missing").integer() == 1, "missing coverage");
        const auto rows = coverage.at("result").at("observations");
        check(rows.at(0).at("entity").text() == "gene:A" &&
                  rows.at(1).at("entity").text() == "gene:B" &&
                  rows.at(2).at("entity").text() == "gene:C",
              "coverage coordinates must retain native cell identities");
        check(session.execute(query(R"({"operation":"coverage","field":"baseline.face"})"))
                      .at("result")
                      .at("observations")
                      .at(0)
                      .at("entity")
                      .text() == "assertion:R201",
              "higher-grade coverage must retain native identity");
        check(rows.at(0).at("reported").boolean() &&
                  rows.at(0).at("value").at("rational").text() == "0" &&
                  !rows.at(2).at("reported").boolean() &&
                  json_object_is_type(rows.at(2).at("value").raw(), json_type_null),
              "reported zero differs from absence");
        const auto summary = session.execute(query(
            R"({"operation":"field_summary","field":"baseline.partial","support":"reported_only"})"));
        check(summary.at("result").at("mean").at("rational").text() == "1" &&
                  summary.at("result").at("sum_of_squares").at("rational").text() == "4" &&
                  summary.at("result").at("reported").integer() == 2,
              "reported-only exact summary changed missing support or zero");
        rejects(
            [&] {
                session.execute(query(
                    R"({"operation":"field_summary","field":"baseline.batch","support":"reported_only"})"));
            },
            "categorical summary admitted");
        rejects(
            [&] {
                session.execute(query(
                    R"({"operation":"hodge","context":"baseline","field":"baseline.partial"})"));
            },
            "missing observations silently zero-extended");
        rejects(
            [&] {
                session.execute(query(
                    R"({"operation":"hodge","context":"baseline","field":"baseline.batch"})"));
            },
            "categorical covariate treated as linear amplitude");
        const auto boundary = session.execute(
            query(R"({"operation":"boundary","context":"baseline","field":"baseline.face"})"));
        check(boundary.at("output_grade").integer() == 1 && boundary.at("result").size() == 2 &&
                  boundary.at("result").at(1).at("entity").text() == "assertion:R103" &&
                  boundary.at("result").at(1).at("value").at("rational").text() == "4",
              "declared grade-two boundary lost cancelling lower assertion");
        const auto upper_hodge = session.execute(
            query(R"({"operation":"hodge","context":"baseline","field":"baseline.face"})"));
        check(upper_hodge.at("result").at(0).at("value").at("rational").text() == "20",
              "non-unit metric higher Hodge");
        const auto evidence =
            session.execute(query(R"({"operation":"provenance","entity":"gene:C"})"));
        check(evidence.at("result").at("occurrences").size() == 2,
              "inactive occurrence missing from provenance");
        const auto cached = query(
            R"({"operation":"native_boundary","field_set":"synthetic.assays","backend":"cpu_cached"})");
        const auto first = session.execute(cached);
        const auto second = session.execute(cached);
        check(!first.at("execution").at("result_reused").boolean() &&
                  second.at("execution").at("result_reused").boolean(),
              "canonical P3 result cache not adopted");
        check(first.at("execution").at("retained_occurrences").integer() == 9 &&
                  first.at("execution").at("active_occurrences").integer() == 8,
              "retained/inactive schedule changed");
        check(first.at("result").at(0).at("values").at(0).at("rational").text() == "-2" &&
                  first.at("result").at(1).at("values").at(0).at("rational").text() == "-3" &&
                  first.at("result").at(2).at("values").at(1).at("rational").text() == "7",
              "native two-coordinate boundary");
        auto tiled = query(
            R"({"operation":"native_boundary","field_set":"synthetic.assays","backend":"cpu_tiles","workers":4,"tile_rows":1})");
        const auto tiled_workers =
            std::min<std::size_t>(4, CpuResources::local()->coordinator().capacity());
        tiled.set("workers", Json::integer(tiled_workers));
        auto parallel = session.execute(tiled);
        check(parallel.at("result").dump() == first.at("result").dump() &&
                  parallel.at("execution").at("launched_threads").integer() ==
                      static_cast<std::int64_t>(tiled_workers - 1) &&
                  parallel.at("execution").at("completed_tiles").integer() == 6,
              "canonical P4 tiles disagree with serial action");
        ScientificSession cache_tiled(study);
        auto cache_query = Json::parse(tiled.dump());
        cache_query.set("backend", Json::string("cpu_cached_tiles"));
        const auto cold_cache = cache_tiled.execute(cache_query);
        const auto warm_cache = cache_tiled.execute(cache_query);
        check(cold_cache.at("result").dump() == first.at("result").dump() &&
                  !cold_cache.at("execution").at("result_reused").boolean() &&
                  cold_cache.at("execution")
                          .at("producer_execution")
                          .at("completed_tiles")
                          .integer() == 6 &&
                  warm_cache.at("result").dump() == first.at("result").dump() &&
                  warm_cache.at("execution").at("result_reused").boolean() &&
                  cache_tiled.native_cache_statistics().at("applications").integer() == 1,
              "parallel cache failed exact result/producer reuse");
        check(cache_tiled.execute(cached).at("execution").at("result_reused").boolean(),
              "serial request did not share canonical parallel result");
        for (std::size_t workers = 1; workers <= std::min<std::size_t>(2, tiled_workers);
             ++workers) {
            auto variant = Json::parse(tiled.dump());
            variant.set("workers", Json::integer(workers));
            check(session.execute(variant).at("result").dump() == first.at("result").dump(),
                  "worker count changed exact native values");
        }
        auto observed = std::make_shared<rfc::runtime::HardwareInventory>();
        observed->platform = "injected CPU contract";
        observed->cpu_model = "two declared test cores";
        observed->affinity_observed = true;
        observed->allowed_cpus = {{0, 0, 0, {}}, {1, 0, 1, {}}};
        rfc::runtime::CpuResourcePolicy resource_policy;
        resource_policy.worker_limit = 2;
        resource_policy.tile_limit = 20000;
        auto resources = std::make_shared<CpuResources>(*observed, resource_policy,
                                                        [observed] { return *observed; });
        auto single_observed = std::make_shared<rfc::runtime::HardwareInventory>(*observed);
        single_observed->allowed_cpus.resize(1);
        auto single_policy = resource_policy;
        single_policy.worker_limit = 1;
        auto single_resources = std::make_shared<CpuResources>(
            *single_observed, single_policy, [single_observed] { return *single_observed; });
        ScientificSession single_worker(study, single_resources);
        check(single_worker.execute(structure).at("result").at("total").dump() ==
                  structural.at("result").at("total").dump() &&
                  single_worker.execute(factored_structure).at("result").at("total").dump() ==
                      structural.at("result").at("total").dump() &&
                  single_resources->coordinator().available_workers() == 1,
              "structure display reserved its caller twice or retained a worker");
        ScientificSession coordinated(study, resources);
        auto parallel_controlled = Json::parse(controlled.dump());
        parallel_controlled.set("workers", Json::integer(2));
        check(coordinated.execute(parallel_controlled).at("result").at("fields").dump() ==
                  conditional.at("result").at("fields").dump(),
              "assigned control workers changed the exact field response");
        auto unassigned = Json::parse(parallel_controlled.dump());
        unassigned.set("workers", Json::integer(3));
        rejects([&] { coordinated.validate(unassigned); },
                "control preflight ignored assigned worker capacity");
        {
            auto occupied = resources->reserve(2);
            rejects([&] { coordinated.execute(parallel_controlled); },
                    "busy control request bypassed native worker admission");
        }
        check(resources->coordinator().available_workers() == 2,
              "prepared control request leaked workers");
        // A larger observed allocation may request more than the former 32
        // workers. The one field request launches only one useful tile.
        auto wide_observed = std::make_shared<rfc::runtime::HardwareInventory>(*observed);
        wide_observed->allowed_cpus.clear();
        for (unsigned i = 0; i < 48; ++i)
            wide_observed->allowed_cpus.push_back({i, 0, i, {}});
        auto wide_policy = resource_policy;
        wide_policy.worker_limit = wide_observed->allowed_cpus.size();
        auto wide_resources = std::make_shared<CpuResources>(
            *wide_observed, wide_policy, [wide_observed] { return *wide_observed; });
        ScientificSession wide_session(study, wide_resources);
        auto wide_query = Json::parse(parallel_controlled.dump());
        wide_query.set("workers", Json::integer(40));
        const auto wide_response = wide_session.execute(wide_query);
        check(wide_response.at("result").at("fields").dump() ==
                      conditional.at("result").at("fields").dump() &&
                  wide_response.at("result").at("execution").at("requested_workers").integer() ==
                      40 &&
                  wide_resources->coordinator().available_workers() == 48,
              "scientific execution retained a fixed worker ceiling");
        auto wide_boundary = Json::parse(tiled.dump());
        wide_boundary.set("workers", Json::integer(40));
        check(wide_session.execute(wide_boundary).at("result").dump() == first.at("result").dump(),
              "native field execution retained a fixed worker ceiling");
        auto wide_workflow = query(
            R"({"schema":"gadget.scientific.workflow.v1","id":"large-observed-allocation","limits":{"depth":2,"nodes":4,"model_calls":0,"parallel_workers":1},"definitions":[],"entry":{"id":"root","type":"parallel","workers":1,"steps":[{"id":"controlled","type":"query","query":{}}]}})");
        wide_workflow.at("entry").at("steps").at(0).set("query", wide_query);
        wide_workflow.set("source_digest", Json::string(study->digest()));
        check(execute_workflow(study, wide_workflow, {}, wide_resources).at("status").text() ==
                  "complete",
              "workflow demand retained a fixed worker ceiling");
        const auto automatic = query(
            R"({"operation":"native_boundary","field_set":"synthetic.assays","backend":"cpu_scheduled","tile_rows":1})");
        const auto scheduled = coordinated.execute(automatic);
        check(scheduled.at("result").dump() == first.at("result").dump() &&
                  !scheduled.at("execution").at("calibration_reused").boolean(),
              "canonical P5 calibration changed exact output");
        check(coordinated.execute(automatic).at("execution").at("calibration_reused").boolean(),
              "P5 calibration lost authentic prepared owner");
        {
            auto occupied = resources->coordinator().reserve(1);
            const auto fitting = coordinated.execute(automatic);
            check(fitting.at("execution").at("selected_workers").integer() == 1 &&
                      fitting.at("result").dump() == first.at("result").dump(),
                  "P5 did not select a measured fitting candidate");
        }
        {
            auto occupied = resources->coordinator().reserve(2);
            auto proposal = Json::object();
            proposal.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
            proposal.set("action", Json::string("query"));
            proposal.set("query", automatic);
            const auto busy = apply_scientific_proposal(coordinated, proposal);
            check(busy.at("status").text() == "busy" && !busy.has("certificate"),
                  "overbooked P5 call published a certificate");
        }
        check(resources->coordinator().available_workers() == 2, "P5 reservation leaked");
        auto overbooked = Json::parse(tiled.dump());
        overbooked.set("workers", Json::integer(3));
        rejects([&] { coordinated.execute(overbooked); },
                "unadmitted fixed CPU worker count executed");
        observed->cpu_model = "changed observations";
        rejects([&] { coordinated.execute(automatic); }, "stale P5 hardware/calibration admitted");
        observed->cpu_model = "two declared test cores";
        check(coordinated.execute(automatic).at("result").dump() == first.at("result").dump(),
              "failed P5 admission changed source or prior calibration");
        auto gpu = Json::parse(cached.dump());
        gpu.set("backend", Json::string("gpu"));
        rejects([&] { session.execute(gpu); }, "unqualified GPU backend admitted");
        std::vector<std::future<void>> readers;
        auto reader_query = Json::parse(tiled.dump());
        reader_query.set("workers", Json::integer(1));
        for (unsigned i = 0; i < 4; ++i) {
            if (readers.size() == tiled_workers) {
                for (auto &reader : readers)
                    reader.get();
                readers.clear();
            }
            readers.push_back(std::async(
                std::launch::async, [study, reader_query, expected = first.at("result").dump()] {
                    ScientificSession independent(study);
                    for (unsigned repeat = 0; repeat < 3; ++repeat)
                        check(independent.execute(reader_query).at("result").dump() == expected,
                              "concurrent immutable sessions changed native result");
                }));
        }
        for (auto &reader : readers)
            reader.get();
        auto plan = Json::object();
        plan.set("schema", Json::string("gadget.scientific.plan.v1"));
        plan.set("id", Json::string("two-visits"));
        plan.set("question", Json::string("Exact GET change on the declared population?"));
        plan.set("source_digest", Json::string(study->digest()));
        auto steps = Json::array();
        for (const auto &[id, q] : {std::pair{"baseline", old_moment},
                                    std::pair{"change", temporal}, std::pair{"native", cached}}) {
            auto step = Json::object();
            step.set("id", Json::string(id));
            step.set("query", q);
            steps.push(step);
        }
        plan.set("steps", steps);
        const auto complete = session.execute_plan(plan);
        check(complete.at("status").text() == "complete" && complete.at("steps").size() == 3,
              "automatic plan did not finish");
        auto bad = Json::parse(plan.dump());
        bad.at("steps").at(0).set("query", cached);
        bad.at("steps").at(2).set(
            "query",
            query(R"({"operation":"hodge","context":"baseline","field":"baseline.partial"})"));
        const auto before = session.native_cache_statistics().dump();
        rejects([&] { session.execute_plan(bad); }, "invalid late plan step admitted");
        check(before == session.native_cache_statistics().dump(),
              "plan executed a first step before admitting all steps");
        bad = Json::parse(plan.dump());
        bad.set("source_digest", Json::string("wrong"));
        rejects([&] { session.execute_plan(bad); }, "stale plan source binding admitted");
        auto mentions = Json::object();
        mentions.set("schema", Json::string("gadget.semantic.proposal.v1"));
        mentions.set("source_digest", Json::string(study->digest()));
        auto mention = Json::object();
        mention.set("entity", Json::string("gene:A"));
        mention.set("kind", Json::string("gene"));
        mention.set("polarity", Json::string("asserted"));
        const auto note = read_text(root / "design.txt");
        const auto start = note.find("A, B");
        auto spans = Json::array();
        auto span = Json::object();
        span.set("document", Json::string("design"));
        span.set("start", Json::integer(start));
        span.set("end", Json::integer(start + 1));
        span.set("text", Json::string("A"));
        spans.push(span);
        mention.set("evidence", spans);
        auto list = Json::array();
        list.push(mention);
        mentions.set("mentions", list);
        const auto receipt = session.propose_annotations(mentions);
        check(receipt.at("proposal").dump() == mentions.dump() &&
                  receipt.at("mentions").at(0).at("evidence").at(0).has("document_sha256"),
              "annotation changed input or lacks evidence hash");
        const auto immutable = study->digest();
        check(immutable == study->digest() &&
                  session.execute(old_moment).at("result").dump() == old.at("result").dump(),
              "semantic proposal changed mathematical source");
        auto invalid = Json::parse(mentions.dump());
        invalid.at("mentions").at(0).at("evidence").at(0).set("text", Json::string("B"));
        rejects([&] { session.propose_annotations(invalid); }, "fabricated source quote admitted");
        invalid = Json::parse(mentions.dump());
        invalid.at("mentions").at(0).set("entity", Json::string("gene:invented"));
        rejects([&] { session.propose_annotations(invalid); }, "undeclared NER identity admitted");
        auto detached = session.execute(old_moment);
        detached.at("design").set("population", Json::string("changed externally"));
        auto catalog = study->catalog();
        catalog.at("design").set("selection", Json::string("changed externally"));
        check(study->declaration().at("design").at("population").text() ==
                      "same synthetic population" &&
                  session.execute(old_moment).dump() == old.dump(),
              "returned JSON mutates immutable declaration");
        auto proposed = Json::object();
        proposed.set("schema", Json::string("gadget.scientific.agent.proposal.v1"));
        proposed.set("action", Json::string("plan"));
        proposed.set("plan", plan);
        check(apply_scientific_proposal(session, proposed).at("status").text() == "executed",
              "model plan bridge failed");
        proposed.set("code", Json::string("arbitrary"));
        check(apply_scientific_proposal(session, proposed).at("status").text() == "rejected",
              "unknown model member admitted");
        const auto original = read_text(root / "study.json");
        auto changed = Json::parse(original);
        changed.at("contexts")
            .at(1)
            .at("time")
            .set("population", Json::string("different population"));
        write_text(root / "study.json", changed.dump());
        ScientificSession incompatible(
            std::make_shared<const ScientificStudy>(root / "study.json"));
        rejects([&] { incompatible.execute(structure); }, "structure population mismatch admitted");
        rejects([&] { incompatible.execute(temporal); },
                "different populations treated as chronology");
        changed = Json::parse(original);
        changed.at("fields").at(4).set("units", Json::string("unrelated"));
        write_text(root / "study.json", changed.dump());
        ScientificSession units(std::make_shared<const ScientificStudy>(root / "study.json"));
        rejects([&] { units.execute(temporal); },
                "different units treated as identical observation type");
        changed = Json::parse(original);
        changed.at("cells").at(0).set("key", Json::string("999999"));
        write_text(root / "study.json", changed.dump());
        rejects([&] { ScientificStudy missing(root / "study.json"); },
                "absent native identity admitted by matching dimensions");
        write_text(root / "study.json", original);
        write_text(root / "design.txt", note + "edited");
        rejects([&] { ScientificStudy changed_document(root / "study.json"); },
                "edited evidence document admitted");
        write_text(root / "design.txt", note);
        auto output = Json::object();
        output.set("old", old);
        output.set("new", next);
        output.set("change", change);
        output.set("structure", structural);
        output.set("coverage", coverage);
        output.set("native", first);
        output.set("parallel", parallel);
        output.set("plan", complete);
        output.set("annotations", receipt);
        {
            rfc::core::RelationSignature signature;
            const auto type = signature.add_type("parallel.sources");
            rfc::core::RelationalComplexBuilder builder(signature);
            for (std::uint64_t e = 0; e < 2; ++e)
                builder.add_grade_one(
                    rfc::core::RelationKey{100 + e}, type,
                    {{1, true, true, Rational{1}}, {2, true, false, Rational{1}}});
            NativeStudy weighted(builder.freeze(),
                                 {rfc::field::ExactMetric::identity(2),
                                  rfc::field::ExactMetric::diagonal({Rational{1}, Rational{3}})});
            const auto witness = analyze_chain(
                weighted,
                query(
                    R"({"operation":"chain_certificate","context":"declared","grade":1,"coefficient_domain":"rational","values":["1","0"]})"));
            check(witness.at("sectors").at("gradient").at(0).text() == "3/4" &&
                      witness.at("sectors").at("gradient").at(1).text() == "1/4",
                  "weighted gradient coordinates");
            check(witness.at("sectors").at("harmonic").at(0).text() == "1/4" &&
                      witness.at("sectors").at("harmonic").at(1).text() == "-1/4",
                  "weighted harmonic coordinates");
            check(witness.at("quadrance").at("gradient").text() == "3/4" &&
                      witness.at("quadrance").at("harmonic").text() == "1/4",
                  "weighted exact quadrance partition");
            output.set("weighted_chain_witness", witness);
        }
        write_text(root / "results.json", output.dump() + "\n");
        std::cout
            << "PASS: general GET, grades, native P3/P4, time, plans, evidence and admission\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
