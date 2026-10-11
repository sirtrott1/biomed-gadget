#include "../src/genomic_detail.hpp"
#include <atomic>
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <iostream>

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        if (argc != 2)
            throw std::invalid_argument("expected output directory");
        const auto root = fs::absolute(argv[1]);
        fs::remove_all(root);
        fs::create_directories(root);
        const auto check = [](bool c, const char *m) {
            if (!c)
                throw std::logic_error(m);
        };
        const auto refuse = [&](auto action) {
            bool rejected = false;
            try {
                action();
            } catch (const std::exception &) {
                rejected = true;
            }
            check(rejected, "invalid request accepted");
        };
        const auto record = [&](const std::string &name, const std::string &bytes) {
            write_text(root / name, bytes);
            auto r = Json::object();
            r.set("path", Json::string(name));
            r.set("sha256", Json::string("sha256:" + sha256(bytes)));
            return r;
        };
        std::string axis = "id\n", values = "id\ta\tb\tc\td\n";
        for (std::size_t i = 0; i < 300; ++i) {
            const auto id = "g" + std::to_string(i);
            axis += id + "\n";
            values += id + "\t" + (i == 0 ? "2" : "0") + "\t" + (i == 299 ? "3" : "0") + "\t0\t" +
                      (i == 200 ? "NA" : "0") + "\n";
        }
        auto spec = Json::parse(
            R"({"schema":"gadget.assays.import.v1","id":"accounting","design":{"question":"Declared accounting","population":"authored observations","selection":"all","limitations":["fixture"]},"assays":[],"resources":[],"links":[]})");
        auto assay = Json::parse(
            R"({"id":"x","kind":"rna_seq","format":"dense_tsv","source_namespace":"1101","carrier":"features","quantity":"read_count","units":"counts","channel":"observation","method":"authored counts","arithmetic":"exact_rational","implicit_entries":"measured_zero","missing_tokens":["NA"],"partition_size":300})");
        assay.set("source", record("matrix.tsv", values));
        assay.set("features", record("features.tsv", axis));
        assay.set("observations", record("observations.tsv", "id\na\nb\nc\nd\n"));
        spec.at("assays").push(assay);
        write_text(root / "import.json", spec.dump());
        import_assays(root / "import.json", root / "workspace");
        AssayWorkspace workspace(root / "workspace/workspace.json");
        const auto member = [](const std::string &id, const std::string &coefficient) {
            auto m = Json::object();
            m.set("id", Json::string(id));
            m.set("coefficient", Json::string(coefficient));
            return m;
        };
        const auto group = [&](const std::string &id, const Json &members) {
            auto g = Json::object();
            g.set("id", Json::string(id));
            g.set("members", members);
            return g;
        };
        auto all = Json::array(), identities = Json::array();
        for (std::size_t i = 0; i < 300; ++i)
            all.push(member("g" + std::to_string(i), "1"));
        for (const auto *id : {"a", "b", "c", "d"}) {
            auto ms = Json::array();
            ms.push(member(id, "1"));
            identities.push(group(id, ms));
        }
        auto q = Json::parse(
            R"({"schema":"gadget.assays.query.v1","operation":"aggregate_assay","assay":"x","feature_groups":[],"observation_groups":[],"map_method":"authored map","map_evidence":"fixture identities","support_policy":"require_all_nonzero_inputs","workers":1})");
        q.at("feature_groups").push(group("sum", all));
        q.set("observation_groups", identities);
        auto result = workspace.execute(q).at("result");
        check(result.at("rows").dump() == Json::parse(R"([["2","3","0","NA"]])").dump(),
              "support-aware native sum");
        check(result.at("execution").at("reduction_owner").text() ==
                  "rfc::database::NativeTensorWindow",
              "native route missing");
        auto parallel = Json::parse(q.dump());
        parallel.set("workers", Json::integer(2));
        check(workspace.execute(parallel).at("result").at("rows").dump() ==
                  result.at("rows").dump(),
              "parallel parity");
        auto weighted = Json::parse(q.dump());
        auto groups = Json::array(), ms = Json::array();
        ms.push(member("a", "2"));
        ms.push(member("b", "-1"));
        groups.push(group("signed", ms));
        ms = Json::array();
        ms.push(member("c", "1"));
        ms.push(member("d", "0"));
        groups.push(group("zero dependency", ms));
        weighted.set("observation_groups", groups);
        check(workspace.execute(weighted).at("result").at("rows").dump() ==
                  Json::parse(R"([["1","0"]])").dump(),
              "weighted outer map/zero coefficient");
        auto partial = Json::parse(q.dump());
        groups = Json::array();
        ms = Json::array();
        ms.push(member("g0", "2"));
        ms.push(member("g299", "-1"));
        groups.push(group("partial", ms));
        partial.set("feature_groups", groups);
        auto selected = workspace.execute(partial).at("result");
        check(selected.at("rows").dump() == Json::parse(R"([["4","-3","0","0"]])").dump() &&
                  selected.at("execution").at("source_owner_reused").boolean() &&
                  !selected.at("execution").at("raw_matrix_reparsed").boolean() &&
                  !selected.at("execution").has("reduction_owner"),
              "selected native matrix action");
        auto bad = Json::parse(q.dump());
        bad.set("count_policy", Json::string("disjoint_count_sum"));
        refuse([&] { (void)workspace.execute(bad); });
        bad = Json::parse(q.dump());
        bad.at("feature_groups").at(0).at("members").push(member("g0", "1"));
        refuse([&] { (void)workspace.execute(bad); });
        // Actual field pins own their immutable publication after the caller dies.
        auto owner = workspace.study("x");
        const auto pin = owner->field_accession("all.a");
        ScientificSession session(owner);
        auto window = Json::parse(
            R"({"operation":"tensor_window","field":"all.a","lower":[0],"upper":[0],"available_sequence":"0"})");
        for (std::size_t i = 0; i <= 300; ++i) {
            window.set("upper", Json::parse("[" + std::to_string(i) + "]"));
            const auto answer = session.execute(window).at("result");
            check(!answer.at("preparation_reused").boolean() &&
                      answer.at("preparation_retained").boolean(),
                  "distinct retention");
            if (i)
                check(answer.at("sum").at("rational").text() == "2", "window coefficient sum");
        }
        check(
            session.tensor_window_preparation_statistics().at("retained_preparations").integer() ==
                301,
            "ceiling remains");
        const auto retained = session.execute(window);
        check(retained.at("result").at("preparation_reused").boolean(), "prepared reuse");
        auto once = Json::parse(window.dump());
        once.set("preparation", Json::string("one_shot"));
        check(!session.execute(once).at("result").at("preparation_retained").boolean(),
              "one-shot retention");
        check(
            session.tensor_window_preparation_statistics().at("retained_preparations").integer() ==
                301,
            "one-shot changed cache");
        rfc::runtime::ExecutionHostPolicy policy;
        const auto max = std::numeric_limits<std::size_t>::max();
        policy.limits = {max, max, max, max, max};
        policy.limits.scalar_work = 0;
        refuse([&] {
            (void)session.execute_tensor_window(window, rfc::runtime::ExecutionLease(policy));
        });
        policy.limits.scalar_work = max;
        auto cancelled = std::make_shared<std::atomic_bool>(true);
        refuse([&] {
            (void)session.execute_tensor_window(window,
                                                rfc::runtime::ExecutionLease(policy, cancelled));
        });
        session.clear_tensor_window_preparations();
        check(
            session.tensor_window_preparation_statistics().at("retained_preparations").integer() ==
                0,
            "cache not released");
        check(retained.at("result").at("sum").at("rational").text() == "2",
              "clear invalidated certificate");
        check(!session.execute(window).at("result").at("preparation_reused").boolean(),
              "clear failed to reset reuse");
        // The transposed storage profile can reduce observation groups with the
        // same declared-map and support semantics, preserving feature identities.
        spec.at("assays").at(0).set("carrier", Json::string("observations"));
        spec.at("assays").at(0).set("source_namespace", Json::string("2101"));
        write_text(root / "transposed.json", spec.dump());
        import_assays(root / "transposed.json", root / "transposed");
        AssayWorkspace reverse(root / "transposed/workspace.json");
        auto gene = Json::parse(q.dump());
        groups = Json::array();
        for (const auto *id : {"g0", "g200", "g299"}) {
            ms = Json::array();
            ms.push(member(id, "1"));
            groups.push(group(id, ms));
        }
        gene.set("feature_groups", groups);
        groups = Json::array();
        ms = Json::array();
        for (const auto *id : {"a", "b", "c", "d"})
            ms.push(member(id, "1"));
        groups.push(group("library", ms));
        gene.set("observation_groups", groups);
        gene.set("count_policy", Json::string("disjoint_count_sum"));
        result = reverse.execute(gene).at("result");
        check(result.at("rows").dump() == Json::parse(R"([["2"],["NA"],["3"]])").dump() &&
                  result.at("count_contract_verified").boolean(),
              "native count contract");
        // Disjoint stored carriers address the original measurement coordinates.
        // Their native complexes remain separate throughout the actions.
        spec.at("assays").at(0).set("carrier", Json::string("features"));
        spec.at("assays").at(0).set("partition_size", Json::integer(127));
        spec.at("assays").at(0).set("source_namespace", Json::string("3101"));
        write_text(root / "partitioned.json", spec.dump());
        import_assays(root / "partitioned.json", root / "partitioned");
        AssayWorkspace partitioned(root / "partitioned/workspace.json");
        const auto reference = assay_detail::read_matrix(root, spec.at("assays").at(0));
        spec.at("assays").at(0).set("carrier", Json::string("observations"));
        spec.at("assays").at(0).set("partition_size", Json::integer(2));
        spec.at("assays").at(0).set("source_namespace", Json::string("4101"));
        write_text(root / "partitioned-reverse.json", spec.dump());
        import_assays(root / "partitioned-reverse.json", root / "partitioned-reverse");
        AssayWorkspace partitioned_reverse(root / "partitioned-reverse/workspace.json");
        auto rational = Json::parse(partial.dump());
        groups = Json::array();
        ms = Json::array();
        ms.push(member("g0", "1/3"));
        ms.push(member("g299", "-2/5"));
        groups.push(group("rational", ms));
        rational.set("feature_groups", groups);
        groups = Json::array();
        ms = Json::array();
        ms.push(member("a", "5/7"));
        ms.push(member("b", "-3/11"));
        groups.push(group("contrast", ms));
        rational.set("observation_groups", groups);
        check(workspace.execute(rational).at("result").at("rows").dump() ==
                  Json::parse(R"([["928/1155"]])").dump(),
              "exact two-axis rational contrast");
        auto absent = Json::parse(partial.dump());
        groups = Json::array();
        ms = Json::array();
        ms.push(member("g200", "1"));
        groups.push(group("reported", ms));
        absent.set("feature_groups", groups);
        auto zero = Json::parse(absent.dump());
        groups = Json::array();
        ms = Json::array();
        ms.push(member("g200", "0"));
        groups.push(group("zero map", ms));
        zero.set("feature_groups", groups);
        check(workspace.execute(zero).at("result").at("rows").dump() ==
                  Json::parse(R"([["0","0","0","0"]])").dump(),
              "empty nonzero dependency set");
        check(workspace.execute(absent).at("result").at("rows").dump() ==
                  Json::parse(R"([["0","0","0","NA"]])").dump(),
              "supported zero versus absence");
        for (const auto &query : {q, weighted, partial, gene, rational, absent, zero}) {
            auto resources = CpuResources::local();
            const auto expected = [&] {
                auto reservation = resources->reserve(1);
                return assay_detail::aggregate_analysis(reference, spec.at("assays").at(0), query,
                                                        reservation);
            }();
            for (const auto *source : {&workspace, &reverse, &partitioned, &partitioned_reverse}) {
                for (const auto workers : {1, 2}) {
                    auto action = Json::parse(query.dump());
                    action.set("workers", Json::integer(workers));
                    const auto answer = source->execute(action).at("result");
                    for (const auto *name :
                         {"rows", "feature_ids", "observation_ids", "quantity", "units",
                          "feature_map_nonzeros", "observation_map_nonzeros"})
                        check(answer.at(name).dump() == expected.at(name).dump(),
                              "original-axis map parity");
                    const auto execution = answer.at("execution");
                    if (execution.has("sparse_map_owner")) {
                        check(execution.at("sparse_map_owner").text() ==
                                      "rfc::numeric::PreparedSparseExactMap" &&
                                  execution.at("intermediate_coordinates").integer() == 0,
                              "sparse native map execution retained a dense intermediate");
                    } else {
                        check(execution.has("reduction_owner"),
                              "native numerical map did not use canonical sparse execution");
                    }
                    check(execution.at("source_owner").text() ==
                                  "rfc::storage::PinnedNativeAccession" &&
                              !execution.at("raw_matrix_reparsed").boolean(),
                          "native owner receipt");
                    const auto bindings = execution.at("source_bindings");
                    check(bindings.size() == static_cast<std::size_t>(
                                                 execution.at("source_partition_count").integer()),
                          "incomplete source bindings");
                    for (std::size_t p = 0; p < bindings.size(); ++p)
                        check(bindings.at(p).at("source_digest").text() ==
                                  source->study("x", p)->digest(),
                              "receipt names a foreign publication");
                }
            }
        }
        const auto reduced = partitioned.execute(q).at("result").at("execution");
        check(reduced.at("reduction_composition_owner").text() == "rfc::numeric::ExactMatrix" &&
                  reduced.at("visited_sparse_entries").integer() == 2,
              "native packet sum composition");
        check(partitioned.execute(absent)
                      .at("result")
                      .at("execution")
                      .at("required_native_partitions")
                      .integer() == 1,
              "unselected source packet supplied fields");
        check(partitioned.execute(zero)
                      .at("result")
                      .at("execution")
                      .at("required_native_fields")
                      .integer() == 0,
              "zero action retained unused field pins");
        auto design = Json::parse(
            R"({"schema":"gadget.assays.query.v1","operation":"design_aggregation","assay":"x","feature_design":{"mode":"identity"},"observation_design":{"mode":"groups","missing_tokens":[],"unassigned_policy":"error","overlap_policy":"error","groups":[{"id":"library","coefficient":"1","where":[{"column":"id","values":["a","b","c","d"]}]}]},"map_method":"dictionary source identities","map_evidence":"fixture original axis","support_policy":"require_all_nonzero_inputs","count_policy":"disjoint_count_sum","workers":1})");
        auto resources = CpuResources::local();
        const auto expected_design = [&] {
            auto reservation = resources->reserve(1);
            return assay_detail::design_aggregation(reference, spec.at("assays").at(0), design,
                                                    reservation);
        }();
        for (const auto *source : {&workspace, &reverse, &partitioned, &partitioned_reverse}) {
            const auto answer = source->execute(design).at("result");
            check(answer.at("rows").dump() == expected_design.at("rows").dump() &&
                      answer.at("design").dump() == expected_design.at("design").dump() &&
                      !answer.at("execution").at("raw_matrix_reparsed").boolean(),
                  "native dictionary design parity");
        }
        std::size_t derived_index = 0;
        for (const auto *source : {&workspace, &reverse, &partitioned, &partitioned_reverse}) {
            for (const auto &query : {gene, design, zero}) {
                const auto retained = source->execute_record(admit_record_view(query));
                const auto fields = retained.at("result").at("result");
                check(fields.at("schema").text() == "gadget.aggregation.fields.v1" &&
                          !fields.has("rows") &&
                          fields.at("feature_map").kind() == NativeRecordKind::Integer &&
                          fields.at("observation_map").kind() == NativeRecordKind::Integer,
                      "aggregation retained a numerical display instead of native fields and maps");
                const auto map = [&](const char *name) -> const rfc::storage::TemporalTensor & {
                    const auto key = std::stoull(fields.at(name).integer().get_str());
                    const auto endpoint = retained.endpoint();
                    const auto found = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                        [&](const auto &tensor) { return tensor.key == key; });
                    check(found != endpoint->tensors.end() &&
                          found->role == rfc::storage::TemporalTensorRole::Observation &&
                          found->tensor.support().storage() ==
                              rfc::numeric::StructuralSupport::Storage::Full,
                          "aggregation map carrier or full support differs");
                    return *found;
                };
                const auto &feature_map = map("feature_map");
                const auto &observation_map = map("observation_map");
                check(feature_map.input == 2 && feature_map.output == 1 &&
                      observation_map.input == 3 && observation_map.output == 4,
                      "aggregation maps lost their declared source axes");
                const auto compiled = query.at("operation").text() == "design_aggregation"
                    ? expected_design.at("design").at("compiled_request") : query;
                const auto check_map = [&](const auto &tensor, const auto &groups, const auto &ids) {
                    std::map<std::size_t, Rational> expected;
                    for (std::size_t i = 0; i < groups.size(); ++i) {
                        const auto members = groups.at(i).at("members");
                        for (std::size_t j = 0; j < members.size(); ++j) {
                            const auto member = members.at(j);
                            const auto position = std::find(ids.begin(), ids.end(), member.at("id").text());
                            check(position != ids.end(), "oracle group member missing");
                            const Rational value(member.at("coefficient").text());
                            if (value != 0)
                                expected.emplace(i * ids.size() + (position - ids.begin()), value);
                        }
                    }
                    check(tensor.domain().extents() == std::vector<std::size_t>{groups.size(), ids.size()} &&
                          tensor.stored_nonzeros() == expected.size(), "native group map shape differs");
                    for (const auto &entry : tensor.values().entries())
                        check(expected.at(entry.index) == entry.value, "native group coefficient differs");
                };
                check_map(feature_map.tensor, compiled.at("feature_groups"), reference.features.ids);
                check_map(observation_map.tensor, compiled.at("observation_groups"), reference.observations.ids);
                const auto encoded = retained.encoded_tensor(feature_map.key);
                check(encoded && encoded->tensor().support().storage() ==
                      rfc::numeric::StructuralSupport::Storage::Full,
                      "native map encoding expanded structural support");
                const auto image = render_record_view(retained.at("result")).at("result");
                check(image.at("rows").dump() ==
                          source->execute(query).at("result").at("rows").dump(),
                      "aggregation field view differs");
                if (query.at("operation").text() == "aggregate_assay" &&
                    query.at("feature_groups").size() == 1)
                    check(fields.at("feature_groups")
                                      .at(0)
                                      .at("members")
                                      .at(0)
                                      .at("coefficient")
                                      .text() == "0" &&
                              feature_map.tensor.stored_nonzeros() == 0,
                          "declared zero map membership was lost");
            }
            auto publication = Json::object();
            publication.set("schema", Json::string("gadget.assays.materialize.v1"));
            publication.set("id", Json::string("summed"));
            publication.set("source_namespace",
                            Json::string(std::to_string(50000 + derived_index)));
            publication.set("design", spec.at("design"));
            publication.set("recipe", gene);
            publication.set("channel", Json::string("observation"));
            const auto destination = root / ("native-sum-" + std::to_string(derived_index++));
            const auto admitted =
                source->materialize_workspace(admit_record_view(publication), destination);
            const auto owner = admitted->study("summed");
            check(owner->field_accession("all.library").tensor()->state(0) == Rational{2} &&
                      !owner->field_accession("all.library").tensor()->state(1) &&
                      owner->field_accession("all.library").tensor()->state(2) == Rational{3},
                  "count sum support or original feature order differs");
            check(owner->declaration().at("fields").at(0).at("quantity").text() == "read_count" &&
                      owner->declaration().at("fields").at(0).at("units").text() == "counts",
                  "count sum meaning differs");
            check(owner->cell_identifier({0}, 0) == "g0" && owner->cell_identifier({0}, 2) == "g299" &&
                      owner->cell_identity({0}, "g299") == rfc::storage::StableCellKey{{0}, 3},
                  "native feature names differ");
            auto annotation = owner->native_publication().at("study").unbound_value();
            auto &annotation_members = std::get<NativeRecordValue::Object>(annotation.data).members;
            for (auto &[name, value] : annotation_members) {
                if (name == "native") value = NativeRecordValue{};
                if (name != "cell_bindings") continue;
                auto &binding = std::get<NativeRecordValue::Object>(
                    std::get<NativeRecordValue::Array>(value.data).elements[0].data);
                for (auto &[field, ids] : binding.members)
                    if (field == "ids")
                        std::get<NativeRecordValue::Array>(ids.data).elements[1] =
                            NativeRecordValue(std::string("g0"));
            }
            std::erase_if(annotation_members, [](const auto &entry) { return entry.first == "native"; });
            refuse([&] {
                ScientificStudy invalid(ScientificStudy::publish_annotations(owner->producer_publication(),
                    annotation, owner->native_publication().at("documents").unbound_value()));
            });
            check(!fs::exists(destination / "matrix.tsv") &&
                      !fs::exists(destination / "features.tsv"),
                  "native sum used a generated matrix input");
            const auto binary = owner->producer_publication().encode();
            check(NativeRecord::decode(binary).encode() == binary,
                  "aggregation canonical bytes differ");
            auto chained = Json::parse(gene.dump());
            chained.set("assay", Json::string("summed"));
            auto observations = Json::array(), members = Json::array();
            members.push(member("library", "1"));
            observations.push(group("again", members));
            chained.set("observation_groups", observations);
            check(admitted->execute(chained).at("result").at("rows").dump() ==
                      Json::parse(R"([["2"],["NA"],["3"]])").dump(),
                  "chained native sum reparsed or changed a field");
            const auto relocated = root / ("relocated-sum-" + std::to_string(derived_index));
            fs::rename(destination, relocated);
            AssayWorkspace reopened(relocated / "workspace.rfcxs");
            check(reopened.study("summed")->field_accession("all.library").tensor()->state(2) ==
                      Rational{3},
                  "relocated aggregation lost its field owner");
            check(reopened.replay_derivation("summed",
                                             root / ("sum-replay-" + std::to_string(derived_index)))
                          .at("status")
                          .text() == "PASS",
                  "aggregation replay differs");
        }
        std::vector<assay_detail::NativeAssayPartition> packets{
            {0, 127, partitioned.study("x", 0)},
            {127, 254, partitioned.study("x", 1)},
            {254, 300, partitioned.study("x", 2)}};
        const auto reject_packets = [&](const auto &altered) {
            refuse([&] {
                auto reservation = resources->reserve(1);
                (void)assay_detail::native_aggregate_analysis(
                    reference.features, reference.observations, altered, spec.at("assays").at(0), q,
                    reservation);
            });
        };
        // Use the matching declared orientation; wrong ranges, gaps, repeated
        // publications and unrelated local source identities still reject.
        spec.at("assays").at(0).set("carrier", Json::string("features"));
        auto altered = packets;
        altered[1].begin = 126;
        reject_packets(altered);
        altered = packets;
        altered[1].begin = 128;
        reject_packets(altered);
        altered = packets;
        altered[2].end = 301;
        reject_packets(altered);
        altered = packets;
        altered.pop_back();
        reject_packets(altered);
        altered = packets;
        altered[1].owner = altered[0].owner;
        reject_packets(altered);
        altered = packets;
        std::swap(altered[0].owner, altered[1].owner);
        reject_packets(altered);
        // The default follows the original carrier, including sizes beyond the
        // former storage default. An explicit partition request remains valid.
        auto wide = Json::parse(spec.dump());
        auto wide_assay = Json::parse(wide.at("assays").at(0).dump());
        json_object_object_del(wide_assay.raw(), "partition_size");
        wide_assay.set("source_namespace", Json::string("5101"));
        std::string wide_axis = "id\n", wide_matrix = "id\ta\tb\n";
        for (std::size_t i = 0; i < 1025; ++i) {
            const auto id = "original_" + std::to_string(i);
            wide_axis += id + "\n";
            wide_matrix += id + "\t0\t0\n";
        }
        wide_assay.set("source", record("wide.tsv", wide_matrix));
        wide_assay.set("features", record("wide-features.tsv", wide_axis));
        wide_assay.set("observations", record("wide-observations.tsv", "id\na\nb\n"));
        auto wide_assays = Json::array();
        wide_assays.push(wide_assay);
        wide.set("assays", wide_assays);
        write_text(root / "wide-default.json", wide.dump());
        import_assays(root / "wide-default.json", root / "wide-default");
        AssayWorkspace complete_default(root / "wide-default/workspace.json");
        check(complete_default.catalog().at("assays").at(0).at("partitions").size() == 1 &&
                  complete_default.study("x")
                          ->endpoint()
                          ->source.grade(rfc::core::Grade{0})
                          .cell_count() == 1025,
              "default fragmented the original source carrier");
        wide_assay.set("partition_size", Json::integer(256));
        wide_assays = Json::array();
        wide_assays.push(wide_assay);
        wide.set("assays", wide_assays);
        write_text(root / "wide-partitioned.json", wide.dump());
        import_assays(root / "wide-partitioned.json", root / "wide-partitioned");
        AssayWorkspace requested_parts(root / "wide-partitioned/workspace.json");
        check(requested_parts.catalog().at("assays").at(0).at("partitions").size() == 5,
              "explicit packet size was ignored");
        // A zero coefficient still has to name a valid source identity; two
        // authored members cannot silently coalesce even when their sum is zero.
        bad = Json::parse(zero.dump());
        bad.at("feature_groups").at(0).at("members").at(0).set("id", Json::string("foreign"));
        refuse([&] { (void)workspace.execute(bad); });
        bad = Json::parse(partial.dump());
        bad.at("feature_groups").at(0).at("members").push(member("g0", "-2"));
        refuse([&] { (void)workspace.execute(bad); });
        // Detached receipts do not mutate the publication or subsequent answers.
        selected.set("rows", Json::parse(R"([["999"]])"));
        check(workspace.execute(partial).at("result").at("rows").at(0).at(0).text() == "4",
              "caller result altered the original field");
        const auto survivor = [&] {
            ScientificStudy transient(root / "workspace/assays/0/partitions/0/study.json");
            return transient.field_accession("all.a");
        }();
        check(survivor.tensor()->state(0) == Rational{2} && pin.tensor()->state(0) == Rational{2},
              "pin lifetime");
        std::cout << "support, maps, both carrier orientations, retention, leases and import "
                     "contracts passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
