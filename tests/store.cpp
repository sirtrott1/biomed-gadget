#include <gadget/execution.hpp>
#include <gadget/store.hpp>
#include <gadget/native_record.hpp>
#include <iostream>
#include <sys/wait.h>

namespace {
using namespace gadget;
namespace native = rfc::storage;
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&operation, const char *message) {
    bool failed = false;
    try {
        operation();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, message);
}
Json request(const Json &root, std::uint64_t transaction) {
    auto value = Json::object();
    value.set("schema", Json::string("gadget.store.publish.v1"));
    value.set("transaction", Json::string(std::to_string(transaction)));
    value.set("expected", Json::parse(root.dump()));
    return value;
}
std::uint64_t key(const Workspace &workspace, const std::string &name) {
    for (const auto &tensor : workspace.native_snapshot->endpoint()->tensors)
        if (tensor.interpretation == name)
            return tensor.key;
    throw std::runtime_error("missing field key");
}
} // namespace
int main() {
    try {
        const fs::path root = GADGET_STORE_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root);
        fs::copy(GADGET_TEST_FIXTURE, root / "input", fs::copy_options::recursive);
        auto original = load_pager(root / "input");
        import_measurements(original, root / "input/measurements.json");
        save_snapshot(original, root / "original", true);
        auto config = Json::parse(read_text(root / "input/measurements.json"));
        auto added = Json::parse(config.at("fields").at(0).dump());
        added.set("id", Json::string("AD.A"));
        config.at("fields").push(added);
        write_text(root / "input/measurements.json", config.dump());
        write_text(root / "input/effects.csv", "gene,log2FC,padj\nA,0,0.125\nB,NA,0\nC,0.5,NA\n");
        auto updated = load_pager(root / "input");
        import_measurements(updated, root / "input/measurements.json");
        save_snapshot(updated, root / "updated", true);
        fs::rename(root / "input", root / "raw-unavailable");
        // Optional views may use different whitespace without changing the native study.
        fs::copy(root / "original", root / "spaced", fs::copy_options::recursive);
        for (const auto *name : {"snapshot.json", "metadata.json"})
            write_text(root / "spaced" / name, "\n" + read_text(root / "spaced" / name) + "\n");
        check(load_snapshot(root / "spaced").fields == original.fields,
              "valid optional view fixture changed native fields");
        const auto spaced = initialize_store(root / "spaced", root / "spaced-store");
        for (const auto *name : {"snapshot.json", "metadata.json", "study.rfcxs"})
            check(read_text(root / "spaced-store/studies/0000000000000000" / name) ==
                      read_text(root / "spaced" / name),
                  "preparation changed original artifact bytes");
        check(load_store(root / "spaced-store").fields == original.fields &&
                  load_store(root / "spaced-store").field_presence == original.field_presence,
              "optional view whitespace changed stored coefficients or support");
        const auto spaced_request = request(spaced.at("root"), 41);
        rejects([&] {
            publish_store(root / "spaced-store", root / "spaced", spaced_request, {},
                          [](auto stage) {
                              if (stage == StudyPreparationStage::Renamed)
                                  throw std::runtime_error("interrupt prepared provenance");
                          });
        }, "prepared interruption did not occur");
        const auto spaced_completion = read_text(root / "spaced/snapshot.json");
        const auto spaced_metadata = read_text(root / "spaced/metadata.json");
        for (const auto *name : {"snapshot.json", "metadata.json"}) {
            const auto original_bytes = read_text(root / "spaced" / name);
            fs::remove(root / "spaced" / name);
            rejects([&] { publish_store(root / "spaced-store", root / "spaced", spaced_request); },
                    "prepared retry lost an original optional artifact");
            write_text(root / "spaced" / name, original_bytes + "\n");
            rejects([&] { publish_store(root / "spaced-store", root / "spaced", spaced_request); },
                    "prepared retry changed original optional bytes");
            write_text(root / "spaced" / name, original_bytes);
            check(store_status(root / "spaced-store", 41).at("status").text() == "absent",
                  "refused prepared retry published a native transaction");
        }
        publish_store(root / "spaced-store", root / "spaced", spaced_request);
        check(read_text(root / "spaced-store/studies/0000000000000029/snapshot.json") ==
                  spaced_completion &&
                  read_text(root / "spaced-store/studies/0000000000000029/metadata.json") ==
                  spaced_metadata,
              "exact prepared retry changed optional bytes");
        fs::copy(root / "original", root / "native-only", fs::copy_options::recursive);
        fs::remove(root / "native-only/snapshot.json");
        fs::remove(root / "native-only/metadata.json");
        const auto native_only = initialize_store(root / "native-only", root / "native-only-store");
        const auto native_only_request = request(native_only.at("root"), 44);
        publish_store(root / "native-only-store", root / "native-only", native_only_request);
        check(load_store(root / "native-only-store").fields == original.fields,
              "native-only study preparation failed");
        write_text(root / "native-only/metadata.json", read_text(root / "original/metadata.json"));
        rejects([&] {
            publish_store(root / "native-only-store", root / "native-only", native_only_request);
        }, "prepared retry added an optional artifact");
        fs::remove(root / "native-only/metadata.json");
        check(publish_store(root / "native-only-store", root / "native-only", native_only_request)
                      .at("receipt").at("root").at("sequence").text() == "1",
              "exact native-only retry did not retain its receipt");
        fs::create_symlink("absent-view", root / "native-only/snapshot.json");
        rejects([&] { initialize_store(root / "native-only", root / "dangling-view-store"); },
                "dangling optional view bypassed regular-file admission");
        check(!fs::exists(root / "dangling-view-store"),
              "refused snapshot created a store directory");
        const auto initial = initialize_store(root / "original", root / "store");
        const auto initial_root = Json::parse(initial.at("root").dump());
        rejects([&] { initialize_store(root / "updated", root / "store"); }, "genesis overwritten");
        const auto reader = load_store(root / "store");
        auto pinned = reader.native_snapshot->field("AD.log2FC");
        const auto pinned_values = pinned.tensor()->values().dense();
        const auto old_key = key(reader, "AD.log2FC");
        const auto update_request = request(initial_root, 1);
        const auto receipt = publish_store(root / "store", root / "updated", update_request);
        check(receipt.at("root").at("sequence").text() == "1", "publication not retained");
        const auto next = load_store(root / "store");
        check(next.fields == updated.fields && next.field_presence == updated.field_presence &&
                  key(next, "AD.log2FC") == old_key &&
                  key(load_snapshot(root / "updated"), "AD.log2FC") != old_key,
              "durable field keys or exact support changed during portable key remap");
        check(pinned.tensor()->values().dense() == pinned_values &&
                  pinned.endpoint()->sequence == 0 &&
                  next.native_snapshot->endpoint()->sequence == 1,
              "old owning pin changed after publication");
        const auto a = next.gene_index.at("A"), b = next.gene_index.at("B");
        check(next.native_snapshot->field("AD.log2FC").tensor()->support().contains(a) &&
                  next.fields.at("AD.log2FC")[a] == 0 &&
                  !next.native_snapshot->field("AD.log2FC").tensor()->support().contains(b),
              "observed zero conflated with missing");
        const auto query = Json::parse(
            R"({"operation":"typed_moment","context":"AD","realization":"support","field":"AD.log2FC","other_field":"AD.G","order":2})");
        auto result = execute_query(next, query);
        check(result.at("result").dump() == execute_query(updated, query).at("result").dump() &&
                  result.at("store").at("root").dump() == receipt.at("root").dump(),
              "store query or certificate incorrect");
        // Recapture remains portable and does not leak the service publication declarations.
        save_snapshot(next, root / "recaptured", true);
        check(read_text(root / "recaptured/study.rfcxs") == read_text(root / "updated/study.rfcxs"),
              "durable recapture changed portable endpoint bytes");
        check(store_status(root / "store", 1).at("status").text() == "committed" &&
                  store_status(root / "store", 99).at("status").text() == "absent",
              "transaction status incorrect");
        const auto second_request = request(receipt.at("root"), 2);
        publish_store(root / "store", root / "updated", second_request);
        const auto retry = publish_store(root / "store", root / "updated", update_request);
        check(retry.at("root").at("sequence").text() == "2" &&
                  retry.at("receipt").dump() == receipt.at("receipt").dump(),
              "historical exact retry did not return original receipt");
        rejects([&] { publish_store(root / "store", root / "original", update_request); },
                "changed retry admitted");
        rejects([&] { publish_store(root / "store", root / "updated", request(initial_root, 3)); },
                "stale new transaction admitted");
        auto wrong = request(receipt.at("root"), 1);
        rejects([&] { publish_store(root / "store", root / "updated", wrong); },
                "changed retry predecessor admitted");
        const auto latest_root = store_status(root / "store").at("root");
        rejects([&] { publish_store(root / "store", root / "original", request(latest_root, 3)); },
                "field removal admitted");
        auto altered = load_pager(root / "raw-unavailable");
        import_measurements(altered, root / "raw-unavailable/measurements.json");
        auto changed_meaning = Json::parse(altered.field_semantics.at("AD.log2FC").dump());
        changed_meaning.set("units", Json::string("different units"));
        altered.field_semantics["AD.log2FC"] = changed_meaning;
        refresh_source_digest(altered);
        save_snapshot(altered, root / "different-meaning", true);
        rejects(
            [&] {
                publish_store(root / "store", root / "different-meaning", request(latest_root, 3));
            },
            "changed field meaning admitted");
        auto malformed = request(latest_root, 3);
        malformed.set("arbitrary", Json::boolean(true));
        rejects([&] { publish_store(root / "store", root / "updated", malformed); },
                "unknown request member admitted");
        for (const auto *text : {"0", "01", "-1", "18446744073709551616", " 1", "+1"})
            rejects([&] { store_transaction(text); }, "invalid transaction admitted");
        check(store_transaction("18446744073709551615") == UINT64_MAX, "uint64 maximum rejected");
        auto nul_request = request(latest_root, 3);
        nul_request.set("transaction", Json::string(std::string("3\0suffix", 8)));
        rejects([&] { publish_store(root / "store", root / "updated", nul_request); },
                "NUL-truncated transaction admitted");
        // Equal native ordinals and shape do not authorize changed gene labels.
        fs::copy(root / "raw-unavailable", root / "relabeled-input", fs::copy_options::recursive);
        for (const auto &entry : fs::recursive_directory_iterator(root / "relabeled-input")) {
            if (entry.path().extension() != ".csv")
                continue;
            auto table = read_csv(read_text(entry.path()));
            std::string text;
            for (std::size_t i = 0; i < table.header.size(); ++i)
                text += (i ? "," : "") + csv_cell(table.header[i]);
            text += '\n';
            for (const auto &row : table.rows) {
                for (std::size_t i = 0; i < row.cells.size(); ++i)
                    text += (i ? "," : "") + csv_cell(row.cells[i] == "A" ? "AA" : row.cells[i]);
                text += '\n';
            }
            write_text(entry.path(), text);
        }
        auto relabeled = load_pager(root / "relabeled-input");
        import_measurements(relabeled, root / "relabeled-input/measurements.json");
        check(native::encode_rfcx(*relabeled.complex).bytes ==
                  native::encode_rfcx(*updated.complex).bytes,
              "equal-native-geometry label fixture invalid");
        save_snapshot(relabeled, root / "relabeled", true);
        rejects([&] { publish_store(root / "store", root / "relabeled", request(latest_root, 3)); },
                "equal-shaped changed biomedical identities admitted");
        {
            auto coordinator = native::DurableNativeStore::open(root / "store/history");
            rejects([&] { store_status(root / "store"); }, "exclusive coordinator bypassed");
        }
        // Inject exceptions and actual process exits at every canonical append cut.
        std::vector<native::DurableEvent> cuts;
        const auto cuts_initial = initialize_store(root / "original", root / "enumerate");
        publish_store(root / "enumerate", root / "updated", request(cuts_initial.at("root"), 7),
                      [&](auto event) { cuts.push_back(event); });
        check(cuts.size() >= 20, "native publication cuts missing");
        std::size_t unknowns = 0;
        for (std::size_t cut = 0; cut < cuts.size(); ++cut) {
            for (bool process_exit : {false, true}) {
                const auto path =
                    root / ((process_exit ? "exit-" : "exception-") + std::to_string(cut));
                const auto genesis = initialize_store(root / "original", path);
                const auto frozen_request = request(genesis.at("root"), 42);
                auto operation = [&] {
                    std::size_t count = 0;
                    publish_store(path, root / "updated", frozen_request, [&](auto) {
                        if (count++ == cut) {
                            if (process_exit)
                                ::_exit(73);
                            throw std::runtime_error("injected canonical cut");
                        }
                    });
                };
                if (process_exit) {
                    // No open coordinator is inherited across fork.
                    const auto child = ::fork();
                    check(child >= 0, "fork failed");
                    if (child == 0) {
                        try {
                            operation();
                        } catch (...) {
                            ::_exit(74);
                        }
                        ::_exit(75);
                    }
                    int status{};
                    check(::waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                              WEXITSTATUS(status) == 73,
                          "publication child did not reach fault cut");
                } else {
                    bool failed = false;
                    try {
                        operation();
                    } catch (const native::DurableUnknownOutcome &error) {
                        check(error.transaction() == 42, "unknown outcome lost caller identity");
                        ++unknowns;
                        failed = true;
                    } catch (const std::runtime_error &) {
                        failed = true;
                    }
                    check(failed, "exception cut not reached");
                }
                const auto status = store_status(path, 42);
                const auto committed = status.at("status").text() == "committed";
                check(status.at("root").at("sequence").text() == (committed ? "1" : "0"),
                      "recovered partial publication");
                check(load_store(path).fields == (committed ? updated.fields : original.fields),
                      "recovered history without matching biomedical files");
                const auto resolved = publish_store(path, root / "updated", frozen_request);
                check(resolved.at("receipt").at("root").at("sequence").text() == "1",
                      "exact retry after interrupted publication failed");
            }
        }
        check(unknowns > 0, "unknown-outcome path untested");
        // Preparation cuts cannot expose unsynchronized biomedical provenance.
        std::size_t preparation_cuts = 0;
        const auto prep_initial = initialize_store(root / "original", root / "prep-enumerate");
        publish_store(root / "prep-enumerate", root / "updated",
                      request(prep_initial.at("root"), 8), {}, [&](auto) { ++preparation_cuts; });
        for (std::size_t cut = 0; cut < preparation_cuts; ++cut) {
            const auto path = root / ("prep-exit-" + std::to_string(cut));
            const auto genesis = initialize_store(root / "original", path);
            const auto frozen_request = request(genesis.at("root"), 43);
            const auto child = ::fork();
            check(child >= 0, "preparation fork failed");
            if (child == 0) {
                std::size_t count = 0;
                try {
                    publish_store(path, root / "updated", frozen_request, {}, [&](auto) {
                        if (count++ == cut)
                            ::_exit(73);
                    });
                } catch (...) {
                    ::_exit(74);
                }
                ::_exit(75);
            }
            int status{};
            check(::waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                      WEXITSTATUS(status) == 73,
                  "preparation cut not reached");
            check(store_status(path, 43).at("status").text() == "absent" &&
                      load_store(path).fields == original.fields,
                  "provenance preparation published native state");
            publish_store(path, root / "updated", frozen_request);
            check(load_store(path).fields == updated.fields, "preparation retry failed");
        }
        std::size_t genesis_cuts = 0;
        initialize_store(root / "original", root / "genesis-enumerate",
                         [&](auto) { ++genesis_cuts; });
        for (std::size_t cut = 0; cut < genesis_cuts; ++cut) {
            const auto path = root / ("genesis-exit-" + std::to_string(cut));
            const auto child = ::fork();
            check(child >= 0, "genesis fork failed");
            if (child == 0) {
                std::size_t count = 0;
                try {
                    initialize_store(root / "original", path, [&](auto) {
                        if (count++ == cut)
                            ::_exit(73);
                    });
                } catch (...) {
                    ::_exit(74);
                }
                ::_exit(75);
            }
            int status{};
            check(::waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                      WEXITSTATUS(status) == 73,
                  "genesis cut not reached");
            if (fs::exists(path / "history/HEAD"))
                check(load_store(path).fields == original.fields, "recovered partial genesis");
            else
                rejects([&] { store_status(path); }, "HEAD-less genesis opened");
            rejects([&] { initialize_store(root / "original", path); },
                    "interrupted genesis overwritten");
        }
        // A historical snapshot is required even when a later publication is valid.
        fs::copy(root / "store", root / "binding-corrupt", fs::copy_options::recursive);
        const auto bound_manifest = root / "binding-corrupt/studies/0000000000000000/snapshot.json";
        write_text(bound_manifest, read_text(bound_manifest) + "\n");
        load_snapshot(bound_manifest.parent_path()); // Still a valid portable snapshot.
        check(store_status(root / "binding-corrupt").at("root").dump()==store_status(root/"store").at("root").dump(),
              "display-only whitespace changed native history");
        fs::copy(root / "store", root / "endpoint-disagrees", fs::copy_options::recursive);
        fs::copy(root / "store/studies/0000000000000002",
                 root / "endpoint-disagrees/studies/000000000000004d", fs::copy_options::recursive);
        {
            // The generic native store accepts a lawful exact field change, but
            // it cannot verify GADGET's external biomedical snapshot agreement.
            auto native_store =
                native::DurableNativeStore::open(root / "endpoint-disagrees/history");
            const auto before = native_store.snapshot().tail().pin();
            auto after = *before;
            ++after.sequence;
            for (auto &tensor : after.tensors) {
                if (tensor.interpretation != "AD.log2FC")
                    continue;
                auto values = tensor.tensor.values().dense();
                values[0] = 123;
                std::vector<rfc::numeric::SparseRationalEntry> entries;
                for (std::size_t i = 0; i < values.size(); ++i)
                    if (values[i] != 0)
                        entries.push_back({i, values[i]});
                tensor.tensor = rfc::numeric::SparseRationalTensor(
                    tensor.tensor.domain(), tensor.tensor.support(),
                    rfc::numeric::SparseRationalVector(values.size(), std::move(entries)));
            }
            const auto old_publication=NativeRecord::from_endpoint(after,0x4741444745540001ULL);
            auto data=old_publication.materialize_value();
            for(auto &[name,value]:std::get<NativeRecordValue::Object>(data.data).members)
                if(name=="transaction")value=NativeRecordValue(rfc::numeric::Integer{77});
            after=*old_publication.derive(data).endpoint();
            const auto update=native::make_temporal_declaration_update(before->declarations,after.declarations,
                rfc::composite::identity_declaration_maps(before->declarations->registry));
            const std::array frames{
                native::encode_temporal_segment(native::make_temporal_segment(*before, after,{}, {},native::RelationalCodeEncoding::Automatic,{},update))};
            native_store.append(77, native_store.root(), frames);
        }
        rejects([&] { store_status(root / "endpoint-disagrees"); },
                "lawful native endpoint disagreeing with bound biomedical snapshot admitted");
        const auto chain_path = root / "ordered-history";
        auto chain_root = initialize_store(root / "original", chain_path).at("root");
        const auto chain_reader = load_store(chain_path);
        for (std::uint64_t transaction = 1; transaction <= 32; ++transaction) {
            const auto published = publish_store(chain_path, root / "updated", request(chain_root, transaction));
            chain_root = published.at("root");
            check(chain_root.at("sequence").text() == std::to_string(transaction),
                  "ordered recovery omitted a publication");
        }
        const auto complete = load_store(chain_path);
        check(complete.fields == updated.fields && complete.field_presence == updated.field_presence &&
                  chain_reader.fields == original.fields && chain_reader.native_snapshot->endpoint()->sequence == 0,
              "ordered recovery changed a retained source or final field support");
        const auto middle = chain_path / "studies/0000000000000004/study.rfcxs";
        const auto hidden = chain_path / "studies/0000000000000004/retained.rfcxs";
        fs::rename(middle, hidden);
        rejects([&] { load_store(chain_path); }, "ordered recovery skipped an earlier source publication");
        fs::rename(hidden, middle);
        check(load_store(chain_path).fields == complete.fields,
              "restored ordered source history failed recovery");
        write_text(root / "store/studies/0000000000000000/metadata.json", "corrupt");
        rejects([&] { store_status(root / "store"); },
                "corrupt reachable historical provenance admitted");
        std::cout << "PASS: durable study identities, exact fields/support, pins, certificates, "
                     "retries, admission, lock, "
                  << cuts.size() << " exception cuts, " << cuts.size() << " append exits, "
                  << preparation_cuts << " preparation exits, " << genesis_cuts
                  << " genesis exits\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
