#include <gadget/execution.hpp>
#include <gadget/snapshot.hpp>
#include <iostream>
#include <rfc/storage/native_series_codec.hpp>
#include <stdexcept>

namespace {
void check(bool result, const char *message) {
    if (!result)
        throw std::runtime_error(message);
}
template <class F> void rejects(F &&f, const char *message) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, message);
}
// Repair transport hashes in negative fixtures so admission must also check
// canonical source channels and the biomedical coordinate map.
void rehash(const gadget::fs::path &path) {
    using gadget::Json;
    auto completion = Json::parse(gadget::read_text(path / "snapshot.json"));
    for (const auto &name : {"metadata.json", "study.rfcxs"}) {
        auto data = gadget::read_text(path / name);
        auto artifact = completion.at("artifacts").at(name);
        artifact.set("bytes", Json::integer(data.size()));
        artifact.set("sha256", Json::string("sha256:" + gadget::sha256(data)));
    }
    completion.set(
        "study_digest",
        Json::string(
            "sha256:" +
            gadget::sha256("gadget.snapshot.v1\n" +
                           completion.at("artifacts").at("study.rfcxs").at("sha256").text() + "\n" +
                           completion.at("artifacts").at("metadata.json").at("sha256").text() +
                           "\n")));
    gadget::write_text(path / "snapshot.json", completion.dump());
}
} // namespace
int main() {
    try {
        using namespace gadget;
        namespace native = rfc::storage;
        const fs::path root = GADGET_SNAPSHOT_TEST_BUILD;
        fs::remove_all(root);
        fs::create_directories(root);
        fs::copy(GADGET_TEST_FIXTURE, root / "input", fs::copy_options::recursive);
        auto source = load_pager(root / "input");
        for(const auto *context:{"AD","HD"}) {
            const auto prefix=std::string(context);
            check(source.audit.has(prefix+"_files") && source.audit.has(prefix+"_raw_rows") &&
                      source.audit.has(prefix+"_duplicate_rows"),
                  "native source audit preserves collected file and duplicate-row evidence");
        }
        import_measurements(source, root / "input/measurements.json");
        const auto original_manifest = source.manifest().dump();
        auto completion = save_snapshot(source, root / "study", true);
        check(completion.at("source_digest").text() == source.source_digest &&
                  !fs::exists(root / "study.partial"),
              "snapshot completion");
        const auto original_native = gadget::read_text(root / "study/study.rfcxs");
        const auto original_metadata = gadget::read_text(root / "study/metadata.json");
        rejects([&] { save_snapshot(source, root / "study", true); },
                "snapshot overwrote existing directory");
        check(gadget::read_text(root / "study/study.rfcxs") == original_native,
              "existing snapshot changed after rejection");
        fs::create_directory(root / "reserved.partial");
        rejects([&] { save_snapshot(source, root / "reserved", true); },
                "snapshot ignored staging collision");
        check(!fs::exists(root / "reserved"), "failed capture published a destination");
        fs::rename(root / "input", root / "original-input-unavailable");
        check(!fs::exists(source.root), "raw input still available");
        auto restored = load_snapshot(root / "study");
        check(restored.manifest().dump() == original_manifest,
              "source manifest changed during native replay");
        check(native::encode_rfcx(*restored.complex).bytes ==
                  native::encode_rfcx(*source.complex).bytes,
              "native retained source changed during replay");
        check(restored.genes == source.genes && restored.fields == source.fields &&
                  restored.field_presence == source.field_presence,
              "exact fields or missingness changed");
        for (std::size_t i = 0; i < source.genes.size(); ++i)
            check(restored.provenance(static_cast<std::uint32_t>(i)).dump() ==
                      source.provenance(static_cast<std::uint32_t>(i)).dump(),
                  "provenance changed during replay");
        for (const auto &name : {"AD", "HD"})
            for (std::uint16_t grade = 0; grade < 2; ++grade) {
                const auto &before = source.contexts.at(name)->grade(rfc::core::Grade{grade});
                const auto &after = restored.contexts.at(name)->grade(rfc::core::Grade{grade});
                for (std::size_t i = 0; i < before.size(); ++i)
                    check(before.get(i) == after.get(i), "context support changed");
            }
        const auto b = restored.gene_index.at("B");
        const auto observed_zero = restored.native_snapshot->field("AD.log2FC");
        const auto absent = restored.native_snapshot->field("AD.padj");
        check(observed_zero.tensor()->support().contains(b) &&
                  observed_zero.tensor()->values().dense()[b] == 0 &&
                  !absent.tensor()->support().contains(b),
              "native field lost supported zero versus absent");
        for (
            const auto &text :
            {R"({"operation":"typed_moment","context":"AD","realization":"support","field":"AD.log2FC","other_field":"AD.G","order":2})",
             R"({"operation":"green","context":"AD","realization":"rp","field":"AD.G","lambda":"0.5","mode":"exact"})",
             R"({"operation":"context_diff","context":"AD","other_context":"HD","realization":"rp","gene":"A"})"}) {
            const auto query = Json::parse(text);
            const auto before = execute_query(source, query),
                       after = execute_query(restored, query);
            check(before.at("result").dump() == after.at("result").dump() &&
                      before.at("exact").boolean() == after.at("exact").boolean() &&
                      before.at("source_manifest_sha256").text() ==
                          after.at("source_manifest_sha256").text(),
                  "query changed during replay");
            check(after.at("snapshot").at("study_digest").text() ==
                      completion.at("study_digest").text(),
                  "query omitted snapshot provenance");
        }
        rejects(
            [&] {
                import_measurements(restored,
                                    root / "original-input-unavailable/measurements.json");
            },
            "measurement import mutated an immutable native snapshot");
        save_snapshot(restored, root / "recaptured", true);
        check(gadget::read_text(root / "recaptured/study.rfcxs") == original_native &&
                  gadget::read_text(root / "recaptured/metadata.json") == original_metadata,
              "recapture was not byte stable");
        auto pin = [&] {
            auto owner = load_snapshot(root / "study");
            return owner.native_snapshot->field("AD.log2FC");
        }();
        check(pin.tensor()->values().dense() == source.fields.at("AD.log2FC") &&
                  pin.endpoint()->source.grade(rfc::core::Grade{1}).cell_count() == 4,
              "native pin did not survive workspace/directory destruction");
        auto carrier = [&] {
            auto owner = load_snapshot(root / "study");
            return owner.native_snapshot->source_grade(rfc::core::Grade{1});
        }();
        check(carrier.support_count() == 4, "native source carrier omitted retained coordinates");
        // A relocated snapshot uses fixed artifact names and no original-root reads.
        fs::copy(root / "study", root / "relocated", fs::copy_options::recursive);
        check(load_snapshot(root / "relocated").source_digest == source.source_digest,
              "snapshot relocation failed");
        auto bad = root / "bad";
        auto reset_bad = [&] {
            fs::remove_all(bad);
            fs::copy(root / "study", bad, fs::copy_options::recursive);
        };
        reset_bad();
        fs::remove(bad / "snapshot.json");
        check(load_snapshot(bad).source_digest==source.source_digest,"native snapshot incorrectly requires a JSON index");
        fs::remove(bad/"metadata.json");
        check(load_snapshot(bad).raw_sources.size()==source.files.size(),"native-only snapshot lost original input bytes");
        reset_bad();
        auto native_bytes = original_native;
        native_bytes.back() ^= 1;
        gadget::write_text(bad / "study.rfcxs", native_bytes);
        rejects([&] { load_snapshot(bad); }, "native corruption admitted");
        rehash(bad);
        rejects([&] { load_snapshot(bad); }, "native codec validation bypassed by transport hash");
        reset_bad();
        gadget::write_text(bad / "metadata.json", original_metadata + "junk");
        rehash(bad);
        rejects([&] { load_snapshot(bad); }, "metadata trailing content admitted");
        reset_bad();
        auto metadata = Json::parse(original_metadata);
        metadata.at("assertions")
            .at(0)
            .at("participants")
            .at(0)
            .set("occurrence", Json::integer(99999));
        gadget::write_text(bad / "metadata.json", metadata.dump());
        rehash(bad);
        rejects([&] { load_snapshot(bad); }, "foreign occurrence coordinate map admitted");
        reset_bad();
        metadata = Json::parse(original_metadata);
        metadata.at("fields").set("AD.log2FC", metadata.at("fields").at("AD.G"));
        gadget::write_text(bad / "metadata.json", metadata.dump());
        rehash(bad);
        rejects([&] { load_snapshot(bad); }, "aliased native field admitted");
        reset_bad();
        metadata = Json::parse(original_metadata);
        metadata.at("source_manifest").at("audit").set("occurrences", Json::integer(1));
        gadget::write_text(bad / "metadata.json", metadata.dump());
        rehash(bad);
        rejects([&] { load_snapshot(bad); }, "false native source audit admitted");
        reset_bad();
        auto manifest = Json::parse(gadget::read_text(bad / "snapshot.json"));
        manifest.at("artifacts").set("../../escape", Json::object());
        gadget::write_text(bad / "snapshot.json", manifest.dump());
        rejects([&] { load_snapshot(bad); }, "artifact path escape admitted");
        reset_bad();
        fs::remove(bad / "metadata.json");
        fs::create_symlink(root / "study/metadata.json", bad / "metadata.json");
        rejects([&] { load_snapshot(bad); }, "symlink artifact admitted");
        reset_bad();
        gadget::write_text(bad / "metadata.json",std::string(16384,'x'));
        rejects([&] { load_snapshot(bad); }, "unrelated display bytes admitted as native source description");
        fs::remove_all(bad);
        analyze(restored, root / "analysis");
        auto analysis = Json::parse(gadget::read_text(root / "analysis/manifests/analysis.json"));
        check(analysis.at("snapshot").at("study_digest").text() ==
                  completion.at("study_digest").text(),
              "analysis omitted snapshot provenance");
        std::cout << "GADGET native snapshot, exact replay, provenance, lifetime and rejection "
                     "contracts passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
