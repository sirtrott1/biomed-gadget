#include "assay_detail.hpp"
#include <gadget/pipeline.hpp>
#include <atomic>
#include <barrier>
#include <csignal>
#include <exception>
#include <iostream>
#include <sys/resource.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>

namespace {
void io_boundary(const gadget::fs::path &root) {
    using namespace gadget;
    fs::create_directories(root);
    const auto check = [](bool ok, const char *message) {
        if (!ok) throw std::logic_error(message);
    };
    const auto rejects = [&](auto action) {
        bool refused = false;
        try { action(); } catch (const std::exception &) { refused = true; }
        check(refused, "invalid file operation admitted");
    };
    const std::string original("original\0bytes\n", 15);
    const auto source = root / "source";
    write_text(source, original);
    check(read_text(source) == original && sha256_file(source) == sha256(original),
          "all original bytes or complete file hash changed");
    const fs::path nul(source.string() + std::string(1, '\0') + "suffix");
    rejects([&] { (void)read_text(nul); });
    rejects([&] { (void)sha256_file(nul); });
    rejects([&] { write_text(nul, "replacement"); });
    rejects([&] { (void)read_text(fs::path{}); });
    rejects([&] { write_text(fs::path{}, "replacement"); });
    rejects([&] { (void)read_text(root); });
    rejects([&] { (void)read_text(source, original.size() - 1); });
    check(read_text(source, original.size()) == original, "exact caller byte budget refused");
    check(read_text(source) == original, "path refusal modified its lexical prefix");

    const auto link = root / "symbolic";
    fs::create_symlink(source, link);
    rejects([&] { write_text(link, "replacement"); });
    check(fs::is_symlink(link) && read_text(source) == original,
          "symbolic output changed original bytes");
    const auto dangling = root / "dangling";
    fs::create_symlink(root / "absent", dangling);
    rejects([&] { write_text(dangling, "replacement"); });
    check(fs::is_symlink(dangling) && !fs::exists(root / "absent"),
          "dangling output link followed");
    const auto hard = root / "hard";
    fs::create_hard_link(source, hard);
    write_text(hard, "independent replacement");
    check(read_text(source) == original && read_text(hard) == "independent replacement" &&
              !fs::equivalent(source, hard), "replacement mutated another linked original");
    write_text(root / "empty", "");
    check(read_text(root / "empty", 0).empty(), "empty file or zero-byte budget refused");
    const auto previous = fs::current_path();
    fs::current_path(root);
    try {
        write_text("relative", original);
        check(read_text("relative") == original, "relative publication changed original bytes");
    } catch (...) {
        fs::current_path(previous);
        throw;
    }
    fs::current_path(previous);

    for (const auto &destination : {source, root / "new"}) {
        const auto pid = fork();
        if (pid < 0) throw std::runtime_error("cannot start file exhaustion consumer");
        if (pid == 0) {
            std::signal(SIGXFSZ, SIG_IGN);
            const rlimit limit{64, 64};
            if (setrlimit(RLIMIT_FSIZE, &limit)) _exit(2);
            try { write_text(destination, std::string(1u << 20, 'x')); }
            catch (const std::exception &) { _exit(0); }
            _exit(3);
        }
        int status{};
        check(waitpid(pid, &status, 0) == pid && WIFEXITED(status) && WEXITSTATUS(status) == 0,
              "file byte exhaustion did not refuse");
        check(read_text(source) == original && !fs::exists(root / "new"),
              "failed write exposed partial or replaced original bytes");
        for (const auto &entry : fs::directory_iterator(root))
            check(!entry.path().filename().string().starts_with(".rfc-stage-"),
                  "failed staged write left residue");
    }
    write_text(root / "new", original);
    check(read_text(root / "new") == original, "write failure prevented exact retry");

    const auto output = root / "concurrent";
    const std::string first(33027, 'a'), second(571, 'b');
    write_text(output, first);
    std::atomic<bool> done{};
    std::atomic<std::size_t> reads{};
    std::exception_ptr failure;
    std::barrier ready(2);
    std::jthread reader([&] {
        ready.arrive_and_wait();
        try {
            do {
                const auto bytes = read_text(output);
                check(bytes == first || bytes == second, "reader saw partial or mixed replacement");
                ++reads;
            } while (!done.load());
        } catch (...) { failure = std::current_exception(); }
    });
    ready.arrive_and_wait();
    try {
        for (std::size_t i = 0; i < 128; ++i)
            write_text(output, i % 2 ? first : second);
    } catch (...) {
        done = true;
        reader.join();
        throw;
    }
    done = true;
    reader.join();
    if (failure) std::rethrow_exception(failure);
    check(reads > 0 && read_text(output) == first,
          "complete replacement reader/final state not exercised");
}
} // namespace

int main(int argc, char **argv) {
    try {
        using namespace gadget;
        using namespace gadget::assay_detail;
        if (argc != 2)
            throw std::invalid_argument("expected fixture directory");
        const auto root = fs::absolute(argv[1]);
        fs::remove_all(root);
        fs::create_directories(root);
        io_boundary(root / "io-boundary");
        const auto check = [](bool ok, const char *message) {
            if (!ok) throw std::logic_error(message);
        };
        auto spec = Json::parse(R"({"format":"matrix_market","arithmetic":"exact_rational","quantity":"expression","implicit_entries":"unobserved","missing_tokens":[]})");
        spec.set("features", file_record("features.tsv", "id\na\nb\nc\n"));
        spec.set("observations", file_record("observations.tsv", "id\nx\ny\n"));
        write_text(root / "features.tsv", "id\na\nb\nc\n");
        write_text(root / "observations.tsv", "id\nx\ny\n");
        const auto read = [&](const std::string &bytes) {
            write_text(root / "matrix.txt", bytes);
            spec.set("source", file_record("matrix.txt", bytes));
            return read_matrix(root, spec);
        };
        const auto refuse = [&](const std::string &bytes) {
            bool rejected = false;
            try { (void)read(bytes); } catch (const std::invalid_argument &) { rejected = true; }
            check(rejected, "invalid matrix admitted");
        };
        const std::string banner = "%%MatrixMarket matrix coordinate real general\n";
        const std::string source = banner + "% original row order\r\n3 2 4\r\n3\t2\t7/3\r\n1 1 0\r\n3 1 -2\r\n2 2 123456789012345678901234567890";
        const auto matrix = read(source);
        check(matrix.bytes == source && matrix.features.ids[2] == "c", "source bytes or axes changed");
        check(matrix.columns[0].state(0) == Rational{0} &&
                  !matrix.columns[0].state(1) && matrix.columns[0].state(2) == Rational{-2},
              "zero and absence conflated");
        check(matrix.columns[1].state(2) == Rational{7, 3} &&
                  matrix.columns[1].state(1) == Rational{"123456789012345678901234567890"},
              "exact coefficients changed");
        for (const auto &rows : {"1 1 0\n1 1 0\n", "1 1 2\n1 1 -2\n", "1 1 0\n1 1 2\n"})
            refuse(banner + "3 2 2\n" + rows);
        for (const auto &body : {"3 2 1\n", "3 2 0\n1 1 2\n", "3 2 1\n0 1 2\n",
                                 "3 2 1\n4 1 2\n", "3 2 1\n1 3 2\n", "3 2 1\n1 1 2 extra\n",
                                 "3 2 1\n1 1\n", "3 2 0\n\n", "3 2 0 extra\n",
                                 "2 2 0\n", "3 2 1\n18446744073709551616 1 2\n"})
            refuse(banner + body);
        refuse(banner + "3 2 1\n1 1 " + std::string(1, '\0') + "\n");
        refuse("%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 1/2\n");
        for (const auto *quantity : {"read_count", "umi_count", "fragment_count"}) {
            spec.set("quantity", Json::string(quantity));
            refuse(banner + "3 2 1\n1 1 -2\n");
            refuse(banner + "3 2 1\n1 1 1/2\n");
            check(read(banner + "3 2 1\n1 1 2\n").columns[0].state(0) == Rational{2},
                  "integer count measurement changed");
        }
        spec.set("quantity", Json::string("expression"));
        const auto mask = [&](const std::string &bytes) {
            write_text(root / "missing.tsv", bytes);
            spec.set("missing", file_record("missing.tsv", bytes));
        };
        mask("2\t1\n");
        check(!read(source).columns[0].state(1), "explicit missing support changed");
        mask("1\t1\n"); refuse(source); // Reported zero cannot be missing.
        mask("2\t1\n2\t1\n"); refuse(source);
        spec = Json::parse(spec.dump()); json_object_object_del(spec.raw(), "missing");
        spec.set("implicit_entries", Json::string("measured_zero"));
        check(read(source).columns[0].support().count() == 3, "full zero support lost");
        mask("2\t1\n");
        check(read(source).columns[0].support().count() == 2, "full support mask lost");
        json_object_object_del(spec.raw(), "missing");
        spec.set("format", Json::string("dense_tsv"));
        spec.set("missing_tokens", Json::parse("[\"NA\"]"));
        const std::string dense = "id\tx\ty\r\na\t0\t7/3\r\nb\tNA\t-2\r\nc\t5\t0";
        const auto d = read(dense);
        check(d.bytes == dense && !d.columns[0].state(1) && d.columns[1].state(2) == Rational{0},
              "dense missing/zero support changed");
        refuse(dense + "\na\t0\t0\n");
        refuse("id\tx\ty\na\t0\t0\nb\t0\t0\n");
        refuse("id\ty\tx\na\t0\t0\nb\t0\t0\nc\t0\t0\n");
        spec.set("arithmetic", Json::string("exact_decimal"));
        check(read("id\tx\ty\na\t1.25\t0\nb\t0\t0\nc\t0\t0\n").columns[0].state(0) == Rational{5, 4},
              "decimal conversion changed");
        // Complete private columns can transfer ownership. A fully reported
        // sparse support still follows the original full support normalization.
        const std::string complete = banner + "3 2 6\n3 2 0\n1 1 0\n3 1 -2\n2 2 7/3\n2 1 5\n1 2 9\n";
        write_text(root / "matrix.txt", complete);
        auto import = Json::parse(R"({"schema":"gadget.assays.import.v1","id":"intake","design":{"question":"Verify carrier-preserving intake","population":"authored fields","selection":"all","limitations":["fixture"]},"assays":[],"resources":[],"links":[]})");
        auto assay = Json::parse(R"({"id":"x","kind":"expression","format":"matrix_market","source_namespace":"93001","carrier":"features","quantity":"expression","units":"dimensionless","channel":"E","method":"authored exact values","arithmetic":"exact_rational","implicit_entries":"unobserved","missing_tokens":[],"partition_size":3})");
        assay.set("source", file_record("matrix.txt", complete));
        for (const auto *name : {"features", "observations"}) assay.set(name, spec.at(name));
        import.at("assays").push(assay);
        write_text(root / "import.json", import.dump());
        std::vector<std::string> stages;
        const auto imported = import_assay_workspace(root / "import.json", root / "whole",
            [&](std::string_view stage) { stages.emplace_back(stage); });
        check(stages.front() == "specification" && stages.back() == "workspace_validation",
              "import stage boundaries missing");
        check(imported->study("x")->digest() ==
                  AssayWorkspace(root / "whole/workspace.rfcxs").study("x")->digest(),
              "returned import owner differs from independent native replay");
        check(imported->native_publication().encode() ==
                  rfc::storage::read_file_exact(root / "whole/workspace.rfcxs",
                                                native_data_limits().max_archive_bytes),
              "returned native publication differs from its complete output bytes");
        bool aborted = false;
        try {
            import_assays(root / "import.json", root / "aborted", [](std::string_view stage) {
                if (stage == "matrix_read") throw std::runtime_error("observer aborted import");
            });
        } catch (const std::runtime_error &) { aborted = true; }
        check(aborted && !fs::exists(root / "aborted"),
              "observer failure retained an unpublished workspace");
        admit_record_view(import).write(root / "import.rfcxs");
        const auto native_imported = import_assay_workspace(root / "import.rfcxs", root / "native-import");
        check(rfc::storage::read_file_exact(root / "native-import/assays/0/partitions/0/source.rfcxs",
                                       native_data_limits().max_archive_bytes) ==
                  rfc::storage::read_file_exact(root / "whole/assays/0/partitions/0/source.rfcxs",
                                           native_data_limits().max_archive_bytes),
              "native import declaration changed the complete numerical source");
        check(rfc::storage::read_file_exact(root / "native-import/import-spec.rfcxs",
                                       native_data_limits().max_archive_bytes) ==
                  rfc::storage::read_file_exact(root / "import.rfcxs", native_data_limits().max_archive_bytes) &&
                  native_imported->study("x")->endpoint()->tensors[0].tensor.state(2) == Rational{-2},
              "native declaration bytes or exact signed values changed");
        const auto publication=NativeRecord::read(root/"whole/assays/0/partitions/0/study.rfcxs");
        check(publication.at("schema").text()=="gadget.biomedical.native.v2" &&
              publication.at("history").at("kind").text()=="publication_prefix",
              "generated assay embeds a duplicate archive instead of its owning native prefix");
        const auto retained=*publication.publication_history();
        check(retained.segment_count()==1,"generated assay annotation history missing");
        const auto prefix=rfc::storage::NativeSeriesPrefix::capture(retained,retained.first_sequence());
        check(rfc::storage::encode_native_series(prefix.history())==rfc::storage::read_file_exact(
            root/"whole/assays/0/partitions/0/source.rfcxs",native_data_limits().max_archive_bytes),
            "generated source archive bytes changed");
        auto corrupt=retained;const auto old_tail=corrupt.tail();auto endpoint=old_tail.endpoint();++endpoint.sequence;
        endpoint.tensors[0].tensor={rfc::numeric::IndexDomain::linear(3),rfc::numeric::StructuralSupport::full(3),
            rfc::numeric::SparseRationalVector(3,{{0,Rational{999}}})};
        corrupt.append_target(endpoint);
        bool refused=false;
        try {ScientificStudy altered(NativeRecord::from_history(corrupt));}
        catch(const std::exception &) {refused=true;}
        check(refused,"annotation was allowed to change its bound field data");
        const auto whole_study=imported->study("x");
        const auto source_tail=whole_study->native_history()->tail();
        check(whole_study->endpoint()==source_tail.pin(),
              "numerical execution uses a copied or later annotation endpoint");
        check(whole_study->endpoint()->sequence<whole_study->native_publication().endpoint()->sequence,
              "annotation order was confused with the declared numerical source cutoff");
        const auto whole = whole_study->endpoint();
        for (const auto &tensor : whole->tensors)
            check(tensor.tensor.support().storage() == rfc::numeric::StructuralSupport::Storage::Full,
                  "reported full support was not normalized");
        check(whole->tensors[0].tensor.state(2) == Rational{-2} &&
                  whole->tensors[1].tensor.state(1) == Rational{7, 3}, "published coefficients changed");
        import.at("assays").at(0).set("partition_size", Json::integer(2));
        write_text(root / "partitioned.json", import.dump());
        import_assays(root / "partitioned.json", root / "partitioned");
        AssayWorkspace partial(root / "partitioned/workspace.json");
        const auto first = partial.study("x", 0)->endpoint(), last = partial.study("x", 1)->endpoint();
        check(first->tensors[0].tensor.state(1) == Rational{5} &&
                  last->tensors[0].tensor.state(0) == Rational{-2} &&
                  last->tensors[1].tensor.state(0) == Rational{0}, "partition slices changed");
        fs::copy_file(root / "partitioned/workspace.rfcxs", root / "matrix-only.rfcxs");
        fs::remove_all(root / "partitioned");
        auto relocated_matrix = std::make_shared<const AssayWorkspace>(root / "matrix-only.rfcxs");
        const auto matrix_binding = Json::parse(R"({"schema":"gadget.assays.bindings.v1","id":"bound-matrix","design":{"question":"retain partitioned native fields","population":"fixture","selection":"all","limitations":["fixture"]},"assays":[{"input":"source","assay":"x","id":"bound-x"}],"resources":[],"links":[],"method":"Explicit native assay alias","evidence":"Complete original native source owners"})");
        AssayWorkspace::bind_inputs({{"source", relocated_matrix}}, matrix_binding, root / "bound-matrix");
        fs::copy_file(root / "bound-matrix/workspace.rfcxs", root / "bound-matrix-only.rfcxs");
        AssayWorkspace legacy_bound(root / "bound-matrix/workspace.json");
        check(legacy_bound.study("bound-x", 1)->digest() == relocated_matrix->study("x", 1)->digest(),
              "compatibility bound partition changed its owner");
        import_assays(root / "bound-matrix/reimport.json", root / "matrix-reimport");
        fs::remove_all(root / "bound-matrix");
        AssayWorkspace bound_matrix(root / "bound-matrix-only.rfcxs");
        for (std::size_t p = 0; p < 2; ++p) {
            check(bound_matrix.study("bound-x", p)->digest() == relocated_matrix->study("x", p)->digest(),
                  "native-only partition binding changed original source identity");
            const auto &before = relocated_matrix->study("x", p)->endpoint()->tensors;
            const auto &after = bound_matrix.study("bound-x", p)->endpoint()->tensors;
            check(before.size() == after.size(), "bound field count changed");
            for (std::size_t f = 0; f < before.size(); ++f)
                for (std::size_t c = 0; c < before[f].tensor.domain().logical_size(); ++c)
                    check(before[f].tensor.state(c) == after[f].tensor.state(c),
                          "bound numerical value, zero or absence changed");
        }
        // Compressed/original resource artifacts are arbitrary bytes retained
        // by the canonical publication, including after native only relocation.
        const std::string opaque("\0\xff\x1f\x8b\0original", 13);
        write_text(root / "original.gz", opaque);
        write_text(root / "reads.fastq", "@read\nACGT\n+\nIIII\n@short\nNN\n+\n!!\n");
        auto resource = Json::object(), provenance = Json::array();
        resource.set("id", Json::string("reads"));
        resource.set("format", Json::string("fastq"));
        resource.set("method", Json::string("Original compressed bytes with decoded FASTQ view"));
        resource.set("source", file_record("reads.fastq", read_text(root / "reads.fastq")));
        provenance.push(file_record("original.gz", opaque));
        resource.set("provenance", provenance);
        import.set("assays", Json::array());
        import.at("resources").push(resource);
        write_text(root / "resource-import.json", import.dump());
        import_assays(root / "resource-import.json", root / "resources");
        const auto original_path = root / "resources/resources/0/provenance-0";
        const auto refuse_publication = [&] {
            bool rejected = false;
            try {
                AssayWorkspace candidate(root / "resources/workspace.rfcxs");
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            check(rejected, "native admission accepted a changed external artifact");
        };
        auto changed = opaque;
        changed.back() = 'x';
        write_text(original_path, changed);
        refuse_publication();
        write_text(original_path, opaque.substr(0, opaque.size() - 1));
        refuse_publication();
        fs::remove(original_path);
        fs::create_symlink(root / "original.gz", original_path);
        refuse_publication();
        fs::remove(original_path);
        write_text(original_path, opaque);
        {
            AssayWorkspace with_originals(root / "resources/workspace.rfcxs");
            with_originals.write_native(root / "verified-artifacts.rfcxs");
            check(read_text(root / "verified-artifacts.rfcxs") ==
                      read_text(root / "resources/workspace.rfcxs"),
                  "native admission changed the complete original publication");
        }
        fs::copy_file(root / "resources/workspace.rfcxs", root / "portable.rfcxs");
        fs::remove_all(root / "resources");
        fs::remove(root / "original.gz");
        fs::remove(root / "reads.fastq");
        auto portable = std::make_shared<const AssayWorkspace>(root / "portable.rfcxs");
        const auto quality_request = admit_record_view(Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_quality","resource":"reads","quality_thresholds":[0,20,40,93]})"));
        const auto quality_operation = portable->execute_record(quality_request);
        const auto quality = quality_operation.at("result").at("result");
        check(quality.at("cycle_histograms").kind() == NativeRecordKind::Tensor,
              "read-quality numerical result was imported as repeated values");
        const auto histogram_record = quality.at("cycle_histograms");
        const auto &histogram = histogram_record.tensor();
        check(histogram.domain().extents() == std::vector<std::size_t>({4, 94}) &&
                  histogram.support().count() == 376,
              "complete native quality axes or measured-zero support lost");
        for (std::size_t b = 0; b < 4; ++b)
            for (std::size_t q = 0; q < 94; ++q)
                check(histogram.state(b * 94 + q) == Rational{q == 40 || (b < 2 && q == 0) ? 1 : 0},
                      "independent native quality histogram coordinate mismatch");
        check(quality.at("lengths").tensor().state(1) == Rational{2} &&
                  quality.at("gc_lower_bases").tensor().state(1) == Rational{0} &&
                  quality.at("gc_upper_bases").tensor().state(1) == Rational{2} &&
                  quality.at("eligible_reads").tensor().state(0) == Rational{2} &&
                  quality.at("eligible_reads").tensor().state(2) == Rational{1},
              "ambiguous bases, zero support or shorter-read cycle eligibility changed");
        quality_operation.write(root / "quality.rfcxs");
        const auto restored_quality = NativeRecord::read(root / "quality.rfcxs");
        check(restored_quality.at("result").at("result").at("cycle_histograms").kind() == NativeRecordKind::Tensor,
              "native quality field became an array after persistence");
        check(render_record_view(restored_quality.at("result")).dump() ==
                  portable->execute(render_record_view(quality_request)).dump(),
              "native quality view changed the complete compatibility certificate");
        const auto quality_defaults = Json::parse(R"({"schema":"gadget.pipeline.default.v1","profile":"rna_seq","question":"Audit every stored base and eligible cycle","goals":["quality"],"reads":"reads","quality_thresholds":[0,40]})");
        const auto quality_index = execute_pipeline(portable, default_pipeline(*portable, quality_defaults),
                                                   root / "quality-pipeline");
        const auto quality_output = NativeRecord::read(root / "quality-pipeline" /
            quality_index.at("results").at(0).at("result").text());
        check(quality_output.at("result").at("result").at("cycle_histograms").kind() == NativeRecordKind::Tensor,
              "pipeline replaced the typed native quality result with imported arrays");
        bool found = false;
        const auto originals = portable->native_publication().at("originals");
        for (std::size_t i = 0; i < originals.size(); ++i) {
            const auto entry = originals.at(i);
            if (entry.at("path").text() == "resources/0/provenance-0") {
                const auto bytes = entry.at("bytes").bytes();
                found = std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()) == opaque;
            }
        }
        check(found, "native resource provenance bytes changed or were lost");
        const auto bindings = Json::parse(R"({"schema":"gadget.assays.bindings.v1","id":"bound","design":{"question":"retain original bytes","population":"fixture","selection":"all","limitations":["fixture"]},"assays":[],"resources":[{"input":"source","resource":"reads","id":"bound-reads"}],"links":[],"method":"Explicit resource alias","evidence":"Actual original source publication"})");
        AssayWorkspace::bind_inputs({{"source", portable}}, bindings, root / "bound");
        check(read_text(root / "bound/resources/0/provenance-0") == opaque,
              "bound resource lost original provenance");
        import_assays(root / "bound/reimport.json", root / "reimported");
        check(read_text(root / "reimported/resources/0/provenance-0") == opaque,
              "reimport lost original resource provenance");
        const auto derivation = Json::parse(R"({"schema":"gadget.assays.materialize.v1","id":"assembly","source_namespace":"94000","design":{"question":"retain original source through derivation","population":"fixture","selection":"all","limitations":["fixture"]},"recipe":{"schema":"gadget.assays.query.v1","operation":"sequence_assembly_graph","resource":"reads","k":2,"ambiguity":"literal_tokens","evidence":"All source read occurrences"},"channel":"E"})");
        portable->materialize(derivation, root / "derived");
        AssayWorkspace derived(root / "derived/workspace.rfcxs");
        const auto replay = derived.replay_derivation("assembly", root / "replayed");
        check(replay.at("status").text() == "PASS", "derived resource provenance replay failed");
        AssayWorkspace replayed_workspace(root / "replayed/workspace.rfcxs");
        const auto source_archive = NativeRecord::decode(replayed_workspace.study("assembly")
            ->producer_publication().at("result").at("workspace_publication").bytes());
        bool replay_bytes = false;
        for (std::size_t i = 0; i < source_archive.at("originals").size(); ++i) {
            const auto bytes = source_archive.at("originals").at(i).at("bytes").bytes();
            if (bytes.size() == opaque.size() && std::equal(bytes.begin(), bytes.end(),
                    reinterpret_cast<const std::byte *>(opaque.data()))) replay_bytes = true;
        }
        check(replay_bytes, "derived replay omitted original resource bytes");
        write_text(root / "empty.fastq", "");
        auto empty_resource = Json::object();
        empty_resource.set("id", Json::string("empty-reads"));
        empty_resource.set("format", Json::string("fastq"));
        empty_resource.set("method", Json::string("Complete empty source contract"));
        empty_resource.set("source", file_record("empty.fastq", ""));
        auto empty_spec = Json::parse(import.dump()), empty_resources = Json::array();
        empty_resources.push(empty_resource); empty_spec.set("resources", empty_resources);
        write_text(root / "empty-import.json", empty_spec.dump());
        import_assays(root / "empty-import.json", root / "empty-workspace");
        AssayWorkspace empty_workspace(root / "empty-workspace/workspace.rfcxs");
        auto empty_request = render_record_view(quality_request);
        empty_request.set("resource", Json::string("empty-reads"));
        const auto empty_result = empty_workspace.execute_record(admit_record_view(empty_request));
        const auto empty_fields = empty_result.at("result").at("result");
        const auto empty_histogram = empty_fields.at("cycle_histograms");
        check(empty_histogram.tensor().domain().extents() == std::vector<std::size_t>{0, 94} &&
                  empty_histogram.tensor().support().count() == 0 &&
                  empty_fields.at("total_bases").integer() == 0 &&
                  render_record_view(empty_result.at("result")).at("result").at("reads").size() == 0,
              "empty source acquired artificial reads, support or bases");
        std::cout << "source bytes, fields, support and portable resource provenance passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
