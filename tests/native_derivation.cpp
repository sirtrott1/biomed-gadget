#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <iostream>
#include <map>
#include <rfc/storage/file_io.hpp>
#include <rfc/storage/rfcx.hpp>
#include <set>

using namespace gadget;
namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
Json file(const std::string &name, const std::string &bytes) {
    auto record = Json::object();
    record.set("path", Json::string(name));
    record.set("sha256", Json::string("sha256:" + sha256(bytes)));
    return record;
}
Json read_request(const ScientificStudy &study, const std::string &section, std::size_t offset,
                  std::size_t limit) {
    auto q = Json::parse(
        R"({"schema":"gadget.scientific.source.read.v1","publication":"self","stream":{"kind":"catalog"},"offset":0,"limit":1})");
    q.set("source_digest", Json::string(study.digest()));
    q.set("source_sequence", Json::string(std::to_string(study.endpoint()->sequence)));
    q.at("stream").set("section", Json::string(section));
    q.set("offset", Json::integer(offset));
    q.set("limit", Json::integer(limit));
    return q;
}
void refuse(const std::function<void()> &action) {
    bool refused = false;
    try {
        action();
    } catch (const std::exception &) {
        refused = true;
    }
    check(refused, "invalid native evidence or metric was accepted");
}
NativeRecordValue without_native(const ScientificStudy &study) {
    auto value = study.native_publication().at("study").unbound_value();
    auto &members = std::get<NativeRecordValue::Object>(value.data).members;
    std::erase_if(members, [](const auto &member) { return member.first == "native"; });
    return value;
}
void check_artifact_owners(const NativeRecord &workspace) {
    const auto artifacts = workspace.at("originals");
    std::set<std::string> paths;
    std::map<std::string, RawByteSource> owners;
    std::size_t aliases = 0;
    for (std::size_t i = 0; i < artifacts.size(); ++i) {
        const auto artifact = artifacts.at(i);
        check(paths.insert(artifact.at("path").text()).second,
              "original artifact path was duplicated");
        const auto value = artifact.at("bytes");
        const auto bytes = value.bytes();
        const std::string contents(reinterpret_cast<const char *>(bytes.data()), bytes.size());
        const auto owner = value.byte_source();
        const auto [previous, inserted] = owners.emplace(contents, owner);
        if (!inserted) {
            check(previous->second.shares_owner(owner) && previous->second.key() == owner.key(),
                  "equal original artifacts have separate native byte owners");
            ++aliases;
        }
        for (const auto &[other, source] : owners)
            if (other != contents)
                check(!source.shares_owner(owner) && source.key() != owner.key(),
                      "different original artifacts share a byte identity");
    }
    check(aliases == 3, "shared reference and provenance paths were lost");
}
} // namespace
int main(int argc, char **argv) {
    const fs::path root = argc > 1 ? argv[1] : "build/native-derivation-contracts";
    fs::remove_all(root);
    fs::create_directories(root);
    const std::string control = "chr\t0\t8\tbase\t1\t+\nchr\t8\t12\tlast\t2\t+\n",
                      treated =
                          "chr\t0\t4\tfirst\t3\t+\nchr\t4\t8\tzero\t0\t+\nchr\t8\t12\tlast\t2\t+\n";
    write_text(root / "control.bed", control);
    write_text(root / "treated.bed", treated);
    const std::string fasta = ">chr\nACGTACGTACGT\n";
    write_text(root / "reference.fa", fasta);
    write_text(root / "shared.txt", "Original protocol A\n");
    write_text(root / "distinct.txt", "Original protocol B\n");
    auto design = Json::parse(
        R"({"question":"Compare retained interval observations","population":"fixture","selection":"all original records","limitations":["synthetic fixture"]})");
    auto reference = Json::parse(
        R"({"assembly":"fixture","sequence_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000","coordinate_convention":"zero_based_half_open","contigs":[{"id":"chr","length":12}]})");
    reference.set("sequence_sha256", Json::string("sha256:" + sha256(fasta)));
    auto resources = Json::array();
    for (const auto &[id, bytes] :
         std::map<std::string, std::string>{{"control", control}, {"treated", treated}}) {
        auto resource = Json::object();
        resource.set("id", Json::string(id));
        resource.set("source", file(id + ".bed", bytes));
        resource.set("format", Json::string("bed"));
        resource.set("method", Json::string("Fixture original intervals"));
        resource.set("weight_column", Json::integer(4));
        resource.set("reference", reference);
        resource.set("reference_source", file("reference.fa", fasta));
        auto provenance = Json::array();
        provenance.push(file("shared.txt", read_text(root / "shared.txt")));
        provenance.push(file("distinct.txt", read_text(root / "distinct.txt")));
        resource.set("provenance", provenance);
        resources.push(resource);
    }
    auto spec = Json::object();
    spec.set("schema", Json::string("gadget.assays.import.v1"));
    spec.set("id", Json::string("intervals"));
    spec.set("design", design);
    spec.set("assays", Json::array());
    spec.set("resources", resources);
    spec.set("links", Json::array());
    write_text(root / "import.json", spec.dump());
    import_assays(root / "import.json", root / "input");
    AssayWorkspace original(root / "input/workspace.rfcxs");
    check_artifact_owners(original.native_publication());
    check_artifact_owners(NativeRecord::decode(original.native_publication().encode()));
    const auto recipe = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"controlled_coverage","samples":[{"id":"base","resource":"control","scale":"1","scale_provenance":"fixture"},{"id":"treated","resource":"treated","scale":"1","scale_provenance":"fixture"}],"contrasts":[{"id":"change","terms":[{"sample":"treated","coefficient":"1"},{"sample":"base","coefficient":"-1"}],"evidence":"fixture"}],"controls":["base"],"intercept":true,"metric":"reference_length","strand_policy":"separate","thresholds":["1"],"absence_policy":"measured_zero_within_union","units":"count","evidence":"fixture","source_kernel":true})");
    auto request = Json::object();
    request.set("schema", Json::string("gadget.assays.materialize.v1"));
    request.set("id", Json::string("coverage"));
    request.set("source_namespace", Json::string("94000"));
    request.set("design", design);
    request.set("channel", Json::string("observation"));
    request.set("recipe", recipe);
    write_text(root / "derive.json", request.dump());
    original.materialize_record(read_record(root / "derive.json").native(), root / "derived");
    check(!fs::exists(root / "derived/workspace.json"),
          "native derivation generated a JSON authority");
    AssayWorkspace derived(root / "derived/workspace.rfcxs");
    auto study = derived.study("coverage");
    check(derived.native_publication().at("manifest").at("design").at("question").text() ==
                  design.at("question").text() &&
              study->declaration().at("design").at("question").text() ==
                  design.at("question").text(),
          "request text was bound to a foreign byte source");
    check(study->native_history()->source_namespace() == 94000, "derived namespace was relabeled");
    const auto counts = study->source_index().at("catalog_sections");
    check(counts.at("cells").integer() == 8 && counts.at("occurrences").integer() == 6,
          "complete native source was not retained");
    check(study->declaration().at("cell_bindings").size() == 2 &&
              !study->declaration().has("cells") && !study->declaration().has("occurrences"),
          "compact study expanded its evidence");
    check(study->declaration().at("contexts").at(0).at("metrics").at(0).at("kind").text() ==
              "diagonal_field",
          "metric field was copied into an array");
    const auto generating = derived.derivation_certificate("coverage");
    check(generating.at("result").at("schema").text() == "gadget.coverage.fields.reference.v1" &&
              !generating.at("result").has("resources") &&
              generating.at("result").at("source_index").at("source_digest").text() ==
                  study->digest(),
          "generating answer expanded or detached the native source");
    const auto cells = study->read_source(read_request(*study, "cells", 3, 2));
    check(cells.at("records").size() == 2, "source page did not honor requested cell range");
    check(cells.at("records").at(0).at("id").text() == "interval:1",
          "primary source order changed");
    const auto evidence = cells.at("records").at(0).at("evidence").at(0);
    check(evidence.at("start").integer() == 0 &&
              evidence.at("end").integer() == static_cast<std::int64_t>(control.find('\n')),
          "original interval byte span changed");
    const auto occurrences = study->read_source(read_request(*study, "occurrences", 4, 2));
    check(occurrences.at("records").size() == 2 && occurrences.at("complete").boolean(),
          "late occurrence page failed");
    check(occurrences.at("records").at(0).at("key").text() == "5",
          "original occurrence order changed");
    ScientificSession session(study);
    const auto summary = session.execute(Json::parse(
        R"({"operation":"field_summary","field":"all.contrast:0:residual","support":"reported_only"})"));
    check(summary.has("result"), "native study numerical query failed");
    check(
        session
            .execute(Json::parse(R"({"operation":"hodge","context":"all","field":"all.sample:0"})"))
            .has("result"),
        "field-bound metric was not usable by Hodge");
    refuse([&] {
        session.execute(Json::parse(
            R"({"operation":"field_summary","field":"all.interval.sample","support":"reported_only"})"));
    });
    const auto hodge =
        session
            .execute(Json::parse(R"({"operation":"hodge","context":"all","field":"all.sample:0"})"))
            .at("result");
    check(hodge.at(0).at("value").at("rational").text() == "12" &&
              hodge.at(2).at("value").at("rational").text() == "16",
          "native metric field did not enter the Hodge operator");
    const auto green =
        session
            .execute(Json::parse(
                R"({"operation":"green_action","context":"all","field":"all.sample:0","lambda":"1"})"))
            .at("result");
    check(green.at(0).at("value").at("rational").text() == "1/13" &&
              green.at(2).at("value").at("rational").text() == "2/9",
          "native metric field did not enter the Green solve");
    const auto producer = study->producer_publication();
    const auto owner = producer.raw_sources().front();
    check(owner.shares_owner(owner.rekey(owner.key() + 100)) &&
              !owner.shares_owner(RawByteSource(owner.key(), owner.bytes())),
          "equal source keys or contents were treated as shared byte ownership");
    for (const auto &source : producer.raw_sources()) {
        const auto all = study->native_publication().raw_sources();
        const auto found = std::find_if(
            all.begin(), all.end(), [&](const auto &item) { return item.key() == source.key(); });
        check(found != all.end() && found->bytes().data() == source.bytes().data(),
              "unchanged historical byte owners were duplicated");
    }

    for (const auto *replacement : {"interval.source_start", "sample:0"}) {
        auto invalid = without_native(*study);
        auto &items = std::get<NativeRecordValue::Object>(invalid.data).members;
        if (std::string(replacement) == "sample:0") {
            const auto bad =
                study->field_accession("all.contrast:0:residual").tensor_declaration()->key;
            for (auto &[name, value] : items)
                if (name == "contexts") {
                    auto &context = std::get<NativeRecordValue::Object>(
                        std::get<NativeRecordValue::Array>(value.data).elements[0].data);
                    for (auto &[field, metric] : context.members)
                        if (field == "metrics") {
                            auto &m = std::get<NativeRecordValue::Object>(
                                std::get<NativeRecordValue::Array>(metric.data).elements[0].data);
                            for (auto &[key, v] : m.members)
                                if (key == "key")
                                    v = NativeRecordValue(std::to_string(bad));
                        }
                }
        } else {
            const auto bad =
                study->field_accession("all." + std::string(replacement)).tensor_declaration()->key;
            for (auto &[name, value] : items)
                if (name == "cell_bindings") {
                    auto &binding = std::get<NativeRecordValue::Object>(
                        std::get<NativeRecordValue::Array>(value.data).elements[1].data);
                    for (auto &[field, evidence] : binding.members)
                        if (field == "evidence")
                            for (auto &[key, v] :
                                 std::get<NativeRecordValue::Object>(evidence.data).members)
                                if (key == "end_field")
                                    v = NativeRecordValue(std::to_string(bad));
                }
        }
        refuse([&] {
            ScientificStudy invalid_study(ScientificStudy::publish_annotations(
                producer, std::move(invalid),
                study->native_publication().at("documents").unbound_value()));
        });
    }
    fs::rename(root / "input", root / "input-away");
    fs::remove(root / "control.bed");
    fs::remove(root / "treated.bed");
    fs::remove(root / "derive.json");
    check(derived.replay_derivation("coverage", root / "replayed").at("status").text() == "PASS",
          "native-only derivation replay failed");
    auto bind = Json::parse(
        R"({"schema":"gadget.assays.bindings.v1","id":"bound","assays":[{"input":"derived","assay":"coverage","id":"alias"}],"resources":[],"links":[],"method":"Preserve actual study owner","evidence":"Fixture"})");
    bind.set("design", design);
    auto retained = std::make_shared<const AssayWorkspace>(root / "derived/workspace.rfcxs");
    AssayWorkspace::bind_inputs({{"derived", retained}}, bind, root / "bound");
    AssayWorkspace bound(root / "bound/workspace.rfcxs");
    check(bound.study("alias")->digest() == study->digest(),
          "native collection binding changed the original study");
    check(bound.derivation_certificate("alias").at("certificate_sha256").text() ==
              generating.at("certificate_sha256").text(),
          "collection alias changed the original generating answer");
    check(bound.replay_derivation("alias", root / "alias-replayed").at("status").text() == "PASS",
          "bound study replay lost its original generating request");
    for (const auto &[format, bytes] : std::map<std::string, std::string>{
             {"bed", control},
             {"gtf", "chr\tfixture\texon\t1\t4\t.\t+\t.\tgene_id "
                     "\"a\";\nchr\tfixture\texon\t5\t8\t.\t-\t.\tgene_id \"b\";\n"},
             {"gff3", "##gff-version "
                      "3\nchr\tfixture\tgene\t1\t4\t.\t+\t.\tID=a\nchr\tfixture\tgene\t5\t8\t.\t-"
                      "\t.\tID=b\n"},
             {"fasta", ">a desc\r\nAC\r\nGT\r\n>b\nA\n"},
             {"fastq", "@a\nACGT\n+a\n!!!!\n@b\nAA\n+\n##\n"}}) {
        write_text(root / (format + ".input"), bytes);
        auto input = Json::object();
        input.set("schema", Json::string("gadget.bioinformatics.file.v1"));
        input.set("source", file(format + ".input", bytes));
        input.set("format", Json::string(format));
        if (format != "fasta" && format != "fastq") {
            input.set("reference", reference);
            if (format == "bed")
                input.set("weight_column", Json::integer(4));
        }
        write_text(root / (format + ".json"), input.dump());
        const auto native = inspect_bioinformatics_record(root / (format + ".json"));
        native.write(root / (format + ".rfcxs"));
        const auto restored = NativeRecord::read(root / (format + ".rfcxs"));
        check(!restored.accession_directory()->cell({{0}, 1}).declared_values().empty(),
              "original byte source has no native document owner");
        const auto original_bytes = restored.at("result").at("original").bytes();
        check(original_bytes.size() == bytes.size() &&
                  std::equal(original_bytes.begin(), original_bytes.end(),
                             reinterpret_cast<const std::byte *>(bytes.data())),
              "native physical file lost original bytes");
        check(restored.endpoint()->source.grade({1}).cell_count() == 2 &&
                  restored.endpoint()->source.grade({1}).incidence_count() == 2,
              "native physical records lost source containment");
        const auto view = inspect_bioinformatics_file(root / (format + ".json"));
        check(view.at("record_count").integer() == 2, "physical compatibility parser changed");
        if (format == "fasta" || format == "fastq")
            check(view.at("records").at(0).at("length").integer() == 4,
                  "sequence parser length changed");
    }
    std::cout << "Native evidence, metric fields, complete source, derivation, relocation and "
                 "replay PASS\n";
}
