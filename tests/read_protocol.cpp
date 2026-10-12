#include "genomic_detail.hpp"
#include <gadget/record_view.hpp>
#include <iostream>

using namespace gadget;
using namespace gadget::assay_detail;
namespace {
std::size_t checks = 0;
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
template<class F> void refuse(F action) {
    bool rejected = false;
    try { action(); } catch (const std::exception &) { rejected = true; }
    check(rejected, "invalid read protocol admitted");
}
bool same_values(const NativeRecord &a, const NativeRecord &b) {
    if (a.kind() != b.kind()) return false;
    switch (a.kind()) {
    case NativeRecordKind::Null: return true;
    case NativeRecordKind::Boolean: return a.boolean() == b.boolean();
    case NativeRecordKind::Integer: return a.integer() == b.integer();
    case NativeRecordKind::Rational: return a.rational() == b.rational();
    case NativeRecordKind::Text: return a.text() == b.text();
    case NativeRecordKind::Object:
        if (a.names() != b.names()) return false;
        for (const auto &name : a.names())
            if (!same_values(a.at(name), b.at(name))) return false;
        return true;
    case NativeRecordKind::Array:
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (!same_values(a.at(i), b.at(i))) return false;
        return true;
    default: throw std::runtime_error("unexpected value kind in literal read protocol");
    }
}
Json query(const std::string &reads, const std::string &whitelist) {
    auto result = Json::parse(R"({"schema":"gadget.assays.query.v1","operation":"read_protocol_evidence","resource":"reads","resource_sha256":"","protocol":{"id":"authored","evidence":"Complete authored fixed length substitution oracle","count_unit":"read_record_occurrence","segments":[{"id":"cell","role":"barcode","start":0,"end":3,"strand":"+"},{"id":"observed_umi","role":"umi","start":3,"end":5,"strand":"+"}],"barcode":{"segment":"cell","whitelist_resource":"whitelist","whitelist_sha256":"","objective":"minimum_substitutions","maximum_substitutions":1,"symbol_rule":"literal_iupac_tokens"}}})");
    result.set("resource_sha256", Json::string("sha256:" + sha256(reads)));
    result.at("protocol").at("barcode").set("whitelist_sha256",
                                            Json::string("sha256:" + sha256(whitelist)));
    return result;
}
ResourceImages images(const std::string &reads, const std::string &whitelist) {
    return {{"reads", {Json::parse(R"({"format":"fastq"})"), reads, {}}},
            {"whitelist", {Json::parse(R"({"format":"fasta"})"), whitelist, {}}}};
}
std::string reverse(const std::string &value) {
    std::string result;
    for (auto i = value.rbegin(); i != value.rend(); ++i)
        result += *i == 'A' ? 'T' : *i == 'C' ? 'G' : *i == 'G' ? 'C' : *i == 'T' ? 'A' : 'N';
    return result;
}
}
int main(int argc, char **argv) {
    try {
        require(argc == 2, "read protocol test requires an owned fixture directory");
        const fs::path fixtures(argv[1]);
        fs::create_directories(fixtures);
        const auto root = fixtures / ("run-" + std::to_string(new_native_record_namespace()));
        require(fs::create_directory(root), "read protocol fixture creation failed");
        std::vector<std::string> reads, whitelist;
        for (char a : std::string("ACGTN"))
            for (char b : std::string("ACGTN"))
                for (char c : std::string("ACGTN")) reads.push_back(std::string{a, b, c});
        for (char a : std::string("ACGT"))
            for (char b : std::string("ACGT"))
                for (char c : std::string("ACGT")) whitelist.push_back(std::string{a, b, c});
        whitelist.push_back("AAA"); // Distinct original whitelist occurrences remain candidates.
        std::string read_bytes, whitelist_bytes;
        for (std::size_t i = 0; i < reads.size(); ++i)
            read_bytes += "@r" + std::to_string(i) + " retained header\r\n" + reads[i] +
                          "nN\r\n+\r\n!5I~!\r\n";
        read_bytes += "@short\nAC\n+\n!!\n@partial\nACGN\n+\n!!!!\n";
        for (std::size_t i = 0; i < whitelist.size(); ++i)
            whitelist_bytes += ">w" + std::to_string(i) + " original dictionary\n" +
                               whitelist[i].substr(0, 1) + "\n" + whitelist[i].substr(1) + "\n";
        const auto resources = images(read_bytes, whitelist_bytes);
        for (const auto &strand : {"+", "-"})
            for (std::size_t budget = 0; budget <= 3; ++budget) {
                auto request = query(read_bytes, whitelist_bytes);
                request.at("protocol").at("segments").at(0).set("strand", Json::string(strand));
                request.at("protocol").at("barcode").set("maximum_substitutions",
                                                          cardinality(budget));
                const auto record = NativeRecord::create(read_protocol_fields(resources, request),
                                                         new_native_record_namespace());
                check(record.at("record_count").integer() == 127 && record.at("records").size() == 127,
                      "source record count changed");
                for (std::size_t r = 0; r < reads.size(); ++r) {
                    const auto row = record.at("records").at(r);
                    const auto tokens = std::string(strand) == "+" ? reads[r] : reverse(reads[r]);
                    std::vector<std::size_t> costs;
                    std::size_t minimum = 3;
                    for (const auto &word : whitelist) {
                        std::size_t cost = 0;
                        for (std::size_t p = 0; p < 3; ++p) cost += word[p] != tokens[p];
                        costs.push_back(cost); minimum = std::min(minimum, cost);
                    }
                    check(row.at("record_index").integer() == r &&
                              row.at("read_id").text() == "r" + std::to_string(r) &&
                              row.at("minimum_substitutions").integer() == minimum,
                          "original occurrence or independent substitution minimum differs");
                    check(row.at("segments").at(0).at("tokens").text() == tokens &&
                              row.at("segments").at(1).at("tokens").text() == "NN",
                          "oriented segment or literal UMI evidence differs");
                    std::size_t count = 0;
                    for (std::size_t w = 0; w < whitelist.size(); ++w)
                        if (costs[w] == minimum && minimum <= budget) {
                            const auto candidate = row.at("barcode_candidates").at(count++);
                            check(candidate.at("whitelist_record_index").integer() == w &&
                                      candidate.at("id").text() == "w" + std::to_string(w) &&
                                      candidate.at("tokens").text() == whitelist[w] &&
                                      candidate.at("substitutions").integer() == minimum,
                                  "minimum candidate family or original dictionary order differs");
                            const auto start = candidate.at("source_start").integer().get_ui(),
                                       end = candidate.at("source_end").integer().get_ui();
                            check(whitelist_bytes.substr(start, end - start).starts_with(
                                      ">w" + std::to_string(w) + " "),
                                  "candidate byte source span differs");
                        }
                    check(row.at("barcode_candidates").size() == count &&
                              row.at("barcode_status").text() == (count == 0 ? "outside_substitution_budget" :
                                  count == 1 ? "unique_minimum" : "ambiguous_minimum"),
                          "complete ambiguity or explicit budget refusal differs");
                    const auto start = row.at("source_start").integer().get_ui(),
                               end = row.at("source_end").integer().get_ui();
                    check(read_bytes.substr(start, end - start).starts_with("@r" + std::to_string(r) + " "),
                          "original read source span differs");
                }
                const auto short_read = record.at("records").at(125),
                           partial_read = record.at("records").at(126);
                check(short_read.at("segments").at(0).at("tokens").kind() == NativeRecordKind::Null &&
                          short_read.at("minimum_substitutions").kind() == NativeRecordKind::Null &&
                          short_read.at("barcode_candidates").size() == 0 &&
                          short_read.at("barcode_status").text() == "segment_outside_original_read",
                      "short original read became a supported barcode");
                check(partial_read.at("segments").at(0).at("tokens").text().size() == 3 &&
                          partial_read.at("segments").at(1).at("tokens").kind() == NativeRecordKind::Null &&
                          partial_read.at("status").text() == "segment_outside_original_read",
                      "partial segment absence discarded existing barcode evidence");
                check(NativeRecord::decode(record.encode()).encode() == record.encode(),
                      "native protocol roundtrip differs");
            }
        const auto original_request = query(read_bytes, whitelist_bytes);
        for (int variant = 0; variant < 12; ++variant) {
            auto bad = Json::parse(original_request.dump());
            auto protocol = bad.at("protocol");
            switch (variant) {
            case 0: bad.set("resource_sha256", Json::string("sha256:wrong")); break;
            case 1: protocol.set("count_unit", Json::string("deduplicated_molecule")); break;
            case 2: protocol.at("segments").at(0).set("end", cardinality(0)); break;
            case 3: protocol.at("segments").at(1).set("id", Json::string("cell")); break;
            case 4: protocol.at("segments").at(0).set("role", Json::string("true_cell")); break;
            case 5: protocol.at("segments").at(0).set("strand", Json::string("both")); break;
            case 6: protocol.at("barcode").set("segment", Json::string("observed_umi")); break;
            case 7: protocol.at("barcode").set("whitelist_sha256", Json::string("sha256:wrong")); break;
            case 8: protocol.at("barcode").set("objective", Json::string("edit_distance")); break;
            case 9: protocol.at("barcode").set("symbol_rule", Json::string("possible_iupac")); break;
            case 10: protocol.at("barcode").set("maximum_substitutions", Json::integer(-1)); break;
            case 11: protocol.set("hidden_limit", cardinality(5)); break;
            }
            refuse([&] { (void)read_protocol_fields(resources, bad); });
        }
        auto wrong_length = resources;
        wrong_length.at("whitelist").bytes = ">w\nAA\n";
        auto bad_length = query(read_bytes, wrong_length.at("whitelist").bytes);
        refuse([&] { (void)read_protocol_fields(wrong_length, bad_length); });

        {
            const auto long_read = "@long\n" + std::string(128, 'A') + "CNN\n+\n" +
                                   std::string(131, '!') + "\n";
            std::string dictionary;
            for (std::size_t i = 0; i < 130; ++i)
                dictionary += ">d" + std::to_string(i) + "\n" + std::string(129, 'A') + "\n";
            dictionary += ">other\n" + std::string(129, 'C') + "\n";
            auto request = query(long_read, dictionary);
            request.at("protocol").at("segments").at(0).set("end", cardinality(129));
            request.at("protocol").at("segments").at(1).set("start", cardinality(129));
            request.at("protocol").at("segments").at(1).set("end", cardinality(131));
            const auto fields = NativeRecord::create(
                read_protocol_fields(images(long_read, dictionary), request),
                new_native_record_namespace());
            const auto row = fields.at("records").at(0);
            check(row.at("barcode_candidates").size() == 130 &&
                      row.at("minimum_substitutions").integer() == 1 &&
                      row.at("segments").at(1).at("tokens").text() == "NN",
                  "complete long barcode candidates or observed UMI changed");
            for (std::size_t i = 0; i < 130; ++i)
                check(row.at("barcode_candidates").at(i).at("whitelist_record_index").integer() == i &&
                          row.at("barcode_candidates").at(i).at("substitutions").integer() == 1,
                      "packed word seams or complete dictionary occurrences changed");
        }

        auto spec = Json::parse(R"({"schema":"gadget.assays.import.v1","id":"protocol","design":{"question":"Retain source bound protocol evidence","population":"Authored fixtures","selection":"All records","limitations":["No biological molecule inference"]},"assays":[],"resources":[],"links":[]})");
        for (const auto &[id, image] : resources) {
            const auto name = id + (id == "reads" ? ".fastq" : ".fa");
            write_text(root / name, image.bytes);
            auto resource = Json::parse(R"({"id":"","format":"","source":{"path":"","sha256":""},"method":"Authored complete fixture"})");
            resource.set("id", Json::string(id)); resource.set("format", image.declaration.at("format"));
            resource.at("source").set("path", Json::string(name));
            resource.at("source").set("sha256", Json::string("sha256:" + sha256(image.bytes)));
            spec.at("resources").push(resource);
        }
        write_text(root / "import.json", spec.dump());
        import_assays(root / "import.json", root / "workspace");
        auto workspace = std::make_unique<AssayWorkspace>(root / "workspace/workspace.rfcxs");
        const auto request = admit_record_view(original_request);
        const auto result = workspace->execute_record(request);
        const auto before = workspace->native_publication().encode();
        check(result.at("result").at("result").at("schema").text() == "gadget.read.protocol.evidence.v1",
              "workspace did not retain the native protocol declaration");
        check(render_record_view(result.at("result")).dump() == workspace->execute(original_request).dump(),
              "optional display differs from native result");
        const auto archived_request = result.at("request_publication").bytes();
        const auto request_bytes = request.encode();
        check(std::vector<std::byte>(archived_request.begin(), archived_request.end()) == request_bytes,
              "original native request bytes changed");
        check(result.source_namespace() == workspace->native_publication().source_namespace(),
              "actual workspace owner namespace changed");
        result.write(root / "result.rfcxs");
        workspace.reset();
        fs::rename(root / "reads.fastq", root / "reads-away.fastq");
        fs::rename(root / "whitelist.fa", root / "whitelist-away.fa");
        AssayWorkspace portable(root / "workspace/workspace.rfcxs");
        const auto replayed = portable.execute_record(request);
        check(portable.native_publication().encode() == before &&
                  replayed.source_namespace() == result.source_namespace() &&
                  same_values(replayed.at("result"), result.at("result")) &&
                  NativeRecord::read(root / "result.rfcxs").encode() == result.encode(),
              "complete portable protocol replay or native archive changed");
        const auto replayed_request = replayed.at("request_publication").bytes();
        check(replayed_request.size() == request_bytes.size() &&
                  std::equal(replayed_request.begin(), replayed_request.end(), request_bytes.begin()),
              "portable replay changed the original native request owner");
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks
                  << ",\"authored_read_occurrences\":127,\"whitelist_occurrences\":65}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
