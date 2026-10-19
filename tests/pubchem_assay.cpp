#include "../examples/source-formats/pubchem.hpp"
#include <iostream>
#include <stdexcept>

using namespace gadget;
using namespace gadget::examples;
namespace {
constexpr PubChemSourceLimits limits{16 * 1024 * 1024, 4096, 100000, 1000000, 100000, 300000};
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F callback) {
    bool refused{};
    try { callback(); } catch (const std::exception &) { refused = true; }
    check(refused, "invalid PubChem source was accepted");
}
const std::string authored = R"({"PC_AssaySubmit":{"assay":{"descr":{"aid":{"id":100,"version":1},"name":"authored assay", "protocol":["retain\\noriginal"],"results":[{"tid":1,"name":"IC50 Modifier","type":4,"unit":254,"ac_qualifier":true},{"tid":2,"name":"IC50","type":1,"unit":5,"ac":true,"constraints":{"fmin":0}},{"tid":3,"name":"integer","type":2},{"tid":4,"name":"Boolean","type":3},{"tid":5,"name":"text","type":4,"sunit":"arbitrary"}]}},"data":[{"sid":100,"version":0,"outcome":1,"data":[{"tid":1,"value":{"sval":"<"}},{"tid":2,"value":{"fval":-0.000e+00}},{"tid":3,"value":{"ival":-0}},{"tid":4,"value":{"bval":false}},{"tid":5,"value":{"sval":""}}]},{"sid":101,"data":[]},{"sid":100,"outcome":2,"data":[{"tid":1,"value":{"sval":"unknown modifier"}},{"tid":2,"value":{"fval":1e309}}]}],"revoke":[101,101]}})";
Json changed(const std::function<void(Json)> &mutation) {
    auto raw = Json::parse(authored); mutation(raw); return raw;
}
void contracts() {
    const auto parsed = inspect_pubchem_assay(authored, "100", limits);
    check(parsed.bytes == authored, "original source bytes changed");
    check(parsed.original.dump_compact() == Json::parse(authored).dump_compact(), "original document changed");
    const auto projection = parsed.projection, rows = projection.at("substances");
    check(projection.at("description").at("original").at("protocol").at(0).text() == "retain\\noriginal", "protocol dropped");
    check(projection.at("definitions").size() == 5, "non-numeric definitions dropped");
    check(rows.size() == 3 && rows.at(0).at("sid").text() == rows.at(2).at("sid").text(), "repeated substance occurrence coalesced");
    check(rows.at(0).at("outcome").integer() == 1 && rows.at(2).at("outcome").integer() == 2, "reported outcomes changed or filtered");
    check(!rows.at(1).at("outcome_reported").boolean() && !rows.at(1).has("outcome"), "omitted outcome became a reported default");
    check(rows.at(0).at("readouts").at(1).at("number_token").text() == "-0.000e+00", "floating zero spelling lost");
    check(rows.at(0).at("readouts").at(2).at("number_token").text() == "-0", "integer negative zero spelling lost");
    check(!rows.at(0).at("readouts").at(3).at("original").at("value").at("bval").boolean(), "reported false treated as absent");
    check(rows.at(0).at("readouts").at(4).at("support").text() == "reported", "empty reported string treated as absent");
    check(rows.at(2).at("readouts").at(1).at("number_token").text() == "1e309", "large finite decimal rounded to infinity");
    check(rows.at(2).at("readouts").at(0).at("original").at("value").at("sval").text() == "unknown modifier", "unknown modifier changed");
    check(!rows.at(0).at("readouts").at(1).has("comparison"), "modifier paired with endpoint by name");
    for (std::size_t i = 0; i < 5; ++i)
        check(rows.at(1).at("readouts").at(i).at("support").text() == "absent", "missing readout became zero");
    check(projection.at("revocations").size() == 2 && projection.at("revocations").at(0).at("sid").text() == "101", "reported exclusions changed or coalesced");
    auto empty = changed([](auto raw) {
        json_object_object_del(raw.at("PC_AssaySubmit").raw(), "data");
        json_object_object_del(raw.at("PC_AssaySubmit").at("assay").at("descr").raw(), "results");
    });
    check(inspect_pubchem_assay(empty.dump_compact(), "100", limits).projection.at("substances").size() == 0, "metadata-only source acquired rows");
    const auto refuse = [&](const std::function<void(Json)> &mutation) {
        rejects([&] { (void)inspect_pubchem_assay(changed(mutation).dump_compact(), "100", limits); });
    };
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("assay").at("descr").at("results").push(raw.at("PC_AssaySubmit").at("assay").at("descr").at("results").at(0)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("assay").at("descr").at("results").at(0).set("type", Json::integer(5)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").push(raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(0)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(0).set("tid", Json::integer(90)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(0).set("value", Json::object()); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(0).at("value").set("fval", Json::integer(1)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(1).at("value").set("fval", Json::string("0")); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(3).at("value").set("bval", Json::integer(0)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(2).at("value").set("ival", Json::number(1.5)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(4).at("value").set("sval", Json{}); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).at("data").at(1).set("value", Json{}); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("assay").at("descr").at("aid").set("version", Json::integer(-1)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).set("outcome", Json::integer(9)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).set("version", Json::integer(-1)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).set("sid", Json::integer(0)); });
    refuse([](auto raw) { raw.at("PC_AssaySubmit").at("data").at(0).set("sid", Json::string("100")); });
    for (const auto &token : {"9223372036854775808", "-9223372036854775809", "18446744073709551616"}) {
        auto bytes = authored; const auto pos = bytes.find("\"ival\":-0"); bytes.replace(pos, 9, "\"ival\":" + std::string(token));
        rejects([&] { (void)inspect_pubchem_assay(bytes, "100", limits); });
    }
    for (const auto &bytes : {"{}", "{\"PC_AssaySubmit\":[]}", "{\"PC_AssaySubmit\":{\"assay\":{\"aid\":100}}}", "{} trailing"})
        rejects([&] { (void)inspect_pubchem_assay(bytes, "100", limits); });
    rejects([&] { (void)inspect_pubchem_assay(authored, "101", limits); });
    rejects([&] { (void)inspect_pubchem_assay(authored, "0100", limits); });
    auto duplicate = authored; duplicate.insert(duplicate.find("\"sid\":100"), "\"sid\":101,");
    rejects([&] { (void)inspect_pubchem_assay(duplicate, "100", limits); });
    for (const auto &bounded : {PubChemSourceLimits{authored.size()-1,4096,100000,1000000,100000,300000},
         PubChemSourceLimits{authored.size(),4,100000,1000000,100000,300000},
         PubChemSourceLimits{authored.size(),4096,2,1000000,100000,300000},
         PubChemSourceLimits{authored.size(),4096,3,14,100000,300000}})
        rejects([&] { (void)inspect_pubchem_assay(authored, "100", bounded); });
    std::stop_source cancelled; cancelled.request_stop();
    rejects([&] { (void)inspect_pubchem_assay(authored, "100", limits, cancelled.get_token()); });

    const std::string csv = "AID,SID,CID,Target GeneID,Target Accession,Activity Value [uM],Other\r\n100,100,100,1,,999,\"a\r\nb\"\r\n100,100,200,,NP_1.1,1,\r\n100,100,200,,,0,\r\n100,101,,,,,\r\n";
    const auto concise = inspect_pubchem_concise(csv, "100", limits);
    const auto links = concise.projection.at("links");
    check(concise.bytes == csv && concise.original.at("rows").size() == 4, "concise source occurrence lost");
    check(links.size() == 5, "concise repeated/ambiguous links coalesced or missing link fabricated");
    check(links.at(0).at("left").text() == "pubchem:sid:100" && links.at(0).at("right").text() == "pubchem:cid:100", "SID/CID namespaces merged");
    check(links.at(2).at("right").text() == "pubchem:cid:200" && links.at(4).at("right").text() == "pubchem:cid:200", "conflicting or repeated link lost");
    check(links.at(3).at("right").text() == "pubchem:target_accession:NP_1.1" && links.at(3).at("relationship").text() == "reported_assay_target_accession", "target accession given an undeclared NCBI/protein identity");
    check(links.at(0).at("cells").at(6).text() == "a\r\nb", "quoted CSV source cell changed");
    check(!links.at(0).has("value") && !links.at(0).has("number_token"), "concise activity became full assay readout");
    check(inspect_pubchem_concise("AID,SID,CID\n100,100,\n", "100", limits).projection.at("links").size() == 0, "missing CID became SID");
    for (const auto &bad : {"AID,SID\n100,100\n", "AID,SID,CID\n101,100,100\n", "AID,SID,CID\n100,1x,100\n", "AID,SID,CID\n100,100,0\n", "AID,SID,CID\n100,100\n", "AID,SID,CID,CID\n100,1,2,3\n"})
        rejects([&] { (void)inspect_pubchem_concise(bad, "100", limits); });
    auto bounded = limits; bounded.concise_rows = 3;
    rejects([&] { (void)inspect_pubchem_concise(csv, "100", bounded); });
    bounded = limits; bounded.links = 4;
    rejects([&] { (void)inspect_pubchem_concise(csv, "100", bounded); });
    rejects([&] { (void)inspect_pubchem_concise(csv, "100", limits, cancelled.get_token()); });
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 4) {
            const auto source = read_text(argv[2], limits.source_bytes);
            const auto inspected = std::string(argv[1]) == "assay"
                ? inspect_pubchem_assay(source, argv[3], limits)
                : std::string(argv[1]) == "concise" ? inspect_pubchem_concise(source, argv[3], limits)
                : throw std::invalid_argument("mode must be assay or concise");
            auto output = Json::object();
            output.set("source_bytes", Json::string(inspected.bytes));
            output.set("source_sha256", Json::string("sha256:" + sha256(inspected.bytes)));
            output.set("original", inspected.original); output.set("projection", inspected.projection);
            std::cout << output.dump_compact() << '\n'; return 0;
        }
        if (argc != 1) throw std::invalid_argument("usage: pubchem_assay [assay|concise source expected-aid]");
        contracts(); std::cout << checks << " PubChem source-format checks passed\n";
        return 0;
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
