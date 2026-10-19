#include "pathway_ora.hpp"
#include <gadget/pathways.hpp>
#include <gadget/record_view.hpp>
#include "scientific_result_fields.hpp"
#include "pathway_fixture.hpp"
#include "feature_identity_fixture.hpp"
#include <iostream>

using namespace gadget;
using namespace gadget::pathway_detail;
namespace {
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void refuses(F action) {
    bool rejected{};
    try {
        action();
    } catch (const std::exception &) {
        rejected = true;
    }
    check(rejected, "invalid pathway admission accepted");
}
Json overrepresentation(const std::vector<std::string> &identities, const Source &source,
    const Json &query, std::stop_token stop = {}) {
    scientific_detail::ResultFields fields(NativeRecord::create(NativeRecordValue{}, new_native_record_namespace()));
    auto value = pathway_detail::overrepresentation(identities, source, query, fields, stop);
    return scientific_detail::result_fields_view(fields.finish(std::move(value)));
}
Rational exact(const Json &value) {return Rational(value.text());}
Json request() {
    return Json::parse(
        R"({"schema":"gadget.pathways.overrepresentation.v1","operation":"pathway_overrepresentation","analysis_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000","gene_sets_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000","feature_namespace":"ensembl:gene","features_path":["features"],"feature_id_path":["id"],"background_indices":[0,1,2,3,4,5],"selected_indices":[0,1,2],"mapping_policy":{"ambiguous":"exclude","unmapped":"exclude","many_to_one":"count_member_once","duplicate_occurrences":"retain_count_once"},"set_ids":["set-a","set-b","filtered"],"size_filter":{"minimum":1,"maximum":100,"basis":"members_in_background"},"membership_policy":{"duplicates":"retain_count_once","empty":"retain_exclude","outside_background":"retain_exclude"},"correction":"benjamini_hochberg","selection_method":"Explicit authored foreground occurrences","null_assumptions":["Uniform conditional gene selection is declared for this fixture"],"limits":{"source_bytes":1048576,"mapping_rows":10000,"membership_cells":100000,"sets":1000,"features":1000,"field_bytes":4096,"kernel_work":1000000,"archive_bytes":1048576,"declared_buffer_bytes":16777216}})");
}
const std::vector<std::string> ids{"ENSG1.1", "ENSG1.1_PAR_Y", "ENSG2.1",
                                   "ENSG3.1", "ENSG4.1",       "ENSG5.1"};
const std::string gmt =
    "set-a\tReleased description\tG1\tG2\tG2\tOUT\t\nset-b\tOther\tG4\nfiltered\tOutside "
    "universe\tOUT\n";
const std::string mapping =
    "feature_id\tmember\tevidence\nENSG1.1\tG1\tdeclared correspondence\nENSG1.1\tG1\trepeated "
    "evidence\nENSG1.1_PAR_Y\tG2\tdistinct Y copy\nENSG2.1\tG3\tambiguous "
    "alternative\nENSG2.1\tX\tambiguous alternative\nENSG4.1\tG4\tdeclared "
    "correspondence\nENSG5.1\tG1\tmany-to-one occurrence\nnot-an-input\tOTHER\tretained outside "
    "input\n";
Source source(std::string_view genes = gmt, std::string_view map = mapping) {
    return {"authored-gene-sets",
            "authored.v1",
            "declared:gene",
            "sha256:" + sha256(genes),
            "declared-correspondences",
            "authored.v1",
            "explicit_declared",
            "ensembl:gene",
            "sha256:" + sha256(map),
            genes,
            map};
}
void contracts() {
    const auto q = request();
    const auto s = source();
    const auto result = overrepresentation(ids, s, q);
    check(result.at("background_members").dump_compact() == "[\"G1\",\"G2\",\"G4\"]",
          "background extension or undeclared collapse");
    check(result.at("selected_members").dump_compact() == "[\"G1\",\"G2\"]",
          "ambiguous selected feature admitted");
    check(result.at("feature_admission").size() == 6 &&
              result.at("feature_admission").at(1).at("feature_id").text() == "ENSG1.1_PAR_Y",
          "source feature version/Y copy lost");
    check(result.at("feature_admission").at(2).at("status").text() == "ambiguous" &&
              result.at("feature_admission").at(3).at("status").text() == "unmapped",
          "mapping refusals/exclusions lost");
    check(result.at("mapping_occurrences").size() == 8,
          "repeated/unmatched mapping occurrences dropped");
    check(result.at("many_to_one_occurrences").size() == 1 &&
              result.at("many_to_one_occurrences").at(0).at("feature_indices").dump_compact() ==
                  "[0,5]",
          "many-to-one contributing occurrences lost");
    check(result.at("tested_family_size").integer() == 2 && result.at("excluded_sets").size() == 1,
          "tested correction family changed");
    const auto a = result.at("tests").at(0);
    check(a.at("N").integer() == 3 && a.at("K").integer() == 2 && a.at("draws").integer() == 2 &&
              a.at("hits").integer() == 2,
          "overlap counts differ from declared universe");
    check(exact(a.at("p_value")) == Rational(1,3),
          "exact inclusive finite tail");
    check(exact(a.at("adjusted_p_value")) == Rational(2,3),
          "declared family BH correction");
    check(a.at("membership_occurrences").size() == 5 &&
              a.at("membership_occurrences").at(2).at("status").text() == "repeated" &&
              a.at("membership_occurrences").at(4).at("status").text() == "empty",
          "original duplicate/empty GMT members lost");
    check(result.at("mapping_origin").text() == "explicit_declared",
          "declared mapping mislabeled released");
    check(result.at("declared_set_family").size() == 3,
          "declared/excluded family distinction lost");
    auto empty = request();
    empty.set("selected_indices", Json::array());
    const auto no_hits = overrepresentation(ids, s, empty);
    check(no_hits.at("selected_members").size() == 0 &&
              exact(no_hits.at("tests").at(0).at("p_value")) == 1,
          "empty declared selection must remain a zero-draw test");
    for (const auto method : {"benjamini_yekutieli", "bonferroni"}) {
        auto copy = Json::parse(q.dump());
        copy.set("correction", Json::string(method));
        check(overrepresentation(ids, s, copy).at("tests").size() == 2,
              "declared correction refused");
    }
    for (const auto &[key, value] :
         std::vector<std::pair<std::string, std::string>>{{"ambiguous", "refuse"},
                                                          {"unmapped", "refuse"},
                                                          {"many_to_one", "refuse"},
                                                          {"duplicate_occurrences", "drop"}}) {
        auto copy = Json::parse(q.dump());
        copy.at("mapping_policy").set(key, Json::string(value));
        refuses([&] { overrepresentation(ids, s, copy); });
    }
    for (const auto &[key, value] :
         std::vector<std::pair<std::string, std::string>>{{"analysis_sha256", "wrong"},
                                                          {"gene_sets_sha256", "wrong"},
                                                          {"feature_namespace", "gene-symbol"},
                                                          {"correction", "unknown"}}) {
        auto copy = Json::parse(q.dump());
        copy.set(key, Json::string(value));
        refuses([&] { overrepresentation(ids, s, copy); });
    }
    for (const auto &[key, value] :
         std::vector<std::pair<std::string, std::string>>{{"background_indices", "[0,0]"},
                                                          {"background_indices", "[100]"},
                                                          {"background_indices", "[]"},
                                                          {"selected_indices", "[100]"},
                                                          {"selected_indices", "[0,0]"},
                                                          {"set_ids", "[\"absent\"]"},
                                                          {"set_ids", "[\"set-a\",\"set-a\"]"},
                                                          {"set_ids", "[]"}}) {
        auto copy = Json::parse(q.dump());
        copy.set(key, Json::parse(value));
        refuses([&] { overrepresentation(ids, s, copy); });
    }
    for (const auto key :
         {"source_bytes", "mapping_rows", "membership_cells", "sets", "features", "field_bytes",
          "kernel_work", "archive_bytes", "declared_buffer_bytes"}) {
        auto copy = Json::parse(q.dump());
        copy.at("limits").set(key, Json::integer(0));
        refuses([&] { overrepresentation(ids, s, copy); });
    }
    for (const auto &[key, n] : std::vector<std::pair<std::string, int>>{{"source_bytes", 20},
                                                                         {"mapping_rows", 2},
                                                                         {"membership_cells", 4},
                                                                         {"sets", 1},
                                                                         {"features", 3},
                                                                         {"field_bytes", 5},
                                                                         {"kernel_work", 1}}) {
        auto copy = Json::parse(q.dump());
        copy.at("limits").set(key, Json::integer(n));
        refuses([&] { overrepresentation(ids, s, copy); });
    }
    auto wrong = s;
    wrong.sha256 = "sha256:" + std::string(64, '0');
    refuses([&] { overrepresentation(ids, wrong, q); });
    wrong = s;
    wrong.mapping_sha256 = "sha256:" + std::string(64, '0');
    refuses([&] { overrepresentation(ids, wrong, q); });
    wrong = s;
    wrong.mapping_origin = "inferred";
    refuses([&] { overrepresentation(ids, wrong, q); });
    wrong = s;
    wrong.release = "";
    refuses([&] { overrepresentation(ids, wrong, q); });
    for (const auto &bad :
         std::vector<std::string>{gmt + gmt, std::string("x\n"), gmt + std::string(1, '\0'), "\n"})
        refuses([&] { overrepresentation(ids, source(bad, mapping), q); });
    for (const auto &bad : std::vector<std::string>{
             "feature_id\tmember\nA\tG\n", "feature_id\tmember\tevidence\nA\tG\t\n",
             "feature_id\tmember\tevidence\n\tG\tevidence\n"})
        refuses([&] { overrepresentation(ids, source(gmt, bad), q); });
    auto changed = Json::parse(q.dump());
    changed.at("size_filter").set("basis", Json::string("all_source_members"));
    refuses([&] { overrepresentation(ids, s, changed); });
    // This strictly positive tail is smaller than the binary64 normal range.
    std::vector<std::string> extreme_ids;
    std::string extreme_mapping = "feature_id\tmember\tevidence\n";
    std::string extreme_gmt = "extreme\tdeclared underflow fixture";
    auto extreme = request();
    auto background = Json::array(), foreground = Json::array();
    for (std::size_t i = 0; i < 1200; ++i) {
        auto id = "E" + std::to_string(i), member = "G" + std::to_string(i);
        extreme_ids.push_back(id);
        extreme_mapping += id + "\t" + member + "\texplicit test correspondence\n";
        background.push(Json::integer(static_cast<std::int64_t>(i)));
        if (i < 600) {
            extreme_gmt += "\t" + member;
            foreground.push(Json::integer(static_cast<std::int64_t>(i)));
        }
    }
    extreme_gmt += "\n";
    extreme.set("background_indices", background);
    extreme.set("selected_indices", foreground);
    extreme.set("set_ids", Json::parse(R"(["extreme"])"));
    extreme.at("size_filter").set("maximum", Json::integer(1200));
    extreme.at("limits").set("features", Json::integer(1200));
    extreme.at("limits").set("declared_buffer_bytes", Json::integer(67108864));
    const auto tiny = overrepresentation(extreme_ids, source(extreme_gmt, extreme_mapping), extreme)
                          .at("tests")
                          .at(0);
    check(exact(tiny.at("p_value")) > 0 && exact(tiny.at("p_value")) < Rational("1/1" + std::string(308, '0')),
          "strictly positive exact tail lost");
    std::stop_source cancellation;
    cancellation.request_stop();
    refuses([&] { overrepresentation(ids, s, q, cancellation.get_token()); });
    std::cout << "PASS: " << checks << " pathway application admission/family/refusal checks\n";
}
void evaluate(const fs::path &path) {
    const auto input = Json::parse(read_text(path));
    const auto meta = input.at("source");
    const auto genes = input.at("gmt").text(), map = input.at("mapping").text();
    Source s{meta.at("id").text(),
             meta.at("release").text(),
             meta.at("identifier_namespace").text(),
             meta.at("sha256").text(),
             meta.at("mapping_id").text(),
             meta.at("mapping_release").text(),
             meta.at("mapping_origin").text(),
             meta.at("feature_namespace").text(),
             meta.at("mapping_sha256").text(),
             genes,
             map};
    std::vector<std::string> features;
    for (std::size_t i = 0; i < input.at("features").size(); ++i)
        features.push_back(input.at("features").at(i).text());
    std::cout << overrepresentation(features, s, input.at("request")).dump() << '\n';
}
} // namespace
int admission_main(int argc, char **argv) {
    try {
        if (argc == 1)
            contracts();
        else if (argc == 2 && std::string_view(argv[1]) == "--request")
            std::cout << request().dump() << '\n';
        else if (argc == 3 && std::string_view(argv[1]) == "--evaluate")
            evaluate(argv[2]);
        else
            throw std::invalid_argument("pathway test invocation");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

int main(int argc, char **argv) {
    if (argc != 2 || std::string_view(argv[1]).starts_with("--")) return admission_main(argc, argv);
    try {
        const char *arguments[] = {argv[0]};
        if (admission_main(1, const_cast<char **>(arguments))) return 1;
        const auto parent = fs::path(argv[1]); fs::create_directories(parent);
        const auto root = private_directory(parent,"pathway-");
        auto fixture = pathway_test::fixture();
        const auto before = pathway_test::archive_sha256(fixture.analysis);
        const auto publication = pathway_overrepresentation(fixture.analysis, fixture.gene_sets, fixture.request);
        check(pathway_test::archive_sha256(fixture.analysis) == before, "analysis owner changed");
        const auto result = scientific_detail::result_fields_view(publication.at("result"));
        check(exact(result.at("tests").at(0).at("p_value")) == Rational(1,6), "native exact tail");
        const auto field = publication.encoded_tensor(101);
        check(field.has_value(), "original field owner absent");
        const auto retained = NativeRecord::decode(publication.at("result").at("view").at("gene_set_source_archive").bytes());
        check(pathway_test::archive_sha256(retained) == pathway_test::archive_sha256(fixture.gene_sets), "released source owner changed");
        publication.write(root/"result.rfcxs");
        fixture.analysis.write(root/"analysis.rfcxs"); fixture.gene_sets.write(root/"sets.rfcxs");
        fixture.request.write(root/"request.rfcxs");
        const auto restored = NativeRecord::read(root/"result.rfcxs");
        const auto replayed = pathway_overrepresentation(NativeRecord::read(root/"analysis.rfcxs"),
            NativeRecord::read(root/"sets.rfcxs"), NativeRecord::read(root/"request.rfcxs"));
        check(result.dump() == scientific_detail::result_fields_view(restored.at("result")).dump(), "cold restored fields differ");
        check(result.dump() == scientific_detail::result_fields_view(replayed.at("result")).dump(), "native replay fields differ");
        auto wrong = render_record_view(fixture.request); wrong.set("analysis_sha256", Json::string("sha256:"+std::string(64, '0')));
        refuses([&] {pathway_overrepresentation(fixture.analysis, fixture.gene_sets, admit_record_view(wrong));});
        const auto repeated=test::repeated_feature_identity_source();
        auto bounded=render_record_view(fixture.request);
        bounded.set("analysis_sha256",Json::string(pathway_test::archive_sha256(repeated)));
        bounded.at("limits").set("declared_buffer_bytes",Json::integer(1048576));
        bool early_refusal{};
        try {(void)pathway_overrepresentation(repeated,fixture.gene_sets,admit_record_view(bounded));}
        catch(const std::exception &error) {
            early_refusal=std::string_view(error.what()).find("identity working byte admission")!=std::string_view::npos;
        }
        check(early_refusal,"shared native label spans expanded beyond working admission");
        std::cout << "PASS: original native fields, exact tail, retained source and cold replay\n";
        return 0;
    } catch (const std::exception &e) {std::cerr << e.what() << '\n'; return 1;}
}
