#include <gadget/facts.hpp>
#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <iostream>

using namespace gadget;
namespace {
std::size_t checks{};
void check(bool condition) {
    if (!condition)
        throw std::runtime_error("fact API qualification failed");
    ++checks;
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected);
}
Json comparison(const Json &path, const std::string &type, const Json &value) {
    auto out = Json::object();
    out.set("path", path);
    out.set("value_type", Json::string(type));
    out.set("comparison", Json::string("equal"));
    out.set("value", value);
    return out;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 5)
            throw std::invalid_argument("manifest request claims artifact-root");
        const AssayWorkspace workspace(argv[1]);
        auto request = Json::parse(read_text(argv[2]));
        ScientificFacts facts(workspace, read_record(argv[2]).native(), argv[4]);
        const auto original_request=read_text(argv[2]);
        bool retained_request=false;
        for(const auto &source:facts.native_publication().raw_sources()) {
            const auto bytes=source.bytes();
            if(bytes.size()==original_request.size() &&
                std::equal(bytes.begin(),bytes.end(),
                    reinterpret_cast<const std::byte *>(original_request.data())))
                retained_request=true;
        }
        check(retained_request);
        const auto before = facts.catalog().dump();
        request.set("question", Json::string("Changed after admission"));
        auto view = facts.catalog();
        view.at("facts").at(0).set("value", Json::integer(999));
        check(facts.catalog().dump() == before);
        const auto reading_before = facts.model_view().dump();
        auto reading = facts.model_view();
        check(reading.at("schema").text() == "gadget.scientific.facts.model.view.v2");
        check(reading.at("facts").size() == facts.catalog().at("facts").size());
        check(reading.at("fact_count").integer() ==
              static_cast<std::int64_t>(reading.at("facts").size()));
        std::map<std::string, Json> owners;
        for (std::size_t i = 0; i < reading.at("owners").size(); ++i) {
            const auto owner = reading.at("owners").at(i);
            check(owners.emplace(owner.at("owner_sha256").text(), owner).second);
        }
        const auto captured = facts.catalog();
        for (std::size_t i = 0; i < reading.at("facts").size(); ++i) {
            const auto shown = reading.at("facts").at(i), original = captured.at("facts").at(i);
            for (const auto *name :
                 {"id", "meaning", "value_type", "value", "from", "path", "fact_sha256"})
                check(shown.at(name).dump() == original.at(name).dump());
            const auto digest = "sha256:" + sha256(original.at("owner").dump());
            check(shown.at("owner_sha256").text() == digest);
            auto restored = Json::parse(owners.at(digest).at("owner").dump());
            check(owners.at(digest).at("kind").dump() == original.at("from").at("kind").dump());
            if (restored.has("reference")) {
                auto reference = restored.at("reference");
                const auto query = original.at("owner").at("reference").at("query");
                check(reference.at("query_sha256").text() == "sha256:" + sha256(query.dump()));
                check(!reference.has("query"));
                json_object_object_del(reference.raw(), "query_sha256");
                reference.set("query", query);
                // Compare individual members: query has a different insertion
                // position after reconstructing the complete owner.
                const auto original_reference = original.at("owner").at("reference");
                check(json_object_object_length(reference.raw()) ==
                      json_object_object_length(original_reference.raw()));
                json_object_object_foreach(original_reference.raw(), key, value) {
                    check(reference.at(key).dump() == Json(json_object_get(value)).dump());
                }
                json_object_object_del(restored.raw(), "reference");
                auto expected = Json::parse(original.at("owner").dump());
                json_object_object_del(expected.raw(), "reference");
                check(restored.dump() == expected.dump());
            } else
                check(restored.dump() == original.at("owner").dump());
        }
        reading.at("facts").at(0).set("value", Json::integer(999));
        reading.at("owners").at(0).at("owner").set("target", Json::string("Changed"));
        check(facts.model_view().dump() == reading_before);
        check(facts.catalog().dump() == before);
        auto proposal = Json::parse(read_text(argv[3]));
        const auto response = facts.check(proposal);
        check(response.at("status").text() == "matched");
        proposal.at("claims").at(0).set("value", Json::integer(999));
        check(response.at("certificate").at("proposal").at("claims").at(0).at("value").integer() ==
              2);
        check(facts.check(proposal).at("status").text() == "contradicted");

        // The owner captures result bytes; reopening must reject changed files.
        const auto path = fs::path(argv[4]) / "relocated-pipeline/quality/result.rfcxs";
        const auto bytes = read_text(path);
        struct Restore {
            fs::path path;
            std::string bytes;
            ~Restore() { write_text(path, bytes); }
        } restore{path, bytes};
        write_text(path, bytes + "\n");
        rejects(
            [&] { ScientificFacts changed(workspace, Json::parse(read_text(argv[2])), argv[4]); });
        check(facts.check(Json::parse(read_text(argv[3]))).at("status").text() == "matched");
        check(facts.model_view().dump() == reading_before);
        write_text(path, bytes);

        // Shared existing certificate comparisons have no inherited source limits.
        auto certificate = Json::object(), result = Json::object(), values = Json::array();
        for (std::size_t i = 0; i < 50002; ++i)
            values.push(Json::integer(i));
        result.set("many", values);
        auto path_large = Json::array();
        path_large.push(Json::string("many"));
        path_large.push(Json::integer(50001));
        auto nested = Json::object(), cursor = nested, deep = Json::array();
        for (std::size_t i = 0; i < 24; ++i) {
            auto child = Json::object();
            cursor.set("part", child);
            cursor = child;
            deep.push(Json::string("part"));
        }
        cursor.set("value", Json::boolean(false));
        deep.push(Json::string("value"));
        result.set("deep", nested);
        auto full_deep = Json::array();
        full_deep.push(Json::string("deep"));
        for (std::size_t i = 0; i < deep.size(); ++i)
            full_deep.push(deep.at(i));
        const std::string key(5000, 'k'), number(5001, '9');
        result.set(key, Json::string(number));
        auto long_path = Json::array();
        long_path.push(Json::string(key));
        result.set("literal", Json::string("RNA α"));
        result.set("signed", Json::integer(-7));
        auto signed_path = Json::array();
        signed_path.push(Json::string("signed"));
        auto literal_path = Json::array();
        literal_path.push(Json::string("literal"));
        certificate.set("schema", Json::string("gadget.scientific.certificate.v1"));
        certificate.set("source_digest", Json::string(workspace.digest()));
        certificate.set("certificate_sha256", Json::string("sha256:" + std::string(64, '0')));
        certificate.set("result", result);
        auto claims = Json::array();
        for (std::size_t i = 0; i < 301; ++i)
            claims.push(comparison(path_large, "count", Json::integer(50001)));
        claims.push(comparison(full_deep, "boolean", Json::boolean(false)));
        claims.push(comparison(long_path, "rational", Json::string(number)));
        claims.push(comparison(literal_path, "text", Json::string("RNA α")));
        claims.push(comparison(signed_path, "integer", Json::integer(-7)));
        const auto checked = check_certificate_claims(certificate, claims);
        check(checked.at("matched").integer() == 305);
        check(checked.at("contradicted").integer() == 0);
        auto model = Json::object();
        model.set("schema", Json::string("gadget.model.claims.v1"));
        model.set("source_digest", Json::string(workspace.digest()));
        model.set("certificate_sha256", certificate.at("certificate_sha256"));
        model.set("items", claims);
        check(admit_certificate_claims(model, workspace.digest()).at("claims").at("items").size() ==
              305);
        claims.at(0).at("path").at(1); // Bounds below are checked by the shared resolver.
        auto invalid = Json::array(), beyond = Json::array();
        beyond.push(Json::string("many"));
        beyond.push(Json::integer(50002));
        invalid.push(comparison(beyond, "count", Json::integer(50002)));
        rejects([&] { check_certificate_claims(certificate, invalid); });
        invalid = Json::array();
        invalid.push(comparison(full_deep, "count", Json::integer(0)));
        rejects([&] { check_certificate_claims(certificate, invalid); });
        invalid = Json::array();
        invalid.push(comparison(path_large, "count", Json::integer(-1)));
        rejects([&] { check_certificate_claims(certificate, invalid); });
        std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
