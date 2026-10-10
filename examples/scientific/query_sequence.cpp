#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <iostream>
#include <set>
#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

namespace {
using namespace gadget;
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
std::vector<std::pair<std::string, Json>> expand(const ScientificStudy &source, const NativeRecord &request) {
    std::vector<std::pair<std::string, Json>> queries;
    if (request.has("queries")) {
        require(!request.has("context_template"), "choose queries or a context template");
        const auto authored = request.at("queries");
        for (std::size_t i = 0; i < authored.size(); ++i)
            queries.emplace_back(authored.at(i).at("id").text(), render_record_view(authored.at(i).at("query")));
        return queries;
    }
    const auto specification = request.at("context_template");
    require(specification.at("restriction").text() == "scalar_identity" &&
                specification.at("metric").text() == "identity", "context template requires declared scalar identity maps and metric");
    const auto fields = specification.at("fields"), methods = specification.at("methods");
    require(fields.size() > 0 && methods.size() > 0, "context template requires fields and methods");
    std::vector<std::string> contexts;
    const auto declared = source.native_publication().at("study").at("contexts");
    if (specification.has("contexts")) {
        for (std::size_t i = 0; i < specification.at("contexts").size(); ++i)
            contexts.push_back(specification.at("contexts").at(i).text());
    } else {
        for (std::size_t i = 0; i < declared.size(); ++i)
            if (declared.at(i).has("witness_scope")) contexts.push_back(declared.at(i).at("id").text());
    }
    require(!contexts.empty(), "context template requires declared block contexts");
    for (const auto &context : contexts) {
        std::optional<std::uint64_t> key;
        for (std::size_t i = 0; i < declared.size(); ++i)
            if (declared.at(i).at("id").text() == context && declared.at(i).has("witness_scope"))
                key = std::stoull(declared.at(i).at("block").text());
        require(key.has_value(), "template context must use its declared block carrier");
        const auto endpoint = source.endpoint();
        const auto block = std::find_if(endpoint->local_blocks.begin(), endpoint->local_blocks.end(),
                                       [&](const auto &b) { return b.block_key == *key; });
        require(block != endpoint->local_blocks.end(), "template block missing");
        for (std::size_t m = 0; m < methods.size(); ++m) {
            const auto method = methods.at(m).text();
            auto q = Json::object(); q.set("operation", Json::string(method)); q.set("context", Json::string(context));
            if (method == "local_sections") {
                auto stalks = Json::array(), maps = Json::array(), sections = Json::array();
                for (const auto &cell : block->cells) {
                    auto item = Json::object(), axes = Json::array(); axes.push(Json::string("scalar"));
                    item.set("grade", Json::integer(cell.grade.value)); item.set("key", Json::string(std::to_string(cell.key)));
                    item.set("axes", axes); stalks.push(item);
                }
                for (const auto &occurrence : block->occurrences) {
                    auto item = Json::object(); item.set("grade", Json::integer(occurrence.upper_grade.value));
                    item.set("key", Json::string(std::to_string(occurrence.key.value)));
                    item.set("matrix", Json::parse(R"([["1"]])")); maps.push(item);
                }
                for (std::size_t f = 0; f < fields.size(); ++f) {
                    const auto name = fields.at(f).text();
                    const auto meaning = source.field_declaration(name);
                    const auto pin = source.field_accession(name);
                    require(pin.carrier()->grade.value == 0, "template sections require grade-zero source fields");
                    auto item = Json::object(), names = Json::array(); names.push(Json::string(name));
                    item.set("id", Json::string(name)); item.set("grade", Json::integer(0)); item.set("fields", names);
                    auto axes = Json::array(); axes.push(Json::string("scalar")); item.set("axes", axes);
                    item.set("quantity", meaning.at("quantity")); item.set("units", meaning.at("units")); sections.push(item);
                }
                q.set("construction", Json::string("cellular_sheaf")); q.set("metric", Json::string("identity"));
                q.set("stalks", stalks); q.set("maps", maps); q.set("sections", sections);
                queries.emplace_back(context + ":" + method, q);
            } else {
                require(method == "green_action" || method == "hodge" || method == "adjoint", "unsupported context template method");
                if (method == "green_action") q.set("lambda", Json::string(specification.at("green_lambda").text()));
                for (std::size_t f = 0; f < fields.size(); ++f) {
                    q.set("field", Json::string(fields.at(f).text()));
                    queries.emplace_back(context + ":" + method + ":" + fields.at(f).text(), Json::parse(q.dump()));
                }
            }
        }
    }
    return queries;
}
} // namespace
int main(int argc, char **argv) try {
    using namespace gadget;
    using N = NativeRecordValue;
    require(argc == 4, "usage: gadget_scientific_sequence WORKSPACE.rfcxs REQUEST.rfcxs OUTPUT.rfcxs");
    require(!fs::exists(argv[3]), "scientific sequence output already exists");
    const auto request = NativeRecord::read(argv[2]);
    require(request.at("schema").text() == "gadget.scientific.sequence.v1", "scientific sequence schema");
    AssayWorkspace workspace(argv[1]);
    const auto source = workspace.study(request.at("assay").text());
    const auto queries = expand(*source, request);
    require(!queries.empty(), "scientific sequence requires queries");
    ScientificSession session(source);
    std::set<std::string> identities;
    for (const auto &[id, query] : queries) {
        require(identities.insert(id).second, "duplicate sequence query identity");
        session.validate(query);
    }
    N::Array answers;
    for (const auto &[id, query] : queries) {
        const auto result = session.execute(query);
        answers.elements.emplace_back(N::Object{{{"id", N(id)},
                                                {"certificate", import_record_view(result)}}});
    }
    N result(N::Object{{{"schema", N(std::string("gadget.scientific.sequence.result.v1"))},
        {"source_digest", N(source->digest())}, {"answers", N(std::move(answers))}}});
    const auto publication = retain_operation(source->native_publication(), request, std::move(result));
    fs::create_directories(fs::absolute(argv[3]).parent_path());
    publication.write(argv[3]);
    std::cout << "PASS\nqueries\t" << queries.size() << '\n';
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
