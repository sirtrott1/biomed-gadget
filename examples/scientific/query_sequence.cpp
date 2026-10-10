#include <gadget/assay.hpp>
#include <gadget/profile_process.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>

namespace {
using namespace gadget;
using Consumer = std::function<void(const std::string &, const Json &)>;
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
class QuerySource {
    NativeRecord request_;
    struct Context { std::string id; const rfc::storage::LocalBlockDefinition *block; };
    struct Field { std::string id; Json meaning; };
    std::vector<Context> contexts_;
    std::vector<Field> fields_;
    std::vector<std::string> methods_;
    std::string lambda_;
    std::size_t count_{};
  public:
    QuerySource(const ScientificStudy &source, NativeRecord request) : request_(std::move(request)) {
        if (request_.has("queries")) {
            require(!request_.has("context_template"), "choose queries or a context template");
            count_ = request_.at("queries").size();
            return;
        }
        const auto specification = request_.at("context_template");
        require(specification.at("restriction").text() == "scalar_identity" &&
                    specification.at("metric").text() == "identity",
                "context template requires declared scalar identity maps and metric");
        const auto fields = specification.at("fields"), methods = specification.at("methods");
        require(fields.size() > 0 && methods.size() > 0, "context template requires fields and methods");
        for (std::size_t i = 0; i < methods.size(); ++i) {
            const auto method = methods.at(i).text();
            require(method == "local_sections" || method == "green_action" || method == "hodge" || method == "adjoint",
                    "unsupported context template method");
            if (method == "green_action") lambda_ = specification.at("green_lambda").text();
            methods_.push_back(method);
        }
        for (std::size_t i = 0; i < fields.size(); ++i) {
            const auto id = fields.at(i).text();
            if (std::find(methods_.begin(), methods_.end(), "local_sections") != methods_.end())
                require(source.field_accession(id).carrier()->grade.value == 0,
                        "template sections require grade-zero source fields");
            fields_.push_back({id, source.field_declaration(id)});
        }
        const auto declared = source.native_publication().at("study").at("contexts");
        std::map<std::string, std::uint64_t> context_keys;
        std::vector<std::string> ordered;
        for (std::size_t i = 0; i < declared.size(); ++i)
            if (declared.at(i).has("witness_scope")) {
                const auto id = declared.at(i).at("id").text();
                context_keys.emplace(id, std::stoull(declared.at(i).at("block").text()));
                ordered.push_back(id);
            }
        if (specification.has("contexts")) {
            ordered.clear();
            for (std::size_t i = 0; i < specification.at("contexts").size(); ++i)
                ordered.push_back(specification.at("contexts").at(i).text());
        }
        require(!ordered.empty(), "context template requires declared block contexts");
        std::map<std::uint64_t, const rfc::storage::LocalBlockDefinition *> blocks;
        for (const auto &block : source.endpoint()->local_blocks) blocks.emplace(block.block_key, &block);
        for (const auto &id : ordered) {
            require(context_keys.contains(id), "template context must use its declared block carrier");
            const auto key = context_keys.at(id);
            require(blocks.contains(key), "template block missing");
            contexts_.push_back({id, blocks.at(key)});
            for (const auto &method : methods_)
                count_ = rfc::machine::checked_size_add(count_, method == "local_sections" ? 1 : fields_.size());
        }
    }
    std::size_t size() const noexcept { return count_; }
    void visit(const Consumer &consume) const {
        if (request_.has("queries")) {
            const auto authored = request_.at("queries");
            for (std::size_t i = 0; i < authored.size(); ++i)
                consume(authored.at(i).at("id").text(), render_record_view(authored.at(i).at("query")));
            return;
        }
        for (const auto &context : contexts_)
            for (const auto &method : methods_) {
                const auto query = [&] {
                    auto q = Json::object();
                    q.set("operation", Json::string(method)); q.set("context", Json::string(context.id));
                    return q;
                };
                if (method != "local_sections") {
                    for (const auto &field : fields_) {
                        auto q = query();
                        if (method == "green_action") q.set("lambda", Json::string(lambda_));
                        q.set("field", Json::string(field.id));
                        consume(context.id + ":" + method + ":" + field.id, q);
                    }
                    continue;
                }
                auto q = query();
                auto stalks = Json::array(), maps = Json::array(), sections = Json::array();
                for (const auto &cell : context.block->cells) {
                    auto item = Json::object(), axes = Json::array(); axes.push(Json::string("scalar"));
                    item.set("grade", Json::integer(cell.grade.value)); item.set("key", Json::string(std::to_string(cell.key)));
                    item.set("axes", axes); stalks.push(item);
                }
                for (const auto &occurrence : context.block->occurrences) {
                    auto item = Json::object(); item.set("grade", Json::integer(occurrence.upper_grade.value));
                    item.set("key", Json::string(std::to_string(occurrence.key.value)));
                    item.set("matrix", Json::parse(R"([["1"]])")); maps.push(item);
                }
                for (const auto &field : fields_) {
                    auto item = Json::object(), names = Json::array(), axes = Json::array();
                    names.push(Json::string(field.id)); axes.push(Json::string("scalar"));
                    item.set("id", Json::string(field.id)); item.set("grade", Json::integer(0)); item.set("fields", names);
                    item.set("axes", axes); item.set("quantity", field.meaning.at("quantity"));
                    item.set("units", field.meaning.at("units")); sections.push(item);
                }
                q.set("construction", Json::string("cellular_sheaf")); q.set("metric", Json::string("identity"));
                q.set("stalks", stalks); q.set("maps", maps); q.set("sections", sections);
                consume(context.id + ":" + method, q);
            }
    }
};
} // namespace
int main(int argc, char **argv) try {
    using namespace gadget;
    using N = NativeRecordValue;
    require(argc == 4 || (argc == 5 && std::string_view(argv[4]) == "--profile"),
            "usage: gadget_scientific_sequence WORKSPACE.rfcxs REQUEST.rfcxs OUTPUT.rfcxs [--profile]");
    require(!fs::exists(argv[3]), "scientific sequence output already exists");
    const auto start = profiling::Clock::now(); auto stage = start;
    const NativeRecordObserver phase = [&](std::string_view name) {
        if (argc == 5) profiling::phase(name, stage);
        stage = profiling::Clock::now();
    };
    const auto request = NativeRecord::read(argv[2]);
    require(request.at("schema").text() == "gadget.scientific.sequence.v1", "scientific sequence schema");
    phase("request");
    AssayWorkspace workspace(argv[1]); phase("workspace");
    const auto source = workspace.study(request.at("assay").text()); phase("study_admission");
    const QuerySource queries(*source, request); phase("query_planning");
    require(queries.size(), "scientific sequence requires queries");
    ScientificSession session(source);
    {
        std::set<std::string> identities;
        queries.visit([&](const auto &id, const auto &query) {
            require(identities.insert(id).second, "duplicate sequence query identity");
            session.validate(query);
        });
    }
    phase("validation");
    N::PackedArray answers{queries.size(), [&](const auto &emit) {
        queries.visit([&](const auto &id, const auto &query) {
            N::Object answer;
            answer.members.emplace_back("id", N(id));
            answer.members.emplace_back("certificate", import_record_view(session.execute(query)));
            emit(N(std::move(answer)));
        });
        phase("execution_and_encoding");
    }};
    N::Object result;
    result.members.emplace_back("schema", N(std::string("gadget.scientific.sequence.result.v1")));
    result.members.emplace_back("source_digest", N(source->digest()));
    result.members.emplace_back("answers", N(std::move(answers)));
    const auto publication = retain_operation(source->native_publication(), request, N(std::move(result)), phase);
    fs::create_directories(fs::absolute(argv[3]).parent_path());
    publication.write(argv[3]); phase("native_write");
    if (argc == 5) profiling::phase("total", start);
    std::cout << "PASS\nqueries\t" << queries.size() << '\n';
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
