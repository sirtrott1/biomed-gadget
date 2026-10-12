#include <gadget/facts.hpp>
#include <gadget/record_view.hpp>
#include <gadget/workflow.hpp>
#include <iostream>

using namespace gadget;
using N = NativeRecordValue;
namespace {
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
template<class F> void rejects(F action) {
    bool failed = false;
    try { action(); } catch (const std::exception &) { failed = true; }
    check(failed, "invalid native fact selection admitted");
}
NativeRecord authored() {
    rfc::core::RelationalComplexBuilder builder(rfc::core::RelationSignature{});
    for (std::uint64_t i = 1; i <= 3; ++i) builder.declare_witness_identity(i);
    rfc::storage::TemporalEndpoint endpoint(builder.freeze());
    rfc::storage::TemporalCarrier carrier;
    carrier.key = 1; carrier.kind = rfc::storage::TemporalCarrierKind::SourceGrade;
    carrier.grade = {0}; carrier.domain = rfc::numeric::IndexDomain::linear(3);
    endpoint.carriers.push_back(carrier);
    rfc::storage::TemporalTensor field;
    field.key = 17; field.output = 1; field.role = rfc::storage::TemporalTensorRole::Field;
    field.interpretation = "Authored supported zero, nonzero and absent coordinate";
    field.tensor = {carrier.domain, rfc::numeric::StructuralSupport::sparse(3, {0, 1}),
                   rfc::numeric::SparseRationalVector(3, {{1, Rational{7, 3}}})};
    endpoint.tensors.push_back(field);
    const auto descriptor = [](const std::string &display, std::optional<std::size_t> index) {
        N::Object out{{{"schema", N(std::string("gadget.scientific.tensor.view.v1"))},
                       {"tensor", N(rfc::numeric::Integer{17})}, {"display", N(display)}}};
        if (index) out.members.emplace_back("index", N(rfc::numeric::Integer(std::to_string(*index))));
        return N(std::move(out));
    };
    N view(N::Object{{{"zero", descriptor("rational", 0)}, {"nonzero", descriptor("rational_object", 1)},
        {"missing", descriptor("rational", 2)}, {"vector", descriptor("vector", {})},
        {"status", N(false)}, {"label", N(std::string("β authored source"))},
        {"count", N(rfc::numeric::Integer("1393796574908163946345982392040522594123776"))},
        {"numeric_text", N(std::string("0"))}}});
    N fields(N::Object{{{"schema", N(std::string("gadget.scientific.fields.v1"))}, {"view", view}}});
    N result(N::Object{{{"schema", N(std::string("gadget.scientific.result.fields.v1"))},
        {"certificate", N(N::Object{{{"source_digest", N(std::string("sha256:" + sha256("authored")))}}})},
        {"fields", fields}}});
    return NativeRecord::create_from_endpoint(N(N::Object{{{"schema", N(std::string("gadget.operation.native.v1"))},
        {"result", result}}}), 976801, {}, std::move(endpoint));
}
Json request(const NativeRecord &operation) {
    auto input = Json::parse(R"({"schema":"gadget.scientific.facts.request.v2","source_digest":"","source_namespace":"","cutoff":0,"question":"Read the authored native values","facts":[]})");
    input.set("source_digest", Json::string("sha256:" + sha256("authored")));
    input.set("source_namespace", Json::string(std::to_string(operation.source_namespace())));
    input.set("cutoff", Json::integer(operation.publication_history()->last_sequence()));
    for (const auto &[name, type] : std::vector<std::pair<std::string, std::string>>{
        {"zero", "rational"}, {"nonzero", "rational"}, {"status", "boolean"}, {"label", "text"}, {"count", "count"}}) {
        auto fact = Json::parse(R"({"id":"","meaning":"Actual authored native value","from":{"kind":"native","path":[]},"path":[],"value_type":""})");
        fact.set("id", Json::string(name)); fact.at("path").push(Json::string(name));
        fact.set("value_type", Json::string(type)); input.at("facts").push(fact);
    }
    auto vector = Json::parse(input.at("facts").at(0).dump());
    vector.set("id", Json::string("original_coordinate"));
    vector.set("path", Json::parse(R"(["vector"])")); vector.set("coordinate", Json::integer(1));
    input.at("facts").push(vector);
    return input;
}
void contracts(const fs::path &root) {
    const auto operation = authored(); const auto input = request(operation);
    const auto native_request = admit_record_view(input);
    const ScientificFacts facts(operation, native_request);
    const auto catalog = facts.catalog(), reading = facts.model_view();
    check(catalog.at("facts").size() == 6 && reading.at("fact_count").integer() == 6, "native reading lost selections");
    check(catalog.at("facts").at(0).at("value").text() == "0", "supported zero changed");
    check(catalog.at("facts").at(1).at("value").text() == "7/3", "native rational changed");
    check(reading.at("facts").at(5).at("coordinate").integer() == 1, "original coordinate omitted from reading");
    auto proposal = Json::parse(R"({"schema":"gadget.scientific.claims.proposal.v1","action":"claims","source_digest":"","facts_sha256":"","question_sha256":"","claims":[]})");
    for (const auto *name : {"source_digest", "facts_sha256", "question_sha256"}) proposal.set(name, catalog.at(name));
    for (std::size_t i = 0; i < catalog.at("facts").size(); ++i) {
        const auto fact = catalog.at("facts").at(i); auto claim = Json::object();
        claim.set("id", fact.at("id")); claim.set("fact", fact.at("id"));
        claim.set("comparison", Json::string("equal")); claim.set("value", fact.at("value")); proposal.at("claims").push(claim);
    }
    check(facts.check(proposal).at("status").text() == "matched", "typed native claims did not match");
    auto wrong = Json::parse(proposal.dump()); wrong.at("claims").at(0).set("value", Json::string("1"));
    check(facts.check(wrong).at("status").text() == "contradicted", "false native literal matched");
    wrong = Json::parse(proposal.dump()); wrong.at("claims").at(0).set("fact", Json::string("undeclared"));
    rejects([&] { facts.check(wrong); });
    for (const auto key : {"source_digest", "source_namespace", "cutoff"}) {
        auto bad = Json::parse(input.dump()); bad.set(key, std::string(key) == "cutoff" ? Json::integer(999) : Json::string("foreign"));
        rejects([&] { ScientificFacts rejected(operation, admit_record_view(bad)); });
    }
    for (const auto key : {"missing", "numeric_text", "absent"}) {
        auto bad = Json::parse(input.dump()); bad.at("facts").at(0).set("path", Json::parse("[\"" + std::string(key) + "\"]"));
        rejects([&] { ScientificFacts rejected(operation, admit_record_view(bad)); });
    }
    for (const auto index : {2, 3}) {
        auto bad = Json::parse(input.dump()); bad.at("facts").at(5).set("coordinate", Json::integer(index));
        rejects([&] { ScientificFacts rejected(operation, admit_record_view(bad)); });
    }
    auto no_index = Json::parse(input.dump()); json_object_object_del(no_index.at("facts").at(5).raw(), "coordinate");
    rejects([&] { ScientificFacts rejected(operation, admit_record_view(no_index)); });
    auto duplicate = Json::parse(input.dump()); duplicate.at("facts").push(duplicate.at("facts").at(0));
    rejects([&] { ScientificFacts rejected(operation, admit_record_view(duplicate)); });
    const auto before = facts.catalog().dump(); auto detached = facts.catalog();
    detached.at("facts").at(0).set("value", Json::string("999"));
    check(before == facts.catalog().dump(), "detached view changed fact owner");
    facts.native_publication().write(root / "authored-facts.rfcxs");
    const auto reopened = NativeRecord::read(root / "authored-facts.rfcxs");
    check(reopened.encode() == facts.native_publication().encode(), "native fact bytes changed on reopening");
    check(reopened.at("result").at("catalog").at("facts").at(0).at("observed").kind() == NativeRecordKind::Rational,
          "native zero became numerical text");
    check(NativeRecord::decode(reopened.at("request_publication").bytes()).encode() == native_request.encode(),
          "native fact request owner changed");
    check(reopened.historical_publication(operation.publication_history()->last_sequence()).encode() == operation.encode(),
          "native facts changed original operation history");
}
}
int main(int argc, char **argv) {
    try {
        if (argc != 2) throw std::invalid_argument("owned test parent required");
        const fs::path parent(argv[1]); fs::create_directories(parent);
        const auto root = parent / std::to_string(new_native_record_namespace()); fs::create_directory(root);
        contracts(root);
        std::cout << "PASS: " << checks << " native fact support, source, literal, request and history checks\n";
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
