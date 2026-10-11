#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <rfc/storage/rfcx.hpp>
#include <iostream>

using namespace gadget;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
const rfc::storage::TemporalTensor *field(const NativeRecord &owner, NativeRecord descriptor) {
    const auto key = (descriptor.kind() == NativeRecordKind::Object ? descriptor.at("tensor") : descriptor).integer();
    for (const auto &tensor : owner.endpoint()->tensors) if (key == tensor.key) return &tensor;
    throw std::runtime_error("retained field absent");
}
Rational state(const rfc::numeric::SparseRationalTensor &tensor, std::size_t i) {
    check(tensor.support().contains(i), "required coefficient absent");
    return tensor.state(i).value();
}
NativeRecordValue &member(NativeRecordValue &value, const std::string &name) {
    for (auto &[key, item] : std::get<NativeRecordValue::Object>(value.data).members)
        if (key == name) return item;
    throw std::runtime_error("fixture member absent");
}
void verify_model(const NativeRecord &publication, const NativeRecord &original) {
    const auto fields = publication.at("result").at("fields");
    check(fields.at("schema").text() == "gadget.scientific.fields.v1", "model has no native owner");
    check(rfc::storage::encode_rfcx(fields.endpoint()->source).bytes ==
          rfc::storage::encode_rfcx(original.endpoint()->source).bytes, "original complex changed");
    const auto view = fields.at("view");
    for (std::size_t s = 0; s < view.at("strata").size(); ++s) {
        const auto stratum = view.at("strata").at(s);
        if (stratum.at("status").text() != "completed") continue;
        const auto rows = stratum.at("coordinates");
        const auto &metric = *field(fields, stratum.at("metric_diagonal"));
        check(metric.tensor.support().count() == rows.size(), "selected support replaced by ambient support");
        for (std::size_t m = 0; m < stratum.at("models").size(); ++m) {
            const auto model = stratum.at("models").at(m);
            const auto &design = *field(fields, model.at("design"));
            const auto &frame = *field(fields, model.at("independent_frame"));
            const auto &kernel = *field(fields, model.at("coefficient_kernel"));
            const auto columns = design.tensor.domain().extents()[1];
            const auto free = kernel.tensor.domain().extents()[1];
            check(design.output == frame.output && design.input != frame.input &&
                  design.input == kernel.output && kernel.input != kernel.output,
                  "parameter and frame axes merged by dimension");
            check(columns == model.at("columns").size() && free == model.at("coefficient_nullity").integer(),
                  "complete coefficient family missing");
            for (std::size_t i = 0; i < rows.size(); ++i)
                for (std::size_t k = 0; k < free; ++k) {
                    Rational sum{0};
                    for (std::size_t j = 0; j < columns; ++j)
                        sum += state(design.tensor, i * columns + j) * state(kernel.tensor, j * free + k);
                    check(sum == 0, "coefficient family changes represented field");
                }
            for (std::size_t r = 0; r < model.at("responses").size(); ++r) {
                const auto response = model.at("responses").at(r);
                const auto &beta = *field(fields, response.at("particular_coefficients"));
                const auto &fitted = *field(fields, response.at("represented"));
                const auto &residual = *field(fields, response.at("residual"));
                const auto &input = *field(fields, response.at("input"));
                const auto &normal = *field(fields, response.at("normal_equations"));
                check(fitted.output == metric.output && fitted.tensor.domain() == metric.tensor.domain() &&
                      fitted.tensor.support().count() == rows.size(), "fit lost original carrier or support");
                Rational quadrance{0};
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    const auto source = rows.at(i).integer().get_ui();
                    Rational value{0};
                    for (std::size_t j = 0; j < columns; ++j)
                        value += state(design.tensor, i * columns + j) * state(beta.tensor, j);
                    check(value == state(fitted.tensor, source), "native design action differs");
                    const auto error = state(residual.tensor, source);
                    check(value + error == state(input.tensor, source), "field reconstruction differs");
                    quadrance += error * error * state(metric.tensor, source);
                }
                check(quadrance == state(field(fields, response.at("residual_quadrance"))->tensor, 0),
                      "retained quadrance differs");
                for (std::size_t j = 0; j < columns; ++j) {
                    Rational pairing{0};
                    for (std::size_t i = 0; i < rows.size(); ++i) {
                        const auto source = rows.at(i).integer().get_ui();
                        pairing += state(design.tensor, i * columns + j) * state(metric.tensor, source) *
                                   state(residual.tensor, source);
                    }
                    check(pairing == 0 && state(normal.tensor, j) == 0,
                          "orthogonality or supported zero differs");
                }
            }
        }
    }
}
void verify_replay(const NativeRecord &expected, const NativeRecord &actual) {
    const auto a = expected.endpoint(), b = actual.endpoint();
    check(a->tensors.size() == b->tensors.size(), "replay field family differs");
    for (std::size_t i = 0; i < a->tensors.size(); ++i)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(a->tensors[i], b->tensors[i]),
              "replay actual tensor differs");
    check(std::ranges::equal(expected.at("request_publication").bytes(), actual.at("request_publication").bytes()),
          "replay retained a result archive as the request owner");
    check(std::ranges::equal(expected.at("result").at("source_study_publication").bytes(),
                            actual.at("result").at("source_study_publication").bytes()),
          "replay changed or extended the original study archive");
    check(expected.raw_sources().size() == actual.raw_sources().size(), "replay retained extra byte owners");
}
}
int main(int argc, char **argv) try {
    if (argc != 5) throw std::invalid_argument("workspace, queries, output and mode required");
    const fs::path workspace_path(argv[1]), queries_path(argv[2]), output(argv[3]);
    fs::create_directories(output);
    const auto study = std::string(argv[4]) == "study"
        ? std::make_shared<const ScientificStudy>(workspace_path)
        : AssayWorkspace(workspace_path).study("measurements");
    const auto original = study->native_publication();
    ScientificSession session(study);
    const auto queries = admit_record_view(Json::parse(read_text(queries_path)));
    std::vector<NativeRecord> results;
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name);
        const auto result = session.execute_record(query);
        if (std::string(argv[4]) != "finite") verify_model(result, original);
        else {
            const auto fields = result.at("result").at("fields");
            const auto binding = fields.at("bindings").at(0);
            const auto &statistics = *field(fields, binding.at("statistics"));
            check(statistics.input != statistics.output, "reference law axes merged");
            if (name == "undefined") check(statistics.tensor.support().count() < statistics.tensor.domain().logical_size(),
                                           "undefined statistic became zero");
        }
        const auto retained = NativeRecord::decode(result.at("result").at("source_study_publication").bytes());
        check(retained.encode() == original.encode(), "original study bytes changed");
        result.write(output / (name + ".rfcxs"));
        results.push_back(NativeRecord::read(output / (name + ".rfcxs")));
    }
    auto plan = Json::object(), steps = Json::array();
    plan.set("schema", Json::string("gadget.scientific.plan.v1"));
    plan.set("id", Json::string("native-models")); plan.set("question", Json::string("Retain declared model fields"));
    plan.set("source_digest", Json::string(study->digest()));
    for (const auto &name : queries.names()) {
        auto step = Json::object(); step.set("id", Json::string(name)); step.set("query", render_record_view(queries.at(name))); steps.push(step);
    }
    plan.set("steps", steps);
    const auto planned = session.execute_plan_record(admit_record_view(plan));
    planned.write(output / "plan.rfcxs");
    const auto moved = workspace_path.parent_path().string() + ".moved";
    fs::rename(workspace_path.parent_path(), moved);
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto replayed = replay_scientific_record(results[i]);
        verify_replay(results[i], replayed);
        replayed.write(output / ("replay-" + std::to_string(i) + ".rfcxs"));
        const auto repeated = replay_scientific_record(replayed);
        verify_replay(replayed, repeated);
    }
    const auto replayed_plan = replay_scientific_record(NativeRecord::read(output / "plan.rfcxs"));
    check(replayed_plan.at("result").at("steps").size() == queries.size(), "plan replay lost step");
    check(std::ranges::equal(planned.at("request_publication").bytes(), replayed_plan.at("request_publication").bytes()),
          "plan replay retained the enclosing result as request owner");
    check(planned.raw_sources().size() == replayed_plan.raw_sources().size(), "plan replay retained extra byte owners");
    replayed_plan.write(output / "replay-plan.rfcxs");
    auto changed_receipt = results.front().unbound_value();
    member(member(member(changed_receipt, "result"), "certificate"), "gadget_sources_sha256") =
        NativeRecordValue(std::string("different software receipt"));
    const auto relocated = results.front().derive(std::move(changed_receipt));
    verify_replay(results.front(), replay_scientific_record(relocated));
    auto changed_request = results.front().unbound_value();
    member(member(changed_request, "request"), "operation") = NativeRecordValue(std::string("unretained_operation"));
    bool request_refused = false;
    try {(void)replay_scientific_record(results.front().derive(std::move(changed_request)));}
    catch (const std::invalid_argument &error) {
        request_refused = std::string(error.what()).find("request owner") != std::string::npos;
    }
    check(request_refused, "query detached from original request owner was accepted");
    auto endpoint = *results.front().endpoint();
    auto &tensor = endpoint.tensors.back().tensor;
    const auto size = tensor.domain().logical_size();
    check(size > 0 && tensor.support().contains(0), "tamper fixture has no supported field");
    auto entries = tensor.values().entries();
    std::vector<rfc::numeric::SparseRationalEntry> changed(entries.begin(), entries.end());
    const auto found = std::find_if(changed.begin(), changed.end(), [](const auto &entry) {return entry.index == 0;});
    if (found == changed.end()) changed.insert(changed.begin(), {0, Rational{1}});
    else found->value += 1;
    tensor = rfc::numeric::SparseRationalTensor(tensor.domain(), tensor.support(),
        rfc::numeric::SparseRationalVector(size, std::move(changed)));
    endpoint.declarations.reset();
    const auto forged = NativeRecord::create_from_endpoint(results.front().unbound_value(), new_native_record_namespace(),
                                                           results.front().raw_sources(), std::move(endpoint));
    forged.write(output / "changed-field.rfcxs");
    bool refused = false;
    try {(void)replay_scientific_record(forged);} catch (const std::invalid_argument &) {refused = true;}
    check(refused, "changed field replay accepted with unchanged source digest");
    auto absent_endpoint = *results.front().endpoint();
    auto &absent = absent_endpoint.tensors.back().tensor;
    auto stored = absent.values().entries();
    std::vector<rfc::numeric::SparseRationalEntry> retained;
    for (const auto &entry : stored) if (entry.index != 0) retained.push_back(entry);
    std::vector<std::size_t> support;
    absent.support().for_each_active([&](std::size_t i) {if (i != 0) support.push_back(i);});
    absent = rfc::numeric::SparseRationalTensor(absent.domain(),
        rfc::numeric::StructuralSupport::sparse(absent.domain().logical_size(), std::move(support)),
        rfc::numeric::SparseRationalVector(absent.domain().logical_size(), std::move(retained)));
    absent_endpoint.declarations.reset();
    const auto missing = NativeRecord::create_from_endpoint(results.front().unbound_value(), new_native_record_namespace(),
        results.front().raw_sources(), std::move(absent_endpoint));
    refused = false;
    try {(void)replay_scientific_record(missing);} catch (const std::invalid_argument &) {refused = true;}
    check(refused, "changed support replay accepted");
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
    return 0;
} catch (const std::exception &error) {std::cerr << error.what() << '\n'; return 1;}
