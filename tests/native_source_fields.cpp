#include "../src/scientific_result_fields.hpp"
#include <gadget/assay.hpp>
#include <gadget/workflow.hpp>
#include <iostream>
#include <atomic>
#include <limits>
#include <rfc/storage/rfcx.hpp>

namespace {
using namespace gadget;
using namespace rfc::storage;
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
bool mutate_value(NativeRecordValue &value, std::string_view name,
                  const NativeRecordValue &replacement) {
    if (auto *object = std::get_if<NativeRecordValue::Object>(&value.data)) {
        for (auto &[key, item] : object->members) {
            if (key == name) { item = replacement; return true; }
            if (mutate_value(item, name, replacement)) return true;
        }
    } else if (auto *array = std::get_if<NativeRecordValue::Array>(&value.data)) {
        for (auto &item : array->elements) if (mutate_value(item, name, replacement)) return true;
    }
    return false;
}
const TemporalTensor &tensor(const NativeRecord &owner, const NativeRecord &reference) {
    const auto key = reference.at("tensor").integer();
    for (const auto &field : owner.endpoint()->tensors)
        if (field.key == key) return field;
    throw std::runtime_error("source result tensor absent");
}
void verify(const NativeRecord &result, const ScientificStudy &study, const NativeRecord &query) {
    const auto original = study.native_publication();
    const auto fields = result.at("result").at("fields");
    const auto view = fields.at("view"), binding = fields.at("bindings").at(0);
    check(result.at("result").at("schema").text() == "gadget.scientific.result.fields.v1",
          "source result used a rendered publication");
    check(encode_rfcx(fields.endpoint()->source).bytes == encode_rfcx(original.endpoint()->source).bytes,
          "source operation changed original relations or occurrences");
    check(NativeRecord::decode(result.at("result").at("source_study_publication").bytes()).encode() == original.encode(),
          "source operation lost original study bytes");
    for (std::size_t i = 0; i < original.endpoint()->tensors.size(); ++i)
        check(temporal_endpoint_detail::tensor_equal(original.endpoint()->tensors[i], fields.endpoint()->tensors[i]),
              "source field or support changed");
    const auto operation = query.at("operation").text();
    check(binding.at("operation").text() == operation, "source binding operation differs");
    if (operation == "coverage") {
        const auto pin = study.field_accession(query.at("field").text());
        const auto &input = *pin.tensor_declaration();
        check(binding.at("source_namespace").text() == std::to_string(pin.source_namespace()) &&
                  binding.at("cutoff").integer() == pin.sequence() &&
                  binding.at("carrier").integer() == *input.output &&
                  temporal_endpoint_detail::tensor_equal(tensor(fields, binding.at("input")), input),
              "coverage producer identity, carrier or cutoff changed");
        check(view.at("observations").size() == input.tensor.domain().logical_size(),
              "coverage omitted original coordinates");
        for (std::size_t i = 0; i < view.at("observations").size(); ++i) {
            const auto row = view.at("observations").at(i), reference = row.at("value");
            check(row.at("coordinate").integer() == i &&
                      row.at("reported").boolean() == input.tensor.support().contains(i) &&
                      reference.at("index").integer() == i && tensor(fields, reference).key == input.key &&
                      tensor(fields, reference).tensor.state(i) == input.tensor.state(i),
                  "coverage zero, missing support or original coordinate differs");
        }
    } else if (operation == "provenance") {
        const StableCellKey cell{{static_cast<std::uint16_t>(binding.at("cell_grade").integer().get_ui())},
                                 binding.at("cell_key").integer().get_ui()};
        std::size_t selected{};
        for (std::size_t k = 1; k < original.endpoint()->source.grade_count(); ++k) {
            const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};
            const auto &arena = original.endpoint()->source.grade(grade);
            for (std::size_t slot = 0; slot < arena.incidence_count(); ++slot) {
                const auto upper = stable_cell_key(original.endpoint()->source, {grade, arena.upper_index[slot]});
                const auto lower = stable_cell_key(original.endpoint()->source,
                    {{static_cast<std::uint16_t>(k - 1)}, arena.lower_index[slot]});
                if (upper != cell && lower != cell) continue;
                const auto identity = binding.at("occurrences").at(selected);
                const auto occurrence = view.at("occurrences").at(selected++);
                check(identity.at("grade").integer() == k && identity.at("slot").integer() == slot &&
                          identity.at("key").integer() == arena.occurrence_identities[slot].value &&
                          identity.at("upper").integer() == upper.key && identity.at("lower").integer() == lower.key,
                      "provenance occurrence identity or ordered source slot differs");
                const auto coefficient = original.endpoint()->source.incidence_coefficient(grade, slot);
                check(tensor(fields, occurrence.at("coefficient")).tensor.state(0) == coefficient &&
                          tensor(fields, occurrence.at("magnitude")).tensor.state(0) == abs(coefficient),
                      "provenance imported or changed an incidence coefficient");
                for (std::size_t a = 0; a < original.endpoint()->source.signature().aspects().size(); ++a) {
                    const auto aspect = static_cast<rfc::composite::AspectId>(a);
                    check(occurrence.at("aspects").at(original.endpoint()->source.signature().aspects().at(aspect).name).boolean() ==
                              arena.aspects.get(aspect, slot), "provenance lost an independent aspect");
                }
            }
        }
        check(selected == view.at("occurrences").size(), "provenance omitted or invented an occurrence");
    } else {
        check(operation == "inspect", "unexpected source operation");
        check(scientific_detail::result_fields_view(fields).dump() == study.catalog().dump(),
              "native catalog differs from original source declarations");
    }
}
}
int main(int argc, char **argv) try {
    check(argc == 5, "study, queries, output and source mode required");
    const fs::path path(argv[1]), output(argv[3]);
    auto base = std::string(argv[4]) == "study" ? std::make_shared<const ScientificStudy>(path)
                                               : AssayWorkspace(path).study("measurements");
    auto study = std::make_shared<const ScientificStudy>(NativeRecord::annotate(
        base->native_publication().unbound_value(), base->native_publication()));
    base.reset();
    const auto queries = read_record(argv[2]).native();
    ScientificSession session(study);
    fs::create_directories(output);
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name);
        check(has_native_scientific_publisher(query.at("operation").text()), "source operation undeclared native");
        const auto result = session.execute_record(query);
        verify(result, *study, query);
        result.write(output/(name+".rfcxs"));
        const auto reopened = NativeRecord::read(output/(name+".rfcxs"));
        verify(reopened, *study, query);
        verify(replay_scientific_record(reopened), *study, query);
        auto root = reopened.unbound_value();
        check(mutate_value(root, "cutoff", NativeRecordValue(rfc::numeric::Integer(999999))),
              "source binding mutation was not exercised");
        const auto changed = reopened.derive(std::move(root));
        bool refused = false;
        try { (void)replay_scientific_record(changed); }
        catch (const std::exception &) { refused = true; }
        check(refused, "source replay admitted a changed producer cutoff");
    }
    auto plan = Json::object(), steps = Json::array();
    plan.set("schema", Json::string("gadget.scientific.plan.v1")); plan.set("id", Json::string("source-fields"));
    plan.set("question", Json::string("Retain original source observations"));
    plan.set("source_digest", Json::string(study->digest()));
    for (const auto &name : queries.names()) {
        auto step = Json::object(); step.set("id", Json::string(name));
        step.set("query", render_record_view(queries.at(name))); steps.push(step);
    }
    plan.set("steps", steps);
    const auto retained = session.execute_plan_record(admit_record_view(plan));
    retained.write(output/"plan.rfcxs");
    const auto replayed = replay_scientific_record(NativeRecord::read(output/"plan.rfcxs"));
    check(replayed.at("result").at("steps").size() == queries.size(), "source plan lost ordered children");
    const auto maximum = std::numeric_limits<std::size_t>::max();
    const auto coverage = render_record_view(queries.at("coverage_A"));
    for (std::size_t scenario = 0; scenario < 4; ++scenario) {
        rfc::runtime::ExecutionHostPolicy policy;
        policy.name = "original-source-contract";
        policy.limits = {maximum, maximum, maximum, maximum, maximum};
        auto cancelled = std::make_shared<std::atomic_bool>(scenario == 0);
        if (scenario == 1) policy.limits.records = 0;
        if (scenario == 2) policy.limits.scalar_work = 0;
        if (scenario == 3) policy.allow_allocation = false;
        bool refused = false;
        try { session.execute_source(coverage, rfc::runtime::ExecutionLease(policy, cancelled)); }
        catch (const std::exception &) { refused = true; }
        check(refused, "source operation ignored declared resource or cancellation refusal");
    }
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << ",\"requests\":" << queries.size() << "}\n";
    return 0;
} catch (const std::exception &error) {std::cerr << error.what() << '\n'; return 1;}
