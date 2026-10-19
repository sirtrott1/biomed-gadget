#include "../src/research_comparisons.hpp"
#include "../src/scientific_result_fields.hpp"
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <gadget/workflow.hpp>
#include <iostream>
#include <atomic>
#include <limits>
#include <rfc/storage/rfcx.hpp>

using namespace gadget;
using namespace rfc::storage;
using Matrix = rfc::numeric::ExactMatrix;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool value, const char *message) {++checks; if (!value) throw std::runtime_error(message);}
const TemporalTensor *field(const NativeRecord &owner, NativeRecord reference) {
    const auto key = (reference.kind() == NativeRecordKind::Object ? reference.at("tensor") : reference).integer();
    for (const auto &t : owner.endpoint()->tensors) if (rfc::numeric::Integer(std::to_string(t.key)) == key) return &t;
    throw std::runtime_error("native comparison field absent");
}
Vector values(const NativeRecord &owner, NativeRecord reference) {
    const auto &t = field(owner, reference)->tensor;
    Vector result;
    if (reference.kind() == NativeRecordKind::Object && reference.has("coordinates")) {
        for (const auto &entry : field(owner, reference.at("coordinates"))->tensor.values().dense()) {
            const auto coordinate = entry.get_num().get_ui();
            check(entry.get_den() == 1 && t.support().contains(coordinate), "comparison row support differs");
            result.push_back(*t.state(coordinate));
        }
    } else {
        check(t.support().count() == t.domain().logical_size(), "complete comparison field support lost");
        for (std::size_t i = 0; i < t.domain().logical_size(); ++i) result.push_back(*t.state(i));
    }
    return result;
}
Matrix matrix(const NativeRecord &owner, NativeRecord reference) {
    const auto *t = field(owner, reference);
    check(t->tensor.support().count() == t->tensor.domain().logical_size(), "complete comparison matrix support lost");
    return temporal_tensor_matrix(*t);
}
void audit(const NativeRecord &owner, const NativeRecord &value) {
    if (value.kind() == NativeRecordKind::Object) {
        if (value.has("schema") && value.at("schema").kind() == NativeRecordKind::Text &&
            value.at("schema").text() == "gadget.scientific.tensor.view.v1") {
            const auto *t = field(owner, value);
            check(t->input || t->output, "comparison field has no native carrier");
            if (value.has("coordinates")) (void)values(owner, value);
            return;
        }
        for (const auto &name : value.names()) audit(owner, value.at(name));
    } else if (value.kind() == NativeRecordKind::Array) {
        for (std::size_t i = 0; i < value.size(); ++i) audit(owner, value.at(i));
    }
}
void verify(const NativeRecord &result, const ScientificStudy &study) {
    const auto owner = result.at("result").at("fields"), view = owner.at("view"), bindings = owner.at("bindings");
    check(encode_rfcx(owner.endpoint()->source).bytes == encode_rfcx(study.endpoint()->source).bytes,
          "comparison changed original source complex");
    check(NativeRecord::decode(result.at("result").at("source_study_publication").bytes()).encode() ==
          study.native_publication().encode(), "comparison lost full original source history");
    audit(owner, view); audit(owner, bindings);
    const auto pins = bindings.at(0).at("fields");
    for (std::size_t i = 0; i < pins.size(); ++i) {
        const auto pin = pins.at(i);
        const auto actual = study.field_accession(pin.at("field").text());
        const auto *retained = field(owner, pin.at("source_field"));
        check(pin.at("source_namespace").text() == std::to_string(actual.source_namespace()) &&
              pin.at("source_cutoff").integer() == rfc::numeric::Integer(std::to_string(actual.sequence())), "comparison substituted annotation cutoff");
        check(retained->output == actual.carrier()->key &&
              temporal_endpoint_detail::tensor_equal(retained->tensor, *actual.tensor()), "comparison changed original field support or coefficients");
    }
    NativeRecord rows = bindings.at(1);
    if (rows.at("kind").text() != "original_comparison_rows") rows = bindings.at(2);
    const auto lift = matrix(owner, rows.at("source_lift"));
    const auto coordinates = values(owner, rows.at("row_coordinates"));
    check(lift.cols() == coordinates.size() && lift.rows() == study.endpoint()->source.grade({0}).cell_count(),
          "comparison source lift shape differs");
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        for (std::size_t j = 0; j < lift.rows(); ++j)
            check(lift.get(j, i) == (coordinates[i] == Rational(std::to_string(j)) ? 1 : 0), "comparison original row lift differs");
    }
    if (view.has("features")) {
        const auto axes = bindings.at(1);
        const auto design = matrix(owner, axes.at("design")), frame = matrix(owner, axes.at("frame")),
                   coefficients = matrix(owner, axes.at("coordinates")), kernel = matrix(owner, axes.at("kernel"));
        check(frame.multiply(coefficients) == design && design.multiply(kernel).nonzeros() == 0 &&
              kernel.cols() + frame.cols() == design.cols(), "paired complete frame or kernel differs");
        check(field(owner, axes.at("design"))->input != field(owner, axes.at("frame"))->input,
              "paired equal dimension coefficient and frame carriers merged");
        const auto scales_l = values(owner, axes.at("left_scales")), scales_r = values(owner, axes.at("right_scales"));
        const auto records = view.at("features"), pairs = view.at("pairs");
        for (std::size_t i = 0; i < records.size(); ++i) {
            const auto record = records.at(i);
            const auto row = record.at("coordinate").integer().get_ui();
            const auto left = values(owner, record.at("left")), right = values(owner, record.at("right")), delta = values(owner, record.at("delta")),
                       consensus = values(owner, record.at("consensus")), residual = values(owner, record.at("residual"));
            check(left.size() == pairs.size() && right.size() == pairs.size(), "paired unit order or arity lost");
            for (std::size_t j = 0; j < pairs.size(); ++j) {
                const auto pair = pairs.at(j);
                check(left[j] == *study.field_accession(pair.at("left").text()).tensor()->state(row) * scales_l[j] &&
                      right[j] == *study.field_accession(pair.at("right").text()).tensor()->state(row) * scales_r[j],
                      "paired normalized counts differ from original native input");
            }
            check(rfc::numeric::sub(left, right) == delta && rfc::numeric::add(consensus, residual) == delta,
                  "paired complete reconstruction differs");
            check(rfc::numeric::dot(consensus, residual) == 0, "paired residual not orthogonal");
            check(frame.apply(values(owner, record.at("consensus_coefficients"))) == consensus,
                  "paired consensus coordinates differ");
            check(rfc::numeric::dot(residual, residual) == values(owner, record.at("heterogeneity_quadrance"))[0],
                  "paired quadrance differs");
        }
    } else {
        const auto design = matrix(owner, view.at("control_design")), frame = matrix(owner, view.at("control_frame")),
                   coefficients = matrix(owner, view.at("control_coordinates")), kernel = matrix(owner, view.at("control_kernel"));
        check(frame.multiply(coefficients) == design && design.multiply(kernel).nonzeros() == 0 &&
              kernel.cols() + frame.cols() == design.cols(), "controlled complete coefficient family differs");
        check(field(owner, view.at("control_design"))->input != field(owner, view.at("control_frame"))->input,
              "controlled independent frame and declared columns merged");
        const auto y = values(owner, view.at("endpoint_input")), py = values(owner, view.at("endpoint_represented")), ry = values(owner, view.at("endpoint_residual"));
        check(rfc::numeric::add(py, ry) == y && rfc::numeric::dot(py, ry) == 0 &&
              rfc::numeric::is_zero(design.transpose().apply(ry)), "controlled endpoint reconstruction differs");
        check(design.apply(values(owner, view.at("endpoint_control_coefficients"))) == py,
              "controlled endpoint coefficients differ");
        const auto fields = view.at("fields");
        for (std::size_t i = 0; i < fields.size(); ++i) {
            const auto record = fields.at(i);
            const auto x = values(owner, record.at("input")), represented = values(owner, record.at("represented")), residual = values(owner, record.at("residual"));
            check(rfc::numeric::add(represented, residual) == x && rfc::numeric::dot(represented, residual) == 0 &&
                  rfc::numeric::is_zero(design.transpose().apply(residual)), "controlled field reconstruction differs");
            check(design.apply(values(owner, record.at("control_coefficients"))) == represented,
                  "controlled field coefficients differ");
            check(rfc::numeric::dot(residual, ry) == values(owner, record.at("pairing"))[0], "controlled pairing differs");
        }
    }
}
} // namespace
int main(int argc, char **argv) try {
    check(argc == 5, "comparison study, queries, output and mode required");
    const fs::path source(argv[1]), queries_path(argv[2]), output(argv[3]);
    auto base = std::string(argv[4]) == "study" ? std::make_shared<const ScientificStudy>(source)
                                               : AssayWorkspace(source).study("measurements");
    // Later declaration publication must not replace the field producer cutoff.
    auto study = std::make_shared<const ScientificStudy>(NativeRecord::annotate(base->native_publication().unbound_value(), base->native_publication()));
    base.reset(); fs::create_directories(output);
    const auto queries = admit_record_view(Json::parse(read_text(queries_path))); auto session = std::make_unique<ScientificSession>(study);
    std::vector<NativeRecord> retained;
    for (const auto &name : queries.names()) {
        const auto query = queries.at(name); check(has_native_scientific_publisher(query.at("operation").text()), "comparison not declared native");
        const auto result = session->execute_record(query); verify(result, *study);
        result.write(output / (name + ".rfcxs"));
        const auto reopened = NativeRecord::read(output / (name + ".rfcxs")); verify(reopened, *study);
        const auto replayed = replay_scientific_record(reopened);
        check(scientific_detail::result_fields_view(replayed.at("result").at("fields")).dump() ==
              scientific_detail::result_fields_view(reopened.at("result").at("fields")).dump(), "comparison portable replay differs");
        verify(replayed, *study);
        check(replayed.endpoint()->tensors.size() == reopened.endpoint()->tensors.size(), "replay complete tensor family differs");
        for (std::size_t i = 0; i < reopened.endpoint()->tensors.size(); ++i)
            check(temporal_endpoint_detail::tensor_equal(replayed.endpoint()->tensors[i], reopened.endpoint()->tensors[i]), "replayed actual tensor differs");
        retained.push_back(reopened);
    }
    const auto first = render_record_view(queries.at(queries.names().front()));
    const auto maximum = std::numeric_limits<std::size_t>::max();
    for (std::size_t scenario = 0; scenario < 6; ++scenario) {
        rfc::runtime::ExecutionHostPolicy policy; policy.name = "declared-comparison-test"; policy.limits = {maximum, maximum, maximum, maximum, maximum};
        auto cancelled = std::make_shared<std::atomic_bool>(scenario == 0);
        if (scenario == 1) policy.limits.rational_slots = 0;
        if (scenario == 2) policy.limits.records = 0;
        if (scenario == 3) policy.limits.integer_bits = 0;
        if (scenario == 4) policy.limits.scalar_work = 0;
        if (scenario == 5) policy.allow_allocation = false;
        bool refused = false;
        try {session->execute_comparison(first, rfc::runtime::ExecutionLease(policy, cancelled));}
        catch (const std::exception &) {refused = true;}
        check(refused, "comparison resource or cancellation policy ignored");
    }
    rfc::runtime::ExecutionHostPolicy admitted; admitted.name = "declared-comparison-test-admitted";
    admitted.limits = {maximum, maximum, maximum, maximum, maximum};
    session->execute_comparison(first, rfc::runtime::ExecutionLease(admitted));
    session.reset(); study.reset();
    for (const auto &owner : retained) audit(owner.at("result").at("fields"), owner.at("result").at("fields").at("view"));
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << ",\"requests\":" << queries.size() << "}\n";
    return 0;
} catch (const std::exception &e) {std::cerr << e.what() << '\n'; return 1;}
