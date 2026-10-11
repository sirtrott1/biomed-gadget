#pragma once
#include <gadget/record_view.hpp>
#include <gadget/scientific.hpp>
#include <rfc/storage/rfcx.hpp>
#include <iostream>

namespace native_scientific_test {
using namespace gadget;
using Rational = rfc::numeric::Rational;
using Tensor = rfc::storage::TemporalTensor;
using Vector = rfc::numeric::RationalVector;
inline std::size_t checks{};
inline void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
inline const Tensor *tensor(const NativeRecord &owner, NativeRecord descriptor) {
    const auto key = descriptor.at("tensor").integer();
    for (const auto &item : owner.endpoint()->tensors)
        if (key == item.key) return &item;
    throw std::runtime_error("native scientific tensor absent");
}
inline Vector vector(const Tensor &value) {
    Vector result;
    for (std::size_t i = 0; i < value.tensor.domain().logical_size(); ++i) {
        check(value.tensor.support().contains(i), "complete field support absent");
        result.push_back(*value.tensor.state(i));
    }
    return result;
}
inline Rational matrix(const Tensor &value, std::size_t row, std::size_t column) {
    const auto columns = value.tensor.domain().extents().at(1);
    check(value.tensor.support().contains(row * columns + column), "complete map support absent");
    return *value.tensor.state(row * columns + column);
}
inline Vector apply_map(const Tensor &map, const Vector &input) {
    const auto shape = map.tensor.domain().extents();
    check(shape.at(1) == input.size(), "native map domain differs");
    Vector output(shape.at(0));
    for (std::size_t i = 0; i < output.size(); ++i)
        for (std::size_t j = 0; j < input.size(); ++j)
            output[i] += matrix(map, i, j) * input[j];
    return output;
}
inline NativeRecordValue &member(NativeRecordValue &value, std::string_view name) {
    for (auto &[key, item] : std::get<NativeRecordValue::Object>(value.data).members)
        if (key == name) return item;
    throw std::runtime_error("native test member absent");
}
inline void compare(const NativeRecord &a, const NativeRecord &b) {
    check(a.endpoint()->tensors.size() == b.endpoint()->tensors.size(), "replay tensor count differs");
    for (std::size_t i = 0; i < a.endpoint()->tensors.size(); ++i)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(a.endpoint()->tensors[i], b.endpoint()->tensors[i]),
              "replay native tensor differs");
    check(std::ranges::equal(a.at("request_publication").bytes(), b.at("request_publication").bytes()),
          "replay retained result as request owner");
    check(std::ranges::equal(a.at("result").at("source_study_publication").bytes(),
                            b.at("result").at("source_study_publication").bytes()), "replay changed original study bytes");
    check(a.raw_sources().size() == b.raw_sources().size(), "replay accumulated byte owners");
}
inline void hodge(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    const auto frames = view.at("frames"), reductions = view.at("reductions");
    for (const auto *name : {"gradient", "curl", "harmonic"}) {
        const auto reduction = reductions.at(name);
        const auto &original = *tensor(owner, reduction.at("original"));
        const auto &image = *tensor(owner, reduction.at("image"));
        const auto &coefficients = *tensor(owner, reduction.at("coefficients"));
        const auto &kernel = *tensor(owner, reduction.at("kernel"));
        const auto n = original.tensor.domain().extents().at(1);
        const auto r = image.tensor.domain().extents().at(1);
        const auto free = kernel.tensor.domain().extents().at(1);
        check(image.input != original.input && kernel.output == original.input &&
              coefficients.input == original.input && coefficients.output == image.input,
              "reduction carrier bindings lost");
        check(n == r + free, "complete reduction family missing");
        for (std::size_t i = 0; i < original.tensor.domain().extents().at(0); ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                Rational sum{0};
                for (std::size_t k = 0; k < r; ++k) sum += matrix(image, i, k) * matrix(coefficients, k, j);
                check(sum == matrix(original, i, j), "reduction image reconstruction differs");
            }
            for (std::size_t k = 0; k < free; ++k) {
                Rational sum{0};
                for (std::size_t j = 0; j < n; ++j) sum += matrix(original, i, j) * matrix(kernel, j, k);
                check(sum == 0, "reduction kernel incomplete");
            }
        }
    }
    for (std::size_t f = 0; f < view.at("fields").size(); ++f) {
        const auto field = view.at("fields").at(f);
        const auto &input_tensor = *tensor(owner, field.at("input"));
        const auto input = vector(input_tensor), reconstructed = vector(*tensor(owner, field.at("reconstructed")));
        Vector sum(input.size());
        std::vector<Vector> sectors;
        for (const auto *name : {"gradient", "curl", "harmonic"}) {
            const auto &sector = *tensor(owner, field.at("sectors").at(name));
            const auto &frame = *tensor(owner, frames.at(name));
            const auto &coordinates = *tensor(owner, field.at("frame_coordinates").at(name));
            check(sector.output == input_tensor.output && frame.output == input_tensor.output &&
                  coordinates.output == frame.input && coordinates.output != sector.output,
                  "Hodge field or coefficient carrier changed");
            const auto value = vector(sector);
            check(apply_map(frame, vector(coordinates)) == value, "frame coordinates do not reconstruct sector");
            for (std::size_t i = 0; i < sum.size(); ++i) sum[i] += value[i];
            sectors.push_back(value);
        }
        check(sum == input && reconstructed == input, "complete sectors do not reconstruct source field");
        const Tensor *metric = nullptr;
        const auto grade = view.at("grade").integer().get_ui();
        if (view.at("metrics").at("kind").text() != "unit") {
            const auto key = view.at("publication").at("metric_tensors").at(grade).at("key").text();
            for (const auto &item : owner.endpoint()->tensors) if (std::to_string(item.key) == key) metric = &item;
            check(metric != nullptr, "actual declared metric absent");
        }
        for (std::size_t a = 0; a < sectors.size(); ++a)
            for (std::size_t b = a + 1; b < sectors.size(); ++b) {
                Rational pairing{0};
                for (std::size_t i = 0; i < input.size(); ++i)
                    for (std::size_t j = 0; j < input.size(); ++j)
                        pairing += sectors[a][i] * (metric ? matrix(*metric, i, j) : Rational(i == j)) * sectors[b][j];
                check(pairing == 0, "retained sectors are not metric orthogonal");
            }
        const auto &constraints = *tensor(owner, reductions.at("harmonic").at("original"));
        for (const auto &value : apply_map(constraints, sectors.at(2))) check(value == 0, "harmonic constraints fail");
    }
}
inline void observation(const NativeRecord &publication) {
    const auto owner = publication.at("result").at("fields"), view = owner.at("view");
    const auto &operation = *tensor(owner, view.at("operator_map"));
    const auto &accession = *tensor(owner, view.at("accession_map"));
    check(operation.input == operation.output && accession.input == operation.output,
          "observation map carrier binding lost");
    if (view.at("observation").at("kind").text() == "coordinates") {
        const auto coordinates = view.at("observation").at("coordinates");
        for (std::size_t i = 0; i < coordinates.size(); ++i)
            for (std::size_t j = 0; j < operation.tensor.domain().extents().at(0); ++j)
                check(matrix(accession, i, j) == Rational(coordinates.at(i).integer() == j),
                      "ordered selector changed original coordinate");
    }
    for (std::size_t i = 0; i < view.at("fields").size(); ++i) {
        const auto field = view.at("fields").at(i);
        const auto &input = *tensor(owner, field.at("input")), &output = *tensor(owner, field.at("values"));
        check(input.output == operation.input && output.output == accession.output &&
              output.tensor.domain().logical_size() == accession.tensor.domain().extents().at(0),
              "observed field lost actual accession carrier");
        const auto values = vector(output), input_values = vector(input);
        for (std::size_t row = 0; row < values.size(); ++row) {
            bool zero = true;
            for (std::size_t column = 0; column < input_values.size(); ++column)
                zero = zero && matrix(accession, row, column) == 0;
            if (zero) check(values[row] == 0, "zero observation row became absent or nonzero");
        }
    }
}
inline void run(const fs::path &root) {
    checks = 0;
    auto study = std::make_shared<const ScientificStudy>(root / "study/study.json");
    ScientificSession session(study);
    const auto query = admit_record_view(Json::parse(read_text(root / "query.json")));
    const auto original = study->native_publication();
    const auto result = session.execute_record(query);
    const auto fields = result.at("result").at("fields");
    check(fields.at("schema").text() == "gadget.scientific.fields.v1", "scientific result has no typed fields");
    check(rfc::storage::encode_rfcx(fields.endpoint()->source).bytes == rfc::storage::encode_rfcx(original.endpoint()->source).bytes,
          "scientific result replaced original complex");
    for (std::size_t i = 0; i < original.endpoint()->tensors.size(); ++i)
        check(rfc::storage::temporal_endpoint_detail::tensor_equal(original.endpoint()->tensors[i], fields.endpoint()->tensors[i]),
              "source field changed during publication");
    const bool is_hodge = query.at("operation").text() == "hodge_frame_decomposition";
    if (is_hodge) hodge(result); else observation(result);
    const auto archive = NativeRecord::decode(result.at("result").at("source_study_publication").bytes());
    check(archive.encode() == original.encode(), "publication changed original study archive bytes");
    fs::create_directories(root / "native");
    result.write(root / "native/result.rfcxs");
    auto plan = Json::object(), steps = Json::array(), step = Json::object();
    plan.set("schema", Json::string("gadget.scientific.plan.v1")); plan.set("id", Json::string("native-fields"));
    plan.set("question", Json::string("Retain complete source fields")); plan.set("source_digest", Json::string(study->digest()));
    step.set("id", Json::string("first")); step.set("query", render_record_view(query)); steps.push(Json::parse(step.dump()));
    step.set("id", Json::string("second")); steps.push(step); plan.set("steps", steps);
    const auto planned = session.execute_plan_record(admit_record_view(plan)); planned.write(root / "native/plan.rfcxs");
    fs::rename(root / "study", root / "study.offline");
    try {
        const auto replayed = replay_scientific_record(NativeRecord::read(root / "native/result.rfcxs"));
        compare(result, replayed); compare(replayed, replay_scientific_record(replayed));
        if (is_hodge) hodge(replayed); else observation(replayed);
        const auto repeated_plan = replay_scientific_record(NativeRecord::read(root / "native/plan.rfcxs"));
        check(repeated_plan.at("result").at("steps").size() == 2, "plan replay lost native step");
        check(planned.raw_sources().size() == repeated_plan.raw_sources().size(), "plan replay accumulated byte owners");
    } catch (...) {fs::rename(root / "study.offline", root / "study"); throw;}
    fs::rename(root / "study.offline", root / "study");
    const auto descriptor = fields.at("view").at("fields").at(0).at(is_hodge ? "sectors" : "values");
    const auto key = (is_hodge ? descriptor.at("gradient") : descriptor).at("tensor").integer();
    const auto target = std::find_if(result.endpoint()->tensors.begin(), result.endpoint()->tensors.end(),
                                    [&](const auto &item) {return item.key == key;});
    if (target->tensor.domain().logical_size()) {
        for (const bool remove_support : {false, true}) {
            auto endpoint = *result.endpoint();
            auto &changed = endpoint.tensors[static_cast<std::size_t>(target - result.endpoint()->tensors.begin())].tensor;
            std::vector<rfc::numeric::SparseRationalEntry> values;
            for (const auto &entry : changed.values().entries()) if (entry.index != 0) values.push_back(entry);
            if (!remove_support) {values.push_back({0, *changed.state(0) + 1}); std::sort(values.begin(), values.end(), [](const auto &a, const auto &b) {return a.index < b.index;});}
            std::vector<std::size_t> support; changed.support().for_each_active([&](auto i) {if (!remove_support || i) support.push_back(i);});
            changed = rfc::numeric::SparseRationalTensor(changed.domain(), rfc::numeric::StructuralSupport::sparse(changed.domain().logical_size(), support),
                rfc::numeric::SparseRationalVector(changed.domain().logical_size(), values));
            endpoint.declarations.reset();
            const auto forged = NativeRecord::create_from_endpoint(result.unbound_value(), new_native_record_namespace(), result.raw_sources(), std::move(endpoint));
            if (!remove_support) forged.write(root / "native/changed-field.rfcxs");
            bool refused = false;
            try {(void)replay_scientific_record(forged);} catch (const std::invalid_argument &) {refused = true;}
            check(refused, "replay accepted changed actual field or support");
        }
    }
    auto changed = result.unbound_value();
    auto &bindings = member(member(member(changed, "result"), "fields"), "bindings");
    member(std::get<NativeRecordValue::Array>(bindings.data).elements.at(0), "source_sequence") = NativeRecordValue(std::string("unbound cutoff"));
    bool refused = false;
    try {(void)replay_scientific_record(result.derive(std::move(changed)));} catch (const std::invalid_argument &) {refused = true;}
    check(refused, "replay accepted changed source cutoff binding");
    std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
}
} // namespace native_scientific_test
