#include "../examples/scientific/fixture.hpp"
#include <gadget/assay.hpp>
#include <gadget/record_view.hpp>
#include <iostream>
#include <rfc/storage/rfcx.hpp>

using namespace gadget;
using namespace rfc::storage;
using Rational = rfc::numeric::Rational;
namespace {
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition)
        throw std::runtime_error(message);
}
const TemporalTensor *tensor(const NativeRecord &owner, NativeRecord descriptor) {
    const auto key = descriptor.kind() == NativeRecordKind::Object
                         ? descriptor.at("tensor").integer()
                         : descriptor.integer();
    for (const auto &field : owner.endpoint()->tensors)
        if (field.key == key)
            return &field;
    throw std::runtime_error("native field absent");
}
const TemporalCarrier *carrier(const NativeRecord &owner, std::uint64_t key) {
    for (const auto &item : owner.endpoint()->carriers)
        if (item.key == key)
            return &item;
    throw std::runtime_error("native carrier absent");
}
Rational scalar(const NativeRecord &owner, NativeRecord descriptor) {
    const auto &value = tensor(owner, descriptor)->tensor;
    check(value.domain().logical_size() == 1 && value.support().contains(0),
          "scalar support changed");
    return *value.state(0);
}
void original_owner(const NativeRecord &publication, const NativeRecord &original) {
    check(NativeRecord::decode(publication.at("result").at("source_study_publication").bytes())
                  .encode() == original.encode(),
          "original byte owner changed");
}
void boundary(const NativeRecord &publication, const NativeRecord &original) {
    original_owner(publication, original);
    const auto owner = publication.at("result").at("fields");
    check(encode_rfcx(owner.endpoint()->source).bytes ==
              encode_rfcx(original.endpoint()->source).bytes,
          "native boundary changed the source complex");
    const auto view = owner.at("view"), binding = owner.at("bindings").at(0);
    check(view.at("schema").text() == "gadget.scientific.native.boundary.result.v1",
          "native boundary missing");
    check(
        NativeRecord::decode(binding.at("source_publication").at("publication").bytes()).encode() ==
            original.encode(),
        "field action lost actual source publication");
    const auto inputs = binding.at("inputs"), outputs = binding.at("outputs");
    const auto action_grade = view.at("field_set").at("grade").integer().get_ui();
    const auto slot_count = original.endpoint()
                                ->source.grade({static_cast<std::uint16_t>(action_grade)})
                                .incidence_count();
    check(binding.at("occurrences").size() == slot_count &&
              view.at("retained_occurrences").size() == slot_count,
          "inactive or ordered occurrences missing");
    std::map<StableCellKey, const TemporalTensor *> source_fields;
    std::map<StableOccurrenceKey, const TemporalTensor *> source_maps;
    std::map<StableCellKey, rfc::numeric::RationalVector> expected;
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        const auto input = inputs.at(i);
        const auto *field = tensor(owner, input.at("field"));
        const auto *source_field = tensor(original, input.at("field"));
        check(temporal_endpoint_detail::tensor_equal(*field, *source_field) &&
                  field->owner.kind == TemporalOwnerKind::Cell &&
                  field->owner.cell.key == input.at("cell_key").integer() &&
                  field->output == input.at("carrier").integer().get_ui(),
              "input pin or cell carrier changed");
        check(source_fields.emplace(source_field->owner.cell, source_field).second,
              "input cell repeated");
    }
    for (const auto &field : original.endpoint()->tensors)
        if (field.role == TemporalTensorRole::CoefficientIndependentMap ||
            field.role == TemporalTensorRole::RealizedIncidenceMap)
            check(source_maps.emplace(field.owner.occurrence, &field).second,
                  "original map repeated");
    for (std::size_t i = 0; i < outputs.size(); ++i) {
        const auto output = outputs.at(i);
        const auto *field = tensor(owner, view.at("result").at(i).at("values"));
        const auto *source_carrier = carrier(original, output.at("carrier").integer().get_ui());
        check(field->owner.kind == TemporalOwnerKind::Cell &&
                  field->owner.cell == *source_carrier->cell &&
                  field->output == source_carrier->key &&
                  field->tensor.domain() == source_carrier->domain &&
                  field->tensor.support().count() == field->tensor.domain().logical_size() &&
                  tensor(owner, output.at("field"))->key == field->key,
              "result cell, original carrier or complete support changed");
        check(expected
                  .emplace(field->owner.cell,
                           rfc::numeric::RationalVector(field->tensor.domain().logical_size()))
                  .second,
              "output cell repeated");
    }
    const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(action_grade)};
    const auto &arena = original.endpoint()->source.grade(grade);
    for (std::size_t slot = 0; slot < arena.incidence_count(); ++slot) {
        const StableOccurrenceKey occurrence{grade, arena.occurrence_identities[slot]};
        const auto upper =
            stable_cell_key(original.endpoint()->source, {grade, arena.upper_index[slot]});
        const auto lower = stable_cell_key(
            original.endpoint()->source,
            {{static_cast<std::uint16_t>(action_grade - 1)}, arena.lower_index[slot]});
        const auto coefficient = original.endpoint()->source.incidence_coefficient(grade, slot);
        const auto active =
            arena.aspects.get(original.endpoint()->source.signature().existence_aspect(), slot);
        const auto item = binding.at("occurrences").at(slot);
        check(item.at("key").integer() == occurrence.key.value &&
                  item.at("grade").integer() == action_grade &&
                  item.at("upper_key").integer() == upper.key &&
                  item.at("lower_key").integer() == lower.key &&
                  scalar(owner, item.at("coefficient")) == coefficient &&
                  item.at("active").boolean() == active,
              "original occurrence identity, coefficient or activity changed");
        const auto map = source_maps.find(occurrence);
        if (map == source_maps.end())
            check(item.at("map").kind() == NativeRecordKind::Null, "implicit map changed");
        else
            check(temporal_endpoint_detail::tensor_equal(*tensor(owner, item.at("map")),
                                                         *map->second) &&
                      item.at("role").integer() == static_cast<unsigned>(map->second->role),
                  "original occurrence map or role changed");
        if (!active)
            continue;
        const auto &values = source_fields.at(upper)->tensor;
        auto &target = expected.at(lower);
        for (std::size_t row = 0; row < target.size(); ++row)
            if (map == source_maps.end())
                target[row] += coefficient * *values.state(row);
            else
                for (std::size_t column = 0; column < values.domain().logical_size(); ++column)
                    target[row] +=
                        *map->second->tensor.state(row * values.domain().logical_size() + column) *
                        *values.state(column) *
                        (map->second->role == TemporalTensorRole::RealizedIncidenceMap
                             ? Rational{1}
                             : coefficient);
    }
    for (std::size_t i = 0; i < outputs.size(); ++i) {
        const auto *field = tensor(owner, outputs.at(i).at("field"));
        const auto &reference = expected.at(field->owner.cell);
        for (std::size_t j = 0; j < reference.size(); ++j)
            check(field->tensor.state(j) == reference[j],
                  "independent original occurrence action differs");
    }
}
void portable(const NativeRecord &publication, const fs::path &path) {
    publication.write(path);
    const auto reopened = NativeRecord::read(path);
    check(reopened.encode() == publication.encode(), "native archive changed bytes");
    const auto replayed = replay_scientific_record(reopened);
    check(publication.endpoint()->tensors.size() == replayed.endpoint()->tensors.size(),
          "replay lost native fields");
    for (std::size_t i = 0; i < publication.endpoint()->tensors.size(); ++i)
        check(temporal_endpoint_detail::tensor_equal(publication.endpoint()->tensors[i],
                                                     replayed.endpoint()->tensors[i]),
              "portable replay changed native fields or support");
}
void green_structure(const NativeRecord &publication, const NativeRecord &original) {
    original_owner(publication, original);
    const auto owner = publication.at("result").at("fields"), binding = owner.at("bindings").at(0);
    check(NativeRecord::decode(binding.at("original_study_publication").bytes()).encode() ==
              original.encode(),
          "Green context lost original study declarations");
    const auto *output = tensor(owner, binding.at("output"));
    const auto *input = tensor(owner, binding.at("context_input"));
    const auto *source = tensor(original, binding.at("input_tensor"));
    const auto *realization = tensor(owner, binding.at("input_realization"));
    const auto lift = temporal_tensor_matrix(*realization);
    check(realization->output == output->output && input->output == output->output &&
              source->output == binding.at("input_carrier").integer().get_ui() &&
              lift.cols() == source->tensor.domain().logical_size() &&
              lift.rows() == input->tensor.domain().logical_size(),
          "Green original input carrier or lift changed");
    for (std::size_t i = 0; i < lift.rows(); ++i) {
        Rational expected{0};
        for (std::size_t j = 0; j < lift.cols(); ++j)
            if (lift.get(i, j) != 0) {
                check(source->tensor.support().contains(j),
                      "context input acquired missing support");
                expected += lift.get(i, j) * *source->tensor.state(j);
            }
        check(input->tensor.state(i) == expected,
              "context input differs from original tensor lift");
    }
    const auto size = input->tensor.domain().logical_size();
    check(output->tensor.domain().logical_size() == size &&
              output->tensor.support().count() == size,
          "Green omitted complete context support");
    const auto *residual = tensor(owner, binding.at("residual"));
    check(residual->output == output->output && residual->tensor.support().count() == size &&
              residual->tensor.values().nonzeros() == 0,
          "Green lost supported zero residual");
    const auto *interior = tensor(owner, binding.at("interior_lift"));
    const auto *interface = tensor(owner, binding.at("interface_lift"));
    const auto *schur = tensor(owner, binding.at("schur"));
    check(interior->input != interface->input && interior->output == output->output &&
              interface->output == output->output && schur->input == interface->input &&
              schur->output == interface->input,
          "partition axes identified by dimension");
    const auto ii = temporal_tensor_matrix(*interior), ig = temporal_tensor_matrix(*interface);
    for (std::size_t i = 0; i < size; ++i) {
        Rational covered{0};
        for (std::size_t j = 0; j < ii.cols(); ++j)
            covered += ii.get(i, j);
        for (std::size_t j = 0; j < ig.cols(); ++j)
            covered += ig.get(i, j);
        check(covered == 1, "partition lost or repeated a source coordinate");
    }
}
void green(const NativeRecord &publication, const NativeRecord &original, bool expression) {
    green_structure(publication, original);
    const auto owner = publication.at("result").at("fields"), binding = owner.at("bindings").at(0);
    const auto *output = tensor(owner, binding.at("output"));
    const rfc::numeric::RationalVector expected =
        expression ? rfc::numeric::RationalVector{Rational{1, 2}, Rational{3, 2}, Rational{0}}
                   : rfc::numeric::RationalVector{Rational{3, 4}, Rational{1, 4}, Rational{1}};
    check(output->tensor.domain().logical_size() == expected.size() &&
              output->tensor.support().count() == expected.size(),
          "Green omitted complete source coordinates");
    for (std::size_t i = 0; i < expected.size(); ++i)
        check(output->tensor.state(i) == expected[i],
              "independent declared weighted Green differs");
    const auto *residual = tensor(owner, binding.at("residual"));
    check(residual->output == output->output &&
              residual->tensor.support().count() == expected.size() &&
              residual->tensor.values().nonzeros() == 0,
          "Green lost supported zero residual");
    const auto *interior = tensor(owner, binding.at("interior_lift"));
    const auto *interface = tensor(owner, binding.at("interface_lift"));
    const auto *schur = tensor(owner, binding.at("schur"));
    check(interior->input != interface->input && interior->output == output->output &&
              interface->output == output->output && schur->input == interface->input &&
              schur->output == interface->input,
          "partition axes identified by dimension");
}
void variant(const fs::path &root, bool rectangular) {
    example::scientific_study(root);
    const auto bytes = read_text(root / "study.rfcxs");
    const auto source = decode_native_series(std::as_bytes(std::span(bytes.data(), bytes.size())));
    const auto retained = source.at(source.last_sequence());
    TemporalEndpoint endpoint(retained.endpoint());
    if (!rectangular) {
        for (auto &field : endpoint.tensors)
            if (field.owner.kind == TemporalOwnerKind::Cell && field.owner.cell.grade.value == 1)
                field.tensor = {rfc::numeric::IndexDomain::linear(2),
                                rfc::numeric::StructuralSupport::full(2),
                                rfc::numeric::SparseRationalVector(2, {})};
    } else {
        for (auto &item : endpoint.carriers)
            if (item.kind == TemporalCarrierKind::CellStalk) {
                item.coordinate_system = item.key;
                if (item.grade.value == 0)
                    item.domain = rfc::numeric::IndexDomain::linear(3);
            }
        for (std::size_t k = 1; k < endpoint.source.grade_count(); ++k) {
            const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(k)};
            const auto &arena = endpoint.source.grade(grade);
            for (std::size_t slot = 0; slot < arena.incidence_count(); ++slot) {
                const auto upper =
                    stable_cell_key(endpoint.source, {grade, arena.upper_index[slot]});
                const auto lower =
                    stable_cell_key(endpoint.source,
                                    {{static_cast<std::uint16_t>(k - 1)}, arena.lower_index[slot]});
                const auto rows = k == 1 ? 3u : 2u;
                rfc::numeric::ExactMatrix map(rows, 2);
                for (std::size_t i = 0; i < rows; ++i)
                    for (std::size_t j = 0; j < 2; ++j)
                        if ((i + j + slot) % 3 != 0)
                            map.set(
                                i, j,
                                Rational(static_cast<unsigned long>((slot + 1) * (i + 1)), j + 1));
                TemporalTensor field;
                field.key = 61000 + arena.occurrence_identities[slot].value;
                field.role = slot % 2 ? TemporalTensorRole::RealizedIncidenceMap
                                      : TemporalTensorRole::CoefficientIndependentMap;
                field.owner.kind = TemporalOwnerKind::Occurrence;
                field.owner.occurrence = {grade, arena.occurrence_identities[slot]};
                field.input = 10000 + upper.key;
                field.output = 10000 + lower.key;
                field.interpretation = "authored rectangular occurrence map";
                field.tensor = temporal_matrix_tensor(map);
                endpoint.tensors.push_back(std::move(field));
            }
        }
    }
    const auto encoded = encode_native_series(NativeTemporalSeries(std::move(endpoint), 9876));
    const std::string native(reinterpret_cast<const char *>(encoded.data()), encoded.size());
    write_text(root / "study.rfcxs", native);
    auto declaration = Json::parse(read_text(root / "study.json"));
    declaration.at("native").set("sha256", Json::string("sha256:" + sha256(native)));
    if (rectangular) {
        const auto old = declaration.at("field_sets").at(0);
        auto set = Json::object();
        for (const auto *key : {"id", "grade", "fields", "carriers"})
            set.set(key, old.at(key));
        set.set("maps", Json::string("explicit"));
        auto sets = Json::array();
        sets.push(set);
        declaration.set("field_sets", sets);
    }
    write_text(root / "study.json", declaration.dump());
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 5) {
            const auto study = std::string(argv[4]) == "workspace"
                                   ? AssayWorkspace(argv[1]).study("measurements")
                                   : std::make_shared<const ScientificStudy>(fs::path(argv[1]));
            const auto original = study->native_publication();
            ScientificSession session(study);
            const auto queries = admit_record_view(Json::parse(read_text(argv[2])));
            const fs::path output(argv[3]);
            fs::create_directories(output);
            for (const auto &name : queries.names()) {
                const auto query = queries.at(name);
                const auto publication = session.execute_record(query);
                if (query.at("operation").text() == "native_boundary")
                    boundary(publication, original);
                else {
                    check(query.at("operation").text() == "green_response",
                          "unsupported native field action");
                    green_structure(publication, original);
                }
                portable(publication, output / (name + ".rfcxs"));
            }
            std::cout << "{\"status\":\"PASS\",\"checks\":" << checks << "}\n";
            return 0;
        }
        if (argc != 2)
            throw std::invalid_argument("usage: native_field_actions OWNED_PARENT");
        const fs::path parent(argv[1]);
        fs::create_directories(parent);
        const auto root = parent / ("run-" + std::to_string(new_native_record_namespace()));
        example::scientific_study(root / "source");
        const auto study = std::make_shared<const ScientificStudy>(root / "source/study.json");
        const auto original = study->native_publication();
        ScientificSession session(study);
        std::vector<NativeRecord> saved;
        variant(root / "zero", false);
        variant(root / "rectangular", true);
        for (const auto *profile : {"source", "zero", "rectangular"}) {
            const auto selected =
                std::make_shared<const ScientificStudy>(root / profile / "study.json");
            ScientificSession boundary_session(selected);
            for (const auto *backend :
                 {"cpu_cached", "cpu_cached_tiles", "cpu_tiles", "cpu_scheduled"}) {
                auto query = Json::parse(
                    R"({"operation":"native_boundary","field_set":"synthetic.assays","backend":"cpu_cached"})");
                query.set("backend", Json::string(backend));
                if (std::string(backend) != "cpu_cached") {
                    query.set("tile_rows", Json::integer(1));
                    if (std::string(backend) != "cpu_scheduled")
                        query.set("workers", Json::integer(1));
                }
                for (unsigned repeat = 0; repeat < 2; ++repeat) {
                    const auto publication =
                        boundary_session.execute_record(admit_record_view(query));
                    boundary(publication, selected->native_publication());
                    saved.push_back(publication);
                    const auto displayed = render_record_view(publication.at("result"));
                    check(displayed.at("retained_occurrences").size() == 9,
                          "optional boundary display changed original values or full occurrences");
                    if (std::string(profile) == "zero")
                        for (std::size_t i = 0; i < displayed.at("result").size(); ++i)
                            for (std::size_t j = 0;
                                 j < displayed.at("result").at(i).at("values").size(); ++j)
                                check(displayed.at("result")
                                              .at(i)
                                              .at("values")
                                              .at(j)
                                              .at("rational")
                                              .text() == "0",
                                      "supported zero output became absent");
                }
            }
        }
        for (const auto *field : {"baseline.G", "baseline.E"}) {
            for (const auto *partition : {"[]", "[0]", "[2,0]", "[0,1,2]"}) {
                auto query = Json::parse(
                    R"({"operation":"green_response","context":"baseline","field":"baseline.G","lambda":"1","coordinates":{"id":"91000","cells":["10","20","30"]},"interior":[]})");
                query.set("field", Json::string(field));
                query.set("interior", Json::parse(partition));
                const auto publication = session.execute_record(admit_record_view(query));
                green(publication, original, std::string(field) == "baseline.E");
                saved.push_back(publication);
            }
        }
        fs::rename(root / "source", root / "moved-originals");
        fs::rename(root / "zero", root / "moved-zero");
        fs::rename(root / "rectangular", root / "moved-rectangular");
        for (std::size_t i = 0; i < saved.size(); ++i)
            portable(saved[i], root / (std::to_string(i) + ".rfcxs"));
        check(study->native_publication().encode() == original.encode(),
              "execution mutated original source");
        std::cout << "PASS " << checks << " independent native field action checks; " << root
                  << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
