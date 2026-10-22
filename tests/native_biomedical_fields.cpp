#include "../src/pipeline_answers.hpp"
#include "../src/source_grade_fields.hpp"
#include <gadget/assay.hpp>
#include <gadget/pipeline.hpp>
#include <gadget/record_view.hpp>
#include <rfc/storage/rfcx.hpp>
#include <rfc/storage/temporal_finite.hpp>
#include <iostream>
#include <set>

namespace {
using namespace gadget;
using namespace rfc;
std::size_t checks{};
void check(bool value, const char *message) {
    ++checks;
    if (!value) throw std::logic_error(message);
}
const storage::TemporalTensor &tensor(const NativeRecord &record, const char *name) {
    return assay_detail::retained_tensor(record, record.at(name));
}
void finite() {
    auto catalog = Json::parse(R"({"source_digest":"original-source","intent_sha256":"original-intent",
        "recipes":[{"id":"one","method":"coverage","goal_sha256":"original-goal"}]})");
    const auto intent = Json::parse(R"({"answer_contract":{"schema":"gadget.pipeline.answer.contract.v1",
        "question":"Are the declared measured coordinates supported?","requirements":[
        {"id":"first","observable":"field_support","bindings":["one"]},
        {"id":"second","observable":"field_support","bindings":["one"]}]}})");
    for (bool supplied : {false, true}) {
        auto selection = Json::array();
        if (supplied) selection.push(Json::string("one"));
        const auto result = pipeline_detail::answer_admission_record(intent, catalog, &selection);
        const auto reopened = NativeRecord::decode(result.encode());
        const auto &endpoint = *reopened.endpoint();
        check(endpoint.source.grade({0}).cell_count() == 2 && endpoint.finite_systems.size() == 1,
              "finite result lost actual obligation carrier");
        const auto native = storage::temporal_finite_constraint_system<std::int64_t>(endpoint, 1);
        const auto section = storage::temporal_finite_section<std::int64_t>(endpoint, 1, 1);
        check(native.globally_valid(section) == supplied && native.charts().size() == 4,
              "finite native section or charts changed");
        check(render_record_view(result).dump() == render_record_view(reopened).dump(),
              "finite result roundtrip changed display");
        if (!supplied) {
            algebra::SubcomplexMask patch(endpoint.source); patch.select({{0}, 0});
            const auto core = native.exact_patch_unsat_core(patch);
            check(core && core->subset_minimal && core->chart_indices == std::vector<std::size_t>{0, 1},
                  "finite result lost minimal conflicting chart pair");
        }
    }
}
void gluing(const AssayWorkspace &workspace, const NativeRecord &fields, const NativeRecord &request) {
    const auto &source = fields.endpoint()->source;
    for (const auto *role : {"common", "left", "right"}) {
        const auto control = request.at("sources").at(role);
        const auto original = workspace.study(control.at("assay").text(), control.at("partition").integer().get_ui());
        const auto publication = NativeRecord::decode(fields.at("participant_publications").at(role).bytes());
        check(publication.encode() == original->native_publication().encode(), "gluing original participant changed");
        const ScientificStudy admitted(publication);
        check(publication.source_namespace() == original->native_publication().source_namespace() &&
                  admitted.endpoint()->sequence == original->endpoint()->sequence &&
                  fields.at("source_bindings").at(role).at("source_sequence").text() ==
                      std::to_string(admitted.endpoint()->sequence), "gluing cutoff binding changed");
    }
    bool supported_zero = false, missing = false;
    for (std::size_t s = 0; s < fields.at("sections").size(); ++s) {
        const auto lineage = fields.at("sections").at(s).at("lineage");
        for (std::size_t j = 0; j < lineage.size(); ++j) {
            const auto chain = lineage.at(j);
            const auto &field = tensor(chain, "field");
            const auto grade = chain.at("grade").integer().get_ui();
            const auto carrier = std::find_if(fields.endpoint()->carriers.begin(), fields.endpoint()->carriers.end(),
                [&](const auto &c) { return c.key == field.output; });
            check(carrier != fields.endpoint()->carriers.end() && carrier->kind == storage::TemporalCarrierKind::SourceGrade &&
                      carrier->grade.value == grade && field.tensor.domain().logical_size() == source.grade({static_cast<std::uint16_t>(grade)}).cell_count(),
                  "glued field is detached from actual constructed grade");
            std::set<std::size_t> expected_support;
            if (j == 0) {
                for (std::size_t i = 0; i < field.tensor.domain().logical_size(); ++i) expected_support.insert(i);
            } else {
                const auto previous = lineage.at(j - 1);
                const auto &before = tensor(previous, "field");
                const auto &upper = source.grade({static_cast<std::uint16_t>(grade + 1)});
                for (std::size_t o = 0; o < upper.incidence_count(); ++o)
                    if (before.tensor.state(upper.upper_index.at(o))) expected_support.insert(upper.lower_index.at(o));
            }
            for (std::size_t i = 0; i < field.tensor.domain().logical_size(); ++i) {
                const auto state = field.tensor.state(i);
                check(state.has_value() == expected_support.contains(i), "glued boundary lost visited zero or fabricated missing support");
                supported_zero = supported_zero || (state && *state == 0);
                missing = missing || !state;
            }
        }
    }
    check(supported_zero && missing, "gluing fixture must distinguish supported zero and missing support");
}
void transfer(const AssayWorkspace &workspace, const NativeRecord &fields, const NativeRecord &request) {
    for (const auto *role : {"left", "right"}) {
        const auto binding = fields.at(std::string(role) + "_binding");
        const auto control = request.at(role);
        const auto original = workspace.study(control.at("assay").text(), control.at("partition").integer().get_ui());
        const auto publication = NativeRecord::decode(binding.at("publication").bytes());
        check(publication.encode() == original->native_publication().encode(), "transfer original participant changed");
        const ScientificStudy admitted(publication);
        check(rfc::numeric::Integer(std::to_string(publication.source_namespace())) == binding.at("source_namespace").integer() &&
                  rfc::numeric::Integer(std::to_string(admitted.endpoint()->sequence)) == binding.at("source_sequence").integer(), "transfer original binding changed");
        const auto original_field = std::find_if(admitted.endpoint()->tensors.begin(), admitted.endpoint()->tensors.end(),
            [&](const auto &t) { return rfc::numeric::Integer(std::to_string(t.key)) == binding.at("source_field").integer(); });
        check(original_field != admitted.endpoint()->tensors.end() &&
                  original_field->output == binding.at("source_carrier").integer().get_ui(), "transfer source field carrier changed");
        const auto &stored = tensor(fields, role);
        check(stored.tensor.domain() == original_field->tensor.domain(), "transfer input domain changed");
        for (std::size_t i = 0; i < stored.tensor.domain().logical_size(); ++i)
            check(stored.tensor.state(i) == original_field->tensor.state(i), "transfer input values or support changed");
        if (std::string(role) == "left")
            check(storage::encode_rfcx(admitted.endpoint()->source).bytes == storage::encode_rfcx(fields.endpoint()->source).bytes,
                  "transfer fabricated its left source");
    }
    const auto &map = tensor(fields, "map_field");
    const auto &left = tensor(fields, "left"), &right = tensor(fields, "right");
    check(map.input == right.output && map.output == left.output && left.output != right.output,
          "transfer conflated independent carriers");
    const auto right_carrier = std::find_if(fields.endpoint()->carriers.begin(), fields.endpoint()->carriers.end(),
        [&](const auto &c) { return c.key == right.output; });
    check(right_carrier != fields.endpoint()->carriers.end() && right_carrier->kind == storage::TemporalCarrierKind::Explicit,
          "independent right axis masquerades as left source grade");
    numeric::RationalVector l, r;
    for (std::size_t i = 0; i < left.tensor.domain().logical_size(); ++i) l.push_back(left.tensor.state(i).value());
    for (std::size_t i = 0; i < right.tensor.domain().logical_size(); ++i) r.push_back(right.tensor.state(i).value());
    const auto transported = storage::temporal_tensor_matrix(map).apply(r);
    const auto residual = numeric::sub(l, transported);
    for (std::size_t i = 0; i < l.size(); ++i) {
        check(tensor(fields, "transported_right").tensor.state(i) == transported[i] &&
                  tensor(fields, "residual").tensor.state(i) == residual[i], "transfer exact action changed");
        check(tensor(fields, "metric_diagonal").tensor.state(i) == numeric::Rational{1}, "transfer identity metric changed");
    }
    check(tensor(fields, "residual_quadrance").tensor.state(0) == numeric::dot(residual, residual) &&
              tensor(fields, "pairing").tensor.state(0) == numeric::dot(l, transported), "transfer exact scalar witnesses changed");
}
void allocation(const NativeRecord &fields) {
    check(fields.at("schema").text() == "gadget.allocation.fields.v1", "allocation display used as owner");
    const auto &source = fields.endpoint()->source;
    const auto &relations = source.grade({1});
    check(relations.cell_count() == 5 && relations.incidence_count() == 13,
          "allocation lost original classes or explicit grouped target diagrams");
    const std::vector<std::vector<std::uint32_t>> originals{{0, 1, 2}, {1, 2}, {0}};
    std::size_t occurrence = 0;
    for (std::size_t c = 0; c < originals.size(); ++c) {
        check(relations.offsets.at(c + 1) - relations.offsets.at(c) == originals[c].size(),
              "allocation collapsed original candidate arity");
        for (const auto candidate : originals[c]) {
            check(relations.lower_index.at(occurrence) == candidate &&
                      relations.occurrence_identities.at(occurrence).value == occurrence + 1 &&
                      source.incidence_coefficient({1}, occurrence) == 1,
                  "allocation changed an ordered original candidate occurrence");
            ++occurrence;
        }
    }
    const auto &sample = fields.at("samples").at(0);
    const auto &counts = tensor(sample, "class_count");
    check(counts.tensor.state(0) == numeric::Rational{0}, "zero original count lost support");
    for (const auto *name : {"lower", "upper"}) {
        const auto &bound = tensor(sample, name);
        const auto coordinate_field = fields.at("target_coordinates");
        const auto &coordinates = coordinate_field.tensor();
        for (std::size_t t = 0; t < coordinates.domain().logical_size(); ++t)
            check(bound.tensor.state(coordinates.state(t)->get_num().get_ui()).has_value(), "target bound lost support");
    }
    const auto coordinate_record = fields.at("target_coordinates");
    const auto &coordinates = coordinate_record.tensor();
    const auto forged = [&](const numeric::Rational &coordinate) {
        auto declaration = fields.unbound_value();
        auto &members = std::get<NativeRecordValue::Object>(declaration.data).members;
        for (auto &[name, value] : members)
            if (name == "target_coordinates")
                value = NativeRecordValue(assay_detail::SourceGradeFields::tensor(
                    coordinates.domain().logical_size(),
                    [&](auto i) { return i == 0 ? coordinate : coordinates.state(i).value(); },
                    coordinates.support()));
        storage::TemporalEndpoint endpoint(source, fields.endpoint()->sequence);
        endpoint.carriers = fields.endpoint()->carriers;
        endpoint.tensors = fields.endpoint()->tensors;
        return NativeRecord::create_from_endpoint(std::move(declaration),
            new_native_record_namespace(), fields.raw_sources(), std::move(endpoint));
    };
    check(render_record_view(forged(numeric::Rational{0})).dump() == render_record_view(fields).dump(),
          "forged-coordinate fixture changed an unrelated field");
    for (const auto &coordinate : numeric::RationalVector{-1, numeric::Rational{3, 2},
             numeric::Rational{1, 2}, 1, 5, numeric::Rational{"1208925819614629174706176"}}) {
        bool refused = false;
        try { (void)render_record_view(forged(coordinate)); }
        catch (const std::invalid_argument &) { refused = true; }
        check(refused, "allocation renderer accepted an invalid native target coordinate");
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 5 && std::string(argv[1]) == "--reserved-pipeline") {
            const auto hardware = runtime::detect_hardware();
            runtime::CpuResourcePolicy policy;
            policy.worker_limit = hardware.allowed_cpus.size();
            policy.declared_buffer_byte_limit = 64 * 1024 * 1024;
            const auto resources = std::make_shared<CpuResources>(hardware, policy,
                [] { return runtime::detect_hardware(); });
            auto &coordinator = resources->coordinator();
            auto retained = coordinator.reserve(1, 4096);
            coordinator.resize(retained, 0, 4096);
            const auto available = coordinator.available_declared_buffer_bytes();
            const auto workspace = std::make_shared<const AssayWorkspace>(argv[2]);
            const auto result = execute_pipeline(workspace, Json::parse(read_text(argv[3])), argv[4], {}, resources);
            check(result.at("status").text() == "complete", "resident bytes blocked pipeline execution");
            check(coordinator.available_declared_buffer_bytes() == available,
                  "pipeline changed another owner's retained byte charge");
            check(coordinator.available_workers() == coordinator.capacity(),
                  "pipeline leaked its worker reservation");
            std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--consensus-artifact") {
            const auto publication = NativeRecord::read(argv[2]);
            const auto certificate = publication.at("result"), fields = certificate.at("result");
            check(certificate.at("schema").text() == "gadget.assay.result.native.v1" &&
                      fields.at("schema").text() == "gadget.sequence.consensus.fields.v1",
                  "typed consensus answer references a numerical display");
            const auto view = fields.at("view"), families = view.at("families");
            check(view.at("family_count").kind() == NativeRecordKind::Integer &&
                      view.at("family_count").integer() == families.size() && view.at("complete").boolean(),
                  "typed consensus answer lost complete native family count");
            const auto optimizer = NativeRecord::decode(fields.at("optimization_publication").bytes());
            check(optimizer.source_namespace() != fields.source_namespace(),
                  "typed consensus conflated biological and optimizer sources");
            check(!fields.at("source_bytes").bytes().empty() && fields.at("original_tokens").size() > 0,
                  "typed consensus lost original raw tokens");
            const auto &cost = tensor(fields, "edit_costs"), &pairs = tensor(fields, "paired_occurrence_map");
            check(cost.role == storage::TemporalTensorRole::Field && cost.output == 2 &&
                      cost.tensor.domain().logical_size() == fields.endpoint()->source.grade({1}).cell_count() &&
                      pairs.input == 1 && pairs.output == 1,
                  "typed consensus lost actual biological source grade fields");
            if (families.size()) {
                bool full_arity = false;
                const auto &relations = fields.endpoint()->source.grade({1});
                for (std::size_t r = 0; r < relations.cell_count(); ++r)
                    full_arity = full_arity || relations.offsets.at(r + 1) - relations.offsets.at(r) > 2;
                check(full_arity, "typed consensus projected full-arity biology into binary optimization");
                check(view.at("minimum_cost").kind() == NativeRecordKind::Rational,
                      "typed consensus lost native exact objective");
                for (std::size_t f = 0; f < families.size(); ++f) {
                    const auto family = families.at(f);
                    check(family.at("read_alignments").size() == fields.at("original_reads").size() &&
                              family.at("route_index").integer() < optimizer.at("routes").size(),
                          "typed consensus dropped an original read or route selector");
                }
            } else {
                check(view.at("status").text() == "no_admissible_consensus" &&
                          view.at("minimum_length").kind() == NativeRecordKind::Null &&
                          view.at("minimum_edits").kind() == NativeRecordKind::Null &&
                          view.at("minimum_cost").kind() == NativeRecordKind::Null,
                      "empty family evidence became an asserted consensus");
            }
            check(NativeRecord::decode(publication.encode()).encode() == publication.encode(),
                  "typed consensus complete parent roundtrip changed");
            std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
            return 0;
        }
        if (argc == 3 && (std::string(argv[1]) == "--planning" || std::string(argv[1]) == "--admission")) {
            const bool require_section = std::string(argv[1]) == "--planning";
            const auto planning = NativeRecord::read(argv[2]).at("result");
            check(planning.at("schema").text() == "gadget.scientific.fields.v1", "planning retained numerical admission displays");
            for (const auto *name : {"scientific_admission", "available_answers"}) {
                const auto reference = planning.at("view").at(name);
                check(reference.at("schema").text() == "gadget.scientific.publication.view.v1", "planning lost finite publication owner");
                const auto owner = NativeRecord::decode(reference.at("publication").bytes());
                if (!require_section && owner.at("status").text() == "untyped") {
                    check(!owner.has("native"), "untyped admission invented a finite constraint family");
                    continue;
                }
                const auto system = storage::temporal_finite_constraint_system<std::int64_t>(*owner.endpoint(), 1);
                const auto section = storage::temporal_finite_section<std::int64_t>(*owner.endpoint(), 1, 1);
                const bool valid = system.globally_valid(section);
                check(valid == owner.at("native").at("global_section_exists").boolean(),
                      "admission family differs from its actual native section validity");
                check(!require_section || valid, "planning section is not valid on its actual carrier");
            }
            std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--replay") {
            const auto original = NativeRecord::read(argv[2]);
            const auto workspace_path = fs::path(argv[2]).string() + ".replay-workspace.rfcxs";
            NativeRecord::decode(original.at("result").at("workspace_publication").bytes()).write(workspace_path);
            const AssayWorkspace embedded(workspace_path);
            const auto replayed = embedded.execute_record(original.at("request"));
            check(render_record_view(replayed.at("result")).dump() == render_record_view(original.at("result")).dump(),
                  "moved-input embedded replay differs");
            std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
            return 0;
        }
        if (argc != 4) throw std::invalid_argument("workspace, requests, output directory required");
        const AssayWorkspace workspace(argv[1]);
        const auto requests = read_record(argv[2]).native();
        const fs::path output(argv[3]); fs::create_directories(output);
        finite();
        for (const auto &name : requests.names()) {
            const auto request = requests.at(name);
            const auto result = workspace.execute_record(request);
            check(AssayWorkspace::supports_native_publication(request.at("operation").text()), "native publication predicate excludes producer");
            check(NativeRecord::decode(result.at("request_publication").bytes()).encode() == requests.encode(), "original whole request publication changed");
            const auto certificate = result.at("result"), fields = certificate.at("result");
            check(certificate.at("schema").text() == "gadget.assay.result.native.v1", "producer retained only a numerical display");
            check(NativeRecord::decode(certificate.at("workspace_publication").bytes()).encode() == workspace.native_publication().encode(), "original workspace owner changed");
            const auto operation = request.at("operation").text();
            if (operation == "relational_gluing") gluing(workspace, fields, request);
            else if (operation == "field_transfer") transfer(workspace, fields, request);
            else if (operation == "allocation_bounds") allocation(fields);
            else throw std::invalid_argument("unexpected producer");
            const auto path = output / (name + ".rfcxs"); result.write(path);
            const auto reopened = NativeRecord::read(path);
            check(reopened.encode() == result.encode(), "complete native archive changed on roundtrip");
            write_text(output / (name + ".json"), render_record_view(certificate).dump());
            const auto workspace_path = output / (name + ".workspace.rfcxs");
            NativeRecord::decode(reopened.at("result").at("workspace_publication").bytes()).write(workspace_path);
            const AssayWorkspace embedded(workspace_path);
            const auto replay = embedded.execute_record(reopened.at("request"));
            check(render_record_view(replay.at("result")).dump() == render_record_view(certificate).dump(), "embedded replay changed exact result");
        }
        std::cout << "{\"status\":\"PASS\",\"native_checks\":" << checks << "}\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
