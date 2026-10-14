#include "assay_field_detail.hpp"
#include "assay_result_fields.hpp"
#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
#include <charconv>
#include <gadget/resources.hpp>
#include <gadget/native_gluing.hpp>
#include <rfc/inference/native_gluing.hpp>
#include <rfc/runtime/cpu_tiles.hpp>

namespace gadget::assay_detail {
namespace {
using namespace rfc;
using Grade = core::Grade;
using Cell = storage::StableCellKey;
using Occurrence = storage::StableOccurrenceKey;
constexpr auto addressable = std::numeric_limits<std::size_t>::max();
Grade grade(const Json &value) {
    return {static_cast<std::uint16_t>(count(value, std::numeric_limits<std::uint16_t>::max(),
                                           "gluing grade"))};
}
std::uint64_t sequence(const Json &value) {
    const auto word = value.text();
    std::uint64_t n{};
    const auto [end, error] = std::from_chars(word.data(), word.data() + word.size(), n);
    require(!word.empty() && (word == "0" || word.front() != '0') &&
                error == std::errc{} && end == word.data() + word.size(),
            "gluing availability requires a canonical unsigned decimal string");
    return n;
}
Json decimal(std::uint64_t n) { return Json::string(std::to_string(n)); }
Json cell_record(Cell c) {
    auto out = Json::object();
    out.set("grade", Json::integer(c.grade.value));
    out.set("key", decimal(c.key));
    return out;
}
Json occurrence_record(Occurrence o) {
    auto out = Json::object();
    out.set("grade", Json::integer(o.upper_grade.value));
    out.set("key", decimal(o.key.value));
    return out;
}
struct Source {
    std::shared_ptr<const ScientificStudy> owner;
    storage::NativeAtlasSource native;
    Json binding;
    runtime::WorkAvailability availability;
    std::map<std::string, Json> fields;
    std::map<std::uint64_t, const storage::TemporalTensor *> tensors;
    std::map<std::uint64_t, const storage::TemporalCarrier *> carriers;
    Source(const AssayWorkspace &workspace, const Json &record)
        : owner(workspace.study(text(record.at("assay")),
                                count(record.at("partition"), addressable, "gluing partition"))),
          native(owner->native_atlas_source(key(record.at("key")))),
          binding(Json::parse(record.dump())),
          availability{native.name_space(), sequence(record.at("available_sequence"))} {
        members(record, {"key", "assay", "partition", "available_sequence"});
        runtime::admit_work_availability(native.name_space(), native.sequence(), availability);
        require(native.endpoint() == owner->endpoint(), "gluing publication owner changed");
        binding.set("source_namespace", decimal(native.name_space()));
        binding.set("source_sequence", decimal(native.sequence()));
        binding.set("source_digest", Json::string(owner->digest()));
        // Compatibility output names the original byte artifacts. The native
        // complex and its accessions supply all cells, occurrences and values.
        const auto artifacts = owner->source_binding();
        binding.set("native", artifacts.at("native"));
        binding.set("declaration_sha256", artifacts.at("declaration_sha256"));
        binding.set("evidence_profile", Json::string("original-source-bytes-v1"));
        for (const auto &t : native.endpoint()->tensors) tensors.emplace(t.key, &t);
        for (const auto &c : native.endpoint()->carriers) carriers.emplace(c.key, &c);
    }
    Cell cell(Grade g, const Json &id) const {
        return owner->cell_identity(g, text(id));
    }
    std::pair<Grade, storage::PinnedNativeAccession> field(const Json &id,
                                                        const Json &section) {
        const auto name = text(id);
        if (!fields.contains(name)) fields.emplace(name, owner->field_declaration(name));
        const auto meaning = fields.at(name);
        require(meaning.at("value_type").text() == "quantitative",
                "gluing sections require declared quantitative scalar fields");
        for (const auto *member : {"quantity", "units", "channel"})
            require(meaning.at(member).text() == text(section.at(member)),
                    "gluing section field meanings differ");
        const auto tensor_key = key(meaning.at("key"));
        require(tensors.contains(tensor_key), "gluing field tensor missing");
        const auto &t = *tensors.at(tensor_key);
        require(t.output && carriers.contains(*t.output), "gluing field carrier missing");
        const auto &c = *carriers.at(*t.output);
        require(c.kind == storage::TemporalCarrierKind::SourceGrade &&
                    c.domain.extents().size() == 1 &&
                    t.tensor.domain().extents().size() == 1 &&
                    t.tensor.domain() == c.domain &&
                    c.domain.logical_size() == native.endpoint()->source.grade(c.grade).cell_count(),
                "gluing field must use the actual source grade carrier");
        // Native accession supplies the existing tensor as an authentic chain;
        // no host value vector or replacement publication is constructed.
        auto pin = owner->field_accession(name).access(
            storage::NativeAccessionRequest::chain(tensor_key, c.grade));
        return {c.grade, std::move(pin)};
    }
};
storage::NativeSourceEmbedding embedding(const Json &record, const Source &common,
                                          const Source &target) {
    members(record, {"key", "cells", "occurrences"});
    storage::NativeSourceEmbedding out;
    out.key = key(record.at("key")); out.from = common.native.key(); out.to = target.native.key();
    for (const auto &name : {"cells", "occurrences"}) source_array(record.at(name), name);
    const auto cells = record.at("cells"), occurrences = record.at("occurrences");
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto item = cells.at(i); members(item, {"grade", "from", "to"});
        const auto g = grade(item.at("grade"));
        out.cells.push_back({common.cell(g, item.at("from")), target.cell(g, item.at("to"))});
    }
    for (std::size_t i = 0; i < occurrences.size(); ++i) {
        const auto item = occurrences.at(i); members(item, {"grade", "from", "to"});
        const auto g = grade(item.at("grade"));
        out.occurrences.push_back({{g, core::OccurrenceKey{key(item.at("from"))}},
                                    {g, core::OccurrenceKey{key(item.at("to"))}}});
    }
    return out;
}
using N = NativeRecordValue;
N integer(std::uint64_t n) { return N(numeric::Integer(std::to_string(n))); }
N native_cell(Cell c) {
    return N(N::Object{{{"grade", integer(c.grade.value)}, {"key", N(std::to_string(c.key))}}});
}
N native_occurrence(Occurrence o) {
    return N(N::Object{{{"grade", integer(o.upper_grade.value)},
                        {"key", N(std::to_string(o.key.value))}}});
}
N archive(std::vector<RawByteSource> &sources, const NativeRecord &owner) {
    std::uint64_t id = 1;
    for (const auto &source : sources)
        id = std::max(id, machine::checked_size_add(source.key(), 1));
    auto bytes = owner.archive_source(id);
    const auto size = bytes.bytes().size();
    sources.push_back(std::move(bytes));
    return N(N::ByteSpan{id, 0, size, false});
}
N grade_carriers(storage::TemporalEndpoint &endpoint) {
    N::Array out;
    for (std::size_t k = 0; k < endpoint.source.grade_count(); ++k) {
        endpoint.carriers.push_back({k + 1, storage::TemporalCarrierKind::SourceGrade,
            core::Grade{static_cast<std::uint16_t>(k)}, {}, 0,
            numeric::IndexDomain::linear(endpoint.source.grade({static_cast<std::uint16_t>(k)}).cell_count())});
        out.elements.push_back(integer(k + 1));
    }
    return N(std::move(out));
}
} // namespace
Json relational_gluing(const AssayWorkspace &workspace, const Json &request,
                       std::shared_ptr<CpuResources> resources,
                       std::optional<AssayFieldsSource> *retained) {
    members(request, {"schema", "operation", "declaration", "sources", "embeddings", "sections",
                       "boundary_policy"}, {"workers"});
    const auto sources = request.at("sources"), maps = request.at("embeddings"), d = request.at("declaration");
    members(sources, {"common", "left", "right"}); members(maps, {"common_left", "common_right"});
    members(d, {"key", "revision", "pushout", "method", "evidence"});
    (void)text(d.at("method")); (void)text(d.at("evidence"));
    const auto boundary = request.at("boundary_policy").text();
    require(boundary == "none" || boundary == "to_grade_zero", "gluing boundary policy");
    Source common(workspace, sources.at("common")), left(workspace, sources.at("left")),
        right(workspace, sources.at("right"));
    require(std::set<std::uint64_t>{common.native.key(), left.native.key(), right.native.key()}.size() == 3,
            "gluing requires three distinct declared source identities");
    storage::NativeAtlasManifest manifest;
    manifest.key = key(d.at("key")); manifest.revision = key(d.at("revision"));
    manifest.embeddings = {embedding(maps.at("common_left"), common, left),
                           embedding(maps.at("common_right"), common, right)};
    manifest.pushouts.push_back({key(d.at("pushout")), manifest.embeddings[0].key,
                                 manifest.embeddings[1].key, {}, {}});
    auto limits = storage::NativeAtlasLimits{};
    limits.max_sources = 3; limits.max_definitions = 3;
    limits.max_records = limits.max_work = limits.max_dimension = limits.max_archive_bytes =
        limits.max_path_length = addressable;
    limits.max_integer_bytes = addressable / std::numeric_limits<unsigned char>::digits;
    limits.source_declarations = {addressable, addressable};
    // NativeAtlas uses these participants' original accession owners and limits.
    std::size_t dimension = 1;
    for (const auto *participant : {&common, &left, &right}) {
        const auto &c = participant->native.endpoint()->source;
        dimension = std::max(dimension, c.grade_count());
        dimension = std::max(dimension, c.signature().aspects().size());
        dimension = std::max(dimension, c.signature().type_count());
    }
    for (std::size_t k = 0; k < std::max(left.native.endpoint()->source.grade_count(),
                                      right.native.endpoint()->source.grade_count()); ++k) {
        const auto size = [k](const Source &s) {
            const auto &c = s.native.endpoint()->source;
            return k < c.grade_count() ? c.grade(Grade{static_cast<std::uint16_t>(k)}).cell_count() : 0;
        };
        dimension = std::max(dimension, machine::checked_size_add(size(left), size(right)));
    }
    const inference::NativeGluingAvailability available{common.availability, left.availability,
                                                       right.availability};
    auto plan = inference::PreparedNativeGluing::prepare(
        {storage::NativeAtlas({common.native, left.native, right.native}, manifest, limits),
         manifest.pushouts[0].key, dimension}, available,
        source_execution_lease("gadget-relational-gluing-v1"));
    const auto constructed = plan.constructed_source();
    const auto sections = request.at("sections"); source_array(sections, "gluing sections");
    std::vector<NativeGluingSection> inputs;
    std::set<std::string> ids; std::set<std::uint64_t> keys;
    for (std::size_t i = 0; i < sections.size(); ++i) {
        const auto s = sections.at(i);
        members(s, {"id", "key", "left_field", "right_field", "quantity", "units", "channel",
                     "method", "evidence"});
        require(ids.insert(text(s.at("id"))).second && keys.insert(key(s.at("key"))).second,
                "duplicate gluing section identity");
        (void)text(s.at("method")); (void)text(s.at("evidence"));
        auto [lg, lp] = left.field(s.at("left_field"), s);
        auto [rg, rp] = right.field(s.at("right_field"), s);
        require(lg == rg, "gluing section grades differ");
        inputs.push_back({lg, key(s.at("key")), std::move(lp), std::move(rp)});
    }
    const auto wanted = request.has("workers") ? count(request.at("workers"), addressable,
                                                       "gluing workers") : 1;
    require(wanted > 0, "gluing workers must be positive");
    auto native_result = execute_native_gluing(plan, std::move(inputs), available,
        boundary == "to_grade_zero", resources->reserve(wanted),
        source_execution_lease("gadget-glued-section-v1"));
    const auto &results = native_result.sections;
    const auto &execution = native_result.execution;
    storage::TemporalEndpoint endpoint(constructed->source);
    const auto carriers = grade_carriers(endpoint);
    SourceGradeFields fields(endpoint);
    N::Array native_sections, cell_maps, occurrence_maps;
    std::vector<RawByteSource> owners;
    N::Object participants;
    N::Object bindings;
    for (const auto &[name, participant] :
         std::vector<std::pair<std::string, const Source *>>{
             {"common", &common}, {"left", &left}, {"right", &right}}) {
        participants.members.emplace_back(name, archive(owners, participant->owner->native_publication()));
        bindings.members.emplace_back(name, import_record_view(participant->binding));
    }
    const auto cell_map = [&](std::uint64_t from, Cell cell, Cell image) {
        cell_maps.elements.emplace_back(N::Object{{{"source", N(std::to_string(from))},
            {"from", native_cell(cell)}, {"to", native_cell(image)}}});
    };
    const auto occurrence_map = [&](std::uint64_t from, Occurrence occurrence, Occurrence image) {
        occurrence_maps.elements.emplace_back(N::Object{{{"source", N(std::to_string(from))},
            {"from", native_occurrence(occurrence)}, {"to", native_occurrence(image)}}});
    };
    for (const auto &[from, to] : constructed->cells) cell_map(from.source, from.cell, to);
    for (const auto &[from, to] : constructed->occurrences)
        occurrence_map(from.source, from.occurrence, to);
    for (const auto &image : manifest.embeddings[0].cells)
        cell_map(common.native.key(), image.from,
                 constructed->cells.at({left.native.key(), image.to}));
    for (const auto &image : manifest.embeddings[0].occurrences)
        occurrence_map(common.native.key(), image.from,
                       constructed->occurrences.at({left.native.key(), image.to}));
    for (std::size_t i = 0; i < results.size(); ++i) {
        N::Array lineage;
        for (const auto &derived : results[i]) {
            const auto &chain = derived.chain();
            const auto size = endpoint.source.grade(chain.grade()).cell_count();
            std::vector<std::size_t> supported;
            supported.reserve(chain.values().size());
            for (const auto &[coordinate, value] : chain.values()) {
                (void)value; supported.push_back(coordinate);
            }
            const auto tensor = fields.column("gluing:" + std::to_string(i) + ":" +
                std::to_string(chain.grade().value), chain.grade().value, size,
                [&](auto j) { return chain.values().at(j); },
                numeric::StructuralSupport::sparse(size, std::move(supported)));
            lineage.elements.emplace_back(N::Object{{{"grade", integer(chain.grade().value)},
                {"field", tensor},
                {"action", N(std::string(derived.input() ? "boundary" : "same_field_sections"))},
                {"execution_steps", integer(derived.provenance().execution_steps)},
                {"preparation_steps", integer(derived.provenance().preparation_steps)},
                {"execution_policy", N(derived.provenance().execution_policy)}}});
        }
        native_sections.elements.emplace_back(N::Object{{
            {"declaration", import_record_view(sections.at(i))},
            {"left_meaning", import_record_view(left.fields.at(text(sections.at(i).at("left_field"))))},
            {"right_meaning", import_record_view(right.fields.at(text(sections.at(i).at("right_field"))))},
            {"lineage", N(std::move(lineage))}}});
    }
    N declaration(N::Object{{{"schema", N(std::string("gadget.assays.gluing.fields.v1"))},
        {"request", import_record_view(request)}, {"source_bindings", N(std::move(bindings))},
        {"participant_publications", N(std::move(participants))}, {"carriers", carriers},
        {"cell_maps", N(std::move(cell_maps))}, {"occurrence_maps", N(std::move(occurrence_maps))},
        {"sections", N(std::move(native_sections))},
        {"execution", N(N::Object{{{"owner", N(std::string("rfc::inference::PreparedNativeGluing"))},
            {"arithmetic", N(std::string("exact_rational"))}, {"dimension_limit", integer(dimension)},
            {"requested_workers", integer(wanted)}, {"worker_slots", integer(execution.worker_slots)},
            {"completed_tiles", integer(execution.completed_tiles)},
            {"participant_publications_reused", N(true)}}})}}});
    if (retained) {
        retained->emplace(AssayFieldsSource{std::move(endpoint), std::move(declaration), std::move(owners)});
        return Json::object();
    }
    return gluing_fields_view(NativeRecord::create_from_endpoint(std::move(declaration),
        new_native_record_namespace(), std::move(owners), std::move(endpoint)));
}
Json gluing_fields_view(const NativeRecord &record) {
    const auto &source = record.endpoint()->source;
    auto out = Json::object(), inventory = Json::array(), retained_cells = Json::array(),
        retained_occurrences = Json::array(), outputs = Json::array();
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        const Grade g{static_cast<std::uint16_t>(k)}; const auto &a = source.grade(g);
        auto record = Json::object(); record.set("grade", Json::integer(g.value));
        record.set("cells", cardinality(a.cell_count())); record.set("occurrences", cardinality(a.incidence_count()));
        inventory.push(record);
        for (std::size_t i = 0; i < a.cell_count(); ++i) {
            auto cell = cell_record(storage::stable_cell_key(source,
                {g, static_cast<std::uint32_t>(i)}));
            cell.set("index", cardinality(i));
            cell.set("scale", Json::string(a.cell_scales.at(i).get_str()));
            if (k > 0) cell.set("type", cardinality(a.relation_types.at(i)));
            retained_cells.push(cell);
        }
        for (std::size_t i = 0; i < a.incidence_count(); ++i) {
            auto occurrence = occurrence_record({g, a.occurrence_identities.at(i)});
            occurrence.set("parent", cell_record(storage::stable_cell_key(source,
                {g, a.upper_index.at(i)})));
            occurrence.set("participant", cell_record(storage::stable_cell_key(source,
                {Grade{static_cast<std::uint16_t>(k - 1)}, a.lower_index.at(i)})));
            occurrence.set("position", cardinality(i - a.offsets.at(a.upper_index.at(i))));
            occurrence.set("lattice_coefficient", Json::string(a.lattice_coefficient.at(i).get_str()));
            occurrence.set("coefficient", Json::string(source.incidence_coefficient(g, i).get_str()));
            auto aspects = Json::object();
            const auto &schema = source.signature().aspects();
            for (std::size_t j = 0; j < schema.size(); ++j) {
                const auto id = static_cast<composite::AspectId>(j);
                aspects.set(schema.at(id).name, Json::boolean(a.aspects.get(id, i)));
            }
            occurrence.set("aspects", aspects); retained_occurrences.push(occurrence);
        }
    }
    const auto sections = record.at("sections");
    for (std::size_t i = 0; i < sections.size(); ++i) {
        const auto section = sections.at(i);
        auto item = Json::object(), lineage = Json::array();
        for (const auto *name : {"declaration", "left_meaning", "right_meaning"})
            item.set(name, render_record_view(section.at(name)));
        for (std::size_t j = 0; j < section.at("lineage").size(); ++j) {
            const auto chain = section.at("lineage").at(j);
            const auto field_key = chain.at("field");
            const auto &field = retained_tensor(record, field_key);
            const auto carrier = std::find_if(record.endpoint()->carriers.begin(),
                record.endpoint()->carriers.end(), [&](const auto &c) { return c.key == field.output; });
            require(carrier != record.endpoint()->carriers.end() &&
                        carrier->kind == storage::TemporalCarrierKind::SourceGrade &&
                        numeric::Integer(carrier->grade.value) == chain.at("grade").integer() &&
                        field.tensor.domain() == carrier->domain,
                    "glued field must retain its actual constructed grade carrier");
            auto row = Json::object(), values = Json::array();
            row.set("grade", render_record_view(chain.at("grade")));
            field.tensor.support().for_each_active([&](std::size_t index) {
                auto value = cell_record(storage::stable_cell_key(source,
                    {carrier->grade, static_cast<std::uint32_t>(index)}));
                value.set("index", cardinality(index));
                value.set("value", Json::string(field.tensor.state(index).value().get_str()));
                values.push(value);
            });
            row.set("values", values);
            for (const auto *name : {"action", "execution_steps", "preparation_steps", "execution_policy"})
                row.set(name, render_record_view(chain.at(name)));
            lineage.push(row);
        }
        item.set("lineage", lineage); outputs.push(item);
    }
    out.set("declaration", render_record_view(record.at("request").at("declaration")));
    out.set("source_bindings", render_record_view(record.at("source_bindings")));
    out.set("source_inventory", inventory);
    out.set("cell_maps", render_record_view(record.at("cell_maps")));
    out.set("occurrence_maps", render_record_view(record.at("occurrence_maps")));
    out.set("retained_cells", retained_cells); out.set("retained_occurrences", retained_occurrences);
    out.set("sections", outputs);
    out.set("boundary_policy", render_record_view(record.at("request").at("boundary_policy")));
    out.set("constructed_source_status", Json::string("prospective"));
    out.set("execution", render_record_view(record.at("execution")));
    return out;
}
Json field_transfer(const AssayWorkspace &workspace, const Json &request,
                    std::shared_ptr<CpuResources> resources,
                    std::optional<AssayFieldsSource> *retained) {
    members(request, {"schema", "operation", "left", "right", "map", "metric"});
    require(request.at("metric").text() == "identity", "field transfer metric profile");
    const auto load = [&](const Json &record) {
        members(record, {"assay", "partition", "field"});
        auto input = load_assay_field(workspace.study(text(record.at("assay")),
            count(record.at("partition"), addressable, "field transfer partition")),
            text(record.at("field")));
        require(!input.values.empty(), "empty field transfer input");
        require(std::all_of(input.supported.begin(), input.supported.end(), [](auto v) { return v; }),
                "field transfer requires complete support; no implicit missing extension");
        return input;
    };
    const auto left = load(request.at("left")), right = load(request.at("right"));
    const auto map = request.at("map");
    members(map, {"domain_quantity", "domain_units", "codomain_quantity", "codomain_units",
                  "method", "evidence", "entries"});
    require(map.at("domain_quantity").text() == right.meaning.at("quantity").text() &&
                map.at("domain_units").text() == right.meaning.at("units").text() &&
                map.at("codomain_quantity").text() == left.meaning.at("quantity").text() &&
                map.at("codomain_units").text() == left.meaning.at("units").text(),
            "field transfer domain/codomain meanings differ from measurements");
    (void)text(map.at("method")); (void)text(map.at("evidence"), 16384);
    std::map<std::string, std::size_t> li, ri;
    for (std::size_t i = 0; i < left.ids.size(); ++i) li.emplace(left.ids.at(i).text(), i);
    for (std::size_t i = 0; i < right.ids.size(); ++i) ri.emplace(right.ids.at(i).text(), i);
    numeric::ExactMatrix matrix(left.values.size(), right.values.size());
    std::set<std::pair<std::size_t, std::size_t>> seen;
    for (std::size_t i = 0; i < map.at("entries").size(); ++i) {
        const auto e = map.at("entries").at(i);
        members(e, {"left_id", "right_id", "coefficient"});
        const auto l = text(e.at("left_id")), r = text(e.at("right_id"));
        require(li.contains(l) && ri.contains(r) && seen.emplace(li.at(l), ri.at(r)).second,
                "field transfer coordinate identity or duplicate entry");
        matrix.set(li.at(l), ri.at(r), exact_rational(text(e.at("coefficient"))));
    }
    auto reservation = resources->reserve(1);
    const auto transported = matrix.apply(right.values), residual = numeric::sub(left.values, transported);
    const auto original_carrier = [&](const AssayField &field) {
        const auto &endpoint = *field.owner->endpoint();
        const auto id = key(field.meaning.at("key"));
        const auto tensor = std::find_if(endpoint.tensors.begin(), endpoint.tensors.end(),
            [&](const auto &t) { return t.key == id; });
        require(tensor != endpoint.tensors.end() && tensor->output.has_value(),
                "field transfer original carrier missing");
        return *tensor->output;
    };
    const auto left_carrier = original_carrier(left), right_original_carrier = original_carrier(right);
    storage::TemporalEndpoint endpoint(left.owner->endpoint()->source, left.owner->endpoint()->sequence);
    endpoint.carriers = left.owner->endpoint()->carriers;
    SourceGradeFields fields(endpoint);
    // The independently admitted right publication has its own source identity.
    // Its explicit input axis is bound below to the original carrier and keys.
    const auto right_carrier = fields.explicit_carrier(right.values.size());
    const auto scalar = fields.explicit_carrier(1);
    std::vector<RawByteSource> owners;
    auto left_owner = archive(owners, left.owner->native_publication());
    auto right_owner = right.owner->endpoint() == left.owner->endpoint()
        ? left_owner : archive(owners, right.owner->native_publication());
    const auto binding = [&](const AssayField &field, std::uint64_t carrier, N owner) {
        return N(N::Object{{{"publication", std::move(owner)},
            {"source_namespace", integer(field.owner->native_publication().source_namespace())},
            {"source_sequence", integer(field.owner->endpoint()->sequence)},
            {"source_carrier", integer(carrier)}, {"source_field", integer(key(field.meaning.at("key")))},
            {"ids", import_record_view(field.ids)}, {"keys", import_record_view(field.keys)}}});
    };
    const auto rank = matrix.rank();
    N declaration(N::Object{{{"schema", N(std::string("gadget.assays.transfer.fields.v1"))},
        {"request", import_record_view(request)},
        {"left_binding", binding(left, left_carrier, std::move(left_owner))},
        {"right_binding", binding(right, right_original_carrier, std::move(right_owner))},
        {"left", fields.vector("transfer.left", left_carrier, left.values)},
        {"right", fields.vector("transfer.right", right_carrier, right.values)},
        {"map_field", fields.map("transfer.map", right_carrier, left_carrier, matrix)},
        {"metric_diagonal", fields.vector("transfer.identity_metric", left_carrier, Vector(left.values.size(), 1))},
        {"transported_right", fields.vector("transfer.transported_right", left_carrier, transported)},
        {"residual", fields.vector("transfer.residual", left_carrier, residual)},
        {"left_carrier_grade", integer(left.grade)}, {"right_carrier_grade", integer(right.grade)},
        {"map_rank", integer(rank)}, {"map_kernel_dimension", integer(matrix.cols() - rank)},
        {"residual_quadrance", fields.vector("transfer.residual_quadrance", scalar,
            Vector{numeric::dot(residual, residual)})},
        {"pairing", fields.vector("transfer.pairing", scalar, Vector{numeric::dot(left.values, transported)})},
        {"left_source_digest", N(left.owner->digest())}, {"right_source_digest", N(right.owner->digest())},
        {"left_meaning", import_record_view(left.meaning)}, {"right_meaning", import_record_view(right.meaning)}}});
    if (retained) {
        retained->emplace(AssayFieldsSource{std::move(endpoint), std::move(declaration), std::move(owners)});
        return Json::object();
    }
    return transfer_fields_view(NativeRecord::create_from_endpoint(std::move(declaration),
        new_native_record_namespace(), std::move(owners), std::move(endpoint)));
}
Json transfer_fields_view(const NativeRecord &record) {
    auto out = Json::object();
    for (const auto *name : {"left", "right", "transported_right", "residual"}) {
        auto values = Json::array();
        const auto &tensor = retained_tensor(record, record.at(name)).tensor;
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i)
            values.push(Json::string(tensor.state(i).value().get_str()));
        out.set(name, values);
    }
    out.set("left_ids", render_record_view(record.at("left_binding").at("ids")));
    out.set("right_ids", render_record_view(record.at("right_binding").at("ids")));
    for (const auto *name : {"left_carrier_grade", "right_carrier_grade", "map_rank", "map_kernel_dimension"})
        out.set(name, render_record_view(record.at(name)));
    for (const auto *name : {"residual_quadrance", "pairing"})
        out.set(name, Json::string(retained_tensor(record, record.at(name)).tensor.state(0).value().get_str()));
    for (const auto *name : {"left_source_digest", "right_source_digest", "left_meaning", "right_meaning"})
        out.set(name, render_record_view(record.at(name)));
    out.set("map", render_record_view(record.at("request").at("map")));
    out.set("scope", Json::string("Full field comparison under the declared measurement map and metric; the "
                                 "map's biological validity is supplied by its evidence"));
    return out;
}
Json allocation_bounds_fields_view(const NativeRecord &record) {
    auto out = Json::object(), samples = Json::array(), contrasts = Json::array();
    const auto targets = record.at("targets");
    const auto coordinate_field = record.at("target_coordinates");
    const auto &coordinates = coordinate_field.tensor();
    const auto bounds = [&](const NativeRecord &sample, bool contrast) {
        const auto &lower = retained_tensor(record, sample.at("lower")).tensor;
        const auto &upper = retained_tensor(record, sample.at("upper")).tensor;
        auto rows = Json::array();
        for (std::size_t t = 0; t < targets.size(); ++t) {
            const auto coordinate = coordinates.state(t).value().get_num();
            require(mpz_fits_ulong_p(coordinate.get_mpz_t()), "allocation target coordinate addressability");
            const auto lo = lower.state(coordinate.get_ui()).value(), hi = upper.state(coordinate.get_ui()).value();
            auto row = Json::object();
            row.set("id", render_record_view(targets.at(t).at("id")));
            row.set("lower", Json::string(lo.get_str())); row.set("upper", Json::string(hi.get_str()));
            if (contrast) row.set("direction", Json::string(lo > 0 ? "positive_for_every_allocation" :
                hi < 0 ? "negative_for_every_allocation" : lo == 0 && hi == 0 ? "zero_for_every_allocation" : "unresolved"));
            else row.set("identifiable", Json::boolean(lo == hi));
            rows.push(row);
        }
        return rows;
    };
    for (std::size_t s = 0; s < record.at("samples").size(); ++s) {
        const auto sample = record.at("samples").at(s);
        auto row = Json::object(); row.set("id", render_record_view(sample.at("id")));
        row.set("total", Json::string(sample.at("total").rational().get_str()));
        row.set("bounds", bounds(sample, false)); samples.push(row);
    }
    for (std::size_t c = 0; c < record.at("contrasts").size(); ++c) {
        const auto contrast = record.at("contrasts").at(c);
        auto row = Json::object();
        for (const auto *name : {"id", "left", "right"}) row.set(name, render_record_view(contrast.at(name)));
        row.set("bounds", bounds(contrast, true)); contrasts.push(row);
    }
    out.set("samples", samples); out.set("contrasts", contrasts);
    const auto summary = allocation_fields_summary_view(record);
    json_object_object_foreach(summary.raw(), name, value) {
        (void)value; out.set(name, summary.at(name));
    }
    return out;
}
} // namespace gadget::assay_detail
