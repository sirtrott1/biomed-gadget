#include "assay_field_detail.hpp"
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
Json chain_record(const inference::DerivedNativeGluedChain &input,
                  const core::RelationalComplex &source) {
    auto out = Json::object(), values = Json::array();
    const auto &chain = input.chain();
    for (const auto &[index, q] : chain.values()) {
        auto item = cell_record(storage::stable_cell_key(source, {chain.grade(), index}));
        item.set("index", cardinality(index)); item.set("value", Json::string(q.get_str()));
        values.push(item);
    }
    out.set("grade", Json::integer(chain.grade().value)); out.set("values", values);
    out.set("action", Json::string(input.input() ? "boundary" : "same_field_sections"));
    out.set("execution_steps", cardinality(input.provenance().execution_steps));
    out.set("preparation_steps", cardinality(input.provenance().preparation_steps));
    out.set("execution_policy", Json::string(input.provenance().execution_policy));
    return out;
}
} // namespace
Json relational_gluing(const AssayWorkspace &workspace, const Json &request,
                       std::shared_ptr<CpuResources> resources) {
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
    auto out = Json::object(), bindings = Json::object(), cell_maps = Json::array(),
        occurrence_maps = Json::array(), inventory = Json::array(), outputs = Json::array(),
        retained_cells = Json::array(), retained_occurrences = Json::array();
    bindings.set("common", common.binding); bindings.set("left", left.binding); bindings.set("right", right.binding);
    for (const auto &[from, to] : constructed->cells) {
        auto item = Json::object(); item.set("source", decimal(from.source));
        item.set("from", cell_record(from.cell)); item.set("to", cell_record(to)); cell_maps.push(item);
    }
    for (const auto &[from, to] : constructed->occurrences) {
        auto item = Json::object(); item.set("source", decimal(from.source));
        item.set("from", occurrence_record(from.occurrence)); item.set("to", occurrence_record(to));
        occurrence_maps.push(item);
    }
    // The native materialization supplies participant inclusions. Address each
    // common-interface witness through its already certified left image when
    // reporting provenance; mathematical commutation is checked by the native
    // prepared owner, rather than inferred here from labels or coefficients.
    for (const auto &image : manifest.embeddings[0].cells) {
        auto item = Json::object(); item.set("source", decimal(common.native.key()));
        item.set("from", cell_record(image.from));
        item.set("to", cell_record(constructed->cells.at({left.native.key(), image.to})));
        cell_maps.push(item);
    }
    for (const auto &image : manifest.embeddings[0].occurrences) {
        auto item = Json::object(); item.set("source", decimal(common.native.key()));
        item.set("from", occurrence_record(image.from));
        item.set("to", occurrence_record(constructed->occurrences.at({left.native.key(), image.to})));
        occurrence_maps.push(item);
    }
    for (std::size_t k = 0; k < constructed->source.grade_count(); ++k) {
        const Grade g{static_cast<std::uint16_t>(k)}; const auto &a = constructed->source.grade(g);
        auto record = Json::object(); record.set("grade", Json::integer(g.value));
        record.set("cells", cardinality(a.cell_count())); record.set("occurrences", cardinality(a.incidence_count()));
        inventory.push(record);
        for (std::size_t i = 0; i < a.cell_count(); ++i) {
            auto cell = cell_record(storage::stable_cell_key(constructed->source,
                {g, static_cast<std::uint32_t>(i)}));
            cell.set("index", cardinality(i));
            cell.set("scale", Json::string(a.cell_scales.at(i).get_str()));
            if (k > 0) cell.set("type", cardinality(a.relation_types.at(i)));
            retained_cells.push(cell);
        }
        for (std::size_t i = 0; i < a.incidence_count(); ++i) {
            auto occurrence = occurrence_record({g, a.occurrence_identities.at(i)});
            occurrence.set("parent", cell_record(storage::stable_cell_key(constructed->source,
                {g, a.upper_index.at(i)})));
            occurrence.set("participant", cell_record(storage::stable_cell_key(constructed->source,
                {Grade{static_cast<std::uint16_t>(k - 1)}, a.lower_index.at(i)})));
            occurrence.set("position", cardinality(i - a.offsets.at(a.upper_index.at(i))));
            occurrence.set("lattice_coefficient", Json::string(a.lattice_coefficient.at(i).get_str()));
            occurrence.set("coefficient", Json::string(constructed->source.incidence_coefficient(g, i).get_str()));
            auto aspects = Json::object();
            const auto &schema = constructed->source.signature().aspects();
            for (std::size_t j = 0; j < schema.size(); ++j) {
                const auto id = static_cast<composite::AspectId>(j);
                aspects.set(schema.at(id).name, Json::boolean(a.aspects.get(id, i)));
            }
            occurrence.set("aspects", aspects); retained_occurrences.push(occurrence);
        }
    }
    for (std::size_t i = 0; i < results.size(); ++i) {
        auto item = Json::object(), lineage = Json::array();
        item.set("declaration", Json::parse(sections.at(i).dump()));
        item.set("left_meaning", left.fields.at(text(sections.at(i).at("left_field"))));
        item.set("right_meaning", right.fields.at(text(sections.at(i).at("right_field"))));
        for (const auto &chain : results[i]) lineage.push(chain_record(chain, constructed->source));
        item.set("lineage", lineage); outputs.push(item);
    }
    auto run = Json::object(); run.set("owner", Json::string("rfc::inference::PreparedNativeGluing"));
    run.set("arithmetic", Json::string("exact_rational")); run.set("dimension_limit", cardinality(dimension));
    run.set("requested_workers", cardinality(wanted)); run.set("worker_slots", cardinality(execution.worker_slots));
    run.set("completed_tiles", cardinality(execution.completed_tiles));
    run.set("participant_publications_reused", Json::boolean(true));
    out.set("declaration", Json::parse(d.dump()));
    out.set("source_bindings", bindings); out.set("source_inventory", inventory);
    out.set("cell_maps", cell_maps); out.set("occurrence_maps", occurrence_maps);
    out.set("retained_cells", retained_cells); out.set("retained_occurrences", retained_occurrences);
    out.set("sections", outputs); out.set("boundary_policy", Json::string(boundary));
    out.set("constructed_source_status", Json::string("prospective")); out.set("execution", run);
    return out;
}
} // namespace gadget::assay_detail
