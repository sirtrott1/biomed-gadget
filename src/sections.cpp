#include "scientific_detail.hpp"
#include "scientific_projection.hpp"
#include "scientific_sections.hpp"
#include <gadget/sections.hpp>
#include <rfc/local/cellular_sheaf.hpp>
#include <rfc/local/fiber_system.hpp>
#include <rfc/storage/temporal_endpoint.hpp>
#include <rfc/topology/sheaf_cohomology.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc;
using M = numeric::ExactMatrix;
Rational scalar(const Json &v) { return exact_rational(text(v, source_byte_limit())); }
M declared_matrix(const Json &v) {
    require(json_object_is_type(v.raw(), json_type_array) && v.size(),
            "local matrix requires rows");
    const auto cols = v.at(0).size();
    require(cols, "local matrix requires columns");
    M m(v.size(), cols);
    for (std::size_t i = 0; i < v.size(); ++i) {
        require(v.at(i).size() == cols, "ragged local map");
        for (std::size_t j = 0; j < cols; ++j)
            m.set(i, j, scalar(v.at(i).at(j)));
    }
    return m;
}
} // namespace
scientific_detail::LocalSectionAnalysis scientific_detail::analyze_local_section_fields(
    const NativeStudy &study, const Json &query, ResultFields *publication,
    const std::map<std::string, Vector> *operands, bool validate_only, bool include_cochain_maps) {
    members(query,
            {"operation", "context", "construction", "stalks", "maps", "sections", "metric"});
    const auto construction = text(query.at("construction"));
    require(construction == "cellular_sheaf" || construction == "boundary_fibers",
            "local construction");
    require(query.at("metric").text() == "identity", "local section metric profile");
    const auto &source = study.source();
    storage::StableSourceIndex index(source);
    local::LinearStalkBundle stalks(source);
    std::set<storage::StableCellKey> declared;
    const auto records = query.at("stalks");
    require(records.size(), "local stalks must be nonempty");
    std::size_t cell_count = 0;
    for (std::size_t g = 0; g < source.grade_count(); ++g)
        cell_count = machine::checked_size_add(
            cell_count, source.grade(core::Grade{static_cast<std::uint16_t>(g)}).cell_count());
    require(records.size() == cell_count, "every context cell requires a declared stalk");
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto item = records.at(i);
        members(item, {"grade", "key", "axes"});
        storage::StableCellKey k{{static_cast<std::uint16_t>(count(
                                     item.at("grade"), source.grade_count() - 1, "stalk grade"))},
                                 key(item.at("key"))};
        const auto cell = index.cell(k);
        require(cell.has_value() && declared.insert(k).second, "unknown or duplicate stalk");
        const auto axes = item.at("axes");
        require(axes.size(), "stalk requires axes");
        std::set<std::string> names;
        for (std::size_t a = 0; a < axes.size(); ++a)
            require(names.insert(text(axes.at(a), source_byte_limit())).second,
                    "duplicate stalk axis");
        stalks.set_dimension(*cell, axes.size());
    }
    local::LinearRelationalFiberSystem fibers(std::move(stalks));
    local::LinearCellularSheaf sheaf(fibers);
    std::set<storage::StableOccurrenceKey> flags;
    std::size_t flag_count = 0;
    for (std::size_t g = 1; g < source.grade_count(); ++g)
        flag_count = machine::checked_size_add(
            flag_count, source.grade(core::Grade{static_cast<std::uint16_t>(g)}).incidence_count());
    const auto maps = query.at("maps");
    require(maps.size() == flag_count, "every retained flag requires an explicit map");
    for (std::size_t i = 0; i < maps.size(); ++i) {
        const auto item = maps.at(i);
        members(item, {"grade", "key", "matrix"});
        storage::StableOccurrenceKey k{
            {static_cast<std::uint16_t>(
                count(item.at("grade"), source.grade_count() - 1, "flag grade"))},
            {key(item.at("key"))}};
        const auto slot = index.occurrence(k);
        require(slot.has_value() && flags.insert(k).second, "unknown or duplicate flag map");
        auto map = declared_matrix(item.at("matrix"));
        if (construction == "cellular_sheaf")
            sheaf.set_restriction(k.upper_grade, *slot, std::move(map));
        else
            fibers.set_transport(k.upper_grade, *slot, std::move(map));
    }
    if (construction == "cellular_sheaf")
        require(sheaf.is_cochain_complex(), "restriction two-step defect is nonzero");
    const auto sections = query.at("sections");
    require(sections.size(), "sections must be nonempty");
    std::set<std::string> ids;
    std::vector<Vector> values;
    for (std::size_t i = 0; i < sections.size(); ++i) {
        const auto item = sections.at(i);
        members(item, {"id", "grade", "quantity", "units"}, {"values", "fields", "axes"});
        require(ids.insert(text(item.at("id"))).second, "duplicate section identity");
        (void)text(item.at("quantity"));
        (void)text(item.at("units"));
        const auto grade = core::Grade{static_cast<std::uint16_t>(
            count(item.at("grade"), source.grade_count() - 1, "section grade"))};
        Vector v;
        if (operands && item.has("fields")) {
            require(operands->contains(item.at("id").text()) && !item.has("values"),
                    "source section requires its admitted native operand");
            v = operands->at(item.at("id").text());
        } else {
            require(item.has("values"), "declared section values absent");
            for (std::size_t j = 0; j < item.at("values").size(); ++j)
                v.push_back(scalar(item.at("values").at(j)));
        }
        require(v.size() == fibers.stalks().total_dimension(grade), "section direct-sum dimension");
        values.push_back(std::move(v));
    }
    if (validate_only) return {};
    require(publication != nullptr, "local sections require native result ownership");
    using Value = ResultValue;
    using Publication = ResultFields;
    const auto matrix = [](const M &value, std::uint64_t input, std::uint64_t output) {
        auto result = Publication::stage_matrix(value, input, output);
        result.set("display", Value::string("dense_matrix"));
        return result;
    };
    LocalSectionAnalysis analysis;
    std::map<storage::StableCellKey, std::uint64_t> stalk_axes;
    std::vector<std::vector<std::size_t>> offsets;
    auto binding = Value::object(), direct_sums = Value::array();
    binding.set("kind", Value::string("local_section_stalks"));
    binding.set("construction", query.at("construction"));
    for (std::size_t g = 0; g < source.grade_count(); ++g) {
        const auto grade = core::Grade{static_cast<std::uint16_t>(g)};
        offsets.push_back(fibers.stalks().offsets(grade));
        const auto axis = publication->axis(offsets.back().back());
        analysis.grade_carriers.push_back(axis);
        auto sum = Value::object();
        sum.set("grade", Value::count(g));
        sum.set("carrier", Value::count(axis));
        Vector positions;
        for (auto offset : offsets.back()) positions.emplace_back(std::to_string(offset));
        sum.set("offsets", Publication::stage_vector(positions, publication->axis(positions.size()), "indices"));
        sum.set("metric", Publication::stage_matrix(M::identity(offsets.back().back()), axis, axis,
            storage::TemporalTensorRole::Metric));
        direct_sums.push(std::move(sum));
    }
    const auto empty_axis = publication->axis(0);
    auto declared_stalks = Value::array();
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto item = records.at(i);
        const storage::StableCellKey key_{{static_cast<std::uint16_t>(item.at("grade").integer())}, key(item.at("key"))};
        const auto cell = *index.cell(key_);
        const auto size = fibers.stalks().dimension(cell), axis = publication->stalk(key_, size);
        stalk_axes.emplace(key_, axis);
        auto declared_stalk = Value(item);
        declared_stalk.set("carrier", Value::count(axis));
        M inclusion(offsets[cell.grade.value].back(), size);
        for (std::size_t a = 0; a < size; ++a) inclusion.set(offsets[cell.grade.value][cell.index] + a, a, Rational{1});
        declared_stalk.set("direct_sum_inclusion", Publication::stage_matrix(inclusion, axis, analysis.grade_carriers[cell.grade.value]));
        declared_stalks.push(std::move(declared_stalk));
    }
    binding.set("direct_sums", std::move(direct_sums));
    binding.set("stalks", std::move(declared_stalks));
    publication->bind(std::move(binding));
    auto native_maps = Value::array();
    for (std::size_t i = 0; i < maps.size(); ++i) {
        const auto item = maps.at(i);
        const core::Grade grade{static_cast<std::uint16_t>(item.at("grade").integer())};
        const auto slot = *index.occurrence({grade, {key(item.at("key"))}});
        const auto &arena = source.grade(grade);
        const auto upper = storage::stable_cell_key(source, {grade, arena.upper_index[slot]});
        const auto lower = storage::stable_cell_key(source, {{static_cast<std::uint16_t>(grade.value - 1u)}, arena.lower_index[slot]});
        const bool upward = construction == "cellular_sheaf";
        auto record = Value(item);
        auto map = Publication::stage_matrix(upward ? sheaf.restriction(grade, slot) : fibers.flag_map(grade, slot).transport,
            stalk_axes.at(upward ? lower : upper), stalk_axes.at(upward ? upper : lower));
        map.set("display", Value::string("matrix_values"));
        map.set("occurrence_grade", Value::count(grade.value));
        map.set("occurrence_key", Value::count(key(item.at("key"))));
        map.set("interpretation", Value::string(upward ? "declared_cellular_sheaf_restriction" : "declared_realized_boundary_transport"));
        record.set("matrix", std::move(map));
        native_maps.push(std::move(record));
    }
    if (construction == "cellular_sheaf")
        for (std::size_t g = 0; g < source.grade_count(); ++g)
            analysis.coboundaries.push_back(sheaf.coboundary(core::Grade{static_cast<std::uint16_t>(g)}));
    auto out = Value::object(), results = Value::array(), groups = Value::array();
    std::size_t incompatible = 0;
    if (construction == "cellular_sheaf") {
        for (const auto &group : topology::sheaf_cohomology(sheaf)) {
            auto r = Value::object();
            r.set("grade", Value::count(group.grade.value));
            r.set("dimension", Value::count(group.dimension));
            const auto axis = analysis.grade_carriers.at(group.grade.value);
            const auto frame = [&](const M &value) {
                return matrix(value, publication->axis(value.cols()), axis);
            };
            r.set("cocycles", frame(group.cocycle_frame));
            r.set("coboundaries", frame(group.coboundary_frame));
            r.set("representatives", frame(group.representative_frame));
            groups.push(r);
        }
    }
    std::map<std::uint16_t, std::pair<M, std::unique_ptr<PreparedProjection>>> prepared;
    for (std::size_t i = 0; i < sections.size(); ++i) {
        auto r = Value(sections.at(i));
        r.set("input_origin",
              Value::string(sections.at(i).has("fields") ? "admitted source fields"
                                                        : "declared query operand"));
        const auto grade =
            core::Grade{static_cast<std::uint16_t>(sections.at(i).at("grade").integer())};
        const auto &v = values[i];
        const auto axis = analysis.grade_carriers.at(grade.value);
        r.set("values", Publication::stage_vector(v, axis));
        if (construction == "cellular_sheaf") {
            if (!prepared.contains(grade.value)) {
                auto delta = analysis.coboundaries.at(grade.value);
                auto projection = std::make_unique<PreparedProjection>(delta.kernel_frame());
                prepared.emplace(grade.value,
                                 std::make_pair(std::move(delta), std::move(projection)));
            }
            const auto &delta = prepared.at(grade.value).first;
            const auto glued = prepared.at(grade.value).second->apply(v);
            const auto residual = numeric::sub(v, glued), defect = delta.apply(v);
            require(numeric::add(glued, residual) == v && numeric::is_zero(delta.apply(glued)) &&
                        numeric::dot(glued, residual) == 0,
                    "section gluing reconstruction");
            const bool compatible = numeric::is_zero(defect);
            if (!compatible)
                ++incompatible;
            r.set("compatible", Value::boolean(compatible));
            analysis.glued.push_back(glued);
            r.set("glued", Publication::stage_vector(glued, axis));
            r.set("residual", Publication::stage_vector(residual, axis));
            const auto target = grade.value + 1u < analysis.grade_carriers.size()
                ? analysis.grade_carriers.at(grade.value + 1u) : empty_axis;
            r.set("coboundary_defect", Publication::stage_vector(defect, target));
            r.set("coboundary", matrix(delta, axis, target));
            r.set("residual_quadrance", publication->stage_scalar(numeric::dot(residual, residual)));
        } else {
            require(grade.value > 0, "fiber relation sections require a positive grade");
            auto branches = Value::array();
            const auto offsets = fibers.stalks().offsets(grade);
            for (std::size_t cell = 0; cell < source.grade(grade).cell_count(); ++cell) {
                Vector local_value(v.begin() + offsets[cell], v.begin() + offsets[cell + 1]);
                const auto section = fibers.realize_relation(grade, cell, local_value);
                require(section.relation.grade == grade && section.relation.index == cell,
                        "relation section source identity");
                const auto upper_key = storage::stable_cell_key(source, section.relation);
                auto branch = Value::object();
                branch.set("key", Value::string(std::to_string(upper_key.key)));
                branch.set("upper", Publication::stage_vector(local_value, stalk_axes.at(upper_key)));
                auto legs = Value::array();
                for (const auto &leg : section.boundary_values) {
                    auto item = Value::object();
                    item.set("slot", Value::count(leg.incidence_slot));
                    item.set("lower_key", Value::string(std::to_string(
                                              storage::stable_cell_key(source, leg.lower).key)));
                    item.set("value", Publication::stage_vector(leg.value, stalk_axes.at(storage::stable_cell_key(source, leg.lower))));
                    legs.push(item);
                }
                branch.set("legs", legs);
                branches.push(branch);
            }
            r.set("relations", branches);
        }
        results.push(r);
    }
    out.set("construction", query.at("construction"));
    if (include_cochain_maps && construction == "cellular_sheaf") {
        auto maps_out = Value::array();
        for (std::size_t k = 0; k < source.grade_count(); ++k)
            maps_out.push(matrix(analysis.coboundaries[k], analysis.grade_carriers[k],
                k + 1u < analysis.grade_carriers.size() ? analysis.grade_carriers[k + 1u] : empty_axis));
        out.set("cochain_maps", maps_out);
    }
    out.set("stalks", records);
    out.set("maps", std::move(native_maps));
    out.set("cohomology", groups);
    out.set("sections", results);
    out.set("incompatible_sections", Value::count(incompatible));
    out.set("reconstruction_exact", Value::boolean(true));
    analysis.view = std::move(out);
    return analysis;
}
Json analyze_local_sections(const NativeStudy &study, const Json &query, bool validate_only,
                            bool include_cochain_maps) {
    if (validate_only) {
        (void)scientific_detail::analyze_local_section_fields(study, query, nullptr, nullptr, true);
        return Json{};
    }
    storage::TemporalEndpoint endpoint(study.source());
    scientific_detail::ResultFields publication(NativeRecord::create_from_endpoint(
        NativeRecordValue{}, new_native_record_namespace(), {}, std::move(endpoint)));
    auto analysis = scientific_detail::analyze_local_section_fields(study, query, &publication,
        nullptr, false, include_cochain_maps);
    return scientific_detail::result_fields_view(publication.finish(std::move(analysis.view)));
}
} // namespace gadget
