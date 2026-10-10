#include "scientific_detail.hpp"
#include "scientific_projection.hpp"
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
Json vector(const Vector &v) {
    auto out = Json::array();
    for (auto q : v) {
        q.canonicalize();
        out.push(Json::string(q.get_str()));
    }
    return out;
}
Json matrix(const M &m) {
    auto out = Json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        Vector row(m.cols());
        for (std::size_t j = 0; j < m.cols(); ++j)
            row[j] = m.get(i, j);
        out.push(vector(row));
    }
    auto result = Json::object();
    result.set("rows", cardinality(m.rows()));
    result.set("columns", cardinality(m.cols()));
    result.set("values", out);
    return result;
}
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
Json analyze_local_sections(const NativeStudy &study, const Json &query, bool validate_only,
                            bool include_cochain_maps) {
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
        members(item, {"id", "grade", "values", "quantity", "units"}, {"fields", "axes"});
        require(ids.insert(text(item.at("id"))).second, "duplicate section identity");
        (void)text(item.at("quantity"));
        (void)text(item.at("units"));
        const auto grade = core::Grade{static_cast<std::uint16_t>(
            count(item.at("grade"), source.grade_count() - 1, "section grade"))};
        require(item.at("values").size() == fibers.stalks().total_dimension(grade),
                "section direct-sum dimension");
        Vector v;
        for (std::size_t j = 0; j < item.at("values").size(); ++j)
            v.push_back(scalar(item.at("values").at(j)));
        values.push_back(std::move(v));
    }
    if (validate_only)
        return Json{};
    auto out = Json::object(), results = Json::array(), groups = Json::array();
    std::size_t incompatible = 0;
    if (construction == "cellular_sheaf") {
        for (const auto &group : topology::sheaf_cohomology(sheaf)) {
            auto r = Json::object();
            r.set("grade", cardinality(group.grade.value));
            r.set("dimension", cardinality(group.dimension));
            r.set("cocycles", matrix(group.cocycle_frame));
            r.set("coboundaries", matrix(group.coboundary_frame));
            r.set("representatives", matrix(group.representative_frame));
            groups.push(r);
        }
    }
    std::map<std::uint16_t, std::pair<M, std::unique_ptr<PreparedProjection>>> prepared;
    for (std::size_t i = 0; i < sections.size(); ++i) {
        auto r = Json::parse(sections.at(i).dump());
        r.set("input_origin",
              Json::string(sections.at(i).has("fields") ? "admitted source fields"
                                                        : "declared query operand"));
        const auto grade =
            core::Grade{static_cast<std::uint16_t>(sections.at(i).at("grade").integer())};
        const auto &v = values[i];
        if (construction == "cellular_sheaf") {
            if (!prepared.contains(grade.value)) {
                auto delta = sheaf.coboundary(grade);
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
            r.set("compatible", Json::boolean(compatible));
            r.set("glued", vector(glued));
            r.set("residual", vector(residual));
            r.set("coboundary_defect", vector(defect));
            r.set("coboundary", matrix(delta));
            r.set("residual_quadrance", Json::string(numeric::dot(residual, residual).get_str()));
        } else {
            require(grade.value > 0, "fiber relation sections require a positive grade");
            auto branches = Json::array();
            const auto offsets = fibers.stalks().offsets(grade);
            for (std::size_t cell = 0; cell < source.grade(grade).cell_count(); ++cell) {
                Vector local_value(v.begin() + offsets[cell], v.begin() + offsets[cell + 1]);
                const auto section = fibers.realize_relation(grade, cell, local_value);
                auto branch = Json::object();
                branch.set(
                    "key",
                    Json::string(std::to_string(
                        storage::stable_cell_key(source, {grade, static_cast<std::uint32_t>(cell)})
                            .key)));
                branch.set("upper", vector(local_value));
                auto legs = Json::array();
                for (const auto &leg : section.boundary_values) {
                    auto item = Json::object();
                    item.set("slot", cardinality(leg.incidence_slot));
                    item.set("lower_key", Json::string(std::to_string(
                                              storage::stable_cell_key(source, leg.lower).key)));
                    item.set("value", vector(leg.value));
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
        auto maps_out = Json::array();
        for (std::size_t k = 0; k < source.grade_count(); ++k)
            maps_out.push(matrix(sheaf.coboundary(core::Grade{static_cast<std::uint16_t>(k)})));
        out.set("cochain_maps", maps_out);
    }
    out.set("stalks", records);
    out.set("maps", maps);
    out.set("cohomology", groups);
    out.set("sections", results);
    out.set("incompatible_sections", cardinality(incompatible));
    out.set("reconstruction_exact", Json::boolean(true));
    return out;
}
} // namespace gadget
