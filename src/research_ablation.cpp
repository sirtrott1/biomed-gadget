#include "research_ablation.hpp"
#include "scientific_ablation_fields.hpp"
#include "scientific_lattice_snapshot.hpp"
#include "scientific_sections.hpp"
#include "scientific_relational_field.hpp"
#include <gadget/resources.hpp>
#include <gadget/sections.hpp>
#include <rfc/algebra/materialize.hpp>
#include <rfc/algebra/relative.hpp>
#include <rfc/local/section.hpp>
#include <rfc/runtime/cpu_tiles.hpp>
#include <rfc/storage/local_blocks.hpp>
#include <rfc/temporal/transfer.hpp>
#include <rfc/topology/integral_homology.hpp>

namespace gadget::research_detail {
namespace {
using namespace rfc;
using namespace scientific_detail;
using M = numeric::ExactMatrix;
constexpr std::array<const char *, 4> sectors{"full", "gradient", "curl", "harmonic"};
M boundary(const core::RelationalComplex &c, std::size_t k) {
    if (!k)
        return M(0, c.grade(core::Grade{0}).cell_count());
    if (k == c.grade_count())
        return M(c.grade(core::Grade{static_cast<std::uint16_t>(k - 1)}).cell_count(), 0);
    return field::materialize_boundary_matrix(c, core::Grade{static_cast<std::uint16_t>(k)});
}
struct SectionResult {
    LocalSectionAnalysis analysis;
    NativeRecord publication;
};
std::shared_ptr<SectionResult> section_result(const NativeStudy &study, const Json &query,
    const std::map<std::string, Vector> &operands, const NativeRecord &source_publication) {
    storage::TemporalEndpoint endpoint(study.source());
    using N = NativeRecordValue;
    const auto source_owner = NativeRecord::create_from_endpoint(
        N(N::Object{{{"schema", N(std::string("gadget.scientific.context.v1"))}}}),
        new_native_record_namespace(), source_publication.raw_sources(), std::move(endpoint));
    ResultFields fields(source_owner);
    auto analysis = analyze_local_section_fields(study, query, &fields, &operands, false, true);
    auto publication = fields.finish(std::move(analysis.view));
    return std::make_shared<SectionResult>(SectionResult{std::move(analysis), std::move(publication)});
}
struct Variant {
    algebra::SubcomplexMask mask;
    std::vector<M> inclusions;
    std::shared_ptr<const NativeStudy> study;
    std::shared_ptr<const PreparedRelationalField> prepared;
    std::shared_ptr<const PreparedProjection> removed_span;
    std::set<std::string> span_ids;
    ResultValue topology, sheaf;
    Json declaration;
    ResultValue coordinates, removed, role_witnesses;
    std::shared_ptr<SectionResult> sections;
    explicit Variant(const core::RelationalComplex &source) : mask(source) {}
};
Json restrict_sheaf(const Json &declared, const core::RelationalComplex &ambient,
                    const core::RelationalComplex &retained, const std::vector<M> &inclusions,
                    const std::map<std::string, Vector> &original_operands,
                    std::map<std::string, Vector> &restricted_operands) {
    auto out = Json::parse(declared.dump()), stalks = Json::array(), maps = Json::array();
    const storage::StableSourceIndex original(ambient), target(retained);
    local::LinearStalkBundle bundle(ambient);
    for (std::size_t i = 0; i < declared.at("stalks").size(); ++i) {
        const auto item = declared.at("stalks").at(i);
        const storage::StableCellKey key{{static_cast<std::uint16_t>(item.at("grade").integer())},
                                         scientific_detail::key(item.at("key"))};
        bundle.set_dimension(original.cell(key).value(), item.at("axes").size());
        if (target.cell(key))
            stalks.push(item);
    }
    for (std::size_t i = 0; i < declared.at("maps").size(); ++i) {
        const auto item = declared.at("maps").at(i);
        const storage::StableOccurrenceKey key{
            {static_cast<std::uint16_t>(item.at("grade").integer())},
            {scientific_detail::key(item.at("key"))}};
        if (target.occurrence(key))
            maps.push(item);
    }
    auto sections = Json::array();
    for (std::size_t s = 0; s < declared.at("sections").size(); ++s) {
        auto item = Json::parse(declared.at("sections").at(s).dump());
        const auto k = static_cast<std::uint16_t>(item.at("grade").integer());
        const auto offsets = bundle.offsets(core::Grade{k});
        Vector values;
        // Columns of the certified inclusion retain the ambient cell order.
        const auto &original = original_operands.at(item.at("id").text());
        for (const auto &[ij, coefficient] : inclusions[k].entries()) {
            require(coefficient == 1, "sheaf restriction requires the source inclusion");
            for (auto a = offsets[ij.first]; a < offsets[ij.first + 1]; ++a)
                values.push_back(original.at(a));
        }
        require(restricted_operands.emplace(item.at("id").text(), std::move(values)).second,
                "duplicate restricted native section operand");
        json_object_object_del(item.raw(), "fields");
        json_object_object_del(item.raw(), "values");
        sections.push(item);
    }
    out.set("stalks", stalks);
    out.set("maps", maps);
    out.set("sections", sections);
    return out;
}
std::vector<M> stalk_inclusions(const Json &declared, const core::RelationalComplex &ambient,
                                const std::vector<M> &cell_inclusions) {
    const storage::StableSourceIndex index(ambient);
    local::LinearStalkBundle bundle(ambient);
    for (std::size_t i = 0; i < declared.at("stalks").size(); ++i) {
        const auto s = declared.at("stalks").at(i);
        bundle.set_dimension(
            index.cell({{static_cast<std::uint16_t>(s.at("grade").integer())}, key(s.at("key"))})
                .value(),
            s.at("axes").size());
    }
    std::vector<M> result;
    for (std::size_t k = 0; k < ambient.grade_count(); ++k) {
        const auto offsets = bundle.offsets(core::Grade{static_cast<std::uint16_t>(k)});
        std::size_t dimension = 0;
        for (const auto &[ij, coefficient] : cell_inclusions[k].entries()) {
            (void)coefficient;
            dimension =
                machine::checked_size_add(dimension, offsets[ij.first + 1] - offsets[ij.first]);
        }
        M inclusion(offsets.back(), dimension);
        std::size_t col = 0;
        for (const auto &[ij, coefficient] : cell_inclusions[k].entries()) {
            (void)coefficient;
            for (auto a = offsets[ij.first]; a < offsets[ij.first + 1]; ++a)
                inclusion.set(a, col++, Rational{1});
        }
        result.push_back(std::move(inclusion));
    }
    return result;
}
ResultValue check_sheaf(const SectionResult &base, const SectionResult &next,
    const Json &declaration, const std::vector<M> &inclusions, AblationFields &publication) {
    const auto &old_axes = publication.stalks(base.publication, base.analysis.grade_carriers);
    const auto &new_axes = publication.stalks(next.publication, next.analysis.grade_carriers);
    const auto archive = publication.fields.publication(next.publication);
    auto certificates = ResultValue::array();
    for (std::size_t k = 0; k + 1 < inclusions.size(); ++k) {
        const auto &delta = base.analysis.coboundaries.at(k), &delta_next = next.analysis.coboundaries.at(k);
        const auto defect = inclusions[k + 1].transpose().multiply(delta).add(
            delta_next.multiply(inclusions[k].transpose()).scale(Rational{-1}));
        require(defect.nonzeros() == 0, "full sheaf restriction cochain family");
        auto item = ResultValue::object();
        item.set("grade", ResultValue::count(k));
        item.set("restriction", AblationFields::matrix(inclusions[k].transpose(), old_axes[k], new_axes[k]));
        item.set("defect", AblationFields::matrix(defect, old_axes[k], new_axes[k + 1]));
        certificates.push(std::move(item));
    }
    auto sections = ResultValue::array();
    for (std::size_t s = 0; s < declaration.at("sections").size(); ++s) {
        const auto k = static_cast<std::size_t>(declaration.at("sections").at(s).at("grade").integer());
        const auto observed = inclusions.at(k).transpose().apply(base.analysis.glued.at(s));
        const auto change = numeric::sub(next.analysis.glued.at(s), observed);
        auto additions = ResultValue::object();
        additions.set("restricted_original_compatible", AblationFields::vector(observed, new_axes[k]));
        additions.set("projection_order_defect", AblationFields::vector(change, new_axes[k]));
        additions.set("projection_order_quadrance", publication.scalar(numeric::dot(change, change)));
        auto selected = ResultFields::member(archive, {"view", "sections"});
        selected = ResultFields::member(std::move(selected), s);
        sections.push(ResultFields::overlay(std::move(selected), std::move(additions)));
    }
    auto additions = ResultValue::object();
    additions.set("cochain_restriction_certificates", std::move(certificates));
    additions.set("sections", std::move(sections));
    return ResultFields::overlay(ResultFields::member(archive, {"view"}), std::move(additions));
}

} // namespace

ResultValue source_structure_ablation(const NativeStudy &base, const Json &query, core::Grade grade,
                               const std::map<std::string, Vector> &fields,
                               const std::map<std::string, SourceAblationRole> &roles,
                               const ResultValue &witness, const Json &sheaf, bool validate,
                               const std::shared_ptr<CpuResources> &resources,
                               const NativeRecord &source_publication,
                               const std::map<std::string, Vector> &sheaf_operands,
                               AblationFields *publication) {
    members(query,
            {"operation", "context", "fields", "roles", "ablations", "homology_output", "lambda",
             "workers", "interpretation"},
            {"sheaf", "lattice_snapshots", "comparisons"});
    require(query.at("homology_output").text() == "ranks" ||
                query.at("homology_output").text() == "representatives",
            "homology output profile");
    const bool reps = query.at("homology_output").text() == "representatives";
    require(!query.has("lattice_snapshots") ||
                json_object_is_type(query.at("lattice_snapshots").raw(), json_type_boolean),
            "lattice snapshots must be boolean");
    const bool integral = query.has("lattice_snapshots") && query.at("lattice_snapshots").boolean();
    const auto lambda = exact_rational(text(query.at("lambda"), source_byte_limit()));
    require(lambda > 0, "ablation Green parameter must be positive");
    const auto interpretation = query.at("interpretation");
    members(interpretation, {"population", "meaning", "evidence"});
    (void)text(interpretation.at("population"), source_byte_limit());
    (void)text(interpretation.at("meaning"), source_byte_limit());
    const auto &source = base.source();
    const auto native_cells = std::move(witness.at("cells")).release();
    const auto cell_record = [&](const ResultValue &, core::CellId cell) {
        const auto &grades = std::get<NativeRecordValue::Array>(native_cells.data).elements;
        const auto &cells = std::get<NativeRecordValue::Array>(grades.at(cell.grade.value).data).elements;
        return ResultValue(cells.at(cell.index));
    };
    (void)PreparedRelationalField::chain(base);
    auto cpu = resources ? resources : CpuResources::local();
    const auto workers =
        count(query.at("workers"), cpu->coordinator().capacity(), "ablation workers");
    require(workers > 0, "ablation workers must be positive");
    source_array(query.at("ablations"), "structural ablations");
    require(query.at("ablations").size(), "ablation family must be nonempty");
    std::vector<Variant> variants;
    std::set<std::string> ids, used_roles;
    std::map<std::vector<std::vector<std::uint32_t>>, std::shared_ptr<const NativeStudy>> owners;
    std::vector<std::vector<std::uint32_t>> all_coordinates(source.grade_count());
    for (std::size_t k = 0; k < source.grade_count(); ++k)
        for (std::size_t i = 0;
             i < source.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count(); ++i)
            all_coordinates[k].push_back(static_cast<std::uint32_t>(i));
    owners.emplace(all_coordinates,
                   std::shared_ptr<const NativeStudy>(&base, [](const NativeStudy *) {}));
    for (std::size_t a = 0; a < query.at("ablations").size(); ++a) {
        const auto declared = query.at("ablations").at(a);
        members(declared, {"id", "remove_roles"});
        require(ids.insert(text(declared.at("id"), source_byte_limit())).second,
                "duplicate ablation identity");
        variants.emplace_back(source);
        auto &v = variants.back();
        v.declaration = declared;
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const auto g = core::Grade{static_cast<std::uint16_t>(k)};
            for (std::size_t i = 0; i < source.grade(g).cell_count(); ++i)
                v.mask.select({g, static_cast<std::uint32_t>(i)});
        }
        std::set<std::string> selected_roles, span_ids;
        source_array(declared.at("remove_roles"), "ablation role selections");
        require(declared.at("remove_roles").size(), "ablation requires a selected role");
        v.role_witnesses = ResultValue::array();
        for (std::size_t r = 0; r < declared.at("remove_roles").size(); ++r) {
            const auto id = declared.at("remove_roles").at(r).text();
            require(roles.contains(id) && selected_roles.insert(id).second,
                    "unknown or repeated ablation role");
            used_roles.insert(id);
            const auto &role = roles.at(id);
            v.role_witnesses.push(role.witness);
            for (const auto cell : role.cells)
                v.mask.select(cell, false);
            span_ids.insert(role.fields.begin(), role.fields.end());
        }
        require(v.mask.boundary_closed(source),
                "source deletion is not boundary closed; explicitly declare the additional cofaces "
                "to remove, including inactive and cancelled attachments");
        std::vector<std::vector<std::uint32_t>> selected(source.grade_count());
        for (std::size_t k = 0; k < source.grade_count(); ++k)
            v.mask.grade(core::Grade{static_cast<std::uint16_t>(k)})
                .for_each_set(
                    [&](std::size_t i) { selected[k].push_back(static_cast<std::uint32_t>(i)); });
        std::vector<field::ExactMetric> metrics;
        v.coordinates = ResultValue::array();
        v.removed = ResultValue::array();
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            M map(source.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count(),
                  selected[k].size());
            auto coordinates = ResultValue::array(), removed = ResultValue::array();
            Vector diagonal;
            const auto metric = base.tower().metrics()[k].materialize();
            require(base.tower().metrics()[k].kind() == field::ExactMetricKind::Identity ||
                        base.tower().metrics()[k].kind() == field::ExactMetricKind::Diagonal,
                    "source restriction metric profile requires identity or diagonal metrics");
            for (std::size_t i = 0; i < selected[k].size(); ++i) {
                const auto c = selected[k][i];
                map.set(c, i, Rational{1});
                diagonal.push_back(metric.get(c, c));
                coordinates.push(
                    cell_record(witness, {core::Grade{static_cast<std::uint16_t>(k)}, c}));
            }
            for (std::size_t i = 0; i < map.rows(); ++i)
                if (!v.mask.contains({core::Grade{static_cast<std::uint16_t>(k)},
                                      static_cast<std::uint32_t>(i)}))
                    removed.push(cell_record(witness, {core::Grade{static_cast<std::uint16_t>(k)},
                                                       static_cast<std::uint32_t>(i)}));
            metrics.push_back(field::ExactMetric::diagonal(std::move(diagonal)));
            v.inclusions.push_back(std::move(map));
            v.coordinates.push(coordinates);
            v.removed.push(removed);
        }
        if (!owners.contains(selected)) {
            auto materialized = algebra::materialize_subcomplex(source, v.mask);
            require(materialized.to_ambient == selected,
                    "native materialization coordinate witness");
            owners.emplace(selected, std::make_shared<const NativeStudy>(
                                         std::move(materialized.complex), std::move(metrics)));
        }
        v.study = owners.at(selected);
        const algebra::RelationalMorphism inclusion(v.study->source(), source, selected);
        require(inclusion.is_strict() && inclusion.is_chain_map(),
                "strict retained source inclusion");
        v.span_ids = std::move(span_ids);
    }
    require(used_roles.size() == roles.size(), "unused scientific role declaration");
    std::set<std::string> comparison_ids;
    if (query.has("comparisons")) {
        source_array(query.at("comparisons"), "ablation comparisons");
        for (std::size_t c = 0; c < query.at("comparisons").size(); ++c) {
            const auto item = query.at("comparisons").at(c);
            members(item, {"id", "left", "right"});
            require(comparison_ids.insert(text(item.at("id"), source_byte_limit())).second &&
                        ids.contains(item.at("left").text()) &&
                        ids.contains(item.at("right").text()),
                    "unknown or repeated ablation comparison");
        }
    }
    auto original_operands = sheaf_operands;
    if (sheaf.raw() != nullptr) {
        for (std::size_t i = 0; i < sheaf.at("sections").size(); ++i) {
            const auto section = sheaf.at("sections").at(i);
            if (original_operands.contains(section.at("id").text())) continue;
            require(section.has("values") && !section.has("fields"),
                    "original section native source operand absent");
            Vector values;
            for (std::size_t j = 0; j < section.at("values").size(); ++j)
                values.push_back(exact_rational(section.at("values").at(j).text()));
            original_operands.emplace(section.at("id").text(), std::move(values));
        }
        std::set<std::size_t> section_grades;
        for (std::size_t s = 0; s < sheaf.at("sections").size(); ++s)
            require(section_grades
                        .insert(count(sheaf.at("sections").at(s).at("grade"),
                                      source.grade_count() - 1, "sheaf section grade"))
                        .second,
                    "ablation sheaf profile requires one section per grade");
        require(section_grades.size() == source.grade_count(),
                "ablation sheaf gluing requires sections on every grade");
    }
    if (sheaf.raw() != nullptr)
        (void)analyze_local_section_fields(base, sheaf, nullptr, &original_operands, true);
    if (validate)
        return ResultValue::object();
    require(publication != nullptr, "ablation requires native result field ownership");
    const auto &original_axes = publication->source(base);
    for (const auto &v : variants) (void)publication->source(*v.study);
    const auto q = [&](const Rational &value) {return publication->scalar(value);};
    const auto vector = [](const Vector &value, std::uint64_t axis) {return AblationFields::vector(value, axis);};
    const auto matrix = [](const M &value, std::uint64_t input, std::uint64_t output) {
        return AblationFields::matrix(value, input, output);
    };
    auto reservation = cpu->reserve(workers);
    std::map<std::set<std::string>, std::shared_ptr<const PreparedProjection>> span_owners;
    for (auto &v : variants)
        if (!v.span_ids.empty()) {
            if (!span_owners.contains(v.span_ids)) {
                M frame(base.tower().metrics()[grade.value].size(), v.span_ids.size());
                std::size_t c = 0;
                for (const auto &id : v.span_ids) {
                    for (std::size_t i = 0; i < frame.rows(); ++i)
                        frame.set(i, c, fields.at(id).at(i));
                    ++c;
                }
                span_owners.emplace(v.span_ids,
                                    std::make_shared<const PreparedProjection>(
                                        frame.image_frame(), base.tower().metrics()[grade.value]));
            }
            v.removed_span = span_owners.at(v.span_ids);
        }
    auto span_bindings = ResultValue::array();
    for (const auto &[names, preparation] : span_owners) {
        auto binding = ResultValue::object(), members = ResultValue::array();
        for (const auto &name : names) members.push(ResultValue::string(name));
        const auto parameters = publication->fields.axis(preparation->frame.cols());
        binding.set("fields", std::move(members));
        binding.set("frame", matrix(preparation->frame, parameters, original_axes.rows[grade.value]));
        binding.set("rank", ResultValue::count(preparation->frame.cols()));
        span_bindings.push(std::move(binding));
    }
    auto spans = ResultValue::object(); spans.set("kind", ResultValue::string("admitted_removed_field_spans"));
    spans.set("preparations", std::move(span_bindings)); publication->fields.bind(std::move(spans));
    const auto original_owner = std::make_shared<const PreparedRelationalField>(base, grade);
    const auto &original = *original_owner;
    std::map<const NativeStudy *, std::shared_ptr<const PreparedRelationalField>> prepared_owners{
        {&base, original_owner}};
    struct SourceResult { ResultValue topology, sheaf; std::shared_ptr<SectionResult> sections; };
    std::map<const NativeStudy *, SourceResult> source_results;
    std::vector<M> original_boundaries, original_images;
    for (std::size_t k = 0; k <= source.grade_count(); ++k) {
        original_boundaries.push_back(boundary(source, k));
        original_images.push_back(original_boundaries.back().image_frame());
    }
    auto base_grades = ResultValue::array();
    std::vector<std::size_t> original_betti;
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        auto g = ResultValue::object();
        const auto n = original_boundaries[k].cols();
        g.set("grade", ResultValue::count(k));
        g.set("cells", ResultValue::count(n));
        original_betti.push_back(n - original_images[k].cols() - original_images[k + 1].cols());
        g.set("betti", ResultValue::count(original_betti.back()));
        const auto empty = publication->fields.axis(0);
        g.set("boundary", matrix(original_boundaries[k], original_axes.rows[k], k ? original_axes.rows[k - 1] : empty));
        const auto image_axis = publication->fields.axis(original_images[k + 1].cols());
        g.set("boundary_image_frame", matrix(original_images[k + 1], image_axis, original_axes.rows[k]));
        if (integral)
            g.set("lattice", lattice_snapshot(source, core::Grade{static_cast<std::uint16_t>(k)}, publication->fields, original_axes.lattice));
        base_grades.push(g);
    }
    ResultValue original_sheaf;
    std::shared_ptr<SectionResult> original_sections;
    if (sheaf.raw() != nullptr) {
        original_sections = section_result(base, sheaf, original_operands, source_publication);
        original_sheaf = ResultFields::member(publication->fields.publication(original_sections->publication), {"view"});
    }
    for (auto &v : variants) {
        if (!prepared_owners.contains(v.study.get()))
            prepared_owners.emplace(
                v.study.get(), std::make_shared<const PreparedRelationalField>(*v.study, grade));
        v.prepared = prepared_owners.at(v.study.get());
        if (source_results.contains(v.study.get())) {
            v.topology = source_results.at(v.study.get()).topology;
            v.sheaf = source_results.at(v.study.get()).sheaf;
            v.sections = source_results.at(v.study.get()).sections;
            continue;
        }
        auto grades = ResultValue::array();
        const auto &retained_axes = publication->source(*v.study);
        algebra::SubcomplexMask whole(source);
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const auto g = core::Grade{static_cast<std::uint16_t>(k)};
            for (std::size_t i = 0; i < source.grade(g).cell_count(); ++i)
                whole.select({g, static_cast<std::uint32_t>(i)});
        }
        const algebra::RelativeChainComplex relative(source, whole, v.mask);
        std::vector<M> relative_boundaries{M(0, relative.cell_count(core::Grade{0}))};
        for (std::size_t k = 1; k < source.grade_count(); ++k) {
            const auto g = core::Grade{static_cast<std::uint16_t>(k)};
            M b(relative.cell_count(core::Grade{static_cast<std::uint16_t>(k - 1)}),
                relative.cell_count(g));
            for (std::size_t i = 0; i < b.cols(); ++i) {
                Vector e(b.cols());
                e[i] = 1;
                const auto column = relative.boundary_apply(g, e);
                for (std::size_t j = 0; j < b.rows(); ++j)
                    b.set(j, i, column[j]);
            }
            relative_boundaries.push_back(std::move(b));
        }
        relative_boundaries.emplace_back(relative.cell_count(source.max_grade()), 0);
        std::vector<std::uint64_t> relative_axes;
        auto relative_binding = ResultValue::object(), relative_grades = ResultValue::array();
        relative_binding.set("kind", ResultValue::string("actual_relative_quotient"));
        relative_binding.set("ambient_source", original_axes.reference);
        relative_binding.set("subcomplex_source", retained_axes.reference);
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const core::Grade g{static_cast<std::uint16_t>(k)};
            const auto axis = publication->fields.axis(relative.cell_count(g));
            relative_axes.push_back(axis);
            std::vector<std::size_t> coordinates;
            for (std::size_t i = 0; i < source.grade(g).cell_count(); ++i)
                if (!v.mask.contains({g, static_cast<std::uint32_t>(i)})) coordinates.push_back(i);
            require(coordinates.size() == relative.cell_count(g), "complete quotient coordinate owner");
            M projection(coordinates.size(), source.grade(g).cell_count());
            for (std::size_t i = 0; i < coordinates.size(); ++i) projection.set(i, coordinates[i], Rational{1});
            auto record = ResultValue::object(); record.set("grade", ResultValue::count(k));
            record.set("carrier", ResultValue::count(axis));
            record.set("ambient_coordinates", publication->fields.indices(coordinates, axis));
            record.set("projection", matrix(projection, original_axes.rows[k], axis));
            relative_grades.push(std::move(record));
        }
        relative_binding.set("grades", std::move(relative_grades));
        publication->fields.bind(std::move(relative_binding));
        const auto empty_axis = publication->fields.axis(0);
        std::vector<M> retained_boundaries;
        std::vector<std::size_t> retained_ranks, relative_ranks;
        for (std::size_t k = 0; k <= source.grade_count(); ++k) {
            retained_boundaries.push_back(boundary(v.study->source(), k));
            retained_ranks.push_back(retained_boundaries.back().rank());
            relative_ranks.push_back(relative_boundaries[k].rank());
        }
        for (std::size_t k = 0; k < source.grade_count(); ++k) {
            const auto &low = retained_boundaries[k], &upper = retained_boundaries[k + 1];
            const auto &j = v.inclusions[k];
            const auto n = j.cols();
            if (k)
                require(original_boundaries[k]
                                .multiply(j)
                                .add(v.inclusions[k - 1].multiply(low).scale(Rational{-1}))
                                .nonzeros() == 0,
                        "complete inclusion chain family");
            if (k + 1 < source.grade_count())
                require(
                    low.multiply(upper).nonzeros() == 0 &&
                        relative_boundaries[k].multiply(relative_boundaries[k + 1]).nonzeros() == 0,
                    "retained and relative chain families");
            const auto pullback =
                j.transpose().multiply(base.tower().metrics()[k].materialize()).multiply(j);
            require(temporal::metric_pullback_defect(v.study->tower().metrics()[k],
                                                     base.tower().metrics()[k], j)
                            .nonzeros() == 0,
                    "native metric pullback");
            const auto &f = original_images[k + 1];
            M outside(f.rows() - n, f.cols());
            std::size_t row = 0;
            for (std::size_t i = 0; i < f.rows(); ++i)
                if (!v.mask.contains({core::Grade{static_cast<std::uint16_t>(k)},
                                      static_cast<std::uint32_t>(i)})) {
                    for (std::size_t c = 0; c < f.cols(); ++c)
                        outside.set(row, c, f.get(i, c));
                    ++row;
                }
            // Independent F identifies ambient boundaries supported in the
            // subcomplex by ker(outside F). Strict inclusion makes them cycles.
            const auto intersection = f.cols() - outside.rank(), cycles = n - retained_ranks[k],
                       betti = cycles - retained_ranks[k + 1];
            require(intersection >= retained_ranks[k + 1] && intersection <= cycles,
                    "inclusion boundary intersection");
            const auto image = cycles - intersection;
            require(image <= betti &&
                        image <= original_betti.at(k),
                    "homology image rank");
            auto g = ResultValue::object();
            g.set("grade", ResultValue::count(k));
            g.set("cells", ResultValue::count(n));
            g.set("boundary", matrix(low, retained_axes.rows[k], k ? retained_axes.rows[k - 1] : empty_axis));
            g.set("inclusion", matrix(j, retained_axes.rows[k], original_axes.rows[k]));
            g.set("cochain_restriction", matrix(j.transpose(), original_axes.rows[k], retained_axes.rows[k]));
            g.set("pullback_metric", matrix(pullback, retained_axes.rows[k], retained_axes.rows[k]));
            g.set("betti", ResultValue::count(betti));
            g.set("homology_image_rank", ResultValue::count(image));
            g.set("killed_by_inclusion", ResultValue::count(betti - image));
            g.set("ambient_classes_outside_image",
                  ResultValue::count(original_betti.at(k) - image));
            g.set("relative_boundary", matrix(relative_boundaries[k], relative_axes[k], k ? relative_axes[k - 1] : empty_axis));
            g.set("relative_betti", ResultValue::count(relative_boundaries[k].cols() - relative_ranks[k] -
                                                relative_ranks[k + 1]));
            if (integral)
                g.set("lattice",
                      lattice_snapshot(v.study->source(), core::Grade{static_cast<std::uint16_t>(k)}, publication->fields, retained_axes.lattice));
            if (reps && k == grade.value) {
                const auto h = v.study->hodge(grade).harmonic_frame(), included = j.multiply(h);
                M mapped(j.rows(), h.cols());
                for (std::size_t c = 0; c < h.cols(); ++c) {
                    const auto value = original.apply(field::matrix_column(included, c))[3];
                    for (std::size_t i = 0; i < value.size(); ++i)
                        mapped.set(i, c, value[i]);
                }
                const auto kernel = mapped.kernel_frame(), killed = included.multiply(kernel);
                require(h.cols() == betti && mapped.rank() == image &&
                            kernel.cols() == betti - image,
                        "homology representatives match rank-nullity");
                std::vector<std::size_t> pivots;
                (void)original_boundaries[k + 1].rref(pivots);
                PreparedProjection filling(original_boundaries[k + 1].select_columns(pivots),
                                           original.metric);
                M preimages(original_boundaries[k + 1].cols(), killed.cols());
                for (std::size_t c = 0; c < killed.cols(); ++c) {
                    const auto value = field::matrix_column(killed, c),
                               coefficients = filling.coefficients(value);
                    require(filling.frame.apply(coefficients) == value,
                            "killed class has an ambient filling");
                    for (std::size_t i = 0; i < coefficients.size(); ++i)
                        preimages.set(pivots[i], c, coefficients[i]);
                }
                require(original_boundaries[k + 1].multiply(preimages).entries() ==
                            killed.entries(),
                        "full ambient filling witness");
                const auto harmonic_axis = publication->fields.axis(h.cols());
                const auto killed_axis = publication->fields.axis(kernel.cols());
                g.set("harmonic_representatives", matrix(h, harmonic_axis, retained_axes.rows[k]));
                g.set("included_representatives", matrix(included, harmonic_axis, original_axes.rows[k]));
                g.set("ambient_harmonic_images", matrix(mapped, harmonic_axis, original_axes.rows[k]));
                g.set("killed_class_coefficients", matrix(kernel, killed_axis, harmonic_axis));
                g.set("killed_representatives", matrix(killed, killed_axis, original_axes.rows[k]));
                g.set("ambient_fillings", matrix(preimages, killed_axis, k + 1 < original_axes.rows.size() ? original_axes.rows[k + 1] : empty_axis));
            }
            grades.push(g);
        }
        v.topology = grades;
        if (sheaf.raw() != nullptr) {
            if (v.study.get() == &base) v.sections = original_sections;
            else {
                std::map<std::string, Vector> operands;
                const auto restricted = restrict_sheaf(sheaf, source, v.study->source(),
                    v.inclusions, original_operands, operands);
                v.sections = section_result(*v.study, restricted, operands, source_publication);
            }
            v.sheaf = check_sheaf(*original_sections, *v.sections, sheaf,
                stalk_inclusions(sheaf, source, v.inclusions), *publication);
        }
        source_results.emplace(v.study.get(), SourceResult{v.topology, v.sheaf, v.sections});
    }
    const auto ids_query = query.at("fields");
    const auto jobs = machine::checked_size_mul(variants.size(), ids_query.size());
    std::vector<ResultValue> readings(jobs);
    std::vector<std::array<Vector, 4>> retained_sectors(jobs);
    std::vector<std::array<Vector, 4>> originals;
    std::vector<Vector> original_green;
    for (std::size_t f = 0; f < ids_query.size(); ++f) {
        const auto &x = fields.at(ids_query.at(f).text());
        originals.push_back(original.apply(x));
        const auto solved = base.green_action(grade, lambda, x);
        require(solved.exact, "original Green solve");
        original_green.push_back(solved.value);
    }
    runtime::CpuTileExecutor executor({workers, cpu->coordinator().capacity(), jobs});
    const auto report = executor.run_tiles(jobs, [&](std::size_t job, std::size_t,
                                                     const runtime::CpuTileStop &stop) {
        const auto a = job / ids_query.size(), f = job % ids_query.size();
        const auto &v = variants[a];
        const auto &axes = publication->source(*v.study);
        const auto &x = originals[f][0];
        const auto removed =
            v.removed_span ? v.removed_span->apply(x, [&] { stop.check(); }) : Vector(x.size());
        const auto y = numeric::sub(x, removed),
                   restricted = v.inclusions[grade.value].transpose().apply(y);
        const auto fields_next = v.prepared->apply(restricted, [&] { stop.check(); });
        auto record = ResultValue::object(), quadrances = ResultValue::object(), fields_json = ResultValue::object(),
             transitions = ResultValue::array(), defects = ResultValue::object();
        const auto restriction = v.inclusions[grade.value].transpose();
        const auto innovation = restriction.apply(numeric::sub(y, x));
        record.set("field", ids_query.at(f));
        record.set("removed_span_field", vector(removed, original_axes.rows[grade.value]));
        record.set("retained_span_field", vector(y, original_axes.rows[grade.value]));
        record.set("restricted_input", vector(restricted, axes.rows[grade.value]));
        record.set("discarded_source_field",
                   vector(numeric::sub(y, v.inclusions[grade.value].apply(restricted)), original_axes.rows[grade.value]));
        record.set("field_innovation", vector(innovation, axes.rows[grade.value]));
        std::array<Vector, 3> reconstructed{Vector(restricted.size()), Vector(restricted.size()),
                                            Vector(restricted.size())};
        for (std::size_t b = 1; b < sectors.size(); ++b) {
            const auto action =
                v.prepared->apply(restriction.apply(originals[f][b]), [&] { stop.check(); });
            for (std::size_t a_sector = 1; a_sector < sectors.size(); ++a_sector) {
                auto t = ResultValue::object();
                t.set("from", ResultValue::string(sectors[b]));
                t.set("to", ResultValue::string(sectors[a_sector]));
                t.set("value", vector(action[a_sector], axes.rows[grade.value]));
                transitions.push(t);
                reconstructed[a_sector - 1] =
                    numeric::add(reconstructed[a_sector - 1], action[a_sector]);
            }
        }
        const auto innovative = v.prepared->apply(innovation, [&] { stop.check(); });
        for (std::size_t a_sector = 0; a_sector < sectors.size(); ++a_sector) {
            fields_json.set(sectors[a_sector], vector(fields_next[a_sector], axes.rows[grade.value]));
            quadrances.set(sectors[a_sector], q(v.prepared->metric.inner(fields_next[a_sector],
                                                                         fields_next[a_sector])));
            if (a_sector) {
                require(numeric::add(reconstructed[a_sector - 1], innovative[a_sector]) ==
                            fields_next[a_sector],
                        "full sector transition and innovation reconstruction");
                defects.set(sectors[a_sector],
                            vector(numeric::sub(fields_next[a_sector],
                                                restriction.apply(originals[f][a_sector])), axes.rows[grade.value]));
            }
        }
        const auto next = v.study->green_action(grade, lambda, restricted);
        require(next.exact, "retained Green solve");
        const auto transported = restriction.apply(original_green[f]);
        const auto operator_defect =
            numeric::sub(restriction.apply(base.hodge(grade).apply(original_green[f])),
                         v.study->hodge(grade).apply(transported));
        const auto structure =
            v.study->green_action(grade, lambda, numeric::scale(lambda, operator_defect));
        const auto changed = v.study->green_action(grade, lambda, innovation);
        require(structure.exact && changed.exact &&
                    numeric::add(structure.value, changed.value) ==
                        numeric::sub(next.value, transported),
                "action-based Green defect and innovation reconstruction");
        auto green = ResultValue::object();
        green.set("value", vector(next.value, axes.rows[grade.value]));
        green.set("restricted_original", vector(transported, axes.rows[grade.value]));
        green.set("structure_term", vector(structure.value, axes.rows[grade.value]));
        green.set("field_term", vector(changed.value, axes.rows[grade.value]));
        green.set("change", vector(numeric::sub(next.value, transported), axes.rows[grade.value]));
        record.set("sectors", fields_json);
        record.set("quadrances", quadrances);
        record.set("sector_transport", transitions);
        record.set("sector_change", defects);
        record.set("green", green);
        retained_sectors[job] = fields_next;
        readings[job] = record;
    });
    auto ablations = ResultValue::array(), base_fields = ResultValue::array();
    for (std::size_t f = 0; f < ids_query.size(); ++f) {
        auto r = ResultValue::object(), quadrances = ResultValue::object();
        r.set("field", ids_query.at(f));
        for (std::size_t s = 0; s < sectors.size(); ++s) {
            r.set(sectors[s], vector(originals[f][s], original_axes.rows[grade.value]));
            quadrances.set(sectors[s], q(original.metric.inner(originals[f][s], originals[f][s])));
        }
        r.set("quadrances", quadrances);
        r.set("green", vector(original_green[f], original_axes.rows[grade.value]));
        base_fields.push(r);
    }
    for (std::size_t a = 0; a < variants.size(); ++a) {
        const auto &v = variants[a];
        auto r = ResultValue(v.declaration), values = ResultValue::array();
        r.set("roles", v.role_witnesses);
        r.set("coordinates", v.coordinates);
        r.set("removed_cells", v.removed);
        r.set("grades", v.topology);
        r.set("removed_span_rank", ResultValue::count(v.removed_span ? v.removed_span->frame.cols() : 0));
        for (std::size_t f = 0; f < ids_query.size(); ++f)
            values.push(readings[a * ids_query.size() + f]);
        r.set("fields", values);
        if (sheaf.raw() != nullptr)
            r.set("sheaf", v.sheaf);
        ablations.push(r);
    }
    auto comparisons = ResultValue::array();
    std::vector<std::shared_ptr<const NativeStudy>> common_owners;
    if (query.has("comparisons")) {
        std::map<std::string, std::size_t> by_id;
        for (std::size_t a = 0; a < variants.size(); ++a)
            by_id.emplace(variants[a].declaration.at("id").text(), a);
        for (std::size_t c = 0; c < query.at("comparisons").size(); ++c) {
            auto item = ResultValue(query.at("comparisons").at(c));
            const auto left = by_id.at(item.at("left").text()),
                       right = by_id.at(item.at("right").text());
            const auto &a = variants[left], &b = variants[right];
            const auto overlap = meet(a.mask, b.mask), combined = join(a.mask, b.mask);
            require(overlap.boundary_closed(source) && combined.boundary_closed(source),
                    "native subcomplex meet and join");
            const auto common = algebra::materialize_subcomplex(source, overlap);
            std::array<std::vector<std::vector<std::uint32_t>>, 2> common_maps;
            M inclusion(source.grade(grade).cell_count(), common.to_ambient[grade.value].size());
            auto common_cells = ResultValue::array();
            bool covers = true;
            for (std::size_t k = 0; k < source.grade_count(); ++k) {
                auto records = ResultValue::array();
                std::array<std::map<std::size_t, std::uint32_t>, 2> inverse;
                for (const auto &[ij, value] : a.inclusions[k].entries()) {
                    (void)value;
                    inverse[0].emplace(ij.first, static_cast<std::uint32_t>(ij.second));
                }
                for (const auto &[ij, value] : b.inclusions[k].entries()) {
                    (void)value;
                    inverse[1].emplace(ij.first, static_cast<std::uint32_t>(ij.second));
                }
                for (std::size_t i = 0; i < common.to_ambient[k].size(); ++i) {
                    const auto ambient = common.to_ambient[k][i];
                    records.push(cell_record(
                        witness, {core::Grade{static_cast<std::uint16_t>(k)}, ambient}));
                    if (k == grade.value)
                        inclusion.set(ambient, i, Rational{1});
                }
                for (std::size_t side = 0; side < 2; ++side) {
                    std::vector<std::uint32_t> map;
                    for (const auto ambient : common.to_ambient[k])
                        map.push_back(inverse[side].at(ambient));
                    common_maps[side].push_back(std::move(map));
                }
                common_cells.push(records);
                covers = covers &&
                         combined.grade(core::Grade{static_cast<std::uint16_t>(k)}).count() ==
                             source.grade(core::Grade{static_cast<std::uint16_t>(k)}).cell_count();
            }
            require(
                algebra::RelationalMorphism(common.complex, a.study->source(), common_maps[0])
                        .is_strict() &&
                    algebra::RelationalMorphism(common.complex, b.study->source(), common_maps[1])
                        .is_strict(),
                "common carrier strict inclusions");
            std::vector<field::ExactMetric> common_metrics;
            for (std::size_t k = 0; k < source.grade_count(); ++k) {
                const auto original_metric = base.tower().metrics()[k].materialize();
                Vector diagonal;
                for (const auto coordinate : common.to_ambient[k]) diagonal.push_back(original_metric.get(coordinate, coordinate));
                common_metrics.push_back(field::ExactMetric::diagonal(std::move(diagonal)));
            }
            common_owners.push_back(std::make_shared<const NativeStudy>(common.complex, std::move(common_metrics)));
            const auto &common_axes = publication->source(*common_owners.back());
            item.set("common_cells", common_cells);
            item.set("union_covers_ambient", ResultValue::boolean(covers));
            const auto restrict_left = inclusion.transpose().multiply(a.inclusions[grade.value]),
                       restrict_right = inclusion.transpose().multiply(b.inclusions[grade.value]);
            auto comparisons_fields = ResultValue::array();
            for (std::size_t f = 0; f < ids_query.size(); ++f) {
                auto record = ResultValue::object(), differences = ResultValue::object(),
                     quadrances = ResultValue::object();
                record.set("field", ids_query.at(f));
                for (std::size_t sector = 0; sector < sectors.size(); ++sector) {
                    const auto &av = retained_sectors[left * ids_query.size() + f][sector];
                    const auto &bv = retained_sectors[right * ids_query.size() + f][sector];
                    const auto delta =
                        numeric::sub(restrict_right.apply(bv), restrict_left.apply(av));
                    differences.set(sectors[sector], vector(delta, common_axes.rows[grade.value]));
                    const auto included = inclusion.apply(delta);
                    quadrances.set(sectors[sector], q(original.metric.inner(included, included)));
                }
                record.set("difference_on_common_carrier", differences);
                record.set("quadrances", quadrances);
                comparisons_fields.push(record);
            }
            item.set("fields", comparisons_fields);
            if (sheaf.raw() != nullptr) {
                const storage::StableSourceIndex index(source);
                local::LinearStalkBundle bundle(source);
                for (std::size_t s = 0; s < sheaf.at("stalks").size(); ++s) {
                    const auto stalk = sheaf.at("stalks").at(s);
                    bundle.set_dimension(
                        index
                            .cell({{static_cast<std::uint16_t>(stalk.at("grade").integer())},
                                   key(stalk.at("key"))})
                            .value(),
                        stalk.at("axes").size());
                }
                const auto make_section = [&](const Variant &v) {
                    local::LinearSection section;
                    const auto maps = stalk_inclusions(sheaf, source, v.inclusions);
                    for (std::size_t s = 0; s < sheaf.at("sections").size(); ++s) {
                        const auto record = sheaf.at("sections").at(s);
                        const auto k = static_cast<std::uint16_t>(record.at("grade").integer());
                        const auto included = maps[k].apply(v.sections->analysis.glued.at(s));
                        const auto offsets = bundle.offsets(core::Grade{k});
                        v.mask.grade(core::Grade{k}).for_each_set([&](std::size_t i) {
                            section.emplace(
                                core::CellId{core::Grade{k}, static_cast<std::uint32_t>(i)},
                                Vector(included.begin() + offsets[i],
                                       included.begin() + offsets[i + 1]));
                        });
                    }
                    return section;
                };
                const auto ls = make_section(a), rs = make_section(b);
                const local::RelationalLocalCarrier carrier(source);
                const bool compatible = local::sections_compatible_on_overlap(
                    carrier, ls, a.mask, rs, b.mask, local::OverlapCoverage::RequireComplete);
                item.set("sheaf_overlap_compatible", ResultValue::boolean(compatible));
                auto disagreements = ResultValue::array();
                for (std::size_t k = 0; k < source.grade_count(); ++k) {
                    const auto g = core::Grade{static_cast<std::uint16_t>(k)};
                    overlap.grade(g).for_each_set([&](std::size_t i) {
                        const core::CellId cell{g, static_cast<std::uint32_t>(i)};
                        if (ls.at(cell) != rs.at(cell)) {
                            auto record = ResultValue::object();
                            record.set("cell", cell_record(witness, cell));
                            const auto key = storage::stable_cell_key(source, cell);
                            const auto axis = publication->fields.stalk(key, bundle.dimension(cell));
                            record.set("left", vector(ls.at(cell), axis));
                            record.set("right", vector(rs.at(cell), axis));
                            record.set("difference", vector(numeric::sub(rs.at(cell), ls.at(cell)), axis));
                            disagreements.push(record);
                        }
                    });
                }
                require(compatible == (disagreements.size() == 0),
                        "native overlap obstruction witness");
                item.set("sheaf_overlap_disagreements", disagreements);
                if (compatible) {
                    const auto glued = local::glue_sections(
                        carrier, ls, a.mask, rs, b.mask, local::OverlapCoverage::RequireComplete);
                    auto glued_cells = ResultValue::array();
                    for (std::size_t k = 0; k < source.grade_count(); ++k) {
                        const auto g = core::Grade{static_cast<std::uint16_t>(k)};
                        const auto offsets = bundle.offsets(g);
                        Vector global(offsets.back());
                        combined.grade(g).for_each_set([&](std::size_t i) {
                            const core::CellId cell{g, static_cast<std::uint32_t>(i)};
                            auto record = ResultValue::object();
                            record.set("cell", cell_record(witness, cell));
                            const auto axis = publication->fields.stalk(storage::stable_cell_key(source, cell), bundle.dimension(cell));
                            record.set("value", vector(glued.at(cell), axis));
                            glued_cells.push(record);
                            for (std::size_t d = 0; d < glued.at(cell).size(); ++d)
                                global[offsets[i] + d] = glued.at(cell)[d];
                        });
                        if (covers)
                            require(numeric::is_zero(
                                        original_sections->analysis.coboundaries.at(k)
                                            .apply(global)),
                                    "glued global section cochain compatibility");
                    }
                    item.set("glued_sections", glued_cells);
                    item.set("global_cochain_verified", ResultValue::boolean(covers));
                }
            }
            comparisons.push(item);
        }
    }
    auto result = ResultValue::object(), execution = ResultValue::object();
    execution.set("backend", ResultValue::string("native_cpu_tiles"));
    execution.set("requested_workers", ResultValue::count(workers));
    execution.set("workers", ResultValue::count(report.worker_slots));
    require(report.completed_tiles == jobs, "complete joined field family");
    execution.set("completed_tiles", ResultValue::count(report.completed_tiles));
    execution.set("launched_threads", ResultValue::count(report.launched_threads));
    execution.set("field_actions", ResultValue::count(jobs));
    execution.set("harmonic_basis_requested", ResultValue::boolean(reps));
    execution.set("source_preparations", ResultValue::count(prepared_owners.size()));
    execution.set("span_preparations", ResultValue::count(span_owners.size()));
    execution.set("joined", ResultValue::boolean(true));
    result.set("execution", execution);
    result.set("grade", ResultValue::count(grade.value));
    result.set("source", witness);
    result.set("original_grades", base_grades);
    result.set("original_fields", base_fields);
    result.set("ablations", ablations);
    if (sheaf.raw() != nullptr)
        result.set("original_sheaf", original_sheaf);
    result.set("comparisons", comparisons);
    result.set("coefficient_domain", ResultValue::string("Q"));
    result.set("identities_verified", ResultValue::boolean(true));
    return result;
}
} // namespace gadget::research_detail
