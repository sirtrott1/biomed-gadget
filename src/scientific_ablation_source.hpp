#pragma once
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <gadget/study.hpp>
#include <set>

namespace gadget::research_detail {
// Actual source and metric snapshots for original, retained and overlap
// owners. Quotient carriers remain separate from primary relations.
struct AblationSourcePublication {
    NativeRecord publication;
    std::vector<std::uint64_t> grade_axes;
};
inline AblationSourcePublication ablation_source_publication(
    const NativeStudy &study, const NativeRecord &original,
    const std::vector<std::vector<std::uint32_t>> &to_original) {
    using namespace scientific_detail;
    using N = NativeRecordValue;
    const auto &source = study.source();
    require(to_original.size() == source.grade_count(), "source coordinate lift grade count");
    auto sources = original.raw_sources();
    std::uint64_t archive_key = 1;
    for (const auto &bytes : sources)
        archive_key = std::max(archive_key, rfc::machine::checked_size_add(bytes.key(), 1));
    auto archive = original.archive_source(archive_key);
    const auto archive_size = archive.bytes().size();
    sources.push_back(std::move(archive));
    const auto seed = NativeRecord::create_from_endpoint(
        N(N::Object{{{"schema", N(std::string("gadget.scientific.context.v1"))},
            {"original_study_publication", N(N::ByteSpan{archive_key, 0, archive_size, false})}}}),
        new_native_record_namespace(), std::move(sources),
        rfc::storage::TemporalEndpoint(source));
    ResultFields fields(seed);
    std::vector<std::uint64_t> grade_axes;
    auto view = ResultValue::object(), grades = ResultValue::array();
    auto binding = ResultValue::object();
    binding.set("kind", ResultValue::string("declared_context_source"));
    binding.set("original_study_publication", ResultValue(N(
        N::ByteSpan{archive_key, 0, archive_size, false})));
    auto lifts = ResultValue::array();
    const auto original_index = rfc::storage::StableSourceIndex(original.endpoint()->source);
    for (std::size_t k = 0; k < source.grade_count(); ++k) {
        const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};
        const auto size = source.grade(grade).cell_count();
        require(to_original.at(k).size() == size, "source coordinate lift size");
        const auto axis = fields.source_grade(k);
        grade_axes.push_back(axis);
        auto declaration = ResultValue::object(), identities = ResultValue::array();
        declaration.set("grade", ResultValue::count(k));
        declaration.set("carrier", ResultValue::count(axis));
        declaration.set("metric", ResultFields::stage_matrix(
            study.tower().metrics().at(k).materialize(), axis, axis,
            rfc::storage::TemporalTensorRole::Metric));
        std::vector<std::size_t> coordinates;
        std::set<std::size_t> unique;
        for (std::size_t i = 0; i < size; ++i) {
            const auto cell = rfc::storage::stable_cell_key(source,
                {grade, static_cast<std::uint32_t>(i)});
            const auto original_cell = original_index.cell(cell);
            require(original_cell && original_cell->index == to_original.at(k).at(i) &&
                        unique.insert(original_cell->index).second,
                    "source lift must name distinct actual original identities");
            identities.push(ResultValue::count(cell.key));
            coordinates.push_back(original_cell->index);
        }
        declaration.set("identities", std::move(identities));
        auto lift = ResultValue::object();
        lift.set("grade", ResultValue::count(k)); lift.set("carrier", ResultValue::count(axis));
        lift.set("original_coordinates", fields.indices(coordinates, axis));
        lifts.push(std::move(lift)); grades.push(std::move(declaration));
    }
    binding.set("lifts", std::move(lifts)); fields.bind(std::move(binding));
    view.set("grades", std::move(grades));
    return {fields.finish(std::move(view)), std::move(grade_axes)};
}
} // namespace gadget::research_detail
