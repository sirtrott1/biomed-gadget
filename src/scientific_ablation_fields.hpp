#pragma once
#include "scientific_ablation_source.hpp"
#include "scientific_relational_field.hpp"

namespace gadget::research_detail {
struct AblationSourceAxes {
    NativeRecord snapshot;
    std::vector<std::uint64_t> rows;
    std::vector<std::uint64_t> lattice;
    scientific_detail::ResultValue reference;
};
class AblationFields {
    const NativeRecord &original_;
    std::map<const NativeStudy *, AblationSourceAxes> owners_;
    std::map<std::pair<const void *, std::uint64_t>, std::vector<std::uint64_t>> stalk_owners_;
  public:
    scientific_detail::ResultFields &fields;
    explicit AblationFields(scientific_detail::ResultFields &publication, const NativeRecord &original)
        : original_(original), fields(publication) {}
    const AblationSourceAxes &source(const NativeStudy &study) {
        auto found = owners_.find(&study);
        if (found != owners_.end()) return found->second;
        std::vector<std::vector<std::uint32_t>> coordinates;
        const rfc::storage::StableSourceIndex index(original_.endpoint()->source);
        for (std::size_t k = 0; k < study.source().grade_count(); ++k) {
            std::vector<std::uint32_t> row;
            const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};
            for (std::size_t i = 0; i < study.source().grade(grade).cell_count(); ++i) {
                const auto key = rfc::storage::stable_cell_key(study.source(),
                    {grade, static_cast<std::uint32_t>(i)});
                const auto cell = index.cell(key);
                scientific_detail::require(cell.has_value(), "ablation source identity absent from original owner");
                row.push_back(cell->index);
            }
            coordinates.push_back(std::move(row));
        }
        auto snapshot = ablation_source_publication(study, original_, coordinates);
        std::vector<std::uint64_t> rows, lattice;
        auto binding = scientific_detail::ResultValue::object();
        binding.set("kind", scientific_detail::ResultValue::string("actual_ablation_source_axes"));
        auto reference = fields.publication(snapshot.publication);
        binding.set("source", reference);
        auto grades = scientific_detail::ResultValue::array();
        for (std::size_t k = 0; k < snapshot.grade_axes.size(); ++k) {
            const auto size = coordinates[k].size();
            rows.push_back(fields.axis(size));
            lattice.push_back(fields.axis(size));
            auto grade = scientific_detail::ResultValue::object();
            grade.set("grade", scientific_detail::ResultValue::count(k));
            grade.set("source_carrier", scientific_detail::ResultValue::count(snapshot.grade_axes[k]));
            grade.set("coordinate_carrier", scientific_detail::ResultValue::count(rows.back()));
            grade.set("lattice_carrier", scientific_detail::ResultValue::count(lattice.back()));
            std::vector<std::size_t> lift(coordinates[k].begin(), coordinates[k].end());
            grade.set("original_coordinates", fields.indices(lift, rows.back()));
            rfc::numeric::ExactMatrix realization(size, size);
            const auto &scales = study.source().grade(rfc::core::Grade{static_cast<std::uint16_t>(k)}).cell_scales;
            for (std::size_t i = 0; i < size; ++i) realization.set(i, i, scales.at(i));
            grade.set("lattice_realization", scientific_detail::ResultFields::stage_matrix(
                realization, lattice.back(), rows.back()));
            grades.push(std::move(grade));
        }
        binding.set("grades", std::move(grades)); fields.bind(std::move(binding));
        return owners_.emplace(&study, AblationSourceAxes{std::move(snapshot.publication),
            std::move(rows), std::move(lattice), std::move(reference)}).first->second;
    }
    const std::vector<std::uint64_t> &stalks(const NativeRecord &owner,
        const std::vector<std::uint64_t> &source_axes) {
        const auto id = std::pair<const void *, std::uint64_t>{owner.endpoint().get(),
            owner.publication_history()->last_sequence()};
        auto found = stalk_owners_.find(id);
        if (found != stalk_owners_.end()) return found->second;
        std::vector<std::uint64_t> axes;
        auto binding = scientific_detail::ResultValue::object(), grades = scientific_detail::ResultValue::array();
        binding.set("kind", scientific_detail::ResultValue::string("actual_section_coordinate_axes"));
        binding.set("source", fields.publication(owner));
        for (std::size_t k = 0; k < source_axes.size(); ++k) {
            const auto carrier = std::find_if(owner.endpoint()->carriers.begin(), owner.endpoint()->carriers.end(),
                [&](const auto &value) {return value.key == source_axes[k];});
            scientific_detail::require(carrier != owner.endpoint()->carriers.end(), "section source carrier absent");
            const auto axis = fields.axis(carrier->domain.logical_size()); axes.push_back(axis);
            auto grade = scientific_detail::ResultValue::object();
            grade.set("grade", scientific_detail::ResultValue::count(k));
            grade.set("source_carrier", scientific_detail::ResultValue::count(source_axes[k]));
            grade.set("coordinate_carrier", scientific_detail::ResultValue::count(axis));
            grades.push(std::move(grade));
        }
        binding.set("grades", std::move(grades)); fields.bind(std::move(binding));
        return stalk_owners_.emplace(id, std::move(axes)).first->second;
    }
    scientific_detail::ResultValue scalar(const Rational &value) const {
        return fields.stage_scalar(value);
    }
    static scientific_detail::ResultValue vector(const Vector &value, std::uint64_t axis) {
        return scientific_detail::ResultFields::stage_vector(value, axis);
    }
    static scientific_detail::ResultValue matrix(const rfc::numeric::ExactMatrix &value,
                                                std::uint64_t input, std::uint64_t output) {
        return scientific_detail::ResultFields::stage_matrix(value, input, output);
    }
};
} // namespace gadget::research_detail
