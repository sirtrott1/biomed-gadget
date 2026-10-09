#include "assay_field_detail.hpp"
#include <algorithm>

namespace gadget::assay_detail {
std::vector<AssayField> load_assay_fields(std::shared_ptr<const ScientificStudy> owner,
                                          const std::vector<std::string> &selected) {
    const auto declaration = owner->declaration(), fields = declaration.at("fields");
    const auto endpoint = owner->endpoint();
    std::map<std::string, Json> meanings;
    std::map<std::uint64_t, std::size_t> tensors, carriers;
    std::map<std::uint16_t, std::map<std::uint64_t, std::string>> labels;
    for (std::size_t i = 0; i < fields.size(); ++i)
        meanings.emplace(fields.at(i).at("id").text(), fields.at(i));
    for (std::size_t i = 0; i < endpoint->tensors.size(); ++i)
        tensors.emplace(endpoint->tensors[i].key, i);
    for (std::size_t i = 0; i < endpoint->carriers.size(); ++i)
        carriers.emplace(endpoint->carriers[i].key, i);
    const auto cells = declaration.at("cells");
    for (std::size_t i = 0; i < cells.size(); ++i)
        labels[static_cast<std::uint16_t>(cells.at(i).at("grade").integer())].emplace(
            key(cells.at(i).at("key")), cells.at(i).at("id").text());
    std::map<std::uint16_t, std::pair<Json, Json>> coordinates;
    std::vector<AssayField> result;
    for (const auto &field : selected) {
        require(meanings.contains(field), "assay field not admitted");
        AssayField input;
        input.owner = owner;
        input.meaning = meanings.at(field);
        require(input.meaning.at("value_type").text() != "categorical",
                "assay field requires quantitative or explicitly declared indicator values");
        const auto tensor = key(input.meaning.at("key"));
        require(tensors.contains(tensor), "admitted assay field has no native tensor");
        const auto &t = endpoint->tensors[tensors.at(tensor)];
        require(t.output.has_value() && carriers.contains(*t.output),
                "assay field requires a native output carrier");
        const auto &carrier = endpoint->carriers[carriers.at(*t.output)];
        require(carrier.kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                    carrier.domain.extents().size() == 1,
                "assay field requires native SourceGrade coordinates");
        input.grade = carrier.grade.value;
        const auto &arena = endpoint->source.grade(carrier.grade);
        require(arena.cell_count() == t.tensor.domain().logical_size(),
                "assay field native coordinate count");
        if (!coordinates.contains(input.grade)) {
            auto ids = Json::array(), keys = Json::array();
            for (std::size_t i = 0; i < arena.cell_count(); ++i) {
                const auto native_key = input.grade == 0 ? arena.identities.at(i)
                                                         : arena.relation_identities.at(i).value;
                require(labels[input.grade].contains(native_key),
                        "assay field native identity label");
                ids.push(Json::string(labels[input.grade].at(native_key)));
                keys.push(Json::string(std::to_string(native_key)));
            }
            coordinates.emplace(input.grade, std::make_pair(std::move(ids), std::move(keys)));
        }
        input.ids = coordinates.at(input.grade).first;
        input.keys = coordinates.at(input.grade).second;
        for (std::size_t i = 0; i < arena.cell_count(); ++i) {
            const auto value = t.tensor.state(i);
            input.supported.push_back(value.has_value());
            input.values.push_back(value ? *value : Rational{0});
        }
        result.push_back(std::move(input));
    }
    return result;
}
AssayField load_assay_field(std::shared_ptr<const ScientificStudy> owner,
                            const std::string &field) {
    auto fields = load_assay_fields(std::move(owner), {field});
    return std::move(fields.front());
}
} // namespace gadget::assay_detail
