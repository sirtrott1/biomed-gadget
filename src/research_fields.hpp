#pragma once
#include "scientific_detail.hpp"
#include <gadget/scientific.hpp>
#include <map>
#include <optional>
#include <set>

namespace gadget::research_detail {
using namespace rfc;
using namespace scientific_detail;
using Indices = std::vector<std::size_t>;
inline std::vector<std::string> names(const Json &a, std::size_t limit, bool empty = false) {
    require(json_object_is_type(a.raw(), json_type_array) && a.size() <= limit &&
                (empty || a.size()),
            "research field-list budget");
    std::vector<std::string> out;
    std::set<std::string> seen;
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto s = text(a.at(i), source_byte_limit());
        require(seen.insert(s).second, "duplicate research field");
        out.push_back(s);
    }
    return out;
}
struct Inputs {
    std::shared_ptr<const storage::TemporalEndpoint> endpoint;
    std::map<std::string, const storage::TemporalTensor *> fields;
    std::map<std::string, Json> meaning;
    std::vector<std::string> entities;
    std::optional<std::uint64_t> carrier;
    core::Grade grade{0};
    explicit Inputs(const ScientificStudy &study, const std::set<std::string> &requested,
                    bool finite_topology_profile,
                    std::optional<core::Grade> required_grade = core::Grade{0})
        : endpoint(study.endpoint()) {
        const auto declaration = study.declaration();
        const auto tensors = storage::temporal_endpoint_detail::tensors(*endpoint);
        const auto carriers = storage::temporal_endpoint_detail::carriers(*endpoint);
        const auto df = declaration.at("fields");
        for (std::size_t i = 0; i < df.size(); ++i) {
            auto f = df.at(i);
            const auto id = f.at("id").text();
            if (!requested.contains(id))
                continue;
            const auto *t = tensors.at(key(f.at("key")));
            require(t->output.has_value(), "research field carrier missing");
            const auto *c = carriers.at(*t->output);
            require(c->kind == storage::TemporalCarrierKind::SourceGrade &&
                        (!required_grade || c->grade == *required_grade) &&
                        c->domain.extents().size() == 1,
                    "research requires scalar fields on the declared source grade");
            require(!carrier || *carrier == *t->output,
                    "research fields require the identical native carrier");
            carrier = *t->output;
            grade = c->grade;
            require(f.at("value_type").text() != "categorical",
                    "categorical codes require declared indicator fields");
            fields.emplace(id, t);
            meaning.emplace(id, f);
        }
        require(fields.size() == requested.size(), "unknown research field");
        const auto n = endpoint->source.grade(grade).cell_count();
        require(n > 0 && (!finite_topology_profile || n <= 4096), "research population budget");
        for (std::size_t i = 0; i < n; ++i)
            entities.push_back(study.cell_identifier(grade, i));
        for (const auto &[id, t] : fields) {
            (void)id;
            require(t->tensor.domain().logical_size() == n, "research field size");
        }
    }
    Rational value(const std::string &id, std::size_t row) const {
        return *fields.at(id)->tensor.state(row);
    }
    void indicator(const std::string &id) const {
        require(meaning.at(id).at("value_type").text() == "indicator",
                "selection must name a declared indicator field");
    }
    Indices select(const std::vector<std::string> &complete,
                   const std::vector<std::string> &selection, bool allow_empty = false) const {
        for (const auto &s : selection)
            indicator(s);
        Indices rows;
        for (std::size_t i = 0; i < entities.size(); ++i) {
            bool admitted = true;
            for (const auto &s : selection)
                admitted =
                    admitted && fields.at(s)->tensor.support().contains(i) && value(s, i) == 1;
            for (const auto &s : complete)
                admitted = admitted && fields.at(s)->tensor.support().contains(i);
            if (admitted)
                rows.push_back(i);
        }
        require(allow_empty || !rows.empty(), "research selection has no jointly observed rows");
        return rows;
    }
    Vector values(const std::string &id, const Indices &rows) const {
        Vector x;
        for (auto i : rows)
            x.push_back(value(id, i));
        return x;
    }
    Json identity(const Indices &rows) const {
        auto out = Json::array();
        for (auto i : rows)
            out.push(Json::string(entities[i]));
        return out;
    }
};
} // namespace gadget::research_detail
