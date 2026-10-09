#pragma once
#include "scientific_temporal_prepare.hpp"
#include <tuple>

namespace gadget::research_detail {
// A certificate representation of complete maps. Numerical execution remains
// in the canonical SDK; this registry introduces no evaluator or solve kernel.
class FactoredTemporalMaps {
    using Matrix = rfc::numeric::ExactMatrix;
    using Grade = rfc::core::Grade;
    std::vector<scientific_detail::PreparedTemporalStudy *> endpoints_;
    Json nodes_ = Json::array();
    std::map<std::tuple<std::size_t, std::uint16_t, std::string>, Json> bound_;
    Json append(Json node, std::size_t rows, std::size_t columns);
    Json bound(std::size_t endpoint, Grade grade, const std::string &role, Json map);

  public:
    explicit FactoredTemporalMaps(
        std::vector<scientific_detail::PreparedTemporalStudy *> endpoints);
    static bool requested(const Json &query);
    Json matrix(const Matrix &value);
    Json identity(std::size_t dimension);
    Json transpose(const Json &map);
    Json compose(const std::vector<Json> &application_order);
    Json sum(const std::vector<Json> &maps);
    Json difference(const Json &left, const Json &right);
    Json scale(const Json &map, Rational coefficient);
    Json solve(const Json &system);
    Json metric(std::size_t endpoint, Grade grade);
    Json boundary(std::size_t endpoint, Grade upper);
    Json adjoint(std::size_t endpoint, Grade upper);
    Json hodge(std::size_t endpoint, Grade grade);
    Json green(std::size_t endpoint, Grade grade, Rational lambda);
    Json sector(std::size_t endpoint, Grade grade, std::size_t sector);
    Json defect(const Json &old_map, const Json &new_map, const Json &input_map,
                const Json &output_map);
    Json archive() const;
};
} // namespace gadget::research_detail
