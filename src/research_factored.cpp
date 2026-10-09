#include "research_factored.hpp"

namespace gadget::research_detail {
using namespace scientific_detail;
namespace {
std::size_t dimension(const Json &map, const char *key) {
    return static_cast<std::size_t>(map.at(key).integer());
}
Json rational(Rational value) {
    value.canonicalize();
    return Json::string(value.get_str());
}
} // namespace
FactoredTemporalMaps::FactoredTemporalMaps(std::vector<PreparedTemporalStudy *> endpoints)
    : endpoints_(std::move(endpoints)) {
    for (const auto *owner : endpoints_)
        require(owner, "factored map requires an admitted endpoint owner");
}
bool FactoredTemporalMaps::requested(const Json &query) {
    if (!query.has("map_output"))
        return false;
    require(json_object_is_type(query.at("map_output").raw(), json_type_string),
            "map_output must be a string");
    const auto mode = query.at("map_output").text();
    require(mode == "factored" || mode == "enumerated", "unknown complete-map output");
    return mode == "factored";
}
Json FactoredTemporalMaps::append(Json node, std::size_t rows, std::size_t columns) {
    node.set("rows", cardinality(rows));
    node.set("columns", cardinality(columns));
    auto map = Json::object();
    map.set("rows", cardinality(rows));
    map.set("columns", cardinality(columns));
    map.set("expression", cardinality(nodes_.size()));
    nodes_.push(node);
    return map;
}
Json FactoredTemporalMaps::bound(std::size_t endpoint, Grade grade, const std::string &role,
                                 Json map) {
    auto node = nodes_.at(dimension(map, "expression"));
    node.set("endpoint", cardinality(endpoint));
    node.set("grade", cardinality(grade.value));
    node.set("role", Json::string(role));
    bound_.emplace(std::make_tuple(endpoint, grade.value, role), map);
    return map;
}
Json FactoredTemporalMaps::matrix(const Matrix &value) {
    auto node = Json::object(), entries = Json::array();
    node.set("kind", Json::string("matrix"));
    for (const auto &[p, v] : value.entries()) {
        auto e = Json::object();
        e.set("row", cardinality(p.first));
        e.set("column", cardinality(p.second));
        e.set("value", rational(v));
        entries.push(e);
    }
    node.set("entries", entries);
    return append(node, value.rows(), value.cols());
}
Json FactoredTemporalMaps::identity(std::size_t n) {
    auto node = Json::object();
    node.set("kind", Json::string("identity"));
    return append(node, n, n);
}
Json FactoredTemporalMaps::transpose(const Json &map) {
    auto node = Json::object();
    node.set("kind", Json::string("transpose"));
    node.set("input", map.at("expression"));
    return append(node, dimension(map, "columns"), dimension(map, "rows"));
}
Json FactoredTemporalMaps::compose(const std::vector<Json> &order) {
    require(!order.empty(), "factored composition requires a typed map");
    auto node = Json::object(), refs = Json::array();
    node.set("kind", Json::string("compose"));
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (i)
            require(dimension(order[i - 1], "rows") == dimension(order[i], "columns"),
                    "factored composition carrier mismatch");
        refs.push(order[i].at("expression"));
    }
    node.set("application_order", refs);
    return append(node, dimension(order.back(), "rows"), dimension(order.front(), "columns"));
}
Json FactoredTemporalMaps::sum(const std::vector<Json> &maps) {
    require(!maps.empty(), "factored sum requires a typed map");
    auto node = Json::object(), refs = Json::array();
    node.set("kind", Json::string("sum"));
    const auto rows = dimension(maps.front(), "rows"), cols = dimension(maps.front(), "columns");
    for (const auto &map : maps) {
        require(dimension(map, "rows") == rows && dimension(map, "columns") == cols,
                "factored sum carrier mismatch");
        refs.push(map.at("expression"));
    }
    node.set("inputs", refs);
    return append(node, rows, cols);
}
Json FactoredTemporalMaps::scale(const Json &map, Rational coefficient) {
    auto node = Json::object();
    node.set("kind", Json::string("scale"));
    node.set("input", map.at("expression"));
    node.set("coefficient", rational(coefficient));
    return append(node, dimension(map, "rows"), dimension(map, "columns"));
}
Json FactoredTemporalMaps::difference(const Json &left, const Json &right) {
    return sum({left, scale(right, Rational{-1})});
}
Json FactoredTemporalMaps::solve(const Json &system) {
    require(dimension(system, "rows") == dimension(system, "columns"),
            "factored solve requires a square complete system");
    auto node = Json::object();
    node.set("kind", Json::string("solve"));
    node.set("system", system.at("expression"));
    return append(node, dimension(system, "rows"), dimension(system, "columns"));
}
Json FactoredTemporalMaps::metric(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("metric"));
    if (bound_.contains(key))
        return bound_.at(key);
    return bound(t, k, "metric",
                 matrix(endpoints_.at(t)->study().tower().metrics().at(k.value).materialize()));
}
Json FactoredTemporalMaps::boundary(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("boundary"));
    if (bound_.contains(key))
        return bound_.at(key);
    require(k.value, "factored boundary at grade zero");
    return bound(
        t, k, "boundary",
        matrix(rfc::field::materialize_boundary_matrix(endpoints_.at(t)->study().source(), k)));
}
Json FactoredTemporalMaps::adjoint(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("adjoint"));
    if (bound_.contains(key))
        return bound_.at(key);
    require(k.value, "factored adjoint requires an upper grade");
    const Grade lower{static_cast<std::uint16_t>(k.value - 1u)};
    return bound(t, k, "adjoint",
                 compose({metric(t, lower), transpose(boundary(t, k)), solve(metric(t, k))}));
}
Json FactoredTemporalMaps::hodge(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("hodge"));
    if (bound_.contains(key))
        return bound_.at(key);
    const auto &study = endpoints_.at(t)->study();
    std::vector<Json> excursions;
    if (k.value)
        excursions.push_back(compose({boundary(t, k), adjoint(t, k)}));
    if (k.value < study.source().max_grade().value) {
        const Grade upper{static_cast<std::uint16_t>(k.value + 1u)};
        excursions.push_back(compose({adjoint(t, upper), boundary(t, upper)}));
    }
    const auto n = study.tower().metrics().at(k.value).size();
    return bound(t, k, "hodge", excursions.empty() ? matrix(Matrix(n, n)) : sum(excursions));
}
Json FactoredTemporalMaps::green(std::size_t t, Grade k, Rational lambda) {
    lambda.canonicalize();
    require(lambda > 0, "factored Green parameter must be positive");
    const auto role = "green:" + lambda.get_str();
    const auto key = std::make_tuple(t, k.value, role);
    if (bound_.contains(key))
        return bound_.at(key);
    const auto n = endpoints_.at(t)->study().tower().metrics().at(k.value).size();
    return bound(t, k, role, solve(sum({identity(n), scale(hodge(t, k), lambda)})));
}
Json FactoredTemporalMaps::sector(std::size_t t, Grade k, std::size_t s) {
    require(s < 3, "unknown factored Hodge sector");
    const auto role = std::string("sector:") + std::to_string(s);
    const auto key = std::make_tuple(t, k.value, role);
    if (bound_.contains(key))
        return bound_.at(key);
    const auto &prepared = endpoints_.at(t)->sectors(k);
    Json map;
    if (s == 2)
        map = difference(identity(prepared.metric.size()), sum({sector(t, k, 0), sector(t, k, 1)}));
    else {
        const auto frame = matrix(prepared.projections[s].frame), dual = transpose(frame);
        const auto gram = compose({frame, metric(t, k), dual});
        map = compose({metric(t, k), dual, solve(gram), frame});
    }
    return bound(t, k, role, map);
}
Json FactoredTemporalMaps::defect(const Json &old_map, const Json &new_map, const Json &input_map,
                                  const Json &output_map) {
    return difference(compose({input_map, new_map}), compose({old_map, output_map}));
}
Json FactoredTemporalMaps::archive() const {
    auto out = Json::object();
    out.set("format", Json::string("gadget-factored-maps-v1"));
    out.set("nodes", nodes_);
    return out;
}
} // namespace gadget::research_detail
