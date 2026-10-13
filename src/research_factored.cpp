#include "research_factored.hpp"

namespace gadget::research_detail {
using namespace scientific_detail;
template<class Value> Value BasicFactoredTemporalMaps<Value>::count(std::size_t value) {
    if constexpr (std::is_same_v<Value, ResultValue>) return ResultValue::count(value);
    else return cardinality(value);
}
template<class Value> std::size_t BasicFactoredTemporalMaps<Value>::dimension(const Value &map, const char *key) {
    if constexpr (std::is_same_v<Value, ResultValue>) return map.at(key).count_value();
    else return static_cast<std::size_t>(map.at(key).integer());
}
template<class Value> void BasicFactoredTemporalMaps<Value>::node_set(Value &node, std::string name, Value value) {
    if constexpr (std::is_same_v<Value, ResultValue>) {
        try {
            if (node.at("schema").text() == "gadget.scientific.view.overlay.v1") {
                auto additions = node.at("additions"); additions.set(std::move(name), std::move(value));
                node.set("additions", std::move(additions)); return;
            }
        } catch (const std::invalid_argument &) {}
    }
    node.set(std::move(name), std::move(value));
}
template<class Value> std::uint64_t BasicFactoredTemporalMaps<Value>::axis(std::size_t endpoint, Grade grade) const {
    if constexpr (std::is_same_v<Value, ResultValue>) return axes_.at(endpoint).at(grade.value);
    else { (void)endpoint; (void)grade; return 0; }
}
template<class Value> Value BasicFactoredTemporalMaps<Value>::typed(Value map, std::uint64_t input, std::uint64_t output) {
    if constexpr (std::is_same_v<Value, ResultValue>) {
        map.set("input_carrier", count(input)); map.set("output_carrier", count(output));
        auto &node = nodes_.at(dimension(map, "expression"));
        node_set(node, "input_carrier", count(input)); node_set(node, "output_carrier", count(output));
    }
    return map;
}

template<class Value>
BasicFactoredTemporalMaps<Value>::BasicFactoredTemporalMaps(std::vector<PreparedTemporalStudy *> endpoints, ResultFields *publication,
    std::vector<std::vector<std::uint64_t>> axes, std::function<void(const Matrix &)> admit_matrix,
    std::function<void(const Rational &)> admit_scalar)
    : endpoints_(std::move(endpoints)), publication_(publication), axes_(std::move(axes)),
      admit_matrix_(std::move(admit_matrix)), admit_scalar_(std::move(admit_scalar)) {
    if constexpr (std::is_same_v<Value, ResultValue>)
        require(publication_ && axes_.size() == endpoints_.size(), "native factored maps require actual endpoint axes");
    for (const auto *owner : endpoints_)
        require(owner, "factored map requires an admitted endpoint owner");
}
template<class Value>
bool BasicFactoredTemporalMaps<Value>::requested(const Json &query) {
    if (!query.has("map_output"))
        return false;
    require(json_object_is_type(query.at("map_output").raw(), json_type_string),
            "map_output must be a string");
    const auto mode = query.at("map_output").text();
    require(mode == "factored" || mode == "enumerated", "unknown complete-map output");
    return mode == "factored";
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::append(Value node, std::size_t rows, std::size_t columns) {
    node_set(node, "rows", count(rows));
    node_set(node, "columns", count(columns));
    auto map = Value::object();
    map.set("rows", count(rows));
    map.set("columns", count(columns));
    map.set("expression", count(nodes_.size()));
    nodes_.push_back(std::move(node));
    return map;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::bound(std::size_t endpoint, Grade grade, const std::string &role,
                                 Value map) {
    auto &node = nodes_.at(dimension(map, "expression"));
    node_set(node, "endpoint", count(endpoint));
    node_set(node, "grade", count(grade.value));
    node_set(node, "role", Value::string(role));
    bound_.emplace(std::make_tuple(endpoint, grade.value, role), map);
    return map;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::matrix(const Matrix &value, std::optional<std::uint64_t> input,
    std::optional<std::uint64_t> output) {
    if constexpr (std::is_same_v<Value, ResultValue>) {
        if (admit_matrix_) admit_matrix_(value);
        const auto domain = input ? *input : publication_->axis(value.cols());
        const auto codomain = output ? *output : publication_->axis(value.rows());
        auto field = ResultFields::stage_matrix(value, domain, codomain);
        field.set("display", ResultValue::string("sparse_matrix"));
        auto properties = ResultValue::object();
        properties.set("kind", ResultValue::string("matrix"));
        auto node = ResultFields::overlay(std::move(field), std::move(properties));
        auto map = append(std::move(node), value.rows(), value.cols());
        return typed(std::move(map), domain, codomain);
    } else {
    auto node = Value::object(), entries = Value::array();
    node_set(node, "kind", Value::string("matrix"));
    for (const auto &[p, v] : value.entries()) {
        auto e = Value::object();
        e.set("row", count(p.first));
        e.set("column", count(p.second));
        e.set("value", Value::string(v.get_str()));
        entries.push(e);
    }
    node_set(node, "entries", entries);
    return append(node, value.rows(), value.cols());
    }
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::identity(std::size_t n, std::optional<std::uint64_t> carrier) {
    auto node = Value::object();
    node_set(node, "kind", Value::string("identity"));
    auto map = append(std::move(node), n, n);
    if constexpr (std::is_same_v<Value, ResultValue>) {
        const auto key = carrier ? *carrier : publication_->axis(n);
        return typed(std::move(map), key, key);
    }
    return map;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::transpose(const Value &map) {
    auto node = Value::object();
    node_set(node, "kind", Value::string("transpose"));
    node_set(node, "input", map.at("expression"));
    auto result = append(node, dimension(map, "columns"), dimension(map, "rows"));
    if constexpr (std::is_same_v<Value, ResultValue>)
        return typed(std::move(result), dimension(map, "output_carrier"), dimension(map, "input_carrier"));
    return result;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::compose(const std::vector<Value> &order) {
    require(!order.empty(), "factored composition requires a typed map");
    auto node = Value::object(), refs = Value::array();
    node_set(node, "kind", Value::string("compose"));
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (i)
            require(dimension(order[i - 1], "rows") == dimension(order[i], "columns"),
                    "factored composition carrier mismatch");
        if constexpr (std::is_same_v<Value, ResultValue>)
            if (i) require(dimension(order[i - 1], "output_carrier") == dimension(order[i], "input_carrier"),
                           "native factored composition requires identical actual carriers");
        refs.push(order[i].at("expression"));
    }
    node_set(node, "application_order", refs);
    auto result = append(node, dimension(order.back(), "rows"), dimension(order.front(), "columns"));
    if constexpr (std::is_same_v<Value, ResultValue>)
        return typed(std::move(result), dimension(order.front(), "input_carrier"), dimension(order.back(), "output_carrier"));
    return result;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::sum(const std::vector<Value> &maps) {
    require(!maps.empty(), "factored sum requires a typed map");
    auto node = Value::object(), refs = Value::array();
    node_set(node, "kind", Value::string("sum"));
    const auto rows = dimension(maps.front(), "rows"), cols = dimension(maps.front(), "columns");
    for (const auto &map : maps) {
        require(dimension(map, "rows") == rows && dimension(map, "columns") == cols,
                "factored sum carrier mismatch");
        if constexpr (std::is_same_v<Value, ResultValue>)
            require(dimension(map, "input_carrier") == dimension(maps.front(), "input_carrier") &&
                    dimension(map, "output_carrier") == dimension(maps.front(), "output_carrier"),
                    "native factored sum requires identical actual carriers");
        refs.push(map.at("expression"));
    }
    node_set(node, "inputs", refs);
    auto result = append(node, rows, cols);
    if constexpr (std::is_same_v<Value, ResultValue>)
        return typed(std::move(result), dimension(maps.front(), "input_carrier"), dimension(maps.front(), "output_carrier"));
    return result;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::scale(const Value &map, Rational coefficient) {
    auto node = Value::object();
    node_set(node, "kind", Value::string("scale"));
    node_set(node, "input", map.at("expression"));
    if constexpr (std::is_same_v<Value, ResultValue>) {
        if (admit_scalar_) admit_scalar_(coefficient);
        node_set(node, "coefficient", publication_->stage_scalar(coefficient));
    }
    else node_set(node, "coefficient", Value::string(coefficient.get_str()));
    auto result = append(node, dimension(map, "rows"), dimension(map, "columns"));
    if constexpr (std::is_same_v<Value, ResultValue>)
        return typed(std::move(result), dimension(map, "input_carrier"), dimension(map, "output_carrier"));
    return result;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::difference(const Value &left, const Value &right) {
    return sum({left, scale(right, Rational{-1})});
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::solve(const Value &system) {
    require(dimension(system, "rows") == dimension(system, "columns"),
            "factored solve requires a square complete system");
    auto node = Value::object();
    node_set(node, "kind", Value::string("solve"));
    node_set(node, "system", system.at("expression"));
    auto result = append(node, dimension(system, "rows"), dimension(system, "columns"));
    if constexpr (std::is_same_v<Value, ResultValue>)
        return typed(std::move(result), dimension(system, "output_carrier"), dimension(system, "input_carrier"));
    return result;
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::metric(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("metric"));
    if (bound_.contains(key))
        return bound_.at(key);
    return bound(t, k, "metric",
                 matrix(endpoints_.at(t)->study().tower().metrics().at(k.value).materialize(), axis(t,k), axis(t,k)));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::boundary(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("boundary"));
    if (bound_.contains(key))
        return bound_.at(key);
    require(k.value, "factored boundary at grade zero");
    return bound(
        t, k, "boundary",
        matrix(rfc::field::materialize_boundary_matrix(endpoints_.at(t)->study().source(), k),
               axis(t,k), axis(t,Grade{static_cast<std::uint16_t>(k.value-1u)})));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::adjoint(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("adjoint"));
    if (bound_.contains(key))
        return bound_.at(key);
    require(k.value, "factored adjoint requires an upper grade");
    const Grade lower{static_cast<std::uint16_t>(k.value - 1u)};
    return bound(t, k, "adjoint",
                 compose({metric(t, lower), transpose(boundary(t, k)), solve(metric(t, k))}));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::hodge(std::size_t t, Grade k) {
    const auto key = std::make_tuple(t, k.value, std::string("hodge"));
    if (bound_.contains(key))
        return bound_.at(key);
    const auto &study = endpoints_.at(t)->study();
    std::vector<Value> excursions;
    if (k.value)
        excursions.push_back(compose({boundary(t, k), adjoint(t, k)}));
    if (k.value < study.source().max_grade().value) {
        const Grade upper{static_cast<std::uint16_t>(k.value + 1u)};
        excursions.push_back(compose({adjoint(t, upper), boundary(t, upper)}));
    }
    const auto n = study.tower().metrics().at(k.value).size();
    return bound(t, k, "hodge", excursions.empty() ? matrix(Matrix(n, n),axis(t,k),axis(t,k)) : sum(excursions));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::green(std::size_t t, Grade k, Rational lambda) {
    lambda.canonicalize();
    require(lambda > 0, "factored Green parameter must be positive");
    const auto role = "green:" + lambda.get_str();
    const auto key = std::make_tuple(t, k.value, role);
    if (bound_.contains(key))
        return bound_.at(key);
    const auto n = endpoints_.at(t)->study().tower().metrics().at(k.value).size();
    return bound(t, k, role, solve(sum({identity(n,axis(t,k)), scale(hodge(t, k), lambda)})));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::sector(std::size_t t, Grade k, std::size_t s) {
    require(s < 3, "unknown factored Hodge sector");
    const auto role = std::string("sector:") + std::to_string(s);
    const auto key = std::make_tuple(t, k.value, role);
    if (bound_.contains(key))
        return bound_.at(key);
    const auto &prepared = endpoints_.at(t)->sectors(k);
    Value map;
    if (s == 2)
        map = difference(identity(prepared.metric.size(),axis(t,k)), sum({sector(t, k, 0), sector(t, k, 1)}));
    else {
        std::optional<std::uint64_t> coefficient;
        if constexpr (std::is_same_v<Value, ResultValue>) coefficient=publication_->axis(prepared.projections[s].frame.cols());
        const auto frame = matrix(prepared.projections[s].frame,coefficient,axis(t,k)), dual = transpose(frame);
        const auto gram = compose({frame, metric(t, k), dual});
        map = compose({metric(t, k), dual, solve(gram), frame});
    }
    return bound(t, k, role, map);
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::defect(const Value &old_map, const Value &new_map, const Value &input_map,
                                  const Value &output_map) {
    return difference(compose({input_map, new_map}), compose({old_map, output_map}));
}
template<class Value>
Value BasicFactoredTemporalMaps<Value>::archive() const {
    auto out = Value::object();
    out.set("format", Value::string("gadget-factored-maps-v1"));
    auto nodes = Value::array();
    for (const auto &node : nodes_) nodes.push(node);
    out.set("nodes", std::move(nodes));
    return out;
}
template class BasicFactoredTemporalMaps<Json>;
template class BasicFactoredTemporalMaps<ResultValue>;
} // namespace gadget::research_detail
