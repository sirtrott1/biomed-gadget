#pragma once
#include "scientific_temporal_prepare.hpp"
#include <tuple>
#include "scientific_result_fields.hpp"

namespace gadget::research_detail {
// A certificate representation of complete maps. Numerical execution remains
// in the canonical SDK; this registry introduces no evaluator or solve kernel.
template<class Value> class BasicFactoredTemporalMaps {
    using Matrix = rfc::numeric::ExactMatrix;
    using Grade = rfc::core::Grade;
    std::vector<scientific_detail::PreparedTemporalStudy *> endpoints_;
    std::vector<Value> nodes_;
    scientific_detail::ResultFields *publication_{};
    std::vector<std::vector<std::uint64_t>> axes_;
    std::function<void(const Matrix &)> admit_matrix_;
    std::function<void(const Rational &)> admit_scalar_;
    static Value count(std::size_t value);
    static std::size_t dimension(const Value &map, const char *key);
    static void node_set(Value &node, std::string name, Value value);
    std::uint64_t axis(std::size_t endpoint, Grade grade) const;
    Value typed(Value map, std::uint64_t input, std::uint64_t output);
    std::map<std::tuple<std::size_t, std::uint16_t, std::string>, Value> bound_;
    Value append(Value node, std::size_t rows, std::size_t columns);
    Value bound(std::size_t endpoint, Grade grade, const std::string &role, Value map);

  public:
    explicit BasicFactoredTemporalMaps(
        std::vector<scientific_detail::PreparedTemporalStudy *> endpoints,
        scientific_detail::ResultFields *publication = nullptr,
        std::vector<std::vector<std::uint64_t>> axes = {},
        std::function<void(const Matrix &)> admit_matrix = {},
        std::function<void(const Rational &)> admit_scalar = {});
    static bool requested(const Json &query);
    Value matrix(const Matrix &value, std::optional<std::uint64_t> input = {},
                 std::optional<std::uint64_t> output = {});
    Value identity(std::size_t dimension, std::optional<std::uint64_t> carrier = {});
    Value transpose(const Value &map);
    Value compose(const std::vector<Value> &application_order);
    Value sum(const std::vector<Value> &maps);
    Value difference(const Value &left, const Value &right);
    Value scale(const Value &map, Rational coefficient);
    Value solve(const Value &system);
    Value metric(std::size_t endpoint, Grade grade);
    Value boundary(std::size_t endpoint, Grade upper);
    Value adjoint(std::size_t endpoint, Grade upper);
    Value hodge(std::size_t endpoint, Grade grade);
    Value green(std::size_t endpoint, Grade grade, Rational lambda);
    Value sector(std::size_t endpoint, Grade grade, std::size_t sector);
    Value defect(const Value &old_map, const Value &new_map, const Value &input_map,
                const Value &output_map);
    Value archive() const;
};
using FactoredTemporalMaps = BasicFactoredTemporalMaps<Json>;
using NativeFactoredTemporalMaps = BasicFactoredTemporalMaps<scientific_detail::ResultValue>;
} // namespace gadget::research_detail
