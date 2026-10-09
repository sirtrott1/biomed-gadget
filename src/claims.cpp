#include "claim_detail.hpp"
#include "scientific_detail.hpp"
#include <limits>

namespace gadget::claim_detail {
using namespace scientific_detail;
Json member_path(Json value, const Json &path) {
    source_array(path, "certificate path");
    require(path.size(), "certificate path requires components");
    for (std::size_t i = 0; i < path.size(); ++i) {
        const auto part = path.at(i);
        if (json_object_is_type(part.raw(), json_type_string))
            value = value.at(text(part, source_byte_limit()));
        else
            value = value.at(
                count(part, value.size() ? value.size() - 1 : 0, "certificate array coordinate"));
    }
    return value;
}
Rational exact_claim_scalar(const Json &v) {
    if (json_object_is_type(v.raw(), json_type_object))
        return exact_claim_scalar(v.at("rational"));
    const auto q = exact_rational(text(v, source_byte_limit()));
    return q;
}
Json literal_value(const Json &value, const std::string &type) {
    if (type == "rational")
        return Json::string(exact_claim_scalar(value).get_str());
    (void)compare_claim(value, value, type, "equal");
    if (type == "boolean")
        return Json::boolean(value.boolean());
    if (type == "count" || type == "integer")
        return Json::integer(value.integer());
    return Json::string(value.text());
}
bool compare_claim(const Json &actual, const Json &expected, const std::string &type,
                   const std::string &comparison) {
    int ordering = 0;
    if (type == "boolean") {
        require(json_object_is_type(actual.raw(), json_type_boolean) &&
                    json_object_is_type(expected.raw(), json_type_boolean),
                "boolean certificate claim type");
        ordering = int(actual.boolean()) - int(expected.boolean());
        require(comparison == "equal" || comparison == "not_equal", "boolean comparison");
    } else if (type == "count" || type == "integer") {
        require(json_object_is_type(actual.raw(), json_type_int) &&
                    json_object_is_type(expected.raw(), json_type_int),
                "count certificate claim type");
        const auto a = actual.integer(), b = expected.integer();
        require(type != "count" || (a >= 0 && b >= 0), "count claims require nonnegative integers");
        ordering = (a > b) - (a < b);
    } else if (type == "text") {
        const auto a = actual.text(), b = expected.text();
        require(a.find('\0') == std::string::npos && b.find('\0') == std::string::npos,
                "literal text contains NUL");
        ordering = (a > b) - (a < b);
        require(comparison == "equal" || comparison == "not_equal", "text comparison");
    } else {
        require(type == "rational", "certificate claim value type");
        const auto a = exact_claim_scalar(actual), b = exact_claim_scalar(expected);
        ordering = (a > b) - (a < b);
    }
    if (comparison == "equal")
        return ordering == 0;
    if (comparison == "not_equal")
        return ordering != 0;
    if (comparison == "greater")
        return ordering > 0;
    if (comparison == "greater_equal")
        return ordering >= 0;
    if (comparison == "less")
        return ordering < 0;
    require(comparison == "less_equal", "certificate comparison operator");
    return ordering <= 0;
}
void admit_condition(const Json &condition) {
    members(condition, {"path", "value_type", "comparison", "value"});
    const auto path = condition.at("path");
    source_array(path, "condition path");
    require(path.size(), "condition path requires components");
    for (std::size_t i = 0; i < path.size(); ++i) {
        const auto part = path.at(i);
        if (json_object_is_type(part.raw(), json_type_string))
            (void)text(part, source_byte_limit());
        else
            (void)count(part, std::numeric_limits<std::size_t>::max(),
                        "condition array coordinate");
    }
    (void)compare_claim(condition.at("value"), condition.at("value"),
                        condition.at("value_type").text(), condition.at("comparison").text());
}
Json check_literal_claim(const Json &value, const Json &claim) {
    admit_condition(claim);
    const auto actual = member_path(value, claim.at("path"));
    const bool agrees = compare_claim(actual, claim.at("value"), claim.at("value_type").text(),
                                      claim.at("comparison").text());
    auto out = Json::object();
    out.set("claim", claim);
    out.set("observed", actual);
    out.set("agrees", Json::boolean(agrees));
    return seal(out, "check_sha256");
}
void verify_seal(const Json &value, const std::string &name) {
    auto body = Json::parse(value.dump());
    const auto supplied = body.at(name).text();
    json_object_object_del(body.raw(), name.c_str());
    require(supplied == "sha256:" + sha256(body.dump()), "invalid retained seal: " + name);
}
} // namespace gadget::claim_detail
