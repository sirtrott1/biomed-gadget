#include "claim_detail.hpp"
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
#include <gadget/workflow.hpp>
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
    members(condition, {"path", "value_type", "comparison", "value"}, {"coordinate"});
    if (condition.has("coordinate"))
        (void)count(condition.at("coordinate"), std::numeric_limits<std::size_t>::max(),
                    "native claim original coordinate");
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
    require(!claim.has("coordinate"), "tensor coordinates require an actual native claim owner");
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

namespace gadget {
namespace {
using namespace scientific_detail;
using N = NativeRecordValue;
using Literal = std::variant<bool, rfc::numeric::Integer, Rational, std::string>;

std::uint64_t native_count(const NativeRecord &value, const char *name) {
    require(value.kind() == NativeRecordKind::Integer, std::string(name) + " requires an integer");
    const auto n = value.integer();
    require(n >= 0 && mpz_fits_ulong_p(n.get_mpz_t()), std::string(name) + " is outside native coordinates");
    return n.get_ui();
}
Literal native_literal(const NativeRecord &value, const std::string &type, bool expected) {
    if (type == "boolean") {
        require(value.kind() == NativeRecordKind::Boolean, "native Boolean claim type");
        return value.boolean();
    }
    if (type == "integer" || type == "count") {
        require(value.kind() == NativeRecordKind::Integer, "native integer claim type");
        const auto n = value.integer();
        require(type != "count" || n >= 0, "native count must be nonnegative");
        return n;
    }
    if (type == "text") {
        require(value.kind() == NativeRecordKind::Text, "native text claim type");
        const auto s = value.text();
        require(s.find('\0') == std::string::npos, "native text claim contains NUL");
        return s;
    }
    require(type == "rational", "native claim value type");
    if (expected && value.kind() == NativeRecordKind::Text)
        return exact_rational(value.text());
    if (value.kind() == NativeRecordKind::Integer) return Rational(value.integer());
    require(value.kind() == NativeRecordKind::Rational || value.kind() == NativeRecordKind::Decimal ||
                value.kind() == NativeRecordKind::ExactNumber,
            "native rational observation requires a typed value, not numerical display text");
    return value.rational();
}
bool compare_native(const Literal &actual, const Literal &expected, const std::string &comparison) {
    require(actual.index() == expected.index(), "native claim literal types differ");
    const auto order = std::visit([&](const auto &a) {
        const auto &b = std::get<std::decay_t<decltype(a)>>(expected);
        return (a > b) - (a < b);
    }, actual);
    if (comparison == "equal") return order == 0;
    if (comparison == "not_equal") return order != 0;
    require(!std::holds_alternative<bool>(actual) && !std::holds_alternative<std::string>(actual),
            "native Boolean/text claims require equality or inequality");
    if (comparison == "greater") return order > 0;
    if (comparison == "greater_equal") return order >= 0;
    if (comparison == "less") return order < 0;
    require(comparison == "less_equal", "native claim comparison operator");
    return order <= 0;
}
ResultValue literal_result(const Literal &literal) {
    return std::visit([](const auto &value) { return ResultValue(N(value)); }, literal);
}
NativeRecord native_path(NativeRecord value, const NativeRecord &path) {
    return select_scientific_result(value, path);
}

} // namespace

claim_detail::NativeLiteralObservation claim_detail::read_native_literal(
    const NativeRecord &selected, const std::string &type,
    std::optional<std::size_t> coordinate) {
    auto binding = ResultValue::object();
    binding.set("source_namespace", ResultValue::string(std::to_string(selected.source_namespace())));
    const auto history = selected.publication_history();
    require(history.has_value(), "native literal requires its actual publication");
    binding.set("cutoff", ResultValue::count(history->last_sequence()));
    if (selected.kind() == NativeRecordKind::Object && selected.has("schema") &&
        selected.at("schema").text() == "gadget.scientific.tensor.view.v1") {
        require(type == "rational", "native tensor facts require rational values");
        const auto key = native_count(selected.at("tensor"), "native fact tensor");
        rfc::storage::NativeAccessionRequest request;
        request.subject = key;
        const auto pin = selected.accession_directory()->pin(request);
        require(pin.source_namespace() == selected.source_namespace() &&
                    pin.sequence() == history->last_sequence(), "native fact accession owner mismatch");
        const auto tensor = pin.tensor();
        const auto descriptor = pin.tensor_declaration();
        const auto index = coordinate ? *coordinate :
            selected.has("index") ? native_count(selected.at("index"), "native fact scalar index") : 0;
        if (!coordinate)
            require((selected.at("display").text() == "rational" ||
                         selected.at("display").text() == "rational_object") &&
                        (selected.has("index") || tensor->domain().logical_size() == 1),
                    "nonscalar native facts require an original flat coordinate");
        require(index < tensor->domain().logical_size(), "native fact coordinate outside field");
        const auto value = tensor->state(index);
        require(value.has_value(), "native fact coordinate is missing");
        binding.set("kind", ResultValue::string("native_tensor"));
        binding.set("tensor", ResultValue::count(key));
        binding.set("coordinate", ResultValue::count(index));
        binding.set("supported", ResultValue::boolean(true));
        binding.set("input_carrier", descriptor->input ? ResultValue::count(*descriptor->input) : ResultValue());
        binding.set("output_carrier", descriptor->output ? ResultValue::count(*descriptor->output) : ResultValue());
        binding.set("interpretation", ResultValue::string(descriptor->interpretation));
        return {NativeRecordValue(*value), std::move(binding).release()};
    }
    require(!coordinate, "primitive native fact cannot name a tensor coordinate");
    binding.set("kind", ResultValue::string("native_literal"));
    return {literal_result(native_literal(selected, type, false)).release(), std::move(binding).release()};
}

bool claim_detail::compare_native_literal(const NativeRecord &actual, const Json &expected,
    const std::string &type, const std::string &comparison) {
    auto proposed = admit_record_view(expected);
    if ((type == "count" || type == "integer") && proposed.kind() == NativeRecordKind::Text) {
        const auto text = proposed.text();
        require(!text.empty() && text.find_first_not_of("-0123456789") == std::string::npos,
                "native integer fact proposal requires decimal integer text");
        proposed = NativeRecord::create(N(rfc::numeric::Integer(text)) , new_native_record_namespace());
    }
    return compare_native(native_literal(actual, type, false), native_literal(proposed, type, true), comparison);
}

NativeRecord check_certificate_claims_record(const NativeRecord &certificate, const NativeRecord &request) {
    require(certificate.is_root() && certificate.at("schema").text() == "gadget.operation.native.v1",
            "native claims require the actual owning operation publication");
    require(request.kind() == NativeRecordKind::Object, "native claims request must be an object");
    const auto names = request.names();
    require(std::set<std::string>(names.begin(), names.end()) ==
                std::set<std::string>({"cutoff", "items", "schema", "source_namespace"}),
            "native claims request members");
    require(request.at("schema").text() == "gadget.claims.native.v1", "native claims request schema");
    require(request.at("source_namespace").kind() == NativeRecordKind::Text &&
                request.at("source_namespace").text() == std::to_string(certificate.source_namespace()),
            "native claims publication namespace mismatch");
    const auto history = certificate.publication_history();
    require(history.has_value(), "native claims require a retained publication cutoff");
    const auto cutoff = native_count(request.at("cutoff"), "native claims cutoff");
    require(cutoff == history->last_sequence(), "native claims publication cutoff mismatch");
    const auto result = certificate.at("result");
    const auto schema = result.at("schema").text();
    require(schema == "gadget.scientific.result.fields.v1" || schema == "gadget.biomedical.certificate.v1" ||
                schema == "gadget.claims.native.receipt.v1",
            "native claims refuse an unmigrated numerical certificate");
    const bool scientific = schema == "gadget.scientific.result.fields.v1";
    const auto metadata = scientific ? result.at("certificate") : result;
    const auto values = scientific ? result.at("fields").at("view") :
        schema == "gadget.claims.native.receipt.v1" ? result : result.at("result");
    if (scientific)
        require(result.at("fields").at("schema").text() == "gadget.scientific.fields.v1",
                "native claims require actual scientific fields");
    const auto items = request.at("items");
    require(items.kind() == NativeRecordKind::Array && items.size(), "native claims require items");
    auto checks = ResultValue::array();
    std::size_t matched{};
    for (std::size_t i = 0; i < items.size(); ++i) {
        const auto item = items.at(i);
        require(item.kind() == NativeRecordKind::Object, "native claim must be an object");
        for (const auto &name : item.names())
            require(name == "path" || name == "value_type" || name == "comparison" || name == "value" ||
                        name == "coordinate", "unknown native claim member");
        const auto type = item.at("value_type").text(), comparison = item.at("comparison").text();
        const auto expected = native_literal(item.at("value"), type, true);
        (void)compare_native(expected, expected, comparison);
        const auto selected = native_path(values, item.at("path"));
        auto check = ResultValue::object();
        check.set("claim_index", ResultValue::count(i));
        Literal actual;
        if (selected.kind() == NativeRecordKind::Object && selected.has("schema") &&
            selected.at("schema").kind() == NativeRecordKind::Text &&
            selected.at("schema").text() == "gadget.scientific.tensor.view.v1") {
            require(type == "rational", "native tensor observations require rational claims");
            const auto key = native_count(selected.at("tensor"), "native claim tensor");
            const auto field = selected.encoded_tensor(key);
            require(field.has_value(), "native claim tensor has no original encoded owner");
            const auto &tensor = field->tensor();
            const auto coordinate = item.has("coordinate")
                ? native_count(item.at("coordinate"), "native claim field coordinate")
                : selected.has("index") ? native_count(selected.at("index"), "native scalar coordinate") : 0;
            if (!item.has("coordinate"))
                require((selected.at("display").text() == "rational" || selected.at("display").text() == "rational_object") &&
                            (selected.has("index") || tensor.domain().logical_size() == 1),
                        "nonscalar native claims require an original flat coordinate");
            require(coordinate < tensor.domain().logical_size(), "native claim coordinate outside field");
            require(tensor.support().contains(coordinate), "native claim coordinate is missing");
            actual = *tensor.value(coordinate);
            const auto &pin = field->publication();
            const auto &endpoint = pin.endpoint();
            const auto descriptor = std::find_if(endpoint.tensors.begin(), endpoint.tensors.end(),
                [&](const auto &t) { return t.key == key; });
            require(descriptor != endpoint.tensors.end(), "native claim field descriptor missing");
            auto binding = ResultValue::object();
            binding.set("source_namespace", ResultValue::string(std::to_string(pin.source_namespace())));
            binding.set("cutoff", ResultValue::count(pin.sequence()));
            binding.set("encoding_sequence", ResultValue::count(field->encoding_sequence()));
            binding.set("tensor", ResultValue::count(key));
            binding.set("coordinate", ResultValue::count(coordinate));
            binding.set("supported", ResultValue::boolean(true));
            binding.set("input_carrier", descriptor->input ? ResultValue::count(*descriptor->input) : ResultValue());
            binding.set("output_carrier", descriptor->output ? ResultValue::count(*descriptor->output) : ResultValue());
            binding.set("interpretation", ResultValue::string(descriptor->interpretation));
            check.set("binding", std::move(binding));
        } else {
            require(!item.has("coordinate"), "primitive native claim cannot name a tensor coordinate");
            actual = native_literal(selected, type, false);
        }
        const bool agrees = compare_native(actual, expected, comparison);
        matched += agrees;
        check.set("observed", literal_result(actual));
        check.set("agrees", ResultValue::boolean(agrees));
        checks.push(std::move(check));
    }
    auto receipt = ResultValue::object();
    receipt.set("schema", ResultValue::string("gadget.claims.native.receipt.v1"));
    receipt.set("source_digest", ResultValue::string(metadata.at("source_digest").text()));
    receipt.set("source_namespace", ResultValue::string(std::to_string(certificate.source_namespace())));
    receipt.set("cutoff", ResultValue::count(cutoff));
    receipt.set("checks", std::move(checks));
    receipt.set("matched", ResultValue::count(matched));
    receipt.set("contradicted", ResultValue::count(items.size() - matched));
    receipt.set("scope", ResultValue::string("Literal comparisons against the supplied native operation; interpretation and computation replay are separate"));
    std::uint64_t key = 1;
    for (const auto &bytes : certificate.raw_sources())
        key = std::max(key, rfc::machine::checked_size_add(bytes.key(), 1));
    auto archive = certificate.archive_source(key);
    receipt.set("source_publication", ResultValue(N(N::ByteSpan{key, 0, archive.bytes().size(), false})));
    const auto retained = certificate.derive(std::move(receipt).release(), {std::move(archive)});
    return retain_operation(retained, request, retained.unbound_value());
}
} // namespace gadget
