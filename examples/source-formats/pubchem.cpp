#include "pubchem.hpp"
#include <charconv>
#include <cstdint>
#include <map>
#include <utility>

namespace gadget::examples {
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::invalid_argument(message);
}
void check_stop(std::stop_token stop) {
    if (stop.stop_requested()) throw std::runtime_error("PubChem source inspection cancelled");
}
void type(const Json &value, json_type expected, const char *message) {
    require(json_object_is_type(value.raw(), expected), message);
}
std::int64_t integer_token(std::string_view token) {
    std::int64_t value{};
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    require(result.ec == std::errc{} && result.ptr == token.data() + token.size(),
            "PubChem integer exceeds signed 64-bit source representation");
    return value;
}
std::string public_id(const std::string &token) {
    const auto value = integer_token(token);
    require(value > 0 && std::to_string(value) == token, "PubChem identity must be a positive canonical integer");
    return token;
}
// Json owns grammar, decoding and the tree. This lexical audit only
// retains numeric spellings that JSON serialization otherwise canonicalizes
// (notably integer -0), and refuses integer saturation in the parser. It also refuses
// lost numeric members, such as overwritten duplicate numeric keys.
class NumberTokens {
    std::map<json_object *, std::string> tokens_;
  public:
    NumberTokens(const std::string &bytes, const Json &root, std::stop_token stop) {
        std::vector<std::string_view> literals;
        for (std::size_t i = 0; i < bytes.size();) {
            check_stop(stop);
            if (bytes[i] == '"') {
                ++i;
                while (i < bytes.size()) {
                    const auto c = bytes[i++];
                    if (c == '\\') ++i;
                    else if (c == '"') break;
                }
            } else if (bytes[i] == '-' || (bytes[i] >= '0' && bytes[i] <= '9')) {
                const auto begin = i++;
                while (i < bytes.size() && (bytes[i] == '+' || bytes[i] == '-' || bytes[i] == '.' ||
                       bytes[i] == 'e' || bytes[i] == 'E' || (bytes[i] >= '0' && bytes[i] <= '9'))) ++i;
                literals.emplace_back(bytes.data() + begin, i - begin);
            } else ++i;
        }
        std::vector<json_object *> pending{root.raw()};
        std::size_t next{};
        while (!pending.empty()) {
            check_stop(stop);
            auto *value = pending.back(); pending.pop_back();
            if (json_object_is_type(value, json_type_int) || json_object_is_type(value, json_type_double)) {
                require(next < literals.size(), "PubChem numeric source/tree mismatch");
                const auto token = literals[next++];
                const Json node(json_object_get(value));
                if (json_object_is_type(value, json_type_int))
                    require(integer_token(token) == node.integer(), "PubChem integer source representation lost");
                else
                    require(token == node.dump_compact(), "PubChem floating source spelling lost");
                tokens_.emplace(value, token);
            } else if (json_object_is_type(value, json_type_array)) {
                for (auto i = json_object_array_length(value); i > 0; --i)
                    pending.push_back(json_object_array_get_idx(value, i - 1));
            } else if (json_object_is_type(value, json_type_object)) {
                std::vector<json_object *> children;
                auto current = json_object_iter_begin(value), end = json_object_iter_end(value);
                while (!json_object_iter_equal(&current, &end)) {
                    children.push_back(json_object_iter_peek_value(&current));
                    json_object_iter_next(&current);
                }
                pending.insert(pending.end(), children.rbegin(), children.rend());
            }
        }
        require(next == literals.size(), "PubChem numeric members lost during JSON decoding");
    }
    const std::string &at(const Json &number) const { return tokens_.at(number.raw()); }
};
std::string id(const Json &value, const NumberTokens &tokens) {
    type(value, json_type_int, "PubChem identity must be an integer");
    return public_id(tokens.at(value));
}
Json pointer(const std::string &path, const Json &original) {
    auto value = Json::object();
    value.set("source_pointer", Json::string(path));
    value.set("original", original);
    return value;
}
void limit_source(const std::string &source, const PubChemSourceLimits &limits, std::stop_token stop) {
    check_stop(stop);
    if (source.size() > limits.source_bytes) throw std::length_error("PubChem source byte limit");
}
} // namespace
PubChemSource inspect_pubchem_assay(std::string source, const std::string &expected_aid,
                                   const PubChemSourceLimits &limits, std::stop_token stop) {
    limit_source(source, limits, stop);
    public_id(expected_aid);
    auto raw = Json::parse(source);
    NumberTokens tokens(source, raw, stop);
    type(raw, json_type_object, "PubChem assay root must be an object");
    const auto submit = raw.at("PC_AssaySubmit");
    type(submit, json_type_object, "PubChem requires full PC_AssaySubmit");
    const auto description = submit.at("assay").at("descr");
    type(description, json_type_object, "PubChem requires full assay description");
    const auto aid = id(description.at("aid").at("id"), tokens);
    require(aid == expected_aid, "PubChem assay differs from declared AID");
    if (description.at("aid").has("version"))
        require(description.at("aid").at("version").integer() >= 0, "negative PubChem assay version");
    type(description.at("name"), json_type_string, "PubChem assay name must be text");
    const auto definitions = description.has("results") ? description.at("results") : Json::array();
    type(definitions, json_type_array, "PubChem result definitions must be an array");
    if (definitions.size() > limits.definitions) throw std::length_error("PubChem definition limit");
    std::map<std::string, std::size_t> tids;
    auto described = Json::array();
    static const char *fields[]{"", "fval", "ival", "bval", "sval"};
    for (std::size_t i = 0; i < definitions.size(); ++i) {
        check_stop(stop);
        const auto definition = definitions.at(i);
        type(definition, json_type_object, "PubChem definition must be an object");
        const auto tid = id(definition.at("tid"), tokens);
        require(tids.emplace(tid, i).second, "duplicate PubChem readout definition");
        type(definition.at("name"), json_type_string, "PubChem readout name must be text");
        const auto code = definition.at("type").integer();
        require(code >= 1 && code <= 4, "unsupported PubChem declared readout type");
        if (definition.has("unit")) (void)definition.at("unit").integer();
        auto entry = pointer("/PC_AssaySubmit/assay/descr/results/" + std::to_string(i), definition);
        entry.set("tid", Json::string(tid));
        entry.set("value_member", Json::string(fields[code]));
        described.push(entry);
    }
    const auto rows = submit.has("data") ? submit.at("data") : Json::array();
    type(rows, json_type_array, "PubChem substance data must be an array");
    if (rows.size() > limits.substances) throw std::length_error("PubChem substance row limit");
    if (rows.size() && definitions.size() > limits.readout_slots / rows.size())
        throw std::length_error("PubChem descriptor/substance coordinate limit");
    auto substances = Json::array();
    for (std::size_t i = 0; i < rows.size(); ++i) {
        check_stop(stop);
        const auto row = rows.at(i);
        type(row, json_type_object, "PubChem substance row must be an object");
        auto occurrence = pointer("/PC_AssaySubmit/data/" + std::to_string(i), row);
        occurrence.set("sid", Json::string(id(row.at("sid"), tokens)));
        if (row.has("version")) require(row.at("version").integer() >= 0, "negative PubChem substance version");
        // Preserve omission separately from the ASN.1 default. Never infer
        // active/inactive status from readout values or discard inactive rows.
        occurrence.set("outcome_reported", Json::boolean(row.has("outcome")));
        if (row.has("outcome")) {
            const auto outcome = row.at("outcome").integer();
            require(outcome >= 1 && outcome <= 5, "unsupported PubChem outcome code");
            occurrence.set("outcome", row.at("outcome"));
        }
        const auto data = row.has("data") ? row.at("data") : Json::array();
        type(data, json_type_array, "PubChem row readouts must be an array");
        std::map<std::string, std::size_t> values;
        for (std::size_t j = 0; j < data.size(); ++j) {
            check_stop(stop);
            const auto value = data.at(j);
            type(value, json_type_object, "PubChem readout must be an object");
            const auto tid = id(value.at("tid"), tokens);
            require(tids.contains(tid), "PubChem readout has no declared definition");
            require(values.emplace(tid, j).second, "duplicate PubChem substance/readout coordinate");
        }
        auto readouts = Json::array();
        for (std::size_t d = 0; d < definitions.size(); ++d) {
            check_stop(stop);
            const auto definition = definitions.at(d);
            const auto tid = id(definition.at("tid"), tokens);
            auto item = Json::object();
            item.set("tid", Json::string(tid));
            item.set("definition_pointer", described.at(d).at("source_pointer"));
            const auto found = values.find(tid);
            item.set("support", Json::string(found == values.end() ? "absent" : "reported"));
            if (found != values.end()) {
                const auto value = data.at(found->second).at("value");
                type(value, json_type_object, "PubChem readout value must be a typed choice");
                const auto code = definition.at("type").integer();
                require(json_object_object_length(value.raw()) == 1 && value.has(fields[code]), "PubChem value differs from declared type");
                const auto scalar = value.at(fields[code]);
                if (code == 1)
                    require(json_object_is_type(scalar.raw(), json_type_int) || json_object_is_type(scalar.raw(), json_type_double),
                            "PubChem fval must be numeric");
                else if (code == 2) type(scalar, json_type_int, "PubChem ival must be an integer");
                else if (code == 3) type(scalar, json_type_boolean, "PubChem bval must be Boolean");
                else type(scalar, json_type_string, "PubChem sval must be text");
                item.set("source_pointer", Json::string("/PC_AssaySubmit/data/" + std::to_string(i) + "/data/" + std::to_string(found->second)));
                item.set("original", data.at(found->second));
                item.set("value_member", Json::string(fields[code]));
                if (code <= 2) item.set("number_token", Json::string(tokens.at(scalar)));
            }
            readouts.push(item);
        }
        occurrence.set("readouts", readouts);
        substances.push(occurrence);
    }
    auto revoked = Json::array();
    if (submit.has("revoke")) {
        const auto source_revoked = submit.at("revoke");
        type(source_revoked, json_type_array, "PubChem revocations must be an array");
        for (std::size_t i = 0; i < source_revoked.size(); ++i) {
            check_stop(stop);
            auto item = pointer("/PC_AssaySubmit/revoke/" + std::to_string(i), source_revoked.at(i));
            item.set("sid", Json::string(id(source_revoked.at(i), tokens)));
            revoked.push(item);
        }
    }
    auto projection = Json::object();
    projection.set("aid", Json::string(aid));
    projection.set("description", pointer("/PC_AssaySubmit/assay/descr", description));
    projection.set("definitions", described);
    projection.set("substances", substances);
    projection.set("revocations", revoked);
    check_stop(stop);
    return {std::move(source), raw, projection};
}
PubChemSource inspect_pubchem_concise(std::string source, const std::string &expected_aid,
                                     const PubChemSourceLimits &limits, std::stop_token stop) {
    limit_source(source, limits, stop);
    public_id(expected_aid);
    const auto table = read_csv(source);
    for (const auto *name : {"AID", "SID", "CID"}) {
        const auto column = table.columns.find(name);
        require(column != table.columns.end() && table.header.at(column->second) == name,
                "PubChem concise identity columns missing or renamed");
    }
    if (table.rows.size() > limits.concise_rows) throw std::length_error("PubChem concise row limit");
    auto original = Json::object(), rows = Json::array(), links = Json::array(), header = Json::array();
    for (const auto &name : table.header) header.push(Json::string(name));
    original.set("header", header);
    for (const auto &row : table.rows) {
        check_stop(stop);
        const auto aid = public_id(table.get(row, "AID"));
        require(aid == expected_aid, "PubChem concise row differs from declared AID");
        auto cells = Json::array();
        for (const auto &cell : row.cells) cells.push(Json::string(cell));
        auto occurrence = Json::object();
        occurrence.set("source_record", Json::integer(static_cast<std::int64_t>(row.record)));
        occurrence.set("cells", cells);
        rows.push(occurrence);
        const auto add = [&](const std::string &left, const std::string &right, const char *meaning) {
            if (links.size() >= limits.links) throw std::length_error("PubChem explicit link limit");
            auto link = Json::object();
            link.set("source_record", occurrence.at("source_record"));
            link.set("cells", cells);
            link.set("left", Json::string(left)); link.set("right", Json::string(right));
            link.set("relationship", Json::string(meaning)); links.push(link);
        };
        const auto &sid = table.get(row, "SID"), &cid = table.get(row, "CID");
        if (!sid.empty()) public_id(sid);
        if (!cid.empty()) public_id(cid);
        if (!sid.empty() && !cid.empty())
            add("pubchem:sid:" + sid, "pubchem:cid:" + cid, "reported_substance_compound_link");
        if (table.columns.contains("Target GeneID") && !table.get(row, "Target GeneID").empty())
            add("pubchem:aid:" + aid, "ncbi:gene:" + public_id(table.get(row, "Target GeneID")), "reported_assay_target");
        if (table.columns.contains("Target Accession") && !table.get(row, "Target Accession").empty())
            add("pubchem:aid:" + aid, "pubchem:target_accession:" + table.get(row, "Target Accession"), "reported_assay_target_accession");
    }
    original.set("rows", rows);
    auto projection = Json::object();
    projection.set("aid", Json::string(expected_aid)); projection.set("links", links);
    check_stop(stop);
    return {std::move(source), original, projection};
}
} // namespace gadget::examples
