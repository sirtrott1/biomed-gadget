#include "pipeline_result_binding.hpp"
#include <algorithm>
#include <set>

namespace gadget::pipeline_detail {
namespace {
void require(bool ok, const std::string &message) {
    if (!ok)
        throw std::invalid_argument(message);
}
void keys(const Json &value, const std::set<std::string> &expected) {
    require(json_object_is_type(value.raw(), json_type_object),
            "result operation requires an object");
    for (const auto &key : expected)
        require(value.has(key), "missing result operation member: " + key);
    json_object_object_foreach(value.raw(), key, member) {
        (void)member;
        require(expected.contains(key), "unknown result operation member: " + std::string(key));
    }
}
std::size_t bound(const Json &value) {
    const auto n = value.integer();
    require(n > 0 && static_cast<std::uint64_t>(n) <= source_byte_limit(),
            "result archive byte bound");
    return static_cast<std::size_t>(n);
}
void digest(const std::string &value) {
    require(
        value.size() == 71 && value.starts_with("sha256:") &&
            std::all_of(value.begin() + 7, value.end(),
                        [](char c) { return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'); }),
        "result archive SHA-256 must be canonical");
}
} // namespace
ResultOperation admit_result_operation(const Json &request) {
    keys(request,
         {"schema", "operation", "source_record", "analysis_bytes", "analysis_binding", "query"});
    require(request.at("schema").text() == "gadget.result-operation.v1",
            "result operation schema");
    const auto method = request.at("operation").text();
    require(method == "pathway_overrepresentation" || method == "feature_evidence_context",
            "unsupported result operation operation");
    const auto source = request.at("source_record");
    keys(source, {"path", "sha256", "max_bytes"});
    const auto path = source.at("path").text();
    require(!path.empty() && path.find('\0') == std::string::npos, "result source archive path");
    const auto hash = source.at("sha256").text();
    digest(hash);
    require(json_object_is_type(request.at("query").raw(), json_type_object),
            "result producer query must be an object");
    const auto mode = request.at("analysis_binding").text();
    require(mode == "exact" || mode == "actual_result", "analysis binding mode");
    require(mode == "exact" ? request.at("query").has("analysis_sha256")
                            : !request.at("query").has("analysis_sha256"),
            "analysis pin conflicts with declared binding mode");
    return {method,
            mode,
            fs::path(path),
            hash,
            bound(source.at("max_bytes")),
            bound(request.at("analysis_bytes"))};
}
void admit_result_step(const Json &node, const std::map<std::string, std::string> &declared) {
    require(node.at("type").text() == "result", "result input belongs only to a result step");
    require(node.has("result_input"), "result step requires a result_input");
    const auto input = node.at("result_input").text();
    require(declared.contains(input) && declared.at(input) == "certificate",
            "result_input must name an earlier result-producing step");
    for (const auto *member : {"assay", "partition", "normalization_evidence_field"})
        require(!node.has(member), "result step cannot change an assay or normalization binding");
    (void)admit_result_operation(node.at("request"));
}
std::string verified_record_image(const fs::path &path, const std::string &hash,
                                  std::size_t limit) {
    digest(hash);
    require(limit > 0 && limit <= source_byte_limit(), "result archive byte bound");
    const auto bytes = read_text(path, limit);
    require("sha256:" + sha256(bytes) == hash, "result archive digest mismatch");
    return bytes;
}
} // namespace gadget::pipeline_detail
