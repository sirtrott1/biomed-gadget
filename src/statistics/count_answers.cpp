#include "count_answers.hpp"
#include "../pipeline_answers.hpp"
#include <gadget/record_view.hpp>
#include <gadget/statistics/count_analysis.hpp>

namespace gadget::statistics {
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::invalid_argument(message);
}
std::string digest(std::span<const std::byte> bytes) {
    return "sha256:" +
           sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
}
NativeRecord checked(const fs::path &directory, const Json &reference, std::size_t limit) {
    const auto path = pipeline_detail::checked_answer_path(directory, reference);
    require(fs::file_size(path) <= limit, "count answer archive byte admission");
    require("sha256:" + sha256_file(path) == reference.at("sha256").text(),
            "count answer archive hash differs");
    return NativeRecord::read(path);
}
} // namespace
Json count_answer_view(const fs::path &directory, const Json &reference) {
    const auto selection = reference.at("native_selection");
    require(selection.at("kind").text() == "count_model" &&
                json_object_object_length(selection.raw()) == 2,
            "count answer selection declaration");
    const auto request = checked(directory, selection.at("request"), source_byte_limit());
    require(request.at("schema").text() == "gadget.statistics.count.request.v1" &&
                request.at("operation").text() == "count_model_analysis",
            "count answer request schema/operation");
    const auto raw_limit = request.at("limits").at("declared_buffer_bytes").integer();
    require(raw_limit > 0 &&
                raw_limit <= rfc::numeric::Integer(std::to_string(source_byte_limit())),
            "count answer numerical owner byte admission");
    const auto publication = checked(directory, reference, raw_limit.get_ui());
    require(publication.at("schema").text() == "gadget.operation.native.v1",
            "count answer requires its complete operation");
    require(digest(publication.at("request_publication").bytes()) ==
                selection.at("request").at("sha256").text(),
            "count answer request owner differs");
    require(render_record_view(publication.at("request")).dump() ==
                render_record_view(request).dump(),
            "count answer request value differs");
    const auto fields = publication.at("result"), view = fields.at("view");
    require(view.at("schema").text() == "gadget.statistics.count.result.v1" &&
                view.at("source_digest").text() == request.at("source_digest").text(),
            "count answer source binding differs");
    const ScientificStudy source(NativeRecord::decode(view.at("source_study_publication").bytes()));
    require(source.digest() == view.at("source_digest").text(),
            "count answer original study differs");
    require(reference.at("source_digest").text() == source.digest(),
            "count answer source reference differs");
    auto out = Json::object();
    out.set("schema", Json::string("gadget.statistics.count.answer.v1"));
    out.set("source_digest", Json::string(source.digest()));
    out.set("query", render_record_view(request));
    out.set("result", count_result_view(publication));
    out.set("certificate_sha256", Json::string("sha256:" + sha256(out.dump())));
    if (reference.has("query"))
        require(reference.at("query").dump() == out.at("query").dump(),
                "count answer displayed request differs");
    if (reference.has("certificate_sha256"))
        require(reference.at("certificate_sha256").text() == out.at("certificate_sha256").text(),
                "count answer certificate differs");
    return out;
}
Json count_answer_reference(const fs::path &directory, const Json &result, const Json &request) {
    auto reference = Json::parse(result.dump()), selection = Json::object();
    selection.set("kind", Json::string("count_model"));
    selection.set("request", request);
    reference.set("native_selection", selection);
    reference.set("pointer", Json::string("/result"));
    const auto retained = checked(directory, request, source_byte_limit());
    reference.set("source_digest", Json::string(retained.at("source_digest").text()));
    const auto certificate = count_answer_view(directory, reference);
    for (const auto *name : {"query", "certificate_sha256"})
        reference.set(name, certificate.at(name));
    return reference;
}
} // namespace gadget::statistics
