#include "assay_coverage_detail.hpp"
#include <filesystem>
#include <gadget/record_view.hpp>
#include <iostream>
#include <rfc/field/boundary.hpp>
#include <rfc/storage/file_io.hpp>

using namespace gadget;
using namespace gadget::assay_detail;
namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
std::uint64_t key(NativeRecord selector) { return std::stoull(selector.integer().get_str()); }
const rfc::numeric::SparseRationalTensor &field(const NativeRecord &record,
                                                std::uint64_t identity) {
    for (const auto &tensor : record.endpoint()->tensors)
        if (identity == tensor.key)
            return tensor.tensor;
    throw std::runtime_error("native coverage field not found");
}
} // namespace
int main(int argc, char **argv) {
    const auto root = std::filesystem::path(argc > 1 ? argv[1] : "native-coverage-contracts");
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    auto declaration = Json::parse(
        R"({"reference":{"assembly":"fixture","sequence_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000","coordinate_convention":"zero_based_half_open","contigs":[{"id":"chr","length":12}]},"format":"bed","weight_column":4})");
    ResourceImages images{
        {"control", {declaration, "chr\t0\t8\tbase\t1\t+\nchr\t8\t12\tlast\t2\t+\n", {}}},
        {"treated",
         {declaration,
          "chr\t0\t4\tfirst\t3\t+\nchr\t4\t8\tzero\t0\t+\nchr\t8\t12\tlast\t2\t+\n",
          {}}}};
    const auto request = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"controlled_coverage","samples":[{"id":"base","resource":"control","scale":"1","scale_provenance":"fixture"},{"id":"treated","resource":"treated","scale":"1","scale_provenance":"fixture"}],"contrasts":[{"id":"change","terms":[{"sample":"treated","coefficient":"1"},{"sample":"base","coefficient":"-1"}],"evidence":"fixture"}],"controls":["base"],"intercept":true,"metric":"reference_length","strand_policy":"separate","thresholds":["1"],"absence_policy":"measured_zero_within_union","units":"count","evidence":"fixture","source_kernel":true})");
    const auto record = controlled_coverage_record(images, request);
    auto source_owner = controlled_coverage_source(images, request);
    const auto request_record = admit_record_view(request);
    const auto operation = retain_operation(source_owner.endpoint, new_native_record_namespace(),
                                            std::move(source_owner.sources), request_record,
                                            std::move(source_owner.declaration));
    check(operation.at("result").at("schema").text() == coverage_fields_schema,
          "new source operation lost its typed result");
    check(operation.endpoint()->source.grade({1}).incidence_count() == 6,
          "new source operation substituted an administrative complex");
    check(!operation.parent().has_value(), "unpublished source acquired redundant parent records");
    check(operation.at("request_publication").bytes().size() == request_record.encode().size(),
          "new source operation lost its original request publication");
    check(record.endpoint()->source.grade({0}).cell_count() == 3, "physical union carrier");
    check(record.endpoint()->source.grade({1}).cell_count() == 5, "original interval carrier");
    check(record.endpoint()->source.grade({1}).incidence_count() == 6, "original occurrences");
    check(record.at("schema").text() == coverage_fields_schema, "native result schema");
    check(!record.has("source_records") && !record.has("relation_declaration"),
          "native result admits display rows");
    const auto contrast = record.at("contrasts").at(0);
    check(field(record, key(contrast.at("input"))).values().dense() == Vector({2, -1, 0}),
          "exact coverage contrast");
    check(field(record, key(contrast.at("represented"))).values().dense() ==
              Vector({Rational("1/2"), Rational("1/2"), 0}),
          "control span");
    check(field(record, key(contrast.at("residual"))).values().dense() ==
              Vector({Rational("3/2"), Rational("-3/2"), 0}),
          "control residual");
    check(contrast.at("quadrances").at("residual").rational() == 18, "length weighted quadrance");
    const auto weights_key = record.at("intervals").at("weight");
    const auto &weights = field(record, key(weights_key));
    check(weights.support().count() == 5 && weights.state(3).has_value() &&
              weights.state(3).value() == 0,
          "supported zero interval disappeared");
    const rfc::field::BoundaryOperator boundary(record.endpoint()->source, {1});
    check(boundary.apply(field(record, key(contrast.at("source_coefficients"))).values()).dense() ==
              field(record, key(contrast.at("input"))).values().dense(),
          "source coefficient action and coverage disagree");
    const auto source = record.at("resources").at(0).at("bytes").bytes();
    check(std::string(reinterpret_cast<const char *>(source.data()), source.size()) ==
              images.at("control").bytes,
          "original source bytes changed");
    for (const auto &tensor : record.endpoint()->tensors) {
        check(tensor.output.has_value(), "coverage column lacks an owning carrier");
        const auto grade = tensor.output.value() == 1 ? 0 : 1;
        check(tensor.tensor.domain().logical_size() ==
                  record.endpoint()->source.grade({static_cast<std::uint16_t>(grade)}).cell_count(),
              "field carrier does not address the owning complex");
        check(tensor.tensor.support().count() == tensor.tensor.domain().logical_size(),
              "measured support omitted");
    }
    record.write(root / "coverage.rfcxs");
    const auto restored = NativeRecord::read(root / "coverage.rfcxs");
    check(restored.encode() == record.encode(), "native coverage roundtrip changed bytes");
    check(field(restored, key(restored.at("contrasts").at(0).at("residual"))).values().dense() ==
              Vector({Rational("3/2"), Rational("-3/2"), 0}),
          "restored field differs");
    const auto view = coverage_fields_view(restored);
    check(view.at("source_records").size() == 5 && view.at("segments").size() == 3,
          "explicit coverage view");
    check(view.at("segments").at(0).at("source_record_indices").size() == 2 &&
              view.at("segments").at(2).at("source_record_indices").size() == 2,
          "segment memberships conflated");
    check(view.at("contrasts").at(0).at("source_action_exact").boolean(),
          "display source certificate");
    std::cout << "Native coverage ownership, fields, support, original bytes and replay PASS\n";
}
