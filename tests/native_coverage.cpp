#include "assay_coverage_detail.hpp"
#include <filesystem>
#include <gadget/record_view.hpp>
#include <iostream>
#include <limits>
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
    const auto interval_request = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"interval_coverage","resource":"treated","thresholds":["0","1","4"],"strand_policy":"separate"})");
    auto interval_source = interval_fields_source(images.at("treated"), interval_request);
    const auto intervals = NativeRecord::create_from_endpoint(
        std::move(interval_source.declaration), 740001, std::move(interval_source.sources),
        std::move(interval_source.endpoint));
    check(intervals.endpoint()->source.grade({0}).cell_count() == 3 &&
              intervals.endpoint()->source.grade({1}).incidence_count() == 3,
          "interval publication replaced original membership");
    check(field(intervals, key(intervals.at("coverage"))).values().dense() == Vector({3, 0, 2}) &&
              field(intervals, key(intervals.at("coverage"))).support().count() == 3,
          "interval supported zero changed");
    check(intervals.at("coverage_integral").rational() == 20 &&
              intervals.at("threshold_components").at(0).at("retained_length").rational() == 12 &&
              intervals.at("threshold_components").at(1).at("regions").size() == 2 &&
              intervals.at("threshold_components").at(2).at("regions").size() == 0,
          "interval exact integral or physical threshold components changed");
    const auto interval_view = render_record_view(intervals);
    check(interval_view.at("segments").at(1).at("coverage").text() == "0" &&
              interval_view.at("threshold_components").at(0).at("retained_atoms").integer() == 3,
          "interval native view omitted support");
    intervals.write(root / "interval.rfcxs");
    const auto interval_restored = NativeRecord::read(root / "interval.rfcxs");
    check(interval_restored.encode() == intervals.encode() &&
              render_record_view(interval_restored).dump() == interval_view.dump(),
          "interval source/value roundtrip differs");
    const auto original = interval_restored.at("source_bytes").bytes();
    check(std::string(reinterpret_cast<const char *>(original.data()), original.size()) ==
              images.at("treated").bytes &&
              interval_restored.at("original_records").at(1).at("bytes").bytes().size() ==
                  std::string("chr\t4\t8\tzero\t0\t+").size(),
          "interval original bytes or span changed");
    auto gff_declaration = Json::parse(declaration.dump());
    gff_declaration.set("format", Json::string("gff3"));
    const ResourceImage annotation{gff_declaration,
        "##gff-version 3\nchr\tfixture\texon\t1\t4\t.\t+\t.\tID=e1;Parent=F1,F2\n"
        "chr\tfixture\texon\t3\t6\t.\t+\t.\tID=e2;Parent=F1\n"
        "chr\tfixture\tgene\t1\t12\t.\t+\t.\tID=F1\n", {}};
    const auto annotation_request = Json::parse(
        R"({"schema":"gadget.assays.query.v1","operation":"annotation_features","resource":"annotation","feature_type":"exon","group_attribute":"Parent"})");
    auto annotation_source = annotation_fields_source(annotation, annotation_request);
    const auto features = NativeRecord::create_from_endpoint(
        std::move(annotation_source.declaration), 740002, std::move(annotation_source.sources),
        std::move(annotation_source.endpoint));
    check(features.endpoint()->source.grade({0}).cell_count() == 3 &&
              features.endpoint()->source.grade({1}).cell_count() == 2 &&
              features.endpoint()->source.grade({1}).incidence_count() == 3,
          "annotation filtered/repeated group source occurrences changed");
    const auto feature_view = render_record_view(features).at("features");
    check(feature_view.at(0).at("source_record_indices").size() == 2 &&
              feature_view.at(1).at("source_record_indices").at(0).integer() == 0 &&
              feature_view.at(0).at("union_intervals").at(0).at(1).integer() == 6 &&
              feature_view.at(1).at("union_intervals").at(0).at(1).integer() == 4,
          "annotation physical union or original group order changed");
    check(field(features, key(features.at("record_start"))).values().dense() == Vector({0, 2, 0}) &&
              field(features, key(features.at("record_start"))).support().count() == 3,
          "annotation original coordinates or supported zero changed");
    features.write(root / "annotation.rfcxs");
    check(NativeRecord::read(root / "annotation.rfcxs").encode() == features.encode(),
          "annotation original owner roundtrip differs");
    const ResourceImage wide_intervals{
        Json::parse(R"({"reference":{"assembly":"fixture","sequence_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000","coordinate_convention":"zero_based_half_open","contigs":[{"id":"chr","length":9223372036854775807}]},"format":"bed","weight_column":4})"),
        "chr\t9007199254740993\t9007199254740997\tfirst\t1\t+\n"
        "chr\t9223372036854775803\t9223372036854775807\tlast\t0\t+\n", {}};
    auto wide_source = interval_fields_source(wide_intervals, interval_request);
    const auto wide = NativeRecord::create_from_endpoint(
        std::move(wide_source.declaration), std::numeric_limits<std::uint64_t>::max() - 1,
        std::move(wide_source.sources), std::move(wide_source.endpoint));
    check(field(wide, key(wide.at("record_start"))).values().dense() ==
              Vector({Rational("9007199254740993"), Rational("9223372036854775803")}) &&
              field(wide, key(wide.at("record_end"))).values().dense() ==
              Vector({Rational("9007199254740997"), Rational("9223372036854775807")}),
          "wide original coordinates lost integer precision");
    check(wide.at("coverage_integral").rational() == 4,
          "wide interval coverage integral changed");
    wide.write(root / "wide.rfcxs");
    const auto wide_restored = NativeRecord::read(root / "wide.rfcxs");
    check(wide_restored.source_namespace() == std::numeric_limits<std::uint64_t>::max() - 1 &&
              wide_restored.encode() == wide.encode(),
          "wide native namespace or coordinate roundtrip changed");
    std::cout << "Native coverage ownership, fields, support, original bytes and replay PASS\n";
}
