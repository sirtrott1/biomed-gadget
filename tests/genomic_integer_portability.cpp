#include "genomic_detail.hpp"
#include <gadget/native_record.hpp>
#include <iostream>

using namespace gadget;
using namespace gadget::assay_detail;
namespace {
void check(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
NativeRecord publish(AssayFieldsSource source) {
    return NativeRecord::create_from_endpoint(std::move(source.declaration),
        new_native_record_namespace(), std::move(source.sources), std::move(source.endpoint));
}
void column(const NativeRecord &record, const char *name, const char *expected) {
    const auto key = record.at(name).integer();
    for (const auto &field : record.endpoint()->tensors)
        if (key == rfc::numeric::Integer(std::to_string(field.key))) {
            check(field.tensor.state(0) == rfc::numeric::Rational(expected),
                  "genomic coordinate lost exact integer precision");
            return;
        }
    throw std::runtime_error("coordinate field is absent");
}
}
int main() try {
    // Odd coordinates beyond 2^53 expose a conversion through double.
    const auto reference = Json::parse(R"({"assembly":"synthetic-large-coordinate-frame",
        "sequence_sha256":"sha256:0000000000000000000000000000000000000000000000000000000000000000",
        "coordinate_convention":"zero_based_half_open",
        "contigs":[{"id":"chrLarge","length":9007199254741001}]})");
    ResourceImage resource;
    resource.declaration = Json::parse(R"({"id":"large.bed","format":"bed"})");
    resource.declaration.set("reference", reference);
    resource.bytes = "chrLarge\t9007199254740993\t9007199254741001\tfeature\t1\t+\n";
    const auto request = Json::parse(R"({"schema":"gadget.assays.query.v1",
        "operation":"interval_coverage","resource":"large.bed",
        "thresholds":["1"],"strand_policy":"separate"})");
    const auto coverage = NativeRecord::decode(publish(interval_fields_source(resource, request)).encode());
    for (const auto *name : {"record_start", "start"}) column(coverage, name, "9007199254740993");
    for (const auto *name : {"record_end", "end"}) column(coverage, name, "9007199254741001");
    column(coverage, "lengths", "8");
    const auto region = coverage.at("threshold_components").at(0).at("regions").at(0);
    check(region.at("start").integer() == rfc::numeric::Integer("9007199254740993") &&
              region.at("end").integer() == rfc::numeric::Integer("9007199254741001"),
          "large threshold region endpoints changed");
    const auto annotation_request = Json::parse(R"({"schema":"gadget.assays.query.v1",
        "operation":"annotation_features","resource":"large.bed",
        "feature_type":"interval","group_attribute":"id"})");
    const auto annotation = NativeRecord::decode(publish(annotation_fields_source(resource, annotation_request)).encode());
    column(annotation, "record_start", "9007199254740993");
    column(annotation, "record_end", "9007199254741001");
    const auto interval = annotation.at("features").at(0).at("union_intervals").at(0);
    check(interval.at(0).integer() == rfc::numeric::Integer("9007199254740993") &&
              interval.at(1).integer() == rfc::numeric::Integer("9007199254741001"),
          "large annotation union endpoints changed");
#if GADGET_HAS_HTSLIB
    resource.declaration.set("id", Json::string("large.vcf"));
    resource.declaration.set("format", Json::string("vcf"));
    resource.bytes = "##fileformat=VCFv4.3\n"
        "##contig=<ID=chrLarge,length=9007199254741001>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "chrLarge\t9007199254740994\tvariant\tA\tC\t.\tPASS\t.\tGT\t0/1\n";
    const auto genotype_request = Json::parse(R"({"schema":"gadget.assays.query.v1",
        "operation":"variant_genotypes","resource":"large.vcf","filter_policy":"all"})");
    std::optional<AssayFieldsSource> retained;
    (void)genotype_analysis(resource, genotype_request, &retained);
    check(retained.has_value(), "genotype fields were not retained");
    const auto genotype = NativeRecord::decode(publish(std::move(*retained)).encode());
    column(genotype, "record_start", "9007199254740993");
    column(genotype, "record_end", "9007199254740994");
#endif
    std::cout << "exact genomic coordinates above 2^53 survive native publication and reopening\n";
} catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
