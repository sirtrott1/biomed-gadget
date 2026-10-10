#include "genomic_detail.hpp"
#include <gadget/record_view.hpp>
#include <iostream>
#include <rfc/field/boundary.hpp>

using namespace gadget;
using namespace gadget::assay_detail;
namespace {
void check(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
const rfc::numeric::SparseRationalTensor &field(const NativeRecord &record, std::string_view group,
                                                std::string_view name) {
    const auto key = record.at(group).at(name).integer();
    for (const auto &tensor : record.endpoint()->tensors)
        if (key == rfc::numeric::Integer(std::to_string(tensor.key)))
            return tensor.tensor;
    throw std::runtime_error("field owner missing");
}
} // namespace
int main() {
    try {
        const std::string sam =
            "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chrA\tLN:100\n@RG\tID:lane\tSM:sample\n"
            "a\t0\tchrA\t1\t60\t2M5N3M\t*\t0\t0\tAAAAA\tIIIII\tRG:Z:lane\n"
            "b\t1024\tchrA\t1\t60\t2M5N3M\t*\t0\t0\tAAAAA\tIIIII\n"
            "c\t16\tchrA\t20\t255\t2M\t*\t0\t0\tAA\tII\n"
            "u\t4\t*\t0\t0\t*\t*\t0\t0\tAA\tII\n"
            "pair\t65\tchrA\t40\t0\t2M\t*\t0\t0\tAA\tII\n";
        auto declaration = Json::parse(
            R"({"format":"sam","reference":{"assembly":"fixture","sequence_sha256":"sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","coordinate_convention":"zero_based_half_open","contigs":[{"id":"chrA","length":100}]}})");
        ResourceImage resource{declaration, sam, {}};
        auto request = Json::parse(
            R"({"schema":"gadget.assays.query.v1","operation":"alignment_evidence","resource":"fixture","filter":{"exclude_flags":1024,"minimum_mapping_quality":0,"unknown_mapping_quality":"include"}})");
        auto source = alignment_fields_source(resource, request);
        const auto record = NativeRecord::create_from_endpoint(
            source.declaration, new_native_record_namespace(), std::move(source.sources),
            std::move(source.endpoint));
        const auto restarted = NativeRecord::decode(record.encode());
        check(restarted.encode() == record.encode(), "complete native replay changed");
        check(restarted.endpoint()->source.grade_count() == 2, "unasserted higher topology");
        check(restarted.endpoint()->source.grade({0}).cell_count() == 5 &&
                  restarted.endpoint()->source.grade({1}).cell_count() == 7,
              "record, block or junction identity lost");
        const auto &starts = field(restarted, "record_fields", "start");
        check(starts.state(0) == Rational{0} && !starts.state(3),
              "mapped zero and unmapped absence confused");
        const auto &selected = field(restarted, "record_fields", "selected");
        check(selected.state(1) == Rational{0} && selected.state(2) == Rational{1},
              "duplicate filter or unknown MAPQ policy");
        const auto &counts = field(restarted, "relation_fields", "selected_record_count");
        check(counts.state(6) == Rational{1}, "junction count changed");
        const auto &grade = restarted.endpoint()->source.grade({1});
        const auto [begin, end] = grade.incidence_range(6);
        check(end - begin == 2 && grade.lower_index[begin] == 0 &&
                  grade.lower_index[begin + 1] == 1,
              "excluded junction occurrence erased");
        const auto raw = restarted.at("original").bytes();
        check(std::string_view(reinterpret_cast<const char *>(raw.data()), raw.size()) == sam,
              "original SAM bytes changed");
        check(render_record_view(restarted).dump() ==
                  genomic_analysis({{"fixture", resource}}, request, CpuResources::local()).dump(),
              "native alignment view changed the complete result");
        const auto before = decode_alignments(resource),
                   compact = decode_alignments(resource, false);
        check(before.records.size() == compact.records.size() &&
                  compact.records[0].decoded.empty() &&
                  before.records[0].blocks == compact.records[0].blocks &&
                  before.records[0].junctions == compact.records[0].junctions,
              "compact decoder changed physical evidence");
        const auto complete = inspect_hts(resource), summary = inspect_hts(resource, true);
        check(complete.at("record_count").integer() == summary.at("record_count").integer() &&
                  !summary.has("records") &&
                  summary.at("source_sha256").text() == complete.at("source_sha256").text(),
              "compact inspection changed validation or source identity");
        request.at("filter").set("unknown_mapping_quality", Json::string("exclude"));
        auto other = alignment_fields_source(resource, request);
        const auto filtered =
            NativeRecord::create_from_endpoint(other.declaration, new_native_record_namespace(),
                                               std::move(other.sources), std::move(other.endpoint));
        check(field(filtered, "record_fields", "selected").state(2) == Rational{0},
              "unknown MAPQ exclusion ignored");
        auto invalid = resource;
        invalid.bytes = sam + "truncated\n";
        bool rejected = false;
        try {
            (void)inspect_hts(invalid, true);
        } catch (const std::exception &) {
            rejected = true;
        }
        check(rejected, "summary inspection accepted a truncated record");
        const std::string gtf =
            "chrA\tfixture\texon\t1\t2\t.\t+\t.\tgene_id \"A\"; transcript_id \"A1\";\n"
            "chrA\tfixture\texon\t8\t10\t.\t+\t.\tgene_id \"A\"; transcript_id \"A1\";\n"
            "chrA\tfixture\texon\t1\t4\t.\t+\t.\tgene_id \"B\"; transcript_id \"B1\";\n";
        auto annotation = Json::parse(declaration.dump());
        annotation.set("format", Json::string("gtf"));
        ResourceImages images{{"fixture", resource}, {"annotation", {annotation, gtf, {}}}};
        auto feature_request = Json::parse(
            R"({"schema":"gadget.assays.query.v1","operation":"feature_counts","resource":"fixture","annotation":"annotation","feature_type":"exon","group_attribute":"gene_id","unit":"sequenced_segment","overlap":"any_aligned_base","strand_policy":"ignore","filter":{"exclude_flags":1024,"minimum_mapping_quality":0,"unknown_mapping_quality":"include"},"candidate_search":{"method":"Authored mappings","completeness":"Complete fixture","reference_sha256":"sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"}})");
        const auto local = CpuResources::local();
        rfc::runtime::CpuResourcePolicy policy;
        policy.worker_limit = 1;
        rfc::runtime::CpuResourceCoordinator coordinator(local->readmit(), policy);
        auto one_worker = local->assigned(coordinator);
        auto feature_source = feature_count_fields_source(images, feature_request, one_worker);
        const auto feature = NativeRecord::create_from_endpoint(
            feature_source.declaration, new_native_record_namespace(),
            std::move(feature_source.sources), std::move(feature_source.endpoint));
        const auto resumed = NativeRecord::decode(feature.encode());
        check(resumed.endpoint()->source.grade({0}).cell_count() == 3 &&
                  resumed.endpoint()->source.grade({1}).cell_count() == 2,
              "allocation class complex not retained");
        const auto &lo = field(resumed, "fields", "lower"), &hi = field(resumed, "fields", "upper");
        check(lo.state(0) == Rational{0} && hi.state(0) == Rational{1} &&
                  lo.state(1) == Rational{0} && hi.state(1) == Rational{1} && !lo.state(2),
              "assignment ambiguity, measured zero or unrequested support changed");
        check(resumed.at("segments").at("class_index").tensor().domain().logical_size() == 4 &&
                  resumed.at("segments").at("source_records").tensor().domain().logical_size() == 4,
              "original sequenced segments lost");
        const auto legacy_view = genomic_analysis(images, feature_request, one_worker);
        check(render_record_view(resumed).dump() == legacy_view.dump(),
              "native assignment view changed the complete result");
        const auto legacy = legacy_view.at("allocation_certificate").at("result");
        check(resumed.at("allocation_summary").at("target_ambiguity_rank").integer() ==
                  legacy.at("target_ambiguity_rank").integer(),
              "allocation rank changed");
        std::cout << "Alignment fields and complete source checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
