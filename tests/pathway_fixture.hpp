#pragma once
#include <gadget/pathways.hpp>
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>

namespace gadget::pathway_test {
struct Fixture {
    NativeRecord analysis, gene_sets, request;
};
inline std::string archive_sha256(const NativeRecord &record) {
    const auto bytes = record.owning_publication().encode();
    return "sha256:" +
           sha256(std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
}
inline Fixture fixture() {
    using namespace rfc::storage;
    rfc::core::RelationSignature signature;
    rfc::core::RelationalComplexBuilder builder(signature);
    for (std::uint64_t i = 1; i <= 4; ++i)
        builder.declare_witness_identity(i);
    TemporalEndpoint endpoint(builder.freeze());
    endpoint.carriers.push_back(
        {1, TemporalCarrierKind::SourceGrade, {0}, {}, 0, rfc::numeric::IndexDomain::linear(4)});
    TemporalTensor tensor;
    tensor.key = 101;
    tensor.role = TemporalTensorRole::Field;
    tensor.output = 1;
    tensor.interpretation = "retained molecular counts";
    tensor.tensor = {rfc::numeric::IndexDomain::linear(4), rfc::numeric::StructuralSupport::full(4),
                     rfc::numeric::SparseRationalVector(
                         4, {{0, Rational(7)}, {2, Rational(9)}, {3, Rational(2)}})};
    endpoint.tensors.push_back(std::move(tensor));
    const auto identities = Json::parse(
        R"({"schema":"gadget.pathways.features.v1","features":[{"id":"ENSG1.1","source_cell":"1"},{"id":"ENSG1.1_PAR_Y","source_cell":"2"},{"id":"ENSG2.1","source_cell":"3"},{"id":"ENSG3.1","source_cell":"4"}],"selection_scope":"Explicit source-bound feature declarations; no inferred molecular significance"})");
    auto analysis = NativeRecord::create(import_record_view(identities), 26001, {}, &endpoint);
    const std::string gmt = "set-a\tauthored release\tG1\tG2\tG2\t\nset-b\tother\tG3\tG4\n";
    const std::string mapping =
        "feature_id\tmember\tevidence\nENSG1.1\tG1\tdeclared exact "
        "correspondence\nENSG1.1_PAR_Y\tG2\tdistinct Y feature\nENSG2.1\tG3\tdeclared exact "
        "correspondence\nENSG3.1\tG4\tdeclared exact correspondence\n";
    auto source = Json::parse(
        R"({"schema":"gadget.pathways.source.v1","id":"authored-sets","release":"authored.v1","identifier_namespace":"declared:gene","format":"gmt","mapping":{"id":"authored-correspondences","release":"authored.v1","origin":"explicit_declared","feature_namespace":"ensembl:gene","identifier_namespace":"declared:gene","format":"tsv"}})");
    source.set("data", Json::string(gmt));
    source.set("sha256", Json::string("sha256:" + sha256(gmt)));
    source.at("mapping").set("data", Json::string(mapping));
    source.at("mapping").set("sha256", Json::string("sha256:" + sha256(mapping)));
    auto gene_sets = NativeRecord::create(import_record_view(source), 26002);
    auto query = Json::parse(
        R"({"schema":"gadget.pathways.overrepresentation.v1","operation":"pathway_overrepresentation","feature_namespace":"ensembl:gene","features_path":["features"],"feature_id_path":["id"],"background_indices":[0,1,2,3],"selected_indices":[0,1],"mapping_policy":{"ambiguous":"refuse","unmapped":"refuse","many_to_one":"refuse","duplicate_occurrences":"retain_count_once"},"set_ids":["set-a","set-b"],"size_filter":{"minimum":1,"maximum":100,"basis":"members_in_background"},"membership_policy":{"duplicates":"retain_count_once","empty":"retain_exclude","outside_background":"retain_exclude"},"correction":"benjamini_yekutieli","selection_method":"Explicit retained native feature occurrences","null_assumptions":["Uniform conditional subset law for this authored fixture"],"limits":{"source_bytes":1048576,"mapping_rows":10000,"membership_cells":100000,"sets":1000,"features":1000,"field_bytes":4096,"kernel_work":1000000,"archive_bytes":1048576,"declared_buffer_bytes":16777216}})");
    query.set("analysis_sha256", Json::string(archive_sha256(analysis)));
    query.set("gene_sets_sha256", Json::string(archive_sha256(gene_sets)));
    auto request = NativeRecord::create(import_record_view(query), 26003);
    return {std::move(analysis), std::move(gene_sets), std::move(request)};
}
} // namespace gadget::pathway_test
