#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
#include <algorithm>
#include <rfc/core/builder.hpp>
#include <rfc/field/boundary.hpp>

namespace gadget::assay_detail {
namespace {
using N = NativeRecordValue;
N integer(std::size_t value) { return N(rfc::numeric::Integer(std::to_string(value))); }
struct Dictionary {
    std::map<std::string, std::size_t> indices;
    N::Array labels;
    std::size_t add(const std::string &label) {
        const auto [it, added] = indices.emplace(label, indices.size());
        if (added)
            labels.elements.emplace_back(label);
        return it->second;
    }
};
struct Relation {
    std::size_t contig, strand, kind;
    std::uint64_t start, end;
};
} // namespace

AssayFieldsSource alignment_fields_source(const ResourceImage &resource, const Json &request) {
    members(request, {"schema", "operation", "resource", "filter"});
    require(request.at("schema").text() == "gadget.assays.query.v1" &&
                request.at("operation").text() == "alignment_evidence",
            "alignment query schema");
    const auto policy = request.at("filter");
    members(policy, {"exclude_flags", "minimum_mapping_quality", "unknown_mapping_quality"});
    const auto excluded = count(policy.at("exclude_flags"), 4095, "alignment excluded flags");
    const auto minimum = count(policy.at("minimum_mapping_quality"), 254, "alignment minimum MAPQ");
    const auto unknown = policy.at("unknown_mapping_quality").text();
    require(unknown == "include" || unknown == "exclude", "MAPQ 255 handling");
    auto decoded = decode_alignments(resource, false);
    const auto size = decoded.records.size();
    Dictionary contigs, strands, groups;
    std::vector<std::size_t> record_contig, record_strand, record_group, name_start, name_length,
        mapped;
    Vector selected;
    std::string names;
    using Key = std::tuple<std::string, std::string, std::uint64_t, std::uint64_t>;
    std::map<Key, std::vector<std::size_t>> junctions;
    rfc::core::RelationSignature signature;
    const auto block_type = signature.add_type("alignment.aligned.block");
    const auto junction_type = signature.add_type("alignment.splice.junction.support");
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade({0});
    builder.declare_grade({1});
    for (std::size_t i = 0; i < size; ++i)
        builder.declare_witness_identity(i + 1);
    std::vector<Relation> relations;
    std::size_t occurrence = 0;
    for (std::size_t i = 0; i < size; ++i) {
        const auto &record = decoded.records[i];
        const auto contig = contigs.add(record.contig), strand = strands.add(record.strand);
        record_contig.push_back(contig);
        record_strand.push_back(strand);
        record_group.push_back(groups.add(record.group));
        name_start.push_back(names.size());
        name_length.push_back(record.query.size());
        require(record.query.size() <= names.max_size() - names.size(),
                "alignment query bytes overflow");
        names.append(record.query);
        selected.emplace_back(!(record.flags & excluded) &&
                              (record.mapq == 255 ? unknown == "include" : record.mapq >= minimum));
        if (!(record.flags & 4))
            mapped.push_back(i);
        for (auto [start, end] : record.blocks) {
            occurrence = rfc::machine::checked_size_add(occurrence, 1);
            builder.add_grade_one(
                rfc::core::RelationKey{relations.size() + 1}, block_type,
                {{i + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{occurrence}}});
            relations.push_back({contig, strand, 0, start, end});
        }
        for (auto [start, end] : record.junctions)
            junctions[{record.contig, record.strand, start, end}].push_back(i);
    }
    const auto blocks = relations.size();
    for (const auto &[key, records] : junctions) {
        std::vector<rfc::core::BoundaryOccurrenceInput> occurrences;
        occurrences.reserve(records.size());
        for (auto record : records) {
            occurrence = rfc::machine::checked_size_add(occurrence, 1);
            occurrences.push_back(
                {record + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{occurrence}});
        }
        builder.add_grade_one(rfc::core::RelationKey{relations.size() + 1}, junction_type,
                              occurrences);
        relations.push_back({contigs.indices.at(std::get<0>(key)),
                             strands.indices.at(std::get<1>(key)), 1, std::get<2>(key),
                             std::get<3>(key)});
    }
    rfc::storage::TemporalEndpoint endpoint(builder.freeze());
    const rfc::field::BoundaryOperator boundary(endpoint.source, {1});
    const auto support_count = boundary.transpose_apply(selected);
    endpoint.carriers.push_back({1,
                                 rfc::storage::TemporalCarrierKind::SourceGrade,
                                 {0},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(size)});
    endpoint.carriers.push_back({2,
                                 rfc::storage::TemporalCarrierKind::SourceGrade,
                                 {1},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(relations.size())});
    SourceGradeFields fields{endpoint};
    N::Object records, pieces;
    const auto record_column = [&](const std::string &name, auto value) {
        records.members.emplace_back(name,
                                     fields.column("alignment.record." + name, 0, size, value));
    };
    record_column("flags", [&](auto i) { return Rational(decoded.records[i].flags); });
    record_column("mapping_quality", [&](auto i) { return Rational(decoded.records[i].mapq); });
    record_column("query_length", [&](auto i) {
        return Rational(std::to_string(decoded.records[i].query_length));
    });
    record_column("selected", [&](auto i) { return selected[i]; });
    record_column("contig", [&](auto i) { return Rational(std::to_string(record_contig[i])); });
    record_column("strand", [&](auto i) { return Rational(std::to_string(record_strand[i])); });
    record_column("read_group", [&](auto i) { return Rational(std::to_string(record_group[i])); });
    record_column("query_name_offset",
                  [&](auto i) { return Rational(std::to_string(name_start[i])); });
    record_column("query_name_length",
                  [&](auto i) { return Rational(std::to_string(name_length[i])); });
    const auto mapped_support =
        mapped.size() == size ? rfc::numeric::StructuralSupport::full(size)
                              : (mapped.size() <= (size / 8 + (size % 8 != 0)) / sizeof(std::size_t)
                                     ? rfc::numeric::StructuralSupport::sparse(size, mapped)
                                     : rfc::numeric::StructuralSupport::packed(size, mapped));
    for (const auto *name : {"start", "end"}) {
        records.members.emplace_back(
            name, fields.column(
                      "alignment.record." + std::string(name), 0, size,
                      [&](auto i) {
                          return Rational(std::to_string(std::string_view(name) == "start"
                                                             ? decoded.records[i].start
                                                             : decoded.records[i].end));
                      },
                      mapped_support));
    }
    const auto piece_column = [&](const std::string &name, auto value) {
        pieces.members.emplace_back(
            name, fields.column("alignment.relation." + name, 1, relations.size(), value));
    };
    piece_column("start", [&](auto i) { return Rational(std::to_string(relations[i].start)); });
    piece_column("end", [&](auto i) { return Rational(std::to_string(relations[i].end)); });
    piece_column("kind", [&](auto i) { return Rational(std::to_string(relations[i].kind)); });
    piece_column("contig", [&](auto i) { return Rational(std::to_string(relations[i].contig)); });
    piece_column("strand", [&](auto i) { return Rational(std::to_string(relations[i].strand)); });
    piece_column("selected_record_count", [&](auto i) { return support_count[i]; });
    std::vector<RawByteSource> originals;
    originals.emplace_back(1,
                           std::as_bytes(std::span(resource.bytes.data(), resource.bytes.size())));
    originals.emplace_back(2, std::as_bytes(std::span(names.data(), names.size())));
    N result(N::Object{
        {{"schema", N(std::string("gadget.alignment.fields.v1"))},
         {"format", N(resource.declaration.at("format").text())},
         {"htslib_version", N(hts_codec_version())},
         {"reference", import_record_view(resource.declaration.at("reference"))},
         {"header", N(decoded.header)},
         {"filter", import_record_view(policy)},
         {"record_count", integer(size)},
         {"block_relation_count", integer(blocks)},
         {"junction_relation_count", integer(junctions.size())},
         {"record_fields", N(std::move(records))},
         {"relation_fields", N(std::move(pieces))},
         {"contigs", N(std::move(contigs.labels))},
         {"strands", N(std::move(strands.labels))},
         {"read_groups", N(std::move(groups.labels))},
         {"decoded_query_names", N(N::ByteSpan{2, 0, names.size(), false})},
         {"original", N(N::ByteSpan{1, 0, resource.bytes.size(), false})},
         {"original_sha256", N("sha256:" + sha256(resource.bytes))},
         {"record_identity", N(std::string("Original file and ordered decoded record index"))},
         {"relation_order", N(std::string("Aligned blocks in record and CIGAR order, then "
                                          "junctions in contig, strand, start, end order"))},
         {"junction_quantity",
          N(std::string(
              "CIGAR N occurrences; selected_record_count applies the declared record filter"))}}});
    if (!resource.reference_bytes.empty()) {
        originals.emplace_back(3, std::as_bytes(std::span(resource.reference_bytes.data(),
                                                          resource.reference_bytes.size())));
        std::get<N::Object>(result.data)
            .members.emplace_back("original_reference",
                                  N(N::ByteSpan{3, 0, resource.reference_bytes.size(), false}));
    }
    return {std::move(endpoint), std::move(result), std::move(originals)};
}
AssayFieldsSource feature_count_fields_source(const ResourceImages &images, const Json &request,
                                                  std::shared_ptr<CpuResources> resources) {
    if (!resources)
        resources = CpuResources::local();
    Json evidence;
    std::vector<Vector> native_counts(1);
    {
        auto reservation = resources->reserve(1);
        evidence = feature_allocation_evidence(images, request, &native_counts.front());
    }
    const auto allocation = evidence.at("allocation_request");
    std::optional<AssayFieldsSource> source;
    const auto certificate = allocation_analysis(allocation, resources, &source, &native_counts);
    auto reservation = resources->reserve(1);
    require(source.has_value(), "allocation source missing");
    auto endpoint = std::move(source->endpoint);
    const auto candidates = allocation.at("candidates");
    const auto classes = allocation.at("classes");
    const auto &declaration = std::get<N::Object>(source->declaration.data).members;
    const auto samples = std::find_if(declaration.begin(), declaration.end(),
        [](const auto &entry) { return entry.first == "samples"; });
    require(samples != declaration.end(), "allocation sample declaration missing");
    const auto &samples_value = samples->second;
    const auto &sample = std::get<N::Object>(
        std::get<N::Array>(samples_value.data).elements.front().data).members;
    N::Object quantitative;
    N total;
    for (const auto &[name, value] : sample) {
        if (name == "class_count" || name == "lower" || name == "upper")
            quantitative.members.emplace_back(name, value);
        if (name == "total") total = value;
    }
    std::map<std::string, std::size_t> candidate_ids;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        candidate_ids.emplace(candidates.at(i).text(), i);
    const auto segments = evidence.at("segments");
    N::Object segment_fields;
    const auto column = [&](std::size_t size, auto value) {
        return N(
            SourceGradeFields::tensor(size, value, rfc::numeric::StructuralSupport::full(size)));
    };
    for (const auto *name : {"mate", "class_index"})
        segment_fields.members.emplace_back(name, column(segments.size(), [&](auto i) {
                                                return Rational(segments.at(i).at(name).integer());
                                            }));
    std::vector<std::size_t> starts, lengths, record_starts, record_lengths;
    Vector source_records, group_ids;
    Dictionary groups;
    std::string names;
    for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto segment = segments.at(i);
        starts.push_back(names.size());
        const auto name = segment.at("query").text();
        lengths.push_back(name.size());
        require(name.size() <= names.max_size() - names.size(), "segment query bytes overflow");
        names.append(name);
        group_ids.emplace_back(std::to_string(groups.add(segment.at("read_group").text())));
        const auto records = segment.at("source_record_indices");
        record_starts.push_back(source_records.size());
        record_lengths.push_back(records.size());
        for (std::size_t j = 0; j < records.size(); ++j)
            source_records.emplace_back(records.at(j).integer());
    }
    segment_fields.members.emplace_back("query_name_offset", column(starts.size(), [&](auto i) {
                                            return Rational(std::to_string(starts[i]));
                                        }));
    segment_fields.members.emplace_back("query_name_length", column(lengths.size(), [&](auto i) {
                                            return Rational(std::to_string(lengths[i]));
                                        }));
    segment_fields.members.emplace_back(
        "read_group", column(group_ids.size(), [&](auto i) { return group_ids[i]; }));
    segment_fields.members.emplace_back("source_record_offset",
                                        column(record_starts.size(), [&](auto i) {
                                            return Rational(std::to_string(record_starts[i]));
                                        }));
    segment_fields.members.emplace_back("source_record_count",
                                        column(record_lengths.size(), [&](auto i) {
                                            return Rational(std::to_string(record_lengths[i]));
                                        }));
    segment_fields.members.emplace_back(
        "source_records", column(source_records.size(), [&](auto i) { return source_records[i]; }));
    const auto witnesses = evidence.at("overlap_witnesses");
    N::Object overlaps;
    for (const auto *name : {"source_record_index", "overlap_bases", "aligned_bases", "eligible"})
        overlaps.members.emplace_back(name, column(witnesses.size(), [&](auto i) {
                                          const auto value = witnesses.at(i).at(name);
                                          if (std::string_view(name) == "eligible")
                                              return Rational(value.boolean());
                                          if (std::string_view(name) == "source_record_index")
                                              return Rational(value.integer());
                                          return exact_rational(value.text());
                                      }));
    overlaps.members.emplace_back(
        "feature", column(witnesses.size(), [&](auto i) {
            return Rational(std::to_string(candidate_ids.at(witnesses.at(i).at("feature").text())));
        }));
    if (request.has("junction_policy"))
        overlaps.members.emplace_back("junction_compatible", column(witnesses.size(), [&](auto i) {
                                          return Rational(
                                              witnesses.at(i).at("junction_compatible").boolean());
                                      }));
    const auto excluded = evidence.at("excluded_record_indices");
    N exclusions =
        column(excluded.size(), [&](auto i) { return Rational(excluded.at(i).integer()); });
    auto identity = Json::object();
    for (const auto *name :
         {"schema", "source_digest", "arithmetic", "gadget_sha256", "rexgraph_sha256"})
        identity.set(name, certificate.at(name));
    const auto &alignments = images.at(request.at("resource").text());
    const auto &annotation = images.at(request.at("annotation").text());
    std::vector<RawByteSource> originals;
    originals.emplace_back(
        1, std::as_bytes(std::span(alignments.bytes.data(), alignments.bytes.size())));
    originals.emplace_back(
        2, std::as_bytes(std::span(annotation.bytes.data(), annotation.bytes.size())));
    originals.emplace_back(3, std::as_bytes(std::span(names.data(), names.size())));
    N result(N::Object{
        {{"schema", N(std::string("gadget.feature.count.fields.v1"))},
         {"reference", import_record_view(evidence.at("reference"))},
         {"request", import_record_view(request)},
         {"unit", N(std::string("sequenced_segment"))},
         {"features", import_record_view(evidence.at("features"))},
         {"candidates", import_record_view(candidates)},
         {"classes", import_record_view(classes)},
         {"fields", N(std::move(quantitative))},
         {"segment_count", integer(segments.size())},
         {"segments", N(std::move(segment_fields))},
         {"read_groups", N(std::move(groups.labels))},
         {"overlap_count", integer(witnesses.size())},
         {"overlaps", N(std::move(overlaps))},
         {"excluded_record_indices", std::move(exclusions)},
         {"allocation_fields", std::move(source->declaration)},
         {"allocation_identity", import_record_view(identity)},
         {"allocation_total",
          std::move(total)},
         {"scope", N(evidence.at("scope").text())},
         {"decoded_query_names", N(N::ByteSpan{3, 0, names.size(), false})},
         {"original_alignments", N(N::ByteSpan{1, 0, alignments.bytes.size(), false})},
         {"original_annotation", N(N::ByteSpan{2, 0, annotation.bytes.size(), false})}}});
    if (evidence.has("execution"))
        std::get<N::Object>(result.data)
            .members.emplace_back("execution", import_record_view(evidence.at("execution")));
    if (!alignments.reference_bytes.empty()) {
        originals.emplace_back(4, std::as_bytes(std::span(alignments.reference_bytes.data(),
                                                         alignments.reference_bytes.size())));
        std::get<N::Object>(result.data).members.emplace_back(
            "original_reference", N(N::ByteSpan{4, 0, alignments.reference_bytes.size(), false}));
    }
    return {std::move(endpoint), std::move(result), std::move(originals)};
}
} // namespace gadget::assay_detail
