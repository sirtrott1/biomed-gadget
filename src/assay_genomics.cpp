#include "genomic_detail.hpp"
#include "interval_detail.hpp"
#include "source_grade_fields.hpp"
#include <algorithm>
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::assay_detail {
namespace {
using rfc::numeric::ExactMatrix;
const ResourceImage *image(const ResourceImages &images, Json id) {
    const auto name = text(id);
    require(images.contains(name), "unknown genomic resource");
    return &images.at(name);
}
void same_reference(const ResourceImage &a, const ResourceImage &b) {
    require(reference_record(a.declaration.at("reference")).dump() ==
                reference_record(b.declaration.at("reference")).dump(),
            "genomic resources require the identical declared assembly, sequence and contigs");
}
struct Feature {
    std::string id, contig, strand;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> intervals;
    std::vector<std::size_t> records;
};
std::vector<Feature> annotation(const ResourceImage &r, const Json &request,
                                const std::vector<Interval> *original = nullptr) {
    const auto format = r.declaration.at("format").text();
    require(format == "gtf" || format == "gff3" || format == "bed", "feature annotation format");
    const auto type = text(request.at("feature_type")),
               attribute = text(request.at("group_attribute"));
    require(format != "bed" || (type == "interval" && attribute == "id"),
            "BED features use interval/id");
    const auto decoded = original ? std::vector<Interval>{} : read_intervals(r.bytes, r.declaration);
    const auto &records = original ? *original : decoded;
    std::vector<Feature> out;
    std::map<std::string, std::size_t> indices;
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto &record = records[i];
        if (format != "bed" && record.feature_type != type)
            continue;
        std::vector<std::string> groups;
        if (format == "bed")
            groups.push_back(record.id);
        else {
            require(record.attributes.contains(attribute),
                    "selected annotation record lacks grouping attribute");
            groups = record.attributes.at(attribute);
        }
        std::set<std::string> seen;
        for (const auto &group : groups) {
            require(!group.empty() && group != "__unassigned__" && seen.insert(group).second &&
                        group.find_first_of("\t\r\n") == std::string::npos,
                    "annotation grouping identity");
            if (!indices.contains(group)) {
                indices.emplace(group, out.size());
                out.push_back({group, record.contig, record.strand, {}, {}});
            }
            auto &feature = out[indices.at(group)];
            require(
                feature.contig == record.contig && feature.strand == record.strand,
                "feature identity spans different contigs or strands; declare distinct features");
            feature.intervals.emplace_back(record.start, record.end);
            feature.records.push_back(i);
        }
    }
    require(!out.empty(), "selected annotation features must be nonempty");
    for (auto &feature : out) {
        std::sort(feature.intervals.begin(), feature.intervals.end());
        std::vector<std::pair<std::uint64_t, std::uint64_t>> merged;
        for (const auto &p : feature.intervals)
            if (!merged.empty() && p.first <= merged.back().second)
                merged.back().second = std::max(merged.back().second, p.second);
            else
                merged.push_back(p);
        feature.intervals = std::move(merged);
    }
    return out;
}
Json feature_json(const Feature &f) {
    auto out = Json::object(), segments = Json::array(), records = Json::array();
    out.set("id", Json::string(f.id));
    out.set("contig", Json::string(f.contig));
    out.set("strand", Json::string(f.strand));
    for (auto [a, b] : f.intervals) {
        auto p = Json::array();
        p.push(cardinality(a));
        p.push(cardinality(b));
        segments.push(p);
    }
    for (auto i : f.records)
        records.push(cardinality(i));
    out.set("union_intervals", segments);
    out.set("source_record_indices", records);
    return out;
}
Json complex_receipt(const rfc::core::RelationalComplex &source) {
    const auto archive = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(source), 1, source_codec_limits(source_byte_limit())));
    auto out = Json::object();
    out.set(
        "sha256",
        Json::string("sha256:" + sha256(std::string(reinterpret_cast<const char *>(archive.data()),
                                                    archive.size()))));
    out.set("archive_bytes", cardinality(archive.size()));
    return out;
}
bool retained(const Alignment &a, const Json &policy) {
    members(policy, {"exclude_flags", "minimum_mapping_quality", "unknown_mapping_quality"});
    const auto exclude = count(policy.at("exclude_flags"), 4095, "alignment excluded flags");
    const auto min = count(policy.at("minimum_mapping_quality"), 254, "alignment minimum MAPQ");
    const auto unknown = policy.at("unknown_mapping_quality").text();
    require(unknown == "include" || unknown == "exclude", "MAPQ 255 handling");
    return !(a.flags & exclude) && (a.mapq == 255 ? unknown == "include" : a.mapq >= min);
}
Json alignment_evidence(const ResourceImage &r, const Json &request) {
    members(request, {"schema", "operation", "resource", "filter"});
    auto out = inspect_hts(r), selected = Json::array(), excluded = Json::array(),
         junctions = Json::array();
    using Key = std::tuple<std::string, std::string, std::uint64_t, std::uint64_t>;
    std::map<Key, std::vector<std::size_t>> support;
    for (std::size_t i = 0; i < out.at("records").size(); ++i) {
        const auto record = out.at("records").at(i);
        Alignment a;
        a.flags = static_cast<std::uint16_t>(record.at("flags").integer());
        a.mapq = static_cast<std::uint8_t>(record.at("mapping_quality").integer());
        if (!retained(a, request.at("filter"))) {
            excluded.push(cardinality(i));
            continue;
        }
        selected.push(cardinality(i));
        for (std::size_t j = 0; j < record.at("splice_junctions").size(); ++j) {
            const auto p = record.at("splice_junctions").at(j);
            support[{record.at("contig").text(), record.at("strand").text(),
                     static_cast<std::uint64_t>(p.at(0).integer()),
                     static_cast<std::uint64_t>(p.at(1).integer())}]
                .push_back(i);
        }
    }
    ExactMatrix map(support.size(), out.at("records").size());
    std::size_t row = 0;
    for (const auto &[key, indices] : support) {
        auto item = Json::object(), ids = Json::array();
        item.set("contig", Json::string(std::get<0>(key)));
        item.set("alignment_strand", Json::string(std::get<1>(key)));
        item.set("intron_start", cardinality(std::get<2>(key)));
        item.set("intron_end", cardinality(std::get<3>(key)));
        for (auto i : indices) {
            map.set(row, i, map.get(row, i) + 1);
            ids.push(cardinality(i));
        }
        item.set("source_record_indices", ids);
        junctions.push(item);
        ++row;
    }
    const auto counts = map.apply(Vector(out.at("records").size(), 1));
    for (std::size_t i = 0; i < junctions.size(); ++i)
        junctions.at(i).set("alignment_record_count", Json::string(counts[i].get_str()));
    out.set("selected_record_indices", selected);
    out.set("excluded_record_indices", excluded);
    out.set("junctions", junctions);
    out.set("junction_quantity", Json::string("Retained CIGAR N occurrences by alignment strand; "
                                              "multiple mappings remain separate records"));
    return out;
}
Json feature_evidence(const ResourceImages &images, const Json &request,
                      Vector *native_counts = nullptr) {
    members(request,
            {"schema", "operation", "resource", "annotation", "feature_type", "group_attribute",
             "unit", "overlap", "strand_policy", "filter", "candidate_search"},
            {"junction_policy"});
    require(request.at("unit").text() == "sequenced_segment",
            "counting currently admits sequenced_segment units; paired templates need a declared "
            "mate model");
    const auto overlap = request.at("overlap").text(), strand = request.at("strand_policy").text();
    require(overlap == "any_aligned_base" || overlap == "all_aligned_bases",
            "feature overlap rule");
    require(strand == "ignore" || strand == "same" || strand == "opposite", "counting strand rule");
    const auto &r = *image(images, request.at("resource")),
               &ann = *image(images, request.at("annotation"));
    same_reference(r, ann);
    const auto features = annotation(ann, request);
    const auto junction_policy =
        request.has("junction_policy") ? text(request.at("junction_policy")) : "ignore";
    require(junction_policy == "ignore" || junction_policy == "declared_exon_chain",
            "feature junction policy");
    require(junction_policy != "declared_exon_chain" || overlap == "all_aligned_bases",
            "exon-chain compatibility requires all aligned bases");
    struct IndexedInterval {
        std::uint64_t start, end;
        std::size_t feature;
    };
    struct Index {
        std::vector<IndexedInterval> intervals;
        std::vector<std::uint64_t> prefix_end;
    };
    std::map<std::string, Index> indices;
    for (std::size_t f = 0; f < features.size(); ++f)
        for (auto [start, end] : features[f].intervals)
            indices[features[f].contig].intervals.push_back({start, end, f});
    for (auto &[contig, index] : indices) {
        (void)contig;
        std::stable_sort(index.intervals.begin(), index.intervals.end(),
                         [](const auto &a, const auto &b) { return a.start < b.start; });
        std::uint64_t end = 0;
        for (const auto &interval : index.intervals) {
            end = std::max(end, interval.end);
            index.prefix_end.push_back(end);
        }
    }
    const auto alignments = decode_alignments(r, false);
    struct Segment {
        std::set<std::size_t> candidates;
        std::vector<std::size_t> records;
    };
    using Key = std::tuple<std::string, std::string, int>;
    std::map<Key, Segment> segments;
    auto excluded = Json::array(), overlap_records = Json::array();
    std::size_t comparisons = 0, incidences = 0;
    for (std::size_t i = 0; i < alignments.records.size(); ++i) {
        const auto &a = alignments.records[i];
        if (!retained(a, request.at("filter"))) {
            excluded.push(cardinality(i));
            continue;
        }
        require(a.query != "*", "anonymous query names cannot define segment identity");
        require(!(a.flags & 1) || (a.flags & (64 | 128)),
                "paired segment lacks first/second identity");
        auto &segment = segments[{a.group, a.query, (a.flags & 64) ? 1 : (a.flags & 128) ? 2 : 0}];
        segment.records.push_back(i);
        bool assigned = false;
        std::set<std::size_t> possible;
        if (indices.contains(a.contig)) {
            const auto &index = indices.at(a.contig);
            for (auto [start, end] : a.blocks) {
                const auto first =
                    std::upper_bound(index.prefix_end.begin(), index.prefix_end.end(), start);
                const auto last = std::lower_bound(
                    index.intervals.begin(), index.intervals.end(), end,
                    [](const auto &interval, auto position) { return interval.start < position; });
                for (auto k = static_cast<std::size_t>(first - index.prefix_end.begin());
                     k < static_cast<std::size_t>(last - index.intervals.begin()); ++k)
                    if (index.intervals[k].end > start)
                        possible.insert(index.intervals[k].feature);
            }
        }
        for (auto f : possible) {
            comparisons = rfc::machine::checked_size_add(comparisons, 1);
            const auto &feature = features[f];
            if (a.contig != feature.contig ||
                (strand != "ignore" &&
                 (feature.strand == "." || ((strand == "same") != (a.strand == feature.strand)))))
                continue;
            Vector lengths, block_lengths;
            for (auto [begin, end] : a.blocks) {
                block_lengths.push_back(exact_rational(std::to_string(end - begin)));
                for (auto [left, right] : feature.intervals) {
                    const auto lo = std::max(begin, left), hi = std::min(end, right);
                    if (lo < hi)
                        lengths.push_back(exact_rational(std::to_string(hi - lo)));
                }
            }
            const auto common = rfc::numeric::dot(lengths, Vector(lengths.size(), 1));
            const auto bases = rfc::numeric::dot(block_lengths, Vector(block_lengths.size(), 1));
            if (common > 0) {
                auto witness = Json::object();
                witness.set("source_record_index", cardinality(i));
                witness.set("feature", Json::string(feature.id));
                witness.set("overlap_bases", Json::string(common.get_str()));
                witness.set("aligned_bases", Json::string(bases.get_str()));
                bool junction_compatible = true;
                if (junction_policy == "declared_exon_chain") {
                    std::set<std::pair<std::uint64_t, std::uint64_t>> gaps;
                    for (std::size_t k = 1; k < feature.intervals.size(); ++k)
                        gaps.emplace(feature.intervals[k - 1].second, feature.intervals[k].first);
                    for (auto gap : a.junctions)
                        junction_compatible = junction_compatible && gaps.contains(gap);
                    const std::set<std::pair<std::uint64_t, std::uint64_t>> observed(
                        a.junctions.begin(), a.junctions.end());
                    for (auto gap : gaps)
                        if (a.start < gap.second && a.end > gap.first)
                            junction_compatible = junction_compatible && observed.contains(gap);
                    // Terminal deletions may not extend beyond the declared exon span.
                    const auto contains = [&](auto base) {
                        return std::any_of(
                            feature.intervals.begin(), feature.intervals.end(), [&](auto interval) {
                                return interval.first <= base && base < interval.second;
                            });
                    };
                    junction_compatible = junction_compatible && a.end > a.start &&
                                          contains(a.start) && contains(a.end - 1);
                }
                const bool eligible =
                    (overlap == "any_aligned_base" || common == bases) && junction_compatible;
                witness.set("eligible", Json::boolean(eligible));
                if (request.has("junction_policy"))
                    witness.set("junction_compatible", Json::boolean(junction_compatible));
                overlap_records.push(witness);
                incidences = rfc::machine::checked_size_add(incidences, 1);
                if (eligible) {
                    segment.candidates.insert(f);
                    assigned = true;
                }
            }
        }
        // An alternative retained mapping with no annotation overlap is a real
        // possible unassigned outcome, even if another mapping overlaps a gene.
        if (!assigned)
            segment.candidates.insert(features.size());
    }
    require(!segments.empty(), "no segment observations remain after the declared filter");
    std::map<std::vector<std::size_t>, std::size_t> classes;
    std::vector<std::size_t> counts;
    auto segment_records = Json::array();
    for (const auto &[key, s] : segments) {
        const std::vector<std::size_t> participants(s.candidates.begin(), s.candidates.end());
        if (!classes.contains(participants)) {
            classes.emplace(participants, counts.size());
            counts.push_back(0);
        }
        const auto c = classes.at(participants);
        ++counts[c];
        auto record = Json::object(), ids = Json::array();
        record.set("read_group", Json::string(std::get<0>(key)));
        record.set("query", Json::string(std::get<1>(key)));
        record.set("mate", cardinality(std::get<2>(key)));
        record.set("class_index", cardinality(c));
        for (auto i : s.records)
            ids.push(cardinality(i));
        record.set("source_record_indices", ids);
        segment_records.push(record);
    }
    auto allocation = Json::object(), candidates = Json::array(), targets = Json::array(),
         class_json = Json::array();
    auto feature_records = Json::array();
    for (const auto &f : features) {
        candidates.push(Json::string(f.id));
        feature_records.push(feature_json(f));
        auto target = Json::object(), ids = Json::array();
        ids.push(Json::string(f.id));
        target.set("id", Json::string(f.id));
        target.set("candidates", ids);
        targets.push(target);
    }
    candidates.push(Json::string("__unassigned__"));
    std::vector<Json> ordered(classes.size());
    for (const auto &[p, index] : classes) {
        auto c = Json::object(), ids = Json::array();
        c.set("id", Json::string("class:" + std::to_string(index)));
        for (auto f : p)
            ids.push(candidates.at(f));
        c.set("candidates", ids);
        ordered[index] = c;
    }
    for (const auto &c : ordered)
        class_json.push(c);
    auto sample = Json::object(), sample_counts = Json::array(), samples = Json::array(),
         evidence = Json::array();
    for (auto n : counts)
        sample_counts.push(Json::string(std::to_string(n)));
    if (native_counts) {
        native_counts->clear();
        for (auto n : counts) native_counts->emplace_back(std::to_string(n));
    }
    sample.set("id", request.at("resource"));
    sample.set("counts", sample_counts);
    sample.set("scale", Json::string("1"));
    sample.set("scale_provenance",
               Json::string("Declared identity scale for retained sequenced-segment counts"));
    samples.push(sample);
    evidence.push(Json::string("sha256:" + sha256(r.bytes)));
    evidence.push(Json::string("sha256:" + sha256(ann.bytes)));
    allocation.set("schema", Json::string("gadget.assay.allocations.v1"));
    allocation.set("operation", Json::string("allocation_bounds"));
    allocation.set("candidates", candidates);
    allocation.set("classes", class_json);
    allocation.set("samples", samples);
    allocation.set("targets", targets);
    allocation.set("contrasts", Json::array());
    allocation.set("evidence", evidence);
    allocation.set("units", Json::string("assigned_segment_count"));
    allocation.set("class_partition", Json::string("disjoint_observations"));
    const auto search = request.at("candidate_search");
    require(search.at("reference_sha256").text() ==
                r.declaration.at("reference").at("sequence_sha256").text(),
            "counting candidate reference witness");
    allocation.set("candidate_search", search);
    auto out = Json::object();
    out.set("features", feature_records);
    out.set("segments", segment_records);
    out.set("excluded_record_indices", excluded);
    out.set("overlap_witnesses", overlap_records);
    out.set("allocation_request", allocation);
    out.set("reference", r.declaration.at("reference"));
    out.set("unit", request.at("unit"));
    out.set("scope",
            Json::string("Sharp marginal bounds for one assignment per retained RG/query/mate "
                         "segment over the retained alignment/annotation candidates, including "
                         "unassigned alternatives. Paired ends are separate segments. Upstream "
                         "alignment completeness and filters are declared evidence."));
    if (request.has("junction_policy")) {
        out.set("junction_policy", request.at("junction_policy"));
        auto execution = Json::object();
        execution.set("candidate_feature_comparisons", cardinality(comparisons));
        execution.set("overlap_witnesses", cardinality(incidences));
        execution.set("interval_index",
                      Json::string("Complete contig/start intervals with exact prefix-end "
                                   "exclusion; original feature order"));
        out.set("execution", execution);
    }
    return out;
}

Json feature_counts(const ResourceImages &images, const Json &request,
                    std::shared_ptr<CpuResources> budget) {
    Json out;
    {
        auto reservation = budget->reserve(1);
        out = feature_evidence(images, request);
    }
    out.set("allocation_certificate", analyze_allocations(out.at("allocation_request"), budget));
    return out;
}

Json chromatin_contrast(const ResourceImages &images, const Json &request,
                        std::optional<AssayFieldsSource> *retained) {
    members(request, {"schema", "operation", "samples", "contrasts", "strand_policy", "thresholds",
                      "absence_policy"});
    require(request.at("absence_policy").text() == "measured_zero_within_union",
            "chromatin absence policy must explicitly declare measured coverage zeros");
    const auto strand = request.at("strand_policy").text();
    require(strand == "ignore" || strand == "separate", "chromatin strand policy");
    const auto samples = request.at("samples");
    require(samples.size() > 0, "coverage samples must be nonempty");
    using interval_detail::Record;
    std::vector<Record> records;
    std::map<std::string, std::size_t> ids;
    Vector scales;
    const ResourceImage *reference = nullptr;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto sample = samples.at(s);
        members(sample, {"id", "resource", "scale", "scale_provenance"});
        require(ids.emplace(text(sample.at("id"), source_byte_limit()), s).second,
                "duplicate coverage sample");
        const auto q = exact_rational(text(sample.at("scale"), source_byte_limit()));
        require(q > 0, "coverage scale must be positive");
        (void)text(sample.at("scale_provenance"), source_byte_limit());
        scales.push_back(q);
        const auto &r = *image(images, sample.at("resource"));
        if (reference)
            same_reference(*reference, r);
        else
            reference = &r;
        const auto intervals = read_intervals(r.bytes, r.declaration);
        for (std::size_t i = 0; i < intervals.size(); ++i) {
            require(intervals[i].weight >= 0, "chromatin coverage weights must be nonnegative");
            records.push_back({intervals[i], s, i});
        }
    }
    const interval_detail::Overlay overlay(records, strand);
    const auto &atoms = overlay.atoms;
    const auto coverage = overlay.coverage(records.size());
    ExactMatrix weights(records.size(), samples.size());
    Vector lengths;
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("joint.physical.coverage.membership");
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(rfc::core::Grade{1});
    for (std::size_t i = 0; i < records.size(); ++i) {
        builder.declare_witness_identity(i + 1);
        weights.set(i, records[i].sample, records[i].interval.weight);
    }
    std::uint64_t occurrence = 0;
    auto segments = Json::array();
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto &a = atoms[i];
        std::vector<rfc::core::BoundaryOccurrenceInput> participants;
        auto source = Json::array();
        for (auto j : a.records) {
            participants.push_back(
                {j + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{++occurrence}});
            if (!retained) {
                auto item = Json::object();
                item.set("sample", samples.at(records[j].sample).at("id"));
                item.set("source_record_index", cardinality(records[j].index));
                source.push(item);
            }
        }
        builder.add_grade_one(rfc::core::RelationKey{rfc::machine::checked_size_add(1000000, i)},
                              type, participants);
        lengths.push_back(exact_rational(std::to_string(a.end - a.start)));
        if (!retained) {
            auto item = Json::object();
            item.set("contig", Json::string(a.frame.first));
            item.set("strand", Json::string(a.frame.second));
            item.set("start", cardinality(a.start));
            item.set("end", cardinality(a.end));
            item.set("memberships", source);
            segments.push(item);
        }
    }
    const auto measured = coverage.multiply(weights);
    using N = NativeRecordValue;
    rfc::storage::TemporalEndpoint endpoint(builder.freeze());
    for (std::size_t k = 0; k < endpoint.source.grade_count(); ++k)
        endpoint.carriers.push_back(
            {k + 1,
             rfc::storage::TemporalCarrierKind::SourceGrade,
             {static_cast<std::uint16_t>(k)},
             {},
             0,
             rfc::numeric::IndexDomain::linear(
                 endpoint.source.grade({static_cast<std::uint16_t>(k)}).cell_count())});
    SourceGradeFields store{endpoint};
    N::Array native_samples, native_contrasts, field_declarations;
    const auto field = [&](std::string id, const Vector &values, bool contrast = true) {
        auto key = store.field(id, 1, values);
        field_declarations.elements.emplace_back(N::Object{
            {{"id", N(id)},
             {"key", key},
             {"quantity",
              N(std::string(contrast ? "scaled_coverage_contrast" : "scaled_coverage"))},
             {"units",
              N(std::string("declared interval weights divided by declared sample scale"))},
             {"method", N(std::string("Declared coverage overlay and sample scales"))}}});
        return key;
    };
    auto sample_fields = Json::array();
    std::vector<Vector> fields;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        Vector v;
        for (std::size_t i = 0; i < atoms.size(); ++i)
            v.push_back(measured.get(i, s));
        v = rfc::numeric::scale(Rational{1} / scales[s], v);
        fields.push_back(v);
        auto item = Json::object(), values = Json::array();
        if (!retained)
            for (const auto &x : v)
                values.push(Json::string(x.get_str()));
        item.set("sample", samples.at(s).at("id"));
        if (!retained)
            item.set("values", values);
        item.set("integral", Json::string(rfc::numeric::dot(lengths, v).get_str()));
        if (retained)
            native_samples.elements.emplace_back(
                N::Object{{{"field", field("sample:" + std::to_string(s), v, false)},
                           {"summary", import_record_view(item)}}});
        sample_fields.push(item);
    }
    const auto contrasts = request.at("contrasts");
    require(contrasts.size() > 0, "chromatin contrasts must be nonempty");
    const auto thresholds = request.at("thresholds");
    require(thresholds.size() > 0, "chromatin thresholds must be nonempty");
    auto results = Json::array();
    std::set<std::string> contrast_ids;
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto contrast = contrasts.at(c);
        members(contrast, {"id", "treated", "control"});
        require(contrast_ids.insert(text(contrast.at("id"), source_byte_limit())).second &&
                    ids.contains(text(contrast.at("treated"), source_byte_limit())) &&
                    ids.contains(text(contrast.at("control"), source_byte_limit())),
                "coverage contrast identities");
        const auto v = rfc::numeric::sub(fields.at(ids.at(contrast.at("treated").text())),
                                         fields.at(ids.at(contrast.at("control").text())));
        auto item = Json::object(), values = Json::array(), components = Json::array();
        if (!retained)
            for (const auto &x : v)
                values.push(Json::string(x.get_str()));
        std::set<Rational> unique;
        for (std::size_t t = 0; t < thresholds.size(); ++t) {
            const auto threshold = exact_rational(text(thresholds.at(t), source_byte_limit()));
            require(unique.insert(threshold).second, "duplicate contrast threshold");
            auto component = interval_detail::threshold(atoms, lengths, v, threshold, false);
            components.push(component);
        }
        item.set("id", contrast.at("id"));
        if (!retained)
            item.set("values", values);
        item.set("integral", Json::string(rfc::numeric::dot(lengths, v).get_str()));
        item.set("quadrance", Json::string(rfc::numeric::dot(v, v).get_str()));
        item.set("threshold_components", components);
        if (retained)
            native_contrasts.elements.emplace_back(
                N::Object{{{"field", field(contrast.at("id").text(), v)},
                           {"summary", import_record_view(item)}}});
        results.push(item);
    }
    auto out = Json::object();
    if (!retained)
        out.set("segments", segments);
    out.set("sample_fields", sample_fields);
    out.set("contrasts", results);
    out.set("reference", reference->declaration.at("reference"));
    out.set("native_complex", complex_receipt(endpoint.source));
    out.set("quantity", Json::string("scaled_coverage_contrast"));
    out.set("units", Json::string("declared interval weights divided by declared sample scale"));
    out.set(
        "scope",
        Json::string("Complete coordinate fields on the union of retained physical intervals, with "
                     "caller-declared measured zeros. Thresholds are observed coverage criteria; "
                     "no statistical peak null or region outside the union is supplied."));
    if (retained) {
        endpoint.carriers.push_back({3,
                                     rfc::storage::TemporalCarrierKind::Explicit,
                                     {},
                                     {},
                                     1,
                                     rfc::numeric::IndexDomain::linear(samples.size())});
        const auto length_key = store.field("segment.length", 1, lengths);
        const auto coverage_key = store.map("coverage_map", 1, 2, coverage);
        const auto weights_key = store.map("source_weights", 3, 1, weights);
        const auto scales_key =
            store.on_carrier("sample.scales", 3,
                             SourceGradeFields::tensor(
                                 scales.size(), [&](auto i) { return scales[i]; },
                                 rfc::numeric::StructuralSupport::full(scales.size())),
                             rfc::storage::TemporalTensorRole::Field);
        N::Object records_value, segments_value;
        const auto record_column = [&](std::string name, auto value) {
            records_value.members.emplace_back(
                name, store.column("record." + name, 0, records.size(), value));
        };
        record_column("sample", [&](auto i) { return Rational(records[i].sample); });
        record_column("record", [&](auto i) { return Rational(records[i].index); });
        record_column("source_start",
                      [&](auto i) { return Rational(records[i].interval.source_start); });
        record_column("source_end",
                      [&](auto i) { return Rational(records[i].interval.source_end); });
        std::map<std::pair<std::string, std::string>, std::size_t> frames;
        N::Array frame_values;
        for (const auto &atom : atoms)
            if (!frames.contains(atom.frame)) {
                frames.emplace(atom.frame, frame_values.elements.size());
                frame_values.elements.emplace_back(
                    N::Object{{{"contig", N(atom.frame.first)}, {"strand", N(atom.frame.second)}}});
            }
        const auto segment_column = [&](std::string name, auto value) {
            segments_value.members.emplace_back(
                name, store.column("segment." + name, 1, atoms.size(), value));
        };
        segment_column("frame", [&](auto i) { return Rational(frames.at(atoms[i].frame)); });
        segment_column("start", [&](auto i) { return Rational(atoms[i].start); });
        segment_column("end", [&](auto i) { return Rational(atoms[i].end); });
        N declaration(N::Object{{{"schema", N(std::string("gadget.chromatin.fields.v1"))},
                                 {"fields", N(std::move(field_declarations))},
                                 {"sample_fields", N(std::move(native_samples))},
                                 {"contrasts", N(std::move(native_contrasts))},
                                 {"coverage_map", coverage_key},
                                 {"source_weights", weights_key},
                                 {"sample_scales", scales_key},
                                 {"segment_lengths", length_key},
                                 {"segments", N(std::move(segments_value))},
                                 {"frames", N(std::move(frame_values))},
                                 {"records", N(std::move(records_value))},
                                 {"request", import_record_view(request)},
                                 {"summary", import_record_view(out)}}});
        retained->emplace(AssayFieldsSource{std::move(endpoint), std::move(declaration), {}});
    }
    return out;
}
} // namespace
Json feature_allocation_evidence(const ResourceImages &images, const Json &request,
                                 Vector *native_counts) {
    return feature_evidence(images, request, native_counts);
}
AssayFieldsSource annotation_fields_source(const ResourceImage &resource, const Json &request) {
    using N = NativeRecordValue;
    members(request, {"schema", "operation", "resource", "feature_type", "group_attribute"});
    const auto records = read_intervals(resource.bytes, resource.declaration);
    const auto features = annotation(resource, request, &records);
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("annotation.feature.record.membership");
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(rfc::core::Grade{1});
    for (std::size_t i = 0; i < records.size(); ++i)
        builder.declare_witness_identity(i + 1);
    std::uint64_t occurrence = 1;
    N::Array groups, original;
    for (std::size_t i = 0; i < features.size(); ++i) {
        const auto &feature = features[i];
        std::vector<rfc::core::BoundaryOccurrenceInput> members;
        for (auto index : feature.records)
            members.push_back({index + 1, true, false, Rational{1}, {},
                               rfc::core::OccurrenceKey{occurrence++}});
        builder.add_grade_one(rfc::core::RelationKey{i + 1}, type, members);
        N::Array intervals;
        for (const auto &[start, end] : feature.intervals)
            intervals.elements.emplace_back(N::Array{{
                N(rfc::numeric::Integer(std::to_string(start))),
                N(rfc::numeric::Integer(std::to_string(end)))}});
        groups.elements.emplace_back(N::Object{{
            {"id", N(feature.id)}, {"contig", N(feature.contig)}, {"strand", N(feature.strand)},
            {"union_intervals", N(std::move(intervals))}}});
    }
    for (const auto &record : records)
        original.elements.emplace_back(N::Object{{
            {"id", N(record.id)}, {"contig", N(record.contig)}, {"strand", N(record.strand)},
            {"feature_type", N(record.feature_type)},
            {"bytes", N(N::ByteSpan{1, record.source_start,
                                    record.source_end - record.source_start, false})}}});
    rfc::storage::TemporalEndpoint endpoint(builder.freeze());
    endpoint.carriers = {
        {1, rfc::storage::TemporalCarrierKind::SourceGrade, {0}, {}, 0,
         rfc::numeric::IndexDomain::linear(records.size())},
        {2, rfc::storage::TemporalCarrierKind::SourceGrade, {1}, {}, 0,
         rfc::numeric::IndexDomain::linear(features.size())}};
    SourceGradeFields store(endpoint);
    N declaration(N::Object{{
        {"schema", N(std::string("gadget.annotation.features.fields.v1"))},
        {"features", N(std::move(groups))}, {"original_records", N(std::move(original))},
        {"reference", import_record_view(resource.declaration.at("reference"))},
        {"source_bytes", N(N::ByteSpan{1, 0, resource.bytes.size(), false})},
        {"record_start", store.column("annotation.start", 0, records.size(), [&](auto i) { return Rational(records[i].start); })},
        {"record_end", store.column("annotation.end", 0, records.size(), [&](auto i) { return Rational(records[i].end); })},
        {"group_record_count", store.column("feature.record_count", 1, features.size(), [&](auto i) { return Rational(features[i].records.size()); })}}});
    return {std::move(endpoint), std::move(declaration),
            {RawByteSource(1, std::as_bytes(std::span(resource.bytes)))}};
}
Json annotation_fields_view(const NativeRecord &record) {
    auto out = Json::object(), features = Json::array();
    const auto declarations = record.at("features");
    const auto &source = record.endpoint()->source.grade(rfc::core::Grade{1});
    for (std::size_t i = 0; i < source.cell_count(); ++i) {
        auto feature = render_record_view(declarations.at(i)), indices = Json::array();
        const auto range = source.incidence_range(i);
        for (auto j = range.first; j < range.second; ++j)
            indices.push(cardinality(source.lower_index.at(j)));
        feature.set("source_record_indices", indices);
        features.push(feature);
    }
    out.set("features", features);
    out.set("reference", render_record_view(record.at("reference")));
    return out;
}
Json genomic_analysis(const ResourceImages &images, const Json &request,
                      std::shared_ptr<CpuResources> budget,
                      std::optional<AssayFieldsSource> *retained) {
    const auto op = request.at("operation").text();
    if (op == "feature_counts")
        return feature_counts(images, request, budget);
    if (op == "expression_allocations")
        return expression_allocations(images, request, budget);
    auto reservation = budget->reserve(1);
    if (op == "indexed_sequence_alignment")
        return indexed_alignment(images, request);
    if (op == "sequence_splice_discovery")
        return splice_search(images, request);
    if (op == "read_overlap_reconstruction")
        return read_reconstruction(*image(images, request.at("resource")), request);
    if (op == "reference_sequence_matches")
        return reference_matches(images, request);
    if (op == "sequence_assembly_graph")
        return assembly_graph(*image(images, request.at("resource")), request, retained);
    if (op == "exact_sequence_alignment")
        return alignment_analysis(images, request, retained);
    if (op == "alignment_evidence")
        return alignment_evidence(*image(images, request.at("resource")), request);
    if (op == "variant_genotypes")
        return genotype_analysis(*image(images, request.at("resource")), request, retained);
    if (op == "read_quality")
        return quality_analysis(*image(images, request.at("resource")), request);
    if (op == "controlled_coverage")
        return controlled_coverage(images, request);
    if (op == "chromatin_contrast")
        return chromatin_contrast(images, request, retained);
    if (op == "annotation_features") {
        const auto &r = *image(images, request.at("resource"));
        auto fields = annotation_fields_source(r, request);
        if (retained) {
            retained->emplace(std::move(fields));
            return Json{};
        }
        return annotation_fields_view(NativeRecord::create_from_endpoint(
            std::move(fields.declaration), new_native_record_namespace(),
            std::move(fields.sources), std::move(fields.endpoint)));
    }
    throw std::invalid_argument("unsupported genomic query");
}
Json chromatin_fields_view(const NativeRecord &record) {
    auto out = render_record_view(record.at("summary"));
    const auto coordinate = [&](const NativeRecord &key, std::size_t i) {
        const auto value = retained_tensor(record, key).tensor.state(i).value();
        require(value >= 0 && value.get_den() == 1, "chromatin source coordinate column");
        return value.get_num().get_ui();
    };
    auto segments = Json::array();
    const auto columns = record.at("segments"), rows = record.at("records");
    const auto endpoint = record.endpoint();
    const auto &arena = endpoint->source.grade(rfc::core::Grade{1});
    for (std::size_t i = 0; i < arena.cell_count(); ++i) {
        auto segment =
            render_record_view(record.at("frames").at(coordinate(columns.at("frame"), i)));
        segment.set("start", cardinality(coordinate(columns.at("start"), i)));
        segment.set("end", cardinality(coordinate(columns.at("end"), i)));
        auto memberships = Json::array();
        const auto range = arena.incidence_range(i);
        for (std::size_t slot = range.first; slot < range.second; ++slot) {
            const auto lower = arena.lower_index.at(slot);
            auto member = Json::object();
            member.set("sample", render_record_view(record.at("request")
                                                        .at("samples")
                                                        .at(coordinate(rows.at("sample"), lower))
                                                        .at("id")));
            member.set("source_record_index", cardinality(coordinate(rows.at("record"), lower)));
            memberships.push(member);
        }
        segment.set("memberships", memberships);
        segments.push(segment);
    }
    out.set("segments", segments);
    for (const auto *name : {"sample_fields", "contrasts"}) {
        auto array = Json::array();
        const auto entries = record.at(name);
        for (std::size_t i = 0; i < entries.size(); ++i) {
            auto value = render_record_view(entries.at(i).at("summary"));
            const auto key = std::stoull(entries.at(i).at("field").integer().get_str());
            const auto endpoint = record.endpoint();
            const auto tensor = std::find_if(endpoint->tensors.begin(), endpoint->tensors.end(),
                                             [&](const auto &field) { return field.key == key; });
            require(tensor != endpoint->tensors.end(), "chromatin field owner absent");
            auto values = Json::array();
            for (std::size_t j = 0; j < tensor->tensor.domain().logical_size(); ++j)
                values.push(Json::string(tensor->tensor.state(j).value().get_str()));
            value.set("values", values);
            array.push(value);
        }
        out.set(name, array);
    }
    return out;
}
} // namespace gadget::assay_detail
