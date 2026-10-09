#include "genomic_detail.hpp"
#include "assay_quality_detail.hpp"
#include <algorithm>
#include <array>
#include <numeric>
#include <rfc/core/builder.hpp>
#include <rfc/field/exact_metric.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::assay_detail {
Json spatial_analysis(const Json &assay, const Json &request) {
    members(request,
            {"schema", "operation", "assay", "coordinates", "frame", "units", "metric_diagonal",
             "radius_squared", "relation_type", "evidence"},
            {"execution_certificate"});
    const auto columns = request.at("coordinates"), diagonal = request.at("metric_diagonal");
    require(columns.size() > 0 && columns.size() == diagonal.size(),
            "spatial coordinate dimensions");
    (void)text(request.at("frame"), source_byte_limit());
    (void)text(request.at("units"), source_byte_limit());
    (void)text(request.at("relation_type"), source_byte_limit());
    (void)text(request.at("evidence"), source_byte_limit());
    Vector weights;
    std::set<std::string> unique;
    for (std::size_t i = 0; i < columns.size(); ++i) {
        require(unique.insert(text(columns.at(i), source_byte_limit())).second,
                "duplicate spatial coordinate column");
        weights.push_back(exact_rational(text(diagonal.at(i), source_byte_limit())));
    }
    const auto metric = rfc::field::ExactMetric::diagonal(weights);
    const auto radius = exact_rational(text(request.at("radius_squared"), source_byte_limit()));
    require(radius >= 0, "spatial radius squared");
    const auto observations = assay.at("observations");
    require(observations.size() > 0, "spatial observations must be nonempty");
    std::vector<Vector> positions;
    for (std::size_t i = 0; i < observations.size(); ++i) {
        const auto row = observations.at(i);
        Vector v;
        require(row.has("coordinate_frame") &&
                    row.at("coordinate_frame").text() == request.at("frame").text(),
                "spatial observation lacks the selected coordinate frame");
        require(row.has("coordinate_units") &&
                    row.at("coordinate_units").text() == request.at("units").text(),
                "spatial observation coordinate units");
        for (std::size_t c = 0; c < columns.size(); ++c)
            v.push_back(exact_rational(row.at(columns.at(c).text()).text()));
        positions.push_back(std::move(v));
    }
    // A positive diagonal form bounds every coordinate contribution by the
    // complete quadrance. Sorted windows therefore exclude only impossible pairs.
    // Choose the widest weighted axis without square roots or fitted coordinates.
    std::size_t axis = 0;
    Rational extent = -1;
    for (std::size_t d = 0; d < columns.size(); ++d) {
        auto low = positions.front()[d], high = low;
        for (const auto &p : positions) {
            low = std::min(low, p[d]);
            high = std::max(high, p[d]);
        }
        const Vector delta{high - low};
        const auto q = rfc::field::ExactMetric::diagonal(Vector{weights[d]}).inner(delta, delta);
        if (q > extent) {
            axis = d;
            extent = q;
        }
    }
    const auto axis_metric = rfc::field::ExactMetric::diagonal(Vector{weights[axis]});
    std::vector<std::size_t> order(positions.size()), inverse(positions.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](auto a, auto b) { return positions[a][axis] < positions[b][axis]; });
    for (std::size_t i = 0; i < order.size(); ++i)
        inverse[order[i]] = i;
    std::size_t comparisons = 0;
    auto relations = Json::array(), distances = Json::array();
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type(request.at("relation_type").text());
    rfc::core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < positions.size(); ++i)
        builder.declare_witness_identity(i + 1);
    std::uint64_t occurrence = 0;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const auto anchor = order.begin() + inverse[i];
        const auto within = [&](auto j) {
            const Vector delta{positions[j][axis] - positions[i][axis]};
            return axis_metric.inner(delta, delta) <= radius;
        };
        const auto first =
            std::partition_point(order.begin(), anchor, [&](auto j) { return !within(j); });
        const auto last = std::partition_point(anchor, order.end(), within);
        std::vector<std::size_t> candidates;
        for (auto it = first; it != last; ++it)
            if (*it > i)
                candidates.push_back(*it);
        // Retain the previous source order, occurrence order and relation IDs.
        std::sort(candidates.begin(), candidates.end());
        for (const auto j : candidates) {
            comparisons = rfc::machine::checked_size_add(comparisons, 1);
            const auto delta = rfc::numeric::sub(positions[j], positions[i]);
            const auto distance = metric.inner(delta, delta);
            if (distance > radius)
                continue;
            const auto k = relations.size();
            auto relation = Json::object(), participants = Json::array();
            relation.set("id", Json::string("physical-neighbor:" + std::to_string(k)));
            relation.set("grade", cardinality(1));
            relation.set("type", request.at("relation_type"));
            relation.set("meaning", Json::string("Declared physical neighborhood in " +
                                                 request.at("frame").text() +
                                                 "; coordinate order supplies orientation"));
            std::vector<rfc::core::BoundaryOccurrenceInput> native;
            for (auto index : {i, j}) {
                auto p = Json::object();
                p.set("cell", observations.at(index).at("id"));
                p.set("exists", Json::boolean(true));
                p.set("negative", Json::boolean(index == i));
                p.set("magnitude", Json::string("1"));
                participants.push(p);
                native.push_back({index + 1,
                                  true,
                                  index == i,
                                  Rational{1},
                                  {},
                                  rfc::core::OccurrenceKey{++occurrence}});
            }
            relation.set("participants", participants);
            relations.push(relation);
            builder.add_grade_one(rfc::core::RelationKey{1000000 + k}, type, native);
            auto record = Json::object();
            record.set("relation", relation.at("id"));
            record.set("squared_distance", Json::string(distance.get_str()));
            distances.push(record);
        }
    }
    auto declaration = Json::object();
    declaration.set("relations", relations);
    const auto bytes = rfc::storage::encode_native_series(
        rfc::storage::NativeTemporalSeries(rfc::storage::TemporalEndpoint(builder.freeze()), 1,
                                           source_codec_limits(source_byte_limit())));
    auto out = Json::object();
    out.set("relations", declaration);
    out.set("distances", distances);
    out.set("frame", request.at("frame"));
    out.set("units", request.at("units"));
    out.set("native_complex_sha256",
            Json::string(
                "sha256:" +
                sha256(std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    if (request.has("execution_certificate") && request.at("execution_certificate").boolean()) {
        auto receipt = Json::object();
        receipt.set("method", Json::string("exact_sorted_coordinate_window"));
        receipt.set("coordinate", columns.at(axis));
        receipt.set("native_metric_comparisons", cardinality(comparisons));
        receipt.set("source_observations", cardinality(positions.size()));
        receipt.set("complete", Json::boolean(true));
        out.set("execution", receipt);
    }
    out.set("scope",
            Json::string("Native exact metric neighborhoods over retained physical coordinates. "
                         "This relation describes the declared geometry; it does not infer a "
                         "biochemical interaction or a higher attachment."));
    return out;
}
NativeRecordValue quality_fields(const ResourceImage &resource, const Json &request) {
    using N = NativeRecordValue;
    using rfc::numeric::SparseRationalEntry;
    using rfc::numeric::SparseRationalTensor;
    using rfc::numeric::SparseRationalVector;
    using rfc::numeric::IndexDomain;
    using rfc::numeric::StructuralSupport;
    members(request, {"schema", "operation", "resource", "quality_thresholds"});
    require(resource.declaration.at("format").text() == "fastq",
            "read quality requires four-line FASTQ");
    const auto inspected = inspect_resource(resource.bytes, resource.declaration),
               thresholds = request.at("quality_thresholds");
    require(thresholds.size() > 0 && thresholds.size() <= 94, "quality threshold budget");
    std::vector<std::size_t> cuts;
    std::set<std::size_t> unique;
    for (std::size_t i = 0; i < thresholds.size(); ++i) {
        const auto n = count(thresholds.at(i), 93, "Phred+33 threshold");
        require(unique.insert(n).second, "duplicate quality threshold");
        cuts.push_back(n);
    }
    const auto ls = lines(resource.bytes);
    N::Array ids;
    std::vector<SparseRationalEntry> lengths, lower, upper, sums, below;
    std::vector<std::array<std::size_t, 94>> cycles;
    std::vector<std::size_t> eligible;
    const auto retain = [](auto &entries, std::size_t index, Rational value) {
        if (value != 0) entries.push_back({index, std::move(value)});
    };
    const auto exact_count = [](std::size_t n) { return Rational(std::to_string(n)); };
    std::size_t bases = 0;
    for (std::size_t i = 0; i < inspected.at("records").size(); ++i) {
        const auto sequence = resource.bytes.substr(
            ls.at(4 * i + 1).first, ls.at(4 * i + 1).second - ls.at(4 * i + 1).first);
        const auto qualities = resource.bytes.substr(
            ls.at(4 * i + 3).first, ls.at(4 * i + 3).second - ls.at(4 * i + 3).first);
        bases = rfc::machine::checked_size_add(bases, sequence.size());
        if (cycles.size() < sequence.size()) {
            cycles.resize(sequence.size());
            eligible.resize(sequence.size());
        }
        Vector definite, possible, scores;
        for (std::size_t b = 0; b < sequence.size(); ++b) {
            const auto mask = iupac_mask(sequence[b]),
                       q = static_cast<unsigned char>(qualities[b]) - 33u;
            definite.emplace_back((mask & ~6u) == 0 ? 1 : 0);
            possible.emplace_back((mask & 6u) != 0 ? 1 : 0);
            scores.emplace_back(q);
            cycles[b][q] = rfc::machine::checked_size_add(cycles[b][q], 1);
            eligible[b] = rfc::machine::checked_size_add(eligible[b], 1);
        }
        ids.elements.emplace_back(inspected.at("records").at(i).at("id").text());
        retain(lengths, i, exact_count(sequence.size()));
        retain(lower, i, rfc::numeric::dot(definite, Vector(definite.size(), 1)));
        retain(upper, i, rfc::numeric::dot(possible, Vector(possible.size(), 1)));
        retain(sums, i, rfc::numeric::dot(scores, Vector(scores.size(), 1)));
        for (std::size_t j = 0; j < cuts.size(); ++j) {
            const auto cut = cuts[j];
            Vector indicator;
            for (const auto &q : scores)
                indicator.emplace_back(q < cut ? 1 : 0);
            retain(below, rfc::machine::checked_size_add(rfc::machine::checked_size_mul(i, cuts.size()), j),
                   rfc::numeric::dot(indicator, Vector(indicator.size(), 1)));
        }
    }
    std::vector<SparseRationalEntry> histograms, exposures, thresholds_values;
    std::size_t retained_bases = 0;
    for (std::size_t b = 0; b < cycles.size(); ++b) {
        std::size_t total = 0;
        for (std::size_t q = 0; q < 94; ++q) {
            total = rfc::machine::checked_size_add(total, cycles[b][q]);
            if (cycles[b][q] != 0)
                retain(histograms,
                       rfc::machine::checked_size_add(rfc::machine::checked_size_mul(b, 94), q),
                       exact_count(cycles[b][q]));
        }
        require(total == eligible[b], "quality histogram/exposure conservation");
        retained_bases = rfc::machine::checked_size_add(retained_bases, total);
        retain(exposures, b, exact_count(eligible[b]));
    }
    require(retained_bases == bases, "quality total-base conservation");
    for (std::size_t i = 0; i < cuts.size(); ++i) retain(thresholds_values, i, exact_count(cuts[i]));
    const auto tensor = [](std::vector<std::size_t> shape, std::vector<SparseRationalEntry> values) {
        IndexDomain domain(std::move(shape));
        const auto size = domain.logical_size();
        return N(SparseRationalTensor(std::move(domain), StructuralSupport::full(size),
                                     SparseRationalVector(size, std::move(values))));
    };
    const auto reads = ids.elements.size();
    return N(N::Object{{
        {"schema", N(std::string(quality_fields_schema))}, {"read_ids", N(std::move(ids))},
        {"lengths", tensor({reads}, std::move(lengths))},
        {"gc_lower_bases", tensor({reads}, std::move(lower))},
        {"gc_upper_bases", tensor({reads}, std::move(upper))},
        {"quality_score_sums", tensor({reads}, std::move(sums))},
        {"quality_thresholds", tensor({cuts.size()}, std::move(thresholds_values))},
        {"bases_below", tensor({reads, cuts.size()}, std::move(below))},
        {"cycle_histograms", tensor({cycles.size(), 94}, std::move(histograms))},
        {"eligible_reads", tensor({cycles.size()}, std::move(exposures))},
        {"total_bases", N(rfc::numeric::Integer(std::to_string(bases)))},
        {"scope", N(std::string("Stored Phred+33 scores and exact marginal GC bounds at "
            "ambiguous IUPAC bases. Score histograms describe the input; no "
            "recalibration, trimming or error-rate inference is performed."))}
    }});
}
Json quality_analysis(const ResourceImage &resource, const Json &request) {
    return quality_fields_view(NativeRecord::create(quality_fields(resource, request),
                                                   new_native_record_namespace()));
}
Json quality_fields_view(const NativeRecord &record) {
    require(record.at("schema").text() == quality_fields_schema, "unknown read-quality declaration");
    const auto ids = record.at("read_ids");
    const std::array owners{record.at("lengths"), record.at("gc_lower_bases"),
        record.at("gc_upper_bases"), record.at("quality_score_sums"),
        record.at("quality_thresholds"), record.at("bases_below"),
        record.at("cycle_histograms"), record.at("eligible_reads")};
    const auto &lengths = owners[0].tensor(), &lower = owners[1].tensor(),
               &upper = owners[2].tensor(), &sums = owners[3].tensor(),
               &thresholds = owners[4].tensor(), &below = owners[5].tensor(),
               &histogram = owners[6].tensor(), &eligible = owners[7].tensor();
    const auto shape = [](const auto &field, const std::vector<std::size_t> &expected) {
        require(field.domain().extents() == expected &&
                    field.support().count() == field.domain().logical_size(),
                "read-quality field domain or support");
        for (const auto &entry : field.values().entries())
            require(entry.value >= 0 && entry.value.get_den() == 1,
                    "read-quality count field requires nonnegative integers");
    };
    const auto n = ids.size(), cuts = thresholds.domain().logical_size();
    require(histogram.domain().rank() == 2 && histogram.domain().extents()[1] == 94,
            "read-quality score axis");
    const auto cycles = histogram.domain().extents()[0];
    shape(lengths, {n}); shape(lower, {n}); shape(upper, {n}); shape(sums, {n});
    shape(thresholds, {cuts}); shape(below, {n, cuts});
    shape(histogram, {cycles, 94}); shape(eligible, {cycles});
    auto reads = Json::array(), histograms = Json::array();
    for (std::size_t i = 0; i < n; ++i) {
        auto item = Json::object(), counts = Json::array();
        item.set("id", Json::string(ids.at(i).text()));
        item.set("length", Json::parse(lengths.state(i)->get_str() + " "));
        item.set("gc_lower_bases", Json::string(lower.state(i)->get_str()));
        item.set("gc_upper_bases", Json::string(upper.state(i)->get_str()));
        item.set("quality_score_sum", Json::string(sums.state(i)->get_str()));
        for (std::size_t j = 0; j < cuts; ++j) {
            auto count = Json::object();
            count.set("threshold", Json::parse(thresholds.state(j)->get_str() + " "));
            count.set("bases_below", Json::string(below.state(i * cuts + j)->get_str()));
            counts.push(count);
        }
        item.set("quality_thresholds", counts); reads.push(item);
    }
    // Walk sparse entries once while rendering the complete requested view.
    // Numerical and structural ownership remains in the original native field.
    const auto &entries = histogram.values().entries();
    auto entry = entries.begin();
    for (std::size_t b = 0; b < cycles; ++b) {
        auto item = Json::object(), counts = Json::array();
        item.set("cycle", cardinality(b + 1));
        for (std::size_t q = 0; q < 94; ++q) {
            const auto index = b * 94 + q;
            if (entry != entries.end() && entry->index == index) {
                counts.push(Json::parse(entry->value.get_str() + " ")); ++entry;
            } else counts.push(Json::integer(0));
        }
        item.set("phred33_counts", counts); histograms.push(item);
    }
    auto out = Json::object();
    out.set("reads", reads); out.set("cycle_histograms", histograms);
    out.set("total_bases", Json::parse(record.at("total_bases").integer().get_str() + " "));
    out.set("scope", Json::string(record.at("scope").text()));
    return out;
}
} // namespace gadget::assay_detail
