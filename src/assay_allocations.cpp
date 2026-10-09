#include "assay_detail.hpp"
#include <gadget/resources.hpp>
#include <limits>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget {
Json analyze_allocations(const Json &request, std::shared_ptr<CpuResources> resources) {
    using namespace assay_detail;
    using rfc::numeric::ExactMatrix;
    if (!resources)
        resources = CpuResources::local();
    auto reservation = resources->reserve(1);
    members(request,
            {"schema", "operation", "candidates", "classes", "samples", "targets", "contrasts",
             "class_partition", "candidate_search", "evidence", "units"},
            {"kernel_certificate"});
    const bool retain_factors =
        request.has("kernel_certificate") && request.at("kernel_certificate").boolean();
    require(request.at("schema").text() == "gadget.assay.allocations.v1" &&
                request.at("operation").text() == "allocation_bounds",
            "allocation schema/operation");
    require(request.at("class_partition").text() == "disjoint_observations",
            "allocation constraints require disjoint observation classes");
    require(request.at("units").text() == "assigned_fragment_count" ||
                request.at("units").text() == "assigned_segment_count" ||
                request.at("units").text() == "assigned_observation_count",
            "allocation measurement units");
    const auto search = request.at("candidate_search");
    members(search, {"method", "completeness", "reference_sha256"});
    for (const auto *m : {"method", "completeness", "reference_sha256"})
        (void)text(search.at(m), source_byte_limit());
    require(request.at("evidence").size() > 0, "allocation source evidence");
    for (std::size_t i = 0; i < request.at("evidence").size(); ++i)
        (void)text(request.at("evidence").at(i), source_byte_limit());
    const auto candidates = request.at("candidates"), classes = request.at("classes"),
               samples = request.at("samples"), targets = request.at("targets");
    require(candidates.size() > 0 && classes.size() > 0 && samples.size() > 0 && targets.size() > 0,
            "allocation shape");
    require(candidates.size() <= std::numeric_limits<std::uint32_t>::max() &&
                classes.size() <= std::numeric_limits<std::uint32_t>::max(),
            "allocation native cell coordinate representation");
    std::map<std::string, std::size_t> ids;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        require(ids.emplace(text(candidates.at(i), source_byte_limit()), i).second,
                "duplicate allocation candidate");
    std::vector<std::vector<std::size_t>> participants;
    std::size_t slots = 0;
    std::set<std::string> class_ids, target_ids;
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("reported.assignment.candidate.class");
    rfc::core::RelationalComplexBuilder builder(signature);
    for (std::size_t i = 0; i < candidates.size(); ++i)
        builder.declare_witness_identity(i + 1);
    for (std::size_t i = 0; i < classes.size(); ++i) {
        const auto c = classes.at(i);
        members(c, {"id", "candidates"});
        require(class_ids.insert(text(c.at("id"), source_byte_limit())).second,
                "duplicate allocation class");
        const auto p = c.at("candidates");
        require(p.size() > 0, "allocation class arity");
        std::set<std::string> seen;
        std::vector<std::size_t> indices;
        std::vector<rfc::core::BoundaryOccurrenceInput> occurrences;
        for (std::size_t j = 0; j < p.size(); ++j) {
            const auto id = text(p.at(j), source_byte_limit());
            require(ids.contains(id) && seen.insert(id).second,
                    "unknown/duplicate class candidate");
            indices.push_back(ids.at(id));
            slots = rfc::machine::checked_size_add(slots, 1);
            occurrences.push_back(
                {ids.at(id) + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{slots}});
        }
        participants.push_back(indices);
        builder.add_grade_one(rfc::core::RelationKey{rfc::machine::checked_size_add(1000000, i)},
                              type, occurrences);
    }
    const auto complex = builder.freeze();
    ExactMatrix lower(targets.size(), classes.size()), upper(targets.size(), classes.size());
    std::vector<std::vector<std::size_t>> candidate_targets(candidates.size());
    for (std::size_t t = 0; t < targets.size(); ++t) {
        const auto target = targets.at(t);
        members(target, {"id", "candidates"});
        require(target_ids.insert(text(target.at("id"), source_byte_limit())).second,
                "duplicate allocation target");
        std::set<std::size_t> set;
        require(target.at("candidates").size() > 0, "empty allocation target");
        for (std::size_t j = 0; j < target.at("candidates").size(); ++j) {
            const auto id = text(target.at("candidates").at(j), source_byte_limit());
            require(ids.contains(id) && set.insert(ids.at(id)).second, "target candidate identity");
        }
        for (auto p : set)
            candidate_targets[p].push_back(t);
    }
    // Construct incidence maps from set membership, then delegate every count
    // action to RFC. A class wholly inside a target contributes to both bounds;
    // a mixed class contributes only to its upper bound. These endpoints are
    // attained by assigning that class to an inside/outside candidate respectively.
    for (std::size_t c = 0; c < classes.size(); ++c) {
        std::map<std::size_t, std::size_t> hits;
        for (auto p : participants[c])
            for (auto t : candidate_targets[p])
                ++hits[t];
        for (const auto &[t, n] : hits) {
            upper.set(t, c, 1);
            if (n == participants[c].size()) {
                lower.set(t, c, 1);
            }
        }
    }
    std::map<std::string, std::size_t> sample_ids;
    std::vector<Vector> lows, highs;
    std::vector<Rational> scales;
    auto analyses = Json::array();
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto sample = samples.at(s);
        members(sample, {"id", "counts", "scale", "scale_provenance"});
        require(sample_ids.emplace(text(sample.at("id"), source_byte_limit()), s).second,
                "duplicate allocation sample");
        (void)text(sample.at("scale_provenance"), source_byte_limit());
        const auto scale = exact_rational(text(sample.at("scale"), source_byte_limit()));
        require(scale > 0, "allocation sample scale");
        scales.push_back(scale);
        require(sample.at("counts").size() == classes.size(), "class counts/coordinate shape");
        Vector counts;
        for (std::size_t c = 0; c < classes.size(); ++c) {
            const auto q = exact_rational(sample.at("counts").at(c).text());
            require(q >= 0 && q.get_den() == 1, "class count must be nonnegative integer");
            counts.push_back(q);
        }
        const auto lo = lower.apply(counts), hi = upper.apply(counts);
        lows.push_back(lo);
        highs.push_back(hi);
        auto record = Json::object(), bounds = Json::array();
        record.set("id", sample.at("id"));
        record.set("total",
                   Json::string(rfc::numeric::dot(Vector(counts.size(), 1), counts).get_str()));
        for (std::size_t t = 0; t < targets.size(); ++t) {
            auto b = Json::object();
            b.set("id", targets.at(t).at("id"));
            b.set("lower", Json::string(lo[t].get_str()));
            b.set("upper", Json::string(hi[t].get_str()));
            b.set("identifiable", Json::boolean(lo[t] == hi[t]));
            bounds.push(b);
        }
        record.set("bounds", bounds);
        analyses.push(record);
    }
    auto contrast_results = Json::array();
    for (std::size_t i = 0; i < request.at("contrasts").size(); ++i) {
        const auto c = request.at("contrasts").at(i);
        members(c, {"id", "left", "right"});
        const auto l = text(c.at("left"), source_byte_limit()),
                   r = text(c.at("right"), source_byte_limit());
        require(l != r && sample_ids.contains(l) && sample_ids.contains(r),
                "allocation contrast sample identity");
        const auto li = sample_ids.at(l), ri = sample_ids.at(r);
        const auto lo = rfc::numeric::sub(rfc::numeric::scale(scales[li], lows[li]),
                                          rfc::numeric::scale(scales[ri], highs[ri]));
        const auto hi = rfc::numeric::sub(rfc::numeric::scale(scales[li], highs[li]),
                                          rfc::numeric::scale(scales[ri], lows[ri]));
        auto record = Json::object(), bounds = Json::array();
        record.set("id", c.at("id"));
        record.set("left", c.at("left"));
        record.set("right", c.at("right"));
        for (std::size_t t = 0; t < targets.size(); ++t) {
            auto b = Json::object();
            b.set("id", targets.at(t).at("id"));
            b.set("lower", Json::string(lo[t].get_str()));
            b.set("upper", Json::string(hi[t].get_str()));
            b.set("direction",
                  Json::string(lo[t] > 0 ? "positive_for_every_allocation"
                                         : (hi[t] < 0 ? "negative_for_every_allocation"
                                                      : (lo[t] == 0 && hi[t] == 0
                                                             ? "zero_for_every_allocation"
                                                             : "unresolved"))));
            bounds.push(b);
        }
        record.set("bounds", bounds);
        contrast_results.push(record);
    }
    auto result = Json::object();
    result.set("samples", analyses);
    result.set("contrasts", contrast_results);
    result.set("bound_scope",
               Json::string("Sharp marginal bounds over independent nonnegative class allocations; "
                            "endpoints for different targets need not be jointly attainable"));
    result.set("native_relation_count", cardinality(complex.grade({1}).cell_count()));
    result.set("native_occurrence_count", cardinality(slots));
    rfc::storage::TemporalEndpoint endpoint(complex);
    const auto archive = rfc::storage::encode_native_series(
        rfc::storage::NativeTemporalSeries(endpoint, 1, source_codec_limits(source_byte_limit())));
    result.set(
        "native_complex_sha256",
        Json::string("sha256:" + sha256(std::string(reinterpret_cast<const char *>(archive.data()),
                                                    archive.size()))));
    // Disjoint observation classes have independent occurrence coordinates.
    // Compute their exact kernels natively, then retain common target rows when
    // measuring the image. Image ranks cannot be added across classes: target
    // directions from different classes may coincide.
    ExactMatrix target_image(targets.size(), slots - classes.size());
    std::size_t slot_offset = 0, kernel_offset = 0, constraint_rank = 0;
    auto factors = Json::array();
    const auto entries = [](const ExactMatrix &matrix) {
        auto out = Json::array();
        for (const auto &[coordinate, q] : matrix.entries()) {
            auto row = Json::array();
            row.push(cardinality(coordinate.first));
            row.push(cardinality(coordinate.second));
            row.push(Json::string(q.get_str()));
            out.push(row);
        }
        return out;
    };
    for (std::size_t c = 0; c < classes.size(); ++c) {
        const auto &p = participants[c];
        ExactMatrix retained(1, p.size()), aggregate(targets.size(), p.size());
        for (std::size_t j = 0; j < p.size(); ++j) {
            retained.set(0, j, 1);
            for (auto t : candidate_targets[p[j]])
                aggregate.set(t, j, 1);
        }
        const auto kernel = retained.kernel_frame();
        require(retained.multiply(kernel).nonzeros() == 0, "native allocation kernel residual");
        const auto image = aggregate.multiply(kernel);
        for (const auto &[coordinate, q] : image.entries())
            target_image.set(coordinate.first,
                             rfc::machine::checked_size_add(kernel_offset, coordinate.second), q);
        constraint_rank = rfc::machine::checked_size_add(constraint_rank, retained.rank());
        if (retain_factors) {
            auto factor = Json::object();
            factor.set("class_id", classes.at(c).at("id"));
            factor.set("native_relation_identity",
                       cardinality(rfc::machine::checked_size_add(1000000, c)));
            factor.set("candidates", classes.at(c).at("candidates"));
            factor.set("occurrence_offset", cardinality(slot_offset));
            factor.set("kernel_offset", cardinality(kernel_offset));
            factor.set("kernel_dimension", cardinality(kernel.cols()));
            factor.set("kernel_frame", entries(kernel));
            factor.set("target_image", entries(image));
            factors.push(factor);
        }
        slot_offset = rfc::machine::checked_size_add(slot_offset, p.size());
        kernel_offset = rfc::machine::checked_size_add(kernel_offset, kernel.cols());
    }
    require(slot_offset == slots && kernel_offset == target_image.cols() &&
                constraint_rank + kernel_offset == slots,
            "allocation factor coordinate shape");
    result.set("constraint_rank", cardinality(constraint_rank));
    result.set("allocation_nullity", cardinality(kernel_offset));
    result.set("target_ambiguity_rank", cardinality(target_image.rank()));
    result.set("kernel_residual_zero", Json::boolean(true));
    if (retain_factors) {
        auto certificate = Json::object();
        certificate.set("factors", factors);
        certificate.set("target_coordinate_order", targets);
        certificate.set(
            "coordinate_rule",
            Json::string(
                "Class order then candidate order; occurrence key is occurrence_offset+j+1; "
                "target rows follow target_coordinate_order; global kernel column is "
                "kernel_offset+k"));
        certificate.set("bound_map_nonzeros", cardinality(rfc::machine::checked_size_add(
                                                  lower.nonzeros(), upper.nonzeros())));
        certificate.set("target_image_nonzeros", cardinality(target_image.nonzeros()));
        result.set("kernel_certificate", certificate);
    }
    auto output = Json::object();
    output.set("schema", Json::string("gadget.assay.result.v1"));
    output.set("query", request);
    output.set("source_digest", Json::string("sha256:" + sha256(request.dump())));
    output.set("result", result);
    output.set("arithmetic", Json::string("rational"));
    output.set("rexgraph_sha256", Json::string(GADGET_REXGRAPH_SHA256));
    output.set("gadget_sha256", Json::string(GADGET_SOURCE_SHA256));
    return seal(output, "certificate_sha256");
}
} // namespace gadget
