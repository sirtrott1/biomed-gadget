#include "interval_detail.hpp"
#include "genomic_detail.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/numeric/exact_matrix.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::assay_detail {
namespace interval_detail {
Overlay::Overlay(const std::vector<Record> &records, const std::string &policy) {
    require(policy == "ignore" || policy == "separate", "interval strand policy");
    struct Event {
        std::vector<std::size_t> enter, leave;
    };
    std::map<Frame, std::map<std::uint64_t, Event>> frames;
    for (std::size_t i = 0; i < records.size(); ++i) {
        const auto &r = records[i].interval;
        auto &frame = frames[{r.contig, policy == "ignore" ? "ignored" : r.strand}];
        frame[r.start].enter.push_back(i);
        frame[r.end].leave.push_back(i);
    }
    for (const auto &[frame, events] : frames) {
        std::set<std::size_t> active;
        for (auto it = events.begin(); it != events.end(); ++it) {
            for (auto i : it->second.leave)
                require(active.erase(i) == 1, "interval leave event");
            for (auto i : it->second.enter)
                require(active.insert(i).second, "interval enter event");
            const auto next = std::next(it);
            if (next == events.end() || active.empty())
                continue;
            memberships = rfc::machine::checked_size_add(memberships, active.size());
            atoms.push_back({frame, it->first, next->first, {active.begin(), active.end()}});
        }
        require(active.empty(), "complete interval event closure");
    }
}
rfc::numeric::ExactMatrix Overlay::coverage(std::size_t records) const {
    rfc::numeric::ExactMatrix map(atoms.size(), records);
    for (std::size_t i = 0; i < atoms.size(); ++i)
        for (auto j : atoms[i].records)
            map.set(i, j, 1);
    return map;
}
Vector Overlay::lengths() const {
    Vector out;
    for (const auto &a : atoms)
        out.push_back(exact_rational(std::to_string(a.end - a.start)));
    return out;
}
Json threshold(const std::vector<Atom> &atoms, const Vector &lengths, const Vector &values,
               const Rational &cut, bool report_atoms, NativeRecordValue *native) {
    require(values.size() == atoms.size() && lengths.size() == atoms.size(),
            "threshold carrier shape");
    auto item = Json::object(), regions = Json::array();
    using N = NativeRecordValue;
    N::Array native_regions;
    Vector included(atoms.size(), 0);
    std::size_t retained = 0;
    for (std::size_t i = 0; i < atoms.size();) {
        if (values[i] < cut) {
            ++i;
            continue;
        }
        const auto begin = i;
        const auto frame = atoms[i].frame;
        auto end = atoms[i].end;
        included[i] = lengths[i];
        ++retained;
        ++i;
        while (i < atoms.size() && atoms[i].frame == frame && atoms[i].start == end &&
               values[i] >= cut) {
            end = atoms[i].end;
            included[i] = lengths[i];
            ++retained;
            ++i;
        }
        if (native)
            native_regions.elements.emplace_back(N::Object{{
                {"contig", N(frame.first)}, {"strand", N(frame.second)},
                {"start", N(rfc::numeric::Integer(std::to_string(atoms[begin].start)))},
                {"end", N(rfc::numeric::Integer(std::to_string(end)))}}});
        else {
            auto region = Json::object();
            region.set("contig", Json::string(frame.first));
            region.set("strand", Json::string(frame.second));
            region.set("start", cardinality(atoms[begin].start));
            region.set("end", cardinality(end));
            regions.push(region);
        }
    }
    const auto length = rfc::numeric::dot(included, Vector(included.size(), 1));
    if (native) {
        N::Object value{{{"threshold", N(cut)}, {"regions", N(std::move(native_regions))}}};
        if (report_atoms)
            value.members.emplace_back("retained_atoms", N(rfc::numeric::Integer(std::to_string(retained))));
        value.members.emplace_back(report_atoms ? "retained_length" : "length", N(length));
        *native = N(std::move(value));
        return Json{};
    }
    item.set("threshold", Json::string(cut.get_str()));
    item.set("regions", regions);
    if (report_atoms)
        item.set("retained_atoms", cardinality(retained));
    item.set(report_atoms ? "retained_length" : "length",
             Json::string(length.get_str()));
    return item;
}
} // namespace interval_detail
AssayFieldsSource interval_fields_source(const ResourceImage &resource, const Json &request) {
    using N = NativeRecordValue;
    members(request, {"schema", "operation", "resource", "thresholds", "strand_policy"});
    const auto records = read_intervals(resource.bytes, resource.declaration);
    const auto policy = request.at("strand_policy").text();
    require(policy == "separate" || policy == "ignore", "coverage strand policy");
    for (const auto &r : records)
        require(
            r.weight >= 0,
            "coverage weights must be nonnegative; signed observations need a different quantity");
    std::vector<interval_detail::Record> source_records;
    for (std::size_t i = 0; i < records.size(); ++i)
        source_records.push_back({records[i], 0, i});
    const interval_detail::Overlay overlay(source_records, policy);
    const auto &atoms = overlay.atoms;
    const auto memberships = overlay.memberships;
    const auto coverage = overlay.coverage(records.size());
    Vector weights, lengths;
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("physical.interval.coverage.membership");
    rfc::core::RelationalComplexBuilder builder(signature);
    builder.declare_grade(rfc::core::Grade{1});
    for (std::size_t i = 0; i < records.size(); ++i) {
        builder.declare_witness_identity(i + 1);
        weights.push_back(records[i].weight);
    }
    std::uint64_t occurrence = 1;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        std::vector<rfc::core::BoundaryOccurrenceInput> p;
        for (auto index : atoms[i].records) {
            p.push_back(
                {index + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{occurrence++}});
        }
        builder.add_grade_one(rfc::core::RelationKey{1000000 + i}, type, p);
        lengths.push_back(exact_rational(std::to_string(atoms[i].end - atoms[i].start)));
    }
    const auto source = builder.freeze();
    const auto values = coverage.apply(weights);
    N::Array components;
    const auto thresholds = request.at("thresholds");
    require(thresholds.size() > 0, "coverage thresholds must be nonempty");
    std::set<Rational> unique;
    for (std::size_t t = 0; t < thresholds.size(); ++t) {
        const auto threshold = exact_rational(text(thresholds.at(t), source_byte_limit()));
        require(threshold >= 0 && unique.insert(threshold).second,
                "coverage threshold nonnegative/unique");
        N item;
        interval_detail::threshold(atoms, lengths, values, threshold, true, &item);
        components.elements.push_back(std::move(item));
    }
    rfc::storage::TemporalEndpoint endpoint(source);
    const auto archive = rfc::storage::encode_native_series(
        rfc::storage::NativeTemporalSeries(endpoint, 1, source_codec_limits(source_byte_limit())));
    auto result = Json::object();
    result.set("reference", resource.declaration.at("reference"));
    result.set("strand_policy", request.at("strand_policy"));
    result.set("integral_units", Json::string("declared interval weight times reference bases"));
    result.set("domain", Json::string("Union of retained interval observations; regions outside "
                                      "this support remain unobserved"));
    result.set(
        "native_complex_sha256",
        Json::string("sha256:" + sha256(std::string(reinterpret_cast<const char *>(archive.data()),
                                                    archive.size()))));
    result.set("native_relation_count", cardinality(atoms.size()));
    result.set("native_occurrence_count", cardinality(memberships));
    result.set(
        "construction",
        Json::string("Exact physical intersection segments; every interval-to-segment occurrence "
                     "retained as a grade-one relation. Threshold regions are physical connected "
                     "components, not significance calls or inferred biochemical cycles."));
    endpoint.carriers = {
        {1, rfc::storage::TemporalCarrierKind::SourceGrade, {0}, {}, 0,
         rfc::numeric::IndexDomain::linear(records.size())},
        {2, rfc::storage::TemporalCarrierKind::SourceGrade, {1}, {}, 0,
         rfc::numeric::IndexDomain::linear(atoms.size())}};
    SourceGradeFields fields(endpoint);
    N::Array frames, original;
    std::map<interval_detail::Frame, std::size_t> frame_indices;
    for (const auto &atom : atoms)
        if (!frame_indices.contains(atom.frame)) {
            frame_indices.emplace(atom.frame, frames.elements.size());
            frames.elements.emplace_back(N::Object{{{"contig", N(atom.frame.first)},
                                                    {"strand", N(atom.frame.second)}}});
        }
    for (const auto &record : records)
        original.elements.emplace_back(N::Object{{
            {"id", N(record.id)}, {"contig", N(record.contig)}, {"strand", N(record.strand)},
            {"bytes", N(N::ByteSpan{1, record.source_start,
                                    record.source_end - record.source_start, false})}}});
    N declaration(N::Object{{
        {"schema", N(std::string("gadget.intervals.coverage.fields.v1"))},
        {"summary", import_record_view(result)}, {"frames", N(std::move(frames))},
        {"original_records", N(std::move(original))},
        {"source_bytes", N(N::ByteSpan{1, 0, resource.bytes.size(), false})},
        {"record_start", fields.column("record.start", 0, records.size(), [&](auto i) { return Rational(records[i].start); })},
        {"record_end", fields.column("record.end", 0, records.size(), [&](auto i) { return Rational(records[i].end); })},
        {"weights", fields.field("record.weight", 0, weights)},
        {"coverage_map", fields.map("interval.coverage", 1, 2, coverage)},
        {"lengths", fields.field("segment.length", 1, lengths)},
        {"coverage", fields.field("segment.coverage", 1, values)},
        {"frame", fields.column("segment.frame", 1, atoms.size(), [&](auto i) { return Rational(frame_indices.at(atoms[i].frame)); })},
        {"start", fields.column("segment.start", 1, atoms.size(), [&](auto i) { return Rational(atoms[i].start); })},
        {"end", fields.column("segment.end", 1, atoms.size(), [&](auto i) { return Rational(atoms[i].end); })},
        {"threshold_components", N(std::move(components))},
        {"coverage_integral", N(rfc::numeric::dot(lengths, values))}}});
    return {std::move(endpoint), std::move(declaration),
            {RawByteSource(1, std::as_bytes(std::span(resource.bytes)))}};
}
Json interval_fields_view(const NativeRecord &record) {
    auto out = render_record_view(record.at("summary")), segments = Json::array();
    const auto &source = record.endpoint()->source.grade(rfc::core::Grade{1});
    const auto &frame = retained_tensor(record, record.at("frame")).tensor;
    const auto &start = retained_tensor(record, record.at("start")).tensor;
    const auto &end = retained_tensor(record, record.at("end")).tensor;
    const auto &values = retained_tensor(record, record.at("coverage")).tensor;
    for (std::size_t i = 0; i < source.cell_count(); ++i) {
        auto atom = render_record_view(record.at("frames").at(frame.state(i).value().get_num().get_ui()));
        atom.set("start", cardinality(start.state(i).value().get_num().get_ui()));
        atom.set("end", cardinality(end.state(i).value().get_num().get_ui()));
        atom.set("coverage", Json::string(values.state(i).value().get_str()));
        auto indices = Json::array();
        const auto range = source.incidence_range(i);
        for (auto j = range.first; j < range.second; ++j)
            indices.push(cardinality(source.lower_index.at(j)));
        atom.set("source_record_indices", indices);
        segments.push(atom);
    }
    out.set("segments", segments);
    out.set("threshold_components", render_record_view(record.at("threshold_components")));
    out.set("coverage_integral", Json::string(record.at("coverage_integral").rational().get_str()));
    return out;
}
} // namespace gadget::assay_detail
