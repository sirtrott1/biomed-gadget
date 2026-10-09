#include "interval_detail.hpp"
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
               const Rational &cut, bool report_atoms) {
    require(values.size() == atoms.size() && lengths.size() == atoms.size(),
            "threshold carrier shape");
    auto item = Json::object(), regions = Json::array();
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
        auto region = Json::object();
        region.set("contig", Json::string(frame.first));
        region.set("strand", Json::string(frame.second));
        region.set("start", cardinality(atoms[begin].start));
        region.set("end", cardinality(end));
        regions.push(region);
    }
    item.set("threshold", Json::string(cut.get_str()));
    item.set("regions", regions);
    if (report_atoms)
        item.set("retained_atoms", cardinality(retained));
    item.set(report_atoms ? "retained_length" : "length",
             Json::string(rfc::numeric::dot(included, Vector(included.size(), 1)).get_str()));
    return item;
}
} // namespace interval_detail
Json interval_analysis(const std::string &bytes, const Json &resource, const Json &request) {
    const auto records = read_intervals(bytes, resource);
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
    auto components = Json::array(), segments = Json::array();
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto &a = atoms[i];
        auto atom = Json::object(), support = Json::array();
        atom.set("contig", Json::string(a.frame.first));
        atom.set("strand", Json::string(a.frame.second));
        atom.set("start", cardinality(a.start));
        atom.set("end", cardinality(a.end));
        atom.set("coverage", Json::string(values[i].get_str()));
        for (auto j : a.records)
            support.push(cardinality(j));
        atom.set("source_record_indices", support);
        segments.push(atom);
    }
    const auto thresholds = request.at("thresholds");
    require(thresholds.size() > 0, "coverage thresholds must be nonempty");
    std::set<Rational> unique;
    for (std::size_t t = 0; t < thresholds.size(); ++t) {
        const auto threshold = exact_rational(text(thresholds.at(t), source_byte_limit()));
        require(threshold >= 0 && unique.insert(threshold).second,
                "coverage threshold nonnegative/unique");
        auto item = interval_detail::threshold(atoms, lengths, values, threshold, true);
        components.push(item);
    }
    rfc::storage::TemporalEndpoint endpoint(source);
    const auto archive = rfc::storage::encode_native_series(
        rfc::storage::NativeTemporalSeries(endpoint, 1, source_codec_limits(source_byte_limit())));
    auto result = Json::object();
    result.set("reference", resource.at("reference"));
    result.set("strand_policy", request.at("strand_policy"));
    result.set("segments", segments);
    result.set("threshold_components", components);
    result.set("coverage_integral", Json::string(rfc::numeric::dot(lengths, values).get_str()));
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
    return result;
}
} // namespace gadget::assay_detail
