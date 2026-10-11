#include "genomic_detail.hpp"
#include <gadget/record_view.hpp>

namespace gadget::assay_detail {
namespace {
const rfc::numeric::SparseRationalTensor &
field(const NativeRecord &record, const NativeRecord &columns, std::string_view name) {
    const auto key = columns.at(name).integer();
    for (const auto &tensor : record.endpoint()->tensors)
        if (key == rfc::numeric::Integer(std::to_string(tensor.key)))
            return tensor.tensor;
    throw std::invalid_argument("genomic field is not retained by its publication");
}
std::size_t index(const Rational &value) {
    require(value >= 0 && value.get_den() == 1, "genomic table index must be nonnegative integer");
    return std::stoull(value.get_str());
}
std::string image(const NativeRecord &bytes) {
    const auto span = bytes.bytes();
    return std::string(reinterpret_cast<const char *>(span.data()), span.size());
}
Json integer_array(const rfc::numeric::SparseRationalTensor &values) {
    auto out = Json::array();
    for (std::size_t i = 0; i < values.domain().logical_size(); ++i)
        out.push(cardinality(index(values.state(i).value())));
    return out;
}
} // namespace

Json alignment_fields_view(const NativeRecord &record) {
    auto declaration = Json::object();
    declaration.set("format", Json::string(record.at("format").text()));
    declaration.set("reference", render_record_view(record.at("reference")));
    auto out =
        inspect_hts({declaration, image(record.at("original")),
                     record.has("original_reference") ? image(record.at("original_reference"))
                                                      : std::string{}});
    out.set("htslib_version", Json::string(record.at("htslib_version").text()));
    const auto records = out.at("records");
    const auto columns = record.at("record_fields");
    const auto relations = record.at("relation_fields");
    const auto &selection = field(record, columns, "selected");
    auto selected = Json::array(), excluded = Json::array(), junctions = Json::array();
    for (std::size_t i = 0; i < records.size(); ++i) {
        auto row = records.at(i);
        for (const auto *name : {"flags", "mapping_quality", "query_length"})
            row.set(name, cardinality(index(field(record, columns, name).state(i).value())));
        for (const auto *name : {"start", "end"}) {
            const auto value = field(record, columns, name).state(i);
            row.set(name, cardinality(value ? index(*value) : 0));
        }
        row.set("aligned_blocks", Json::array());
        row.set("splice_junctions", Json::array());
        const auto keep = selection.state(i).value();
        require(keep == 0 || keep == 1, "alignment selection field is not binary");
        (keep == 1 ? selected : excluded).push(cardinality(i));
    }
    const auto &grade = record.endpoint()->source.grade({1});
    const auto &kind = field(record, relations, "kind"),
               &starts = field(record, relations, "start"), &ends = field(record, relations, "end"),
               &counts = field(record, relations, "selected_record_count"),
               &contigs = field(record, relations, "contig"),
               &strands = field(record, relations, "strand");
    for (std::size_t i = 0; i < grade.cell_count(); ++i) {
        const auto type = index(kind.state(i).value());
        require(type <= 1, "alignment relation kind");
        const auto start = index(starts.state(i).value()), end = index(ends.state(i).value());
        auto span = Json::array(), support = Json::array();
        span.push(cardinality(start));
        span.push(cardinality(end));
        const auto [begin, stop] = grade.incidence_range(i);
        for (auto j = begin; j < stop; ++j) {
            const auto lower = grade.lower_index[j];
            records.at(lower).at(type == 0 ? "aligned_blocks" : "splice_junctions").push(span);
            if (selection.state(lower) == Rational{1})
                support.push(cardinality(lower));
        }
        const auto count = counts.state(i).value();
        if (type == 1 && count != 0) {
            auto junction = Json::object();
            junction.set(
                "contig",
                Json::string(record.at("contigs").at(index(contigs.state(i).value())).text()));
            junction.set(
                "alignment_strand",
                Json::string(record.at("strands").at(index(strands.state(i).value())).text()));
            junction.set("intron_start", cardinality(start));
            junction.set("intron_end", cardinality(end));
            junction.set("source_record_indices", support);
            junction.set("alignment_record_count", Json::string(count.get_str()));
            junctions.push(junction);
        }
    }
    out.set("selected_record_indices", selected);
    out.set("excluded_record_indices", excluded);
    out.set("junctions", junctions);
    out.set("junction_quantity", Json::string("Retained CIGAR N occurrences by alignment strand; "
                                              "multiple mappings remain separate records"));
    return out;
}

Json feature_count_fields_view(const NativeRecord &record) {
    const auto columns = record.at("segments"), groups = record.at("read_groups");
    const auto names = record.at("decoded_query_names").bytes();
    const auto number = [&](const char *name, std::size_t i) {
        return index(columns.at(name).tensor().state(i).value());
    };
    auto segments = Json::array();
    for (std::size_t i = 0; i < std::stoull(record.at("segment_count").integer().get_str()); ++i) {
        auto row = Json::object(), originals = Json::array();
        const auto offset = number("query_name_offset", i), length = number("query_name_length", i);
        require(offset <= names.size() && length <= names.size() - offset,
                "segment query byte span");
        row.set("read_group", Json::string(groups.at(number("read_group", i)).text()));
        row.set("query", Json::string(std::string(
                             reinterpret_cast<const char *>(names.data() + offset), length)));
        row.set("mate", cardinality(number("mate", i)));
        row.set("class_index", cardinality(number("class_index", i)));
        const auto begin = number("source_record_offset", i),
                   count = number("source_record_count", i);
        require(begin <= columns.at("source_records").tensor().domain().logical_size() &&
                    count <= columns.at("source_records").tensor().domain().logical_size() - begin,
                "segment source record range");
        for (auto j = begin; j < begin + count; ++j)
            originals.push(cardinality(number("source_records", j)));
        row.set("source_record_indices", originals);
        segments.push(row);
    }
    const auto candidates = record.at("candidates"), overlaps = record.at("overlaps");
    auto witnesses = Json::array();
    for (std::size_t i = 0; i < std::stoull(record.at("overlap_count").integer().get_str()); ++i) {
        auto row = Json::object();
        row.set("source_record_index",
                cardinality(index(overlaps.at("source_record_index").tensor().state(i).value())));
        row.set("feature",
                Json::string(
                    candidates.at(index(overlaps.at("feature").tensor().state(i).value())).text()));
        for (const auto *name : {"overlap_bases", "aligned_bases"})
            row.set(name, Json::string(overlaps.at(name).tensor().state(i)->get_str()));
        row.set("eligible",
                Json::boolean(overlaps.at("eligible").tensor().state(i) == Rational{1}));
        if (overlaps.has("junction_compatible"))
            row.set(
                "junction_compatible",
                Json::boolean(overlaps.at("junction_compatible").tensor().state(i) == Rational{1}));
        witnesses.push(row);
    }
    const auto fields = record.at("fields");
    const auto &counts = field(record, fields, "class_count"),
               &lower = field(record, fields, "lower"), &upper = field(record, fields, "upper");
    const auto features = render_record_view(record.at("features")),
               request = render_record_view(record.at("request"));
    auto allocation = Json::object(), targets = Json::array(), sample_counts = Json::array(),
         bounds = Json::array();
    std::map<std::string, std::size_t> ids;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        ids.emplace(candidates.at(i).text(), i);
    for (std::size_t i = 0; i < features.size(); ++i) {
        const auto id = features.at(i).at("id");
        auto target = Json::object(), member = Json::array(), bound = Json::object();
        member.push(id);
        target.set("id", id);
        target.set("candidates", member);
        targets.push(target);
        const auto coordinate = ids.at(id.text());
        const auto lo = lower.state(coordinate).value(), hi = upper.state(coordinate).value();
        bound.set("id", id);
        bound.set("lower", Json::string(lo.get_str()));
        bound.set("upper", Json::string(hi.get_str()));
        bound.set("identifiable", Json::boolean(lo == hi));
        bounds.push(bound);
    }
    for (std::size_t i = 0; i < counts.domain().logical_size(); ++i)
        sample_counts.push(Json::string(counts.state(i)->get_str()));
    auto sample = Json::object(), samples = Json::array(), evidence = Json::array();
    sample.set("id", request.at("resource"));
    sample.set("counts", sample_counts);
    sample.set("scale", Json::string("1"));
    sample.set("scale_provenance",
               Json::string("Declared identity scale for retained sequenced-segment counts"));
    samples.push(sample);
    evidence.push(Json::string("sha256:" + sha256(image(record.at("original_alignments")))));
    evidence.push(Json::string("sha256:" + sha256(image(record.at("original_annotation")))));
    allocation.set("schema", Json::string("gadget.assay.allocations.v1"));
    allocation.set("operation", Json::string("allocation_bounds"));
    allocation.set("candidates", render_record_view(candidates));
    allocation.set("classes", render_record_view(record.at("classes")));
    allocation.set("samples", samples);
    allocation.set("targets", targets);
    allocation.set("contrasts", Json::array());
    allocation.set("evidence", evidence);
    allocation.set("units", Json::string("assigned_segment_count"));
    allocation.set("class_partition", Json::string("disjoint_observations"));
    allocation.set("candidate_search", request.at("candidate_search"));
    auto analysis = Json::object(), analyses = Json::array();
    analysis.set("id", request.at("resource"));
    analysis.set("total", Json::string(record.at("allocation_total").rational().get_str()));
    analysis.set("bounds", bounds);
    analyses.push(analysis);
    auto result = Json::object();
    result.set("samples", analyses);
    result.set("contrasts", Json::array());
    const auto summary = record.has("allocation_fields")
        ? allocation_fields_summary_view(record.at("allocation_fields"))
        : render_record_view(record.at("allocation_summary"));
    json_object_object_foreach(summary.raw(), name, value) {
        (void)value;
        result.set(name, summary.at(name));
    }
    const auto identity = record.at("allocation_identity");
    auto certificate = Json::object();
    certificate.set("schema", Json::string(identity.at("schema").text()));
    certificate.set("query", allocation);
    certificate.set("source_digest", Json::string(identity.at("source_digest").text()));
    certificate.set("result", result);
    for (const auto *name : {"arithmetic", "rexgraph_sha256", "gadget_sha256"})
        certificate.set(name, Json::string(identity.at(name).text()));
    certificate = seal(certificate, "certificate_sha256");
    auto out = Json::object();
    out.set("features", features);
    out.set("segments", segments);
    out.set("excluded_record_indices",
            integer_array(record.at("excluded_record_indices").tensor()));
    out.set("overlap_witnesses", witnesses);
    out.set("allocation_request", allocation);
    out.set("reference", render_record_view(record.at("reference")));
    out.set("unit", Json::string(record.at("unit").text()));
    out.set("scope", Json::string(record.at("scope").text()));
    if (request.has("junction_policy")) {
        out.set("junction_policy", request.at("junction_policy"));
        out.set("execution", render_record_view(record.at("execution")));
    }
    out.set("allocation_certificate", certificate);
    return out;
}
} // namespace gadget::assay_detail
