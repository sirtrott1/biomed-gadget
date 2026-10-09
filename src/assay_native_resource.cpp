#include "assay_detail.hpp"
#include "source_grade_fields.hpp"
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>

namespace gadget {
NativeRecord inspect_bioinformatics_record(const fs::path &specification) {
    using namespace assay_detail;
    using namespace rfc::storage;
    using N = NativeRecordValue;
    const auto path = fs::canonical(specification);
    const auto request = read_record(path).native();
    const auto spec = render_record_view(request);
    members(spec, {"schema", "source", "format"},
            {"reference", "weight_column", "reference_source"});
    require(spec.at("schema").text() == "gadget.bioinformatics.file.v1", "physical file schema");
    const auto format = spec.at("format").text();
    require(format == "bed" || format == "gtf" || format == "gff3" || format == "fasta" ||
                format == "fastq",
            "native physical records support BED, GTF, GFF3, FASTA and FASTQ");
    const auto bytes = local_artifact(path.parent_path(), spec.at("source"));
    const bool sequence = format == "fasta" || format == "fastq";
    const auto intervals = sequence ? std::vector<Interval>{} : read_intervals(bytes, spec);
    const auto sequences =
        sequence ? read_sequence_spans(bytes, format) : std::vector<SequenceSpan>{};
    const auto count = sequence ? sequences.size() : intervals.size();
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("raw.source.record");
    auto source = [&] {
        rfc::core::RelationalComplexBuilder builder(signature);
        builder.declare_witness_identity(1);
        builder.declare_grade({1});
        for (std::size_t i = 0; i < count; ++i)
            builder.add_grade_one(
                rfc::core::RelationKey{i + 1}, type,
                {{1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{i + 1}}});
        return builder.freeze();
    }();
    TemporalEndpoint endpoint(std::move(source));
    endpoint.carriers.push_back(
        {1, TemporalCarrierKind::SourceGrade, {0}, {}, 0, rfc::numeric::IndexDomain::linear(1)});
    endpoint.carriers.push_back({2,
                                 TemporalCarrierKind::SourceGrade,
                                 {1},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(count)});
    SourceGradeFields fields{endpoint};
    N::Object columns;
    N::Object dictionaries;
    const auto store = [&](const std::string &name, auto value) {
        columns.members.emplace_back(name, fields.column(name, 1, count, value));
    };
    store("source_start", [&](auto i) {
        return Rational(std::to_string(sequence ? sequences[i].start : intervals[i].source_start));
    });
    store("source_end", [&](auto i) {
        return Rational(std::to_string(sequence ? sequences[i].end : intervals[i].source_end));
    });
    if (sequence) {
        store("length", [&](auto i) { return Rational(std::to_string(sequences[i].length)); });
        store("id_start", [&](auto i) { return Rational(std::to_string(sequences[i].id_start)); });
        store("id_length",
              [&](auto i) { return Rational(std::to_string(sequences[i].id_length)); });
    } else {
        store("start", [&](auto i) { return Rational(std::to_string(intervals[i].start)); });
        store("end", [&](auto i) { return Rational(std::to_string(intervals[i].end)); });
        store("weight", [&](auto i) { return intervals[i].weight; });
        const auto dictionary = [&](const std::string &name, auto value) {
            std::map<std::string, std::size_t> indices;
            N::Array labels;
            store(name, [&](auto i) {
                const auto label = value(i);
                const auto [it, added] = indices.emplace(label, indices.size());
                if (added)
                    labels.elements.emplace_back(label);
                return Rational(std::to_string(it->second));
            });
            dictionaries.members.emplace_back(name, N(std::move(labels)));
        };
        dictionary("contig", [&](auto i) { return intervals[i].contig; });
        dictionary("strand", [&](auto i) { return intervals[i].strand; });
        dictionary("feature_type", [&](auto i) { return intervals[i].feature_type; });
    }
    std::vector<RawByteSource> originals;
    originals.emplace_back(1, std::as_bytes(std::span(bytes.data(), bytes.size())));
    N::Object declaration{
        {{"schema", N(std::string("gadget.physical.records.v1"))},
         {"format", N(format)},
         {"records", N(rfc::numeric::Integer(std::to_string(count)))},
         {"columns", N(std::move(columns))},
         {"dictionaries", N(std::move(dictionaries))},
         {"original",
          N(N::Owned{0,
                     StableCellKey{{0}, 1},
                     {},
                     std::make_shared<const N>(N(N::ByteSpan{1, 0, bytes.size(), false}))})},
         {"source_semantics",
          N(std::string("Original record to original byte source containment"))}}};
    if (spec.has("reference"))
        declaration.members.emplace_back(
            "reference", import_record_view(reference_record(spec.at("reference"))));
    if (spec.has("reference_source")) {
        const auto reference = local_artifact(path.parent_path(), spec.at("reference_source"));
        originals.emplace_back(2, std::as_bytes(std::span(reference.data(), reference.size())));
        declaration.members.emplace_back("reference_original",
                                         N(N::ByteSpan{2, 0, reference.size(), false}));
    }
    return retain_operation(std::move(endpoint), new_native_record_namespace(),
                            std::move(originals), request, N(std::move(declaration)));
}
} // namespace gadget
