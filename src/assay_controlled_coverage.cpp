#include "assay_coverage_detail.hpp"
#include "interval_detail.hpp"
#include "scientific_projection.hpp"
#include <algorithm>
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/field/boundary.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget::assay_detail {
namespace {
using N = NativeRecordValue;
using rfc::numeric::ExactMatrix;
using rfc::numeric::SparseRationalTensor;
N integer(std::uint64_t value) { return N(rfc::numeric::Integer(std::to_string(value))); }
N vector_value(const Vector &values) {
    std::vector<rfc::numeric::SparseRationalEntry> entries;
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] != 0)
            entries.push_back({i, values[i]});
    return N(SparseRationalTensor(
        rfc::numeric::IndexDomain::linear(values.size()),
        rfc::numeric::StructuralSupport::full(values.size()),
        rfc::numeric::SparseRationalVector(values.size(), std::move(entries))));
}
N matrix_value(const ExactMatrix &matrix) {
    const auto domain = rfc::numeric::IndexDomain({matrix.rows(), matrix.cols()});
    std::vector<rfc::numeric::SparseRationalEntry> entries;
    for (const auto &[coordinate, value] : matrix.entries())
        entries.push_back({coordinate.first * matrix.cols() + coordinate.second, value});
    return N(SparseRationalTensor(
        domain, rfc::numeric::StructuralSupport::full(domain.logical_size()),
        rfc::numeric::SparseRationalVector(domain.logical_size(), std::move(entries))));
}
bool zero(const Vector &values) {
    return std::all_of(values.begin(), values.end(), [](const auto &x) { return x == 0; });
}
struct Fields {
    rfc::storage::TemporalEndpoint &endpoint;
    template <class Value>
    N column(std::string name, std::uint16_t grade, std::size_t size, Value value) {
        std::vector<rfc::numeric::SparseRationalEntry> entries;
        for (std::size_t i = 0; i < size; ++i) {
            Rational x = value(i);
            if (x != 0)
                entries.push_back({i, std::move(x)});
        }
        rfc::storage::TemporalTensor tensor;
        tensor.key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
        tensor.role = rfc::storage::TemporalTensorRole::Field;
        tensor.output = grade + 1;
        tensor.interpretation = std::move(name);
        tensor.tensor = SparseRationalTensor(
            rfc::numeric::IndexDomain::linear(size), rfc::numeric::StructuralSupport::full(size),
            rfc::numeric::SparseRationalVector(size, std::move(entries)));
        const auto key = tensor.key;
        endpoint.tensors.push_back(std::move(tensor));
        return integer(key);
    }
    N field(std::string name, std::uint16_t grade, const Vector &values) {
        return column(std::move(name), grade, values.size(), [&](auto i) { return values[i]; });
    }
};
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
std::uint64_t field_key(NativeRecord key) {
    const auto value = key.integer();
    require(value > 0 && value <= rfc::numeric::Integer(
                                      std::to_string(std::numeric_limits<std::uint64_t>::max())),
            "coverage field identity is outside its declared native domain");
    return std::stoull(value.get_str());
}
const SparseRationalTensor &field(const NativeRecord &record, std::uint64_t identity) {
    for (const auto &tensor : record.endpoint()->tensors)
        if (identity == tensor.key)
            return tensor.tensor;
    throw std::invalid_argument("coverage field is not retained by its owning endpoint");
}
std::size_t index(const Rational &value) {
    require(value >= 0 && value.get_den() == 1 &&
                value <=
                    rfc::numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())),
            "coverage coordinate is outside the host address domain");
    return std::stoull(value.get_num().get_str());
}
Json values(const SparseRationalTensor &tensor) {
    auto out = Json::array();
    for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i)
        out.push(Json::string(tensor.state(i).value().get_str()));
    return out;
}
Json sparse_matrix(const NativeRecord &record) {
    const auto &tensor = record.tensor();
    const auto &shape = tensor.domain().extents();
    require(shape.size() == 2, "coverage matrix declaration shape");
    auto out = Json::object(), entries = Json::array();
    for (const auto &entry : tensor.values().entries()) {
        auto row = Json::object();
        row.set("row", cardinality(entry.index / shape[1]));
        row.set("column", cardinality(entry.index % shape[1]));
        row.set("value", Json::string(entry.value.get_str()));
        entries.push(row);
    }
    out.set("rows", cardinality(shape[0]));
    out.set("columns", cardinality(shape[1]));
    out.set("entries", entries);
    return out;
}
} // namespace

CoverageSource controlled_coverage_source(const ResourceImages &images, const Json &request) {
    members(request,
            {"schema", "operation", "samples", "contrasts", "controls", "intercept", "metric",
             "strand_policy", "thresholds", "absence_policy", "units", "evidence"},
            {"source_kernel"});
    require(request.at("metric").text() == "reference_length",
            "coverage metric must be reference_length");
    require(request.at("absence_policy").text() == "measured_zero_within_union",
            "coverage requires explicitly declared measured zeros within the observed union");
    (void)text(request.at("units"), source_byte_limit());
    (void)text(request.at("evidence"), source_byte_limit());
    const auto samples = request.at("samples");
    require(samples.size() > 0, "coverage samples must be nonempty");
    std::map<std::string, std::size_t> sample_ids;
    std::map<std::string, std::uint64_t> resource_keys;
    std::vector<RawByteSource> raw;
    N::Array resources;
    std::vector<interval_detail::Record> records;
    Vector scales;
    Json reference;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto sample = samples.at(s);
        members(sample, {"id", "resource", "scale", "scale_provenance"});
        require(sample_ids.emplace(text(sample.at("id"), source_byte_limit()), s).second,
                "duplicate coverage sample");
        const auto scale = exact_rational(text(sample.at("scale"), source_byte_limit()));
        require(scale > 0, "coverage scale divisor must be positive");
        (void)text(sample.at("scale_provenance"), source_byte_limit());
        scales.push_back(scale);
        const auto name = text(sample.at("resource"), source_byte_limit());
        require(images.contains(name), "unknown coverage resource");
        const auto &image = images.at(name);
        const auto frame = reference_record(image.declaration.at("reference"));
        if (s == 0)
            reference = frame;
        else
            require(reference.dump() == frame.dump(),
                    "coverage reference assemblies, sequences and contigs must agree");
        if (!resource_keys.contains(name)) {
            const auto key = rfc::machine::checked_size_add(raw.size(), 1);
            resource_keys.emplace(name, key);
            raw.emplace_back(key, std::as_bytes(std::span(image.bytes.data(), image.bytes.size())));
            resources.elements.emplace_back(
                N::Object{{{"id", N(name)},
                           {"declaration", import_record_view(image.declaration)},
                           {"bytes", N(N::ByteSpan{key, 0, image.bytes.size(), false})}}});
        }
        auto intervals = read_intervals(image.bytes, image.declaration);
        for (std::size_t j = 0; j < intervals.size(); ++j) {
            require(intervals[j].weight >= 0, "coverage weights must be nonnegative");
            records.push_back({std::move(intervals[j]), s, j});
        }
    }
    interval_detail::Overlay overlay(records, request.at("strand_policy").text());
    require(!overlay.atoms.empty(), "coverage requires an observed physical union");
    const auto lengths = overlay.lengths();
    const auto metric = rfc::field::ExactMetric::diagonal(lengths);
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("physical.interval.segment.membership");
    auto source = [&] {
        rfc::core::RelationalComplexBuilder builder(signature);
        builder.declare_grade({1});
        std::vector<std::vector<std::size_t>> memberships(records.size());
        std::vector<std::size_t> counts(records.size());
        for (const auto &atom : overlay.atoms)
            for (const auto j : atom.records)
                ++counts[j];
        for (std::size_t j = 0; j < records.size(); ++j)
            memberships[j].reserve(counts[j]);
        for (std::size_t i = 0; i < overlay.atoms.size(); ++i) {
            builder.declare_witness_identity(i + 1);
            for (const auto j : overlay.atoms[i].records)
                memberships[j].push_back(i);
            std::vector<std::size_t>().swap(overlay.atoms[i].records);
        }
        std::uint64_t occurrence = 0;
        for (std::size_t j = 0; j < records.size(); ++j) {
            std::vector<rfc::core::BoundaryOccurrenceInput> native;
            native.reserve(memberships[j].size());
            for (const auto i : memberships[j]) {
                require(occurrence != std::numeric_limits<std::uint64_t>::max(),
                        "interval occurrence identity space exhausted");
                native.push_back(
                    {i + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{++occurrence}});
            }
            builder.add_grade_one(rfc::core::RelationKey{j + 1}, type, native);
            std::vector<std::size_t>().swap(memberships[j]);
        }
        return builder.freeze();
    }();
    rfc::storage::TemporalEndpoint endpoint(std::move(source));
    std::vector<std::vector<rfc::numeric::SparseRationalEntry>> sample_weights(samples.size());
    ExactMatrix original_weights(records.size(), samples.size());
    for (std::size_t j = 0; j < records.size(); ++j) {
        sample_weights[records[j].sample].push_back({j, records[j].interval.weight});
        original_weights.set(j, records[j].sample, records[j].interval.weight);
    }
    endpoint.carriers.push_back({1,
                                 rfc::storage::TemporalCarrierKind::SourceGrade,
                                 {0},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(overlay.atoms.size())});
    endpoint.carriers.push_back({2,
                                 rfc::storage::TemporalCarrierKind::SourceGrade,
                                 {1},
                                 {},
                                 0,
                                 rfc::numeric::IndexDomain::linear(records.size())});
    Fields store{endpoint};
    const rfc::field::BoundaryOperator boundary(endpoint.source, {1});
    Dictionary contigs, strands;
    for (const auto &record : records) {
        contigs.add(record.interval.contig);
        strands.add(record.interval.strand);
    }
    for (const auto &atom : overlay.atoms) {
        contigs.add(atom.frame.first);
        strands.add(atom.frame.second);
    }
    N::Object interval_columns, segment_columns;
    const auto intervals_column = [&](const std::string &name, auto value) {
        interval_columns.members.emplace_back(
            name, store.column("interval." + name, 1, records.size(), value));
    };
    intervals_column("sample", [&](auto i) { return Rational(std::to_string(records[i].sample)); });
    intervals_column("source_record_index",
                     [&](auto i) { return Rational(std::to_string(records[i].index)); });
    intervals_column("contig", [&](auto i) {
        return Rational(std::to_string(contigs.indices.at(records[i].interval.contig)));
    });
    intervals_column("strand", [&](auto i) {
        return Rational(std::to_string(strands.indices.at(records[i].interval.strand)));
    });
    intervals_column("start",
                     [&](auto i) { return Rational(std::to_string(records[i].interval.start)); });
    intervals_column("end",
                     [&](auto i) { return Rational(std::to_string(records[i].interval.end)); });
    intervals_column("weight", [&](auto i) { return records[i].interval.weight; });
    intervals_column("source_start", [&](auto i) {
        return Rational(std::to_string(records[i].interval.source_start));
    });
    intervals_column("source_end", [&](auto i) {
        return Rational(std::to_string(records[i].interval.source_end));
    });
    std::string names;
    std::vector<std::size_t> name_offsets;
    for (const auto &record : records) {
        name_offsets.push_back(names.size());
        names += record.interval.id;
    }
    intervals_column("id_offset",
                     [&](auto i) { return Rational(std::to_string(name_offsets[i])); });
    intervals_column("id_length", [&](auto i) {
        return Rational(std::to_string(records[i].interval.id.size()));
    });
    const auto names_key = rfc::machine::checked_size_add(raw.size(), 1);
    raw.emplace_back(names_key, std::as_bytes(std::span(names.data(), names.size())));
    for (const auto &name : {"contig", "strand", "start", "end"}) {
        segment_columns.members.emplace_back(
            name,
            store.column("segment." + std::string(name), 0, overlay.atoms.size(), [&](auto i) {
                const auto &atom = overlay.atoms[i];
                if (std::string_view(name) == "contig")
                    return Rational(std::to_string(contigs.indices.at(atom.frame.first)));
                if (std::string_view(name) == "strand")
                    return Rational(std::to_string(strands.indices.at(atom.frame.second)));
                return Rational(
                    std::to_string(std::string_view(name) == "start" ? atom.start : atom.end));
            }));
    }
    const auto metric_key = store.field("reference_length", 0, lengths);
    std::vector<Vector> fields;
    N::Array sample_fields;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const rfc::numeric::SparseRationalVector weights(records.size(),
                                                         std::move(sample_weights[s]));
        fields.push_back(
            rfc::numeric::scale(Rational{1} / scales[s], boundary.apply(weights).dense()));
        sample_fields.elements.emplace_back(
            N::Object{{{"sample", N(samples.at(s).at("id").text())},
                       {"scale_divisor", N(scales[s])},
                       {"scale_provenance", N(samples.at(s).at("scale_provenance").text())},
                       {"field", store.field("sample:" + std::to_string(s), 0, fields.back())},
                       {"integral", N(rfc::numeric::dot(lengths, fields.back()))}}});
    }
    const auto controls = request.at("controls");
    const auto intercept = request.at("intercept").boolean();
    ExactMatrix declared(overlay.atoms.size(),
                         rfc::machine::checked_size_add(controls.size(), intercept ? 1 : 0));
    std::set<std::string> seen;
    for (std::size_t c = 0; c < controls.size(); ++c) {
        const auto id = text(controls.at(c), source_byte_limit());
        require(sample_ids.contains(id) && seen.insert(id).second,
                "unknown or duplicate coverage control");
        const auto &v = fields.at(sample_ids.at(id));
        for (std::size_t i = 0; i < v.size(); ++i)
            declared.set(i, c, v[i]);
    }
    if (intercept)
        for (std::size_t i = 0; i < declared.rows(); ++i)
            declared.set(i, controls.size(), 1);
    const PreparedProjection projection(declared.image_frame(), metric);
    std::set<Rational> thresholds;
    N::Array cuts;
    for (std::size_t i = 0; i < request.at("thresholds").size(); ++i) {
        const auto cut = exact_rational(text(request.at("thresholds").at(i), source_byte_limit()));
        require(thresholds.insert(cut).second, "duplicate coverage threshold");
        cuts.elements.emplace_back(cut);
    }
    const auto contrasts = request.at("contrasts");
    require(contrasts.size() > 0, "coverage contrasts must be nonempty");
    N::Array results;
    ExactMatrix residuals(overlay.atoms.size(), contrasts.size());
    seen.clear();
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto contrast = contrasts.at(c);
        members(contrast, {"id", "terms", "evidence"});
        require(seen.insert(text(contrast.at("id"), source_byte_limit())).second,
                "duplicate coverage contrast");
        (void)text(contrast.at("evidence"), source_byte_limit());
        const auto terms = contrast.at("terms");
        require(terms.size() > 0, "coverage contrast terms must be nonempty");
        Vector input(overlay.atoms.size(), 0), source_coefficients(samples.size(), 0);
        std::set<std::string> term_ids;
        for (std::size_t t = 0; t < terms.size(); ++t) {
            const auto term = terms.at(t);
            members(term, {"sample", "coefficient"});
            const auto id = text(term.at("sample"), source_byte_limit());
            require(sample_ids.contains(id) && term_ids.insert(id).second,
                    "unknown or duplicate contrast sample");
            const auto coefficient =
                exact_rational(text(term.at("coefficient"), source_byte_limit()));
            source_coefficients[sample_ids.at(id)] = coefficient / scales.at(sample_ids.at(id));
            input = rfc::numeric::add(
                input, rfc::numeric::scale(coefficient, fields.at(sample_ids.at(id))));
        }
        const auto source_contrast = original_weights.apply(source_coefficients);
        require(boundary.apply(source_contrast) == input,
                "source contrast and physical coverage action must commute");
        const auto coordinates = projection.coefficients(input),
                   represented = projection.frame.apply(coordinates);
        const auto residual = rfc::numeric::sub(input, represented);
        require(rfc::numeric::add(represented, residual) == input &&
                    zero(projection.dual.apply(metric.apply(residual))) &&
                    metric.inner(represented, residual) == 0,
                "native controlled coverage residual certificate");
        for (std::size_t i = 0; i < residual.size(); ++i)
            residuals.set(i, c, residual[i]);
        results.elements.emplace_back(N::Object{
            {{"id", N(contrast.at("id").text())},
             {"terms", import_record_view(terms)},
             {"evidence", N(contrast.at("evidence").text())},
             {"input", store.field("contrast:" + std::to_string(c) + ":input", 0, input)},
             {"represented",
              store.field("contrast:" + std::to_string(c) + ":represented", 0, represented)},
             {"residual", store.field("contrast:" + std::to_string(c) + ":residual", 0, residual)},
             {"source_coefficients",
              store.field("interval-contrast:" + std::to_string(c), 1, source_contrast)},
             {"frame_coordinates", vector_value(coordinates)},
             {"quadrances", N(N::Object{{{"input", N(metric.inner(input, input))},
                                         {"represented", N(metric.inner(represented, represented))},
                                         {"residual", N(metric.inner(residual, residual))}}})},
             {"residual_integral", N(rfc::numeric::dot(lengths, residual))}}});
    }
    N::Object root{
        {{"schema", N(std::string(coverage_fields_schema))},
         {"reference", import_record_view(reference)},
         {"request", import_record_view(request)},
         {"resources", N(std::move(resources))},
         {"contigs", N(std::move(contigs.labels))},
         {"strands", N(std::move(strands.labels))},
         {"intervals", N(std::move(interval_columns))},
         {"segments", N(std::move(segment_columns))},
         {"source_ids", N(N::ByteSpan{names_key, 0, names.size(), false})},
         {"metric_diagonal", metric_key},
         {"sample_fields", N(std::move(sample_fields))},
         {"contrasts", N(std::move(results))},
         {"thresholds", N(std::move(cuts))},
         {"controls", N(N::Object{{{"samples", import_record_view(controls)},
                                   {"intercept", N(intercept)},
                                   {"declared_frame", matrix_value(declared)},
                                   {"independent_frame", matrix_value(projection.frame)},
                                   {"rank", integer(projection.frame.cols())},
                                   {"gram", matrix_value(projection.solve.matrix())},
                                   {"factorizations", integer(1)}}})},
         {"residual_gram", matrix_value(rfc::field::metric_gram(residuals, metric, residuals))}}};
    if (request.has("source_kernel") && request.at("source_kernel").boolean()) {
        const auto matrix = rfc::field::materialize_boundary_matrix(endpoint.source, {1});
        const auto frame = matrix.kernel_frame();
        require(matrix.multiply(frame).nonzeros() == 0, "native interval source kernel residual");
        root.members.emplace_back(
            "source_kernel",
            N(N::Object{{{"frame", matrix_value(frame)},
                         {"rank", integer(records.size() - frame.cols())},
                         {"nullity", integer(frame.cols())},
                         {"boundary_residual_zero", N(true)},
                         {"meaning", N(std::string("Exact signed combinations of original interval "
                                                   "variables with zero physical coverage"))}}}));
    }
    return {std::move(endpoint), N(std::move(root)), std::move(raw)};
}

NativeRecord controlled_coverage_record(const ResourceImages &images, const Json &request) {
    auto source = controlled_coverage_source(images, request);
    return NativeRecord::create_from_endpoint(source.declaration, new_native_record_namespace(),
                                             std::move(source.sources), std::move(source.endpoint));
}

Json controlled_coverage(const ResourceImages &images, const Json &request) {
    return coverage_fields_view(controlled_coverage_record(images, request));
}

Json coverage_fields_view(const NativeRecord &record) {
    require(record.at("schema").text() == coverage_fields_schema, "unknown coverage declaration");
    const auto request = render_record_view(record.at("request"));
    const auto intervals = record.at("intervals"), segments = record.at("segments");
    const auto column = [&](const NativeRecord &columns,
                            std::string_view name) -> const SparseRationalTensor & {
        return field(record, field_key(columns.at(name)));
    };
    const auto &sample_index = column(intervals, "sample"), &starts = column(intervals, "start"),
               &ends = column(intervals, "end"), &weights = column(intervals, "weight"),
               &record_index = column(intervals, "source_record_index"),
               &contig_index = column(intervals, "contig"),
               &strand_index = column(intervals, "strand"),
               &source_start = column(intervals, "source_start"),
               &source_end = column(intervals, "source_end"),
               &id_offset = column(intervals, "id_offset"),
               &id_length = column(intervals, "id_length");
    const auto &segment_start = column(segments, "start"), &segment_end = column(segments, "end"),
               &segment_contig = column(segments, "contig"),
               &segment_strand = column(segments, "strand");
    const auto ids = record.at("source_ids").bytes();
    std::map<std::string, std::string> raw;
    const auto resources = record.at("resources");
    for (std::size_t i = 0; i < resources.size(); ++i) {
        const auto bytes = resources.at(i).at("bytes").bytes();
        raw.emplace(resources.at(i).at("id").text(),
                    std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
    }
    auto original = Json::array();
    for (std::size_t j = 0; j < sample_index.domain().logical_size(); ++j) {
        const auto sample = request.at("samples").at(index(sample_index.state(j).value()));
        const auto name = sample.at("resource").text();
        const auto offset = index(id_offset.state(j).value()),
                   length = index(id_length.state(j).value());
        require(offset <= ids.size() && length <= ids.size() - offset,
                "coverage source name byte span");
        auto row = Json::object();
        row.set("id", Json::string("interval:" + std::to_string(j)));
        row.set("sample", sample.at("id"));
        row.set("resource", sample.at("resource"));
        row.set("source_record_index", cardinality(index(record_index.state(j).value())));
        row.set("source_id", Json::string(std::string(
                                 reinterpret_cast<const char *>(ids.data() + offset), length)));
        row.set("contig",
                Json::string(record.at("contigs").at(index(contig_index.state(j).value())).text()));
        row.set("strand",
                Json::string(record.at("strands").at(index(strand_index.state(j).value())).text()));
        row.set("start", cardinality(index(starts.state(j).value())));
        row.set("end", cardinality(index(ends.state(j).value())));
        row.set("weight", Json::string(weights.state(j).value().get_str()));
        row.set("source", span(name, raw.at(name), index(source_start.state(j).value()),
                               index(source_end.state(j).value())));
        original.push(row);
    }
    std::vector<interval_detail::Atom> atoms;
    std::vector<Json> support;
    support.reserve(segment_start.domain().logical_size());
    for (std::size_t i = 0; i < segment_start.domain().logical_size(); ++i)
        support.push_back(Json::array());
    auto relations = Json::array();
    const auto &grade = record.endpoint()->source.grade({1});
    for (std::size_t j = 0; j < grade.cell_count(); ++j) {
        auto row = Json::object(), participants = Json::array();
        row.set("id", Json::string("interval:" + std::to_string(j)));
        row.set("grade", cardinality(1));
        row.set("type", Json::string("physical.interval.segment.membership"));
        row.set(
            "meaning",
            Json::string("Every physical segment covered by this retained interval observation"));
        const auto [begin, end] = grade.incidence_range(j);
        for (auto slot = begin; slot < end; ++slot) {
            auto p = Json::object();
            const auto i = grade.lower_index[slot];
            p.set("cell", Json::string("segment:" + std::to_string(i)));
            p.set("exists", Json::boolean(grade.aspects.get(
                                record.endpoint()->source.signature().existence_aspect(), slot)));
            p.set("negative",
                  Json::boolean(grade.aspects.get(
                      record.endpoint()->source.signature().orientation_aspect(), slot)));
            auto magnitude = record.endpoint()->source.incidence_coefficient({1}, slot);
            if (magnitude < 0)
                magnitude = -magnitude;
            p.set("magnitude", Json::string(magnitude.get_str()));
            participants.push(p);
            support.at(i).push(cardinality(j));
        }
        row.set("participants", participants);
        relations.push(row);
    }
    auto segment_rows = Json::array();
    for (std::size_t i = 0; i < segment_start.domain().logical_size(); ++i) {
        const auto contig = record.at("contigs").at(index(segment_contig.state(i).value())).text();
        const auto strand = record.at("strands").at(index(segment_strand.state(i).value())).text();
        const auto start = index(segment_start.state(i).value()),
                   end = index(segment_end.state(i).value());
        atoms.push_back({{contig, strand}, start, end, {}});
        auto row = Json::object();
        row.set("id", Json::string("segment:" + std::to_string(i)));
        row.set("contig", Json::string(contig));
        row.set("strand", Json::string(strand));
        row.set("start", cardinality(start));
        row.set("end", cardinality(end));
        row.set("source_record_indices", support[i]);
        segment_rows.push(row);
    }
    auto sample_fields = Json::array();
    const auto samples = record.at("sample_fields");
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto sample = samples.at(s);
        auto row = Json::object();
        row.set("sample", Json::string(sample.at("sample").text()));
        row.set("scale_divisor", request.at("samples").at(s).at("scale"));
        row.set("scale_provenance", Json::string(sample.at("scale_provenance").text()));
        row.set("values", values(field(record, field_key(sample.at("field")))));
        row.set("integral", Json::string(sample.at("integral").rational().get_str()));
        sample_fields.push(row);
    }
    auto results = Json::array(), relation_fields = Json::array();
    const auto contrasts = record.at("contrasts");
    const auto lengths = field(record, field_key(record.at("metric_diagonal"))).values().dense();
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto contrast = contrasts.at(c);
        auto row = Json::object();
        row.set("id", Json::string(contrast.at("id").text()));
        row.set("terms", render_record_view(contrast.at("terms")));
        row.set("evidence", Json::string(contrast.at("evidence").text()));
        for (const auto *name : {"input", "represented", "residual", "source_coefficients"})
            row.set(name, values(field(record, field_key(contrast.at(name)))));
        row.set("frame_coordinates", values(contrast.at("frame_coordinates").tensor()));
        row.set("source_action_exact", Json::boolean(true));
        row.set("quadrances", render_record_view(contrast.at("quadrances")));
        row.set("residual_integral",
                Json::string(contrast.at("residual_integral").rational().get_str()));
        row.set("reconstruction_exact", Json::boolean(true));
        row.set("control_orthogonal", Json::boolean(true));
        auto components = Json::array();
        const auto residual = field(record, field_key(contrast.at("residual"))).values().dense();
        for (std::size_t t = 0; t < record.at("thresholds").size(); ++t)
            components.push(interval_detail::threshold(
                atoms, lengths, residual, record.at("thresholds").at(t).rational(), true));
        row.set("residual_threshold_components", components);
        results.push(row);
        auto source_field = Json::object(), entries = Json::array();
        source_field.set("id", Json::string("interval-contrast:" + std::to_string(c)));
        source_field.set("grade", cardinality(1));
        source_field.set("channel", Json::string("observation"));
        source_field.set("quantity", Json::string("signed_interval_contrast_coefficient"));
        source_field.set("units", request.at("units"));
        source_field.set("method",
                         Json::string("Declared contrast coefficients on original interval "
                                      "observations, divided by declared sample scales"));
        const auto source_key = contrast.at("source_coefficients");
        const auto &tensor = field(record, field_key(source_key));
        for (std::size_t j = 0; j < tensor.domain().logical_size(); ++j) {
            auto entry = Json::object();
            entry.set("relation", Json::string("interval:" + std::to_string(j)));
            entry.set("value", Json::string(tensor.state(j).value().get_str()));
            entries.push(entry);
        }
        source_field.set("values", entries);
        relation_fields.push(source_field);
    }
    auto controls = Json::object();
    const auto control = record.at("controls");
    for (const auto *name : {"samples", "intercept", "rank", "factorizations"})
        controls.set(name, render_record_view(control.at(name)));
    for (const auto *name : {"declared_frame", "independent_frame", "gram"})
        controls.set(name, sparse_matrix(control.at(name)));
    auto out = Json::object(), declaration = Json::object();
    declaration.set("relations", relations);
    declaration.set("fields", relation_fields);
    out.set("segments", segment_rows);
    out.set("source_records", original);
    out.set("sample_fields", sample_fields);
    out.set("contrasts", results);
    out.set("controls", controls);
    out.set("metric_diagonal", values(field(record, field_key(record.at("metric_diagonal")))));
    out.set("residual_gram", sparse_matrix(record.at("residual_gram")));
    out.set("relation_declaration", declaration);
    out.set("reference", render_record_view(record.at("reference")));
    const auto bytes = rfc::storage::encode_native_series(rfc::storage::NativeTemporalSeries(
        rfc::storage::TemporalEndpoint(record.endpoint()->source), 1,
        source_codec_limits(source_byte_limit())));
    out.set("native_complex_sha256",
            Json::string(
                "sha256:" +
                sha256(std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()))));
    const auto occurrences = grade.incidence_count();
    out.set("native_occurrence_count", cardinality(occurrences));
    if (record.has("source_kernel")) {
        const auto source_kernel = record.at("source_kernel");
        auto kernel = Json::object();
        kernel.set("frame", sparse_matrix(source_kernel.at("frame")));
        for (const auto *name : {"rank", "nullity", "boundary_residual_zero", "meaning"})
            kernel.set(name, render_record_view(source_kernel.at(name)));
        out.set("source_kernel", kernel);
    }
    out.set("quantity", Json::string("controlled_coverage_field"));
    out.set("units", request.at("units"));
    out.set("evidence", request.at("evidence"));
    out.set(
        "scope",
        Json::string(
            "Control projection on the complete physical interval union. Reference lengths weight "
            "the shared metric; source intervals remain distinct grade-one cells. Zero-boundary "
            "combinations describe coverage redundancy, not biochemical feedback."));
    return out;
}
} // namespace gadget::assay_detail
