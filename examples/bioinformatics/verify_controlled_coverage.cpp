#include <algorithm>
#include <charconv>
#include <gadget/native_record.hpp>
#include <iostream>
#include <map>
#include <rfc/storage/mapped_file.hpp>
#include <string_view>

using gadget::NativeRecord;
using rfc::numeric::Rational;
using Vector = std::vector<Rational>;
using Matrix = std::vector<Vector>;
namespace {
void check(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
std::uint64_t integer(std::string_view text) {
    std::uint64_t value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    check(!text.empty() && error == std::errc{} && end == text.data() + text.size(),
          "invalid unsigned source coordinate");
    return value;
}
std::uint64_t integer(const Rational &value) {
    check(value >= 0 && value.get_den() == 1, "nonintegral source coordinate");
    return integer(value.get_num().get_str());
}
struct Column {
    const rfc::numeric::SparseRationalTensor *tensor;
    Rational at(std::size_t i) const {
        const auto value = tensor->state(i);
        check(value.has_value(), "missing measured coordinate");
        return *value;
    }
};
Column field(const NativeRecord &record, const NativeRecord &key) {
    const auto identity = integer(key.integer().get_str());
    for (const auto &item : record.endpoint()->tensors)
        if (item.key == identity) return {&item.tensor};
    throw std::runtime_error("missing native field");
}
Rational entry(const NativeRecord &matrix, std::size_t row, std::size_t column) {
    const auto &tensor = matrix.tensor();
    const auto &shape = tensor.domain().extents();
    check(shape.size() == 2 && row < shape[0] && column < shape[1], "matrix coordinate");
    return tensor.state(row * shape[1] + column).value();
}
std::size_t rank(Matrix matrix) {
    if (matrix.empty()) return 0;
    const auto columns = matrix.front().size();
    std::size_t pivot = 0;
    for (std::size_t column = 0; column < columns && pivot < matrix.size(); ++column) {
        auto next = pivot;
        while (next < matrix.size() && matrix[next][column] == 0) ++next;
        if (next == matrix.size()) continue;
        std::swap(matrix[pivot], matrix[next]);
        const Rational scale = matrix[pivot][column];
        for (auto &value : matrix[pivot]) value /= scale;
        for (std::size_t row = pivot + 1; row < matrix.size(); ++row) {
            const Rational coefficient = matrix[row][column];
            for (std::size_t c = column; c < columns; ++c)
                matrix[row][c] -= coefficient * matrix[pivot][c];
        }
        ++pivot;
    }
    return pivot;
}
std::vector<std::string_view> split(std::string_view row) {
    std::vector<std::string_view> result;
    for (std::size_t p = 0;;) {
        auto q = row.find('\t', p);
        if (q == std::string_view::npos) q = row.size();
        result.push_back(row.substr(p, q - p));
        if (q == row.size()) return result;
        p = q + 1;
    }
}
using Frame = std::pair<std::string, std::string>;
struct Interval {
    Frame frame;
    std::uint64_t start, end;
};
struct Event {
    std::uint64_t position;
    std::size_t sample;
    Rational weight;
    int count;
};
void same_bytes(std::span<const std::byte> a, std::span<const std::byte> b) {
    check(a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin()),
          "retained original bytes differ");
}
NativeRecord select(const NativeRecord &archive, const std::string &assay) {
    auto owner = archive;
    if (!assay.empty()) {
        std::optional<NativeRecord> selected;
        const auto publications = archive.at("publications");
        for (std::size_t i = 0; i < publications.size(); ++i) {
            const auto item = publications.at(i);
            if (item.at("assay").text() == assay) {
                check(!selected, "select a standalone publication for partitioned assays");
                selected = NativeRecord::decode(item.at("bytes").bytes());
            }
        }
        check(selected.has_value(), "assay publication is absent");
        owner = *selected;
    }
    return owner.has("history")
        ? owner.historical_publication(integer(owner.at("history").at("cutoff").integer().get_str()))
        : owner;
}
} // namespace
int verify(int argc, char **argv) {
    std::string assay;
    if (argc >= 3 && std::string_view(argv[1]) == "--assay") {
        assay = argv[2]; argc -= 2; argv += 2;
    }
    check(argc >= 5, "usage: gadget_verify_controlled_coverage [--assay ID] RECORD WORKSPACE REQUEST BED... (sample order)");
    const auto archive = NativeRecord::read(argv[1]);
    const auto operation = select(archive, assay);
    const auto record = operation.at("result").at("result");
    check(record.at("schema").text() == "gadget.coverage.fields.v1", "unsupported coverage source");
    const rfc::storage::ReadOnlyMappedFile workspace(argv[2]), request_bytes(argv[3]);
    same_bytes(operation.at("result").at("workspace_publication").bytes(), workspace.bytes());
    same_bytes(operation.at("request_publication").bytes(), request_bytes.bytes());
    const auto request = record.at("request"), samples = request.at("samples");
    check(samples.size() == static_cast<std::size_t>(argc - 4), "one original BED is required per declared sample");
    check(request.at("metric").text() == "reference_length" &&
              request.at("absence_policy").text() == "measured_zero_within_union",
          "unsupported physical metric or absence policy");
    const auto strand_policy = request.at("strand_policy").text();
    check(strand_policy == "ignore" || strand_policy == "separate", "strand policy");
    const auto source = record.endpoint();
    check(source->source.grade_count() == 2, "coverage requires its retained two grade source");
    const auto &relations = source->source.grade({1});
    const auto atoms = source->source.grade({0}).cell_count();
    const auto intervals = record.at("intervals"), segments = record.at("segments");
    const auto starts = field(record, intervals.at("start")), ends = field(record, intervals.at("end")),
        weights = field(record, intervals.at("weight")), row_indices = field(record, intervals.at("source_record_index")),
        byte_starts = field(record, intervals.at("source_start")), byte_ends = field(record, intervals.at("source_end")),
        sample_indices = field(record, intervals.at("sample")), contigs = field(record, intervals.at("contig")),
        strands = field(record, intervals.at("strand")), id_offsets = field(record, intervals.at("id_offset")),
        id_lengths = field(record, intervals.at("id_length")), atom_starts = field(record, segments.at("start")),
        atom_ends = field(record, segments.at("end")), atom_contigs = field(record, segments.at("contig")),
        atom_strands = field(record, segments.at("strand")), metric = field(record, record.at("metric_diagonal"));
    const auto names = record.at("source_ids").bytes();
    std::map<std::string, std::uint64_t> reference;
    const auto contig_records = record.at("reference").at("contigs");
    for (std::size_t i = 0; i < contig_records.size(); ++i)
        check(reference.emplace(contig_records.at(i).at("id").text(),
                  integer(contig_records.at(i).at("length").integer().get_str())).second,
              "duplicate reference contig");
    std::map<std::string, std::size_t> sample_ids;
    std::vector<Column> sample_fields;
    Vector scales, supports(samples.size(), 0), integrals(samples.size(), 0);
    std::vector<std::size_t> sample_rows(samples.size());
    std::vector<Interval> rows;
    std::map<Frame, std::vector<Event>> frames;
    const auto contrasts = record.at("contrasts");
    std::vector<Column> inputs, represented, residuals, coefficients;
    std::vector<NativeRecord> frame_coordinates;
    Matrix terms(contrasts.size(), Vector(samples.size(), 0));
    for (std::size_t s = 0; s < samples.size(); ++s) {
        check(sample_ids.emplace(samples.at(s).at("id").text(), s).second, "duplicate sample");
        scales.emplace_back(samples.at(s).at("scale").text());
        check(scales.back() > 0, "invalid sample divisor");
        const auto item = record.at("sample_fields").at(s);
        check(item.at("sample").text() == samples.at(s).at("id").text() &&
                  item.at("scale_divisor").rational() == scales.back(), "sample field order or divisor");
        sample_fields.push_back(field(record, item.at("field")));
    }
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        const auto item = contrasts.at(c);
        check(item.at("id").text() == request.at("contrasts").at(c).at("id").text(), "contrast order");
        inputs.push_back(field(record, item.at("input")));
        represented.push_back(field(record, item.at("represented")));
        residuals.push_back(field(record, item.at("residual")));
        coefficients.push_back(field(record, item.at("source_coefficients")));
        frame_coordinates.push_back(item.at("frame_coordinates"));
        const auto declared = request.at("contrasts").at(c).at("terms");
        for (std::size_t t = 0; t < declared.size(); ++t)
            terms[c][sample_ids.at(declared.at(t).at("sample").text())] = Rational(declared.at(t).at("coefficient").text());
    }
    for (std::size_t s = 0; s < samples.size(); ++s) {
        const auto resource_name = samples.at(s).at("resource").text();
        std::optional<NativeRecord> resource;
        const auto resources = record.at("resources");
        for (std::size_t i = 0; i < resources.size(); ++i)
            if (resources.at(i).at("id").text() == resource_name) resource = resources.at(i);
        check(resource.has_value() && resource->at("declaration").at("format").text() == "bed",
              "verifier requires declared BED resources");
        const rfc::storage::ReadOnlyMappedFile original(argv[s + 4]);
        same_bytes(resource->at("bytes").bytes(), original.bytes());
        const auto raw = original.bytes();
        const std::string_view bytes(reinterpret_cast<const char *>(raw.data()), raw.size());
        std::size_t sample_row = 0;
        for (std::size_t start = 0; start < bytes.size();) {
            auto end = bytes.find('\n', start);
            if (end == std::string_view::npos) end = bytes.size();
            auto row = bytes.substr(start, end - start);
            if (!row.empty() && row.back() == '\r') row.remove_suffix(1);
            const auto next = end == bytes.size() ? end : end + 1;
            if (row.empty() || row.starts_with('#') || row.starts_with("track ") || row.starts_with("browser ")) {
                start = next; continue;
            }
            const auto tokens = split(row);
            check(tokens.size() >= 3 && tokens.size() <= 6, "BED profile requires three to six columns");
            const std::string contig(tokens[0]);
            const auto a = integer(tokens[1]), b = integer(tokens[2]);
            check(reference.contains(contig) && a < b && b <= reference.at(contig), "original interval bounds");
            const auto id = tokens.size() > 3 ? std::string(tokens[3]) : "record:" + std::to_string(sample_row + 1);
            const auto strand = tokens.size() > 5 ? std::string(tokens[5]) : ".";
            Rational weight = 1;
            if (resource->at("declaration").has("weight_column")) {
                const auto column = integer(resource->at("declaration").at("weight_column").integer().get_str());
                check(column >= 3 && column < tokens.size(), "weight column");
                weight = Rational(std::string(tokens[column]));
            }
            const auto j = rows.size();
            check(j < relations.cell_count() && integer(starts.at(j)) == a && integer(ends.at(j)) == b &&
                      weights.at(j) == weight && integer(row_indices.at(j)) == sample_row &&
                      integer(byte_starts.at(j)) == start && integer(byte_ends.at(j)) == end &&
                      integer(sample_indices.at(j)) == s &&
                      record.at("contigs").at(integer(contigs.at(j))).text() == contig &&
                      record.at("strands").at(integer(strands.at(j))).text() == strand,
                  "original row, weight, sample or byte span differs");
            const auto name_start = integer(id_offsets.at(j)), name_size = integer(id_lengths.at(j));
            check(name_start <= names.size() && name_size <= names.size() - name_start &&
                      std::string_view(reinterpret_cast<const char *>(names.data() + name_start), name_size) == id,
                  "original interval identity bytes differ");
            for (std::size_t c = 0; c < contrasts.size(); ++c)
                check(coefficients[c].at(j) == weight * terms[c][s] / scales[s], "original contrast coefficient differs");
            Frame frame{contig, strand_policy == "ignore" ? "ignored" : strand};
            rows.push_back({frame, a, b});
            frames[frame].push_back({a, s, weight, 1});
            frames[frame].push_back({b, s, -weight, -1});
            supports[s] += weight;
            integrals[s] += Rational(std::to_string(b - a)) * weight / scales[s];
            ++sample_row; start = next;
        }
        sample_rows[s] = sample_row;
    }
    check(rows.size() == relations.cell_count(), "extra or missing original intervals");
    const auto controls = record.at("controls"), basis = controls.at("independent_frame"), declared_frame = controls.at("declared_frame");
    const auto rank_count = integer(controls.at("rank").integer().get_str());
    std::vector<std::size_t> control_samples;
    for (std::size_t i = 0; i < request.at("controls").size(); ++i)
        control_samples.push_back(sample_ids.at(request.at("controls").at(i).text()));
    const auto control_count = control_samples.size() + (request.at("intercept").boolean() ? 1 : 0);
    check(basis.tensor().domain().extents() == std::vector<std::size_t>{atoms, rank_count} &&
              declared_frame.tensor().domain().extents() == std::vector<std::size_t>{atoms, control_count}, "control frame shape");
    Matrix control_gram(control_count, Vector(control_count, 0)), basis_gram(rank_count, Vector(rank_count, 0)),
        orthogonality(contrasts.size(), Vector(control_count, 0)), residual_gram(contrasts.size(), Vector(contrasts.size(), 0));
    std::vector<std::vector<bool>> candidates(rank_count, std::vector<bool>(control_count, true));
    Vector input_quadrance(contrasts.size(), 0), represented_quadrance(contrasts.size(), 0), residual_integral(contrasts.size(), 0), cross(contrasts.size(), 0), sample_integral(samples.size(), 0);
    std::map<Frame, std::pair<std::size_t, std::size_t>> ranges;
    std::size_t atom = 0;
    for (auto &[frame, events] : frames) {
        std::sort(events.begin(), events.end(), [](const auto &a, const auto &b) { return a.position < b.position; });
        Vector values(samples.size(), 0);
        std::int64_t active = 0;
        const auto begin = atom;
        for (std::size_t e = 0; e < events.size();) {
            const auto position = events[e].position;
            while (e < events.size() && events[e].position == position) {
                values[events[e].sample] += events[e].weight; active += events[e].count; ++e;
            }
            check(active >= 0, "negative original support count");
            if (e == events.size() || !active) continue;
            const auto end = events[e].position;
            check(atom < atoms && integer(atom_starts.at(atom)) == position && integer(atom_ends.at(atom)) == end &&
                      record.at("contigs").at(integer(atom_contigs.at(atom))).text() == frame.first &&
                      record.at("strands").at(integer(atom_strands.at(atom))).text() == frame.second,
                  "physical atom differs from original endpoint union");
            const Rational length(std::to_string(end - position));
            check(metric.at(atom) == length, "physical length metric differs");
            Vector scaled(samples.size(), 0), declared(control_count, 0), independent(rank_count, 0), observed_residuals(contrasts.size(), 0);
            for (std::size_t s = 0; s < samples.size(); ++s) {
                scaled[s] = values[s] / scales[s];
                check(sample_fields[s].at(atom) == scaled[s], "sample coverage differs from complete original sweep");
                sample_integral[s] += length * scaled[s];
            }
            for (std::size_t c = 0; c < control_samples.size(); ++c) declared[c] = scaled[control_samples[c]];
            if (request.at("intercept").boolean()) declared.back() = 1;
            for (std::size_t c = 0; c < control_count; ++c)
                check(entry(declared_frame, atom, c) == declared[c], "declared control field differs");
            for (std::size_t k = 0; k < rank_count; ++k) {
                independent[k] = entry(basis, atom, k);
                for (std::size_t c = 0; c < control_count; ++c)
                    candidates[k][c] = candidates[k][c] && independent[k] == declared[c];
            }
            for (std::size_t c = 0; c < control_count; ++c)
                for (std::size_t d = 0; d < control_count; ++d) control_gram[c][d] += length * declared[c] * declared[d];
            for (std::size_t c = 0; c < rank_count; ++c)
                for (std::size_t d = 0; d < rank_count; ++d) basis_gram[c][d] += length * independent[c] * independent[d];
            for (std::size_t c = 0; c < contrasts.size(); ++c) {
                Rational input = 0, component = 0;
                for (std::size_t s = 0; s < samples.size(); ++s) input += terms[c][s] * scaled[s];
                const auto &coordinates = frame_coordinates[c].tensor();
                check(coordinates.domain().logical_size() == rank_count, "projection coordinate shape");
                for (std::size_t k = 0; k < rank_count; ++k) component += independent[k] * coordinates.state(k).value();
                const Rational residual = input - component;
                check(inputs[c].at(atom) == input && represented[c].at(atom) == component && residuals[c].at(atom) == residual,
                      "contrast, represented component or residual differs");
                input_quadrance[c] += length * input * input;
                represented_quadrance[c] += length * component * component;
                residual_integral[c] += length * residual;
                cross[c] += length * component * residual;
                observed_residuals[c] = residual;
                for (std::size_t k = 0; k < control_count; ++k) orthogonality[c][k] += length * declared[k] * residual;
            }
            for (std::size_t c = 0; c < contrasts.size(); ++c)
                for (std::size_t d = 0; d < contrasts.size(); ++d) residual_gram[c][d] += length * observed_residuals[c] * observed_residuals[d];
            ++atom;
        }
        check(active == 0 && std::all_of(values.begin(), values.end(), [](const auto &value) { return value == 0; }), "incomplete original sweep");
        ranges.emplace(frame, std::pair{begin, atom});
    }
    check(atom == atoms && sample_integral == integrals, "complete atom count or source mass conservation");
    for (std::size_t s = 0; s < samples.size(); ++s)
        check(record.at("sample_fields").at(s).at("integral").rational() == integrals[s], "sample integral certificate differs");
    for (const auto &candidate : candidates)
        check(std::any_of(candidate.begin(), candidate.end(), [](bool value) { return value; }), "independent control column is outside the declared source frame");
    check(rank(control_gram) == rank_count && rank(basis_gram) == rank_count, "control span rank differs");
    for (std::size_t c = 0; c < rank_count; ++c)
        for (std::size_t d = 0; d < rank_count; ++d)
            check(entry(controls.at("gram"), c, d) == basis_gram[c][d], "control Gram differs");
    std::uint64_t occurrence = 0;
    for (std::size_t j = 0; j < rows.size(); ++j) {
        const auto [first, last] = ranges.at(rows[j].frame);
        std::size_t lo = first, hi = last;
        while (lo < hi) { const auto mid = lo + (hi - lo) / 2; if (integer(atom_starts.at(mid)) < rows[j].start) lo = mid + 1; else hi = mid; }
        const auto [begin, end] = relations.incidence_range(j);
        auto slot = begin;
        check(relations.relation_identities[j].value == j + 1, "original relation identity differs");
        for (auto i = lo; i < last && integer(atom_starts.at(i)) < rows[j].end; ++i) {
            check(slot < end && integer(atom_ends.at(i)) <= rows[j].end && relations.lower_index[slot] == i &&
                      relations.occurrence_identities[slot].value == ++occurrence &&
                      relations.aspects.get(source->source.signature().existence_aspect(), slot) &&
                      !relations.aspects.get(source->source.signature().orientation_aspect(), slot) &&
                      source->source.incidence_coefficient({1}, slot) == 1, "ordered original occurrence differs");
            ++slot;
        }
        check(slot == end, "extra original membership");
    }
    check(occurrence == relations.incidence_count(), "complete occurrence count");
    for (const auto &tensor : source->tensors) {
        check(tensor.output.has_value(), "missing field carrier");
        const auto carrier = std::find_if(source->carriers.begin(), source->carriers.end(), [&](const auto &item) { return item.key == *tensor.output; });
        check(carrier != source->carriers.end() && carrier->kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                  tensor.tensor.domain().logical_size() == source->source.grade(carrier->grade).cell_count() &&
                  tensor.tensor.support().count() == tensor.tensor.domain().logical_size(), "field carrier or measured support differs");
    }
    for (std::size_t c = 0; c < contrasts.size(); ++c) {
        check(std::all_of(orthogonality[c].begin(), orthogonality[c].end(), [](const auto &value) { return value == 0; }) && cross[c] == 0,
              "residual is not orthogonal to the original declared control span");
        const auto q = contrasts.at(c).at("quadrances");
        check(q.at("input").rational() == input_quadrance[c] && q.at("represented").rational() == represented_quadrance[c] &&
                  q.at("residual").rational() == residual_gram[c][c] && input_quadrance[c] == represented_quadrance[c] + residual_gram[c][c] &&
                  contrasts.at(c).at("residual_integral").rational() == residual_integral[c], "projection metric certificate differs");
        for (std::size_t d = 0; d < contrasts.size(); ++d)
            check(entry(record.at("residual_gram"), c, d) == residual_gram[c][d], "complete residual Gram differs");
    }
    check(record.at("thresholds").size() == request.at("thresholds").size(), "threshold count differs");
    for (std::size_t i = 0; i < record.at("thresholds").size(); ++i)
        check(record.at("thresholds").at(i).rational() == Rational(request.at("thresholds").at(i).text()), "threshold value differs");
    const rfc::storage::ReadOnlyMappedFile original_archive(argv[1]);
    const auto encoded = archive.encode();
    same_bytes(encoded, original_archive.bytes());
    std::cout << "sample\toriginal_records\treported_weight\tweighted_base_integral\n";
    for (std::size_t s = 0; s < samples.size(); ++s)
        std::cout << samples.at(s).at("id").text() << '\t' << sample_rows[s] << '\t' << supports[s] << '\t' << integrals[s] << '\n';
    std::cout << "contrast\tinput_quadrance\trepresented_quadrance\tresidual_quadrance\tresidual_integral\n";
    for (std::size_t c = 0; c < contrasts.size(); ++c)
        std::cout << contrasts.at(c).at("id").text() << '\t' << input_quadrance[c] << '\t' << represented_quadrance[c] << '\t' << residual_gram[c][c] << '\t' << residual_integral[c] << '\n';
    std::cout << "status\tPASS\noriginal_records\t" << rows.size() << "\nphysical_atoms\t" << atoms << "\noccurrences\t" << occurrence << "\ncontrol_rank\t" << rank_count << '\n';
    return 0;
}
int main(int argc, char **argv) {
    try { return verify(argc, argv); }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
