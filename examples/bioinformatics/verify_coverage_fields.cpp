#include <algorithm>
#include <charconv>
#include <fstream>
#include <gadget/native_record.hpp>
#include <iostream>
#include <map>
#include <rfc/core/grade_arena.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <set>
#include <string_view>
using namespace gadget;
using rfc::numeric::Rational;
namespace {
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
struct Column {
    const rfc::numeric::SparseRationalTensor *tensor;
    Rational at(std::size_t i) const {
        const auto value = tensor->state(i);
        check(value.has_value(), "missing measured coordinate");
        return *value;
    }
};
Column field(const NativeRecord &record, NativeRecord selector) {
    const auto key = std::stoull(selector.integer().get_str());
    for (const auto &value : record.endpoint()->tensors)
        if (value.key == key)
            return {&value.tensor};
    throw std::runtime_error("missing field owner");
}
std::uint64_t integer(std::string_view text) {
    std::uint64_t value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    check(error == std::errc{} && end == text.data() + text.size(), "invalid original coordinate");
    return value;
}
std::uint64_t integer(const Rational &value) {
    check(value.get_den() == 1 && value >= 0, "nonintegral native coordinate");
    return integer(value.get_num().get_str());
}
struct Interval {
    std::string contig;
    std::uint64_t start, end;
    Rational weight;
};
struct Summary {
    std::size_t rows{}, atoms{}, occurrences{};
    Rational support{}, integral{}, quadrance{};
};
struct Event {
    std::uint64_t position;
    Rational delta;
    int count;
};
} // namespace
int verify(int argc, char **argv) {
    std::string assay;
    if (argc >= 3 && std::string_view(argv[1]) == "--assay") {
        assay = argv[2];
        argc -= 2;
        argv += 2;
    }
    check(argc == 3 || argc == 5, "usage: verify_coverage_fields [--assay ID] RECORD.rfcxs "
                                  "ORIGINAL.bed [WORKSPACE.rfcxs REQUEST.rfcxs]");
    const auto archive = NativeRecord::read(argv[1]);
    auto owning = archive;
    if (!assay.empty()) {
        const auto schema = archive.at("schema").text();
        check(schema == "gadget.assay.workspace.native.v1" ||
                  schema == "gadget.assay.workspace.native.v2",
              "assay selection requires a native workspace");
        std::optional<NativeRecord> selected;
        const auto publications = archive.at("publications");
        for (std::size_t i = 0; i < publications.size(); ++i)
            if (publications.at(i).at("assay").text() == assay) {
                check(!selected.has_value(),
                      "select a standalone publication for partitioned assays");
                selected = NativeRecord::decode(publications.at(i).at("bytes").bytes());
            }
        check(selected.has_value(), "selected assay is absent from this workspace");
        owning = std::move(*selected);
    }
    const auto operation = owning.has("history")
                               ? owning.historical_publication(std::stoull(
                                     owning.at("history").at("cutoff").integer().get_str()))
                               : owning;
    const auto record = operation.at("result").at("result");
    check(record.at("schema").text() == "gadget.coverage.fields.v1", "wrong coverage schema");
    if (argc == 5) {
        const rfc::storage::ReadOnlyMappedFile workspace(argv[3]), original_request(argv[4]);
        const auto retained_workspace = operation.at("result").at("workspace_publication").bytes();
        const auto retained_request = operation.at("request_publication").bytes();
        check(retained_workspace.size() == workspace.size() &&
                  std::equal(retained_workspace.begin(), retained_workspace.end(),
                             workspace.bytes().begin()),
              "input workspace publication bytes changed");
        check(retained_request.size() == original_request.size() &&
                  std::equal(retained_request.begin(), retained_request.end(),
                             original_request.bytes().begin()),
              "input request publication bytes changed");
    }
    const auto request = record.at("request");
    check(request.at("samples").size() == 1 && request.at("contrasts").size() == 1 &&
              request.at("controls").size() == 0 && !request.at("intercept").boolean(),
          "oracle requires one complete identity observation");
    check(request.at("samples").at(0).at("scale").text() == "1" &&
              request.at("contrasts").at(0).at("terms").size() == 1 &&
              request.at("contrasts").at(0).at("terms").at(0).at("coefficient").text() == "1",
          "oracle identity coefficients");
    const auto original = record.at("resources").at(0).at("bytes").bytes();
    std::ifstream input(argv[2], std::ios::binary);
    check(bool(input), "original BED unavailable");
    const std::string bytes((std::istreambuf_iterator<char>(input)), {});
    check(original.size() == bytes.size() &&
              std::equal(original.begin(), original.end(),
                         reinterpret_cast<const std::byte *>(bytes.data())),
          "original byte owner changed");
    const auto intervals = record.at("intervals"), segments = record.at("segments");
    const auto starts = field(record, intervals.at("start")),
               ends = field(record, intervals.at("end")),
               weights = field(record, intervals.at("weight")),
               row_indices = field(record, intervals.at("source_record_index")),
               source_start = field(record, intervals.at("source_start")),
               source_end = field(record, intervals.at("source_end"));
    const auto segment_start = field(record, segments.at("start")),
               segment_end = field(record, segments.at("end"));
    const auto sample_indices = field(record, intervals.at("sample")),
               contigs = field(record, intervals.at("contig")),
               strands = field(record, intervals.at("strand")),
               id_offsets = field(record, intervals.at("id_offset")),
               id_lengths = field(record, intervals.at("id_length")),
               segment_contigs = field(record, segments.at("contig")),
               segment_strands = field(record, segments.at("strand"));
    const auto names = record.at("source_ids").bytes();

    const auto contrast = record.at("contrasts").at(0);
    const auto sample = field(record, record.at("sample_fields").at(0).at("field")),
               coverage = field(record, contrast.at("input")),
               residual = field(record, contrast.at("residual")),
               represented = field(record, contrast.at("represented")),
               coefficients = field(record, contrast.at("source_coefficients")),
               metric = field(record, record.at("metric_diagonal"));
    check(request.at("metric").text() == "reference_length" &&
              request.at("strand_policy").text() == "ignore" &&
              request.at("absence_policy").text() == "measured_zero_within_union",
          "oracle physical frame contract");
    check(request.at("contrasts").at(0).at("terms").at(0).at("sample").text() ==
              request.at("samples").at(0).at("id").text(),
          "identity contrast sample differs");
    std::map<std::string, std::uint64_t> lengths;
    const auto reference = record.at("reference").at("contigs");
    for (std::size_t i = 0; i < reference.size(); ++i)
        check(lengths
                  .emplace(reference.at(i).at("id").text(),
                           integer(reference.at(i).at("length").integer().get_str()))
                  .second,
              "duplicate declared contig");
    std::vector<Interval> rows;
    std::map<std::string, std::vector<Event>> frames;
    std::map<std::string, Summary> summaries;
    Rational support = 0, integral = 0;
    for (std::size_t start = 0; start < bytes.size();) {
        auto end = bytes.find('\n', start);
        if (end == std::string::npos)
            end = bytes.size();
        const std::string_view row(bytes.data() + start, end - start);
        std::vector<std::string_view> tokens;
        for (std::size_t p = 0; p <= row.size();) {
            auto q = row.find('\t', p);
            if (q == std::string_view::npos)
                q = row.size();
            tokens.push_back(row.substr(p, q - p));
            if (q == row.size())
                break;
            p = q + 1;
        }
        check(tokens.size() == 5, "oracle requires five column fragments BED");
        const std::string contig(tokens[0]);
        const auto a = integer(tokens[1]), b = integer(tokens[2]);
        const Rational weight{std::string(tokens[4])};
        check(lengths.contains(contig) && a < b && b <= lengths.at(contig) && weight >= 0 &&
                  weight.get_den() == 1,
              "source fragment bounds or multiplicity");
        const auto j = rows.size();
        check(integer(starts.at(j)) == a && integer(ends.at(j)) == b && weights.at(j) == weight &&
                  coefficients.at(j) == weight,
              "original interval field differs");
        check(integer(row_indices.at(j)) == j && integer(source_start.at(j)) == start &&
                  integer(source_end.at(j)) == end,
              "original row span/order differs");
        check(integer(sample_indices.at(j)) == 0 &&
                  record.at("contigs").at(integer(contigs.at(j))).text() ==
                      std::string(tokens[0]) &&
                  record.at("strands").at(integer(strands.at(j))).text() == ".",
              "original source sample/contig/strand differs");
        const auto name_start = integer(id_offsets.at(j)), name_length = integer(id_lengths.at(j));
        check(name_start <= names.size() && name_length <= names.size() - name_start &&
                  std::string_view(reinterpret_cast<const char *>(names.data() + name_start),
                                   name_length) == tokens[3],
              "original source id bytes differ");
        rows.push_back({contig, a, b, weight});
        frames[contig].push_back({a, weight, 1});
        frames[contig].push_back({b, -weight, -1});
        auto &summary = summaries[contig];
        ++summary.rows;
        summary.support += weight;
        summary.integral += Rational(std::to_string(b - a)) * weight;
        support += weight;
        integral += Rational(std::to_string(b - a)) * weight;
        start = end == bytes.size() ? end : end + 1;
    }
    check(rows.size() == record.endpoint()->source.grade({1}).cell_count(),
          "complete interval count");
    Rational quadrance = 0;
    std::size_t atom = 0;
    std::vector<std::uint64_t> atom_starts, atom_ends;
    std::map<std::string, std::pair<std::size_t, std::size_t>> frame_ranges;
    for (auto &[contig, events] : frames) {
        std::sort(events.begin(), events.end(),
                  [](const auto &a, const auto &b) { return a.position < b.position; });
        Rational value = 0;
        std::int64_t active = 0;
        const auto frame_begin = atom;
        for (std::size_t i = 0; i < events.size();) {
            const auto position = events[i].position;
            while (i < events.size() && events[i].position == position) {
                value += events[i].delta;
                active += events[i].count;
                ++i;
            }
            check(active >= 0, "source sweep support");
            if (i < events.size() && active > 0) {
                const auto end = events[i].position;
                check(record.at("contigs").at(integer(segment_contigs.at(atom))).text() == contig &&
                          record.at("strands").at(integer(segment_strands.at(atom))).text() ==
                              "ignored",
                      "physical segment reference frame differs");
                check(integer(segment_start.at(atom)) == position &&
                          integer(segment_end.at(atom)) == end,
                      "native physical segment differs");
                check(sample.at(atom) == value && coverage.at(atom) == value &&
                          residual.at(atom) == value && represented.at(atom) == 0,
                      "complete native field differs from independent original sweep");
                check(integer(metric.at(atom)) == end - position, "physical metric differs");
                const Rational q = Rational(std::to_string(end - position)) * value * value;
                quadrance += q;
                summaries[contig].quadrance += q;
                ++summaries[contig].atoms;
                atom_starts.push_back(position);
                atom_ends.push_back(end);
                ++atom;
            }
        }
        check(active == 0 && value == 0, "complete source sweep closure");
        frame_ranges.emplace(contig, std::pair{frame_begin, atom});
    }
    check(atom == record.endpoint()->source.grade({0}).cell_count(),
          "complete physical atom count");
    for (const auto &tensor : record.endpoint()->tensors) {
        check(tensor.output.has_value(), "field has no original source carrier");
        const auto &carriers = record.endpoint()->carriers;
        const auto carrier = std::find_if(carriers.begin(), carriers.end(),
                                          [&](const auto &c) { return c.key == *tensor.output; });
        check(carrier != carriers.end() &&
                  carrier->kind == rfc::storage::TemporalCarrierKind::SourceGrade,
              "field uses a foreign carrier");
        check(tensor.tensor.domain().logical_size() ==
                      record.endpoint()->source.grade(carrier->grade).cell_count() &&
                  tensor.tensor.support().count() == tensor.tensor.domain().logical_size(),
              "field domain or complete measured support differs");
    }
    const auto &grade = record.endpoint()->source.grade({1});
    std::uint64_t occurrence = 0;
    for (std::size_t j = 0; j < rows.size(); ++j) {
        const auto [begin, end] = grade.incidence_range(j);
        auto slot = begin;
        const auto [frame_begin, frame_end] = frame_ranges.at(rows[j].contig);
        auto i = std::lower_bound(atom_starts.begin() + frame_begin,
                                  atom_starts.begin() + frame_end, rows[j].start) -
                 atom_starts.begin();
        check(grade.relation_identities[j].value == j + 1, "original relation identity");
        for (; static_cast<std::size_t>(i) < frame_end && atom_starts[i] < rows[j].end; ++i) {
            check(slot < end && atom_ends[i] <= rows[j].end &&
                      grade.lower_index[slot] == static_cast<std::size_t>(i),
                  "complete occurrence membership differs");
            check(grade.occurrence_identities[slot].value == ++occurrence &&
                      grade.aspects.get(record.endpoint()->source.signature().existence_aspect(),
                                        slot) &&
                      !grade.aspects.get(record.endpoint()->source.signature().orientation_aspect(),
                                         slot) &&
                      record.endpoint()->source.incidence_coefficient({1}, slot) == 1,
                  "original occurrence identity/flags/coefficient differs");
            ++slot;
            ++summaries[rows[j].contig].occurrences;
        }
        check(slot == end, "extra original occurrence");
    }
    check(occurrence == grade.incidence_count(), "complete native occurrence count");
    check(contrast.at("residual_integral").rational() == integral &&
              contrast.at("quadrances").at("input").rational() == quadrance &&
              contrast.at("quadrances").at("residual").rational() == quadrance,
          "complete exact integral/quadrance");
    check(contrast.at("quadrances").at("represented").rational() == 0 &&
              record.at("residual_gram").tensor().state(0).value() == quadrance,
          "native zero component or residual Gram differs");
    check(record.at("controls").at("rank").integer() == 0 &&
              record.at("controls").at("independent_frame").tensor().domain().extents() ==
                  std::vector<std::size_t>({atom, 0}),
          "native empty control span differs");
    check(contrast.at("frame_coordinates").tensor().domain().logical_size() == 0,
          "native empty frame coordinate domain");
    check(record.at("thresholds").size() == request.at("thresholds").size(),
          "declared threshold count differs");
    for (std::size_t i = 0; i < request.at("thresholds").size(); ++i)
        check(record.at("thresholds").at(i).rational() ==
                  Rational(request.at("thresholds").at(i).text()),
              "declared threshold differs");
    const auto encoded = archive.encode();
    const rfc::storage::ReadOnlyMappedFile archived(argv[1]);
    check(encoded.size() == archived.size() &&
              std::equal(encoded.begin(), encoded.end(), archived.bytes().begin()),
          "complete native reencoding differs");
    std::cout << "contig\trecords\tatoms\toccurrences\treported_support\tweighted_base_"
                 "integral\tquadrance\n";
    for (const auto &[contig, summary] : summaries)
        std::cout << contig << '\t' << summary.rows << '\t' << summary.atoms << '\t'
                  << summary.occurrences << '\t' << summary.support << '\t' << summary.integral
                  << '\t' << summary.quadrance << '\n';
    std::size_t raw_bytes = 0;
    for (const auto &owner : operation.raw_sources()) {
        raw_bytes = rfc::machine::checked_size_add(raw_bytes, owner.bytes().size());
    }
    std::size_t entry_slots = 0, entry_objects = 0, entry_limbs = 0;
    for (const auto &field : operation.endpoint()->tensors) {
        const auto &entries = field.tensor.values().entries();
        entry_slots += entries.size();
        entry_objects += entries.capacity() * sizeof(rfc::numeric::SparseRationalEntry);
        for (const auto &entry : entries)
            entry_limbs += rfc::numeric::rational_limb_bytes(entry.value);
    }
    std::size_t lattice_objects = 0, lattice_limbs = 0;
    for (std::size_t g = 0; g < operation.endpoint()->source.grade_count(); ++g) {
        const auto &lane =
            operation.endpoint()->source.grade({static_cast<std::uint16_t>(g)}).lattice_coefficient;
        lattice_objects += lane.capacity() * sizeof(rfc::numeric::Integer);
        for (const auto &value : lane)
            lattice_limbs += rfc::numeric::integer_limb_bytes(value);
    }
    std::cout << "retained_raw_owner_bytes\t" << raw_bytes << "\nfield_nonzero_entries\t"
              << entry_slots << "\nfield_entry_object_capacity_bytes\t" << entry_objects
              << "\nfield_gmp_limb_capacity_bytes\t" << entry_limbs
              << "\nlattice_object_capacity_bytes\t" << lattice_objects
              << "\nlattice_gmp_limb_capacity_bytes\t" << lattice_limbs << '\n';
    std::cout << "status\tPASS\noriginal_records\t" << rows.size() << "\nphysical_segments\t"
              << atom << "\noccurrences\t" << occurrence << "\nreported_support\t"
              << support.get_str() << "\nweighted_base_integral\t" << integral.get_str()
              << "\nquadrance\t" << quadrance.get_str() << '\n';
    return 0;
}
int main(int argc, char **argv) {
    try {
        return verify(argc, argv);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
