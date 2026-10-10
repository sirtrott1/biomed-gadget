#include <algorithm>
#include <fstream>
#include <gadget/native_record.hpp>
#include <iostream>
#include <map>
#include <numeric>
#include <rfc/storage/file_io.hpp>
#include <sstream>

using gadget::NativeRecord;
using rfc::numeric::Rational;
namespace {
void check(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
std::vector<std::string> split(const std::string &line) {
    std::vector<std::string> parts;
    std::istringstream stream(line);
    for (std::string part; std::getline(stream, part, '\t');) parts.push_back(part);
    return parts;
}
struct Axis {
    std::string bytes;
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;
    std::vector<std::pair<std::size_t, std::size_t>> spans;
};
Axis axis(const char *path) {
    std::ifstream file(path, std::ios::binary);
    check(bool(file), "original axis is missing");
    Axis result;
    result.bytes.assign(std::istreambuf_iterator<char>(file), {});
    std::size_t start = 0;
    while (start < result.bytes.size()) {
        const auto end = result.bytes.find('\n', start);
        const auto stop = end == std::string::npos ? result.bytes.size() : end;
        auto text = result.bytes.substr(start, stop - start);
        if (!text.empty() && text.back() == '\r') text.pop_back();
        auto row = split(text);
        if (start == 0) {
            check(!row.empty() && row[0] == "id", "original axis requires an id header");
            result.columns = std::move(row);
        } else {
            check(row.size() == result.columns.size(), "original axis width differs");
            result.rows.push_back(std::move(row));
            result.spans.emplace_back(start, stop);
        }
        start = stop + (end == std::string::npos ? 0 : 1);
    }
    return result;
}
std::size_t column(const Axis &axis, const std::string &name) {
    const auto p = std::find(axis.columns.begin(), axis.columns.end(), name);
    check(p != axis.columns.end(), "original coordinate column missing");
    return std::size_t(p - axis.columns.begin());
}
const rfc::numeric::SparseRationalTensor &field(const NativeRecord &publication, std::uint64_t key) {
    for (const auto &tensor : publication.endpoint()->tensors)
        if (tensor.key == key) {
            check(tensor.role == rfc::storage::TemporalTensorRole::Field && tensor.output && !tensor.input,
                  "spatial field must use its declared source carrier");
            const auto &carriers = publication.endpoint()->carriers;
            const auto c = std::find_if(carriers.begin(), carriers.end(), [&](const auto &x) {
                return x.key == *tensor.output;
            });
            check(c != carriers.end() && c->kind == rfc::storage::TemporalCarrierKind::SourceGrade &&
                      c->domain == tensor.tensor.domain() && c->domain.rank() == 1 &&
                      c->domain.logical_size() == publication.endpoint()->source.grade(c->grade).cell_count(),
                  "spatial field carrier mismatch");
            return tensor.tensor;
        }
    throw std::runtime_error("spatial field is missing");
}
std::uint64_t key(const NativeRecord &record) { return std::stoull(record.text()); }
void original_bytes(const NativeRecord &span, const char *path) {
    const auto expected = rfc::storage::read_file_exact(path, std::filesystem::file_size(path));
    const auto actual = span.bytes();
    check(expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin()),
          "generating publication differs from original bytes");
}
} // namespace
int main(int argc, char **argv) try {
    std::string assay;
    if (argc > 2 && std::string_view(argv[1]) == "--assay") {
        assay = argv[2]; argc -= 2; argv += 2;
    }
    check(argc == 7, "usage: verify_spatial_fields [--assay ID] RECORD.rfcxs MATRIX.mtx "
                     "OBSERVATIONS.tsv FEATURES.tsv WORKSPACE.rfcxs REQUEST.rfcxs");
    const auto archive = NativeRecord::read(argv[1]);
    auto publication = archive;
    if (!assay.empty()) {
        bool found = false;
        const auto publications = archive.at("publications");
        for (std::size_t i = 0; i < publications.size(); ++i)
            if (publications.at(i).at("assay").text() == assay) {
                check(!found, "spatial assay must have one complete carrier"); found = true;
                publication = NativeRecord::decode(publications.at(i).at("bytes").bytes());
            }
        check(found, "spatial assay missing");
    }
    if (publication.has("history"))
        publication = publication.historical_publication(
            std::stoull(publication.at("history").at("cutoff").integer().get_str()));
    const auto generated = publication.at("result");
    const auto result = generated.at("result"), request = generated.at("query");
    check(result.at("schema").text() == "gadget.spatial.fields.v1", "spatial producer schema");
    original_bytes(generated.at("workspace_publication"), argv[5]);
    original_bytes(publication.at("request_publication"), argv[6]);
    const auto workspace = NativeRecord::decode(generated.at("workspace_publication").bytes());
    const auto assays = workspace.at("manifest").at("assays");
    std::optional<NativeRecord> original;
    for (std::size_t i = 0; i < assays.size(); ++i)
        if (assays.at(i).at("declaration").at("id").text() == request.at("assay").text()) original = assays.at(i);
    check(original.has_value(), "generating source assay missing");
    const auto retained = workspace.at("originals");
    for (const auto &[name, path] : std::vector<std::pair<std::string, const char *>>{
             {"matrix", argv[2]}, {"observations", argv[3]}, {"features", argv[4]}}) {
        const auto wanted = original->at("originals").at(name).at("path").text();
        bool found = false;
        for (std::size_t i = 0; i < retained.size(); ++i)
            if (retained.at(i).at("path").text() == wanted) {
                original_bytes(retained.at(i).at("bytes"), path); found = true;
            }
        check(found, "original scientific artifact missing");
    }
    const auto observations = axis(argv[3]), features = axis(argv[4]);
    const auto &source = publication.endpoint()->source;
    std::optional<NativeRecord> original_publication;
    const auto inputs = workspace.at("publications");
    for (std::size_t i = 0; i < inputs.size(); ++i)
        if (inputs.at(i).at("assay").text() == request.at("assay").text()) {
            check(!original_publication, "original assay carrier is partitioned");
            original_publication = NativeRecord::decode(inputs.at(i).at("bytes").bytes());
        }
    check(original_publication.has_value(), "original measured publication missing");
    const auto &old = original_publication->endpoint()->source;
    check(old.grade_count() <= source.grade_count() && old.semantic_scale() == source.semantic_scale(),
          "original graded carrier or semantic scale changed");
    for (std::size_t k = 0; k < old.grade_count(); ++k) {
        const rfc::core::Grade grade{static_cast<std::uint16_t>(k)};
        const auto &before = old.grade(grade), &after = source.grade(grade);
        check(before.cell_count() <= after.cell_count(), "original source cell missing");
        for (std::size_t i = 0; i < before.cell_count(); ++i) {
            check((k ? before.relation_identities[i].value : before.identities[i]) ==
                      (k ? after.relation_identities[i].value : after.identities[i]),
                  "original primary cell identity changed");
            if (!k) continue;
            check(before.relation_types[i] == after.relation_types[i], "original relation type changed");
            const auto [a, b] = before.incidence_range(i);
            const auto [c, d] = after.incidence_range(i);
            check(a == c && b == d, "original attaching order changed");
            for (auto p = a; p < b; ++p) {
                check(before.occurrence_identities[p] == after.occurrence_identities[p] &&
                      before.lower_index[p] == after.lower_index[p] &&
                      old.incidence_coefficient(grade, p) == source.incidence_coefficient(grade, p),
                      "original occurrence identity, participant or exact coefficient changed");
                for (std::size_t bit = 0; bit < old.signature().aspects().size(); ++bit)
                    check(before.aspects.get(static_cast<rfc::composite::AspectId>(bit), p) ==
                          after.aspects.get(static_cast<rfc::composite::AspectId>(bit), p),
                          "original occurrence aspect changed");
            }
        }
    }
    const auto n = observations.rows.size(), genes = features.rows.size();
    check(source.grade({0}).cell_count() == n, "complete original spot carrier differs");
    std::map<std::string, std::uint64_t> keys;
    const auto additions = result.at("fields");
    for (std::size_t i = 0; i < additions.size(); ++i)
        check(keys.emplace(additions.at(i).at("id").text(), key(additions.at(i).at("key"))).second,
              "duplicate spatial field identity");
    const auto &starts = field(publication, keys.at("spatial.source_start"));
    const auto &ends = field(publication, keys.at("spatial.source_end"));
    const auto &ids = field(publication, keys.at("spatial.source_document"));
    const auto raw = result.at("source_documents");
    bool axis_bound = false;
    for (std::size_t i = 0; i < raw.size(); ++i)
        if (raw.at(i).at("id").text() == "carrier_axis") {
            const auto bytes = raw.at(i).at("bytes").bytes();
            check(bytes.size() == observations.bytes.size() &&
                std::equal(bytes.begin(), bytes.end(), reinterpret_cast<const std::byte *>(observations.bytes.data())),
                "original observation bytes differ"); axis_bound = true;
        }
    check(axis_bound, "original observation byte owner missing");
    const auto dimensions = request.at("coordinates").size();
    std::vector<std::vector<Rational>> positions(n);
    std::vector<Rational> metric;
    for (std::size_t d = 0; d < dimensions; ++d) {
        const auto name = request.at("coordinates").at(d).text();
        const auto c = column(observations, name);
        const auto &coordinates = field(publication, keys.at("spatial.coordinate:" + name));
        metric.emplace_back(request.at("metric_diagonal").at(d).text());
        check(metric.back() > 0 && coordinates.support().count() == n, "coordinate support or metric");
        for (std::size_t i = 0; i < n; ++i) {
            positions[i].emplace_back(observations.rows[i][c]);
            check(coordinates.state(i).value() == positions[i].back(), "original physical coordinate differs");
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        check(starts.state(i).value() == Rational(std::to_string(observations.spans[i].first)) &&
              ends.state(i).value() == Rational(std::to_string(observations.spans[i].second)) &&
              ids.state(i).value() == 0, "original observation span differs");
        check(observations.rows[i][column(observations, "coordinate_frame")] == request.at("frame").text() &&
              observations.rows[i][column(observations, "coordinate_units")] == request.at("units").text(),
              "original coordinate frame or units differ");
    }
    const auto begin = std::stoull(result.at("neighbor_begin").integer().get_str());
    const auto &arena = source.grade({1});
    const auto &distances = field(publication, keys.at("spatial.squared_distance"));
    const Rational radius(request.at("radius_squared").text());
    std::size_t accepted = 0, pairs = 0;
    std::vector<std::size_t> parent(n); std::iota(parent.begin(), parent.end(), 0);
    const auto root = [&](std::size_t v) { while (parent[v] != v) v = parent[v]; return v; };
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            ++pairs;
            Rational q{0};
            for (std::size_t d = 0; d < dimensions; ++d) {
                const Rational delta = positions[j][d] - positions[i][d];
                q += metric[d] * delta * delta;
            }
            if (q > radius) continue;
            const auto index = begin + accepted++;
            check(index < arena.cell_count(), "original admissible neighbor missing");
            const auto [a, b] = arena.incidence_range(index);
            check(b == a + 2 && arena.lower_index[a] == i && arena.lower_index[a + 1] == j &&
                  source.incidence_coefficient({1}, a) == -1 && source.incidence_coefficient({1}, a + 1) == 1,
                  "ordered physical boundary differs");
            check(arena.aspects.get(source.signature().existence_aspect(), a) &&
                  arena.aspects.get(source.signature().existence_aspect(), a + 1) &&
                  arena.aspects.get(source.signature().orientation_aspect(), a) &&
                  !arena.aspects.get(source.signature().orientation_aspect(), a + 1),
                  "physical occurrence existence or orientation differs");
            check(source.signature().type(arena.relation_types[index]).name == request.at("relation_type").text(),
                  "physical relation type differs");
            check(distances.state(index).value() == q, "native distance field differs");
            parent[root(j)] = root(i);
        }
    check(begin + accepted == arena.cell_count() && distances.support().count() == accepted &&
          accepted == std::stoull(result.at("neighbor_count").integer().get_str()), "complete neighbor support differs");
    const auto meanings = result.at("source_study").at("fields");
    std::map<std::string, std::uint64_t> measured;
    for (std::size_t i = 0; i < meanings.size(); ++i)
        measured.emplace(meanings.at(i).at("id").text(), key(meanings.at(i).at("key")));
    using Entry = std::pair<std::size_t, Rational>;
    std::vector<std::vector<Entry>> counts(genes);
    std::ifstream matrix(argv[2]); check(bool(matrix), "original matrix missing");
    std::string line; bool shape = false; std::size_t records = 0, expected_records = 0;
    Rational total{0};
    while (std::getline(matrix, line)) {
        if (line.empty() || line.front() == '%') continue;
        std::istringstream row(line);
        if (!shape) {
            std::size_t rows{}, cols{};
            check(bool(row >> rows >> cols >> expected_records) && rows == genes && cols == n, "original matrix shape differs");
            shape = true; continue;
        }
        std::size_t gene{}, spot{}; std::string value;
        check(bool(row >> gene >> spot >> value) && gene > 0 && gene <= genes && spot > 0 && spot <= n,
              "original matrix coordinate invalid");
        Rational q(value); check(q >= 0 && q.get_den() == 1, "original UMI count invalid");
        counts[gene - 1].emplace_back(spot - 1, q); total += q; ++records;
    }
    check(shape && records == expected_records, "original matrix records incomplete");
    std::size_t nonzero = 0;
    for (std::size_t gene = 0; gene < genes; ++gene) {
        auto &entries = counts[gene];
        std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
        const auto &tensor = field(publication, measured.at("all." + features.rows[gene][0]));
        check(tensor.domain().logical_size() == n && tensor.support().count() == n, "complete measured-zero support differs");
        const auto &actual = tensor.values().entries();
        std::size_t index = 0;
        for (std::size_t p = 0; p < entries.size();) {
            const auto spot = entries[p].first; Rational q{0};
            do { q += entries[p++].second; } while (p < entries.size() && entries[p].first == spot);
            if (q == 0) continue;
            check(index < actual.size() && actual[index].index == spot && actual[index].value == q,
                  "original measured expression differs"); ++index;
        }
        check(index == actual.size(), "unexpected measured expression coordinate"); nonzero += index;
        std::vector<Entry>().swap(entries);
    }
    const auto charts = result.at("charts");
    if (charts.size()) {
        check(charts.size() == n, "complete original witness cover differs");
        for (std::size_t i = 0; i < n; ++i) {
            const auto b = std::find_if(publication.endpoint()->local_blocks.begin(), publication.endpoint()->local_blocks.end(),
                [&](const auto &block) { return block.block_key == key(charts.at(i).at("block")); });
            check(b != publication.endpoint()->local_blocks.end() && b->kind == rfc::storage::LocalBlockKind::StrictSubcomplex,
                  "native strict chart missing");
            std::vector<rfc::storage::StableCellKey> cells{{{0}, source.grade({0}).identities[i]}};
            std::vector<rfc::storage::StableOccurrenceKey> occurrences;
            for (std::size_t edge = 0; edge < arena.cell_count(); ++edge) {
                const auto [a, z] = arena.incidence_range(edge);
                bool incident = false;
                for (auto p = a; p < z; ++p) incident |= arena.lower_index[p] == i;
                if (!incident) continue;
                cells.push_back({{1}, arena.relation_identities[edge].value});
                for (auto p = a; p < z; ++p) {
                    cells.push_back({{0}, source.grade({0}).identities[arena.lower_index[p]]});
                    occurrences.push_back({{1}, arena.occurrence_identities[p]});
                }
            }
            std::sort(cells.begin(), cells.end()); cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
            std::sort(occurrences.begin(), occurrences.end());
            check(cells == b->cells && occurrences == b->occurrences, "complete original chart closure differs");
        }
    }
    const auto original_archive = rfc::storage::read_file_exact(argv[1], std::filesystem::file_size(argv[1]));
    const auto encoded = archive.encode();
    check(original_archive == encoded, "canonical archive round trip differs");
    std::size_t components = 0;
    for (std::size_t i = 0; i < n; ++i) components += root(i) == i;
    std::cout << "PASS\nspots\t" << n << "\ngenes\t" << genes << "\nmatrix_records\t" << records
        << "\nmeasured_nonzero\t" << nonzero << "\numi_total\t" << total.get_str()
        << "\npairs_checked\t" << pairs << "\nneighbors\t" << accepted
        << "\nordered_occurrences\t" << arena.incidence_count() << "\ncharts\t" << charts.size()
        << "\nphysical_components\t" << components << '\n';
    if (begin == 0) std::cout << "physical_cycle_rank\t" << accepted + components - n << '\n';
} catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
