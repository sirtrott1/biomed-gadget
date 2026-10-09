#include "assay_detail.hpp"
#include "genomic_detail.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <sstream>
#include <string_view>

namespace gadget::assay_detail {
namespace {
// Retain the original bytes without retaining a separate index of every row.
class MatrixRecords {
  public:
    explicit MatrixRecords(std::string_view bytes) : bytes_(bytes) {
        require(bytes.find('\0') == std::string_view::npos, "NUL in text source");
    }
    bool next(std::string_view &record) {
        if (position_ == bytes_.size())
            return false;
        auto end = bytes_.find('\n', position_);
        if (end == std::string_view::npos)
            end = bytes_.size();
        auto stop = end;
        if (stop > position_ && bytes_[stop - 1] == '\r')
            --stop;
        require(stop > position_, "empty text record");
        record = bytes_.substr(position_, stop - position_);
        position_ = end == bytes_.size() ? end : end + 1;
        return true;
    }

  private:
    std::string_view bytes_;
    std::size_t position_{};
};
std::array<std::string_view, 3> matrix_tokens(std::string_view record) {
    std::array<std::string_view, 3> tokens;
    std::size_t position = 0;
    const auto space = [](char c) { return std::isspace(static_cast<unsigned char>(c)); };
    for (auto &token : tokens) {
        while (position < record.size() && space(record[position]))
            ++position;
        const auto begin = position;
        while (position < record.size() && !space(record[position]))
            ++position;
        require(position > begin, "Matrix Market record requires three tokens");
        token = record.substr(begin, position - begin);
    }
    while (position < record.size() && space(record[position]))
        ++position;
    require(position == record.size(), "Matrix Market record requires three tokens");
    return tokens;
}
std::uint64_t matrix_unsigned(std::string_view token) {
    std::uint64_t value{};
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
    require(!token.empty() && error == std::errc{} && end == token.data() + token.size(),
            "expected an unsigned integer token");
    return value;
}
} // namespace
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> out;
    std::size_t p = 0;
    for (;;) {
        auto q = s.find(delimiter, p);
        if (q == std::string::npos)
            q = s.size();
        out.push_back(s.substr(p, q - p));
        if (q == s.size())
            return out;
        p = q + 1;
    }
}
std::vector<std::pair<std::size_t, std::size_t>> lines(const std::string &bytes) {
    require(bytes.find('\0') == std::string::npos, "NUL in text source");
    std::vector<std::pair<std::size_t, std::size_t>> out;
    for (std::size_t p = 0; p < bytes.size();) {
        auto q = bytes.find('\n', p);
        if (q == std::string::npos)
            q = bytes.size();
        auto end = q;
        if (end > p && bytes[end - 1] == '\r')
            --end;
        require(end > p, "empty text record");
        out.emplace_back(p, end);
        p = q == bytes.size() ? q : q + 1;
    }
    return out;
}
std::uint64_t unsigned_token(const std::string &s) {
    std::uint64_t n{};
    const auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), n);
    require(!s.empty() && ec == std::errc{} && end == s.data() + s.size(),
            "expected an unsigned integer token");
    return n;
}
Json file_record(const std::string &path, const std::string &bytes) {
    auto out = Json::object();
    out.set("path", Json::string(path));
    out.set("sha256", Json::string("sha256:" + sha256(bytes)));
    return out;
}
Json span(const std::string &document, const std::string &bytes, std::size_t start,
          std::size_t end) {
    auto out = Json::array(), record = Json::object();
    record.set("document", Json::string(document));
    record.set("start", cardinality(start));
    record.set("end", cardinality(end));
    record.set("text", Json::string(bytes.substr(start, end - start)));
    out.push(record);
    return out;
}
Axis read_axis(const std::string &bytes) {
    Axis out;
    out.bytes = bytes;
    const auto records = lines(bytes);
    require(records.size() >= 2, "axis requires a header and identities");
    const auto header =
        split(bytes.substr(records[0].first, records[0].second - records[0].first), '\t');
    require(header.front() == "id", "axis first column must be id");
    std::set<std::string> headers, ids;
    for (const auto &name : header)
        require(!name.empty() && headers.insert(name).second, "axis header");
    for (std::size_t i = 1; i < records.size(); ++i) {
        const auto [p, q] = records[i];
        const auto cells = split(bytes.substr(p, q - p), '\t');
        require(cells.size() == header.size(), "axis metadata columns");
        require(!cells[0].empty() && ids.insert(cells[0]).second,
                "axis identities must be nonempty and unique");
        auto record = Json::object();
        for (std::size_t c = 0; c < cells.size(); ++c)
            record.set(header[c], Json::string(cells[c]));
        out.ids.push_back(cells[0]);
        out.records.push(record);
        out.spans.push_back(records[i]);
    }
    return out;
}
Json reference_record(const Json &value) {
    members(value, {"assembly", "sequence_sha256", "coordinate_convention", "contigs"});
    (void)text(value.at("assembly"));
    const auto hash = text(value.at("sequence_sha256"));
    require(hash.size() == 71 && hash.starts_with("sha256:"), "reference sequence hash");
    for (auto ch : hash.substr(7))
        require(std::isxdigit(static_cast<unsigned char>(ch)) && !(ch >= 'A' && ch <= 'F'),
                "reference hash encoding");
    require(value.at("coordinate_convention").text() == "zero_based_half_open",
            "canonical resource coordinates must be zero-based half-open");
    const auto contigs = value.at("contigs");
    require(contigs.size() > 0, "reference requires contigs");
    std::set<std::string> seen;
    for (std::size_t i = 0; i < contigs.size(); ++i) {
        const auto c = contigs.at(i);
        members(c, {"id", "length"});
        require(seen.insert(text(c.at("id"))).second, "duplicate reference contig");
        require(c.at("length").integer() > 0, "reference length");
    }
    return Json::parse(value.dump());
}
Matrix read_matrix(const fs::path &root, const Json &spec) {
    const auto format = spec.at("format").text();
    if (format == "h5ad")
        return read_h5ad(root, spec);
    require(format == "matrix_market" || format == "dense_tsv", "unsupported assay matrix format");
    Matrix out;
    out.features = read_axis(local_artifact(root, spec.at("features")));
    out.observations = read_axis(local_artifact(root, spec.at("observations")));
    out.bytes = local_artifact(root, spec.at("source"));
    const auto nr = out.features.ids.size(), nc = out.observations.ids.size();
    (void)rfc::numeric::IndexDomain({nr, nc});
    const auto arithmetic = spec.at("arithmetic").text();
    require(arithmetic == "exact_decimal" || arithmetic == "exact_rational", "matrix arithmetic");
    const bool counts = is_count_quantity(spec.at("quantity").text());
    const auto value = [&](const std::string &s) {
        auto q = arithmetic == "exact_decimal" ? exact_decimal(s) : exact_rational(s);
        if (counts)
            require(q >= 0 && q.get_den() == 1, "raw counts require nonnegative integer values");
        return q;
    };
    std::vector<std::vector<rfc::numeric::SparseRationalEntry>> entries(nc);
    std::vector<std::vector<std::size_t>> reported(nc);
    std::vector<std::set<std::size_t>> missing(nc);
    MatrixRecords records(out.bytes);
    std::string_view record;
    const auto implicit = spec.at("implicit_entries").text();
    require(implicit == "measured_zero" || implicit == "unobserved",
            "matrix implicit-entry policy");
    std::set<std::string> tokens;
    for (std::size_t i = 0; i < spec.at("missing_tokens").size(); ++i)
        require(tokens.insert(spec.at("missing_tokens").at(i).text()).second,
                "duplicate missing token");
    if (format == "dense_tsv") {
        require(!spec.has("missing"),
                "dense matrices express missingness through their declared tokens");
        require(records.next(record), "dense matrix header missing");
        const auto header = split(std::string(record), '\t');
        require(header.size() == nc + 1 && header[0] == "id", "dense matrix header");
        for (std::size_t c = 0; c < nc; ++c)
            require(header[c + 1] == out.observations.ids[c],
                    "dense observation order differs from axis");
        for (std::size_t r = 0; r < nr; ++r) {
            require(records.next(record), "dense matrix row count");
            const auto row = split(std::string(record), '\t');
            require(row.size() == nc + 1 && row[0] == out.features.ids[r],
                    "dense feature order/shape");
            for (std::size_t c = 0; c < nc; ++c) {
                if (tokens.contains(row[c + 1])) {
                    missing[c].insert(r);
                    continue;
                }
                auto v = value(row[c + 1]);
                if (v != 0) {
                    entries[c].push_back({r, std::move(v)});
                }
            }
        }
        require(!records.next(record), "dense matrix row count");
    } else {
        require(tokens.empty(),
                "Matrix Market missingness requires an explicit missing-coordinate mask");
        require(records.next(record), "empty Matrix Market input");
        const auto banner = record;
        require(banner == "%%MatrixMarket matrix coordinate integer general" ||
                    banner == "%%MatrixMarket matrix coordinate real general",
                "Matrix Market requires coordinate integer/real general encoding");
        const bool integer = banner.find("integer") != std::string_view::npos;
        do {
            require(records.next(record), "Matrix Market dimensions missing");
        } while (record.front() == '%');
        const auto shape = matrix_tokens(record);
        const auto nnz = matrix_unsigned(shape[2]);
        require(matrix_unsigned(shape[0]) == nr && matrix_unsigned(shape[1]) == nc,
                "Matrix Market shape differs from axes");
        for (std::uint64_t n = 0; n < nnz; ++n) {
            require(records.next(record), "Matrix Market entry count");
            const auto f = matrix_tokens(record);
            const auto ri = matrix_unsigned(f[0]), ci = matrix_unsigned(f[1]);
            require(ri && ri <= nr && ci && ci <= nc, "Matrix Market one-based coordinate range");
            const auto r = ri - 1, c = ci - 1;
            reported[c].push_back(r);
            auto v = value(std::string(f[2]));
            require(!integer || v.get_den() == 1, "integer matrix token");
            if (v != 0) {
                entries[c].push_back({r, std::move(v)});
            }
        }
        require(!records.next(record), "Matrix Market entry count");
        // Explicit zeros also participate in duplicate and missing-mask checks.
        // Sorting evidence indices does not change the retained source bytes.
        for (auto &column : reported) {
            std::sort(column.begin(), column.end());
            require(std::adjacent_find(column.begin(), column.end()) == column.end(),
                    "duplicate matrix coordinate; aggregation must be declared upstream");
        }
    }
    if (spec.has("missing")) {
        const auto mask = local_artifact(root, spec.at("missing"));
        for (const auto [p, q] : lines(mask)) {
            const auto f = split(mask.substr(p, q - p), '\t');
            require(f.size() == 2, "missing-mask columns");
            const auto r = unsigned_token(f[0]), c = unsigned_token(f[1]);
            require(r && r <= nr && c && c <= nc, "missing-mask one-based coordinates");
            require(!std::binary_search(reported[c - 1].begin(), reported[c - 1].end(), r - 1) &&
                        missing[c - 1].insert(r - 1).second,
                    "missing mask overlaps a reported value or repeats a coordinate");
        }
    }
    for (std::size_t c = 0; c < nc; ++c) {
        std::vector<std::size_t> support;
        const bool full =
            (format == "dense_tsv" || implicit == "measured_zero") && missing[c].empty();
        if (!full && format == "matrix_market" && implicit == "unobserved")
            support = std::move(reported[c]);
        else if (!full)
            for (std::size_t r = 0; r < nr; ++r)
                if (!missing[c].contains(r) &&
                    (format == "dense_tsv" || implicit == "measured_zero" ||
                     std::binary_search(reported[c].begin(), reported[c].end(), r)))
                    support.push_back(r);
        out.columns.emplace_back(rfc::numeric::IndexDomain::linear(nr),
                                 full ? rfc::numeric::StructuralSupport::full(nr)
                                      : rfc::numeric::StructuralSupport::sparse(nr, std::move(support)),
                                 rfc::numeric::SparseRationalVector(nr, std::move(entries[c])));
        std::vector<std::size_t>().swap(reported[c]);
    }
    return out;
}
std::vector<Interval> read_intervals(const std::string &bytes, const Json &resource) {
    const auto format = resource.at("format").text();
    require(format == "bed" || format == "gtf" || format == "gff3", "interval resource format");
    const auto reference = reference_record(resource.at("reference"));
    std::map<std::string, std::uint64_t> contigs;
    for (std::size_t i = 0; i < reference.at("contigs").size(); ++i) {
        const auto c = reference.at("contigs").at(i);
        contigs.emplace(c.at("id").text(), c.at("length").integer());
    }
    std::vector<Interval> out;
    for (const auto [p, q] : lines(bytes)) {
        const auto line = bytes.substr(p, q - p);
        require(line != "##FASTA", "GFF3 embedded FASTA requires a separate reference resource");
        if (line.starts_with("#") || line.starts_with("track ") || line.starts_with("browser "))
            continue;
        const auto f = split(line, '\t');
        Interval record;
        record.source_start = p;
        record.source_end = q;
        if (format == "bed") {
            require(f.size() >= 3 && f.size() <= 6,
                    "BED profile admits 3..6 columns; blocked BED needs an explicit adapter");
            record.contig = f[0];
            record.start = unsigned_token(f[1]);
            record.end = unsigned_token(f[2]);
            record.id = f.size() > 3 ? f[3] : "record:" + std::to_string(out.size() + 1);
            record.strand = f.size() > 5 ? f[5] : ".";
            if (resource.has("weight_column")) {
                const auto c = count(resource.at("weight_column"), 5, "BED weight column");
                require(c >= 3 && c < f.size(), "BED weight column is absent or is a coordinate");
                record.weight = exact_rational(f[c]);
            }
        } else {
            require(format == "gtf" || format == "gff3", "annotation format");
            require(f.size() == 9, "annotation requires nine columns");
            record.contig = f[0];
            const auto first = unsigned_token(f[3]);
            require(first > 0, "GTF coordinates are one-based closed");
            record.start = first - 1;
            record.end = unsigned_token(f[4]);
            record.strand = f[6];
            record.feature_type = f[2];
            record.id = "record:" + std::to_string(out.size() + 1);
            require(f[7] == "." || f[7] == "0" || f[7] == "1" || f[7] == "2", "GTF frame");
            if (format == "gtf" && f[8] != ".") {
                std::size_t p = 0;
                while (p < f[8].size()) {
                    while (p < f[8].size() && f[8][p] == ' ')
                        ++p;
                    if (p == f[8].size())
                        break;
                    const auto end = f[8].find(' ', p);
                    require(end != std::string::npos && end > p, "GTF attribute key");
                    const auto key = f[8].substr(p, end - p);
                    p = end + 1;
                    while (p < f[8].size() && f[8][p] == ' ')
                        ++p;
                    require(p < f[8].size() && f[8][p++] == '"', "GTF attribute must be quoted");
                    const auto close = f[8].find('"', p);
                    require(close != std::string::npos, "GTF attribute closing quote");
                    const auto value = f[8].substr(p, close - p);
                    require(value.find('\\') == std::string::npos,
                            "unsupported escaped GTF attribute");
                    p = close + 1;
                    require(p < f[8].size() && f[8][p++] == ';', "GTF attribute semicolon");
                    record.attributes[key].push_back(value);
                }
            } else if (f[8] != ".") {
                const auto decode = [](const std::string &token) {
                    std::string result;
                    for (std::size_t i = 0; i < token.size(); ++i) {
                        if (token[i] != '%') {
                            result += token[i];
                            continue;
                        }
                        require(i + 2 < token.size(), "GFF3 percent escape");
                        unsigned n{};
                        const auto [end, ec] =
                            std::from_chars(token.data() + i + 1, token.data() + i + 3, n, 16);
                        require(ec == std::errc{} && end == token.data() + i + 3 && n != 0,
                                "GFF3 percent escape");
                        result += static_cast<char>(n);
                        i += 2;
                    }
                    return result;
                };
                for (const auto &attribute : split(f[8], ';')) {
                    const auto eq = attribute.find('=');
                    require(eq != std::string::npos && eq > 0, "GFF3 attribute key/value");
                    const auto key = decode(attribute.substr(0, eq));
                    require(!record.attributes.contains(key), "duplicate GFF3 attribute key");
                    for (const auto &v : split(attribute.substr(eq + 1), ','))
                        record.attributes[key].push_back(decode(v));
                }
                if (record.attributes.contains("ID")) {
                    require(record.attributes.at("ID").size() == 1, "GFF3 ID multiplicity");
                    record.id = record.attributes.at("ID").front();
                }
            }
        }
        require(!record.id.empty() &&
                    (record.strand == "+" || record.strand == "-" || record.strand == "."),
                "interval identity/strand");
        require(contigs.contains(record.contig) && record.start < record.end &&
                    record.end <= contigs.at(record.contig),
                "interval outside declared reference; zero-length intervals are not coverage");
        out.push_back(record);
    }
    require(!out.empty(), "interval resource contains no observations");
    return out;
}
Json inspect_resource(const std::string &bytes, const Json &resource,
                      const std::string &reference_bytes) {
    auto out = Json::object();
    out.set("format", resource.at("format"));
    out.set("source_sha256", Json::string("sha256:" + sha256(bytes)));
    const auto format = resource.at("format").text();
    if (format == "sam" || format == "bam" || format == "cram" || format == "vcf" ||
        format == "bcf")
        return inspect_hts({resource, bytes, reference_bytes});
    if (format == "bed" || format == "gtf" || format == "gff3") {
        const auto intervals = read_intervals(bytes, resource);
        auto records = Json::array();
        for (const auto &r : intervals) {
            auto item = Json::object();
            item.set("id", Json::string(r.id));
            item.set("contig", Json::string(r.contig));
            item.set("start", cardinality(r.start));
            item.set("end", cardinality(r.end));
            item.set("strand", Json::string(r.strand));
            item.set("weight", Json::string(r.weight.get_str()));
            item.set("source_start", cardinality(r.source_start));
            item.set("source_end", cardinality(r.source_end));
            if (!r.feature_type.empty()) {
                item.set("feature_type", Json::string(r.feature_type));
                auto attributes = Json::object();
                for (const auto &[key, values] : r.attributes) {
                    auto entries = Json::array();
                    for (const auto &v : values)
                        entries.push(Json::string(v));
                    attributes.set(key, entries);
                }
                item.set("attributes", attributes);
            }
            records.push(item);
        }
        out.set("reference", resource.at("reference"));
        out.set("records", records);
        out.set("record_count", cardinality(intervals.size()));
    } else if (format == "fasta" || format == "fastq") {
        const auto ls = lines(bytes);
        auto records = Json::array();
        std::set<std::string> ids;
        const auto line = [&](std::size_t i) {
            return bytes.substr(ls.at(i).first, ls.at(i).second - ls.at(i).first);
        };
        const auto sequence = [&](const std::string &s) {
            require(!s.empty(), "empty sequence");
            for (auto c : s)
                require(std::string("ACGTURYSWKMBDHVNacgturyswkmbdhvn").find(c) !=
                            std::string::npos,
                        "invalid IUPAC base");
        };
        for (std::size_t i = 0; i < ls.size();) {
            const auto start = i;
            const auto header = line(i++);
            require(header.size() > 1 && header[0] == (format == "fasta" ? '>' : '@'),
                    "sequence header");
            const auto id = header.substr(1, header.find_first_of(" \t", 1) - 1);
            require(!id.empty() && ids.insert(id).second,
                    "sequence identities must be nonempty and unique");
            std::string bases;
            if (format == "fasta") {
                while (i < ls.size() && line(i)[0] != '>') {
                    const auto s = line(i++);
                    sequence(s);
                    bases += s;
                }
                require(!bases.empty(), "FASTA record has no sequence");
            } else {
                require(i + 2 < ls.size(), "truncated four-line FASTQ record");
                bases = line(i++);
                sequence(bases);
                const auto plus = line(i++), qualities = line(i++);
                require(plus.starts_with("+"), "FASTQ separator");
                require(plus == "+" || plus.substr(1) == header.substr(1),
                        "FASTQ repeated header differs");
                require(qualities.size() == bases.size(), "FASTQ quality/sequence length");
                for (auto c : qualities)
                    require(c >= 33 && c <= 126,
                            "FASTQ quality encoding requires printable Phred+33 tokens");
            }
            auto record = Json::object();
            record.set("id", Json::string(id));
            record.set("length", cardinality(bases.size()));
            record.set("source_start", cardinality(ls[start].first));
            record.set("source_end", cardinality(ls[i - 1].second));
            records.push(record);
        }
        out.set("records", records);
        out.set("record_count", cardinality(records.size()));
        out.set("quality_interpretation",
                Json::string(format == "fastq"
                                 ? "Phred+33 encoding; scores are retained, not recalibrated"
                                 : "not applicable"));
    } else if (format == "salmon_eq") {
        const auto parsed = salmon_classes(bytes);
        json_object_object_foreach(parsed.raw(), name, value) {
            (void)value;
            out.set(name, parsed.at(name));
        }
    } else
        throw std::invalid_argument("unsupported physical file adapter: " + format);
    return out;
}
} // namespace gadget::assay_detail

namespace gadget {
Json inspect_bioinformatics_file(const fs::path &specification) {
    const auto path = fs::canonical(specification);
    const auto spec = Json::parse(read_text(path));
    assay_detail::members(spec, {"schema", "source", "format"},
                          {"reference", "weight_column", "reference_source"});
    assay_detail::require(spec.at("schema").text() == "gadget.bioinformatics.file.v1",
                          "physical file schema");
    auto result = assay_detail::inspect_resource(
        assay_detail::local_artifact(path.parent_path(), spec.at("source")), spec,
        spec.has("reference_source")
            ? assay_detail::local_artifact(path.parent_path(), spec.at("reference_source"))
            : std::string{});
    result.set("schema", Json::string("gadget.bioinformatics.file.result.v1"));
    result.set("request", spec);
    result.set("source_digest", Json::string("sha256:" + sha256(spec.dump())));
    return assay_detail::seal(result, "certificate_sha256");
}
} // namespace gadget
