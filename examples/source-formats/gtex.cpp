#include "gtex.hpp"
#include <algorithm>
#include <charconv>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

namespace gadget::examples {
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
void check_stop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("GTEx source visit cancelled");
}
class Records {
  public:
    Records(std::string_view bytes, const GtexLimits &limits, std::stop_token stop)
        : bytes_(bytes), limits_(limits), stop_(stop) {
        require(limits.source_bytes && limits.rows && limits.columns && limits.record_bytes &&
                    limits.field_bytes && limits.cells,
                "GTEx source limits must be positive");
        require(bytes.size() <= limits.source_bytes, "GTEx source byte limit");
        check_stop(stop_);
        require(bytes.find('\0') == std::string_view::npos, "NUL in GTEx source");
        check_stop(stop_);
    }
    bool next(GtexSourceRecord &source, std::string_view &content) {
        check_stop(stop_);
        if (offset_ == bytes_.size())
            return false;
        const auto start = offset_;
        auto end = bytes_.find('\n', start);
        const auto terminated = end != std::string_view::npos;
        if (!terminated)
            end = bytes_.size();
        offset_ = end + static_cast<std::size_t>(terminated);
        require(offset_ - start <= limits_.record_bytes, "GTEx record byte limit");
        source = {++line_, start, offset_ - start, bytes_.substr(start, offset_ - start)};
        if (end > start && bytes_[end - 1] == '\r')
            --end;
        require(end > start, "empty GTEx record");
        content = bytes_.substr(start, end - start);
        require(content.find('\r') == std::string_view::npos, "embedded CR in GTEx source");
        return true;
    }

  private:
    std::string_view bytes_;
    const GtexLimits &limits_;
    std::stop_token stop_;
    std::size_t offset_{}, line_{};
};
std::vector<std::string_view> fields(std::string_view line, const GtexLimits &limits,
                                     std::size_t maximum) {
    std::vector<std::string_view> out;
    std::size_t position{};
    for (;;) {
        require(out.size() < maximum, "GTEx column limit or row width");
        auto end = line.find('\t', position);
        if (end == std::string_view::npos)
            end = line.size();
        require(end - position <= limits.field_bytes, "GTEx field byte limit");
        out.push_back(line.substr(position, end - position));
        if (end == line.size())
            return out;
        position = end + 1;
    }
}
std::size_t natural(std::string_view token) {
    std::size_t value{};
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
    require(!token.empty() && error == std::errc{} && end == token.data() + token.size(),
            "GTEx dimension or ordinal must be an unsigned integer");
    return value;
}
// Validate exact nonnegative decimal syntax and integer support without a
// floating conversion, rounding or a bounded integer representation of counts.
struct Decimal {
    bool zero{}, integer{};
};
Decimal decimal(std::string_view token) {
    require(!token.empty(), "empty GTEx value");
    std::size_t at = token.front() == '+' ? 1 : 0;
    std::size_t digits{}, fractional{}, trailing_zeros{};
    bool dot{}, nonzero{};
    for (; at < token.size() && token[at] != 'e' && token[at] != 'E'; ++at) {
        const auto ch = token[at];
        if (ch == '.' && !dot) {
            dot = true;
            continue;
        }
        require(ch >= '0' && ch <= '9', "GTEx value is not a nonnegative decimal");
        ++digits;
        fractional += static_cast<std::size_t>(dot);
        nonzero = nonzero || ch != '0';
        trailing_zeros = ch == '0' ? trailing_zeros + 1 : 0;
    }
    require(digits > 0, "GTEx decimal has no digits");
    bool negative_exponent{};
    std::size_t exponent{};
    if (at < token.size()) {
        ++at;
        if (at < token.size() && (token[at] == '+' || token[at] == '-'))
            negative_exponent = token[at++] == '-';
        require(at < token.size(), "GTEx decimal exponent has no digits");
        // Only the exponent's comparison to the significand length matters.
        // Saturate at token.size()+1 rather than overflowing on long exponents.
        for (; at < token.size(); ++at) {
            const auto ch = token[at];
            require(ch >= '0' && ch <= '9', "GTEx decimal exponent syntax");
            if (exponent <= token.size()) {
                const auto digit = static_cast<std::size_t>(ch - '0');
                if (exponent > (token.size() - std::min(token.size(), digit)) / 10)
                    exponent = token.size() + 1;
                else
                    exponent = exponent * 10 + digit;
            }
        }
    }
    bool integer = !nonzero;
    if (nonzero && negative_exponent)
        integer = fractional <= trailing_zeros && exponent <= trailing_zeros - fractional;
    else if (nonzero)
        integer = fractional <= trailing_zeros || exponent >= fractional - trailing_zeros;
    return {!nonzero, integer};
}
bool versioned_ensembl(std::string_view id, GtexFeature feature) {
    const auto prefix =
        feature == GtexFeature::Gene ? std::string_view("ENSG") : std::string_view("ENST");
    if (!id.starts_with(prefix))
        return false;
    id.remove_prefix(prefix.size());
    // GENCODE v26 retains Y pseudoautosomal copies as distinct source IDs.
    // Validate the documented suffix without stripping it from the row view.
    if (id.ends_with("_PAR_Y"))
        id.remove_suffix(6);
    const auto dot = id.find('.');
    if (dot == std::string_view::npos || !dot || dot + 1 == id.size())
        return false;
    for (std::size_t i = 0; i < id.size(); ++i)
        if (i != dot && (id[i] < '0' || id[i] > '9'))
            return false;
    return true;
}
std::size_t column(const GtexMetadataTable &table, std::string_view name) {
    const auto found = std::find(table.columns.begin(), table.columns.end(), name);
    require(found != table.columns.end(), "GTEx required metadata column absent");
    return static_cast<std::size_t>(found - table.columns.begin());
}
} // namespace
GtexSummary visit_gtex_gct(std::string_view source, GtexDeclaration declaration,
                           const GtexLimits &limits,
                           const std::function<void(const GtexRowView &)> &visitor,
                           std::stop_token stop) {
    require(static_cast<bool>(visitor), "GTEx source visitor is empty");
    require(declaration.layout == GtexLayout::Gct12 ||
                declaration.layout == GtexLayout::PortalTissue13,
            "GTEx layout declaration");
    require(declaration.feature == GtexFeature::Gene ||
                declaration.feature == GtexFeature::Transcript,
            "GTEx feature declaration");
    require(declaration.measurement == GtexMeasurement::IntegerCounts ||
                declaration.measurement == GtexMeasurement::EstimatedCounts ||
                declaration.measurement == GtexMeasurement::Tpm,
            "GTEx measurement declaration");
    require((declaration.feature == GtexFeature::Gene &&
             declaration.measurement != GtexMeasurement::EstimatedCounts) ||
                (declaration.feature == GtexFeature::Transcript &&
                 declaration.measurement != GtexMeasurement::IntegerCounts),
            "GTEx gene reads and estimated transcript counts are distinct source families");
    require(declaration.layout != GtexLayout::PortalTissue13 ||
                declaration.feature == GtexFeature::Gene,
            "GTEx portal tissue profile supports gene rows only");
    Records records(source, limits, stop);
    GtexHeader header;
    header.declaration = declaration;
    GtexSourceRecord source_record;
    std::string_view record;
    require(records.next(header.version_record, record), "GTEx GCT version absent");
    const bool portal = declaration.layout == GtexLayout::PortalTissue13;
    require(record == (portal ? "#1.3" : "#1.2"), "GTEx GCT version differs from declared layout");
    require(records.next(header.dimensions_record, record), "GTEx GCT dimensions absent");
    const auto shape = fields(record, limits, 4);
    require(shape.size() == (portal ? 4u : 2u), "GTEx GCT dimensions width");
    header.features = natural(shape[0]);
    header.declared_columns = natural(shape[1]);
    require(header.features && header.features <= limits.rows, "GTEx feature dimension limit");
    require(!portal || (shape[2] == "0" && shape[3] == "0" && header.declared_columns > 2),
            "GTEx portal tissue export annotation convention");
    const auto samples = portal ? header.declared_columns - 2 : header.declared_columns;
    require(samples && samples <= limits.columns, "GTEx sample dimension limit");
    require(samples <= limits.cells / header.features, "GTEx matrix cell limit");
    require(samples <= std::numeric_limits<std::size_t>::max() - 3, "GTEx row width overflow");
    const auto prefix = portal ? 3u : 2u;
    require(records.next(header.labels_record, record), "GTEx GCT labels absent");
    const auto labels = fields(record, limits, samples + prefix);
    require(labels.size() == samples + prefix, "GTEx GCT header width");
    if (portal)
        require(labels[0] == "id" && labels[1] == "Name" && labels[2] == "Description",
                "GTEx portal tissue header");
    else if (declaration.feature == GtexFeature::Transcript)
        require(labels[0] == "transcript_id" && labels[1] == "gene_id", "GTEx transcript header");
    else
        require((labels[0] == "Name" || labels[0] == "name") &&
                    (labels[1] == "Description" || labels[1] == "description"),
                "GTEx GCT gene header");
    std::set<std::string_view> seen_samples;
    for (std::size_t i = prefix; i < labels.size(); ++i) {
        require(!labels[i].empty() && seen_samples.insert(labels[i]).second,
                "GTEx sample identities must be nonempty and unique");
        header.samples.push_back(labels[i]);
    }
    GtexSummary summary;
    summary.samples = samples;
    std::set<std::string_view> seen_features;
    std::vector<GtexValue> values(samples);
    for (std::size_t row = 0; row < header.features; ++row) {
        require(records.next(source_record, record), "GTEx GCT feature count");
        const auto cells = fields(record, limits, samples + prefix);
        require(cells.size() == samples + prefix, "GTEx GCT value width");
        require(!portal || natural(cells[0]) == row, "GTEx portal row ordinal differs");
        const auto feature = cells[prefix - 2], description = cells[prefix - 1];
        require(versioned_ensembl(feature, declaration.feature) &&
                    seen_features.insert(feature).second,
                "GTEx feature must be a unique versioned Ensembl identity");
        require(declaration.feature != GtexFeature::Transcript ||
                    versioned_ensembl(description, GtexFeature::Gene),
                "GTEx transcript gene link must retain its versioned Ensembl identity");
        for (std::size_t c = 0; c < samples; ++c) {
            const auto token = cells[prefix + c];
            const auto observed = token != "NA";
            values[c] = {token, observed};
            if (!observed) {
                ++summary.missing;
                continue;
            }
            const auto number = decimal(token);
            require(declaration.measurement != GtexMeasurement::IntegerCounts || number.integer,
                    "GTEx integer counts require nonnegative integral values");
            ++summary.observed;
            summary.observed_zeros += static_cast<std::size_t>(number.zero);
        }
        visitor({header, row, source_record, feature, description, values});
        ++summary.rows;
    }
    require(!records.next(source_record, record), "extra GTEx GCT feature row");
    check_stop(stop);
    return summary;
}
GtexMetadataTable read_gtex_metadata(std::string_view source, GtexMetadataKind kind,
                                     const GtexLimits &limits, std::stop_token stop) {
    require(kind == GtexMetadataKind::SampleAttributes ||
                kind == GtexMetadataKind::SubjectPhenotypes,
            "GTEx metadata kind");
    Records records(source, limits, stop);
    GtexMetadataTable table{kind, {}, {}, {}};
    GtexSourceRecord source_record;
    std::string_view record;
    require(records.next(table.header_record, record), "GTEx metadata header absent");
    table.columns = fields(record, limits, limits.columns);
    require(table.columns.front() ==
                (kind == GtexMetadataKind::SampleAttributes ? "SAMPID" : "SUBJID"),
            "GTEx metadata identifier column differs");
    std::set<std::string_view> names;
    for (const auto name : table.columns)
        require(!name.empty() && names.insert(name).second, "GTEx metadata columns must be unique");
    while (records.next(source_record, record)) {
        require(table.rows.size() < limits.rows, "GTEx metadata row limit");
        require(table.columns.size() <= limits.cells / (table.rows.size() + 1),
                "GTEx metadata cell limit");
        auto cells = fields(record, limits, table.columns.size());
        require(cells.size() == table.columns.size(), "GTEx metadata row width");
        require(!cells.front().empty() && cells.front() != "NA", "GTEx metadata identity absent");
        const auto id = cells.front();
        require(table.rows.emplace(id, GtexMetadataRow{source_record, std::move(cells)}).second,
                "duplicate GTEx metadata identity");
    }
    check_stop(stop);
    return table;
}
GtexSampleSet join_gtex_samples(std::span<const std::string_view> samples,
                                const GtexMetadataTable &sample_attributes,
                                const GtexMetadataTable &subject_phenotypes, std::stop_token stop) {
    require(sample_attributes.kind == GtexMetadataKind::SampleAttributes &&
                subject_phenotypes.kind == GtexMetadataKind::SubjectPhenotypes,
            "GTEx metadata tables have different roles");
    const auto tissue_column = column(sample_attributes, "SMTSD");
    // The tables remain separate: no same named sample/subject field is silently
    // overwritten, and source categorical metadata is never coerced to numbers.
    GtexSampleSet out;
    std::set<std::string_view> seen;
    std::map<std::string_view, std::size_t> donors;
    for (const auto id : samples) {
        check_stop(stop);
        require(seen.insert(id).second, "duplicate GTEx selected sample");
        require(id.starts_with("GTEX-"), "GTEx sample identifier prefix");
        const auto end = id.find('-', 5);
        require(end != std::string_view::npos && end > 5 && end + 1 < id.size(),
                "GTEx sample cannot resolve its donor");
        const auto donor = id.substr(0, end);
        const auto sample = sample_attributes.rows.find(id);
        const auto subject = subject_phenotypes.rows.find(donor);
        require(sample != sample_attributes.rows.end() && subject != subject_phenotypes.rows.end(),
                "GTEx sample or donor absent from source metadata");
        const auto tissue = sample->second.fields[tissue_column];
        require(!tissue.empty() && tissue != "NA", "GTEx tissue detail absent");
        out.samples.push_back({id, donor, tissue, &sample->second, &subject->second, false});
        ++donors[donor];
    }
    out.distinct_donors = donors.size();
    for (const auto &[id, count] : donors) {
        (void)id;
        out.repeated_donors += static_cast<std::size_t>(count > 1);
    }
    for (auto &sample : out.samples)
        sample.repeated_donor = donors.at(sample.donor) > 1;
    check_stop(stop);
    return out;
}
void GtexSampleSet::require_independent_donors() const {
    require(!repeated_donors, "repeated GTEx donor cannot be an independent replicate");
}
} // namespace gadget::examples
