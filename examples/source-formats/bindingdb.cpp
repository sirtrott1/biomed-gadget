#include "bindingdb.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace gadget::examples {
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::invalid_argument(message);
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested()) throw std::runtime_error("BindingDB bulk read cancelled");
}
struct Line {
    std::string_view content, original;
    std::size_t next;
};
Line line_at(std::string_view source, std::size_t offset) {
    const auto newline = source.find('\n', offset);
    const auto end = newline == source.npos ? source.size() : newline;
    auto content_end = end;
    if (content_end > offset && source[content_end - 1] == '\r') --content_end;
    const auto next = newline == source.npos ? end : end + 1;
    return {source.substr(offset, content_end - offset),
            source.substr(offset, next - offset), next};
}
std::vector<std::string_view> fields(std::string_view content,
                                    const BindingDBBulkLimits &limits) {
    std::vector<std::string_view> values;
    std::size_t offset = 0;
    while (true) {
        const auto tab = content.find('\t', offset);
        const auto end = tab == content.npos ? content.size() : tab;
        require(values.size() < limits.columns, "BindingDB column budget exceeded");
        const auto value = content.substr(offset, end - offset);
        require(value.size() <= limits.field_bytes, "BindingDB field byte budget exceeded");
        require(value.find('\0') == value.npos && value.find('\r') == value.npos,
                "BindingDB field contains an unsupported control byte");
        values.push_back(value);
        if (tab == content.npos) return values;
        offset = tab + 1;
    }
}
std::optional<std::size_t> column(const BindingDBBulkHeader &header,
                                std::string_view name, bool required = false) {
    std::optional<std::size_t> found;
    for (std::size_t i = 0; i < header.columns.size(); ++i) {
        if (header.columns[i] != name) continue;
        require(!found, "BindingDB semantic column is duplicated");
        found = i;
    }
    require(!required || found.has_value(), "BindingDB required column is absent");
    return found;
}
std::optional<std::size_t> chain_start(const BindingDBBulkHeader &header) {
    constexpr std::array<std::string_view, 12> names{
        "BindingDB Target Chain Sequence", "PDB ID(s) of Target Chain",
        "UniProt (SwissProt) Recommended Name of Target Chain",
        "UniProt (SwissProt) Entry Name of Target Chain",
        "UniProt (SwissProt) Primary ID of Target Chain",
        "UniProt (SwissProt) Secondary ID(s) of Target Chain",
        "UniProt (SwissProt) Alternative ID(s) of Target Chain",
        "UniProt (TrEMBL) Submitted Name of Target Chain",
        "UniProt (TrEMBL) Entry Name of Target Chain",
        "UniProt (TrEMBL) Primary ID of Target Chain",
        "UniProt (TrEMBL) Secondary ID(s) of Target Chain",
        "UniProt (TrEMBL) Alternative ID(s) of Target Chain"};
    std::size_t first = 0;
    while (first < header.columns.size() &&
           !header.columns[first].starts_with(names[0])) ++first;
    if (first == header.columns.size() || (header.columns.size() - first) % 12) return {};
    for (std::size_t i = first; i < header.columns.size(); ++i) {
        const auto base = names[(i - first) % 12];
        const auto numbered = std::string(base) + " " + std::to_string((i - first) / 12 + 1);
        if (header.columns[i] != base && header.columns[i] != numbered) return {};
    }
    return first;
}
std::string_view trim(std::string_view value) {
    while (!value.empty() && value.front() == ' ') value.remove_prefix(1);
    while (!value.empty() && value.back() == ' ') value.remove_suffix(1);
    return value;
}
bool digit(char value) { return value >= '0' && value <= '9'; }
bool decimal(std::string_view value) {
    std::size_t i = 0, digits = 0;
    if (i < value.size() && (value[i] == '+' || value[i] == '-')) ++i;
    while (i < value.size() && digit(value[i])) { ++i; ++digits; }
    if (i < value.size() && value[i] == '.') {
        ++i;
        while (i < value.size() && digit(value[i])) { ++i; ++digits; }
    }
    if (!digits) return false;
    if (i < value.size() && (value[i] == 'e' || value[i] == 'E')) {
        ++i;
        if (i < value.size() && (value[i] == '+' || value[i] == '-')) ++i;
        const auto start = i;
        while (i < value.size() && digit(value[i])) ++i;
        if (i == start) return false;
    }
    return i == value.size();
}
BindingDBAffinityView affinity(std::size_t index, std::string_view endpoint,
                              std::string_view reported) {
    auto value = trim(reported);
    BindingDBAffinityView result{index, endpoint, reported, {}, {}};
    if (value.empty()) return result;
    result.relation = BindingDBRelation::Equal;
    for (const auto &[prefix, relation] :
         std::array<std::pair<std::string_view, BindingDBRelation>, 6>{{
             {"<=", BindingDBRelation::LessEqual}, {">=", BindingDBRelation::GreaterEqual},
             {"<", BindingDBRelation::Less}, {">", BindingDBRelation::Greater},
             {"=", BindingDBRelation::Equal}, {"~", BindingDBRelation::Approximate}}}) {
        if (!value.starts_with(prefix)) continue;
        result.relation = relation;
        value = trim(value.substr(prefix.size()));
        break;
    }
    require(decimal(value), "BindingDB affinity requires a complete decimal token");
    result.number_token = value;
    return result;
}
} // namespace

BindingDBBulkSummary visit_bindingdb_bulk(
    std::string_view source, const BindingDBBulkLimits &limits,
    const std::function<void(const BindingDBBulkRowView &)> &consume,
    std::stop_token stop) {
    require(limits.source_bytes && limits.rows && limits.columns && limits.field_bytes,
            "BindingDB budgets must be positive");
    require(source.size() <= limits.source_bytes, "BindingDB source byte budget exceeded");
    require(!source.empty() && static_cast<bool>(consume), "BindingDB source or consumer is empty");
    cancelled(stop);
    const auto first = line_at(source, 0);
    BindingDBBulkHeader header{fields(first.content, limits), {}};
    for (const auto name : header.columns) require(!name.empty(), "BindingDB header name is empty");
    const auto reaction = column(header, "BindingDB Reactant_set_id", true).value();
    const auto monomer = column(header, "BindingDB MonomerID", true).value();
    (void)column(header, "Target Name", true);
    (void)column(header, "Target Source Organism According to Curator or DataSource", true);
    for (const auto name : {"pH", "Temp (C)", "Curation/DataSource", "Article DOI",
                           "BindingDB Entry DOI", "PMID", "PubChem AID", "PubChem SID",
                           "PubChem CID", "Number of Protein Chains in Target (>1 implies a multichain complex)"})
        (void)column(header, name);
    header.chain_start = chain_start(header);
    constexpr std::array<std::string_view, 4> endpoints{"Ki", "Kd", "IC50", "EC50"};
    std::vector<std::pair<std::size_t, std::string_view>> affinity_columns;
    for (const auto endpoint : endpoints) {
        const auto index = column(header, std::string(endpoint) + " (nM)");
        if (index) affinity_columns.emplace_back(*index, endpoint);
    }
    require(!affinity_columns.empty(), "BindingDB affinity columns are absent");
    std::sort(affinity_columns.begin(), affinity_columns.end());
    BindingDBBulkSummary result;
    for (auto offset = first.next; offset < source.size();) {
        cancelled(stop);
        require(result.rows < limits.rows, "BindingDB row budget exceeded");
        const auto line = line_at(source, offset);
        require(!line.content.empty(), "BindingDB blank row is not an observation");
        const auto cells = fields(line.content, limits);
        if (cells.size() != header.columns.size()) {
            require(header.chain_start && cells.size() >= *header.chain_start + 12 &&
                        (cells.size() - *header.chain_start) % 12 == 0,
                    "BindingDB row has an incomplete or undeclared column layout");
            require(reaction < cells.size() && monomer < cells.size(),
                    "BindingDB row lacks its declared identity columns");
        }
        require(!cells[reaction].empty() && !cells[monomer].empty(),
                "BindingDB reaction or monomer identity is empty");
        std::vector<BindingDBAffinityView> values;
        for (const auto &[index, endpoint] : affinity_columns) {
            require(index < cells.size(), "BindingDB row lacks its declared affinity column");
            values.push_back(affinity(index, endpoint, cells[index]));
        }
        std::size_t reported = 0;
        for (const auto &value : values) reported += value.relation.has_value() ? 1 : 0;
        require(reported <= std::numeric_limits<std::size_t>::max() - result.reported_affinities,
                "BindingDB affinity count overflow");
        require(result.rows <= std::numeric_limits<std::size_t>::max() - 2,
                "BindingDB line number overflow");
        const BindingDBBulkRowView row{header, result.rows, result.rows + 2, offset,
            line.original.size(), line.original, cells, values};
        consume(row);
        cancelled(stop);
        ++result.rows;
        result.reported_affinities += reported;
        offset = line.next;
    }
    return result;
}
} // namespace gadget::examples
