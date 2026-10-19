#include "../examples/source-formats/bindingdb.hpp"
#include <json-c/json.h>
#include <array>
#include <gadget/io.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

using namespace gadget::examples;
namespace {
constexpr BindingDBBulkLimits limits{16 * 1024 * 1024, 100000, 4096, 1024 * 1024};
std::size_t checks{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F callback) {
    bool rejected = false;
    try { callback(); } catch (const std::exception &) { rejected = true; }
    check(rejected, "invalid BindingDB source or policy was admitted");
}
std::string tsv(const std::vector<std::string> &cells) {
    std::string result;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (i) result += '\t';
        result += cells[i];
    }
    return result;
}
std::vector<std::string> header() {
    return {"BindingDB Reactant_set_id", "BindingDB MonomerID", "Target Name",
        "Target Source Organism According to Curator or DataSource",
        "Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)", "Original extra"};
}
std::vector<std::string> row() {
    return {"reaction-01", "compound-01", "declared construct", "", "0",
        " <= 0.50 ", "", ">1e309", "\"literal quotes\""};
}
const std::array<std::string, 12> chain{
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
void contracts() {
    const auto h = header();
    const auto r = row();
    const auto first = tsv(r) + "\r\n";
    const auto source = tsv(h) + "\r\n" + first + tsv(r);
    std::size_t visits{};
    const auto summary = visit_bindingdb_bulk(source, limits, [&](const auto &value) {
        check(value.ordinal == visits && value.line == visits + 2, "original row order changed");
        check(value.header.columns.size() == h.size(), "original header changed");
        check(value.fields.size() == r.size(), "original field count changed");
        for (std::size_t i = 0; i < r.size(); ++i)
            check(value.fields[i] == r[i], "original field token changed");
        check(source.substr(value.byte_offset, value.byte_length) == value.original_bytes,
              "row does not name its complete original byte interval");
        check(value.original_bytes == (visits ? tsv(r) : first), "line ending or final row changed");
        check(value.affinities.size() == 4, "missing affinity column was discarded");
        check(value.affinities[0].number_token == "0" &&
              value.affinities[0].relation == BindingDBRelation::Equal, "reported zero lost");
        check(value.affinities[1].reported_token == " <= 0.50 " &&
              value.affinities[1].number_token == "0.50" &&
              value.affinities[1].relation == BindingDBRelation::LessEqual, "censored token changed");
        check(!value.affinities[2].relation && value.affinities[2].number_token.empty(),
              "missing affinity became a measured zero");
        check(value.affinities[3].number_token == "1e309" &&
              value.affinities[3].relation == BindingDBRelation::Greater,
              "lexical reading performed an overflowing floating conversion");
        ++visits;
    });
    check(summary.rows == 2 && summary.reported_affinities == 6, "repeated source rows coalesced");
    check(source == tsv(h) + "\r\n" + first + tsv(r), "source bytes mutated");
    check(visit_bindingdb_bulk(tsv(h), limits, [](const auto &) {}).rows == 0,
          "empty declared table acquired an observation");

    auto hc = h;
    for (const auto &name : chain) hc.push_back(name + " 1");
    auto rc = r;
    for (std::size_t i = 0; i < 24; ++i) rc.push_back("chain-field-" + std::to_string(i));
    const auto variable = tsv(hc) + "\n" + tsv(rc) + "\n";
    visit_bindingdb_bulk(variable, limits, [&](const auto &value) {
        check(value.header.chain_start == h.size(), "complete chain layout was not recognized");
        check(value.fields.size() == r.size() + 24, "additional chain fields were lost");
        for (std::size_t i = 0; i < 24; ++i)
            check(value.fields[r.size() + i] == "chain-field-" + std::to_string(i),
                  "chain occurrences or order changed");
    });
    auto hc2 = hc;
    for (const auto &name : chain) hc2.push_back(name + " 2");
    rc.resize(r.size() + 12);
    check(visit_bindingdb_bulk(tsv(hc2) + "\n" + tsv(rc), limits,
          [](const auto &) {}).rows == 1, "complete shorter chain row was padded or refused");
    rc.push_back("partial-chain");
    rejects([&] { visit_bindingdb_bulk(tsv(hc) + "\n" + tsv(rc), limits, [](const auto &) {}); });
    auto unknown = h;
    unknown.push_back("Original extra");
    auto extra = r; extra.push_back("second original value");
    visit_bindingdb_bulk(tsv(unknown) + "\n" + tsv(extra), limits, [&](const auto &value) {
        check(value.fields[8] != value.fields[9], "unknown repeated columns were coalesced");
    });
    for (const auto &token : {"1e", "<", "NaN", "inf", "1,000", "1.2.3", "1x", "<>2", "==1"}) {
        auto bad = r; bad[4] = token;
        rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(bad), limits, [](const auto &) {}); });
    }
    for (const auto &token : {"-0.00", "+.5", "~ 5.", "1E-999", "= 2", "> 3", "< 4", ">=5"}) {
        auto good = r; good[4] = token;
        check(visit_bindingdb_bulk(tsv(h) + "\n" + tsv(good), limits,
              [](const auto &) {}).rows == 1, "complete decimal/comparison token refused");
    }
    auto duplicate = h; duplicate.push_back(h[0]);
    rejects([&] { visit_bindingdb_bulk(tsv(duplicate), limits, [](const auto &) {}); });
    duplicate = h; duplicate.push_back(h[4]);
    rejects([&] { visit_bindingdb_bulk(tsv(duplicate), limits, [](const auto &) {}); });
    auto absent = h; absent[0] = "unrelated";
    rejects([&] { visit_bindingdb_bulk(tsv(absent), limits, [](const auto &) {}); });
    auto bad = r; bad[0].clear();
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(bad), limits, [](const auto &) {}); });
    bad = r; bad[1].clear();
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(bad), limits, [](const auto &) {}); });
    bad = r; bad.pop_back();
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(bad), limits, [](const auto &) {}); });
    bad = r; bad[2] += '\0';
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(bad), limits, [](const auto &) {}); });
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n\n", limits, [](const auto &) {}); });
    rejects([&] { visit_bindingdb_bulk({}, limits, [](const auto &) {}); });
    rejects([&] { visit_bindingdb_bulk(source, limits, {}); });
    for (const auto &bounded : {BindingDBBulkLimits{source.size() - 1, 100, 100, 1000},
         BindingDBBulkLimits{source.size(), 1, 100, 1000},
         BindingDBBulkLimits{source.size(), 100, h.size() - 1, 1000},
         BindingDBBulkLimits{source.size(), 100, 100, 4},
         BindingDBBulkLimits{0, 1, 1, 1}}) {
        rejects([&] { visit_bindingdb_bulk(source, bounded, [](const auto &) {}); });
    }
    std::stop_source cancellation;
    cancellation.request_stop();
    rejects([&] { visit_bindingdb_bulk(source, limits, [](const auto &) {}, cancellation.get_token()); });
    std::stop_source during;
    visits = 0;
    rejects([&] { visit_bindingdb_bulk(source, limits, [&](const auto &) {
        ++visits; during.request_stop();
    }, during.get_token()); });
    check(visits == 1, "cancellation continued into another source row");
    std::stop_source final;
    visits = 0;
    rejects([&] { visit_bindingdb_bulk(tsv(h) + "\n" + tsv(r), limits, [&](const auto &) {
        ++visits; final.request_stop();
    }, final.get_token()); });
    check(visits == 1, "final callback was not visited before its cancellation");
    rejects([&] { visit_bindingdb_bulk(source, limits, [](const auto &) {
        throw std::runtime_error("consumer refusal");
    }); });
}
const char *relation(std::optional<BindingDBRelation> value) {
    if (!value) return "";
    switch (*value) {
    case BindingDBRelation::Equal: return "=";
    case BindingDBRelation::Less: return "<";
    case BindingDBRelation::LessEqual: return "<=";
    case BindingDBRelation::Greater: return ">";
    case BindingDBRelation::GreaterEqual: return ">=";
    case BindingDBRelation::Approximate: return "~";
    }
    throw std::logic_error("unknown relation");
}
json_object *text(std::string_view value) {
    return json_object_new_string_len(value.empty() ? "" : value.data(), static_cast<int>(value.size()));
}
void count(json_object *object, const char *name, std::size_t value) {
    json_object_object_add(object, name, json_object_new_uint64(value));
}
void emit(json_object *object) {
    std::cout << json_object_to_json_string_ext(object, JSON_C_TO_STRING_PLAIN) << '\n';
    json_object_put(object);
}
void inspect(const char *path) {
    const auto source = gadget::read_text(path, limits.source_bytes);
    const auto summary = visit_bindingdb_bulk(source, limits, [&](const auto &row_view) {
        if (!row_view.ordinal) {
            auto *record = json_object_new_object(), *columns = json_object_new_array();
            for (const auto value : row_view.header.columns)
                json_object_array_add(columns, text(value));
            json_object_object_add(record, "header", columns);
            if (row_view.header.chain_start) count(record, "chain_start", *row_view.header.chain_start);
            emit(record);
        }
        auto *record = json_object_new_object(), *cells = json_object_new_array(),
             *affinities = json_object_new_array();
        count(record, "ordinal", row_view.ordinal); count(record, "line", row_view.line);
        count(record, "byte_offset", row_view.byte_offset);
        count(record, "byte_length", row_view.byte_length);
        for (const auto value : row_view.fields) json_object_array_add(cells, text(value));
        for (const auto &value : row_view.affinities) {
            auto *item = json_object_new_object(); count(item, "column", value.column);
            json_object_object_add(item, "endpoint", text(value.endpoint));
            json_object_object_add(item, "reported", text(value.reported_token));
            json_object_object_add(item, "number", text(value.number_token));
            json_object_object_add(item, "relation", json_object_new_string(relation(value.relation)));
            json_object_array_add(affinities, item);
        }
        json_object_object_add(record, "cells", cells);
        json_object_object_add(record, "affinities", affinities);
        emit(record);
    });
    auto *record = json_object_new_object(); count(record, "rows", summary.rows);
    count(record, "reported_affinities", summary.reported_affinities); emit(record);
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 2) inspect(argv[1]);
        else {
            if (argc != 1) throw std::invalid_argument("expected an optional source TSV path");
            contracts();
            std::cout << "PASS: " << checks << " BindingDB bulk format and refusal checks\n";
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
