#include "../examples/source-formats/bindingdb.hpp"
#include <json-c/json.h>
#include <algorithm>
#include <array>
#include <charconv>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

using namespace gadget::examples;
namespace {
constexpr BindingDBBulkLimits limits{16 * 1024 * 1024, 100000, 4096, 1024 * 1024};
std::size_t checks{};
bool streaming{};
void check(bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F callback, std::string_view message) {
    try { callback(); }
    catch (const std::invalid_argument &error) {
        ++checks;
        if (error.what() != message)
            throw std::runtime_error("expected refusal '" + std::string(message) + "', got '" + error.what() + "'");
        return;
    }
    check(false, "invalid BindingDB source or policy was admitted");
}
template<class F> void cancelled(F callback) {
    try { callback(); }
    catch (const std::system_error &error) {
        check(error.code() == std::errc::operation_canceled, "cancellation error code changed");
        return;
    }
    check(false, "cancellation was ignored");
}
BindingDBBulkSummary read(std::string_view source, const BindingDBBulkLimits &budget,
    const std::function<void(const BindingDBBulkRowView &)> &consume,
    std::stop_token stop = {},
    const std::function<void(const BindingDBBulkHeader &)> &consume_header = {}) {
    if (!streaming) return visit_bindingdb_bulk(source, budget, consume, stop, consume_header);
    std::istringstream input{std::string(source)};
    return visit_bindingdb_bulk(input, budget, consume, stop, consume_header);
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
    const auto h = header(), r = row();
    auto second = r;
    second[0] = "reaction-02-distinct";
    second[2] = "another construct";
    second[4] = "= 123.4500";
    const auto first = tsv(r) + "\r\n", last = tsv(second) + "\n";
    const auto header_bytes = tsv(h) + "\r\n";
    const auto source = header_bytes + first + last;
    std::size_t visits{}, headers{};
    const std::array<std::string_view, 4> endpoints{"Ki", "IC50", "Kd", "EC50"};
    const auto summary = read(source, limits, [&](const auto &value) {
        check(headers == 1, "header callback did not precede the first row");
        const auto &expected = visits ? second : r;
        check(value.ordinal == visits && value.line == visits + 2, "original row order changed");
        check(value.header.columns.size() == h.size(), "original header changed");
        check(value.fields.size() == expected.size(), "original field count changed");
        for (std::size_t i = 0; i < expected.size(); ++i)
            check(value.fields[i] == expected[i], "original field token changed");
        check(value.byte_offset == header_bytes.size() + (visits ? first.size() : 0),
              "row byte offset changed or was reused");
        check(value.byte_length == (visits ? last.size() : first.size()), "row byte length changed");
        check(source.substr(value.byte_offset, value.byte_length) == value.original_bytes,
              "row does not name its complete original byte interval");
        check(value.original_bytes == (visits ? last : first), "line ending or row bytes changed");
        check(value.affinities.size() == 4, "missing affinity column was discarded");
        for (std::size_t i = 0; i < endpoints.size(); ++i) {
            check(value.affinities[i].column == i + 4, "affinity column order changed");
            check(value.affinities[i].endpoint == endpoints[i], "affinity endpoint mislabeled");
            check(value.header.columns[value.affinities[i].column] == std::string(endpoints[i]) + " (nM)",
                  "endpoint and column label disagree");
        }
        check(value.affinities[0].number_token == (visits ? "123.4500" : "0") &&
              value.affinities[0].relation == BindingDBRelation::Equal, "exact value or zero lost");
        check(value.affinities[1].reported_token == " <= 0.50 " &&
              value.affinities[1].number_token == "0.50" &&
              value.affinities[1].relation == BindingDBRelation::LessEqual, "censored token changed");
        check(!value.affinities[2].relation && value.affinities[2].number_token.empty(),
              "missing affinity became a measured zero");
        check(value.affinities[3].number_token == "1e309" &&
              value.affinities[3].relation == BindingDBRelation::Greater,
              "lexical reading performed an overflowing floating conversion");
        ++visits;
    }, {}, [&](const auto &value) {
        ++headers;
        check(value.columns.size() == h.size(), "header callback lost columns");
    });
    check(summary.rows == 2 && summary.reported_affinities == 6 && summary.ends_with_newline &&
          summary.bytes_consumed == source.size(),
          "source summary changed");
    check(source == header_bytes + first + last, "source bytes mutated");
    check(read(tsv(h) + "\n" + first + first, limits, [](const auto &) {}).rows == 2,
          "repeated source rows coalesced");
    headers = 0;
    check(read(tsv(h) + "\n", limits, [](const auto &) {
        check(false, "header-only table acquired an observation");
    }, {}, [&](const auto &) { ++headers; }).rows == 0 && headers == 1,
          "header-only table did not publish its header");
    rejects([&] { read(tsv(h), limits, [](const auto &) {}); },
            "BindingDB source must end with a newline");
    for (const auto ending : {"", "\r"}) {
        const auto incomplete = tsv(h) + "\n" + tsv(r) + ending;
        rejects([&] { read(incomplete, limits, [](const auto &) {}); },
                "BindingDB source must end with a newline");
        auto opt_in = limits;
        opt_in.allow_unterminated_final_line = true;
        check(!read(incomplete, opt_in, [&](const auto &value) {
            check(value.original_bytes == tsv(r) + ending, "opt-in final byte range changed");
        }).ends_with_newline, "unterminated source was reported as terminated");
    }
    auto cut = r; cut[8] = "123.";
    auto final_numeric_header = h; final_numeric_header[4] = "Original numeric field";
    final_numeric_header[8] = "Ki (nM)";
    rejects([&] { read(tsv(final_numeric_header) + "\n" + tsv(cut), limits, [](const auto &) {}); },
            "BindingDB source must end with a newline");
    for (bool numbered : {false, true}) {
        auto hc = h;
        for (const auto &name : chain) hc.push_back(name + (numbered ? " 1" : ""));
        auto rc = r;
        for (std::size_t i = 0; i < 24; ++i) rc.push_back("chain-field-" + std::to_string(i));
        read(tsv(hc) + "\n" + tsv(rc) + "\n", limits, [&](const auto &value) {
            check(value.header.chain_start == h.size(), "complete chain layout not recognized");
            check(value.fields.size() == r.size() + 24, "additional chain fields lost");
            for (std::size_t i = 0; i < 24; ++i)
                check(value.fields[r.size() + i] == "chain-field-" + std::to_string(i),
                      "chain occurrences or order changed");
        });
        auto hc2 = hc;
        for (const auto &name : chain) hc2.push_back(name + (numbered ? " 2" : ""));
        rc.resize(r.size() + 12);
        check(read(tsv(hc2) + "\n" + tsv(rc) + "\n", limits, [](const auto &) {}).rows == 1,
              "complete shorter chain row was padded or refused");
        rc.push_back("partial-chain");
        rejects([&] { read(tsv(hc) + "\n" + tsv(rc) + "\n", limits, [](const auto &) {}); },
                "BindingDB row has an incomplete or undeclared column layout");
        hc[hc.size() - 1] = "unrecognized chain label";
        rc.pop_back(); rc.insert(rc.end(), 12, "additional chain");
        rejects([&] { read(tsv(hc) + "\n" + tsv(rc) + "\n", limits, [](const auto &) {}); },
                "BindingDB row has an incomplete or undeclared column layout");
    }
    auto unknown = h; unknown.push_back("Original extra");
    auto extra = r; extra.push_back("second original value");
    read(tsv(unknown) + "\n" + tsv(extra) + "\n", limits, [&](const auto &value) {
        check(value.fields[8] != value.fields[9], "unknown repeated columns coalesced");
    });
    for (const auto token : {"1e", "<", "NaN", "inf", "1,000", "1.2.3", "1x", "<>2", "==1", ".", "+", "-", "e5"}) {
        auto bad = r; bad[4] = token;
        rejects([&] { read(tsv(h) + "\n" + tsv(bad) + "\n", limits, [](const auto &) {}); },
                "BindingDB affinity requires a complete decimal token");
    }
    struct Token { const char *reported, *number; BindingDBRelation relation; };
    for (const auto &token : std::array<Token, 11>{{
         {"-0.00", "-0.00", BindingDBRelation::Equal}, {"-5", "-5", BindingDBRelation::Equal},
         {"+.5", "+.5", BindingDBRelation::Equal}, {"~ 5.", "5.", BindingDBRelation::Approximate},
         {"1E-999", "1E-999", BindingDBRelation::Equal}, {"= 2", "2", BindingDBRelation::Equal},
         {"2", "2", BindingDBRelation::Equal}, {"> 3", "3", BindingDBRelation::Greater},
         {"< 4", "4", BindingDBRelation::Less}, {">=5", "5", BindingDBRelation::GreaterEqual},
         {" <= 0.50 ", "0.50", BindingDBRelation::LessEqual}}}) {
        auto good = r; good[4] = token.reported;
        check(read(tsv(h) + "\n" + tsv(good) + "\n", limits, [&](const auto &value) {
            const auto &actual = value.affinities[0];
            check(actual.relation == token.relation, "comparison relation changed");
            check(actual.number_token == token.number, "exact number token changed");
            check(actual.reported_token == token.reported, "reported token changed");
        }).rows == 1, "complete decimal/comparison token refused");
    }
    auto refuse = [&](const auto &cells, std::string_view message) {
        rejects([&] { read(tsv(h) + "\n" + tsv(cells) + "\n", limits, [](const auto &) {}); }, message);
    };
    auto duplicate = h; duplicate.push_back(h[0]);
    rejects([&] { read(tsv(duplicate) + "\n", limits, [](const auto &) {}); }, "BindingDB semantic column is duplicated");
    duplicate = h; duplicate.push_back(h[4]);
    rejects([&] { read(tsv(duplicate) + "\n", limits, [](const auto &) {}); }, "BindingDB semantic column is duplicated");
    auto absent = h; absent[0] = "unrelated";
    rejects([&] { read(tsv(absent) + "\n", limits, [](const auto &) {}); }, "BindingDB required column is absent");
    rejects([&] { read("\xef\xbb\xbf" + source, limits, [](const auto &) {}); }, "BindingDB header contains a UTF-8 BOM");
    auto bad = r; bad[0].clear(); refuse(bad, "BindingDB reaction or monomer identity is empty");
    bad = r; bad[1].clear(); refuse(bad, "BindingDB reaction or monomer identity is empty");
    bad = r; bad.pop_back(); refuse(bad, "BindingDB row has an incomplete or undeclared column layout");
    bad = r; bad[2] += '\0'; refuse(bad, "BindingDB field contains an unsupported control byte");
    bad = r; bad[2] += "\rinside"; refuse(bad, "BindingDB field contains an unsupported control byte");
    rejects([&] { read(tsv(h) + "\n\n", limits, [](const auto &) {}); }, "BindingDB blank row is not an observation");
    rejects([&] { read({}, limits, [](const auto &) {}); }, "BindingDB source is empty");
    rejects([&] { read(source, limits, {}); }, "BindingDB consumer is empty");
    BindingDBBulkLimits uninitialized;
    check(!uninitialized.source_bytes && !uninitialized.rows && !uninitialized.columns && !uninitialized.field_bytes,
          "default policy contains uninitialized budgets");
    rejects([&] { read(source, uninitialized, [](const auto &) {}); }, "BindingDB budgets must be positive");
    std::size_t longest{};
    for (const auto &record : {h, r, second})
        for (const auto &cell : record) longest = std::max(longest, cell.size());
    const BindingDBBulkLimits exact{source.size(), 2, h.size(), longest};
    check(read(source, exact, [](const auto &) {}).rows == 2, "exact policy boundary was refused");
    auto bounded = exact; --bounded.source_bytes;
    rejects([&] { read(source, bounded, [](const auto &) {}); }, "BindingDB source byte budget exceeded");
    bounded = exact; --bounded.rows;
    rejects([&] { read(source, bounded, [](const auto &) {}); }, "BindingDB row budget exceeded");
    bounded = exact; --bounded.columns;
    rejects([&] { read(source, bounded, [](const auto &) {}); }, "BindingDB column budget exceeded");
    bounded = exact; --bounded.field_bytes;
    rejects([&] { read(source, bounded, [](const auto &) {}); }, "BindingDB field byte budget exceeded");
    bounded = exact; bounded.expected_source_bytes = source.size();
    check(read(source, bounded, [](const auto &) {}).bytes_consumed == source.size(),
          "independently known source length was not retained");
    bounded.expected_source_bytes = source.size() - last.size();
    rejects([&] { read(source, bounded, [](const auto &) {}); },
            "BindingDB source byte count differs from expected size");
    bounded = limits; bounded.expected_source_bytes = source.size();
    rejects([&] { read(header_bytes + first, bounded, [](const auto &) {}); },
            "BindingDB source byte count differs from expected size");
    std::stop_source before; before.request_stop();
    cancelled([&] { read(source, limits, [](const auto &) {}, before.get_token()); });
    std::stop_source during; visits = 0;
    cancelled([&] { read(source, limits, [&](const auto &) {
        ++visits; during.request_stop();
    }, during.get_token()); });
    check(visits == 1, "cancellation continued into another source row");
    std::stop_source final; visits = 0;
    cancelled([&] { read(tsv(h) + "\n" + first, limits, [&](const auto &) {
        ++visits; final.request_stop();
    }, final.get_token()); });
    check(visits == 1, "final callback was not visited before its cancellation");
    std::stop_source after_header;
    cancelled([&] { read(tsv(h) + "\n", limits, [](const auto &) {}, after_header.get_token(),
        [&](const auto &) { after_header.request_stop(); }); });
    struct ConsumerError : std::runtime_error { using std::runtime_error::runtime_error; };
    try {
        read(source, limits, [](const auto &) { throw ConsumerError("consumer refusal"); });
        check(false, "consumer exception was swallowed");
    } catch (const ConsumerError &error) {
        check(std::string_view(error.what()) == "consumer refusal", "consumer exception changed");
    }
    auto long_row = r; long_row[8] = std::string(130000, 'x');
    check(read(tsv(h) + "\r\n" + tsv(long_row) + "\r\n", limits, [&](const auto &value) {
        check(value.fields[8] == long_row[8], "field split across read buffers changed");
    }).rows == 1, "buffer-spanning source lost a row");
}
void stream_contracts() {
    const auto h = tsv(header()) + "\n";
    // Pin CRLF immediately before, across, and after the 64 KiB refill edge.
    for (std::size_t cr_offset : {65534, 65535, 65536}) {
        auto values = row();
        values[8].append(cr_offset - h.size() - tsv(values).size(), 'x');
        const auto source = h + tsv(values) + "\r\n";
        std::istringstream input(source);
        const auto summary = visit_bindingdb_bulk(input, limits, [&](const auto &value) {
            check(value.original_bytes == source.substr(h.size()), "refill split changed CRLF bytes");
            check(value.fields[8] == values[8], "refill split changed last field");
        });
        check(summary.bytes_consumed == source.size(), "refill split changed consumed length");
    }
    const auto source = h + tsv(row()) + "\n";
    for (auto mask : {std::ios::failbit, std::ios::eofbit,
                      std::ios::failbit | std::ios::eofbit | std::ios::badbit}) {
        std::istringstream input(source);
        input.exceptions(mask);
        check(visit_bindingdb_bulk(input, limits, [](const auto &) {}).bytes_consumed == source.size(),
              "normal EOF with stream exceptions enabled was refused");
    }
    class CutBuffer : public std::stringbuf {
        std::size_t remaining;
        std::streamsize xsgetn(char *output, std::streamsize count) override {
            const auto available = std::min<std::size_t>(remaining, static_cast<std::size_t>(count));
            const auto taken = std::stringbuf::xsgetn(output, static_cast<std::streamsize>(available));
            remaining -= static_cast<std::size_t>(taken);
            return taken;
        }
    public:
        CutBuffer(const std::string &source, std::size_t cutoff) : std::stringbuf(source), remaining(cutoff) {}
    } cut(source + tsv(row()) + "\n", source.size());
    std::istream faulty(&cut);
    auto expected = limits;
    expected.expected_source_bytes = source.size() + tsv(row()).size() + 1;
    rejects([&] { visit_bindingdb_bulk(faulty, expected, [](const auto &) {}); },
            "BindingDB source byte count differs from expected size");
    check(faulty.eof() && !faulty.bad(), "fault injection did not reproduce EOF-only failure");

    std::stop_source stop;
    class CancelBuffer : public std::stringbuf {
        std::stop_source &stop;
        std::size_t calls{};
        std::streamsize xsgetn(char *output, std::streamsize count) override {
            const auto taken = std::stringbuf::xsgetn(output, count);
            if (++calls == 2) stop.request_stop();
            return taken;
        }
    public:
        CancelBuffer(const std::string &source, std::stop_source &stop) : std::stringbuf(source), stop(stop) {}
    };
    auto large = row(); large[8] = std::string(130000, 'x');
    CancelBuffer buffer(h + tsv(large) + "\n", stop);
    std::istream input(&buffer);
    std::size_t visits{};
    cancelled([&] { visit_bindingdb_bulk(input, limits, [&](const auto &) { ++visits; }, stop.get_token()); });
    check(!visits, "cancellation during refill still emitted its row");
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
json_object *allocated(json_object *value) {
    if (!value) throw std::bad_alloc();
    return value;
}
json_object *object() { return allocated(json_object_new_object()); }
json_object *array() { return allocated(json_object_new_array()); }
json_object *text(std::string_view value) {
    if (value.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("JSON string is too large");
    return allocated(json_object_new_string_len(value.empty() ? "" : value.data(), static_cast<int>(value.size())));
}
void count(json_object *object, const char *name, std::size_t value) {
    json_object_object_add(object, name, allocated(json_object_new_uint64(value)));
}
void emit(json_object *object) {
    const std::unique_ptr<json_object, decltype(&json_object_put)> owner(object, json_object_put);
    std::cout << json_object_to_json_string_ext(object, JSON_C_TO_STRING_PLAIN) << '\n';
}
std::size_t positive(std::string_view value) {
    std::size_t number{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), number);
    if (error != std::errc{} || end != value.data() + value.size() || !number)
        throw std::invalid_argument("limits must be positive size_t integers");
    return number;
}
void inspect(int argc, char **argv) {
    const auto size = std::filesystem::file_size(argv[1]);
    if (size > std::numeric_limits<std::size_t>::max()) throw std::length_error("source size exceeds size_t");
    const auto bytes = std::max<std::size_t>(1, static_cast<std::size_t>(size));
    BindingDBBulkLimits budget{bytes, bytes, 4096, 1024 * 1024};
    budget.expected_source_bytes = static_cast<std::size_t>(size);
    std::set<std::string> seen;
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 == argc) throw std::invalid_argument("limit option requires a value");
        const std::string option = argv[i];
        if (!seen.insert(option).second) throw std::invalid_argument("duplicate limit option");
        const auto number = positive(argv[i + 1]);
        if (option == "--source-bytes") budget.source_bytes = number;
        else if (option == "--rows") budget.rows = number;
        else if (option == "--columns") budget.columns = number;
        else if (option == "--field-bytes") budget.field_bytes = number;
        else throw std::invalid_argument("unknown limit option");
    }
    std::ifstream input(argv[1], std::ios::binary);
    if (!input) throw std::runtime_error("source file could not be opened");
    const auto summary = visit_bindingdb_bulk(input, budget, [&](const auto &row_view) {
        auto *record = object(), *cells = array(), *affinities = array();
        count(record, "ordinal", row_view.ordinal); count(record, "line", row_view.line);
        count(record, "byte_offset", row_view.byte_offset);
        count(record, "byte_length", row_view.byte_length);
        for (const auto value : row_view.fields) json_object_array_add(cells, text(value));
        for (const auto &value : row_view.affinities) {
            auto *item = object(); count(item, "column", value.column);
            json_object_object_add(item, "endpoint", text(value.endpoint));
            json_object_object_add(item, "reported", text(value.reported_token));
            json_object_object_add(item, "number", text(value.number_token));
            json_object_object_add(item, "relation", text(relation(value.relation)));
            json_object_array_add(affinities, item);
        }
        json_object_object_add(record, "cells", cells);
        json_object_object_add(record, "affinities", affinities);
        emit(record);
    }, {}, [&](const auto &header) {
        auto *record = object(), *columns = array(), *declared = object();
        for (const auto value : header.columns) json_object_array_add(columns, text(value));
        json_object_object_add(record, "header", columns);
        if (header.chain_start) count(record, "chain_start", *header.chain_start);
        count(declared, "source_bytes", budget.source_bytes); count(declared, "rows", budget.rows);
        count(declared, "columns", budget.columns); count(declared, "field_bytes", budget.field_bytes);
        count(declared, "expected_source_bytes", *budget.expected_source_bytes);
        json_object_object_add(declared, "allow_unterminated_final_line", allocated(json_object_new_boolean(false)));
        json_object_object_add(record, "limits", declared);
        emit(record);
    });
    if (summary.bytes_consumed != size)
        throw std::runtime_error("BindingDB source byte count differs from file size");
    auto *record = object(); count(record, "rows", summary.rows);
    count(record, "reported_affinities", summary.reported_affinities);
    count(record, "bytes_consumed", summary.bytes_consumed);
    json_object_object_add(record, "ends_with_newline", allocated(json_object_new_boolean(summary.ends_with_newline)));
    emit(record);
}
} // namespace
int main(int argc, char **argv) {
    try {
        // Exercise ostream's error reporting rather than dying on a closed pipe.
        if (std::signal(SIGPIPE, SIG_IGN) == SIG_ERR)
            throw std::runtime_error("could not configure output error handling");
        std::cout.exceptions(std::ios::badbit | std::ios::failbit);
        if (argc >= 2) inspect(argc, argv);
        else {
            contracts();
            streaming = true;
            contracts();
            stream_contracts();
            class BrokenBuffer : public std::streambuf {
                std::streamsize xsgetn(char *, std::streamsize) override {
                    throw std::runtime_error("injected input failure");
                }
            } buffer;
            std::istream input(&buffer);
            try {
                visit_bindingdb_bulk(input, limits, [](const auto &) {});
                check(false, "input read failure was swallowed");
            } catch (const std::ios_base::failure &) { check(true, "input error propagated"); }
            std::cout << "PASS: " << checks << " BindingDB bulk format and refusal checks\n";
        }
        std::cout.flush();
        return 0;
    } catch (const std::exception &error) {
        std::cout.exceptions(std::ios::goodbit);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
