#pragma once
#include <filesystem>
#include <functional>
#include <span>
#include <json-c/json.h>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <rfc/numeric/exact.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace gadget {
namespace fs = std::filesystem;
using Rational = rfc::numeric::Rational;
using Vector = rfc::numeric::RationalVector;

class Json {
    std::shared_ptr<json_object> value_;
    // json-c reference counts and serialization caches are not atomic. All
    // access through this wrapper shares a short metadata lock; mathematical
    // execution uses independent canonical workspaces outside this lock.
    static std::recursive_mutex &mutex();

  public:
    explicit Json(json_object *value = nullptr)
        : value_(value, [](json_object *object) {
              std::lock_guard lock(mutex());
              json_object_put(object);
          }) {}
    static Json object() { return Json(json_object_new_object()); }
    static Json array() { return Json(json_object_new_array()); }
    static Json string(const std::string &s) {
        if (s.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("JSON string exceeds json-c integer representation");
        return Json(json_object_new_string_len(s.data(), static_cast<int>(s.size())));
    }
    static Json integer(std::int64_t n) { return Json(json_object_new_int64(n)); }
    static Json number(double n);
    static Json boolean(bool b) { return Json(json_object_new_boolean(b)); }
    // json-c grammar with complete input and representable member names.
    // Escaped NUL string values are retained; NUL member names are refused.
    static Json parse(const std::string &s);
    json_object *raw() const { return value_.get(); }
    bool has(const std::string &key) const;
    Json at(const std::string &key) const;
    Json at(std::size_t index) const;
    std::size_t size() const;
    std::string text() const;
    // Signed access refuses an unsigned value outside this return type.
    std::int64_t integer() const;
    bool boolean() const;
    void set(const std::string &key, const Json &value);
    void push(const Json &value);
    std::string dump() const;
    std::string dump_compact() const;
};
// The complete-file adapter is constrained by its string and stream
// representations. Callers may still declare a narrower byte contract.
std::size_t source_byte_limit();
// Complete bytes from one opened regular file, with size and EOF validation.
std::string read_text(const fs::path &path, std::size_t limit = source_byte_limit());
// Staged SDK replacement of a regular file. Failed writes preserve the former
// bytes. Parents must be trusted; this does not promise fsync durability.
void write_text(const fs::path &path, const std::string &text);
std::string sha256(std::string_view text);
using ByteChunkSink = std::function<void(std::span<const std::byte>)>;
std::string sha256_chunks(const std::function<void(const ByteChunkSink &)> &producer);
// Incremental identity of a complete regular file; does not retain its bytes.
std::string sha256_file(const fs::path &path);
// Application staging shared by table, assay and sequence adapters.
fs::path private_directory(const fs::path &parent, const std::string &prefix);
std::string trim(std::string text);
Rational exact_decimal(const std::string &literal);
Rational exact_rational(const std::string &literal);
struct CsvRow {
    std::size_t record;
    std::vector<std::string> cells;
};
struct CsvTable {
    std::vector<std::string> header;
    std::map<std::string, std::size_t> columns;
    std::vector<CsvRow> rows;
    const std::string &get(const CsvRow &row, const std::string &column) const;
};
CsvTable read_csv(const std::string &text, char delimiter = ',');
std::string csv_cell(const std::string &text);
Json exact_number(const Rational &value);
} // namespace gadget
