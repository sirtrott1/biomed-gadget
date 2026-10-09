#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <gadget/io.hpp>
#include <iomanip>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <sstream>
#include <stdexcept>

namespace gadget {
std::recursive_mutex &Json::mutex() {
    static std::recursive_mutex lock;
    return lock;
}
Json Json::number(double n) {
    if (!std::isfinite(n))
        throw std::invalid_argument("nonfinite JSON result");
    return Json(json_object_new_double(n));
}
Json Json::parse(const std::string &s) {
    if (s.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("JSON input exceeds json-c integer representation");
    // Grow only when the codec reports insufficient stack depth. Every nested
    // container in valid JSON needs two source bytes; no second parser or
    // application nesting ceiling is needed to size the complete input.
    std::size_t depth = JSON_TOKENER_DEFAULT_DEPTH;
    const auto source_depth = s.size() / 2 + 1;
    for (;;) {
        auto tok = std::unique_ptr<json_tokener, decltype(&json_tokener_free)>(
            json_tokener_new_ex(static_cast<int>(depth)), json_tokener_free);
        if (!tok)
            throw std::bad_alloc{};
        json_tokener_set_flags(tok.get(), JSON_TOKENER_STRICT | JSON_TOKENER_VALIDATE_UTF8);
        Json result(json_tokener_parse_ex(tok.get(), s.c_str(), static_cast<int>(s.size())));
        const auto error = json_tokener_get_error(tok.get());
        if (error == json_tokener_error_depth && depth < source_depth) {
            depth = std::min(depth * 2, source_depth);
            continue;
        }
        if (error != json_tokener_success)
            throw std::invalid_argument("invalid JSON");
        const auto consumed = json_tokener_get_parse_end(tok.get());
        if (!trim(s.substr(consumed)).empty())
            throw std::invalid_argument("trailing JSON content");
        return result;
    }
}
bool Json::has(const std::string &key) const {
    std::lock_guard lock(mutex());
    json_object *out = nullptr;
    return json_object_is_type(raw(), json_type_object) &&
           json_object_object_get_ex(raw(), key.c_str(), &out);
}
Json Json::at(const std::string &key) const {
    std::lock_guard lock(mutex());
    json_object *out = nullptr;
    if (!json_object_is_type(raw(), json_type_object) ||
        !json_object_object_get_ex(raw(), key.c_str(), &out))
        throw std::invalid_argument("missing JSON member: " + key);
    return Json(json_object_get(out));
}
Json Json::at(std::size_t index) const {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_array) || index >= size())
        throw std::invalid_argument("JSON array index");
    return Json(json_object_get(json_object_array_get_idx(raw(), index)));
}
std::size_t Json::size() const {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_array))
        throw std::invalid_argument("expected JSON array");
    return json_object_array_length(raw());
}
std::string Json::text() const {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_string))
        throw std::invalid_argument("expected JSON string");
    return {json_object_get_string(raw()),
            static_cast<std::size_t>(json_object_get_string_len(raw()))};
}
std::int64_t Json::integer() const {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_int))
        throw std::invalid_argument("expected JSON integer");
    return json_object_get_int64(raw());
}
bool Json::boolean() const {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_boolean))
        throw std::invalid_argument("expected JSON boolean");
    return json_object_get_boolean(raw());
}
void Json::set(const std::string &key, const Json &value) {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_object))
        throw std::logic_error("JSON object required");
    if (json_object_object_add(raw(), key.c_str(), json_object_get(value.raw())) != 0)
        throw std::runtime_error("JSON allocation");
}
void Json::push(const Json &value) {
    std::lock_guard lock(mutex());
    if (!json_object_is_type(raw(), json_type_array) ||
        json_object_array_add(raw(), json_object_get(value.raw())) != 0)
        throw std::runtime_error("JSON array allocation");
}
std::string Json::dump() const {
    std::lock_guard lock(mutex());
    const auto *bytes = json_object_to_json_string_ext(raw(), JSON_C_TO_STRING_PRETTY |
                                                                  JSON_C_TO_STRING_NOSLASHESCAPE);
    if (!bytes)
        throw std::length_error("JSON serialization exceeds available representation or memory");
    return bytes;
}
std::string Json::dump_compact() const {
    std::lock_guard lock(mutex());
    const auto *bytes = json_object_to_json_string_ext(raw(), JSON_C_TO_STRING_PLAIN |
                                                                  JSON_C_TO_STRING_NOSLASHESCAPE);
    if (!bytes)
        throw std::length_error("JSON serialization exceeds available representation or memory");
    return bytes;
}
std::size_t source_byte_limit() {
    return std::min(std::string{}.max_size(),
                    static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()));
}
std::string read_text(const fs::path &p, std::size_t limit) {
    std::ifstream f(p, std::ios::binary);
    if (!f)
        throw std::runtime_error("cannot open " + p.string());
    auto size = fs::file_size(p);
    if (size > std::min(limit, source_byte_limit()))
        throw std::runtime_error("input size budget exceeded: " + p.string());
    std::string s(size, '\0');
    f.read(s.data(), static_cast<std::streamsize>(size));
    if (!f && size)
        throw std::runtime_error("cannot read " + p.string());
    return s;
}
void write_text(const fs::path &p, const std::string &text) {
    fs::create_directories(p.parent_path());
    std::ofstream f(p, std::ios::binary | std::ios::trunc);
    if (!f || !f.write(text.data(), static_cast<std::streamsize>(text.size())))
        throw std::runtime_error("cannot write " + p.string());
}
std::string sha256(std::string_view text) {
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int n = 0;
    if (EVP_Digest(text.data(), text.size(), digest, &n, EVP_sha256(), nullptr) != 1)
        throw std::runtime_error("SHA-256 failed");
    std::ostringstream s;
    s << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < n; ++i)
        s << std::setw(2) << static_cast<unsigned int>(digest[i]);
    return s.str();
}
std::string sha256_file(const fs::path &path) {
    if (!fs::is_regular_file(path))
        throw std::invalid_argument("file identity requires a regular file");
    std::ifstream stream(path, std::ios::binary);
    auto digest =
        std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)>(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!stream || !digest || EVP_DigestInit_ex(digest.get(), EVP_sha256(), nullptr) != 1)
        throw std::runtime_error("cannot initialize file identity");
    std::array<char, 1024 * 1024> buffer{};
    while (stream) {
        stream.read(buffer.data(), buffer.size());
        if (EVP_DigestUpdate(digest.get(), buffer.data(), stream.gcount()) != 1)
            throw std::runtime_error("file identity hashing failed");
    }
    if (!stream.eof())
        throw std::runtime_error("file identity read failed");
    unsigned char hash[EVP_MAX_MD_SIZE];
    unsigned size{};
    if (EVP_DigestFinal_ex(digest.get(), hash, &size) != 1 || size != 32)
        throw std::runtime_error("file identity finalization failed");
    constexpr char digits[] = "0123456789abcdef";
    std::string result;
    for (unsigned i = 0; i < size; ++i) {
        result += digits[hash[i] >> 4];
        result += digits[hash[i] & 15];
    }
    return result;
}
std::string sha256_chunks(const std::function<void(const ByteChunkSink &)> &producer) {
    auto digest=std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)>(EVP_MD_CTX_new(),EVP_MD_CTX_free);
    if(!digest || EVP_DigestInit_ex(digest.get(),EVP_sha256(),nullptr)!=1)
        throw std::runtime_error("cannot initialize stream identity");
    producer([&](auto bytes) {
        if(EVP_DigestUpdate(digest.get(),bytes.data(),bytes.size())!=1)
            throw std::runtime_error("stream identity hashing failed");
    });
    unsigned char hash[EVP_MAX_MD_SIZE];unsigned size{};
    if(EVP_DigestFinal_ex(digest.get(),hash,&size)!=1 || size!=32)
        throw std::runtime_error("stream identity finalization failed");
    std::ostringstream result;result<<std::hex<<std::setfill('0');
    for(unsigned i=0;i<size;++i)result<<std::setw(2)<<static_cast<unsigned>(hash[i]);
    return result.str();
}
fs::path private_directory(const fs::path &parent, const std::string &prefix) {
    for (;;) {
        unsigned char random[16];
        if (RAND_bytes(random, sizeof(random)) != 1)
            throw std::runtime_error("private source staging identity");
        std::string name = prefix;
        for (auto byte : random) {
            name += "0123456789abcdef"[byte >> 4];
            name += "0123456789abcdef"[byte & 15];
        }
        const auto path = parent / name;
        if (!fs::create_directory(path))
            continue;
        try {
            fs::permissions(path, fs::perms::owner_all, fs::perm_options::replace);
        } catch (...) {
            std::error_code ec;
            fs::remove(path, ec);
            throw;
        }
        return path;
    }
}
std::string trim(std::string s) {
    auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return {};
    auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}
Rational exact_rational(const std::string &literal) {
    const auto token = trim(literal);
    if (token.empty() || token.size() > source_byte_limit() ||
        token.find_first_not_of("0123456789/-") != std::string::npos)
        throw std::invalid_argument("invalid or oversized rational token");
    Rational value(token);
    if (value.get_den() == 0)
        throw std::invalid_argument("zero rational denominator");
    value.canonicalize();
    return value;
}
Rational exact_decimal(const std::string &input) {
    auto s = trim(input);
    if (s.empty() || s.size() > 4096)
        throw std::invalid_argument("missing or oversized exact decimal");
    std::size_t i = 0;
    bool negative = false;
    if (s[i] == '+' || s[i] == '-') {
        negative = s[i] == '-';
        ++i;
    }
    std::string digits;
    long fractional = 0;
    bool point = false;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (c >= '0' && c <= '9') {
            digits += c;
            if (point)
                ++fractional;
        } else if (c == '.' && !point) {
            point = true;
        } else
            break;
    }
    if (digits.empty())
        throw std::invalid_argument("invalid decimal: " + s);
    long exponent = 0;
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        ++i;
        bool eneg = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) {
            eneg = s[i] == '-';
            ++i;
        }
        std::size_t begin = i;
        for (; i < s.size() && s[i] >= '0' && s[i] <= '9'; ++i) {
            exponent = exponent * 10 + s[i] - '0';
            if (exponent > 4096)
                throw std::invalid_argument("decimal exponent budget exceeded");
        }
        if (i == begin)
            throw std::invalid_argument("invalid decimal exponent");
        if (eneg)
            exponent = -exponent;
    }
    if (i != s.size())
        throw std::invalid_argument("invalid decimal: " + s);
    rfc::numeric::Integer numerator(digits, 10), denominator = 1, power = 1;
    const long scale = fractional - exponent;
    if (std::abs(scale) > 4096)
        throw std::invalid_argument("decimal scale budget exceeded");
    mpz_ui_pow_ui(power.get_mpz_t(), 10, static_cast<unsigned long>(std::abs(scale)));
    if (scale >= 0)
        denominator = power;
    else
        numerator *= power;
    if (negative)
        numerator = -numerator;
    Rational value(numerator, denominator);
    value.canonicalize();
    return value;
}
CsvTable read_csv(const std::string &text, char delimiter) {
    if (delimiter != ',' && delimiter != '\t')
        throw std::invalid_argument("unsupported table delimiter");
    CsvTable table;
    std::vector<std::string> cells;
    std::string cell;
    bool quoted = false, closed = false;
    std::size_t record = 0;
    auto emit_cell = [&] {
        cells.push_back(cell);
        cell.clear();
        closed = false;
    };
    auto emit_record = [&] {
        emit_cell();
        ++record;
        if (record == 1) {
            table.header = cells;
            for (std::size_t i = 0; i < cells.size(); ++i) {
                auto name = trim(cells[i]);
                if (name.empty() || !table.columns.emplace(name, i).second)
                    throw std::invalid_argument("empty/duplicate CSV header");
            }
        } else if (!(cells.size() == 1 && trim(cells[0]).empty())) {
            if (cells.size() != table.header.size())
                throw std::invalid_argument("CSV width at record " + std::to_string(record));
            table.rows.push_back({record, cells});
        }
        cells.clear();
    };
    for (std::size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        if (quoted) {
            if (c == '"') {
                if (i + 1 < text.size() && text[i + 1] == '"') {
                    cell += '"';
                    ++i;
                } else {
                    quoted = false;
                    closed = true;
                }
            } else
                cell += c;
        } else if (c == delimiter)
            emit_cell();
        else if (c == '\n' || c == '\r') {
            if (c == '\r' && i + 1 < text.size() && text[i + 1] == '\n')
                ++i;
            emit_record();
        } else if (c == '"' && cell.empty() && !closed)
            quoted = true;
        else {
            if (closed)
                throw std::invalid_argument("content after closing CSV quote");
            if (c == '"')
                throw std::invalid_argument("quote inside unquoted CSV cell");
            cell += c;
        }
    }
    if (quoted)
        throw std::invalid_argument("unterminated CSV quote");
    if (!cell.empty() || !cells.empty() || closed)
        emit_record();
    return table;
}
const std::string &CsvTable::get(const CsvRow &row, const std::string &column) const {
    auto it = columns.find(column);
    if (it == columns.end())
        throw std::invalid_argument("missing CSV column: " + column);
    return row.cells.at(it->second);
}
std::string csv_cell(const std::string &s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos)
        return s;
    std::string out = "\"";
    for (char c : s) {
        out += c;
        if (c == '\"')
            out += c;
    }
    return out + '\"';
}
Json exact_number(const Rational &value) {
    auto j = Json::object();
    j.set("rational", Json::string(value.get_str()));
    auto approximate = value.get_d();
    if (std::isfinite(approximate))
        j.set("approximation", Json::number(approximate));
    return j;
}
} // namespace gadget
