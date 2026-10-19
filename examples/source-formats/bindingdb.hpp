#pragma once
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <vector>

namespace gadget::examples {

// This reader describes original unquoted TSV bytes. It owns no source file or
// native publication. Its views are valid only while the caller retains the
// input bytes; streaming views and row/affinity spans are valid only during
// their callback. A later invalid row or cancellation can follow earlier
// callbacks: the caller must
// keep its output provisional until the entire visit succeeds.
struct BindingDBBulkLimits {
    std::size_t source_bytes{}, rows{}, columns{}, field_bytes{};
    // Opt in only when the caller has independently established completeness.
    // A decimal truncated in transit can remain a valid decimal.
    bool allow_unterminated_final_line{};
    // Independent length metadata catches reads that an underlying streambuf
    // incorrectly reports as ordinary EOF. A budget is only an upper bound.
    std::optional<std::size_t> expected_source_bytes{};
};
enum class BindingDBRelation { Equal, Less, LessEqual, Greater, GreaterEqual, Approximate };
struct BindingDBAffinityView {
    std::size_t column;
    // endpoint is one of Ki, Kd, IC50, EC50. The reader emits relation and
    // number_token together, or neither for a missing cell. Signed decimals
    // (including negative zero) remain lexical source evidence; nonnegative
    // physical interpretation belongs to admission. Bare and explicit '='
    // both mean Equal here; reported_token preserves that source distinction.
    std::string_view endpoint, reported_token, number_token;
    std::optional<BindingDBRelation> relation;
};
struct BindingDBBulkHeader {
    std::vector<std::string_view> columns;
    // Present only for a completely recognized suffix of twelve-column chain
    // blocks. Original header names and all actual row columns remain intact.
    std::optional<std::size_t> chain_start;
};
struct BindingDBBulkRowView {
    const BindingDBBulkHeader &header;
    std::size_t ordinal, line, byte_offset, byte_length;
    // byte_length == original_bytes.size(); byte_offset is the start of that
    // exact interval in the original source, including its line ending.
    std::string_view original_bytes;
    std::span<const std::string_view> fields;
    // One entry per declared affinity column, including an empty/missing cell.
    // Lexical decimal recognition performs no floating conversion, numerical
    // source admission, unit conversion, pooling or assay identity inference.
    std::span<const BindingDBAffinityView> affinities;
};
struct BindingDBBulkSummary {
    std::size_t rows{}, reported_affinities{};
    bool ends_with_newline{};
    std::size_t bytes_consumed{};
};
BindingDBBulkSummary visit_bindingdb_bulk(
    std::string_view source, const BindingDBBulkLimits &limits,
    const std::function<void(const BindingDBBulkRowView &)> &consume,
    std::stop_token stop = {},
    const std::function<void(const BindingDBBulkHeader &)> &consume_header = {});

// Same grammar and budgets, with bounded buffering rather than retaining the
// entire file. Header and row bytes are each bounded by columns * (field_bytes
// + 1), plus line endings. Header storage lasts through the visit; the row byte
// buffer is reused, while field/affinity vectors are allocated per row.
// consume_header runs once after header validation, even for a header-only file.
// Normal EOF is supported with failbit/eofbit exceptions enabled. For streambufs
// that cannot distinguish read failure from EOF, supply expected_source_bytes
// or compare bytes_consumed with independent length metadata before publication.
BindingDBBulkSummary visit_bindingdb_bulk(
    std::istream &source, const BindingDBBulkLimits &limits,
    const std::function<void(const BindingDBBulkRowView &)> &consume,
    std::stop_token stop = {},
    const std::function<void(const BindingDBBulkHeader &)> &consume_header = {});

} // namespace gadget::examples
