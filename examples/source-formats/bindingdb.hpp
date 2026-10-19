#pragma once
#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <stop_token>
#include <string_view>
#include <vector>

namespace gadget::examples {

// This reader describes original unquoted TSV bytes. It owns no source file or
// native publication. Its views are valid only while the caller retains the
// input bytes; row/affinity spans are valid only during the callback. A later
// invalid row or cancellation can follow earlier callbacks: the caller must
// keep its output provisional until the entire visit succeeds.
struct BindingDBBulkLimits {
    std::size_t source_bytes, rows, columns, field_bytes;
};
enum class BindingDBRelation { Equal, Less, LessEqual, Greater, GreaterEqual, Approximate };
struct BindingDBAffinityView {
    std::size_t column;
    std::string_view endpoint, reported_token, number_token;
    std::optional<BindingDBRelation> relation;
};
struct BindingDBBulkHeader {
    std::vector<std::string_view> columns;
    // Present only for a completely recognized suffix of twelve column chain
    // blocks. Original header names and all actual row columns remain intact.
    std::optional<std::size_t> chain_start;
};
struct BindingDBBulkRowView {
    const BindingDBBulkHeader &header;
    std::size_t ordinal, line, byte_offset, byte_length;
    std::string_view original_bytes;
    std::span<const std::string_view> fields;
    // One entry per declared affinity column, including an empty/missing cell.
    // Lexical decimal recognition performs no floating conversion, numerical
    // source admission, unit conversion, pooling or assay identity inference.
    std::span<const BindingDBAffinityView> affinities;
};
struct BindingDBBulkSummary {
    std::size_t rows{}, reported_affinities{};
};
BindingDBBulkSummary visit_bindingdb_bulk(
    std::string_view source, const BindingDBBulkLimits &limits,
    const std::function<void(const BindingDBBulkRowView &)> &consume,
    std::stop_token stop = {});

} // namespace gadget::examples
