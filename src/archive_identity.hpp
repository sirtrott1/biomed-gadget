#pragma once
#include <gadget/record_view.hpp>
#include <stop_token>

namespace gadget {
inline std::string archive_identity(const NativeRecord &record, std::size_t max_bytes,
                                    std::stop_token stop = {}) {
    const auto history = record.publication_history();
    if (!record.is_root() || !history) throw std::invalid_argument("archive identity requires its owning publication");
    std::size_t consumed{};
    return "sha256:" + sha256_chunks([&](const auto &sink) {
        rfc::storage::stream_native_series(*history, [&](std::span<const std::byte> bytes) {
            if (stop.stop_requested()) throw std::runtime_error("native archive hashing cancelled");
            if (bytes.size() > max_bytes - consumed) throw std::length_error("native archive exceeds its byte admission");
            consumed += bytes.size(); sink(bytes);
        });
    });
}
}
