#pragma once
#include <gadget/io.hpp>
#include <gadget/native_record.hpp>

namespace gadget {
struct RecordByteBindings {
    std::vector<RawByteSource> sources;
    std::map<std::string,std::uint64_t> documents;
};
NativeRecord admit_record_view(const Json &input, RecordByteBindings bytes = {},
                              const rfc::storage::TemporalEndpoint *source = nullptr,
                              std::optional<std::uint64_t> source_namespace = {});
NativeRecordValue import_record_view(const Json &input, const RecordByteBindings &bytes = {});
Json render_record_view(const NativeRecord &record);
// Canonical layout/payload and their actual referenced byte spans.
std::string native_declaration_digest(const NativeRecord &record);
// A retained record is always a native owner or a pinned member of one. JSON
// is generated only for the existing external views and adapter calls.
class RetainedRecord {
    std::optional<NativeRecord> record_;
  public:
    RetainedRecord() = default;
    RetainedRecord(const Json &view) :record_(admit_record_view(view)) {}
    RetainedRecord(NativeRecord record) :record_(std::move(record)) {}
    RetainedRecord &operator=(const Json &view) {record_=admit_record_view(view);return *this;}
    operator Json() const {return record_?render_record_view(*record_):Json{};}
    const NativeRecord &native() const {return record_.value();}
    RetainedRecord member(std::string_view name) const {return native().at(name);}
    RetainedRecord element(std::size_t index) const {return native().at(index);}
    bool has(const std::string &name) const {return record_ && native().has(name);}
    RetainedRecord at(const std::string &name) const {return member(name);}
    RetainedRecord at(std::size_t index) const {return element(index);}
    std::size_t size() const {return native().size();}
    std::string text() const {return native().text();}
    std::string dump() const {return static_cast<Json>(*this).dump();}
    std::string dump_compact() const {return static_cast<Json>(*this).dump_compact();}
    bool boolean() const {return native().boolean();}
    std::int64_t integer() const {return std::stoll(native().integer().get_str());}
    void set(const std::string &name,const Json &value);
    void push(const Json &value);
};
// Readers route by bytes. JSON input is an explicit compatibility adapter;
// native output retains typed values and every original source byte owner.
RetainedRecord read_record(const fs::path &path);
RetainedRecord decode_record(std::span<const std::byte> bytes);
void write_record(const fs::path &path, const Json &view,
                  const rfc::storage::TemporalEndpoint *source = nullptr);
} // namespace gadget
