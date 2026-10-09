#pragma once
#include <gadget/io.hpp>
#include <gadget/record_view.hpp>
#include <rfc/algebra/subcomplex.hpp>
#include <rfc/core/relational_complex.hpp>
#include <set>
#include <unordered_map>

namespace gadget {
class NativeSnapshot;
struct Origin {
    std::size_t file{}, record{};
};
struct SourceFile {
    std::string path, hash, parser;
    std::size_t bytes{};
};
struct Participant {
    std::uint32_t gene{};
    Rational rp{0};
    std::string literal;
    std::uint64_t occurrence{};
    std::vector<Origin> origins;
};
struct Assertion {
    std::string external, id;
    std::vector<Participant> participants;
    std::vector<std::size_t> files;
    std::set<std::string> contexts;
    std::map<std::string, RetainedRecord> observations;
};
struct Workspace {
    fs::path root;
    std::string source_digest;
    std::shared_ptr<const rfc::core::RelationalComplex> complex;
    std::shared_ptr<const NativeSnapshot> native_snapshot;
    std::optional<NativeRecord> native_publication;
    RetainedRecord snapshot_origin = Json::object();
    RetainedRecord store_origin = Json::object();
    std::vector<std::string> genes;
    std::unordered_map<std::string, std::uint32_t> gene_index;
    std::vector<Assertion> assertions;
    std::vector<SourceFile> files;
    std::vector<RawByteSource> raw_sources;
    std::map<std::string, std::unique_ptr<rfc::algebra::SubcomplexMask>> contexts;
    std::map<std::string, Vector> fields;
    std::map<std::string, RetainedRecord> field_semantics;
    std::map<std::string, std::vector<bool>> field_presence;
    // Original GET rows retain scores, identifiers and missing observations.
    std::map<std::string, RetainedRecord> measurements;
    RetainedRecord audit = Json::object();
    Json manifest() const;
    Json provenance(std::uint32_t gene) const;
};
Workspace load_pager(const fs::path &root);
void import_measurements(Workspace &source, const fs::path &configuration);
void refresh_source_digest(Workspace &source);
} // namespace gadget
