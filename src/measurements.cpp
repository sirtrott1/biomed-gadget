#include <algorithm>
#include <gadget/source.hpp>
#include <gadget/snapshot.hpp>
#include <set>
#include <stdexcept>

namespace gadget {
namespace {
void keys(const Json &object, const std::set<std::string> &allowed) {
    if (!json_object_is_type(object.raw(), json_type_object))
        throw std::invalid_argument("measurement declaration must be an object");
    json_object_object_foreach(object.raw(), key, value) {
        (void)value;
        if (!allowed.contains(key))
            throw std::invalid_argument(std::string("unknown measurement declaration member: ") +
                                        key);
    }
}
struct Captured {
    SourceFile file;
    std::string bytes;
};
Captured capture(const Workspace &source, const fs::path &path, const std::string &parser) {
    auto absolute = fs::canonical(path);
    auto relative = absolute.lexically_relative(source.root);
    if (relative.empty() || *relative.begin() == "..")
        throw std::invalid_argument("measurement input escapes declared source root");
    auto bytes = read_text(absolute);
    return {{relative.generic_string(), "sha256:" + sha256(bytes), parser, bytes.size()},
            std::move(bytes)};
}
struct StagedField {
    std::string id;
    Vector values;
    std::vector<bool> present;
    Json semantics;
};
} // namespace
void refresh_source_digest(Workspace &source) {
    auto files = source.files;
    std::sort(files.begin(), files.end(),
              [](const auto &a, const auto &b) { return a.path < b.path; });
    std::string digest = "gadget.pager.source.v1\n";
    for (const auto &f : files)
        digest += std::to_string(f.path.size()) + ":" + f.path + f.hash + "\n";
    source.source_digest = "sha256:" + sha256(digest);
}
void import_measurements(Workspace &source, const fs::path &configuration) {
    if (source.snapshot_origin.has("schema") || source.store_origin.has("root"))
        throw std::invalid_argument("snapshot is immutable; import measurements before capture");
    auto config = capture(source, configuration, "gadget.measurement_declaration.v1");
    auto declaration = Json::parse(config.bytes);
    keys(declaration, {"schema", "context", "path", "delimiter", "identity_column", "fields"});
    if (declaration.at("schema").text() != "gadget.measurements.v1")
        throw std::invalid_argument("unsupported measurement declaration schema");
    auto context = declaration.at("context").text();
    if (!source.contexts.contains(context))
        throw std::invalid_argument("unknown measurement context");
    auto relative = fs::path(declaration.at("path").text());
    if (relative.is_absolute())
        throw std::invalid_argument("measurement table path must be relative to source root");
    auto input = capture(source, source.root / relative, "gadget.measurements_table.v1");
    auto delimiter = declaration.has("delimiter") ? declaration.at("delimiter").text() : ",";
    if (delimiter.size() != 1)
        throw std::invalid_argument("measurement delimiter must be comma or tab");
    auto table = read_csv(input.bytes, delimiter[0]);
    auto identity = declaration.at("identity_column").text();
    if (!table.columns.contains(identity))
        throw std::invalid_argument("missing measurement identity column");
    auto specs = declaration.at("fields");
    if (specs.size() == 0 || specs.size() > 128)
        throw std::invalid_argument("measurement field budget (1..128)");
    const std::set<std::string> kinds = {"measurement",    "count",       "effect_size",
                                         "standard_error", "p_value",     "adjusted_p_value",
                                         "score",          "probability", "binary_annotation"};
    std::set<std::string> identifiers;
    std::vector<StagedField> staged;
    for (std::size_t i = 0; i < specs.size(); ++i) {
        auto spec = specs.at(i);
        keys(spec, {"id", "column", "kind", "units", "method", "null_model", "missing_policy",
                    "description"});
        auto id = spec.at("id").text();
        if (!id.starts_with(context + ".") || source.fields.contains(id) ||
            !identifiers.insert(id).second)
            throw std::invalid_argument("measurement field identity conflict: " + id);
        if (!kinds.contains(spec.at("kind").text()))
            throw std::invalid_argument("unknown measurement field kind");
        if (spec.at("method").text().empty())
            throw std::invalid_argument("measurement method must be declared");
        (void)spec.at("units").text();
        (void)spec.at("null_model").text();
        auto policy = spec.at("missing_policy").text();
        if (policy != "reject" && policy != "zero_extension")
            throw std::invalid_argument("missing policy must be reject or zero_extension");
        if (!table.columns.contains(spec.at("column").text()))
            throw std::invalid_argument("measurement column absent");
        auto semantics = Json::parse(spec.dump());
        semantics.set("source_file", Json::integer(source.files.size() + 1));
        semantics.set("declaration_file", Json::integer(source.files.size()));
        semantics.set("grade", Json::integer(0));
        semantics.set("context", Json::string(context));
        semantics.set(
            "representation",
            Json::string(
                "exact reported decimal; observation uncertainty is retained by method metadata"));
        staged.push_back({id, Vector(source.genes.size(), Rational{0}),
                          std::vector<bool>(source.genes.size(), false), semantics});
    }
    auto observations = Json::array();
    std::vector<std::set<std::uint32_t>> seen(staged.size());
    for (const auto &row : table.rows) {
        auto gene = trim(table.get(row, identity));
        auto it = source.gene_index.find(gene);
        if (it == source.gene_index.end())
            throw std::invalid_argument("unmapped measurement gene: " + gene);
        auto index = it->second;
        auto record = Json::object();
        record.set("gene", Json::string(gene));
        record.set("source_record", Json::integer(row.record));
        record.set("source_file", Json::integer(source.files.size() + 1));
        auto original = Json::object();
        for (std::size_t i = 0; i < table.header.size(); ++i)
            original.set(table.header[i], Json::string(row.cells[i]));
        record.set("reported", original);
        observations.push(record);
        for (std::size_t i = 0; i < staged.size(); ++i) {
            auto spec = specs.at(i);
            auto literal = trim(table.get(row, spec.at("column").text()));
            bool present = !literal.empty() && literal != "NA" && literal != "N/A";
            if (!present && spec.at("missing_policy").text() == "reject")
                throw std::invalid_argument("missing measurement observation at record " +
                                            std::to_string(row.record));
            Rational value = present ? exact_decimal(literal) : Rational{0};
            auto kind = spec.at("kind").text();
            if (present &&
                (kind == "p_value" || kind == "adjusted_p_value" || kind == "probability") &&
                (value < 0 || value > 1))
                throw std::invalid_argument("probability field outside [0,1]");
            if (present && kind == "count" && (value < 0 || value.get_den() != 1))
                throw std::invalid_argument("count must be a nonnegative integer");
            if (present && kind == "binary_annotation" && value != 0 && value != 1)
                throw std::invalid_argument("binary annotation must be zero or one");
            if (!seen[i].insert(index).second &&
                (staged[i].values[index] != value || staged[i].present[index] != present))
                throw std::invalid_argument("conflicting repeated measurement identity");
            staged[i].values[index] = value;
            staged[i].present[index] = present;
        }
    }
    for (std::size_t i = 0; i < staged.size(); ++i)
        if (specs.at(i).at("missing_policy").text() == "reject" &&
            seen[i].size() != source.genes.size())
            throw std::invalid_argument("reject missingness requires complete ambient gene "
                                        "coverage; declare zero_extension for a measured subset");
    // Admission is complete before mutation of the workspace's field registry.
    source.raw_sources.emplace_back(source.files.size()+1,std::span<const std::byte>{reinterpret_cast<const std::byte *>(config.bytes.data()),config.bytes.size()});
    source.raw_sources.emplace_back(source.files.size()+2,std::span<const std::byte>{reinterpret_cast<const std::byte *>(input.bytes.data()),input.bytes.size()});
    source.files.push_back(std::move(config.file));
    source.files.push_back(std::move(input.file));
    source.measurements.emplace("import:" + source.files.back().path, std::move(observations));
    for (auto &field : staged) {
        source.field_semantics.emplace(field.id, std::move(field.semantics));
        source.field_presence.emplace(field.id, std::move(field.present));
        source.fields.emplace(std::move(field.id), std::move(field.values));
    }
    refresh_source_digest(source);
    admit_workspace(source);
}
} // namespace gadget
