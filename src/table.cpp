#include "scientific_detail.hpp"
#include <array>
#include <fstream>
#include <gadget/table.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/storage/mapped_file.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <set>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc::storage;
Json span(const std::string &document, std::string_view bytes, std::size_t start, std::size_t end,
          bool quoted) {
    auto item = Json::object(), result = Json::array();
    item.set("document", Json::string(document));
    item.set("start", cardinality(start));
    item.set("end", cardinality(end));
    if (quoted)
        item.set("text", Json::string(std::string(bytes.substr(start, end - start))));
    result.push(item);
    return result;
}
struct Row {
    std::size_t start{}, end{};
};
struct Group {
    std::string column, value;
    std::uint64_t key{};
    std::size_t axis{};
    std::vector<std::size_t> rows;
};
} // namespace
void import_table(const fs::path &specification, const fs::path &destination) {
    require(fs::is_regular_file(fs::symlink_status(specification)),
            "table specification must be regular");
    require(!fs::exists(destination), "table destination already exists");
    const auto spec_path = fs::canonical(specification);
    const auto spec_bytes = read_text(spec_path);
    const auto spec = Json::parse(spec_bytes);
    members(spec,
            {"schema", "id", "source_namespace", "source", "dictionary", "columns",
             "missing_tokens", "relation_columns", "question", "population", "selection",
             "limitations", "source_url", "attribution", "license"},
            {"identity_column", "contrasts", "arithmetic", "evidence_profile", "intake"});
    require(spec.at("schema").text() == "gadget.table.import.v1",
            "unsupported table import schema");
    const auto arithmetic = spec.has("arithmetic") ? text(spec.at("arithmetic")) : "exact_decimal";
    require(arithmetic == "exact_decimal" || arithmetic == "exact_rational",
            "table arithmetic profile");
    const auto id = text(spec.at("id"));
    for (const auto *name :
         {"question", "population", "selection", "source_url", "attribution", "license"})
        text(spec.at(name));
    const auto namespace_key = key(spec.at("source_namespace"));
    const auto input = spec.at("source"), dictionary_record = spec.at("dictionary");
    members(input, {"path", "sha256", "delimiter", "header"});
    members(dictionary_record, {"path", "sha256"});
    const auto profile =
        spec.has("evidence_profile") ? text(spec.at("evidence_profile")) : "quoted_v1";
    require(profile == "quoted_v1" || profile == "document_spans_v2", "table evidence profile");
    const bool quoted = profile == "quoted_v1";
    std::size_t max_source_bytes = source_byte_limit(), max_record_bytes = source_byte_limit(),
                max_token_offsets = std::numeric_limits<std::size_t>::max();
    if (spec.has("intake")) {
        const auto limits = spec.at("intake");
        members(limits, {}, {"max_source_bytes", "max_record_bytes", "max_token_offsets"});
        if (limits.has("max_source_bytes"))
            max_source_bytes =
                count(limits.at("max_source_bytes"), max_source_bytes, "source bytes");
        if (limits.has("max_record_bytes"))
            max_record_bytes =
                count(limits.at("max_record_bytes"), max_record_bytes, "record bytes");
        if (limits.has("max_token_offsets"))
            max_token_offsets =
                count(limits.at("max_token_offsets"), max_token_offsets, "token offsets");
    }
    const auto source_path = local_artifact_path(spec_path.parent_path(), input);
    const auto source_bytes = fs::file_size(source_path);
    require(source_bytes > 0 && source_bytes <= max_source_bytes, "table source byte budget");
    const auto parent = fs::absolute(destination).parent_path();
    fs::create_directories(parent);
    const auto stage = private_directory(parent, ".gadget-table-");
    struct Cleanup {
        fs::path path;
        ~Cleanup() {
            std::error_code error;
            fs::remove_all(path, error);
        }
    } cleanup{stage};
    // Copy once to a private retained image before parsing. mmap is a byte
    // backend, not an immutable view of a caller owned file being edited.
    publish_file_exact(stage / "table.txt", FilePublication::CreateNew, [&](std::ostream &out) {
        std::ifstream stream(source_path, std::ios::binary);
        require(static_cast<bool>(stream), "table source open");
        // A transfer chunk is not a source size limit. Read exactly the admitted
        // image and reject growth/truncation before a larger file can be retained.
        std::array<char, 1024 * 1024> buffer{};
        for (std::uintmax_t copied = 0; copied < source_bytes;) {
            const auto n = static_cast<std::size_t>(
                std::min<std::uintmax_t>(source_bytes - copied, buffer.size()));
            stream.read(buffer.data(), static_cast<std::streamsize>(n));
            require(stream.gcount() == static_cast<std::streamsize>(n), "table source truncated");
            write_stream_exact(out, {reinterpret_cast<const std::byte *>(buffer.data()), n});
            copied += n;
        }
        require(stream.peek() == std::char_traits<char>::eof() && !stream.bad(),
                "table source grew or became unreadable");
    });
    ReadOnlyMappedFile image(stage / "table.txt");
    require(image.size() == source_bytes && image.size() <= max_source_bytes,
            "table source changed size during retention");
    const std::string_view bytes(reinterpret_cast<const char *>(image.data()), image.size());
    require(input.at("sha256").text() == "sha256:" + sha256(bytes),
            "table source SHA-256 mismatch");
    const auto dictionary = local_artifact(spec_path.parent_path(), dictionary_record);
    const auto delimiter = text(input.at("delimiter"), 1)[0];
    require(delimiter == ',' || delimiter == '\t', "table delimiter must be comma or tab");
    const auto columns = spec.at("columns");
    require(columns.size() > 0, "table requires declared columns");
    std::map<std::string, std::size_t> indices;
    for (std::size_t c = 0; c < columns.size(); ++c) {
        const auto column = columns.at(c);
        members(column, {"id", "channel", "quantity", "units", "value_type", "dictionary_start",
                         "dictionary_end"});
        const auto name = text(column.at("id"));
        require(column.at("value_type").text() == "categorical" ||
                    column.at("value_type").text() == "quantitative" ||
                    column.at("value_type").text() == "indicator",
                "table value type");
        require(indices.emplace(name, c).second, "duplicate table column");
        const auto start =
            count(column.at("dictionary_start"), dictionary.size(), "dictionary offset");
        const auto end = count(column.at("dictionary_end"), dictionary.size(), "dictionary offset");
        require(start < end, "dictionary evidence must be nonempty");
    }
    std::set<std::string> missing;
    for (std::size_t i = 0; i < spec.at("missing_tokens").size(); ++i)
        require(missing.insert(spec.at("missing_tokens").at(i).text()).second,
                "duplicate missing token");
    std::vector<Row> rows;
    std::vector<std::size_t> token_ends;
    std::size_t largest_record{};
    bool header = input.at("header").boolean();
    for (std::size_t start = 0; start < bytes.size();) {
        auto end = bytes.find('\n', start);
        if (end == std::string_view::npos)
            end = bytes.size();
        const auto line = bytes.substr(start, end - start);
        require(line.size() <= max_record_bytes, "table record byte budget");
        largest_record = std::max(largest_record, line.size());
        require(!line.empty() && line.find('"') == std::string_view::npos &&
                    line.find('\0') == std::string_view::npos,
                "only nonempty unquoted records are admitted");
        std::size_t position = 0, c = 0;
        for (;;) {
            auto stop = line.find(delimiter, position);
            if (stop == std::string_view::npos)
                stop = line.size();
            require(c < columns.size(), "row has undeclared columns");
            if (header)
                require(trim(std::string(line.substr(position, stop - position))) ==
                            columns.at(c).at("id").text(),
                        "table header differs from declared order");
            else {
                require(token_ends.size() < max_token_offsets, "table token offset budget");
                token_ends.push_back(start + stop);
            }
            ++c;
            if (stop == line.size())
                break;
            position = stop + 1;
        }
        require(c == columns.size(), "row has missing columns");
        if (header)
            header = false;
        else
            rows.push_back({start, end});
        start = end == bytes.size() ? end : end + 1;
    }
    // The only retained lexical index contains original byte offsets, not
    // copies of every token. Each value is decoded when its native field needs it.
    const auto token = [&](std::size_t r, std::size_t c) {
        const auto i =
            rfc::machine::checked_size_add(rfc::machine::checked_size_mul(r, columns.size()), c);
        const auto start = c == 0 ? rows.at(r).start : token_ends.at(i - 1) + 1;
        return trim(std::string(bytes.substr(start, token_ends.at(i) - start)));
    };
    require(!rows.empty(), "table has no observed records");
    std::vector<std::string> row_identities;
    if (spec.has("identity_column")) {
        const auto column = text(spec.at("identity_column"));
        require(indices.contains(column), "unknown identity column");
        require(columns.at(indices.at(column)).at("value_type").text() == "categorical",
                "source identities must be lexical, not linear measurements");
        std::set<std::string> seen;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const auto value = token(r, indices.at(column));
            require(!missing.contains(value) && !value.empty() && seen.insert(value).second,
                    "source identities must be reported and unique; repeated samples need explicit "
                    "sample identities");
            row_identities.push_back(value);
        }
    }
    const auto relation_columns = spec.at("relation_columns");
    require(relation_columns.size() > 0, "table requires declared categorical relation columns");
    std::vector<Group> groups;
    std::set<std::string> relation_names;
    for (std::size_t axis = 0; axis < relation_columns.size(); ++axis) {
        const auto name = relation_columns.at(axis).text();
        require(name != "all", "relation column name all is reserved for the full context");
        require(indices.contains(name) && relation_names.insert(name).second,
                "unknown/duplicate relation column");
        require(columns.at(indices.at(name)).at("value_type").text() == "categorical",
                "relation membership column must be explicitly categorical");
        std::map<std::string, std::vector<std::size_t>> membership;
        for (std::size_t r = 0; r < rows.size(); ++r)
            if (!missing.contains(token(r, indices.at(name))))
                membership[token(r, indices.at(name))].push_back(r);
        for (auto &[value, members] : membership) {
            groups.push_back({name, value, rfc::machine::checked_size_add(groups.size(), 1), axis,
                              std::move(members)});
        }
    }
    require(!groups.empty(), "no reported categorical membership");
    rfc::core::RelationSignature signature;
    const auto type = signature.add_type("reported.lexical.category.membership");
    rfc::core::RelationalComplexBuilder builder(signature);
    for (std::size_t r = 0; r < rows.size(); ++r)
        builder.declare_witness_identity(r + 1);
    std::map<std::uint64_t, std::size_t> occurrence_rows;
    std::uint64_t occurrence = 1;
    for (const auto &group : groups) {
        std::vector<rfc::core::BoundaryOccurrenceInput> memberships;
        for (const auto r : group.rows) {
            occurrence_rows.emplace(occurrence, r);
            memberships.push_back(
                {r + 1, true, false, Rational{1}, {}, rfc::core::OccurrenceKey{occurrence++}});
        }
        builder.add_grade_one(rfc::core::RelationKey{group.key}, type, memberships);
    }
    TemporalEndpoint endpoint(builder.freeze());
    // Keys belong to their declared namespaces. Cell keys are grade relative;
    // carrier and tensor keys each use one complete source ordered namespace.
    std::map<StableCellKey, std::uint64_t> stalk_carriers;
    for (std::size_t k = 0; k < endpoint.source.grade_count(); ++k) {
        const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(k)};
        endpoint.carriers.push_back(
            {rfc::machine::checked_size_add(endpoint.carriers.size(), 1),
             TemporalCarrierKind::SourceGrade,
             grade,
             {},
             0,
             rfc::numeric::IndexDomain::linear(endpoint.source.grade(grade).cell_count())});
    }
    auto cells = Json::array(), occurrences = Json::array(), stalks = Json::array();
    const auto spec_evidence = span("import_spec", spec_bytes, 0, spec_bytes.size(), quoted);
    for (std::size_t k = 0; k < 2; ++k) {
        const auto grade = rfc::core::Grade{static_cast<std::uint16_t>(k)};
        const auto &arena = endpoint.source.grade(grade);
        for (std::size_t i = 0; i < arena.cell_count(); ++i) {
            const auto cell =
                stable_cell_key(endpoint.source, {grade, static_cast<std::uint32_t>(i)});
            auto item = Json::object();
            item.set("grade", cardinality(k));
            item.set("key", Json::string(std::to_string(cell.key)));
            const auto identity =
                k == 0 ? (row_identities.empty() ? "record:" + id + ":" + std::to_string(cell.key)
                                                 : row_identities.at(cell.key - 1))
                       : "category:" + groups.at(i).column + ":" + groups.at(i).value;
            item.set("id", Json::string(identity));
            item.set("label", Json::string(identity));
            item.set("kind", Json::string(k == 0 ? "source_record" : "reported_category"));
            item.set("evidence", k == 0 ? span("table", bytes, rows.at(cell.key - 1).start,
                                               rows.at(cell.key - 1).end, quoted)
                                        : spec_evidence);
            cells.push(item);
            const auto carrier = rfc::machine::checked_size_add(endpoint.carriers.size(), 1);
            stalk_carriers.emplace(cell, carrier);
            endpoint.carriers.push_back(
                {carrier, TemporalCarrierKind::CellStalk, grade, cell, 900,
                 rfc::numeric::IndexDomain::linear(relation_columns.size())});
            stalks.push(Json::string(std::to_string(carrier)));
        }
        for (std::size_t p = 0; p < arena.incidence_count(); ++p) {
            const auto key = arena.occurrence_identities[p].value;
            const auto &row = rows.at(occurrence_rows.at(key));
            auto item = Json::object();
            item.set("grade", cardinality(k));
            item.set("key", Json::string(std::to_string(key)));
            item.set("evidence", span("table", bytes, row.start, row.end, quoted));
            occurrences.push(item);
        }
    }
    auto contexts = Json::array();
    std::vector<std::string> context_names{"all"};
    for (std::size_t i = 0; i < relation_columns.size(); ++i)
        context_names.push_back(relation_columns.at(i).text());
    for (std::size_t c = 0; c < context_names.size(); ++c) {
        auto context = Json::object();
        context.set("id", Json::string(context_names[c]));
        context.set("time", Json{});
        context.set("metrics", Json::parse(R"([{"kind":"identity"},{"kind":"identity"}])"));
        context.set("block", c == 0 ? Json{} : Json::string(std::to_string(c)));
        if (c) {
            LocalBlockDefinition block;
            block.block_key = c;
            block.kind = LocalBlockKind::StrictSubcomplex;
            for (const auto &group : groups)
                if (group.column == context_names[c]) {
                    block.cells.push_back({{1}, group.key});
                    for (const auto row : group.rows)
                        block.cells.push_back({{0}, row + 1});
                }
            std::sort(block.cells.begin(), block.cells.end());
            block.cells.erase(std::unique(block.cells.begin(), block.cells.end()),
                              block.cells.end());
            require(!block.cells.empty(), "relation context has no retained observations");
            endpoint.local_blocks.push_back(std::move(block));
        }
        contexts.push(context);
    }
    auto fields = Json::array();
    const auto declaration_field = [&](const Json &column, std::uint64_t tensor,
                                       const std::string &context, const Json &evidence) {
        auto field = Json::object();
        field.set("id", Json::string(context + "." + column.at("id").text()));
        field.set("key", Json::string(std::to_string(tensor)));
        field.set("context", Json::string(context));
        for (const auto &name : {"channel", "quantity", "units", "value_type"})
            field.set(name, column.at(name));
        field.set(
            "method",
            column.has("method")
                ? column.at("method")
                : Json::string("unchanged reported token -> exact rational; declared missing "
                               "tokens preserve absence; context aliases share one native tensor"));
        field.set("evidence", evidence);
        fields.push(field);
    };
    for (std::size_t c = 0; c < columns.size(); ++c) {
        if (spec.has("identity_column") &&
            columns.at(c).at("id").text() == spec.at("identity_column").text())
            continue;
        std::vector<std::size_t> support;
        std::vector<rfc::numeric::SparseRationalEntry> values;
        std::map<std::string, Rational> codes;
        if (columns.at(c).at("value_type").text() == "categorical") {
            for (std::size_t r = 0; r < rows.size(); ++r)
                if (!missing.contains(token(r, c)))
                    codes.emplace(token(r, c), Rational{0});
            std::size_t code = 0;
            for (auto &[label, value] : codes) {
                (void)label;
                value = Rational{++code};
            }
        }
        for (std::size_t r = 0; r < rows.size(); ++r)
            if (!missing.contains(token(r, c))) {
                const Rational value =
                    codes.empty() ? (arithmetic == "exact_rational" ? exact_rational(token(r, c))
                                                                    : exact_decimal(token(r, c)))
                                  : codes.at(token(r, c));
                support.push_back(r);
                if (value != 0)
                    values.push_back({r, value});
            }
        TemporalTensor tensor;
        const auto tensor_key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
        tensor.key = tensor_key;
        tensor.role = TemporalTensorRole::Field;
        tensor.output = 1;
        tensor.interpretation = columns.at(c).at("id").text();
        tensor.tensor = {rfc::numeric::IndexDomain::linear(rows.size()),
                         rfc::numeric::StructuralSupport::sparse(rows.size(), support),
                         rfc::numeric::SparseRationalVector(rows.size(), std::move(values))};
        endpoint.tensors.push_back(std::move(tensor));
        auto col = Json::parse(columns.at(c).dump());
        if (!codes.empty()) {
            auto codebook = Json::object();
            for (const auto &[label, value] : codes)
                codebook.set(label, Json::string(value.get_str()));
            col.set("method",
                    Json::string("Lexical category codebook (no ordinal or linear meaning): " +
                                 codebook.dump()));
        }
        const auto evidence =
            span("dictionary", dictionary,
                 count(col.at("dictionary_start"), dictionary.size(), "dictionary start"),
                 count(col.at("dictionary_end"), dictionary.size(), "dictionary end"), quoted);
        for (const auto &context : context_names)
            declaration_field(col, tensor_key, context, evidence);
    }
    if (spec.has("contrasts")) {
        const auto contrasts = spec.at("contrasts");
        source_array(contrasts, "contrasts");
        std::set<std::string> contrast_ids;
        for (std::size_t i = 0; i < contrasts.size(); ++i) {
            const auto contrast = contrasts.at(i);
            members(contrast, {"id", "left", "right"});
            const auto id_ = text(contrast.at("id")), left = text(contrast.at("left")),
                       right = text(contrast.at("right"));
            require(!indices.contains(id_) && contrast_ids.insert(id_).second &&
                        indices.contains(left) && indices.contains(right),
                    "contrast field identities");
            auto a = columns.at(indices.at(left)), b = columns.at(indices.at(right));
            require(a.at("value_type").text() == "quantitative" &&
                        b.at("value_type").text() == "quantitative" &&
                        a.at("units").text() == b.at("units").text(),
                    "contrast requires quantitative fields with identical units");
            std::vector<std::size_t> support;
            std::vector<rfc::numeric::SparseRationalEntry> values;
            for (std::size_t r = 0; r < rows.size(); ++r)
                if (!missing.contains(token(r, indices.at(left))) &&
                    !missing.contains(token(r, indices.at(right)))) {
                    support.push_back(r);
                    const auto parse = [&](const std::string &token) {
                        return arithmetic == "exact_rational" ? exact_rational(token)
                                                              : exact_decimal(token);
                    };
                    Rational v =
                        parse(token(r, indices.at(left))) - parse(token(r, indices.at(right)));
                    v.canonicalize();
                    if (v != 0)
                        values.push_back({r, v});
                }
            TemporalTensor tensor;
            const auto tensor_key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
            tensor.key = tensor_key;
            tensor.role = TemporalTensorRole::Field;
            tensor.output = 1;
            tensor.interpretation = id_;
            tensor.tensor = {rfc::numeric::IndexDomain::linear(rows.size()),
                             rfc::numeric::StructuralSupport::sparse(rows.size(), support),
                             rfc::numeric::SparseRationalVector(rows.size(), std::move(values))};
            endpoint.tensors.push_back(std::move(tensor));
            auto col = Json::parse(a.dump());
            col.set("id", Json::string(id_));
            col.set("quantity", Json::string("declared_difference:" + left + "-" + right));
            col.set("method",
                    Json::string("Native exact rational difference on joint reported support: " +
                                 left + " minus " + right));
            for (const auto &context : context_names)
                declaration_field(col, tensor_key, context, spec_evidence);
        }
    }
    auto set_fields = Json::array();
    for (std::size_t i = 0; i < groups.size(); ++i) {
        const auto &group = groups[i];
        TemporalTensor tensor;
        const auto tensor_key = rfc::machine::checked_size_add(endpoint.tensors.size(), 1);
        tensor.key = tensor_key;
        tensor.role = TemporalTensorRole::Field;
        tensor.output = stalk_carriers.at({{1}, group.key});
        tensor.owner.kind = TemporalOwnerKind::Cell;
        tensor.owner.cell = {{1}, group.key};
        tensor.interpretation = "declared categorical membership basis";
        tensor.tensor = {rfc::numeric::IndexDomain::linear(relation_columns.size()),
                         rfc::numeric::StructuralSupport::full(relation_columns.size()),
                         rfc::numeric::SparseRationalVector(relation_columns.size(),
                                                            {{group.axis, Rational{1}}})};
        endpoint.tensors.push_back(std::move(tensor));
        auto column = Json::parse(
            R"({"channel":"observation","quantity":"declared_category_membership_basis","units":"dimensionless","value_type":"indicator"})");
        column.set("id", Json::string("membership." + std::to_string(group.key)));
        declaration_field(column, tensor_key, "all", spec_evidence);
        set_fields.push(Json::string("all.membership." + std::to_string(group.key)));
    }
    auto sets = Json::array(), set = Json::object();
    set.set("id", Json::string("category.memberships"));
    set.set("grade", Json::integer(1));
    set.set("fields", set_fields);
    set.set("carriers", stalks);
    set.set("maps", Json::string("constant"));
    set.set("coordinate_system", Json::string("900"));
    sets.push(set);
    const auto archive = encode_native_series(
        NativeTemporalSeries(endpoint, namespace_key, source_codec_limits(source_byte_limit())));
    const std::string_view native(reinterpret_cast<const char *>(archive.data()), archive.size());
    auto receipt = Json::object();
    receipt.set("schema", Json::string("gadget.table.receipt.v1"));
    receipt.set("producer", Json::string("gadget.table.import"));
    receipt.set("producer_version", Json::integer(3));
    receipt.set("evidence_profile", Json::string(profile));
    auto intake = Json::object();
    intake.set("source_bytes", cardinality(image.size()));
    intake.set("largest_record_bytes", cardinality(largest_record));
    intake.set("token_offsets", cardinality(token_ends.size()));
    intake.set("retained_index_bytes",
               cardinality(rfc::machine::checked_size_add(
                   rfc::machine::checked_size_mul(rows.size(), sizeof(Row)),
                   rfc::machine::checked_size_mul(token_ends.size(), sizeof(std::size_t)))));
    intake.set("limits", spec.has("intake") ? spec.at("intake") : Json::object());
    intake.set("byte_backend", Json::string("rexgraph_read_only_mapped_file"));
    receipt.set("intake", intake);
    receipt.set("key_policy", Json::string("positive_source_order_per_namespace"));
    receipt.set("source_namespace", Json::string(std::to_string(namespace_key)));
    receipt.set("specification_sha256", Json::string("sha256:" + sha256(spec_bytes)));
    receipt.set("table_sha256", input.at("sha256"));
    receipt.set("dictionary_sha256", dictionary_record.at("sha256"));
    receipt.set("native_sha256", Json::string("sha256:" + sha256(native)));
    receipt.set("rows", cardinality(rows.size()));
    receipt.set("columns", cardinality(columns.size()));
    receipt.set("categories", cardinality(groups.size()));
    receipt.set("carriers", cardinality(endpoint.carriers.size()));
    receipt.set("tensors", cardinality(endpoint.tensors.size()));
    receipt = seal(receipt, "receipt_sha256");
    const auto receipt_bytes = receipt.dump() + "\n";
    auto declaration = Json::object();
    declaration.set("schema", Json::string(quoted ? "gadget.scientific.study.v1"
                                                  : "gadget.scientific.study.v2"));
    declaration.set("id", Json::string(id));
    auto native_record = Json::object();
    native_record.set("path", Json::string("source.rfcxs"));
    native_record.set("sha256", Json::string("sha256:" + sha256(native)));
    native_record.set("source_namespace", Json::string(std::to_string(namespace_key)));
    declaration.set("native", native_record);
    auto design = Json::object();
    for (const auto &name : {"question", "population", "selection", "limitations"})
        design.set(name, spec.at(name));
    auto covariates = Json::array();
    for (std::size_t c = 0; c < columns.size(); ++c)
        if (columns.at(c).at("channel").text() == "covariate" &&
            (!spec.has("identity_column") ||
             columns.at(c).at("id").text() != spec.at("identity_column").text())) {
            auto covariate = Json::object();
            covariate.set("name", columns.at(c).at("id"));
            covariate.set("description", columns.at(c).at("quantity"));
            covariate.set("field", Json::string("all." + columns.at(c).at("id").text()));
            covariate.set("evidence", spec_evidence);
            covariates.push(covariate);
        }
    design.set("covariates", covariates);
    declaration.set("design", design);
    auto documents = Json::array();
    const auto doc = [&](const char *id_, const char *path_, std::string_view content) {
        auto d = Json::object();
        d.set("id", Json::string(id_));
        d.set("path", Json::string(path_));
        d.set("format", Json::string("utf8"));
        d.set("sha256", Json::string("sha256:" + sha256(content)));
        documents.push(d);
    };
    doc("table", "table.txt", bytes);
    doc("dictionary", "dictionary.txt", dictionary);
    doc("import_spec", "import-spec.json", spec_bytes);
    doc("import_receipt", "import-receipt.json", receipt_bytes);
    declaration.set("documents", documents);
    declaration.set("cells", cells);
    declaration.set("occurrences", occurrences);
    declaration.set("contexts", contexts);
    declaration.set("fields", fields);
    declaration.set("field_sets", sets);
    write_text(stage / "dictionary.txt", dictionary);
    write_text(stage / "import-spec.json", spec_bytes);
    write_text(stage / "import-receipt.json", receipt_bytes);
    write_file_exact(stage / "source.rfcxs", archive, FilePublication::CreateNew);
    write_text(stage / "study.json", declaration.dump() + "\n");
    // The lexical index is no longer needed once all complete native tensors
    // and semantic spans have been emitted. Its budget is separate from native
    // source admission, which still owns the whole complex and exact fields.
    std::vector<Row>().swap(rows);
    std::vector<std::size_t>().swap(token_ends);
    declaration = Json{};
    cells = Json{};
    occurrences = Json{};
    fields = Json{};
    contexts = Json{};
    stalks = Json{};
    set_fields = Json{};
    set = Json{};
    sets = Json{};
    {
        ScientificStudy admitted(stage / "study.json");
        write_text(stage / "admission.json", admitted.catalog().dump() + "\n");
        admitted.write_native(stage/"study.rfcxs");
    }
    require(fs::create_directory(destination), "table destination creation conflict");
    try {
        for (const auto &entry : fs::directory_iterator(stage))
            fs::rename(entry.path(), destination / entry.path().filename());
    } catch (...) {
        std::error_code ignored;
        fs::remove_all(destination, ignored);
        throw;
    }
}
} // namespace gadget
