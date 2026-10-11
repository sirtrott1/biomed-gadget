#include <gadget/scientific.hpp>
#include <gadget/study.hpp>
#include <iostream>
#include <rfc/core/builder.hpp>
#include <rfc/storage/native_series_codec.hpp>

// A source producer for examples. Mathematical maps use canonical RexGraph
// materializers or explicitly authored coefficients; this is not an inference engine.
namespace {
using namespace gadget;
using namespace rfc;
using namespace rfc::storage;
using Matrix = numeric::ExactMatrix;
std::uint64_t id(const Json &j) { return std::stoull(j.text()); }
std::size_t size(const Json &j) {
    const auto n = j.integer();
    if (n < 0)
        throw std::invalid_argument("negative source extent");
    return static_cast<std::size_t>(n);
}
Json evidence(const std::string &bytes) {
    auto e = Json::object(), out = Json::array();
    e.set("document", Json::string("declarations"));
    e.set("start", Json::integer(0));
    e.set("end", Json::integer(static_cast<std::int64_t>(bytes.size())));
    e.set("text", Json::string(bytes));
    out.push(e);
    return out;
}
void document(Json &docs, const fs::path &root, const std::string &name, const std::string &path,
              const std::string &bytes, const std::string &format) {
    fs::create_directories((root / path).parent_path());
    write_text(root / path, bytes);
    auto d = Json::object();
    d.set("id", Json::string(name));
    d.set("path", Json::string(path));
    d.set("sha256", Json::string("sha256:" + sha256(bytes)));
    d.set("format", Json::string(format));
    docs.push(d);
}
TemporalEndpoint authored_source(const Json &spec) {
    core::RelationSignature signature;
    const auto type = signature.add_type("declared.relation", core::BoundaryPolicy::DeclaredExact);
    core::RelationalComplexBuilder builder(signature);
    const auto s = spec.at("source");
    for (std::size_t i = 0; i < s.at("vertices").size(); ++i)
        builder.declare_witness_identity(core::IdentityKey{id(s.at("vertices").at(i))});
    for (std::size_t i = 0; i < s.at("relations").size(); ++i) {
        const auto r = s.at("relations").at(i);
        std::vector<core::BoundaryOccurrenceInput> occurrences;
        for (std::size_t j = 0; j < r.at("occurrences").size(); ++j) {
            const auto o = r.at("occurrences").at(j);
            occurrences.push_back({core::IdentityKey{id(o.at("cell"))},
                                   o.at("exists").boolean(),
                                   o.at("negative").boolean(),
                                   exact_rational(o.at("magnitude").text()),
                                   {},
                                   core::OccurrenceKey{id(o.at("key"))}});
        }
        builder.add_grade_one(core::RelationKey{id(r.at("key"))}, type, occurrences);
    }
    if (s.has("higher"))
        for (std::size_t i = 0; i < s.at("higher").size(); ++i) {
            const auto r = s.at("higher").at(i);
            const core::Grade grade{static_cast<std::uint16_t>(size(r.at("grade")))};
            const auto current = builder.freeze();
            const StableSourceIndex index(current);
            std::vector<core::HigherIncidenceInput> occurrences;
            for (std::size_t j = 0; j < r.at("occurrences").size(); ++j) {
                const auto o = r.at("occurrences").at(j);
                const auto cell = index.cell(
                    {core::Grade{static_cast<std::uint16_t>(grade.value - 1)}, id(o.at("cell"))});
                if (!cell)
                    throw std::invalid_argument("higher attachment lower cell missing");
                occurrences.push_back({cell->index,
                                       o.at("exists").boolean(),
                                       o.at("negative").boolean(),
                                       exact_rational(o.at("magnitude").text()),
                                       {},
                                       core::OccurrenceKey{id(o.at("key"))}});
            }
            builder.add_higher(grade, core::RelationKey{id(r.at("key"))}, type, occurrences);
        }
    return TemporalEndpoint(builder.freeze(), id(spec.at("sequence")));
}
Matrix coefficients(const Json &j) {
    Matrix out(size(j.at("rows")), size(j.at("columns")));
    for (std::size_t i = 0; i < j.at("entries").size(); ++i) {
        const auto e = j.at("entries").at(i);
        out.set(size(e.at("row")), size(e.at("column")), exact_rational(e.at("value").text()));
    }
    return out;
}
void produce(const fs::path &out, const fs::path &spec_path,
             const std::optional<fs::path> &original) {
    if (fs::exists(out))
        throw std::invalid_argument("destination exists");
    const auto bytes = read_text(spec_path);
    const auto spec = Json::parse(bytes), spans = evidence(bytes);
    std::shared_ptr<const ScientificStudy> retained;
    if (original)
        retained = std::make_shared<const ScientificStudy>(*original);
    auto endpoint = retained ? *retained->endpoint() : authored_source(spec);
    auto d = retained ? retained->declaration() : Json::object();
    if (d.has("publications") && d.at("publications").size())
        throw std::invalid_argument("augment an actual leaf publication");
    fs::create_directories(out);
    auto docs = Json::array();
    if (retained) {
        for (std::size_t i = 0; i < d.at("documents").size(); ++i) {
            const auto doc = d.at("documents").at(i);
            document(docs, out, doc.at("id").text(), "original/document-" + std::to_string(i),
                     read_text(original->parent_path() / doc.at("path").text()),
                     doc.at("format").text());
        }
        document(docs, out, "origin-native", "original/native.rfcxs",
                 read_text(original->parent_path() / d.at("native").at("path").text()), "opaque");
        document(docs, out, "origin-declaration", "original/study.json", read_text(*original),
                 "utf8");
    }
    document(docs, out, "declarations", "declarations.json", bytes, "utf8");
    if (!retained) {
        auto cells = Json::array(), occurrences = Json::array(), contexts = Json::array(),
             metrics = Json::array();
        for (std::size_t g = 0; g < endpoint.source.grade_count(); ++g) {
            const core::Grade grade{static_cast<std::uint16_t>(g)};
            const auto &arena = endpoint.source.grade(grade);
            for (std::size_t i = 0; i < arena.cell_count(); ++i) {
                const auto k =
                    stable_cell_key(endpoint.source, {grade, static_cast<std::uint32_t>(i)});
                auto c = Json::object();
                c.set("grade", Json::integer(static_cast<std::int64_t>(g)));
                c.set("key", Json::string(std::to_string(k.key)));
                c.set("id",
                      Json::string("cell:" + std::to_string(g) + ":" + std::to_string(k.key)));
                c.set("kind", Json::string("authored_source"));
                c.set("label", c.at("id"));
                c.set("evidence", spans);
                cells.push(c);
            }
            for (const auto k : arena.occurrence_identities) {
                auto o = Json::object();
                o.set("grade", Json::integer(static_cast<std::int64_t>(g)));
                o.set("key", Json::string(std::to_string(k.value)));
                o.set("evidence", spans);
                occurrences.push(o);
            }
            auto m = Json::object();
            m.set("kind", Json::string("identity"));
            metrics.push(m);
        }
        auto context = Json::object();
        context.set("id", Json::string("all"));
        context.set("block", Json{});
        context.set("time", Json{});
        context.set("metrics", metrics);
        contexts.push(context);
        d.set("schema", Json::string("gadget.scientific.study.v1"));
        d.set(
            "design",
            Json::parse(
                R"({"question":"How do declared source, realization, operator and observation changes contribute?","population":"authored source","selection":"Every original coordinate and occurrence","covariates":[],"limitations":["Algebra qualification"]})"));
        d.set("cells", cells);
        d.set("occurrences", occurrences);
        d.set("contexts", contexts);
        d.set("fields", Json::array());
        d.set("field_sets", Json::array());
    }
    d.set("id", spec.at("id"));
    d.set("documents", docs);
    if (spec.has("carriers"))
        for (std::size_t i = 0; i < spec.at("carriers").size(); ++i) {
            const auto c = spec.at("carriers").at(i);
            endpoint.carriers.push_back(
                {id(c.at("key")),
                 c.has("grade") ? TemporalCarrierKind::SourceGrade : TemporalCarrierKind::Explicit,
                 core::Grade{static_cast<std::uint16_t>(c.has("grade") ? size(c.at("grade")) : 0)},
                 {},
                 c.has("coordinate_system") ? id(c.at("coordinate_system")) : 0,
                 numeric::IndexDomain::linear(size(c.at("size")))});
        }
    auto declarations = d.has("operators") ? d.at("operators") : Json::array();
    for (std::size_t i = 0; i < spec.at("operators").size(); ++i) {
        const auto op = spec.at("operators").at(i);
        const auto kind = op.at("role").text();
        TemporalTensor t;
        t.key = id(op.at("key"));
        t.input = id(op.at("input"));
        t.output = id(op.at("output"));
        if (kind != "resolving" && kind != "realization" && kind != "accession" && kind != "metric")
            throw std::invalid_argument("unknown declared response role");
        t.role = kind == "resolving"     ? TemporalTensorRole::SourceRealization
                 : kind == "realization" ? TemporalTensorRole::SourceRealization
                 : kind == "accession"   ? TemporalTensorRole::Accession
                                         : TemporalTensorRole::Metric;
        Matrix map;
        if (op.has("hodge_grade")) {
            std::vector<field::ExactMetric> metrics;
            const auto context = d.at("contexts").at(0);
            for (std::size_t g = 0; g < endpoint.source.grade_count(); ++g) {
                const auto m = context.at("metrics").at(g);
                if (m.at("kind").text() == "identity")
                    metrics.push_back(field::ExactMetric::identity(
                        endpoint.source.grade(core::Grade{static_cast<std::uint16_t>(g)})
                            .cell_count()));
                else {
                    Vector diagonal;
                    for (std::size_t j = 0; j < m.at("values").size(); ++j)
                        diagonal.push_back(exact_decimal(m.at("values").at(j).text()));
                    metrics.push_back(field::ExactMetric::diagonal(std::move(diagonal)));
                }
            }
            map = NativeStudy(endpoint.source, std::move(metrics))
                      .hodge(core::Grade{static_cast<std::uint16_t>(size(op.at("hodge_grade")))})
                      .materialize();
            if (op.has("identity_shift"))
                map = map.add(Matrix::identity(map.rows())
                                  .scale(exact_rational(op.at("identity_shift").text())));
        } else if (op.has("boundary_transpose_grade"))
            map = field::materialize_boundary_matrix(endpoint.source,
                                                     core::Grade{static_cast<std::uint16_t>(
                                                         size(op.at("boundary_transpose_grade")))})
                      .transpose();
        else
            map = coefficients(op.at("map"));
        t.interpretation = op.at("method").text();
        t.tensor = temporal_matrix_tensor(map);
        endpoint.tensors.push_back(t);
        auto metadata = Json::object();
        for (const auto *name : {"id", "key", "quantity", "units", "method"})
            metadata.set(name, op.at(name));
        metadata.set("evidence", spans);
        declarations.push(metadata);
    }
    d.set("operators", declarations);
    auto fields = d.at("fields");
    if (spec.has("fields"))
        for (std::size_t i = 0; i < spec.at("fields").size(); ++i) {
            const auto f = spec.at("fields").at(i);
            std::vector<numeric::SparseRationalEntry> values;
            std::vector<std::size_t> support;
            for (std::size_t j = 0; j < f.at("values").size(); ++j)
                if (!json_object_is_type(f.at("values").at(j).raw(), json_type_null)) {
                    support.push_back(j);
                    const auto q = exact_rational(f.at("values").at(j).text());
                    if (q != 0)
                        values.push_back({j, q});
                }
            TemporalTensor t;
            t.key = id(f.at("key"));
            t.role = TemporalTensorRole::Field;
            t.output = id(f.at("carrier"));
            t.interpretation = f.at("quantity").text();
            t.tensor = {
                numeric::IndexDomain::linear(f.at("values").size()),
                numeric::StructuralSupport::sparse(f.at("values").size(), std::move(support)),
                numeric::SparseRationalVector(f.at("values").size(), std::move(values))};
            endpoint.tensors.push_back(t);
            auto metadata = Json::object();
            for (const auto *name : {"id", "key", "quantity", "units", "value_type"})
                metadata.set(name, f.at(name));
            metadata.set("context", Json::string("all"));
            metadata.set("channel", Json::string("observation"));
            metadata.set("method", Json::string("declared source coordinate field"));
            metadata.set("evidence", spans);
            fields.push(metadata);
        }
    d.set("fields", fields);
    const NativeTemporalSeries native_history(endpoint, id(spec.at("namespace")), native_data_limits());
    const auto encoded = encode_native_series(native_history);
    const std::string archive(reinterpret_cast<const char *>(encoded.data()), encoded.size());
    write_text(out / "study.rfcxs", archive);
    auto native = Json::object();
    native.set("path", Json::string("study.rfcxs"));
    native.set("sha256", Json::string("sha256:" + sha256(archive)));
    native.set("source_namespace", spec.at("namespace"));
    native.set("cutoff", Json::string(std::to_string(endpoint.sequence)));
    d.set("native", native);
    write_text(out / "study.json", d.dump() + "\n");
    (void)ScientificStudy(out / "study.json");
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 3 && argc != 4)
            throw std::invalid_argument("Usage: gadget_declared_response_fixture SPEC.json "
                                        "NEW_DIRECTORY [ORIGINAL_STUDY.json]");
        produce(argv[2], argv[1], argc == 4 ? std::optional<fs::path>(argv[3]) : std::nullopt);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
