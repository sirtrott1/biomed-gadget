#include <gadget/scientific.hpp>
#include <iostream>
#include <rfc/core/builder.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace {
using namespace gadget;
using namespace rfc;
using namespace rfc::storage;
Json n(std::size_t value) { return Json::integer(static_cast<std::int64_t>(value)); }
Json span(const std::string &bytes) {
    auto e = Json::object(), out = Json::array();
    e.set("document", Json::string("observations"));
    e.set("start", n(0));
    e.set("end", n(bytes.size()));
    e.set("text", Json::string(bytes));
    out.push(e);
    return out;
}
TemporalEndpoint publication(const fs::path &root, std::size_t t) {
    fs::create_directories(root);
    std::string note =
        "Authored nonclinical publication " + std::to_string(t) +
        ". Original row identities and occurrence keys persist only in namespace 7001. "
        "The compared populations change. Repeated, cancelling and inactive occurrences are "
        "retained. "
        "R200 explicitly attaches to the cancelling R103; no face is inferred. "
        "The field x is authored, not an observed clinical trajectory. "
        "The incomplete field lacks its second coordinate; extension value 4 is a declared "
        "scenario.\n";
    note += t == 0 ? "Source rows: 10=1,20=2,30=3. Relation field: 101=1,102=2,103=3.\n"
                   : (t == 1 ? "Source rows: 20=3,30=4,40=5. Relation field: 101=2,102=2,103=3. "
                               "Partial row 30 is absent.\n"
                             : "Source rows: 20=7,40=9. Relation field: 101=7,103=9.\n");
    write_text(root / "observations.txt", note);
    const auto evidence = span(note);
    core::RelationSignature signature;
    const auto type = signature.add_type("authored.signed.relation");
    core::RelationalComplexBuilder builder(signature);
    const std::uint64_t a = t == 0 ? 10 : 20, b = t == 0 ? 20 : (t == 1 ? 30 : 40);
    builder.add_grade_one(core::RelationKey{101}, type,
                          {{core::IdentityKey{a},
                            true,
                            true,
                            Rational{t == 1 ? 2 : 1},
                            {},
                            core::OccurrenceKey{1001}},
                           {core::IdentityKey{b},
                            true,
                            false,
                            Rational{t == 1 ? 2 : 1},
                            {},
                            core::OccurrenceKey{1002}}});
    if (t != 2)
        builder.add_grade_one(core::RelationKey{102}, type,
                              {{core::IdentityKey{t == 0 ? 20u : 30u},
                                true,
                                t == 0,
                                Rational{1},
                                {},
                                core::OccurrenceKey{1003}},
                               {core::IdentityKey{t == 0 ? 30u : 40u},
                                true,
                                t != 0,
                                Rational{1},
                                {},
                                core::OccurrenceKey{1004}}});
    builder.add_grade_one(core::RelationKey{103}, type,
                          {{core::IdentityKey{20},
                            true,
                            false,
                            Rational{t == 1 ? 2 : 1},
                            {},
                            core::OccurrenceKey{1005}},
                           {core::IdentityKey{20},
                            true,
                            true,
                            Rational{t == 1 ? 2 : 1},
                            {},
                            core::OccurrenceKey{1006}},
                           {core::IdentityKey{t == 2 ? 40u : 30u},
                            t != 0,
                            false,
                            Rational{7},
                            {},
                            core::OccurrenceKey{1007}},
                           {core::IdentityKey{t == 2 ? 40u : 30u},
                            t != 0,
                            true,
                            Rational{7},
                            {},
                            core::OccurrenceKey{1008}}});
    builder.add_higher(core::Grade{2}, core::RelationKey{200}, type,
                       {{static_cast<std::uint32_t>(t == 2 ? 1 : 2),
                         true,
                         false,
                         Rational{1},
                         {},
                         core::OccurrenceKey{2001}}});
    TemporalEndpoint endpoint(builder.freeze(), t);
    auto cells = Json::array(), occurrences = Json::array(), metrics = Json::array();
    for (std::size_t g = 0; g < endpoint.source.grade_count(); ++g) {
        const auto grade = core::Grade{static_cast<std::uint16_t>(g)};
        const auto &arena = endpoint.source.grade(grade);
        endpoint.carriers.push_back({g + 1,
                                     TemporalCarrierKind::SourceGrade,
                                     grade,
                                     {},
                                     0,
                                     numeric::IndexDomain::linear(arena.cell_count())});
        auto metric = Json::object();
        metric.set("kind", Json::string("identity"));
        metrics.push(metric);
        for (std::size_t i = 0; i < arena.cell_count(); ++i) {
            const auto key =
                stable_cell_key(endpoint.source, {grade, static_cast<std::uint32_t>(i)});
            auto c = Json::object();
            c.set("grade", n(g));
            c.set("key", Json::string(std::to_string(key.key)));
            const auto id = std::string(g == 0 ? "row:" : "relation:") + std::to_string(key.key);
            c.set("id", Json::string(id));
            c.set("label", Json::string(id));
            c.set("kind", Json::string(g == 0 ? "authored_record" : "authored_relation"));
            c.set("evidence", evidence);
            cells.push(c);
        }
        for (auto key : arena.occurrence_identities) {
            auto o = Json::object();
            o.set("grade", n(g));
            o.set("key", Json::string(std::to_string(key.value)));
            o.set("evidence", evidence);
            occurrences.push(o);
        }
    }
    auto fields = Json::array();
    for (std::size_t f = 0; f < 3; ++f) {
        const auto grade = core::Grade{static_cast<std::uint16_t>(f == 1 ? 1 : 0)};
        const auto size = endpoint.source.grade(grade).cell_count();
        std::vector<numeric::SparseRationalEntry> entries;
        auto support = numeric::StructuralSupport::full(size);
        if (f == 2 && t == 1)
            support = numeric::StructuralSupport::sparse(size, {0, 2});
        const Vector values =
            grade.value == 1
                ? (t == 2 ? Vector{7, 9}
                          : Vector{Rational{static_cast<unsigned long>(t + 1)}, 2, 3})
                : (t == 0 ? Vector{1, 2, 3} : (t == 1 ? Vector{3, 4, 5} : Vector{7, 9}));
        for (std::size_t i = 0; i < size; ++i)
            if (support.contains(i) && values[i] != 0)
                entries.push_back({i, values[i]});
        TemporalTensor tensor;
        tensor.key = 10 + f;
        tensor.role = TemporalTensorRole::Field;
        tensor.output = grade.value + 1;
        tensor.interpretation = f == 1 ? "relation amplitude" : "authored scalar";
        tensor.tensor = {numeric::IndexDomain::linear(size), support,
                         numeric::SparseRationalVector(size, std::move(entries))};
        endpoint.tensors.push_back(tensor);
        auto d = Json::object();
        d.set("id", Json::string(f == 0 ? "x" : (f == 1 ? "r" : "partial")));
        d.set("key", Json::string(std::to_string(tensor.key)));
        d.set("context", Json::string("all"));
        d.set("channel", Json::string("observation"));
        d.set("quantity", Json::string(f == 1 ? "authored_relation" : "authored_scalar"));
        d.set("units", Json::string("dimensionless"));
        d.set("method", Json::string("authored rational fixture"));
        d.set("value_type", Json::string("quantitative"));
        d.set("evidence", evidence);
        fields.push(d);
    }
    auto contexts = Json::array(), context = Json::object(), time = Json::object();
    time.set("axis", Json::string("authored_visit"));
    time.set("position", Json::string(std::to_string(t)));
    time.set("unit", Json::string("day"));
    time.set("population", Json::string("authored population " + std::to_string(t)));
    context.set("id", Json::string("all"));
    context.set("block", Json{});
    context.set("time", time);
    context.set("metrics", metrics);
    contexts.push(context);
    auto bytes = encode_native_series(NativeTemporalSeries(endpoint, 7001));
    const std::string native(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    write_text(root / "study.rfcxs", native);
    auto archive = Json::object();
    archive.set("path", Json::string("study.rfcxs"));
    archive.set("sha256", Json::string("sha256:" + sha256(native)));
    archive.set("source_namespace", Json::string("7001"));
    auto design = Json::parse(
        R"({"question":"What changes under declared source and field evolution?","population":"authored changing records","selection":"Every retained original record","covariates":[],"limitations":["Constructed algebra qualification"]})");
    auto docs = Json::array(), doc = Json::object();
    doc.set("id", Json::string("observations"));
    doc.set("path", Json::string("observations.txt"));
    doc.set("sha256", Json::string("sha256:" + sha256(note)));
    doc.set("format", Json::string("utf8"));
    docs.push(doc);
    auto d = Json::object();
    d.set("schema", Json::string("gadget.scientific.study.v1"));
    d.set("id", Json::string("authored-publication-" + std::to_string(t)));
    d.set("native", archive);
    d.set("design", design);
    d.set("documents", docs);
    d.set("cells", cells);
    d.set("occurrences", occurrences);
    d.set("contexts", contexts);
    d.set("fields", fields);
    d.set("field_sets", Json::array());
    write_text(root / "study.json", d.dump() + "\n");
    return endpoint;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 2 && !(argc == 3 && std::string(argv[2]) == "--series"))
            throw std::invalid_argument(
                "Usage: gadget_source_history_fixture NEW_DIRECTORY [--series]");
        const gadget::fs::path root = argv[1];
        if (gadget::fs::exists(root))
            throw std::invalid_argument("destination exists");
        std::vector<rfc::storage::TemporalEndpoint> endpoints;
        for (std::size_t t = 0; t < 3; ++t)
            endpoints.push_back(publication(root / ("p" + std::to_string(t)), t));
        if (argc == 3) {
            rfc::storage::NativeTemporalSeries series(endpoints.front(), 7001);
            for (std::size_t t = 1; t < endpoints.size(); ++t)
                series.append_target(endpoints[t]);
            const auto bytes = rfc::storage::encode_native_series(series);
            const std::string archive(reinterpret_cast<const char *>(bytes.data()), bytes.size());
            for (std::size_t t = 0; t < endpoints.size(); ++t) {
                const auto dir = root / ("p" + std::to_string(t));
                gadget::write_text(dir / "history.rfcxs", archive);
                auto d = gadget::Json::parse(gadget::read_text(dir / "study.json"));
                d.at("native").set("path", gadget::Json::string("history.rfcxs"));
                d.at("native").set("sha256",
                                   gadget::Json::string("sha256:" + gadget::sha256(archive)));
                d.at("native").set("cutoff", gadget::Json::string(std::to_string(t)));
                gadget::write_text(dir / "study.json", d.dump() + "\n");
            }
        }
        std::cout << "Created three authored native publications\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
