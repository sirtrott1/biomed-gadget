#pragma once
#include <gadget/scientific.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/storage/native_series_codec.hpp>
#include <tuple>

namespace gadget::example {
// Synthetic arithmetic/evidence fixture. Every higher assertion is explicitly
// declared here; none is inferred from overlap or a biomedical model.
inline void scientific_study(const fs::path &root) {
    using namespace rfc;
    using namespace rfc::storage;
    if (fs::exists(root))
        throw std::invalid_argument("example destination already exists");
    fs::create_directories(root);
    const std::string note =
        "Synthetic nonclinical example. A, B and C are declared gene labels. β is UTF-8.\n"
        "Contrast R101: -A +B. Contrast R102: +A -A -B +C; inactive C with magnitude 7.\n"
        "R103: +A -A retains a cancelling column. R201 explicitly attaches to R103.\n"
        "baseline and followup use the same synthetic population at days 0 and 2.\n"
        "Cell-stalk coordinates 900 name two synthetic assay amplitudes, with declared identity "
        "maps.\n"
        "Batch codes are categorical and cannot be treated as linear measurements.\n";
    const std::string observations =
        "context,gene,G,E,T,batch,partial\n"
        "baseline,A,1,0,0,0,0\nbaseline,B,0,2,0,1,2\nbaseline,C,1,0,1,1,NA\n"
        "followup,A,1,1,1,0,NA\nfollowup,B,1,3,0,1,NA\nfollowup,C,1,2,1,1,NA\n";
    write_text(root / "design.txt", note);
    write_text(root / "observations.csv", observations);
    auto span = [](const std::string &document, const std::string &bytes) {
        auto evidence = Json::array();
        auto item = Json::object();
        item.set("document", Json::string(document));
        item.set("start", Json::integer(0));
        item.set("end", Json::integer(static_cast<std::int64_t>(bytes.size())));
        item.set("text", Json::string(bytes));
        evidence.push(item);
        return evidence;
    };
    const auto source_evidence = span("design", note);
    const auto measured_evidence = span("observations", observations);
    core::RelationSignature signature;
    const auto contrast = signature.add_type("explicit.synthetic.contrast");
    core::RelationalComplexBuilder builder(signature);
    builder.add_grade_one(
        core::RelationKey{101}, contrast,
        {{core::IdentityKey{10}, true, true, Rational{1}, {}, core::OccurrenceKey{1001}},
         {core::IdentityKey{20}, true, false, Rational{1}, {}, core::OccurrenceKey{1002}}});
    builder.add_grade_one(
        core::RelationKey{102}, contrast,
        {{core::IdentityKey{10}, true, false, Rational{1}, {}, core::OccurrenceKey{1003}},
         {core::IdentityKey{10}, true, true, Rational{1}, {}, core::OccurrenceKey{1004}},
         {core::IdentityKey{20}, true, true, Rational{1}, {}, core::OccurrenceKey{1005}},
         {core::IdentityKey{30}, true, false, Rational{1}, {}, core::OccurrenceKey{1006}},
         {core::IdentityKey{30}, false, false, Rational{7}, {}, core::OccurrenceKey{1007}}});
    builder.add_grade_one(
        core::RelationKey{103}, contrast,
        {{core::IdentityKey{10}, true, false, Rational{1}, {}, core::OccurrenceKey{1008}},
         {core::IdentityKey{10}, true, true, Rational{1}, {}, core::OccurrenceKey{1009}}});
    builder.add_higher(core::Grade{2}, core::RelationKey{201}, contrast,
                       {{2, true, false, Rational{1}, {}, core::OccurrenceKey{2001}}});
    TemporalEndpoint endpoint(builder.freeze());
    endpoint.carriers.push_back({1,
                                 TemporalCarrierKind::SourceGrade,
                                 core::Grade{0},
                                 {},
                                 0,
                                 numeric::IndexDomain::linear(3)});
    endpoint.carriers.push_back({2,
                                 TemporalCarrierKind::SourceGrade,
                                 core::Grade{1},
                                 {},
                                 0,
                                 numeric::IndexDomain::linear(3)});
    endpoint.carriers.push_back({3,
                                 TemporalCarrierKind::SourceGrade,
                                 core::Grade{2},
                                 {},
                                 0,
                                 numeric::IndexDomain::linear(1)});
    auto cells = Json::array(), occurrences = Json::array(), stalks = Json::array();
    for (std::size_t k = 0; k < endpoint.source.grade_count(); ++k) {
        const auto grade = core::Grade{static_cast<std::uint16_t>(k)};
        const auto &arena = endpoint.source.grade(grade);
        for (std::size_t i = 0; i < arena.cell_count(); ++i) {
            const auto cell =
                stable_cell_key(endpoint.source, {grade, static_cast<std::uint32_t>(i)});
            const auto id = k == 0 ? "gene:" + std::string(1, "ABC"[i])
                                   : "assertion:R" + std::to_string(cell.key);
            auto item = Json::object();
            item.set("grade", Json::integer(static_cast<std::int64_t>(k)));
            item.set("key", Json::string(std::to_string(cell.key)));
            item.set("id", Json::string(id));
            item.set("label", Json::string(id));
            item.set("kind", Json::string(k == 0 ? "gene" : "synthetic_assertion"));
            item.set("evidence", source_evidence);
            cells.push(item);
            const auto carrier_key = 10000 + cell.key;
            endpoint.carriers.push_back({carrier_key, TemporalCarrierKind::CellStalk, grade, cell,
                                         900, numeric::IndexDomain::linear(2)});
            stalks.push(Json::string(std::to_string(carrier_key)));
        }
        for (std::size_t p = 0; p < arena.incidence_count(); ++p) {
            auto item = Json::object();
            item.set("grade", Json::integer(static_cast<std::int64_t>(k)));
            item.set("key", Json::string(std::to_string(arena.occurrence_identities[p].value)));
            item.set("evidence", source_evidence);
            occurrences.push(item);
        }
    }
    for (unsigned which = 1; which <= 2; ++which) {
        LocalBlockDefinition block;
        block.block_key = which;
        block.kind = LocalBlockKind::StrictSubcomplex;
        block.cells = {{{0}, 10}, {{0}, 20}, {{0}, 30}, {{1}, 101}, {{1}, 103}, {{2}, 201}};
        if (which == 2)
            block.cells.push_back({{1}, 102});
        endpoint.local_blocks.push_back(std::move(block));
    }
    auto fields = Json::array();
    const auto add = [&](std::string id, std::uint64_t tensor_key, std::string context,
                         std::string channel, Vector values, std::uint64_t carrier,
                         TemporalTensorOwner owner = {}, bool partial = false,
                         std::string value_type = "quantitative") {
        TemporalTensor tensor;
        tensor.key = tensor_key;
        tensor.role = TemporalTensorRole::Field;
        tensor.owner = owner;
        tensor.output = carrier;
        tensor.interpretation = id;
        const auto support = partial ? numeric::StructuralSupport::sparse(values.size(), {0, 1})
                                     : numeric::StructuralSupport::full(values.size());
        std::vector<numeric::SparseRationalEntry> nonzero;
        for (std::size_t i = 0; i < values.size(); ++i)
            if (values[i] != 0)
                nonzero.push_back({i, values[i]});
        tensor.tensor = {numeric::IndexDomain::linear(values.size()), support,
                         numeric::SparseRationalVector(values.size(), std::move(nonzero))};
        endpoint.tensors.push_back(std::move(tensor));
        auto item = Json::object();
        item.set("id", Json::string(id));
        item.set("key", Json::string(std::to_string(tensor_key)));
        item.set("context", Json::string(context));
        item.set("channel", Json::string(channel));
        item.set("quantity",
                 Json::string(channel == "E" ? "synthetic_expression" : "synthetic_" + channel));
        item.set("units", Json::string(channel == "E" ? "synthetic_count" : "dimensionless"));
        item.set("method", Json::string("explicit synthetic fixture"));
        item.set("value_type", Json::string(value_type));
        item.set("evidence", measured_evidence);
        fields.push(item);
    };
    add("baseline.G", 11, "baseline", "G", {1, 0, 1}, 1, {}, false, "indicator");
    add("baseline.E", 12, "baseline", "E", {0, 2, 0}, 1);
    add("baseline.T", 13, "baseline", "T", {0, 0, 1}, 1, {}, false, "indicator");
    add("followup.G", 21, "followup", "G", {1, 1, 1}, 1, {}, false, "indicator");
    add("followup.E", 22, "followup", "E", {1, 3, 2}, 1);
    add("followup.T", 23, "followup", "T", {1, 0, 1}, 1, {}, false, "indicator");
    add("baseline.partial", 31, "baseline", "observation", {0, 2, 0}, 1, {}, true);
    add("baseline.batch", 32, "baseline", "covariate", {0, 1, 1}, 1, {}, false, "categorical");
    add("baseline.face", 33, "baseline", "observation", {4}, 3);
    const auto set_fields = Json::parse(R"(["global.R101","global.R102","global.R103"])");
    for (unsigned i = 0; i < 3; ++i) {
        TemporalTensorOwner owner;
        owner.kind = TemporalOwnerKind::Cell;
        owner.cell = {{1}, 101 + i};
        add("global.R" + std::to_string(101 + i), 41 + i, "global", "observation",
            i == 0 ? Vector{2, 3} : (i == 1 ? Vector{5, 7} : Vector{11, 13}), 10101 + i, owner);
    }
    auto contexts = Json::array();
    for (unsigned i = 0; i < 3; ++i) {
        auto item = Json::object();
        item.set("id", Json::string(i == 0 ? "baseline" : (i == 1 ? "followup" : "global")));
        item.set("block", i < 2 ? Json::string(std::to_string(i + 1)) : Json{});
        item.set(
            "metrics",
            Json::parse(
                R"([{"kind":"identity"},{"kind":"diagonal","values":["2","3","5"]},{"kind":"identity"}])"));
        auto time = Json{};
        if (i < 2) {
            time = Json::object();
            time.set("axis", Json::string("synthetic_visits"));
            time.set("position", Json::string(i == 0 ? "0" : "2"));
            time.set("unit", Json::string("day"));
            time.set("population", Json::string("same synthetic population"));
        }
        item.set("time", time);
        contexts.push(item);
    }
    const auto archive = encode_native_series(NativeTemporalSeries(endpoint, 9876));
    const std::string native(reinterpret_cast<const char *>(archive.data()), archive.size());
    write_text(root / "study.rfcxs", native);
    auto declaration = Json::object();
    declaration.set("schema", Json::string("gadget.scientific.study.v1"));
    declaration.set("id", Json::string("synthetic-longitudinal-get"));
    auto source = Json::object();
    source.set("path", Json::string("study.rfcxs"));
    source.set("sha256", Json::string("sha256:" + sha256(native)));
    source.set("source_namespace", Json::string("9876"));
    declaration.set("native", source);
    auto design = Json::object();
    design.set("question", Json::string("How do independent GET fields and their declared "
                                        "relational context change over two synthetic visits?"));
    design.set("population", Json::string("same synthetic population"));
    design.set("selection",
               Json::string("synthetic authored coordinates; no disease-conditioned selection"));
    auto covariates = Json::array();
    auto batch = Json::object();
    batch.set("name", Json::string("assay_batch"));
    batch.set("description",
              Json::string("declared categorical batch code; no adjustment inferred"));
    batch.set("field", Json::string("baseline.batch"));
    batch.set("evidence", source_evidence);
    covariates.push(batch);
    design.set("covariates", covariates);
    design.set(
        "limitations",
        Json::parse(
            R"(["Synthetic arithmetic fixture; no clinical or biological inference.","Declared batch is visible, not a causal adjustment or significance calibration."])"));
    declaration.set("design", design);
    auto documents = Json::array();
    for (const auto &[id, filename, bytes] :
         {std::tuple{"design", "design.txt", note},
          std::tuple{"observations", "observations.csv", observations}}) {
        auto doc = Json::object();
        doc.set("id", Json::string(id));
        doc.set("path", Json::string(filename));
        doc.set("sha256", Json::string("sha256:" + sha256(bytes)));
        doc.set("format", Json::string("utf8"));
        documents.push(doc);
    }
    declaration.set("documents", documents);
    declaration.set("cells", cells);
    declaration.set("occurrences", occurrences);
    declaration.set("contexts", contexts);
    declaration.set("fields", fields);
    auto sets = Json::array();
    auto set = Json::object();
    set.set("id", Json::string("synthetic.assays"));
    set.set("grade", Json::integer(1));
    set.set("fields", set_fields);
    set.set("carriers", stalks);
    set.set("maps", Json::string("constant"));
    set.set("coordinate_system", Json::string("900"));
    sets.push(set);
    declaration.set("field_sets", sets);
    write_text(root / "study.json", declaration.dump() + "\n");
}
} // namespace gadget::example
