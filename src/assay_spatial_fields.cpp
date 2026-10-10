#include "genomic_detail.hpp"
#include <gadget/record_view.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/machine/checked_size.hpp>
#include <rfc/local/carrier.hpp>

namespace gadget::assay_detail {
namespace {
using namespace rfc;
using N = NativeRecordValue;
using Entry = numeric::SparseRationalEntry;
std::uint64_t successor(std::uint64_t key) {
    require(key < std::numeric_limits<std::uint64_t>::max(), "native identity space exhausted");
    return key + 1;
}
void retain_source(core::RelationalComplexBuilder &builder, const core::RelationalComplex &source) {
    const auto &signature = source.signature();
    for (auto identity : source.grade({0}).identities) builder.declare_witness_identity(identity);
    for (std::size_t k = 1; k < source.grade_count(); ++k) {
        const core::Grade grade{static_cast<std::uint16_t>(k)};
        builder.declare_grade(grade);
        const auto &arena = source.grade(grade);
        for (std::size_t i = 0; i < arena.cell_count(); ++i) {
            std::vector<core::BoundaryOccurrenceInput> first;
            std::vector<core::HigherIncidenceInput> higher;
            const auto [begin, end] = arena.incidence_range(i);
            for (auto p = begin; p < end; ++p) {
                const auto exists = arena.aspects.get(signature.existence_aspect(), p);
                const auto negative = arena.aspects.get(signature.orientation_aspect(), p);
                const Rational coefficient = source.incidence_coefficient(grade, p);
                const Rational magnitude = coefficient < 0 ? -coefficient : coefficient;
                std::vector<core::AspectValue> aspects;
                for (std::size_t a = 0; a < signature.aspects().size(); ++a)
                    if (a != signature.existence_aspect() && a != signature.orientation_aspect())
                        aspects.push_back({static_cast<composite::AspectId>(a),
                                           arena.aspects.get(static_cast<composite::AspectId>(a), p)});
                const auto key = arena.occurrence_identities[p];
                if (k == 1)
                    first.emplace_back(source.grade({0}).identities[arena.lower_index[p]], exists,
                                       negative, magnitude, aspects, key);
                else
                    higher.emplace_back(arena.lower_index[p], exists, negative, magnitude, aspects, key);
            }
            if (k == 1) builder.add_grade_one(arena.relation_identities[i], arena.relation_types[i], first);
            else builder.add_higher(grade, arena.relation_identities[i], arena.relation_types[i], higher);
        }
    }
}
numeric::SparseRationalTensor extend(const numeric::SparseRationalTensor &input, std::size_t size) {
    require(input.domain().rank() == 1 && size >= input.domain().logical_size(),
            "spatial extension requires a source-grade scalar field");
    std::vector<std::size_t> support;
    input.support().for_each_active([&](auto i) { support.push_back(i); });
    return {numeric::IndexDomain::linear(size), support.empty()
                ? numeric::StructuralSupport::sparse(size, {})
                : numeric::StructuralSupport::packed(size, support),
            numeric::SparseRationalVector(size, input.values().entries())};
}
} // namespace
AlignmentFieldsSource spatial_fields_source(const ScientificStudy &source, const Json &assay,
                                            const Json &request) {
    const auto original = source.endpoint();
    require(assay.at("declaration").at("carrier").text() == "observations",
            "native spatial fields require an observation-carrier assay");
    require(original->source.grade({0}).cell_count() == assay.at("observations").size(),
            "spatial coordinates must cover the complete observation carrier");
    auto signature = original->source.signature();
    const auto name = request.at("relation_type").text();
    core::RelationTypeId type{};
    bool known = false;
    for (std::size_t t = 0; t < signature.type_count(); ++t)
        if (signature.type(static_cast<core::RelationTypeId>(t)).name == name) {
            type = static_cast<core::RelationTypeId>(t); known = true; break;
        }
    if (!known) type = signature.add_type(name);
    core::RelationalComplexBuilder builder(signature, original->source.semantic_scale());
    retain_source(builder, original->source);
    builder.declare_grade({1});
    std::uint64_t relation = 1, occurrence = 1;
    std::size_t begin = 0;
    for (std::size_t k = 1; k < original->source.grade_count(); ++k) {
        const auto &arena = original->source.grade({static_cast<std::uint16_t>(k)});
        if (k == 1) {
            begin = arena.cell_count();
            for (const auto key : arena.relation_identities) relation = std::max(relation, successor(key.value));
        }
        for (const auto key : arena.occurrence_identities) occurrence = std::max(occurrence, successor(key.value));
    }
    std::vector<Entry> distances;
    const auto receipt = spatial_analysis(assay, request, [&](auto i, auto j, const Rational &q) {
        const auto &identities = original->source.grade({0}).identities;
        builder.add_grade_one(core::RelationKey{relation}, type,
            {{identities[i], true, true, Rational{1}, {}, core::OccurrenceKey{occurrence}},
             {identities[j], true, false, Rational{1}, {}, core::OccurrenceKey{successor(occurrence)}}});
        relation = successor(relation);
        occurrence = successor(successor(occurrence));
        distances.push_back({machine::checked_size_add(begin, distances.size()), q});
    });
    storage::TemporalEndpoint endpoint = *original;
    endpoint.source = builder.freeze();
    endpoint.sequence = 0;
    const auto edges = endpoint.source.grade({1}).cell_count();
    std::set<std::uint64_t> expanded;
    std::uint64_t carrier_key = 1, tensor_key = 1;
    std::optional<std::uint64_t> points;
    for (auto &carrier : endpoint.carriers) {
        carrier_key = std::max(carrier_key, successor(carrier.key));
        if (carrier.kind != storage::TemporalCarrierKind::SourceGrade) continue;
        if (carrier.grade.value == 0) points = carrier.key;
        if (carrier.grade.value == 1 && carrier.domain.logical_size() != edges) {
            require(carrier.domain.rank() == 1, "spatial source carrier must be linear");
            carrier.domain = numeric::IndexDomain::linear(edges);
            expanded.insert(carrier.key);
        }
    }
    for (auto &tensor : endpoint.tensors) {
        tensor_key = std::max(tensor_key, successor(tensor.key));
        if ((tensor.output && expanded.contains(*tensor.output)) ||
            (tensor.input && expanded.contains(*tensor.input))) {
            require(tensor.role == storage::TemporalTensorRole::Field && !tensor.input,
                    "spatial topology extension needs an explicit transport for existing operators");
            tensor.tensor = extend(tensor.tensor, edges);
        }
    }
    require(points.has_value(), "native observation carrier missing");
    const auto edge_carrier = carrier_key;
    endpoint.carriers.push_back({carrier_key, storage::TemporalCarrierKind::SourceGrade, {1}, {}, 0,
                                numeric::IndexDomain::linear(edges)});
    N::Array fields;
    const auto add = [&](const std::string &id, std::uint64_t carrier, std::size_t size,
                         numeric::StructuralSupport support, std::vector<Entry> entries,
                         const std::string &quantity, const std::string &units) {
        storage::TemporalTensor tensor;
        tensor.key = tensor_key; tensor_key = successor(tensor_key);
        tensor.role = storage::TemporalTensorRole::Field; tensor.output = carrier;
        tensor.interpretation = "spatial." + id;
        tensor.tensor = {numeric::IndexDomain::linear(size), std::move(support),
                         numeric::SparseRationalVector(size, std::move(entries))};
        fields.elements.emplace_back(N::Object{{{"id", N("spatial." + id)},
            {"key", N(std::to_string(tensor.key))}, {"quantity", N(quantity)}, {"units", N(units)}}});
        endpoint.tensors.push_back(std::move(tensor));
    };
    const auto coordinates = request.at("coordinates"), observations = assay.at("observations");
    const auto documents = source.native_publication().at("documents");
    std::optional<Axis> axis;
    for (std::size_t i = 0; i < documents.size(); ++i) {
        if (documents.at(i).at("id").text() != "carrier_axis") continue;
        const auto bytes = documents.at(i).at("bytes").bytes();
        axis = read_axis(std::string(reinterpret_cast<const char *>(bytes.data()), bytes.size()));
    }
    require(axis.has_value() && axis->ids.size() == observations.size(),
            "spatial source requires its original observation axis bytes");
    std::vector<Entry> starts, ends;
    for (std::size_t i = 0; i < observations.size(); ++i) {
        require(axis->ids[i] == observations.at(i).at("id").text(), "spatial observation identity/order");
        starts.push_back({i, Rational(numeric::Integer(std::to_string(axis->spans[i].first)))});
        ends.push_back({i, Rational(numeric::Integer(std::to_string(axis->spans[i].second)))});
    }
    add("source_start", *points, observations.size(), numeric::StructuralSupport::full(observations.size()),
        std::move(starts), "source_byte_offset", "byte");
    add("source_end", *points, observations.size(), numeric::StructuralSupport::full(observations.size()),
        std::move(ends), "source_byte_offset", "byte");
    add("source_document", *points, observations.size(), numeric::StructuralSupport::full(observations.size()),
        {}, "source_document_index", "index");
    for (std::size_t d = 0; d < coordinates.size(); ++d) {
        std::vector<Entry> values;
        const auto column = coordinates.at(d).text();
        for (std::size_t i = 0; i < observations.size(); ++i) {
            const Rational q = exact_rational(observations.at(i).at(column).text());
            if (q != 0) values.push_back({i, q});
        }
        add("coordinate:" + column, *points, observations.size(), numeric::StructuralSupport::full(observations.size()),
            std::move(values), "physical_coordinate", request.at("units").text());
    }
    std::vector<std::size_t> support;
    for (const auto &entry : distances) support.push_back(entry.index);
    add("squared_distance", edge_carrier, edges, begin == 0
            ? numeric::StructuralSupport::full(edges) : numeric::StructuralSupport::packed(edges, support),
        std::move(distances), "physical_squared_distance", request.at("units").text() + "^2");
    N::Array charts;
    if (request.has("charts")) {
        require(request.at("charts").text() == "none" || request.at("charts").text() == "grade_one_stars",
                "spatial charts must be none or grade_one_stars");
        if (request.at("charts").text() == "grade_one_stars") {
            std::uint64_t block_key = 1;
            for (const auto &block : endpoint.local_blocks) block_key = std::max(block_key, successor(block.block_key));
            local::RelationalLocalCarrier carrier(endpoint.source);
            const auto &arena = endpoint.source.grade({1});
            const auto &identities = endpoint.source.grade({0}).identities;
            for (std::size_t i = 0; i < identities.size(); ++i) {
                std::vector<core::CellId> seeds{{{0}, static_cast<std::uint32_t>(i)}};
                const auto [first, last] = arena.coface_incidence_range(i);
                for (auto p = first; p < last; ++p)
                    seeds.push_back({{1}, arena.upper_index[arena.incidence_by_lower[p]]});
                const auto patch = carrier.closure(std::move(seeds));
                endpoint.local_blocks.push_back(storage::make_strict_local_block(endpoint.source, block_key, patch));
                charts.elements.emplace_back(N::Object{{{"id", N("star:" + std::to_string(identities[i]))},
                    {"block", N(std::to_string(block_key))}, {"anchor", N(std::to_string(identities[i]))}}});
                block_key = successor(block_key);
            }
        }
    }
    N declaration(N::Object{{{"schema", N(std::string("gadget.spatial.fields.v1"))},
        {"assay", N(request.at("assay").text())}, {"frame", N(request.at("frame").text())},
        {"units", N(request.at("units").text())}, {"neighbor_begin", N(numeric::Integer(std::to_string(begin)))},
        {"neighbor_count", N(numeric::Integer(std::to_string(edges - begin)))},
        {"fields", N(std::move(fields))}, {"charts", N(std::move(charts))}, {"request", import_record_view(request)},
        {"execution", import_record_view(receipt.has("execution") ? receipt.at("execution") : Json::object())}}});
    return {std::move(endpoint), std::move(declaration), {}};
}
} // namespace gadget::assay_detail
