#include "scientific_detail.hpp"
#include <algorithm>
#include <functional>
#include <rfc/machine/checked_size.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
bool boundary(std::string_view bytes, std::size_t at) {
    return at == bytes.size() || (static_cast<unsigned char>(bytes.at(at)) & 0xc0) != 0x80;
}
Json extent_json(const rfc::numeric::IndexDomain &domain) {
    auto result = Json::array();
    for (auto extent : domain.extents())
        result.push(cardinality(extent));
    return result;
}
} // namespace
Json ScientificStudy::State::referenced(const Json &metadata) const {
    auto item = Json::object();
    json_object_object_foreach(metadata.raw(), name, value) {
        (void)value;
        if (std::string(name) != "evidence")
            item.set(name, metadata.at(name));
    }
    if (metadata.has("evidence")) {
        auto references = Json::array();
        const auto evidence = metadata.at("evidence");
        for (std::size_t i = 0; i < evidence.size(); ++i) {
            const auto span = evidence.at(i);
            auto reference = Json::object();
            for (const auto *name : {"document", "start", "end"})
                reference.set(name, span.at(name));
            const auto id = span.at("document").text();
            const auto &document = documents.at(id);
            reference.set("document_sha256", document.metadata.at("sha256"));
            const auto start = count(span.at("start"), document.bytes.size(), "evidence start");
            const auto end = count(span.at("end"), document.bytes.size(), "evidence end");
            const auto key = std::make_tuple(id, start, end);
            std::lock_guard lock(evidence_reference_mutex);
            if (!evidence_reference_hashes.contains(key)) {
                evidence_reference_hashes.emplace(
                    key, start == 0 && end == document.bytes.size()
                             ? document.metadata.at("sha256").text()
                             : "sha256:" + sha256(document.bytes.substr(start, end - start)));
            }
            reference.set("text_sha256", Json::string(evidence_reference_hashes.at(key)));
            references.push(reference);
        }
        item.set("evidence", references);
    }
    return Json::parse(item.dump());
}
Json ScientificStudy::source_index() const {
    const auto &s = *state_;
    const auto endpoint = s.directory->pin_snapshot();
    auto result = Json::object(), counts = Json::object();
    result.set("schema", Json::string("gadget.scientific.source.index.v1"));
    result.set("study_id", s.declaration.at("id"));
    result.set("source_digest", Json::string(s.digest));
    result.set("source_sequence", Json::string(std::to_string(endpoint->sequence)));
    result.set("declaration_sha256", Json::string(s.declaration_digest));
    result.set("native_sha256", Json::string(s.native_digest));
    result.set("design", s.declaration.at("design"));
    counts.set("fields", cardinality(s.fields.size()));
    counts.set("cells", cardinality(s.cell_count()));
    counts.set("occurrences", cardinality(s.occurrence_count()));
    counts.set("contexts", cardinality(s.contexts.size()));
    counts.set("operators", cardinality(s.operators.size()));
    counts.set("carriers", cardinality(endpoint->carriers.size()));
    counts.set("documents", cardinality(s.documents.size()));
    counts.set("field_sets", cardinality(s.field_sets.size()));
    counts.set("publications", cardinality(s.publications.size()));
    result.set("catalog_sections", counts);
    result.set("read_schema", Json::string("gadget.scientific.source.read.v1"));
    result.set("evidence_profile", Json::string("document-bound-byte-spans-v1"));
    return seal(result, "index_sha256");
}
Json ScientificStudy::read_source(const Json &request) const {
    return read_source(request, source_execution_lease("scientific-source-read-v1"));
}
Json ScientificStudy::read_source(const Json &input, rfc::runtime::ExecutionLease lease) const {
    lease.check_cancelled();
    rfc::runtime::admit_execution({1, 0, 0, 0, 0}, lease.policy());
    const auto request = Json::parse(input.dump());
    members(request, {"schema", "source_digest", "publication", "source_sequence", "stream",
                      "offset", "limit"});
    require(request.at("schema").text() == "gadget.scientific.source.read.v1" &&
                request.at("source_digest").text() == digest(),
            "source read schema/owner mismatch");
    const auto publication = text(request.at("publication"));
    const ScientificStudy *selected = this;
    if (publication != "self") {
        require(state_->publications.contains(publication), "source read publication unknown");
        selected = state_->publications.at(publication).get();
    }
    const auto &s = *selected->state_;
    const auto endpoint = s.directory->pin_snapshot();
    require(request.at("source_sequence").text() == std::to_string(endpoint->sequence),
            "source read cutoff mismatch");
    const auto stream = request.at("stream");
    const auto kind = text(stream.at("kind"));
    const auto limit =
        count(request.at("limit"), std::numeric_limits<std::size_t>::max(), "read limit");
    require(limit > 0, "source read needs a positive caller page size");
    std::size_t total{};
    std::function<Json(std::size_t)> record;
    auto descriptor = Json::object();
    std::string bytes;
    const Document *document{};
    if (kind == "catalog") {
        members(stream, {"kind", "section"});
        const auto section = text(stream.at("section"));
        // Immutable key indices avoid rescanning earlier entries for a late page.
        auto bind = [&](const auto &items, const auto &keys, auto transform) {
            total = items.size();
            record = [&items, &keys, transform](std::size_t at) {
                return transform(*items.find(keys.at(at)));
            };
        };
        if (section == "fields")
            bind(s.fields, s.catalog_keys.at(section), [&](const auto &item) {
                auto result = s.field_metadata(item.second);
                result.set("definition", s.referenced(item.second.declaration));
                return result;
            });
        else if (section == "cells") {
            total=s.cell_count();
            record=[&s](std::size_t at){return s.referenced(s.cell_record(s.cell_key_at(at)));};
        }
        else if (section == "occurrences") {
            total=s.occurrence_count();
            record=[&s,endpoint](std::size_t at) {
                const auto key=s.occurrence_key_at(at);
                auto result = s.referenced(s.occurrence_record(key));
                const auto pin = s.directory->occurrence(key);
                auto identity = [](const auto &cell) {
                    auto value = Json::object();
                    value.set("grade", Json::integer(cell.grade.value));
                    value.set("key", Json::string(std::to_string(cell.key)));
                    return value;
                };
                result.set("parent", identity(pin.parent()));
                result.set("participant", identity(pin.participant()));
                result.set("coefficient", Json::string(pin.coefficient().get_str()));
                auto aspects = Json::array();
                const auto &schema = endpoint->source.signature().aspects();
                for (std::size_t a = 0; a < schema.size(); ++a) {
                    const auto &definition = schema.at(static_cast<rfc::composite::AspectId>(a));
                    auto aspect = Json::object();
                    aspect.set("id", cardinality(a));
                    aspect.set("name", Json::string(definition.name));
                    aspect.set("role", Json::integer(static_cast<std::int64_t>(definition.role)));
                    aspect.set("value", Json::boolean(pin.aspect(definition.id)));
                    aspects.push(aspect);
                }
                result.set("aspects", aspects);
                return result;
            };
        } else if (section == "contexts")
            bind(s.contexts, s.catalog_keys.at(section),
                 [&](const auto &item) { return s.referenced(item.second.declaration); });
        else if (section == "operators")
            bind(s.operators, s.catalog_keys.at(section),
                 [&](const auto &item) { return s.referenced(s.operator_metadata(item.first)); });
        else if (section == "field_sets")
            bind(s.field_sets, s.catalog_keys.at(section),
                 [&](const auto &item) { return s.referenced(item.second.declaration); });
        else if (section == "carriers") {
            total = endpoint->carriers.size();
            record = [&, endpoint](std::size_t at) {
                return s.carrier_metadata(endpoint->carriers.at(at));
            };
        } else if (section == "documents")
            bind(s.documents, s.catalog_keys.at(section), [&](const auto &item) {
                auto result = Json::parse(item.second.metadata.dump());
                result.set("bytes", cardinality(item.second.bytes.size()));
                return result;
            });
        else if (section == "publications")
            bind(s.publications, s.catalog_keys.at(section), [&](const auto &item) {
                auto result = s.referenced(s.publication_declarations.at(item.first));
                result.set("source_sequence",
                           Json::string(std::to_string(item.second->endpoint()->sequence)));
                return result;
            });
        else
            throw std::invalid_argument("unknown source catalogue section");
        descriptor.set("ordering",Json::string(section=="carriers"?"native_declaration_order":
            s.compact && (section=="cells" || section=="occurrences")?"native_attaching_order":"stable_key_order"));
    } else if (kind == "document") {
        members(stream, {"kind", "id"});
        const auto id = text(stream.at("id"));
        require(s.documents.contains(id), "source read document unknown");
        document = &s.documents.at(id);
        require(document->metadata.at("format").text() == "utf8",
                "text read requires UTF-8 document");
        total = document->bytes.size();
        descriptor = Json::parse(document->metadata.dump());
    } else if (kind == "field" || kind == "operator") {
        members(stream, {"kind", "id"});
        const auto id = text(stream.at("id"));
        require(kind == "field" ? s.fields.contains(id) : s.operators.contains(id),
                "source read tensor unknown");
        const auto pin =
            kind == "field" ? s.field_pins.at(s.fields.at(id).tensor) : s.operator_pins.at(id);
        const auto tensor = pin.tensor();
        total = tensor->domain().logical_size();
        descriptor = kind == "field" ? s.field_metadata(s.fields.at(id)) : s.operator_metadata(id);
        descriptor = s.referenced(descriptor);
        descriptor.set("definition", s.referenced(kind == "field" ? s.fields.at(id).declaration
                                                                  : s.operators.at(id)));
        descriptor.set("source_namespace", Json::string(std::to_string(pin.source_namespace())));
        record = [&, tensor, pin](std::size_t index) {
            lease.check_cancelled();
            auto item = Json::object();
            item.set("source_coordinate", cardinality(index));
            auto coordinates = Json::array();
            for (auto coordinate : tensor->domain().unflatten(index))
                coordinates.push(cardinality(coordinate));
            item.set("coordinates", coordinates);
            const auto value = tensor->state(index);
            item.set("reported", Json::boolean(value.has_value()));
            item.set("value", value ? Json::string(value->get_str()) : Json{});
            if (value) {
                const auto bits = std::max(mpz_sizeinbase(value->get_num().get_mpz_t(), 2),
                                           mpz_sizeinbase(value->get_den().get_mpz_t(), 2));
                rfc::runtime::admit_execution({0, 1, 0, 0, bits}, lease.policy());
            }
            const auto carrier = pin.carrier();
            if (kind == "field" &&
                carrier->kind == rfc::storage::TemporalCarrierKind::SourceGrade) {
                const auto cell = s.coordinate(carrier->grade.value,index);
                auto identity = Json::object();
                identity.set("grade", Json::integer(cell.grade.value));
                identity.set("key", Json::string(std::to_string(cell.key)));
                item.set("source_cell", identity);
            }
            return item;
        };
        descriptor.set("extents", extent_json(tensor->domain()));
    } else
        throw std::invalid_argument("unknown source stream kind");
    const auto offset = count(request.at("offset"), total, "source read offset");
    auto end = offset + std::min(limit, total - offset);
    if (document) {
        require(boundary(document->bytes, offset), "text read starts inside UTF-8 code point");
        while (end > offset && !boundary(document->bytes, end))
            --end;
        require(end > offset || offset == total, "page size cannot hold the next UTF-8 code point");
    }
    const auto returned = end - offset;
    rfc::runtime::admit_execution({rfc::machine::checked_size_add(1, document ? 1 : returned),
                                   (kind == "field" || kind == "operator") && returned ? 1u : 0u, 0,
                                   returned, 0},
                                  lease.policy());
    auto result = Json::object(), binding = Json::object(), records = Json::array();
    binding.set("source_digest", Json::string(digest()));
    binding.set("publication", request.at("publication"));
    binding.set("publication_digest", Json::string(selected->digest()));
    binding.set("source_sequence", request.at("source_sequence"));
    binding.set("stream", stream);
    result.set("schema", Json::string("gadget.scientific.source.page.v1"));
    result.set("source_digest", Json::string(digest()));
    result.set("request", request);
    result.set("binding", binding);
    auto execution = Json::object();
    execution.set("policy", Json::string(lease.policy().name));
    execution.set("backend", Json::string("exact_cpu"));
    execution.set("fallback", Json::boolean(lease.policy().backend !=
                                            rfc::runtime::ExecutionBackend::ExactCpu));
    result.set("execution", execution);
    result.set("stream_sha256", Json::string("sha256:" + sha256(binding.dump())));
    result.set("descriptor", descriptor);
    result.set("offset", cardinality(offset));
    result.set("end", cardinality(end));
    result.set("total", cardinality(total));
    result.set("complete", Json::boolean(end == total));
    if (document) {
        bytes = document->bytes.substr(offset, returned);
        result.set("text", Json::string(bytes));
        result.set("text_sha256", Json::string("sha256:" + sha256(bytes)));
    } else {
        for (std::size_t i = offset; i < end; ++i) {
            lease.check_cancelled();
            records.push(record(i));
        }
        result.set("records", records);
    }
    if (end < total) {
        auto next = Json::parse(request.dump());
        next.set("offset", cardinality(end));
        result.set("next_request", next);
    } else
        result.set("next_request", Json{});
    result.set("evidence_profile", Json::string("document-bound-byte-spans-v1"));
    lease.check_cancelled();
    return seal(result, "page_sha256");
}
} // namespace gadget
