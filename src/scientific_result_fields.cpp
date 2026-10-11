#include "scientific_result_fields.hpp"
#include "scientific_detail.hpp"
#include <rfc/storage/rfcx.hpp>

namespace gadget::scientific_detail {
Json result_fields_view(const NativeRecord &record) {
    std::map<std::pair<const void *, std::uint64_t>, NativeRecord> owners;
    const auto render = [&](const auto &self, const NativeRecord &value) -> Json {
        if (value.kind() == NativeRecordKind::Object) {
            if (value.has("schema") && value.at("schema").kind() == NativeRecordKind::Text) {
                const auto schema = value.at("schema").text();
                if (schema == "gadget.scientific.tensor.view.v1") return result_tensor_view(value);
                if (schema == "gadget.scientific.fields.v1") return self(self, value.at("view"));
                if (schema == "gadget.scientific.publication.view.v1") {
                    const auto bytes = value.at("publication");
                    const auto id = std::pair<const void *, std::uint64_t>{
                        bytes.endpoint().get(), bytes.byte_source().key()};
                    auto owner = owners.find(id);
                    if (owner == owners.end())
                        owner = owners.emplace(id, NativeRecord::decode(bytes.bytes())).first;
                    auto selected = owner->second;
                    for (std::size_t i = 0; i < value.at("path").size(); ++i)
                        selected = selected.at(value.at("path").at(i).text());
                    return self(self, selected);
                }
            }
            auto out = Json::object();
            for (const auto &name : value.names()) out.set(name, self(self, value.at(name)));
            return out;
        }
        if (value.kind() == NativeRecordKind::Array) {
            auto out = Json::array();
            for (std::size_t i = 0; i < value.size(); ++i) out.push(self(self, value.at(i)));
            return out;
        }
        return render_record_view(value);
    };
    return render(render, record.at("view"));
}
Json result_tensor_view(const NativeRecord &record) {
    const auto &tensor = assay_detail::retained_tensor(record, record.at("tensor")).tensor;
    const auto display = record.at("display").text();
    if (display == "rational") {
        const auto position = record.has("index") ? record.at("index").integer() : rfc::numeric::Integer(0);
        if (position < 0 || !mpz_fits_ulong_p(position.get_mpz_t()) ||
            position >= rfc::numeric::Integer(std::to_string(tensor.domain().logical_size())))
            throw std::invalid_argument("scientific scalar coordinate");
        const auto index = position.get_ui();
        if (!record.has("index") && tensor.domain().logical_size() != 1)
            throw std::invalid_argument("scientific scalar support");
        return tensor.support().contains(index) ? Json::string(tensor.state(index)->get_str()) : Json{};
    }
    if (display == "indices") {
        auto out = Json::array();
        for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i) {
            const auto value = tensor.state(i);
            if (!value || value->get_den() != 1 || value->get_num() < 0)
                throw std::invalid_argument("scientific coordinate support");
            out.push(exact_count(value->get_num()));
        }
        return out;
    }
    if (display == "vector") {
        auto out = Json::array();
        if (record.has("coordinates")) {
            const auto &coordinates = assay_detail::retained_tensor(record, record.at("coordinates")).tensor;
            for (std::size_t i = 0; i < coordinates.domain().logical_size(); ++i) {
                const auto value = coordinates.state(i);
                if (!value || value->get_den() != 1 || value->get_num() < 0 ||
                    !mpz_fits_ulong_p(value->get_num().get_mpz_t()))
                    throw std::invalid_argument("scientific source coordinate");
                const auto index = value->get_num().get_ui();
                if (!tensor.support().contains(index))
                    throw std::invalid_argument("scientific selected coordinate outside support");
                out.push(Json::string(tensor.state(index)->get_str()));
            }
        } else {
            for (std::size_t i = 0; i < tensor.domain().logical_size(); ++i)
                out.push(tensor.support().contains(i) ? Json::string(tensor.state(i)->get_str()) : Json{});
        }
        return out;
    }
    if ((display != "sparse_matrix" && display != "declared_map") || tensor.domain().extents().size() != 2)
        throw std::invalid_argument("scientific tensor display");
    auto out = Json::object(), entries = Json::array();
    const auto rows = tensor.domain().extents()[0], columns = tensor.domain().extents()[1];
    out.set("rows", cardinality(rows));
    out.set("columns", cardinality(columns));
    for (const auto &entry : tensor.values().entries()) {
        auto value = Json::object();
        value.set("row", cardinality(entry.index / columns));
        value.set("column", cardinality(entry.index % columns));
        value.set("value", Json::string(entry.value.get_str()));
        entries.push(value);
    }
    out.set("entries", entries);
    if (display == "declared_map")
        out.set("complete_support", Json::boolean(tensor.support().count() == tensor.domain().logical_size()));
    return out;
}
NativeRecord retain_scientific_plan(const NativeRecord &source, const NativeRecord &plan,
    const std::string &source_digest, const std::vector<std::pair<std::string, NativeRecord>> &steps) {
    using N = NativeRecordValue;
    N::Array results;
    std::vector<RawByteSource> sources;
    std::uint64_t key = 1;
    for (const auto &bytes : source.raw_sources())
        key = std::max(key, rfc::machine::checked_size_add(bytes.key(), 1));
    for (const auto &[id, step] : steps) {
        auto archive = step.archive_source(key);
        results.elements.emplace_back(N::Object{{{"id", N(id)},
            {"publication", N(N::ByteSpan{key, 0, archive.bytes().size(), false})}}});
        sources.push_back(std::move(archive));
        key = rfc::machine::checked_size_add(key, 1);
    }
    N value(N::Object{{{"schema", N(std::string("gadget.scientific.plan.fields.v1"))},
        {"plan", import_record_view(render_record_view(plan))}, {"source_digest", N(source_digest)},
        {"steps", N(std::move(results))}}});
    const auto retained = source.derive(std::move(value), std::move(sources));
    return retain_operation(retained, plan, retained.unbound_value());
}
} // namespace gadget::scientific_detail

namespace gadget {
namespace {
void same_result_fields(const NativeRecord &a, const NativeRecord &b);
void same_result_value(const NativeRecord &a, const NativeRecord &b, bool omit_execution = false,
                       bool execution_receipt = false) {
    scientific_detail::require(a.kind() == b.kind(), "replayed scientific value type differs");
    switch (a.kind()) {
    case NativeRecordKind::Null: return;
    case NativeRecordKind::Boolean:
        scientific_detail::require(a.boolean() == b.boolean(), "replayed scientific boolean differs"); return;
    case NativeRecordKind::Integer:
        scientific_detail::require(a.integer() == b.integer(), "replayed scientific integer differs"); return;
    case NativeRecordKind::Rational:
    case NativeRecordKind::Decimal:
    case NativeRecordKind::ExactNumber:
        scientific_detail::require(a.rational() == b.rational(), "replayed scientific scalar differs"); return;
    case NativeRecordKind::Text:
        scientific_detail::require(a.text() == b.text(), "replayed scientific declaration differs"); return;
    case NativeRecordKind::Bytes:
        scientific_detail::require(std::ranges::equal(a.bytes(), b.bytes()), "replayed scientific bytes differ"); return;
    case NativeRecordKind::Tensor:
        scientific_detail::require(rfc::storage::temporal_endpoint_detail::tensor_equal(a.tensor(), b.tensor()),
                                   "replayed scientific tensor differs"); return;
    case NativeRecordKind::Array:
        scientific_detail::require(a.size() == b.size(), "replayed scientific array differs");
        for (std::size_t i = 0; i < a.size(); ++i) same_result_value(a.at(i), b.at(i), omit_execution);
        return;
    case NativeRecordKind::Object:
        scientific_detail::require(a.names() == b.names(), "replayed scientific object differs");
        if (a.has("schema") && a.at("schema").kind() == NativeRecordKind::Text &&
            a.at("schema").text() == "gadget.scientific.publication.view.v1") {
            same_result_value(a.at("schema"), b.at("schema"));
            same_result_value(a.at("path"), b.at("path"));
            same_result_fields(NativeRecord::decode(a.at("publication").bytes()),
                               NativeRecord::decode(b.at("publication").bytes()));
            return;
        }
        for (const auto &name : a.names()) {
            if (omit_execution && (name == "execution" || name == "native_work")) continue;
            if (execution_receipt && (name == "action_limits" || name == "requirements" ||
                name == "retained_requirements" || name == "steps" || name == "execution_steps" ||
                name == "preparation_steps" || name == "green_preparation_steps" ||
                name == "selected_requirements" || name.ends_with("policy") ||
                name == "fallback" || name == "cpu_fallback")) continue;
            same_result_value(a.at(name), b.at(name), omit_execution,
                omit_execution && (name == "native_execution" || name == "native_preparation"));
        }
        return;
    }
}
std::optional<NativeRecord> find_retained_request(const NativeRecord &owner, const NativeRecord &request) {
    if (owner.kind() == request.kind() &&
        (owner.kind() != NativeRecordKind::Object || owner.names() == request.names())) {
        try {same_result_value(owner, request); return owner;}
        catch (const std::invalid_argument &) {}
    }
    if (owner.kind() == NativeRecordKind::Object)
        for (const auto &name : owner.names()) {
            auto found = find_retained_request(owner.at(name), request);
            if (found) return found;
        }
    if (owner.kind() == NativeRecordKind::Array)
        for (std::size_t i = 0; i < owner.size(); ++i) {
            auto found = find_retained_request(owner.at(i), request);
            if (found) return found;
        }
    return std::nullopt;
}
NativeRecord retained_request(const NativeRecord &publication) {
    const auto owner = NativeRecord::decode(publication.at("request_publication").bytes());
    const auto request = find_retained_request(owner, publication.at("request"));
    scientific_detail::require(request.has_value(), "scientific query differs from its retained request owner");
    return *request;
}
void same_result_fields(const NativeRecord &a, const NativeRecord &b) {
    const auto original = a.raw_sources();
    const auto bytes = b.at("source_bytes");
    for (std::size_t i = 0; i < bytes.size(); ++i) {
        const auto source = bytes.at(i);
        const auto found = std::find_if(original.begin(), original.end(),
            [&](const auto &value) {return value.key() == source.byte_source().key();});
        scientific_detail::require(found != original.end() &&
            std::ranges::equal(source.bytes(), found->bytes()), "replayed scientific original bytes differ");
    }
    if (a.has("source_bytes")) same_result_value(a.at("source_bytes"), bytes);
    const auto x = a.endpoint(), y = b.endpoint();
    scientific_detail::require(rfc::storage::encode_rfcx(x->source).bytes ==
        rfc::storage::encode_rfcx(y->source).bytes, "replayed scientific complex differs");
    scientific_detail::require(x->carriers.size() == y->carriers.size() &&
        x->tensors.size() == y->tensors.size(), "replayed scientific carrier/field family differs");
    scientific_detail::require(x->local_blocks.size() == y->local_blocks.size(),
                               "replayed scientific chart family differs");
    for (std::size_t i = 0; i < x->local_blocks.size(); ++i) {
        const auto &p = x->local_blocks[i], &q = y->local_blocks[i];
        scientific_detail::require(p.block_key == q.block_key && p.kind == q.kind &&
            p.cells == q.cells && p.occurrences == q.occurrences && p.payloads.empty() && q.payloads.empty(),
            "replayed scientific chart differs");
    }
    for (std::size_t i = 0; i < x->carriers.size(); ++i) {
        const auto &p = x->carriers[i], &q = y->carriers[i];
        scientific_detail::require(p.key == q.key && p.kind == q.kind && p.grade == q.grade &&
            p.cell == q.cell && p.domain == q.domain && p.coordinate_system == q.coordinate_system,
            "replayed scientific carrier differs");
    }
    for (std::size_t i = 0; i < x->tensors.size(); ++i)
        scientific_detail::require(rfc::storage::temporal_endpoint_detail::tensor_equal(x->tensors[i], y->tensors[i]),
                                   "replayed scientific field or support differs");
    same_result_value(a.at("view"), b.at("view"), true);
    same_result_value(a.at("bindings"), b.at("bindings"));
}
}
NativeRecord replay_scientific_record(const NativeRecord &publication,
                                      std::shared_ptr<CpuResources> resources) {
    scientific_detail::require(publication.at("schema").text() == "gadget.operation.native.v1",
                               "scientific replay requires an owning operation publication");
    const auto result = publication.at("result");
    const auto schema = result.at("schema").text();
    const auto request = retained_request(publication);
    if (schema == "gadget.scientific.plan.fields.v1") {
        scientific_detail::require(result.at("steps").size() > 0, "scientific replay plan is empty");
        std::optional<NativeRecord> original;
        std::optional<NativeRecord> first_step;
        std::vector<std::pair<std::string, NativeRecord>> replayed;
        scientific_detail::require(request.at("steps").size() == result.at("steps").size(),
                                   "scientific replay plan step count differs");
        same_result_value(request, result.at("plan"));
        for (std::size_t i = 0; i < result.at("steps").size(); ++i) {
            const auto declaration = result.at("steps").at(i);
            const auto step = NativeRecord::decode(declaration.at("publication").bytes());
            scientific_detail::require(declaration.at("id").text() == request.at("steps").at(i).at("id").text(),
                                       "scientific replay plan step identity differs");
            same_result_value(step.at("request"), request.at("steps").at(i).at("query"));
            scientific_detail::require(step.at("result").at("schema").text() == "gadget.scientific.result.fields.v1",
                                       "scientific replay plan contains an unmigrated operation");
            const auto bytes = step.at("result").at("source_study_publication").bytes();
            if (!original) {
                original.emplace(NativeRecord::decode(bytes));
                first_step.emplace(step);
            } else {
                scientific_detail::require(std::ranges::equal(bytes,
                    first_step->at("result").at("source_study_publication").bytes()),
                    "scientific replay plan source owners differ");
            }
            replayed.emplace_back(declaration.at("id").text(), replay_scientific_record(step, resources));
        }
        const ScientificStudy study(*original);
        scientific_detail::require(result.at("source_digest").text() == study.digest() &&
            request.at("source_digest").text() == study.digest(), "scientific replay plan source binding differs");
        return scientific_detail::retain_scientific_plan(*original, request, study.digest(), replayed);
    }
    scientific_detail::require(schema == "gadget.scientific.result.fields.v1",
                               "scientific replay operation is not supported");
    const auto original = NativeRecord::decode(result.at("source_study_publication").bytes());
    same_result_value(publication.at("request"), result.at("certificate").at("query"));
    auto study = std::make_shared<const ScientificStudy>(original);
    scientific_detail::require(result.at("certificate").at("source_digest").text() == study->digest(),
                               "scientific replay source binding differs");
    ScientificSession session(std::move(study), std::move(resources));
    const auto replayed = session.execute_record(request);
    same_result_fields(result.at("fields"), replayed.at("result").at("fields"));
    return replayed;
}
} // namespace gadget
