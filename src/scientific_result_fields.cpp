#include "scientific_result_fields.hpp"
#include "scientific_detail.hpp"
#include <rfc/storage/rfcx.hpp>

namespace gadget::scientific_detail {
std::uint64_t ResultFields::source_grade(std::size_t grade, std::optional<std::uint64_t> preferred) {
    if (grade >= endpoint_.source.grade_count())
        throw std::invalid_argument("scientific source grade absent");
    const auto k = rfc::core::Grade{static_cast<std::uint16_t>(grade)};
    const auto domain = rfc::numeric::IndexDomain::linear(endpoint_.source.grade(k).cell_count());
    const rfc::storage::TemporalCarrier *found = nullptr;
    for (const auto &carrier : endpoint_.carriers) {
        if (carrier.kind != rfc::storage::TemporalCarrierKind::SourceGrade || carrier.grade != k ||
            carrier.domain != domain || (preferred && carrier.key != *preferred)) continue;
        if (found) throw std::invalid_argument("ambiguous scientific source grade carrier");
        found = &carrier;
    }
    if (found) return found->key;
    if (preferred) throw std::invalid_argument("scientific input carrier differs from its source grade");
    const auto key = fields_.explicit_carrier(domain.logical_size());
    auto &carrier = endpoint_.carriers.back();
    carrier.kind = rfc::storage::TemporalCarrierKind::SourceGrade;
    carrier.grade = k;
    carrier.coordinate_system = 0;
    return key;
}
NativeRecordValue ResultFields::consume_staged(NativeRecordValue value) {
    if (auto *object = std::get_if<NativeRecordValue::Object>(&value.data)) {
        const auto member = [&](std::string_view name) -> NativeRecordValue * {
            for (auto &[key, item] : object->members) if (key == name) return &item;
            return nullptr;
        };
        const auto *schema = member("schema");
        if (schema && std::get_if<std::string>(&schema->data) &&
            std::get<std::string>(schema->data) == "gadget.scientific.tensor.stage.v1") {
            const auto carrier = [&](std::string_view name) {
                const auto *item = member(name);
                if (!item) throw std::invalid_argument("staged scientific carrier absent");
                const auto &key = std::get<rfc::numeric::Integer>(item->data);
                if (key <= 0 || !mpz_fits_ulong_p(key.get_mpz_t()))
                    throw std::invalid_argument("staged scientific carrier key");
                return static_cast<std::uint64_t>(key.get_ui());
            };
            rfc::storage::TemporalTensor tensor;
            tensor.output = carrier("output");
            tensor.role = rfc::storage::TemporalTensorRole::Field;
            if (member("input")) {
                tensor.input = carrier("input");
                tensor.role = rfc::storage::TemporalTensorRole::Observation;
            }
            if (auto *role = member("role")) {
                const auto &key = std::get<rfc::numeric::Integer>(role->data);
                if (key < 0 || key > static_cast<unsigned>(rfc::storage::TemporalTensorRole::OpaqueExact))
                    throw std::invalid_argument("staged scientific tensor role");
                tensor.role = static_cast<rfc::storage::TemporalTensorRole>(key.get_ui());
            }
            tensor.interpretation = "scientific_result";
            if (auto *interpretation = member("interpretation"))
                tensor.interpretation = std::get<std::string>(interpretation->data);
            if (member("occurrence_key")) {
                const auto &grade = std::get<rfc::numeric::Integer>(member("occurrence_grade")->data);
                const auto &key = std::get<rfc::numeric::Integer>(member("occurrence_key")->data);
                if (grade <= 0 || grade > std::numeric_limits<std::uint16_t>::max() ||
                    key <= 0 || !mpz_fits_ulong_p(key.get_mpz_t()))
                    throw std::invalid_argument("staged scientific occurrence owner");
                tensor.owner.kind = rfc::storage::TemporalOwnerKind::Occurrence;
                tensor.owner.occurrence = {{static_cast<std::uint16_t>(grade.get_ui())}, {key.get_ui()}};
            }
            tensor.tensor = std::move(std::get<rfc::numeric::SparseRationalTensor>(member("tensor")->data));
            std::optional<NativeRecordValue> coordinates;
            if (auto *item = member("coordinates")) {
                if (auto *tensor_coordinates = std::get_if<rfc::numeric::SparseRationalTensor>(&item->data)) {
                    const auto axis = fields_.explicit_carrier(tensor_coordinates->domain().logical_size());
                    coordinates = fields_.on_carrier("selected_coordinates", axis, std::move(*tensor_coordinates));
                } else coordinates = std::move(*item);
            }
            return std::move(reference(fields_.insert(std::move(tensor)),
                std::get<std::string>(member("display")->data), std::move(coordinates))).release();
        }
        for (auto &[key, item] : object->members) {
            (void)key;
            item = consume_staged(std::move(item));
        }
    } else if (auto *array = std::get_if<NativeRecordValue::Array>(&value.data)) {
        for (auto &item : array->elements) item = consume_staged(std::move(item));
    }
    return value;
}
ResultFieldContext::ResultFieldContext(ResultFields &fields, const NativeRecord &source,
    const std::vector<std::vector<std::uint32_t>> &coordinates,
    std::optional<std::pair<std::size_t, std::uint64_t>> field_carrier) : publication(fields) {
    const auto endpoint = source.endpoint();
    for (std::size_t k = 0; k < coordinates.size(); ++k) {
        const auto carrier = publication.source_grade(k,
            field_carrier && field_carrier->first == k ? std::optional(field_carrier->second) : std::nullopt);
        Grade grade{carrier, publication.axis(coordinates[k].size()),
            endpoint->source.grade({static_cast<std::uint16_t>(k)}).cell_count(), {}, {}};
        grade.coordinates.assign(coordinates[k].begin(), coordinates[k].end());
        grade.coordinate_field = publication.coordinates(grade.coordinates, grade.rows);
        rfc::numeric::ExactMatrix lift(grade.size, grade.coordinates.size());
        for (std::size_t i = 0; i < grade.coordinates.size(); ++i) lift.set(grade.coordinates[i], i, 1);
        auto binding = ResultValue::object();
        binding.set("grade", ResultValue::count(k));
        binding.set("source_carrier", ResultValue::count(grade.source));
        binding.set("row_carrier", ResultValue::count(grade.rows));
        binding.set("row_coordinates", ResultValue(grade.coordinate_field));
        binding.set("source_lift", publication.matrix(lift, grade.rows, grade.source));
        publication.bind(std::move(binding));
        grades.push_back(std::move(grade));
    }
}
ResultValue ResultFieldContext::field(const Vector &value, std::size_t k) const {
    const auto &grade = grades.at(k);
    if (value.size() != grade.coordinates.size())
        throw std::invalid_argument("staged scientific source field shape");
    std::vector<rfc::numeric::SparseRationalEntry> entries;
    for (std::size_t i = 0; i < value.size(); ++i)
        if (value[i] != 0) entries.push_back({grade.coordinates[i], value[i]});
    rfc::numeric::SparseRationalTensor tensor(rfc::numeric::IndexDomain::linear(grade.size),
        rfc::numeric::StructuralSupport::sparse(grade.size, grade.coordinates),
        rfc::numeric::SparseRationalVector(grade.size, std::move(entries)));
    return ResultFields::staged(std::move(tensor), {}, grade.source, "vector", grade.coordinate_field);
}
namespace {
NativeRecord result_member(NativeRecord value, const NativeRecord &component) {
    if (value.kind() == NativeRecordKind::Object && value.has("schema") &&
        value.at("schema").kind() == NativeRecordKind::Text) {
        const auto schema = value.at("schema").text();
        if (schema == "gadget.scientific.publication.view.v1") {
            auto selected = NativeRecord::decode(value.at("publication").bytes());
            const auto path = value.at("path");
            require(path.kind() == NativeRecordKind::Array, "scientific archived path type");
            for (std::size_t i = 0; i < path.size(); ++i) selected = result_member(std::move(selected), path.at(i));
            return result_member(std::move(selected), component);
        }
        if (schema == "gadget.scientific.view.overlay.v1") {
            require(component.kind() == NativeRecordKind::Text, "scientific overlay object path type");
            const auto additions = value.at("additions");
            require(additions.kind() == NativeRecordKind::Object, "scientific overlay additions type");
            if (additions.has(component.text())) return additions.at(component.text());
            return result_member(value.at("base"), component);
        }
    }
    if (component.kind() == NativeRecordKind::Text) {
        require(value.kind() == NativeRecordKind::Object, "scientific native object path type");
        return value.at(component.text());
    }
    require(component.kind() == NativeRecordKind::Integer && value.kind() == NativeRecordKind::Array,
            "scientific native array path type");
    const auto index = component.integer();
    require(index >= 0 && mpz_fits_ulong_p(index.get_mpz_t()) && index < value.size(),
            "scientific native array coordinate outside owner");
    return value.at(index.get_ui());
}
NativeRecord result_leaf(NativeRecord value) {
    while (value.kind() == NativeRecordKind::Object && value.has("schema") &&
           value.at("schema").kind() == NativeRecordKind::Text &&
           value.at("schema").text() == "gadget.scientific.publication.view.v1") {
        const auto path = value.at("path");
        auto selected = NativeRecord::decode(value.at("publication").bytes());
        require(path.kind() == NativeRecordKind::Array, "scientific archived path type");
        for (std::size_t i = 0; i < path.size(); ++i) selected = result_member(std::move(selected), path.at(i));
        value = std::move(selected);
    }
    return value;
}
} // namespace
Json result_fields_view(const NativeRecord &record) {
    std::map<std::pair<const void *, std::uint64_t>, NativeRecord> owners;
    const auto render = [&](const auto &self, const NativeRecord &value) -> Json {
        if (value.kind() == NativeRecordKind::Object) {
            if (value.has("schema") && value.at("schema").kind() == NativeRecordKind::Text) {
                const auto schema = value.at("schema").text();
                if (schema == "gadget.scientific.tensor.view.v1") return result_tensor_view(value);
                if (schema == "gadget.scientific.fields.v1") {
                    auto out = self(self, value.at("view"));
                    if (value.has("view_execution")) {
                        require(json_object_is_type(out.raw(), json_type_object),
                                "scientific execution display requires an object");
                        const auto metadata = value.at("view_execution");
                        for (const auto &name : metadata.names()) {
                            require(!out.has(name), "execution display shadows a scientific field");
                            out.set(name, self(self, metadata.at(name)));
                        }
                    }
                    return out;
                }
                if (schema == "gadget.scientific.publication.view.v1") {
                    const auto bytes = value.at("publication");
                    const auto id = std::pair<const void *, std::uint64_t>{
                        bytes.endpoint().get(), bytes.byte_source().key()};
                    auto owner = owners.find(id);
                    if (owner == owners.end())
                        owner = owners.emplace(id, NativeRecord::decode(bytes.bytes())).first;
                    auto selected = owner->second;
                    for (std::size_t i = 0; i < value.at("path").size(); ++i) {
                        const auto component = value.at("path").at(i);
                        if (component.kind() == NativeRecordKind::Text) {
                            if (selected.kind() != NativeRecordKind::Object)
                                throw std::invalid_argument("scientific publication object path type");
                            selected = selected.at(component.text());
                        } else {
                            if (selected.kind() != NativeRecordKind::Array ||
                                component.kind() != NativeRecordKind::Integer)
                                throw std::invalid_argument("scientific publication array path type");
                            const auto index = component.integer();
                            if (index < 0 || !mpz_fits_ulong_p(index.get_mpz_t()) || index >= selected.size())
                                throw std::invalid_argument("scientific publication array path coordinate");
                            selected = selected.at(index.get_ui());
                        }
                    }
                    return self(self, selected);
                }
                if (schema == "gadget.scientific.view.overlay.v1") {
                    auto base = self(self, value.at("base"));
                    const auto additions = self(self, value.at("additions"));
                    if (!json_object_is_type(base.raw(), json_type_object) ||
                        !json_object_is_type(additions.raw(), json_type_object))
                        throw std::invalid_argument("scientific view overlay requires object views");
                    json_object_object_foreach(additions.raw(), name, child) {
                        base.set(name, Json(json_object_get(child)));
                    }
                    return base;
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
    return render(render, record);
}
Json result_tensor_view(const NativeRecord &record) {
    const auto &tensor = assay_detail::retained_tensor(record, record.at("tensor")).tensor;
    const auto display = record.at("display").text();
    if (display == "rational" || display == "rational_object") {
        const auto position = record.has("index") ? record.at("index").integer() : rfc::numeric::Integer(0);
        if (position < 0 || !mpz_fits_ulong_p(position.get_mpz_t()) ||
            position >= rfc::numeric::Integer(std::to_string(tensor.domain().logical_size())))
            throw std::invalid_argument("scientific scalar coordinate");
        const auto index = position.get_ui();
        if (!record.has("index") && tensor.domain().logical_size() != 1)
            throw std::invalid_argument("scientific scalar support");
        if (!tensor.support().contains(index)) return Json{};
        const auto value = Json::string(tensor.state(index)->get_str());
        if (display == "rational") return value;
        auto out = Json::object(); out.set("rational", value); return out;
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
    if (display == "dense_matrix" || display == "matrix_values") {
        if (tensor.domain().extents().size() != 2) throw std::invalid_argument("scientific matrix display");
        const auto shape = tensor.domain().extents();
        auto values = Json::array();
        for (std::size_t i = 0; i < shape[0]; ++i) {
            auto row = Json::array();
            for (std::size_t j = 0; j < shape[1]; ++j) {
                const auto value = tensor.state(i * shape[1] + j);
                row.push(value ? Json::string(value->get_str()) : Json{});
            }
            values.push(row);
        }
        if (display == "matrix_values") return values;
        auto out = Json::object();
        out.set("rows", cardinality(shape[0]));
        out.set("columns", cardinality(shape[1]));
        out.set("values", values);
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
NativeRecord select_scientific_result(const NativeRecord &value, const NativeRecord &path) {
    scientific_detail::require(path.kind() == NativeRecordKind::Array && path.size(),
                               "scientific native path requires components");
    auto selected = value;
    for (std::size_t i = 0; i < path.size(); ++i)
        selected = scientific_detail::result_member(std::move(selected), path.at(i));
    return scientific_detail::result_leaf(std::move(selected));
}

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
