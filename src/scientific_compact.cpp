#include "scientific_detail.hpp"
#include <algorithm>
#include <charconv>
#include <rfc/machine/checked_size.hpp>
#include <rfc/storage/native_series_codec.hpp>

namespace gadget {
namespace {
using namespace scientific_detail;
using namespace rfc::storage;
using N = NativeRecordValue;
std::size_t address(const Rational &value) {
    require(value >= 0 && value.get_den() == 1 &&
                value <=
                    rfc::numeric::Integer(std::to_string(std::numeric_limits<std::size_t>::max())),
            "native evidence coordinate is outside the address domain");
    return std::stoull(value.get_num().get_str());
}
const TemporalTensor &tensor(const TemporalEndpoint &endpoint, std::uint64_t key) {
    const auto found = std::find_if(endpoint.tensors.begin(), endpoint.tensors.end(),
                                    [&](const auto &t) { return t.key == key; });
    require(found != endpoint.tensors.end(), "native evidence field is missing");
    return *found;
}
const rfc::numeric::SparseRationalTensor &column(const TemporalEndpoint &endpoint,
                                                 std::uint64_t key, std::size_t grade) {
    const auto &field = tensor(endpoint, key);
    require(field.role == TemporalTensorRole::Field && field.output && !field.input,
            "evidence column must name a native source-grade field");
    const auto carriers = temporal_endpoint_detail::carriers(endpoint);
    require(carriers.contains(*field.output), "native evidence carrier is missing");
    const auto &carrier = *carriers.at(*field.output);
    const auto size = endpoint.source.grade({static_cast<std::uint16_t>(grade)}).cell_count();
    require(carrier.kind == TemporalCarrierKind::SourceGrade && carrier.grade.value == grade &&
                carrier.domain == field.tensor.domain() && field.tensor.domain().rank() == 1 &&
                field.tensor.domain().logical_size() == size &&
                field.tensor.support().count() == size,
            "native evidence column must cover the complete applicable grade");
    return field.tensor;
}
Json reference(const std::string &document, std::size_t start, std::size_t end) {
    auto out = Json::object();
    out.set("document", Json::string(document));
    out.set("start", cardinality(start));
    out.set("end", cardinality(end));
    return out;
}
} // namespace
NativeRecord ScientificStudy::publish_annotations(const NativeRecord &producer,
                                                  NativeRecordValue study,
                                                  NativeRecordValue documents) {
    require(producer.publication_history().has_value(), "producer has no immutable native history");
    auto history = *producer.publication_history();
    auto &members = std::get<N::Object>(study.data).members;
    require(std::none_of(members.begin(), members.end(),
                         [](const auto &m) { return m.first == "native"; }),
            "producer binding is supplied by the actual history");
    const auto hash =
        "sha256:" + sha256_chunks([&](const auto &sink) { stream_native_series(history, sink); });
    members.emplace_back(
        "native", N(N::Object{{{"path", N(std::string("publication-prefix"))},
                               {"sha256", N(hash)},
                               {"source_namespace", N(std::to_string(history.source_namespace()))},
                               {"cutoff", N(std::to_string(history.last_sequence()))}}}));
    N root(N::Object{{{"schema", N(std::string("gadget.biomedical.native.v3"))},
                      {"history", N(N::Object{{{"kind", N(std::string("publication_prefix"))},
                                               {"cutoff", N(rfc::numeric::Integer(std::to_string(
                                                              history.last_sequence())))}}})},
                      {"study", std::move(study)},
                      {"documents", std::move(documents)}}});
    return NativeRecord::annotate(root, producer);
}
NativeRecord ScientificStudy::producer_publication() const {
    require(state_->compact, "study does not declare a generating application publication");
    return state_->native_publication->historical_publication(state_->history->cutoff());
}
void ScientificStudy::State::admit_compact(NativeRecord publication) {
    compact = true;
    native_publication = std::move(publication);
    const auto &owner = *native_publication;
    const auto declaration = owner.at("study");
    this->declaration = declaration;
    const auto view = render_record_view(declaration);
    members(view,
            {"schema", "id", "native", "design", "documents", "cell_bindings",
             "occurrence_evidence", "contexts", "fields", "field_sets"},
            {"operators"});
    require(view.at("schema").text() == "gadget.scientific.study.v3",
            "unsupported compact study schema");
    (void)text(view.at("id"));
    const auto binding = owner.at("history");
    require(binding.at("kind").text() == "publication_prefix",
            "compact study needs an actual prefix");
    const auto retained = owner.publication_history();
    require(retained.has_value(), "missing native history");
    const auto cutoff = std::stoull(binding.at("cutoff").integer().get_str());
    history =
        std::make_shared<const NativeSeriesPrefix>(NativeSeriesPrefix::capture(*retained, cutoff));
    require(cutoff < retained->last_sequence(), "study annotation must follow its producer");
    const auto tail = history->tail();
    const auto endpoint = tail.pin();
    for (auto sequence = cutoff; sequence < retained->last_sequence(); ++sequence) {
        const auto update =
            decode_temporal_segment(retained->segment_bytes(static_cast<std::size_t>(
                                        sequence - retained->first_sequence())),
                                    retained->limits().temporal);
        require(update.source.grades.empty() &&
                    update.source.target_grade_count == endpoint->source.grade_count() &&
                    update.removed_blocks.empty() && update.blocks.empty() &&
                    update.removed_carriers.empty() && update.carriers.empty() &&
                    update.removed_tensors.empty() && update.tensors.empty() &&
                    update.correspondences.empty() && update.removed_finite_systems.empty() &&
                    update.finite_systems.empty() && !update.codes && !update.source_migration,
                "study annotation changed its generating mathematical source");
    }
    native_digest = "sha256:" + sha256_chunks([&](const auto &sink) {
                        stream_native_series(history->history(), sink);
                    });
    const auto native = view.at("native");
    members(native, {"path", "sha256", "source_namespace", "cutoff"});
    require(native.at("sha256").text() == native_digest &&
                key(native.at("source_namespace")) == history->source_namespace() &&
                native.at("cutoff").text() == std::to_string(cutoff),
            "compact study producer binding mismatch");
    declaration_digest = native_declaration_digest(declaration);
    digest = "sha256:" +
             sha256("gadget.scientific.study.v3\n" + declaration_digest + "\n" + native_digest);
    directory = std::make_shared<const NativeAccessionDirectory>(
        tail, std::vector<RelationalCodeObservationQuery>{},
        [](auto value) { return rfcx_detail::fnv1a64(value); }, native_accession_limits());
    const auto design = view.at("design");
    members(design, {"question", "population", "selection", "covariates", "limitations"});
    for (const auto *name : {"question", "population", "selection"})
        (void)text(design.at(name));
    source_array(design.at("limitations"), "limitations");
    for (std::size_t i = 0; i < design.at("limitations").size(); ++i)
        (void)text(design.at("limitations").at(i));
    const auto docs = view.at("documents");
    const auto raw = owner.at("documents");
    const auto producer = owner.historical_publication(history->cutoff());
    std::map<std::uint64_t, RawByteSource> original_sources;
    for (auto source : producer.raw_sources())
        original_sources.emplace(source.key(), std::move(source));
    source_array(docs, "documents");
    require(raw.size() == docs.size(), "native document count mismatch");
    for (std::size_t i = 0; i < docs.size(); ++i) {
        const auto item = docs.at(i);
        members(item, {"id", "path", "sha256", "format"});
        const auto id = text(item.at("id"));
        require(raw.at(i).at("id").text() == id, "document identity mismatch");
        const auto retained_source = raw.at(i).at("bytes").byte_source();
        require(original_sources.contains(retained_source.key()),
                "document is not retained by the generating source");
        const auto original_bytes = original_sources.at(retained_source.key()).bytes();
        const auto bytes = raw.at(i).at("bytes").bytes();
        require(original_bytes.size() == bytes.size() &&
                    std::equal(original_bytes.begin(), original_bytes.end(), bytes.begin()),
                "document differs from the original generating bytes");
        const auto span =
            std::string_view(reinterpret_cast<const char *>(bytes.data()), bytes.size());
        require(item.at("sha256").text() == "sha256:" + sha256(span),
                "native document byte digest mismatch");
        require(item.at("format").text() == "utf8" || item.at("format").text() == "binary",
                "document format");
        require(documents.emplace(id, Document{declaration.at("documents").at(i), span}).second,
                "duplicate document");
    }
    source_array(design.at("covariates"), "covariates");
    std::set<std::string> covariates;
    for (std::size_t i = 0; i < design.at("covariates").size(); ++i) {
        const auto item = design.at("covariates").at(i);
        members(item, {"name", "description", "field", "evidence"});
        require(covariates.insert(text(item.at("name"))).second, "duplicate covariate");
        (void)text(item.at("description"));
        (void)evidence(item.at("evidence"));
    }
    require(view.at("occurrence_evidence").text() == "parent_cell",
            "unsupported occurrence evidence binding");
    const auto bindings = view.at("cell_bindings");
    source_array(bindings, "cell bindings");
    require(bindings.size() == endpoint->source.grade_count(),
            "complete cell-grade bindings required");
    std::set<std::string> prefixes;
    for (std::size_t k = 0; k < bindings.size(); ++k) {
        const auto item = bindings.at(k);
        members(item, {"grade", "id_prefix", "label_prefix", "kind", "evidence"});
        require(count(item.at("grade"), bindings.size() - 1, "binding grade") == k,
                "cell bindings must follow native grade order");
        const auto prefix = text(item.at("id_prefix"));
        for (const auto &existing : prefixes)
            require(!prefix.starts_with(existing) && !existing.starts_with(prefix),
                    "overlapping native cell prefixes");
        require(prefix.back() == ':' && prefixes.insert(prefix).second,
                "ambiguous native cell prefix");
        (void)text(item.at("label_prefix"));
        (void)text(item.at("kind"));
        cell_bindings.push_back(declaration.at("cell_bindings").at(k));
    }
    for (std::size_t k = 0; k < bindings.size(); ++k) {
        const auto rule = bindings.at(k).at("evidence");
        const auto kind = text(rule.at("kind"));
        if (kind == "spans") {
            members(rule,
                    {"kind", "document_index_field", "start_field", "end_field", "documents"});
            const auto &ids = column(*endpoint, key(rule.at("document_index_field")), k);
            const auto &starts = column(*endpoint, key(rule.at("start_field")), k);
            const auto &ends = column(*endpoint, key(rule.at("end_field")), k);
            const auto names = rule.at("documents");
            source_array(names, "span document dictionary");
            for (std::size_t j = 0; j < names.size(); ++j)
                require(documents.contains(text(names.at(j))), "unknown span document");
            for (std::size_t j = 0; j < ids.domain().logical_size(); ++j) {
                const auto id = address(ids.state(j).value());
                require(id < names.size(), "span document index outside dictionary");
                const auto start = address(starts.state(j).value()),
                           end = address(ends.state(j).value());
                require(start < end && end <= documents.at(names.at(id).text()).bytes.size(),
                        "invalid complete evidence span");
            }
        } else if (kind == "documents") {
            members(rule, {"kind", "documents"});
            const auto names = rule.at("documents");
            source_array(names, "evidence documents");
            require(names.size() > 0, "evidence needs at least one document");
            for (std::size_t j = 0; j < names.size(); ++j)
                require(documents.contains(text(names.at(j))) &&
                            !documents.at(names.at(j).text()).bytes.empty(),
                        "empty or unknown evidence document");
        } else {
            require(kind == "cofaces", "unsupported compact cell evidence binding");
            members(rule, {"kind", "upper_grade"});
            const auto upper =
                count(rule.at("upper_grade"), bindings.size() - 1, "evidence upper grade");
            require(upper == k + 1 &&
                        bindings.at(upper).at("evidence").at("kind").text() != "cofaces",
                    "coface evidence needs an adjacent directly evidenced grade");
            const auto &arena = endpoint->source.grade({static_cast<std::uint16_t>(upper)});
            for (std::size_t j = 0;
                 j < endpoint->source.grade({static_cast<std::uint16_t>(k)}).cell_count(); ++j) {
                const auto range = arena.coface_incidence_range(j);
                require(range.first < range.second, "cell lacks an original evidence coface");
            }
        }
    }
    const auto declared_contexts = view.at("contexts");
    source_array(declared_contexts, "contexts");
    require(declared_contexts.size() > 0, "study needs a context");
    for (std::size_t i = 0; i < declared_contexts.size(); ++i) {
        const auto item = declared_contexts.at(i);
        members(item, {"id", "block", "metrics", "time"}, {"witness_scope"});
        if (item.has("witness_scope"))
            require(item.at("witness_scope").text() == "block" && item.at("block").raw(),
                    "block witness scope requires a declared strict subcomplex");
        Context context{declaration.at("contexts").at(i), {}};
        if (item.at("block").raw()) {
            context.block = key(item.at("block"));
            require(std::any_of(endpoint->local_blocks.begin(), endpoint->local_blocks.end(),
                                [&](const auto &b) {
                                    return b.block_key == *context.block &&
                                           b.kind == LocalBlockKind::StrictSubcomplex;
                                }),
                    "context block must be a strict subcomplex");
        }
        const auto metrics = item.at("metrics");
        source_array(metrics, "metrics");
        require(metrics.size() == bindings.size(), "one explicit metric per grade required");
        for (std::size_t k = 0; k < metrics.size(); ++k) {
            const auto metric = metrics.at(k);
            const auto kind = text(metric.at("kind"));
            if (kind == "identity")
                members(metric, {"kind"});
            else {
                require(kind == "diagonal_field",
                        "compact diagonal metric must bind a native field");
                members(metric, {"kind", "key"});
                const auto &values = column(*endpoint, key(metric.at("key")), k);
                for (std::size_t j = 0; j < values.domain().logical_size(); ++j)
                    require(values.state(j).value() > 0, "nonpositive metric field");
            }
        }
        const auto time = item.at("time");
        if (time.raw()) {
            members(time, {"axis", "position", "unit", "population"});
            for (const auto *name : {"axis", "unit", "population"})
                (void)text(time.at(name));
            (void)exact_decimal(time.at("position").text());
        }
        require(contexts.emplace(text(item.at("id")), std::move(context)).second,
                "duplicate context");
    }
    const auto fields = view.at("fields");
    source_array(fields, "fields");
    for (std::size_t i = 0; i < fields.size(); ++i) {
        const auto item = fields.at(i);
        members(item, {"id", "key", "context", "channel", "quantity", "units", "method",
                       "value_type", "evidence"});
        const auto identity = key(item.at("key"));
        const auto &field = tensor(*endpoint, identity);
        require(field.role == TemporalTensorRole::Field && field.output,
                "field must name a retained native field");
        require(contexts.contains(text(item.at("context"))), "unknown field context");
        const auto channel = text(item.at("channel"));
        require(channel == "G" || channel == "E" || channel == "T" || channel == "covariate" ||
                    channel == "observation",
                "unknown field channel");
        for (const auto *name : {"quantity", "units", "method"})
            (void)text(item.at(name));
        const auto type = text(item.at("value_type"));
        require(type == "quantitative" || type == "indicator" || type == "categorical",
                "field value type");
        if (type == "indicator")
            for (const auto &entry : field.tensor.values().entries())
                require(entry.value == 0 || entry.value == 1, "indicator field must be zero/one");
        (void)evidence(item.at("evidence"));
        if (!field_pins.contains(identity)) {
            NativeAccessionRequest request;
            request.kind = NativeAccessionKind::Field;
            request.subject = identity;
            field_pins.emplace(identity, directory->pin(request));
        }
        require(this->fields
                    .emplace(text(item.at("id")), Field{declaration.at("fields").at(i), identity})
                    .second,
                "duplicate field id");
    }
    for (std::size_t i = 0; i < design.at("covariates").size(); ++i) {
        const auto field = design.at("covariates").at(i).at("field");
        if (field.raw())
            require(this->fields.contains(field.text()) &&
                        this->fields.at(field.text()).declaration.at("channel").text() ==
                            "covariate",
                    "covariate must name a declared covariate field");
    }
    admit_scoped_actions(view.has("operators") ? view.at("operators") : Json{},
                         view.at("field_sets"));
    for (std::size_t k = 0; k < bindings.size(); ++k)
        grade_pins.push_back(directory->pin(
            NativeAccessionRequest::sparse_chain({static_cast<std::uint16_t>(k)}, {}, {})));
}
std::size_t ScientificStudy::State::grade_count() const {
    return directory->pin_snapshot()->source.grade_count();
}
StableCellKey ScientificStudy::State::coordinate(std::size_t grade, std::size_t index) const {
    if (!compact)
        return coordinate_cells.at(grade).at(index);
    require(grade < grade_count() && index < directory->pin_snapshot()
                                                 ->source.grade({static_cast<std::uint16_t>(grade)})
                                                 .cell_count(),
            "native cell coordinate outside grade");
    return stable_cell_key(directory->pin_snapshot()->source, {{static_cast<std::uint16_t>(grade)},
                                                               static_cast<std::uint32_t>(index)});
}
std::size_t ScientificStudy::State::CoordinateRange::size() const {
    return owner->directory->pin_snapshot()
        ->source.grade({static_cast<std::uint16_t>(grade)})
        .cell_count();
}
ScientificStudy::State::CoordinateRange
ScientificStudy::State::coordinate_keys(std::size_t grade) const {
    require(grade < grade_count(), "source grade absent");
    return {this, grade};
}
std::vector<ScientificStudy::State::CoordinateRange>
ScientificStudy::State::coordinate_grades() const {
    std::vector<CoordinateRange> result;
    result.reserve(grade_count());
    for (std::size_t k = 0; k < grade_count(); ++k)
        result.push_back(coordinate_keys(k));
    return result;
}
std::size_t ScientificStudy::State::cell_count() const {
    if (!compact)
        return cells.size();
    std::size_t count = 0;
    for (std::size_t k = 0; k < grade_count(); ++k)
        count = rfc::machine::checked_size_add(
            count,
            directory->pin_snapshot()->source.grade({static_cast<std::uint16_t>(k)}).cell_count());
    return count;
}
std::size_t ScientificStudy::State::occurrence_count() const {
    if (!compact)
        return occurrences.size();
    std::size_t count = 0;
    for (std::size_t k = 0; k < grade_count(); ++k)
        count = rfc::machine::checked_size_add(count,
                                               directory->pin_snapshot()
                                                   ->source.grade({static_cast<std::uint16_t>(k)})
                                                   .incidence_count());
    return count;
}
StableCellKey ScientificStudy::State::cell_key_at(std::size_t index) const {
    if (!compact)
        return cell_keys.at(index);
    for (std::size_t k = 0; k < grade_count(); ++k) {
        const auto size =
            directory->pin_snapshot()->source.grade({static_cast<std::uint16_t>(k)}).cell_count();
        if (index < size)
            return coordinate(k, index);
        index -= size;
    }
    throw std::out_of_range("native cell page coordinate");
}
StableOccurrenceKey ScientificStudy::State::occurrence_key_at(std::size_t index) const {
    if (!compact)
        return occurrence_keys.at(index);
    for (std::size_t k = 0; k < grade_count(); ++k) {
        const auto &arena =
            directory->pin_snapshot()->source.grade({static_cast<std::uint16_t>(k)});
        if (index < arena.incidence_count())
            return {{static_cast<std::uint16_t>(k)}, arena.occurrence_identities.at(index)};
        index -= arena.incidence_count();
    }
    throw std::out_of_range("native occurrence page coordinate");
}
StableCellKey ScientificStudy::State::resolve_cell(const std::string &id) const {
    if (!compact) {
        const auto found = cell_ids.find(id);
        require(found != cell_ids.end(), "unknown cell identity");
        return found->second;
    }
    for (std::size_t k = 0; k < cell_bindings.size(); ++k) {
        const auto prefix = cell_bindings[k].at("id_prefix").text();
        if (!id.starts_with(prefix))
            continue;
        const auto number = std::string_view(id).substr(prefix.size());
        std::uint64_t key{};
        const auto [end, error] =
            std::from_chars(number.data(), number.data() + number.size(), key);
        require(error == std::errc{} && end == number.data() + number.size() && key &&
                    number == std::to_string(key),
                "noncanonical native cell label");
        const StableCellKey cell{{static_cast<std::uint16_t>(k)}, key};
        (void)directory->cell(cell);
        return cell;
    }
    throw std::invalid_argument("unknown native cell label");
}
bool ScientificStudy::State::has_cell(const std::string &id) const {
    try {
        (void)resolve_cell(id);
        return true;
    } catch (const std::invalid_argument &) {
        return false;
    } catch (const std::out_of_range &) {
        return false;
    }
}
Json ScientificStudy::State::compact_evidence(std::size_t grade, std::size_t coordinate_) const {
    const auto rule = cell_bindings.at(grade).at("evidence");
    const auto kind = rule.at("kind").text();
    auto out = Json::array();
    const auto endpoint = directory->pin_snapshot();
    if (kind == "spans") {
        const auto &ids = column(
            *endpoint, std::stoull(rule.at("document_index_field").text()), grade);
        const auto &starts =
            column(*endpoint, std::stoull(rule.at("start_field").text()), grade);
        const auto &ends =
            column(*endpoint, std::stoull(rule.at("end_field").text()), grade);
        out.push(reference(rule.at("documents").at(address(ids.state(coordinate_).value())).text(),
                           address(starts.state(coordinate_).value()),
                           address(ends.state(coordinate_).value())));
    } else if (kind == "documents") {
        const auto names = rule.at("documents");
        for (std::size_t i = 0; i < names.size(); ++i) {
            const auto id = names.at(i).text();
            out.push(reference(id, 0, documents.at(id).bytes.size()));
        }
    } else {
        const auto upper = std::stoull(rule.at("upper_grade").integer().get_str());
        const auto &arena = endpoint->source.grade({static_cast<std::uint16_t>(upper)});
        const auto [start, end] = arena.coface_incidence_range(coordinate_);
        std::optional<std::uint32_t> previous;
        for (auto at = start; at < end; ++at) {
            const auto parent = arena.upper_index.at(arena.incidence_by_lower.at(at));
            if (previous && *previous == parent)
                continue;
            previous = parent;
            const auto spans = compact_evidence(upper, parent);
            for (std::size_t i = 0; i < spans.size(); ++i)
                out.push(spans.at(i));
        }
    }
    return out;
}
Json ScientificStudy::State::cell_record(StableCellKey cell) const {
    if (!compact)
        return cells.at(cell);
    const auto native = directory->cell(cell);
    const auto grade = cell.grade.value;
    const auto &binding = cell_bindings.at(grade);
    auto out = Json::object();
    out.set("grade", Json::integer(grade));
    out.set("key", Json::string(std::to_string(cell.key)));
    out.set("id", Json::string(binding.at("id_prefix").text() + std::to_string(cell.key)));
    out.set("kind", Json::string(binding.at("kind").text()));
    out.set("label", Json::string(binding.at("label_prefix").text() + std::to_string(cell.key)));
    (void)native;
    const auto &arena = directory->pin_snapshot()->source.grade(cell.grade);
    std::size_t position = 0;
    if (cell.key <= arena.cell_count() && coordinate(grade, cell.key - 1) == cell)
        position = cell.key - 1;
    else {
        while (position < arena.cell_count() && coordinate(grade, position) != cell)
            ++position;
    }
    require(position < arena.cell_count(), "native cell lacks a source coordinate");
    out.set("evidence", compact_evidence(grade, position));
    return out;
}
Json ScientificStudy::State::occurrence_record(StableOccurrenceKey occurrence) const {
    if (!compact)
        return occurrences.at(occurrence);
    const auto native = directory->occurrence(occurrence);
    auto out = Json::object();
    out.set("grade", Json::integer(occurrence.upper_grade.value));
    out.set("key", Json::string(std::to_string(occurrence.key.value)));
    out.set("evidence", cell_record(native.parent()).at("evidence"));
    return out;
}
Rational ScientificStudy::State::metric_value(const RetainedRecord &metric,
                                              std::size_t coordinate_) const {
    if (metric.at("kind").text() == "diagonal")
        return exact_rational(metric.at("values").at(coordinate_).text());
    return tensor(*directory->pin_snapshot(), std::stoull(metric.at("key").text()))
        .tensor.state(coordinate_)
        .value();
}
} // namespace gadget
