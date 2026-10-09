#include "scientific_native_novelty.hpp"
#include "scientific_detail.hpp"

namespace gadget::scientific_detail {
namespace {
using namespace rfc;
using Pin = storage::PinnedNativeAccession;
using Plan = inference::PreparedNativeMetricNovelty;
Json vector(const numeric::RationalVector &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(Json::string(value.get_str()));
    return out;
}
Json matrix(const numeric::ExactMatrix &value) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(value.rows()));
    out.set("columns", cardinality(value.cols()));
    out.set("complete_support", Json::boolean(true));
    for (const auto &[position, coefficient] : value.entries()) {
        auto item = Json::object();
        item.set("row", cardinality(position.first));
        item.set("column", cardinality(position.second));
        item.set("value", Json::string(coefficient.get_str()));
        entries.push(item);
    }
    out.set("entries", entries);
    return out;
}
Json requirements(const runtime::ExecutionRequirements &value) {
    auto out = Json::object();
    out.set("records", Json::string(std::to_string(value.records)));
    out.set("rational_slots", Json::string(std::to_string(value.rational_slots)));
    out.set("sparse_entries", Json::string(std::to_string(value.sparse_entries)));
    out.set("scalar_work", Json::string(std::to_string(value.scalar_work)));
    out.set("integer_bits", Json::string(std::to_string(value.integer_bits)));
    return out;
}
void complete_field(const Pin &field) {
    const auto carrier = field.carrier();
    require(field.kind() == storage::NativeAccessionKind::Field &&
                field.tensor_declaration()->role == storage::TemporalTensorRole::Field &&
                carrier->kind == storage::TemporalCarrierKind::SourceGrade &&
                carrier->domain == numeric::IndexDomain::linear(carrier->domain.logical_size()) &&
                field.tensor()->domain() == carrier->domain &&
                field.tensor()->support().count() == carrier->domain.logical_size(),
            "novelty requires a complete Field on its actual linear SourceGrade carrier");
}
void endpoint(const Json &declaration, const Pin &field, bool current) {
    members(declaration, {"publication", "source_sequence", "field"},
            current ? std::set<std::string>{"metric"} : std::set<std::string>{});
    (void)text(declaration.at("publication"));
    (void)text(declaration.at("field"));
    require(declaration.at("source_sequence").text() == std::to_string(field.sequence()),
            "novelty source_sequence differs from actual publication");
    if (current)
        (void)text(declaration.at("metric"));
    complete_field(field);
}
inference::NativeMetricNoveltyRequest request(const Json &query, const NativeNoveltyBinding &b) {
    members(query, {"operation", "old", "new", "comparison"}, {"causes", "preparation"});
    require(query.at("operation").text() == "metric_relative_novelty", "wrong novelty operation");
    if (query.has("preparation"))
        require(query.at("preparation").text() == "reuse" ||
                    query.at("preparation").text() == "ephemeral",
                "novelty preparation must be reuse or ephemeral");
    endpoint(query.at("old"), b.old_field, false);
    endpoint(query.at("new"), b.new_field, true);
    const auto old = b.old_field.carrier(), current = b.new_field.carrier();
    require(old->grade == current->grade, "novelty carrier grades differ");
    require(b.metric.same_execution_authority(b.new_field) &&
                b.metric.endpoint().get() == b.new_field.endpoint().get(),
            "novelty metric requires the actual new publication");
    const auto metric = b.metric.tensor_declaration();
    const auto n = old->domain.logical_size(), m = current->domain.logical_size();
    require(b.metric.kind() == storage::NativeAccessionKind::Tensor &&
                metric->role == storage::TemporalTensorRole::Metric &&
                metric->input == current->key && metric->output == current->key &&
                b.metric.carrier().get() == current.get() &&
                metric->tensor.domain().extents() == std::vector<std::size_t>{m, m} &&
                metric->tensor.support().count() == metric->tensor.domain().logical_size(),
            "novelty metric must be a complete bilinear form on the actual new carrier");
    const auto comparison = query.at("comparison");
    members(comparison, {"old_carrier", "new_carrier", "map", "meaning", "evidence"});
    require(key(comparison.at("old_carrier")) == old->key &&
                key(comparison.at("new_carrier")) == current->key,
            "novelty comparison names different actual carriers");
    (void)text(comparison.at("meaning"));
    std::size_t channels = 0;
    if (query.has("causes")) {
        source_array(query.at("causes"), "novelty cause word");
        channels = query.at("causes").size();
        require(channels == b.causes.size(), "novelty cause bindings differ from ordered word");
        for (std::size_t i = 0; i < channels; ++i) {
            (void)text(query.at("causes").at(i));
            complete_field(b.causes[i]);
            require(b.causes[i].same_execution_authority(b.new_field) &&
                        b.causes[i].endpoint().get() == b.new_field.endpoint().get() &&
                        b.causes[i].carrier().get() == current.get(),
                    "novelty cause requires the actual complete new carrier");
        }
    }
    // Mathematical dimensions are the complete source and requested cause axes.
    // Machine addressability permits reuse for later words; the current lease
    // admits their complete native workspace, bit growth and scalar work.
    return {b.old_field,
            b.new_field,
            old->key,
            current->key,
            declared_comparison(comparison.at("map"), m, n),
            b.metric,
            std::numeric_limits<std::size_t>::max()};
}
inference::NativeMetricAvailability availability(const NativeNoveltyBinding &b) {
    return {{b.old_field.source_namespace(), b.old_field.sequence()},
            {b.new_field.source_namespace(), b.new_field.sequence()}};
}
std::string cache_key(const inference::NativeMetricNoveltyRequest &r) {
    // The digest indexes a candidate only; complete map bytes and actual owners
    // are checked before reuse. Do not retain a second comparison/quotation copy.
    std::string declaration =
        std::to_string(r.old_publication.source_namespace()) + ":" +
        std::to_string(r.old_publication.sequence()) + ":" + std::to_string(r.old_carrier) + ":" +
        std::to_string(r.new_publication.source_namespace()) + ":" +
        std::to_string(r.new_publication.sequence()) + ":" + std::to_string(r.new_carrier) + ":" +
        std::to_string(r.new_metric.tensor_declaration()->key) + ":";
    // Private retained-owner addresses separate independent publications. They
    // are cache addresses, never mathematical identities or public receipts.
    declaration +=
        std::to_string(reinterpret_cast<std::uintptr_t>(r.old_publication.endpoint().get())) + ":";
    declaration +=
        std::to_string(reinterpret_cast<std::uintptr_t>(r.new_publication.endpoint().get())) + ":";
    for (const auto extent : r.comparison.domain().extents())
        declaration += std::to_string(extent) + ":";
    for (const auto &e : r.comparison.values().entries())
        declaration += std::to_string(e.index) + "=" + e.value.get_str() + ";";
    return sha256(declaration);
}
bool same_request(const inference::NativeMetricNoveltyRequest &a,
                  const inference::NativeMetricNoveltyRequest &b) {
    const auto same = [](const Pin &x, const Pin &y) {
        return x.same_execution_authority(y) && x.endpoint().get() == y.endpoint().get();
    };
    if (!same(a.old_publication, b.old_publication) ||
        !same(a.new_publication, b.new_publication) || !same(a.new_metric, b.new_metric) ||
        a.old_carrier != b.old_carrier || a.new_carrier != b.new_carrier ||
        a.new_metric.tensor().get() != b.new_metric.tensor().get() ||
        a.comparison.domain() != b.comparison.domain())
        return false;
    const auto &x = a.comparison.values().entries(), &y = b.comparison.values().entries();
    if (x.size() != y.size())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i].index != y[i].index || x[i].value != y[i].value)
            return false;
    return true;
}
Json witness(const inference::NativeMetricEndpointWitness &w) {
    auto out = Json::object();
    out.set("source_namespace", Json::string(std::to_string(w.source_namespace)));
    out.set("source_sequence", Json::string(std::to_string(w.sequence)));
    out.set("carrier", Json::string(std::to_string(w.carrier)));
    out.set("field", Json::string(std::to_string(w.input_field)));
    out.set("grade", cardinality(w.grade.value));
    return out;
}
Json available(const inference::NativeMetricAvailability &a) {
    auto out = Json::object();
    for (const auto &[name, endpoint] :
         {std::pair{"old", a.old_endpoint}, std::pair{"new", a.new_endpoint}}) {
        auto item = Json::object();
        item.set("source_namespace", Json::string(std::to_string(endpoint.source_namespace)));
        item.set("maximum_sequence", Json::string(std::to_string(endpoint.max_sequence)));
        out.set(name, item);
    }
    return out;
}
} // namespace
void NativeNoveltyPlans::validate(const Json &query, const NativeNoveltyBinding &binding) {
    (void)request(query, binding);
}
Json NativeNoveltyPlans::execute(const Json &query, const NativeNoveltyBinding &b,
                                 runtime::ExecutionLease lease, CpuResources &cpu) {
    lease.check_cancelled();
    auto reservation = cpu.reserve(1);
    auto r = request(query, b);
    const auto a = availability(b);
    const auto id = cache_key(r);
    const bool retain = !query.has("preparation") || query.at("preparation").text() == "reuse";
    auto cached = plans_.find(id);
    const bool reuse = retain && cached != plans_.end() &&
                       same_request(cached->second.request(), r) &&
                       cached->second.request().dimension_limit >= r.dimension_limit;
    // A digest collision or independent equal-key owner cannot borrow this plan.
    auto plan = reuse ? cached->second : Plan::prepare(std::move(r), a, lease);
    auto result = query.has("causes")
                      ? plan.execute_causes(b.old_field, b.new_field, b.causes, a, lease)
                      : plan.execute(b.old_field, b.new_field, a, lease);
    const auto &v = result.values();
    auto out = Json::object(), split = Json::object(), work = Json::object();
    auto publications = Json::object();
    for (const auto *name : {"old", "new", "evidence_profile"})
        publications.set(name, b.declarations.at(name));
    out.set("publications", publications);
    split.set("transported_old", vector(v.transported_old));
    split.set("innovation", vector(v.innovation));
    split.set("persistent_coordinates", vector(v.persistent_coordinates));
    split.set("persistent", vector(v.persistent));
    split.set("complement", vector(v.complement));
    split.set("reconstructed", vector(v.reconstructed));
    split.set("pairings", matrix(v.pairings));
    split.set("quadrance", Json::string(v.quadrance.get_str()));
    out.set("split", split);
    if (result.cause_pairings()) {
        const auto &c = *result.cause_pairings();
        auto causes = Json::object(), dictionary = Json::array(), word = Json::array();
        std::map<std::uint64_t, std::size_t> seen;
        for (std::size_t i = 0; i < b.causes.size(); ++i) {
            const auto key = b.causes[i].tensor_declaration()->key;
            auto found = seen.find(key);
            if (found == seen.end()) {
                auto item = Json::object();
                item.set("field_key", Json::string(std::to_string(key)));
                item.set("values", vector(c.channels[i]));
                found = seen.emplace(key, dictionary.size()).first;
                dictionary.push(item);
            }
            auto slot = Json::object();
            slot.set("field", query.at("causes").at(i));
            slot.set("dictionary_index", cardinality(found->second));
            word.push(slot);
        }
        causes.set("fields", dictionary);
        causes.set("meanings", b.declarations.at("cause_fields"));
        causes.set("word", word);
        causes.set("aggregate", vector(c.aggregate));
        causes.set("pairings", matrix(c.pairings));
        causes.set("quadrance", Json::string(c.quadrance.get_str()));
        out.set("causes", causes);
    }
    const auto &p = result.provenance();
    work.set("profile", Json::string(p.profile));
    work.set("old", witness(p.old_endpoint));
    work.set("new", witness(p.new_endpoint));
    work.set("metric", Json::string(std::to_string(p.new_metric)));
    work.set("declared_causes", Json::boolean(p.declared_causes));
    work.set("selected_requirements", requirements(p.selected_requirements));
    work.set("preparation_availability", available(p.preparation_availability));
    work.set("execution_availability", available(p.execution_availability));
    work.set("preparation_policy", Json::string(p.preparation_policy));
    work.set("execution_policy", Json::string(p.execution_policy));
    work.set("preparation_steps", cardinality(p.preparation_steps));
    work.set("execution_steps", cardinality(p.execution_steps));
    work.set("cpu_fallback", Json::boolean(p.fallback));
    auto limits = Json::object();
    limits.set("dimension", Json::string(std::to_string(p.action_limits.dimension)));
    limits.set("entries", Json::string(std::to_string(p.action_limits.entries)));
    limits.set("work", Json::string(std::to_string(p.action_limits.steps)));
    limits.set("integer_bits", Json::string(std::to_string(p.action_limits.integer_bits)));
    work.set("action_limits", limits);
    out.set("native_execution", work);
    auto storage = Json::object();
    storage.set("preparation", Json::string(retain ? "reuse" : "ephemeral"));
    storage.set("retained_requirements", requirements(plan.retained_requirements()));
    storage.set("numeric_entries", cardinality(plan.numeric_plan().stored_entries()));
    storage.set("factor_index_entries", cardinality(plan.numeric_plan().stored_index_entries()));
    storage.set("dimension_admission",
                Json::string(std::to_string(plan.request().dimension_limit)));
    storage.set("old_dimension", cardinality(plan.numeric_plan().old_dimension()));
    storage.set("new_dimension", cardinality(plan.numeric_plan().new_dimension()));
    storage.set("cause_axes", cardinality(b.causes.size()));
    storage.set("counts_scope",
                Json::string("native named containers; excludes allocator and GMP/RSS"));
    out.set("native_storage", storage);
    auto detached = Json::parse(out.dump());
    lease.check_cancelled();
    const auto next = machine::checked_size_add(applications_, 1);
    const auto preparations = machine::checked_size_add(preparations_, reuse ? 0 : 1);
    if (retain && !reuse)
        plans_.insert_or_assign(id, plan);
    applications_ = next;
    preparations_ = preparations;
    return detached;
}
Json NativeNoveltyPlans::statistics() const {
    auto out = Json::object();
    std::size_t scalars = 0, indices = 0;
    for (const auto &[id, plan] : plans_) {
        (void)id;
        scalars = machine::checked_size_add(scalars, plan.numeric_plan().stored_entries());
        indices = machine::checked_size_add(indices, plan.numeric_plan().stored_index_entries());
    }
    out.set("retained_plans", cardinality(plans_.size()));
    out.set("preparations", cardinality(preparations_));
    out.set("applications", cardinality(applications_));
    out.set("numeric_entries", cardinality(scalars));
    out.set("factor_index_entries", cardinality(indices));
    return out;
}
void NativeNoveltyPlans::clear() { plans_.clear(); }
} // namespace gadget::scientific_detail
