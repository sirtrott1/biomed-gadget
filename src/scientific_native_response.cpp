#include "scientific_native_response.hpp"
#include "scientific_detail.hpp"
#include <algorithm>
#include <charconv>
#include <rfc/numeric/checked_exact_actions.hpp>

namespace gadget::scientific_detail {
namespace {
using namespace rfc;
using Pin = storage::PinnedNativeAccession;
using Matrix = numeric::ExactMatrix;
using Correspondence = inference::NativeGreenCorrespondence;
Json vec(const Vector &v) {
    auto out = Json::array();
    for (const auto &q : v)
        out.push(Json::string(q.get_str()));
    return out;
}
Json mat(const Matrix &m) {
    auto out = Json::object(), entries = Json::array();
    out.set("rows", cardinality(m.rows()));
    out.set("columns", cardinality(m.cols()));
    for (const auto &[p, q] : m.entries()) {
        auto e = Json::object();
        e.set("row", cardinality(p.first));
        e.set("column", cardinality(p.second));
        e.set("value", Json::string(q.get_str()));
        entries.push(e);
    }
    out.set("entries", entries);
    return out;
}
std::uint64_t sequence(const Json &v) {
    const auto s = v.text();
    std::uint64_t n{};
    const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), n);
    require(!s.empty() && (s.size() == 1 || s.front() != '0') && error == std::errc{} &&
                end == s.data() + s.size(),
            "source_sequence must be a canonical unsigned decimal");
    return n;
}
const storage::TemporalCarrier &carrier(const Pin &p, std::uint64_t id) {
    const auto endpoint = p.endpoint();
    for (const auto &c : endpoint->carriers)
        if (c.key == id) {
            require(c.domain.rank() == 1, "response carrier requires declared linear coordinates");
            return c;
        }
    throw std::invalid_argument("response carrier absent from actual publication");
}
void complete(const Pin &p, bool form = false) {
    require(p.kind() == (form ? storage::NativeAccessionKind::Tensor
                              : storage::NativeAccessionKind::Operator) &&
                p.tensor()->domain().rank() == 2 &&
                p.tensor()->support().count() == p.tensor()->domain().logical_size(),
            "response requires a complete stored operator");
}
void same(const Pin &a, const Pin &b) {
    require(a.same_execution_authority(b) && a.endpoint().get() == b.endpoint().get(),
            "response operands must share their actual publication and execution authority");
}
Correspondence comparison(const Json &d, const Pin &a, const Pin &b, std::uint64_t old_id,
                          std::uint64_t new_id) {
    members(d, {"old_carrier", "new_carrier", "map", "meaning", "evidence"});
    require(key(d.at("old_carrier")) == old_id && key(d.at("new_carrier")) == new_id,
            "response comparison binds different actual carriers");
    (void)text(d.at("meaning"), source_byte_limit());
    const auto &ac = carrier(a, old_id), &bc = carrier(b, new_id);
    return {a, b, old_id, new_id,
            declared_comparison(d.at("map"), bc.domain.logical_size(), ac.domain.logical_size())};
}
struct Request {
    inference::NativeGreenResponseRequest green;
    inference::NativeGreenAvailability availability;
    std::optional<Correspondence> output;
};
Request request(const Json &query, const DeclaredResponseBinding &b) {
    members(query, {"operation", "old", "new", "lambda", "comparisons"},
            {"observation", "pairing"});
    require(query.at("operation").text() == "declared_response_transition",
            "wrong response operation");
    const auto endpoint = [](const Json &d, const Pin &L, const Pin &I, const Pin &u) {
        members(d, {"publication", "resolving", "realization", "field", "source_sequence"});
        for (const auto *name : {"publication", "resolving", "realization", "field"})
            (void)text(d.at(name));
        require(sequence(d.at("source_sequence")) == L.sequence(),
                "response source_sequence differs from selected publication");
        same(L, I);
        same(L, u);
        complete(L);
        complete(I);
        const auto l = L.tensor_declaration(), injection = I.tensor_declaration(),
                   field = u.tensor_declaration();
        require(l->input == l->output &&
                    L.carrier()->kind == storage::TemporalCarrierKind::SourceGrade,
                "response resolving map must be a SourceGrade endomorphism");
        const auto &ambient = carrier(L, *l->output);
        require(injection->role == storage::TemporalTensorRole::SourceRealization &&
                    injection->output == ambient.key,
                "response realization must bind its actual ambient carrier");
        const auto &source = carrier(I, *injection->input);
        require(u.kind() == storage::NativeAccessionKind::Field &&
                    field->role == storage::TemporalTensorRole::Field &&
                    field->output == source.key && field->tensor.domain() == source.domain &&
                    field->tensor.support().count() == source.domain.logical_size(),
                "response requires full stored field support on the actual source carrier");
    };
    endpoint(query.at("old"), b.old_operator, b.old_realization, b.old_field);
    endpoint(query.at("new"), b.new_operator, b.new_realization, b.new_field);
    const auto maps = query.at("comparisons");
    members(maps, {"ambient", "source"});
    const auto old_x = *b.old_operator.tensor_declaration()->output,
               new_x = *b.new_operator.tensor_declaration()->output,
               old_u = *b.old_realization.tensor_declaration()->input,
               new_u = *b.new_realization.tensor_declaration()->input;
    std::size_t dimension = 1; // The native contract is positive even for empty carriers.
    for (const auto *p : {&b.old_operator, &b.new_operator})
        for (const auto &c : p->endpoint()->carriers)
            if (c.domain.rank() == 1)
                dimension = std::max(dimension, c.domain.logical_size());
    auto J = comparison(maps.at("ambient"), b.old_operator, b.new_operator, old_x, new_x),
         K = comparison(maps.at("source"), b.old_operator, b.new_operator, old_u, new_u);
    Request out{{b.old_operator, b.new_operator, b.old_realization, b.new_realization, std::move(J),
                 std::move(K), exact_rational(query.at("lambda").text()), dimension},
                {{b.old_operator.source_namespace(), b.old_operator.sequence()},
                 {b.new_operator.source_namespace(), b.new_operator.sequence()}},
                {}};
    auto target = new_x;
    if (query.has("observation")) {
        const auto d = query.at("observation");
        members(d, {"old", "new", "comparison"});
        require(b.old_accession.has_value() && b.new_accession.has_value(),
                "response accession missing");
        for (const auto &p : {*b.old_accession, *b.new_accession}) {
            complete(p);
            require(p.tensor_declaration()->role == storage::TemporalTensorRole::Accession,
                    "response observation requires an actual Accession");
            dimension = std::max(dimension, p.carrier()->domain.logical_size());
        }
        same(*b.old_accession, b.old_operator);
        same(*b.new_accession, b.new_operator);
        require(b.old_accession->tensor_declaration()->input == old_x &&
                    b.new_accession->tensor_declaration()->input == new_x,
                "response accession input differs from actual ambient carrier");
        out.output = comparison(d.at("comparison"), b.old_operator, b.new_operator,
                                *b.old_accession->tensor_declaration()->output,
                                *b.new_accession->tensor_declaration()->output);
        out.green.dimension_limit = dimension;
        target = *b.new_accession->tensor_declaration()->output;
    }
    if (query.has("pairing")) {
        members(query.at("pairing"), {"metric"});
        require(b.metric.has_value(), "response metric missing");
        complete(*b.metric, true);
        same(*b.metric, b.new_operator);
        const auto m = b.metric->tensor_declaration();
        require(m->role == storage::TemporalTensorRole::Metric && m->input == target &&
                    m->output == target,
                "cause pairing requires a declared Metric on the actual output carrier");
    }
    return out;
}
Json endpoint_witness(const inference::NativeGreenEndpointWitness &w) {
    auto out = Json::object();
    for (const auto &[name, value] : std::vector<std::pair<std::string, std::uint64_t>>{
             {"source_namespace", w.source_namespace},
             {"sequence", w.sequence},
             {"ambient_carrier", w.ambient_carrier},
             {"source_carrier", w.source_carrier},
             {"resolving_operator", w.resolving_operator},
             {"realization", w.injection},
             {"input_field", w.input_field}})
        out.set(name, Json::string(std::to_string(value)));
    if (w.output_carrier)
        out.set("output_carrier", Json::string(std::to_string(*w.output_carrier)));
    if (w.accession)
        out.set("accession", Json::string(std::to_string(*w.accession)));
    return out;
}
Json availability(const inference::NativeGreenAvailability &a) {
    auto out = Json::object();
    for (const auto &[name, v] : std::vector<std::pair<std::string, runtime::WorkAvailability>>{
             {"old", a.old_endpoint}, {"new", a.new_endpoint}}) {
        auto e = Json::object();
        e.set("source_namespace", Json::string(std::to_string(v.source_namespace)));
        e.set("source_sequence", Json::string(std::to_string(v.max_sequence)));
        out.set(name, e);
    }
    return out;
}
Json provenance(const inference::NativeGreenProvenance &p) {
    auto out = Json::object();
    out.set("profile", Json::string(p.profile));
    out.set("old", endpoint_witness(p.old_endpoint));
    out.set("new", endpoint_witness(p.new_endpoint));
    out.set("lambda", Json::string(p.lambda.get_str()));
    out.set("green_preparation_policy", Json::string(p.green_preparation_policy));
    out.set("preparation_policy", Json::string(p.preparation_policy));
    out.set("execution_policy", Json::string(p.execution_policy));
    out.set("green_preparation_steps", cardinality(p.green_preparation_steps));
    out.set("preparation_steps", cardinality(p.preparation_steps));
    out.set("execution_steps", cardinality(p.execution_steps));
    out.set("green_preparation_availability", availability(p.green_preparation_availability));
    out.set("preparation_availability", availability(p.preparation_availability));
    out.set("execution_availability", availability(p.execution_availability));
    out.set("observed", Json::boolean(p.observed));
    out.set("reused_green", Json::boolean(p.reused_green));
    out.set("fallback", Json::boolean(p.fallback));
    auto limits = Json::object(), selected = Json::object();
    for (const auto &[name, n] : std::vector<std::pair<std::string, std::size_t>>{
             {"dimension", p.action_limits.dimension},
             {"entries", p.action_limits.entries},
             {"steps", p.action_limits.steps},
             {"integer_bits", p.action_limits.integer_bits}})
        // Host bounds may span the unsigned machine range, beyond json-c's
        // signed integer representation. Preserve them as decimal strings.
        limits.set(name, Json::string(std::to_string(n)));
    for (const auto &[name, n] : std::vector<std::pair<std::string, std::size_t>>{
             {"records", p.selected_requirements.records},
             {"rational_slots", p.selected_requirements.rational_slots},
             {"sparse_entries", p.selected_requirements.sparse_entries},
             {"scalar_work", p.selected_requirements.scalar_work},
             {"integer_bits", p.selected_requirements.integer_bits}})
        selected.set(name, cardinality(n));
    out.set("action_limits", limits);
    out.set("selected_requirements", selected);
    return out;
}
Matrix stored_map(const Pin &pin, numeric::ExactActionGuard &guard) {
    const auto t = pin.tensor();
    Matrix out(t->domain().extents()[0], t->domain().extents()[1]);
    for (const auto &e : t->values().entries()) {
        guard.step();
        guard.value(e.value);
        out.set(e.index / out.cols(), e.index % out.cols(), e.value);
    }
    guard.matrix(out);
    return out;
}
Json metric_pairing(const Pin &metric, const std::vector<Vector> &causes,
                    const runtime::ExecutionLease &lease, std::size_t dimension) {
    const auto &host = lease.policy().limits;
    numeric::ExactActionGuard guard({std::max(dimension, causes.size()),
                                     std::min(host.rational_slots, host.sparse_entries),
                                     host.scalar_work, host.integer_bits},
                                    [&lease] { lease.check_cancelled(); });
    const auto m = stored_map(metric, guard);
    Matrix columns(m.rows(), causes.size());
    for (std::size_t j = 0; j < causes.size(); ++j) {
        guard.vector(causes[j]);
        require(causes[j].size() == m.rows(), "cause metric dimension mismatch");
        for (std::size_t i = 0; i < causes[j].size(); ++i)
            columns.set(i, j, causes[j][i]);
    }
    const auto gram = numeric::checked_multiply(
        columns.transpose(), numeric::checked_multiply(m, columns, guard), guard);
    Rational total{0};
    for (const auto &[p, q] : gram.entries()) {
        (void)p;
        total = guard.add(total, q);
    }
    auto out = Json::object();
    out.set("metric_key", Json::string(std::to_string(metric.tensor_declaration()->key)));
    out.set("gram", mat(gram));
    out.set("aggregate_quadrance", Json::string(total.get_str()));
    out.set("execution_steps", cardinality(guard.steps()));
    return out;
}
} // namespace
runtime::ExecutionLease declared_response_lease() {
    return source_execution_lease("gadget-native-response-machine-addressable-v1");
}
runtime::ExecutionLease source_execution_lease(const std::string &name) {
    runtime::ExecutionHostPolicy policy;
    policy.name = name;
    const auto n = std::numeric_limits<std::size_t>::max();
    policy.limits = {n, n, n, n, n};
    return runtime::ExecutionLease(policy);
}
void DeclaredResponsePlans::validate(const Json &q, const DeclaredResponseBinding &b) {
    (void)request(q, b);
}
Json DeclaredResponsePlans::execute(const Json &q, const DeclaredResponseBinding &b,
                                    runtime::ExecutionLease lease) {
    lease.check_cancelled();
    const auto r = request(q, b);
    auto descriptor = Json::object();
    for (const auto *name : {"old", "new"}) {
        auto e = Json::object();
        for (const auto *k : {"publication", "resolving", "realization", "source_sequence"})
            e.set(k, q.at(name).at(k));
        descriptor.set(name, e);
    }
    descriptor.set("comparisons", q.at("comparisons"));
    descriptor.set("lambda", Json::string(r.green.lambda.get_str()));
    const auto base_key = descriptor.dump_compact();
    if (!green_.contains(base_key))
        green_.emplace(base_key, inference::PreparedNativeGreenResponseDelta::prepare(
                                     r.green, r.availability, lease));
    const auto &green = green_.at(base_key);
    auto out = Json::object(), causes = Json::object(), names = Json::array();
    std::vector<Vector> vectors;
    auto add = [&](const std::string &name, const Vector &v) {
        names.push(Json::string(name));
        causes.set(name, vec(v));
        vectors.push_back(v);
    };
    if (r.output) {
        descriptor.set("observation", q.at("observation"));
        const auto obs_key = descriptor.dump_compact();
        if (!observed_.contains(obs_key))
            observed_.emplace(obs_key, inference::PreparedNativeObservedGreenResponseDelta::prepare(
                                           green, *b.old_accession, *b.new_accession, *r.output,
                                           r.availability, lease));
        const auto result =
            observed_.at(obs_key).execute(b.old_field, b.new_field, r.availability, lease);
        const auto &v = result.values();
        add("observation_change", v.accession_term);
        add("realization_change", v.realization_term);
        add("source_innovation", v.innovation_term);
        add("operator_change", v.operator_term);
        out.set("accession_defect", mat(v.accession_defect));
        out.set("aggregate", vec(v.aggregate));
        out.set("native_execution", provenance(result.provenance()));
    } else {
        const auto result = green.execute(b.old_field, b.new_field, r.availability, lease);
        const auto &v = result.values();
        add("realization_change", v.realization_term);
        add("source_innovation", v.innovation_term);
        add("operator_change", v.operator_term);
        out.set("old_response", vec(v.old_response));
        out.set("new_response", vec(v.new_response));
        out.set("coordinate_innovation", vec(v.coordinate_innovation));
        out.set("aggregate", vec(v.aggregate));
        out.set("native_execution", provenance(result.provenance()));
    }
    out.set("realization_defect", mat(green.numeric_plan().injection_defect()));
    out.set("operator_defect", mat(green.numeric_plan().operator_defect()));
    out.set("cause_order", names);
    out.set("causes", causes);
    out.set("publications", b.publications);
    auto maps = Json::object();
    for (const auto &[name, p] :
         std::vector<std::pair<std::string, Pin>>{{"old_resolving", b.old_operator},
                                                  {"new_resolving", b.new_operator},
                                                  {"old_realization", b.old_realization},
                                                  {"new_realization", b.new_realization}})
        maps.set(name, mat(storage::temporal_tensor_matrix(*p.tensor_declaration())));
    if (b.old_accession) {
        maps.set("old_accession",
                 mat(storage::temporal_tensor_matrix(*b.old_accession->tensor_declaration())));
        maps.set("new_accession",
                 mat(storage::temporal_tensor_matrix(*b.new_accession->tensor_declaration())));
    }
    if (b.metric)
        maps.set("metric", mat(storage::temporal_tensor_matrix(*b.metric->tensor_declaration())));
    out.set("declared_maps", maps);
    out.set("identity_verified", Json::boolean(true));
    if (b.metric)
        out.set("pairing", metric_pairing(*b.metric, vectors, lease, r.green.dimension_limit));
    lease.check_cancelled();
    applications_ = machine::checked_size_add(applications_, 1);
    return out;
}
Json DeclaredResponsePlans::statistics() const {
    auto out = Json::object();
    out.set("green_pair_preparations", cardinality(green_.size()));
    out.set("observation_preparations", cardinality(observed_.size()));
    out.set("applications", cardinality(applications_));
    return out;
}
} // namespace gadget::scientific_detail
