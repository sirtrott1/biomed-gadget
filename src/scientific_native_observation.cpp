#include "scientific_native_observation.hpp"
#include "scientific_detail.hpp"
#include <rfc/runtime/cpu_tiles.hpp>
#include <type_traits>

namespace gadget::scientific_detail {
namespace {
using namespace rfc;
using Pin = storage::PinnedNativeAccession;
using Base = inference::PreparedNativeWoodburySolve::Base;
Json indices(const std::vector<std::size_t> &values) {
    auto out = Json::array();
    for (auto value : values)
        out.push(cardinality(value));
    return out;
}
Json vector(const numeric::RationalVector &values) {
    auto out = Json::array();
    for (const auto &value : values)
        out.push(Json::string(value.get_str()));
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
Json map(const Pin &pin) {
    auto out = Json::object(), entries = Json::array();
    const auto tensor = pin.tensor();
    const auto shape = tensor->domain().extents();
    out.set("rows", cardinality(shape.at(0)));
    out.set("columns", cardinality(shape.at(1)));
    out.set("complete_support", Json::boolean(true));
    for (const auto &e : tensor->values().entries()) {
        auto item = Json::object();
        item.set("row", cardinality(e.index / shape[1]));
        item.set("column", cardinality(e.index % shape[1]));
        item.set("value", Json::string(e.value.get_str()));
        entries.push(item);
    }
    out.set("entries", entries);
    return out;
}
void complete_map(const Pin &pin) {
    require(pin.kind() == storage::NativeAccessionKind::Operator &&
                pin.tensor()->domain().extents().size() == 2 &&
                pin.tensor()->support().count() == pin.tensor()->domain().logical_size(),
            "linear observation requires a complete declared native map");
}
std::vector<std::size_t> coordinates(const Json &value, std::size_t n, bool permutation) {
    source_array(value, "source coordinates");
    require(!permutation || value.size() == n, "ordering must cover the complete source");
    std::vector<std::size_t> result;
    std::set<std::size_t> seen;
    for (std::size_t i = 0; i < value.size(); ++i) {
        const auto index = count(value.at(i), n ? n - 1 : 0, "original source coordinate");
        require(index < n, "coordinate outside source");
        require(!permutation || seen.insert(index).second, "ordering must be a bijection");
        result.push_back(index);
    }
    return result;
}
std::size_t validate_request(const Json &query, const NativeObservationBinding &binding) {
    members(query,
            {"operation", "publication", "source_sequence", "operator", "accession", "fields",
             "ordering", "observation"},
            {"workers"});
    require(query.at("operation").text() == "linear_response_observation",
            "wrong observation operation");
    for (const auto *name : {"publication", "operator", "accession"})
        (void)text(query.at(name));
    require(query.at("source_sequence").text() == std::to_string(binding.operation.sequence()),
            "observation source_sequence differs from actual publication");
    if (query.has("workers"))
        require(count(query.at("workers"), std::numeric_limits<std::size_t>::max(), "workers") > 0,
                "observation workers must be positive");
    complete_map(binding.operation);
    complete_map(binding.accession);
    const auto origin = binding.operation.carrier(), target = binding.accession.carrier();
    const auto n = origin->domain.logical_size(), q = target->domain.logical_size();
    require(binding.operation.same_execution_authority(binding.accession) &&
                binding.operation.endpoint().get() == binding.accession.endpoint().get(),
            "observation requires maps from the actual same publication");
    const auto op = binding.operation.tensor_declaration(),
               access = binding.accession.tensor_declaration();
    require(
        origin->kind == storage::TemporalCarrierKind::SourceGrade &&
            origin->domain == numeric::IndexDomain::linear(n) && op->input == op->output &&
            op->tensor.domain().extents() == std::vector<std::size_t>{n, n},
        "observation operator must be a complete endomorphism on its actual SourceGrade carrier");
    require(access->role == storage::TemporalTensorRole::Accession &&
                access->input == origin->key && target->domain == numeric::IndexDomain::linear(q) &&
                access->tensor.domain().extents() == std::vector<std::size_t>{q, n},
            "observation requires an actual accession with declared input and output carriers");
    source_array(query.at("fields"), "observation fields");
    require(!binding.fields.empty() && binding.fields.size() == query.at("fields").size(),
            "observation requires nonempty ordered fields");
    for (const auto &field : binding.fields) {
        const auto t = field.tensor_declaration();
        require(field.same_execution_authority(binding.operation) &&
                    field.endpoint().get() == binding.operation.endpoint().get() &&
                    field.kind() == storage::NativeAccessionKind::Field &&
                    t->role == storage::TemporalTensorRole::Field && t->output == origin->key &&
                    t->tensor.domain() == numeric::IndexDomain::linear(n) &&
                    t->tensor.support().count() == n,
                "observation requires complete fields on the actual operator carrier");
    }
    const auto order = query.at("ordering");
    const auto kind = order.at("kind").text();
    if (kind == "declared") {
        members(order, {"kind", "rows", "columns"});
        (void)coordinates(order.at("rows"), n, true);
        (void)coordinates(order.at("columns"), n, true);
    } else {
        members(order, {"kind"});
        require(kind == "source" || kind == "minimum_degree", "unknown execution ordering");
    }
    const auto observation = query.at("observation");
    if (observation.at("kind").text() == "coordinates") {
        members(observation, {"kind", "coordinates"});
        const auto selected = coordinates(observation.at("coordinates"), n, false);
        require(selected.size() == q, "selector must cover every declared output row");
        // Validate the entire declared selector, including structural zeros.
        require(binding.accession.tensor()->values().nonzeros() == q,
                "accession differs from complete ordered selector");
        std::size_t row = 0;
        for (const auto &entry : binding.accession.tensor()->values().entries()) {
            require(entry.index == row * n + selected[row] && entry.value == 1,
                    "accession differs from complete ordered selector");
            ++row;
        }
    } else {
        members(observation, {"kind"});
        require(observation.at("kind").text() == "linear", "unknown observation kind");
    }
    return std::max<std::size_t>({1, n, q});
}
numeric::ExactActionGuard guard(const Pin &pin, const runtime::ExecutionLease &lease,
                                std::size_t dimension) {
    const auto &host = lease.policy().limits;
    return numeric::ExactActionGuard({dimension, std::min(host.rational_slots, host.sparse_entries),
                                      std::min(host.scalar_work, pin.limits().max_action_work),
                                      host.integer_bits},
                                     [&lease] { lease.check_cancelled(); });
}
template <class Report> Json receipt(const Report &report) {
    auto out = Json::object(), limits = Json::object();
    out.set("profile", Json::string(report.profile));
    out.set("source_namespace", Json::string(std::to_string(report.source_namespace)));
    out.set("source_sequence", Json::string(std::to_string(report.sequence)));
    out.set("operator_key", Json::string(std::to_string(report.operation)));
    out.set("accession_key", Json::string(std::to_string(report.accession)));
    out.set("input_field_key", Json::string(std::to_string(report.input_field)));
    out.set("input_carrier", Json::string(std::to_string(report.input_carrier)));
    out.set("output_carrier", Json::string(std::to_string(report.output_carrier)));
    out.set("source_grade", cardinality(report.source_grade.value));
    out.set("base_preparation_policy", Json::string(report.base_preparation_policy));
    out.set("preparation_policy", Json::string(report.preparation_policy));
    out.set("execution_policy", Json::string(report.execution_policy));
    out.set("preparation_source_sequence",
            Json::string(std::to_string(report.preparation_availability.max_sequence)));
    out.set("execution_source_sequence",
            Json::string(std::to_string(report.execution_availability.max_sequence)));
    out.set("preparation_steps", cardinality(report.preparation_steps));
    out.set("execution_steps", cardinality(report.execution_steps));
    out.set("uses_full_solve", Json::boolean(report.uses_full_solve));
    out.set("fallback", Json::boolean(report.fallback));
    limits.set("dimension", Json::string(std::to_string(report.action_limits.dimension)));
    limits.set("entries", Json::string(std::to_string(report.action_limits.entries)));
    limits.set("steps", Json::string(std::to_string(report.action_limits.steps)));
    limits.set("integer_bits", Json::string(std::to_string(report.action_limits.integer_bits)));
    out.set("action_limits", limits);
    out.set("requirements", requirements(report.selected_requirements));
    return out;
}
} // namespace
void NativeObservationPlans::validate(const Json &query, const NativeObservationBinding &binding) {
    (void)validate_request(query, binding);
}
Json NativeObservationPlans::execute(const Json &query, const NativeObservationBinding &binding,
                                     runtime::ExecutionLease lease, CpuResources &cpu) {
    lease.check_cancelled();
    const auto dimension = validate_request(query, binding);
    const runtime::WorkAvailability available{binding.operation.source_namespace(),
                                              binding.operation.sequence()};
    const auto capacity = cpu.coordinator().capacity();
    const auto workers = query.has("workers") ? count(query.at("workers"), capacity, "workers")
                                              : std::min(capacity, binding.fields.size());
    auto reservation = cpu.reserve(workers);
    auto descriptor = Json::object();
    for (const auto *name : {"publication", "source_sequence", "operator", "ordering"})
        descriptor.set(name, query.at(name));
    auto declared_order = Json::object();
    declared_order.set("kind", query.at("ordering").at("kind"));
    if (query.at("ordering").at("kind").text() == "declared") {
        declared_order.set("rows", query.at("ordering").at("rows"));
        declared_order.set("columns", query.at("ordering").at("columns"));
    }
    descriptor.set("ordering", declared_order);
    const auto base_key = descriptor.dump_compact();
    auto existing_base = bases_.find(base_key);
    std::optional<Base> prepared;
    auto ordering = Json::parse(query.at("ordering").dump());
    if (existing_base == bases_.end()) {
        const auto kind = ordering.at("kind").text();
        if (kind == "source")
            prepared.emplace(inference::PreparedNativeWoodburySolve::prepare_base(
                binding.operation, lease, available, dimension));
        else {
            numeric::ExactSolveOrdering order;
            const auto n = binding.operation.carrier()->domain.logical_size();
            if (kind == "declared") {
                order.rows = coordinates(ordering.at("rows"), n, true);
                order.columns = coordinates(ordering.at("columns"), n, true);
            } else {
                auto action = guard(binding.operation, lease, dimension);
                action.entries(binding.operation.tensor()->values().nonzeros());
                numeric::ExactMatrix matrix(n, n);
                for (const auto &e : binding.operation.tensor()->values().entries()) {
                    action.step();
                    action.value(e.value);
                    matrix.set(e.index / n, e.index % n, e.value);
                }
                order = numeric::minimum_degree_ordering(std::move(matrix), action);
                ordering.set("rows", indices(order.rows));
                ordering.set("columns", indices(order.columns));
            }
            prepared.emplace(inference::PreparedNativeWoodburySolve::prepare_base_ordered(
                binding.operation, std::move(order), lease, available, dimension));
        }
    } else {
        prepared.emplace(existing_base->second.owner);
        ordering = existing_base->second.ordering;
    }
    descriptor.set("accession", query.at("accession"));
    auto declared_observation = Json::object();
    declared_observation.set("kind", query.at("observation").at("kind"));
    if (query.at("observation").has("coordinates"))
        declared_observation.set("coordinates", query.at("observation").at("coordinates"));
    descriptor.set("observation", declared_observation);
    const auto plan_key = descriptor.dump_compact();
    auto existing_plan = plans_.find(plan_key);
    std::optional<Plan> plan;
    if (existing_plan != plans_.end())
        plan.emplace(existing_plan->second);
    else if (query.at("observation").at("kind").text() == "linear")
        plan.emplace(inference::PreparedNativeLinearObservation::prepare(
            *prepared, {binding.accession, dimension}, available, lease));
    else
        plan.emplace(inference::PreparedNativeSelectedObservation::prepare(
            *prepared,
            {binding.accession,
             coordinates(query.at("observation").at("coordinates"),
                         binding.operation.carrier()->domain.logical_size(), false),
             dimension},
            available, lease));
    // Workers handle immutable native owners and private result slots only.
    // JSON metadata is assembled after joining, avoiding shared codec refcounts.
    using Derived = std::variant<inference::DerivedNativeLinearObservation,
                                 inference::DerivedNativeSelectedObservation>;
    std::vector<std::optional<Derived>> results(binding.fields.size());
    runtime::CpuTileExecutor executor({workers, capacity, binding.fields.size()});
    (void)executor.run_tiles(
        binding.fields.size(), [&](std::size_t i, std::size_t, const runtime::CpuTileStop &stop) {
            stop.check();
            std::visit(
                [&](const auto &owner) {
                    results[i].emplace(owner.execute(binding.fields[i], available, lease));
                },
                *plan);
            stop.check();
        });
    auto out = Json::object(), fields = Json::array(), preparation = Json::object(),
         work = Json::object();
    for (std::size_t i = 0; i < results.size(); ++i) {
        auto item = Json::object();
        item.set("id", query.at("fields").at(i));
        item.set("meaning", binding.declarations.at("fields").at(i));
        numeric::RationalVector input(binding.fields[i].tensor()->domain().logical_size());
        for (const auto &e : binding.fields[i].tensor()->values().entries())
            input[e.index] = e.value;
        item.set("input", vector(input));
        std::visit(
            [&](const auto &derived) {
                item.set("values", vector(derived.values()));
                item.set("native_execution", receipt(derived.provenance()));
            },
            *results[i]);
        fields.push(item);
    }
    std::visit(
        [&](const auto &owner) {
            preparation.set("steps", cardinality(owner.preparation_steps()));
            preparation.set("retained_requirements", requirements(owner.retained_requirements()));
            preparation.set("uses_full_solve",
                            Json::boolean(owner.numeric_plan().uses_full_solve()));
            using T = std::decay_t<decltype(owner)>;
            if constexpr (std::is_same_v<T, inference::PreparedNativeLinearObservation>) {
                preparation.set("row_reuse_word", indices(owner.numeric_plan().row_reuse_word()));
                preparation.set("distinct_coefficient_rows",
                                cardinality(owner.numeric_plan().distinct_rows()));
                preparation.set("certified_nonzero_rows",
                                cardinality(owner.numeric_plan().certified_nonzero_rows()));
            } else {
                auto cells = Json::array();
                for (const auto &cell : owner.source_cells()) {
                    auto c = Json::object();
                    c.set("grade", cardinality(cell.grade.value));
                    c.set("key", Json::string(std::to_string(cell.key)));
                    cells.push(c);
                }
                out.set("selected_source_cells", cells);
            }
        },
        *plan);
    preparation.set("dimension_admission", cardinality(dimension));
    work.set("reserved_workers", cardinality(workers));
    work.set("completed_fields", cardinality(results.size()));
    work.set("executor", Json::string("rfc::runtime::CpuTileExecutor"));
    out.set("publication", binding.declarations);
    out.set("operator_map", map(binding.operation));
    out.set("accession_map", map(binding.accession));
    out.set("ordering", ordering);
    out.set("observation", query.at("observation"));
    out.set("fields", fields);
    out.set("native_preparation", preparation);
    out.set("native_work", work);
    // A rejected action must not publish a new cached preparation or counters.
    auto detached = Json::parse(out.dump());
    lease.check_cancelled();
    const auto next = machine::checked_size_add(applications_, binding.fields.size());
    if (existing_base == bases_.end())
        bases_.emplace(base_key, BaseRecord{*prepared, Json::parse(ordering.dump())});
    if (existing_plan == plans_.end())
        plans_.emplace(plan_key, *plan);
    applications_ = next;
    return detached;
}
Json NativeObservationPlans::statistics() const {
    auto out = Json::object();
    out.set("base_solve_preparations", cardinality(bases_.size()));
    out.set("observation_preparations", cardinality(plans_.size()));
    out.set("field_applications", cardinality(applications_));
    return out;
}
} // namespace gadget::scientific_detail
