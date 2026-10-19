#include <gadget/build.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <gadget/statistics/execution.hpp>
#include <mutex>
#include <numeric>
#include <rfc/field/augmented_green.hpp>
#include <rfc/statistics/field_moments.hpp>
#include <set>
#include <stdexcept>

namespace gadget::statistics {
namespace {
void members(const Json &query, const std::set<std::string> &allowed) {
    if (!json_object_is_type(query.raw(), json_type_object))
        throw std::invalid_argument("query must be an object");
    json_object_object_foreach(query.raw(), key, value) {
        (void)value;
        if (!allowed.contains(key))
            throw std::invalid_argument(std::string("unsupported query member: ") + key);
    }
}
std::size_t bounded_integer(const Json &j, const std::string &key, std::size_t limit,
                            std::size_t fallback) {
    if (!j.has(key))
        return fallback;
    auto n = j.at(key).integer();
    if (n < 0 || static_cast<std::uint64_t>(n) > limit)
        throw std::invalid_argument("query budget/range: " + key);
    return n;
}
Vector input_field(const Workspace &w, const Execution &op, const Json &query, std::uint16_t grade,
                   const std::string &key = "field") {
    if (query.has("values")) {
        if (query.has(key))
            throw std::invalid_argument("choose a named field or explicit values");
        auto values = query.at("values");
        auto size = op.complex().grade(rfc::core::Grade{grade}).cell_count();
        if (values.size() != size)
            throw std::invalid_argument("values carrier dimension mismatch");
        Vector out;
        out.reserve(size);
        for (std::size_t i = 0; i < size; ++i)
            out.push_back(exact_decimal(values.at(i).text()));
        return out;
    }
    if (grade != 0)
        throw std::invalid_argument("named GET field has grade zero");
    auto name = query.at(key).text();
    if (name.starts_with("UNION."))
        return union_field(w, name.substr(6));
    if (!w.fields.contains(name) || !name.starts_with(op.context + "."))
        throw std::invalid_argument("field/context mismatch: " + name);
    return w.fields.at(name);
}
Json ranked(const Workspace &w, const Execution &op, const Vector &value, std::uint16_t grade,
            std::size_t limit) {
    std::vector<std::size_t> order(value.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](auto a, auto b) {
        Rational av = abs(value[a]), bv = abs(value[b]);
        return av == bv ? a < b : av > bv;
    });
    auto result = Json::array();
    for (std::size_t rank = 0; rank < std::min(limit, order.size()); ++rank) {
        auto i = order[rank];
        auto row = Json::object();
        row.set("rank", Json::integer(rank + 1));
        row.set("coordinate", Json::integer(i));
        row.set("value", exact_number(value[i]));
        if (grade == 0) {
            row.set("identity", Json::string(w.genes[i]));
            row.set("provenance", w.provenance(static_cast<std::uint32_t>(i)));
        } else {
            auto id = op.source_relations.at(i);
            row.set("identity", Json::string(w.assertions[id].id));
            row.set("source_relation_key", Json::integer(id + 1));
        }
        result.push(row);
    }
    return result;
}
} // namespace
using ExecutionProvider = std::function<Execution &(const std::string &, const std::string &)>;
static Json execute_query_prepared(const Workspace &w, const Json &query,
                                   const ExecutionProvider &execution,
                                   const std::string &manifest_hash) {
    auto operation = query.at("operation").text();
    const std::map<std::string, std::string> definitions = {
        {"inspect", "retained source manifest and independent field semantics"},
        {"provenance",
         "retained assertions, occurrences and original measurement rows incident to a named gene"},
        {"boundary", "B_k x on declared source and target carriers"},
        {"adjoint", "M_1^-1 B_1^T M_0 x"},
        {"hodge", "L_k x = (B_k^dagger B_k + B_(k+1) B_(k+1)^dagger) x"},
        {"typed_moment", "<x,L_k^n y>_M on the declared context and metric"},
        {"green", "(I+lambda L_k)^-1 j with declared arithmetic and residual check"},
        {"context_diff",
         "||(L_other-L_context)e_gene||^2 on the identical ambient gene identities"}};
    if (!definitions.contains(operation))
        throw std::invalid_argument("unsupported operation: " + operation);
    auto certificate = Json::object();
    certificate.set("schema", Json::string("gadget.certificate.v1"));
    certificate.set("source_digest", Json::string(w.source_digest));
    certificate.set("source_root", Json::string(w.root.string()));
    if (w.snapshot_origin.has("study_digest"))
        certificate.set("snapshot", w.snapshot_origin);
    if (w.store_origin.has("root"))
        certificate.set("store", w.store_origin);
    certificate.set("rexgraph_headers_sha256", Json::string(gadget::build::rexgraph_sha256()));
    certificate.set("gadget_sources_sha256", Json::string(gadget::build::source_sha256()));
    certificate.set("query", query);
    certificate.set("mathematical_definition", Json::string(definitions.at(operation)));
    certificate.set("source_manifest_sha256", Json::string(manifest_hash));
    auto meanings = Json::object();
    for (const auto &key : {"field", "other_field"})
        if (query.has(key)) {
            auto name = query.at(key).text();
            if (w.field_semantics.contains(name))
                meanings.set(name, w.field_semantics.at(name));
        }
    certificate.set("field_semantics", meanings);
    certificate.set("exact", Json::boolean(true));
    if (operation == "inspect") {
        members(query, {"operation"});
        certificate.set("result", w.manifest());
    } else if (operation == "provenance") {
        members(query, {"operation", "gene"});
        auto name = query.at("gene").text();
        if (!w.gene_index.contains(name))
            throw std::invalid_argument("unknown gene: " + name);
        certificate.set("result", w.provenance(w.gene_index.at(name)));
    } else {
        const std::set<std::string> operations = {"boundary",     "adjoint", "hodge",
                                                  "typed_moment", "green",   "context_diff"};
        if (!operations.contains(operation))
            throw std::invalid_argument("unsupported operation: " + operation);
        std::set<std::string> allowed = {"operation", "context", "realization", "grade"};
        if (operation == "context_diff")
            allowed.insert({"gene", "other_context"});
        else if (operation == "typed_moment")
            allowed.insert({"field", "other_field", "order"});
        else
            allowed.insert({"field", "values", "top_k"});
        if (operation == "green")
            allowed.insert({"lambda", "tolerance", "mode"});
        members(query, allowed);
        auto context = query.at("context").text(), realization = query.at("realization").text();
        if (!w.contexts.contains(context))
            throw std::invalid_argument("unknown context");
        if (realization != "support" && realization != "rp")
            throw std::invalid_argument("unknown realization");
        auto grade = static_cast<std::uint16_t>(bounded_integer(query, "grade", 1, 0));
        auto parameters = Json::object();
        parameters.set("grade", Json::integer(grade));
        if (operation == "typed_moment")
            parameters.set("order", Json::integer(bounded_integer(query, "order", 16, 1)));
        certificate.set("effective_parameters", parameters);
        if (operation == "boundary" && grade != 1)
            throw std::invalid_argument("boundary source grade must be one");
        if (operation == "adjoint" && grade != 0)
            throw std::invalid_argument("adjoint source grade must be zero");
        if (operation == "context_diff" && grade != 0)
            throw std::invalid_argument(
                "context differential currently requires common gene carrier");
        auto &op = execution(context, realization);
        certificate.set(
            "metric", Json::string("M0=identity; M1=boundary quadrance; zero-column convention=1"));
        certificate.set("carrier",
                        Json::string(grade == 0 ? "common ambient gene identities"
                                                : "explicit selected source assertion keys"));
        if (operation == "context_diff") {
            auto other = query.at("other_context").text();
            if (other == context || !w.contexts.contains(other))
                throw std::invalid_argument("comparison requires distinct named contexts");
            auto name = query.at("gene").text();
            if (!w.gene_index.contains(name))
                throw std::invalid_argument("unknown gene");
            auto gene = w.gene_index.at(name);
            auto &target = execution(other, realization);
            Vector basis(w.genes.size(), Rational{0});
            basis[gene] = 1;
            auto delta = rfc::numeric::sub(target.hodge().apply(basis), op.hodge().apply(basis));
            auto result = Json::object();
            result.set("delta_L_squared", exact_number(rfc::numeric::dot(delta, delta)));
            result.set("delta_L_norm", Json::number(norm(delta)));
            result.set("gene", Json::string(name));
            result.set("provenance", w.provenance(gene));
            certificate.set("result", result);
        } else if (operation == "typed_moment") {
            auto a = input_field(w, op, query, grade);
            auto b = input_field(w, op, query, grade, "other_field");
            auto order = bounded_integer(query, "order", 16, 1);
            rfc::field::ExactHodgeOperator hodge(op.complex(), rfc::core::Grade{grade},
                                                 op.metrics());
            auto powered = rfc::statistics::hodge_power_apply(hodge, order, std::move(b));
            certificate.set("result", exact_number(hodge.metric().inner(a, powered)));
        } else {
            auto field = input_field(w, op, query, grade);
            Vector value;
            std::uint16_t output_grade = grade;
            auto limit = bounded_integer(query, "top_k", 100, 10);
            if (operation == "boundary") {
                value = op.boundary().apply(field);
                output_grade = 0;
            } else if (operation == "adjoint") {
                value = op.boundary().adjoint(field, op.metrics()[0], op.metrics()[1]);
                output_grade = 1;
            } else {
                rfc::field::ExactHodgeOperator hodge(op.complex(), rfc::core::Grade{grade},
                                                     op.metrics());
                if (operation == "hodge")
                    value = hodge.apply(field);
                else {
                    auto lambda = exact_decimal(query.at("lambda").text());
                    if (lambda <= 0)
                        throw std::invalid_argument("Green lambda must be positive");
                    auto mode = query.has("mode") ? query.at("mode").text() : "exact";
                    if (mode == "exact") {
                        auto result = rfc::field::exact_green_apply(hodge, lambda, field);
                        if (!result.exact)
                            throw std::runtime_error("canonical exact Green did not certify");
                        value = std::move(result.value);
                        certificate.set("green_execution",
                                        Json::string("rfc::field::exact_green_apply; factored "
                                                     "boundary and metric actions"));
                        certificate.set("green_iterations", Json::integer(result.iterations));
                    } else if (mode == "numerical") {
                        if (grade != 0)
                            throw std::invalid_argument(
                                "relation-space numerical Green currently requires grade zero");
                        auto tolerance = query.has("tolerance")
                                             ? exact_decimal(query.at("tolerance").text()).get_d()
                                             : 1e-9;
                        Json diagnostic;
                        value = op.green(field, lambda, tolerance, diagnostic);
                        certificate.set("exact", Json::boolean(false));
                        certificate.set("numerical", diagnostic);
                    } else
                        throw std::invalid_argument("unknown Green mode");
                }
            }
            certificate.set("output_grade", Json::integer(output_grade));
            certificate.set(
                "ranking_observable",
                Json::string("absolute coordinate value with stable coordinate tie-break"));
            certificate.set("result", ranked(w, op, value, output_grade, limit));
        }
    }
    certificate.set("certificate_sha256", Json::string("sha256:" + sha256(certificate.dump())));
    return certificate;
}
Json execute_query(const Workspace &source, const Json &query) {
    std::map<std::pair<std::string, std::string>, std::unique_ptr<Execution>> owners;
    auto provider = [&](const std::string &context, const std::string &realization) -> Execution & {
        auto &owner = owners[{context, realization}];
        if (!owner)
            owner = std::make_unique<Execution>(source, context, realization);
        return *owner;
    };
    return execute_query_prepared(source, query, provider,
                                  "sha256:" + sha256(source.manifest().dump()));
}
struct QuerySession::State {
    std::shared_ptr<const Workspace> source;
    std::shared_ptr<CpuResources> resources;
    std::string manifest_hash;
    mutable std::mutex mutex;
    std::map<std::pair<std::string, std::string>, std::unique_ptr<Execution>> owners;
    State(std::shared_ptr<const Workspace> study, std::shared_ptr<CpuResources> cpu)
        : source(std::move(study)), resources(cpu ? std::move(cpu) : CpuResources::local()) {
        if (!source || !source->native_snapshot)
            throw std::invalid_argument("query sessions require an owning immutable native study");
        manifest_hash = "sha256:" + sha256(source->manifest().dump());
    }
};
QuerySession::QuerySession(std::shared_ptr<const Workspace> source,
                           std::shared_ptr<CpuResources> resources)
    : state_(std::make_unique<State>(std::move(source), std::move(resources))) {}
QuerySession::~QuerySession() = default;
Json QuerySession::execute(const Json &query) {
    std::lock_guard lock(state_->mutex);
    std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> worker;
    const bool numerical_green = query.at("operation").text() == "green" &&
        query.has("mode") && query.at("mode").text() == "numerical";
    // Numerical Green owns its existing dense phase reservation, including the caller.
    if (numerical_green) (void)state_->resources->readmit();
    else worker.emplace(state_->resources->reserve(1));
    auto provider = [&](const std::string &context, const std::string &realization) -> Execution & {
        auto &owner = state_->owners[{context, realization}];
        if (!owner)
            owner = std::make_unique<Execution>(*state_->source, context, realization,
                                                 state_->resources);
        return *owner;
    };
    return execute_query_prepared(*state_->source, query, provider, state_->manifest_hash);
}
std::size_t QuerySession::preparation_count() const {
    std::lock_guard lock(state_->mutex);
    return std::count_if(state_->owners.begin(), state_->owners.end(),
                         [](const auto &item) { return !!item.second; });
}
} // namespace gadget::statistics
