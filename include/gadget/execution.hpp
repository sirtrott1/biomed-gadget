#pragma once
#include <gadget/source.hpp>
#include <gadget/resources.hpp>
#include <optional>
#include <rfc/approximation/dense.hpp>
#include <rfc/field/exact_hodge.hpp>

namespace gadget {
// Owns every object referenced by RexGraph's non-owning operator facades.
class Execution {
    std::string source_digest_;
    std::unique_ptr<rfc::core::RelationalComplex> complex_;
    std::vector<rfc::field::ExactMetric> metrics_;
    std::unique_ptr<rfc::field::ExactHodgeOperator> hodge_;
    std::shared_ptr<CpuResources> resources_;
    std::optional<rfc::runtime::CpuResourceCoordinator::Reservation> gram_reservation_;
    std::optional<rfc::approximation::DenseMatrix> gram_;
    friend rfc::approximation::DenseMatrix cross_relation_gram(const Execution &, const Execution &);
    void prepare_gram();

  public:
    struct GreenWitnesses {
        std::shared_ptr<CpuResources> resource_owner;
        std::unique_ptr<rfc::runtime::CpuResourceCoordinator::Reservation> dense_memory;
        Vector pulled, solved, residual;
        rfc::approximation::DenseMatrix system;
        double absolute{}, relative{}, tolerance{};
    };
    std::string context, realization;
    std::vector<std::uint32_t> source_relations;
    Execution(const Workspace &source, std::string context, std::string realization,
              std::shared_ptr<CpuResources> resources = {});
    Execution(const Execution &) = delete;
    Execution &operator=(const Execution &) = delete;
    const rfc::field::ExactHodgeOperator &hodge() const { return *hodge_; }
    const rfc::core::RelationalComplex &complex() const { return *complex_; }
    const std::string &source_digest() const { return source_digest_; }
    const std::vector<rfc::field::ExactMetric> &metrics() const { return metrics_; }
    const rfc::approximation::DenseMatrix &relation_gram();
    rfc::field::BoundaryOperator boundary() const { return {*complex_, rfc::core::Grade{1}}; }
    Vector green(const Vector &source, const Rational &lambda, double tolerance, Json &diagnostic,
                 std::optional<GreenWitnesses> *retained = nullptr);
};
double norm(const Vector &value);
Vector binary64_values(const rfc::approximation::Vector &value);
rfc::approximation::Vector approximate_values(const Vector &value);
rfc::approximation::DenseMatrix cross_relation_gram(const Execution &left, const Execution &right);
Vector union_field(const Workspace &source, const std::string &channel);
Json execute_query(const Workspace &source, const Json &query);
// One immutable native study; preparation never crosses publication identities.
// Calls on one session serialize. Independent sessions have private preparation.
class QuerySession {
    struct State;
    std::unique_ptr<State> state_;

  public:
    explicit QuerySession(std::shared_ptr<const Workspace> source,
                          std::shared_ptr<CpuResources> resources = {});
    ~QuerySession();
    QuerySession(const QuerySession &) = delete;
    QuerySession &operator=(const QuerySession &) = delete;
    Json execute(const Json &query);
    std::size_t preparation_count() const;
};
void analyze(const Workspace &source, const fs::path &destination);
} // namespace gadget
