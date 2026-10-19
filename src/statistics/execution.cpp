#include <cmath>
#include <gadget/statistics/execution.hpp>
#include <rfc/algebra/materialize.hpp>
#include <rfc/core/builder.hpp>
#include <rfc/numeric/sparse_vector.hpp>
#include <rfc/machine/checked_size.hpp>
#include <stdexcept>

namespace gadget::statistics {
namespace {
std::size_t dense_bytes(std::size_t rows, std::size_t columns) {
    return rfc::machine::checked_size_mul(rfc::machine::checked_size_mul(rows, columns), sizeof(double));
}
Vector coefficients(const rfc::numeric::SparseRationalVector &value) {
    Vector out;
    out.reserve(value.nonzeros());
    for (const auto &e : value.entries())
        out.push_back(e.value);
    return out;
}
} // namespace
Vector binary64_values(const rfc::approximation::Vector &values) {
    Vector out;
    out.reserve(values.size());
    for (double x : values) {
        if (!std::isfinite(x))
            throw std::runtime_error("nonfinite numerical solve");
        out.emplace_back(x);
    }
    return out;
}
rfc::approximation::Vector approximate_values(const Vector &values) {
    rfc::approximation::Vector out;
    out.reserve(values.size());
    for (const auto &x : values) {
        auto d = x.get_d();
        if (!std::isfinite(d))
            throw std::invalid_argument("binary64 range exceeded");
        out.push_back(d);
    }
    return out;
}
double norm(const Vector &values) {
    auto squared = rfc::numeric::dot(values, values).get_d();
    if (!std::isfinite(squared))
        throw std::runtime_error("norm exceeds binary64 output range");
    return std::sqrt(squared);
}
Execution::Execution(const Workspace &source, std::string c, std::string r,
                     std::shared_ptr<CpuResources> resources)
    : source_digest_(source.source_digest), resources_(resources ? std::move(resources) : CpuResources::local()), context(std::move(c)), realization(std::move(r)) {
    if (!source.contexts.contains(context))
        throw std::invalid_argument("unknown context: " + context);
    if (realization != "support" && realization != "rp")
        throw std::invalid_argument("unknown realization: " + realization);
    // The ambient gene carrier and explicit keys are identical for both contexts.
    // Selected assertion coordinates retain a local -> source relation map.
    if (realization == "support") {
        auto carrier = *source.contexts.at(context);
        for (std::size_t i = 0; i < source.genes.size(); ++i)
            carrier.select({rfc::core::Grade{0}, static_cast<std::uint32_t>(i)});
        auto selected = rfc::algebra::materialize_subcomplex(*source.complex, carrier);
        source_relations = std::move(selected.to_ambient[1]);
        complex_ = std::make_unique<rfc::core::RelationalComplex>(std::move(selected.complex));
    } else {
        auto signature = source.complex->signature();
        auto type = source.complex->grade(rfc::core::Grade{1}).relation_types.front();
        rfc::core::RelationalComplexBuilder builder(signature);
        for (std::size_t i = 0; i < source.genes.size(); ++i)
            builder.declare_witness_identity(rfc::core::IdentityKey{i + 1});
        for (std::size_t i = 0; i < source.assertions.size(); ++i) {
            if (!source.contexts.at(context)->contains(
                    {rfc::core::Grade{1}, static_cast<std::uint32_t>(i)}))
                continue;
            source_relations.push_back(static_cast<std::uint32_t>(i));
            std::vector<rfc::core::BoundaryOccurrenceInput> incidences;
            for (const auto &p : source.assertions[i].participants)
                incidences.emplace_back(
                    rfc::core::IdentityKey{p.gene + 1u}, true, false,
                    realization == "support" ? Rational{1} : p.rp,
                    rfc::algebra::declared_aspects_for_incidence(
                        *source.complex, rfc::core::Grade{1},
                        source.complex->grade(rfc::core::Grade{1}).incidence_range(i).first +
                            incidences.size()),
                    rfc::core::OccurrenceKey{p.occurrence});
            builder.add_grade_one(rfc::core::RelationKey{i + 1}, type, incidences);
        }
        complex_ = std::make_unique<rfc::core::RelationalComplex>(builder.freeze());
    }
    auto B = boundary();
    Vector diagonal(B.domain_size(), Rational{1});
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
        auto column = B.apply(rfc::numeric::SparseRationalVector::basis(diagonal.size(), i));
        auto values = coefficients(column);
        auto q = rfc::numeric::dot(values, values);
        if (q != 0)
            diagonal[i] = q;
    }
    metrics_.push_back(rfc::field::ExactMetric::identity(source.genes.size()));
    metrics_.push_back(rfc::field::ExactMetric::diagonal(std::move(diagonal)));
    hodge_ =
        std::make_unique<rfc::field::ExactHodgeOperator>(*complex_, rfc::core::Grade{0}, metrics_);
}
const rfc::approximation::DenseMatrix &Execution::relation_gram() {
    prepare_gram();
    return *gram_;
}
rfc::approximation::DenseMatrix cross_relation_gram(const Execution &left, const Execution &right) {
    auto a = left.boundary(), b = right.boundary();
    if (left.source_digest() != right.source_digest())
        throw std::invalid_argument(
            "cross Gram requires the same retained source and declared identity carrier");
    if (a.codomain_size() != b.codomain_size())
        throw std::invalid_argument("cross Gram carrier shape mismatch");
    (void)left.resources_->readmit();
    auto reservation = left.resources_->coordinator().reserve(1, dense_bytes(a.domain_size(), b.domain_size()));
    // Executions are built from the same declared ambient identities; no
    // dimension only identification is used by the calling biomedical recipe.
    rfc::approximation::DenseMatrix result(a.domain_size(), b.domain_size());
    for (std::size_t j = 0; j < b.domain_size(); ++j) {
        auto column = b.apply(rfc::numeric::SparseRationalVector::basis(b.domain_size(), j));
        auto pulled = a.transpose_apply(column);
        for (const auto &e : pulled.entries()) {
            const auto value = e.value.get_d();
            if (!std::isfinite(value)) throw std::invalid_argument("cross Gram binary64 range exceeded");
            result(e.index, j) = value;
        }
    }
    return result;
}
void Execution::prepare_gram() {
    (void)resources_->readmit();
    if (gram_)
        return;
    auto B = boundary();
    auto n = B.domain_size();
    auto reservation = resources_->coordinator().reserve(1, dense_bytes(n, n));
    rfc::approximation::DenseMatrix prepared(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        auto column = B.apply(rfc::numeric::SparseRationalVector::basis(n, i));
        auto gram_column = B.transpose_apply(column);
        for (const auto &e : gram_column.entries()) {
            auto value = e.value.get_d();
            if (!std::isfinite(value))
                throw std::invalid_argument("relation Gram binary64 range exceeded");
            prepared(e.index, i) = value;
        }
    }
    resources_->coordinator().resize(reservation, 0, dense_bytes(n, n));
    gram_reservation_.emplace(std::move(reservation));
    gram_.emplace(std::move(prepared));
}
Vector Execution::green(const Vector &source, const Rational &lambda, double tolerance,
                        Json &diagnostic, std::optional<GreenWitnesses> *retained) {
    if (lambda <= 0)
        throw std::invalid_argument("Green lambda must be positive");
    if (!std::isfinite(tolerance) || tolerance <= 0 || tolerance > 1e-3)
        throw std::invalid_argument("Green tolerance must lie in (0,1e-3]");
    if (source.size() != hodge().size())
        throw std::invalid_argument("Green source carrier mismatch");
    double parameter = lambda.get_d();
    if (!std::isfinite(parameter) || parameter <= 0)
        throw std::invalid_argument("lambda exceeds numerical recipe range");
    prepare_gram();
    auto B = boundary();
    const auto matrix_bytes = dense_bytes(B.domain_size(), B.domain_size());
    const auto vector_bytes = rfc::machine::checked_size_mul(B.domain_size(), sizeof(double));
    const auto peak_bytes = rfc::machine::checked_size_add(rfc::machine::checked_size_mul(matrix_bytes, 2),
        rfc::machine::checked_size_mul(vector_bytes, 3));
    auto reservation = resources_->coordinator().reserve(1, peak_bytes);
    auto system = rfc::approximation::DenseMatrix::from_exact(metrics_[1].materialize());
    for (std::size_t i = 0; i < B.domain_size(); ++i)
        for (std::size_t j = 0; j < B.domain_size(); ++j) {
            // Preserve the separate binary64 multiply and add of scale().add().
            volatile double scaled = (*gram_)(i, j) * parameter;
            system(i, j) += scaled;
            if (!std::isfinite(system(i, j))) throw std::invalid_argument("Green system binary64 range exceeded");
        }
    auto pulled = B.transpose_apply(source);
    auto solved = binary64_values(system.solve(approximate_values(pulled)));
    // All exact carrier actions and vector algebra remain canonical. The only
    // approximation occurs inside RexGraph's declared binary64 linear solve.
    auto value =
        rfc::numeric::sub(source, rfc::numeric::scale(lambda, B.apply(solved)));
    auto residual = rfc::numeric::sub(
        rfc::numeric::add(value, rfc::numeric::scale(lambda, hodge().apply(value))), source);
    double absolute = norm(residual), relative = absolute / std::max(1.0, norm(source));
    diagnostic = Json::object();
    diagnostic.set("backend",
                   Json::string("rfc.approximation.DenseMatrix.solve / relation-space Green"));
    diagnostic.set("exact", Json::boolean(false));
    diagnostic.set("relative_residual", Json::number(relative));
    diagnostic.set("absolute_residual", Json::number(absolute));
    diagnostic.set("tolerance", Json::number(tolerance));
    diagnostic.set("relation_dimension", Json::integer(B.domain_size()));
    diagnostic.set(
        "metric",
        Json::string("M0=identity; M1=exact boundary quadrance; zero-column convention=1"));
    diagnostic.set("lambda", exact_number(lambda));
    auto allocation = Json::object();
    allocation.set("retained_gram_bytes", Json::string(std::to_string(matrix_bytes)));
    allocation.set("solve_peak_bytes", Json::string(std::to_string(peak_bytes)));
    allocation.set("retained_system_bytes", Json::string(std::to_string(retained ? matrix_bytes : 0)));
    allocation.set("scope", Json::string("Cooperating caller and dense double matrices/vectors; exact source, rational vectors, GMP limbs and returned cross Gram buffers are outside this charge"));
    diagnostic.set("allocation", allocation);
    if (relative > tolerance)
        throw std::runtime_error("numerical Green residual exceeds requested tolerance: " +
                                 std::to_string(relative));
    if (retained) {
        resources_->coordinator().resize(reservation, 0, matrix_bytes);
        retained->emplace(GreenWitnesses{resources_,
            std::make_unique<rfc::runtime::CpuResourceCoordinator::Reservation>(std::move(reservation)),
            std::move(pulled), std::move(solved), std::move(residual), std::move(system), absolute, relative, tolerance});
    }
    return value;
}
Vector union_field(const Workspace &source, const std::string &channel) {
    if (channel.size() != 1 || std::string("GET").find(channel) == std::string::npos)
        throw std::invalid_argument("unknown GET channel");
    auto a = source.fields.at("AD." + channel), b = source.fields.at("HD." + channel);
    for (std::size_t i = 0; i < a.size(); ++i)
        if (b[i] != 0)
            a[i] = 1;
    return a;
}
} // namespace gadget::statistics
