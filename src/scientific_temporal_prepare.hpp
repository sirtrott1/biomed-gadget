#pragma once
#include "scientific_relational_field.hpp"
#include <map>
#include <memory>
#include <optional>
#include <rfc/field/low_rank_update.hpp>

namespace gadget::scientific_detail {
// Preparations belong to one admitted source and metric tower. A session owns
// the study; the borrowed constructor is confined to synchronous comparisons.
class PreparedTemporalStudy {
    using Matrix = rfc::numeric::ExactMatrix;
    using Grade = rfc::core::Grade;
    std::optional<rfc::numeric::ExactActionLimits> action_limits_;
    std::function<void()> checkpoint_;
    rfc::numeric::ExactActionLimits limits() const {
        if (action_limits_) return *action_limits_;
        const auto addressable = std::numeric_limits<std::size_t>::max();
        return {addressable, addressable, addressable, addressable};
    }

  public:
    class Green {
        const NativeStudy &study_;
        Grade grade_;
        Rational lambda_;
        rfc::numeric::ExactActionLimits limits_;
        std::function<void()> checkpoint_;
        std::optional<rfc::field::PreparedWoodburySolve> woodbury_;
        std::optional<rfc::numeric::PreparedExactSolve> direct_;
        std::optional<Matrix> matrix_;
        static Matrix adjoint(const NativeStudy &study, Grade upper) {
            const rfc::field::BoundaryOperator boundary(study.source(), upper);
            const auto &metrics = study.tower().metrics();
            Matrix result(boundary.domain_size(), boundary.codomain_size());
            for (std::size_t c = 0; c < result.cols(); ++c) {
                Vector e(result.cols());
                e[c] = 1;
                const auto column = boundary.adjoint(e, metrics.at(upper.value - 1),
                                                     metrics.at(upper.value));
                for (std::size_t r = 0; r < column.size(); ++r)
                    result.set(r, c, column[r]);
            }
            return result;
        }

      public:
        Green(PreparedTemporalStudy &prepared, Grade grade, Rational lambda)
            : study_(prepared.study()), grade_(grade), lambda_(std::move(lambda)),
              limits_(prepared.limits()), checkpoint_(prepared.checkpoint_) {
            using namespace rfc;
            const auto &study = study_;
            require(lambda_ > 0, "prepared temporal Green parameter must be positive");
            lambda_.canonicalize();
            const auto n = study.tower().metrics().at(grade.value).size();
            const auto lower = grade.value
                                   ? study.source().grade(Grade{static_cast<std::uint16_t>(
                                         grade.value - 1u)}).cell_count()
                                   : 0;
            const auto upper = grade.value < study.source().max_grade().value
                                   ? study.source().grade(Grade{static_cast<std::uint16_t>(
                                         grade.value + 1u)}).cell_count()
                                   : 0;
            const auto m = machine::checked_size_add(lower, upper);
            if (m < n) {
                // L = [B^dagger C] [B; C^dagger]. Retain every adjacent
                // coordinate, including dependent, cancelled and zero columns.
                // The native Woodbury correction has m coordinates, rather
                // than n. This comparison is source derived, not a size cap.
                Matrix u(n, m), v(n, m);
                const auto append = [&](const Matrix &a, const Matrix &b, std::size_t offset) {
                    for (const auto &[p, value] : a.entries())
                        u.set(p.first, machine::checked_size_add(offset, p.second), value);
                    for (const auto &[p, value] : b.entries())
                        v.set(p.first, machine::checked_size_add(offset, p.second), value);
                };
                if (lower)
                    append(adjoint(study, grade),
                           field::materialize_boundary_matrix(study.source(), grade).transpose(), 0);
                if (upper) {
                    const auto next = Grade{static_cast<std::uint16_t>(grade.value + 1u)};
                    append(field::materialize_boundary_matrix(study.source(), next),
                           adjoint(study, next).transpose(), lower);
                }
                woodbury_.emplace(Matrix::identity(n), std::move(u),
                                  Matrix::diagonal(Vector(m, lambda_)), std::move(v), limits_, checkpoint_);
            } else {
                // No smaller adjacent coordinate system: reuse the native
                // exact factorization on the original field coordinates.
                direct_.emplace(Matrix::identity(n).add(prepared.hodge(grade).scale(lambda_)),
                                limits_, checkpoint_);
            }
        }
        Vector apply(const Vector &x) const {
            rfc::numeric::ExactActionGuard guard(limits_,checkpoint_);
            const auto input = guard.canonical_input(x);
            auto value = woodbury_ ? woodbury_->apply(input, guard) : direct_->apply(input, guard);
            require(rfc::numeric::add(value, rfc::numeric::scale(
                        lambda_, study_.hodge(grade_).apply(value))) == input,
                    "prepared temporal Green original-source residual");
            return value;
        }
        const Matrix &matrix() {
            if (!matrix_)
                matrix_ = rfc::field::materialize_action_matrix(
                    study_.tower().metrics().at(grade_.value),
                    [&](const Vector &x) { return apply(x); });
            return *matrix_;
        }
        bool uses_woodbury() const { return woodbury_.has_value(); }
        bool has_matrix() const { return matrix_.has_value(); }
        std::size_t stored_entries() const {
            return woodbury_ ? woodbury_->stored_entries() : direct_->stored_entries();
        }
        std::size_t preparation_steps() const {
            return woodbury_ ? woodbury_->preparation_steps() : direct_->preparation_steps();
        }
    };

  private:
    std::shared_ptr<const NativeStudy> owner_;
    const NativeStudy &study_;
    std::map<std::uint16_t, std::unique_ptr<PreparedRelationalField>> sectors_;
    std::map<std::uint16_t, Matrix> hodge_;
    std::map<std::pair<std::uint16_t, Rational>, std::unique_ptr<Green>> green_;
    static const NativeStudy &checked(const std::shared_ptr<const NativeStudy> &study) {
        require(static_cast<bool>(study), "temporal preparation requires an admitted study owner");
        return *study;
    }

  public:
    explicit PreparedTemporalStudy(const NativeStudy &study,
        std::optional<rfc::numeric::ExactActionLimits> limits = {}, std::function<void()> checkpoint = {})
        : action_limits_(limits), checkpoint_(std::move(checkpoint)), study_(study) {}
    explicit PreparedTemporalStudy(std::shared_ptr<const NativeStudy> study)
        : owner_(std::move(study)), study_(checked(owner_)) {}
    const NativeStudy &study() const { return study_; }
    const std::function<void()> &checkpoint() const { return checkpoint_; }
    PreparedRelationalField &sectors(Grade grade) {
        auto it = sectors_.find(grade.value);
        if (it == sectors_.end())
            it = sectors_.emplace(grade.value,
                    std::make_unique<PreparedRelationalField>(study_, grade, action_limits_,checkpoint_)).first;
        return *it->second;
    }
    const Matrix &hodge(Grade grade) {
        auto it = hodge_.find(grade.value);
        if (it == hodge_.end())
            it = hodge_.emplace(grade.value, study_.hodge(grade).materialize()).first;
        return it->second;
    }
    Green &green(Grade grade, Rational lambda) {
        lambda.canonicalize();
        const auto key = std::make_pair(grade.value, lambda);
        auto it = green_.find(key);
        if (it == green_.end())
            it = green_.emplace(key, std::make_unique<Green>(*this, grade, lambda)).first;
        return *it->second;
    }
    std::array<std::size_t, 6> statistics() const {
        std::array<std::size_t, 6> counts{sectors_.size(), hodge_.size(), green_.size(), 0, 0, 0};
        for (const auto &[key, plan] : green_) {
            (void)key;
            counts[3] += plan->uses_woodbury();
            counts[4] += plan->has_matrix();
            counts[5] = rfc::machine::checked_size_add(counts[5], plan->stored_entries());
        }
        return counts;
    }
};
} // namespace gadget::scientific_detail
