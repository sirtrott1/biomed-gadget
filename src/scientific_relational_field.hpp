#pragma once
#include "scientific_detail.hpp"
#include "scientific_projection.hpp"
#include <array>
#include <gadget/study.hpp>
#include <rfc/core/axioms.hpp>

namespace gadget::scientific_detail {
// Source bound orthogonal actions; native SDK owns frames, solves and arithmetic.
struct PreparedRelationalField {
    const NativeStudy &study;
    rfc::core::Grade grade;
    const rfc::field::ExactMetric &metric;
    rfc::core::BianchiReport bianchi;
    std::array<PreparedProjection, 2> projections;
    static rfc::core::BianchiReport chain(const NativeStudy &owner) {
        auto report = rfc::core::bianchi_through(owner.source(), owner.source().max_grade());
        require(report.satisfied(), "complete chain family required for Hodge decomposition");
        return report;
    }
    PreparedRelationalField(const NativeStudy &owner, rfc::core::Grade k,
                           std::optional<rfc::numeric::ExactActionLimits> limits = {},
                           std::function<void()> checkpoint = {})
        : study(owner), grade(k), metric(owner.tower().metrics().at(k.value)),
          bianchi(chain(owner)),
          projections{PreparedProjection(owner.hodge(k).gradient_frame(), metric, limits, checkpoint),
                      PreparedProjection(owner.hodge(k).curl_frame(), metric, limits, checkpoint)} {
        require(rfc::machine::checked_size_add(projections[0].frame.cols(),
                                               projections[1].frame.cols()) <= metric.size(),
                "Hodge complementary dimension");
    }
    std::array<Vector, 4> apply(const Vector &x, std::function<void()> checkpoint = {}) const {
        using namespace rfc;
        std::array<Vector, 4> fields{
            x, projections[0].apply(x, checkpoint), projections[1].apply(x, checkpoint), {}};
        fields[3] = numeric::sub(x, numeric::add(fields[1], fields[2]));
        if (grade.value)
            require(numeric::is_zero(study.source().boundary_apply(grade, fields[3])),
                    "harmonic lower boundary");
        if (grade.value < study.source().max_grade().value) {
            const auto upper = core::Grade{static_cast<std::uint16_t>(grade.value + 1u)};
            const field::BoundaryOperator boundary(study.source(), upper);
            require(numeric::is_zero(
                        boundary.adjoint(fields[3], metric, study.tower().metrics()[upper.value])),
                    "harmonic upper adjoint");
        }
        require(numeric::add(numeric::add(fields[1], fields[2]), fields[3]) == x,
                "full field reconstruction");
        Rational quadrance{0};
        for (std::size_t a = 1; a < fields.size(); ++a) {
            quadrance += metric.inner(fields[a], fields[a]);
            for (std::size_t b = a + 1; b < fields.size(); ++b)
                require(metric.inner(fields[a], fields[b]) == 0, "sector orthogonality");
        }
        require(quadrance == metric.inner(x, x), "sector quadrance budget");
        return fields;
    }
};
} // namespace gadget::scientific_detail
