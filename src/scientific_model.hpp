#pragma once
#include "scientific_projection.hpp"
#include <stdexcept>

namespace gadget::scientific_detail {
// One canonical projection owner for original coordinate observed models.
// SDK actions supply every rank, kernel, Gram preparation and field application.
struct PreparedObservedModel {
    rfc::numeric::ExactMatrix design, kernel;
    std::vector<std::size_t> pivots;
    PreparedProjection projection;
    PreparedObservedModel(rfc::numeric::ExactMatrix map, const rfc::field::ExactMetric &metric)
        : design(std::move(map)), kernel(design.kernel_frame()),
          projection(design.image_frame(), metric) {
        (void)design.rref(pivots);
        if (pivots.size() != projection.frame.cols() || !design.multiply(kernel).entries().empty())
            throw std::runtime_error("observed model frame/kernel certificate");
    }
    struct Fit {
        rfc::numeric::RationalVector coefficients, represented, residual;
        rfc::numeric::Rational quadrance;
    };
    Fit apply(const rfc::numeric::RationalVector &input,
              std::function<void()> checkpoint = {}) const {
        using namespace rfc::numeric;
        const auto frame_coefficients = projection.coefficients(input, std::move(checkpoint));
        Fit out;
        out.coefficients.assign(design.cols(), Rational{0});
        for (std::size_t j = 0; j < pivots.size(); ++j)
            out.coefficients[pivots[j]] = frame_coefficients[j];
        out.represented = projection.frame.apply(frame_coefficients);
        out.residual = sub(input, out.represented);
        out.quadrance = projection.metric.inner(out.residual, out.residual);
        if (!is_zero(design.transpose().apply(projection.metric.apply(out.residual))) ||
            design.apply(out.coefficients) != out.represented ||
            add(out.represented, out.residual) != input ||
            projection.metric.inner(out.represented, out.residual) != 0)
            throw std::runtime_error("observed model reconstruction/normal equations");
        return out;
    }
};
} // namespace gadget::scientific_detail
