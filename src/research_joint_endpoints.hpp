#pragma once
#include "scientific_detail.hpp"
#include <gadget/study.hpp>

namespace gadget::research_detail {
Json joint_endpoint_family(const NativeStudy &study, const Json &declaration, const Vector &lower,
                           const Vector &upper, const Vector &lower_closed,
                           const Vector &upper_closed, const Vector &upper_finite,
                           const Vector &selected, const std::vector<Rational> &times,
                           const Rational &origin, const rfc::numeric::ExactMatrix &readout,
                           const rfc::numeric::ExactMatrix &contrasts, const Json &coordinates,
                           bool validate);
}
