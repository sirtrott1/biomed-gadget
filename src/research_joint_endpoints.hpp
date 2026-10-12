#pragma once
#include "research_temporal.hpp"
#include <gadget/study.hpp>

namespace gadget::research_detail {
scientific_detail::ResultValue joint_endpoint_family(const NativeStudy &study, const Json &declaration, const Vector &lower,
                           const Vector &upper, const Vector &lower_closed,
                           const Vector &upper_closed, const Vector &upper_finite,
                           const Vector &selected, const std::vector<Rational> &times,
                           const Rational &origin, const rfc::numeric::ExactMatrix &readout,
                           const rfc::numeric::ExactMatrix &contrasts, const Json &coordinates,
                           bool validate, CensoredResultFields *native);
}
