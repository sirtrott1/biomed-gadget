#pragma once
#include <gadget/statistics/counts.hpp>
#include <rfc/numeric/exact_matrix.hpp>
namespace gadget::statistics {
bool finite_count_coefficients(const rfc::numeric::ExactMatrix &design,
                               std::span<const double> counts, CountWork &work,
                               std::stop_token stop);
}
