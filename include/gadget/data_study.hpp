#pragma once
#include <gadget/scientific.hpp>
#include <stop_token>

namespace gadget {
// Bind selected native matrix coordinates to an explicit observation chart.
// Values and support come from the actual retained tensor; the source complex,
// coordinate declarations and original record dependencies remain archived.
NativeRecord bind_measurement_study(const NativeRecord &data, const NativeRecord &request,
                                   std::stop_token stop = {});
}
