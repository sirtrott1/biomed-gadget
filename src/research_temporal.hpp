#pragma once
#include "scientific_detail.hpp"
#include <gadget/study.hpp>

namespace gadget::scientific_detail { class PreparedTemporalStudy; }
namespace gadget::research_detail {
struct TemporalReading {
    const NativeStudy *study;
    rfc::core::Grade grade;
    Vector values;
    Json declaration, coordinates;
    std::map<std::size_t, Vector> diagonal_fields;
    bool shared_identity = true;
    scientific_detail::PreparedTemporalStudy *prepared = nullptr;
};
using TemporalMaps = std::vector<rfc::numeric::ExactMatrix>;
Json temporal_field_history(const Json &query, const std::vector<TemporalReading> &readings,
                            const std::vector<TemporalMaps> &maps, bool validate);
Json censored_endpoint_bounds(const NativeStudy &study, const Json &query,
                              const std::map<std::string, Vector> &fields,
                              const std::vector<bool> &finite_upper_support,
                              const std::map<std::string, Json> &meanings, const Json &coordinates,
                              bool validate);
} // namespace gadget::research_detail
