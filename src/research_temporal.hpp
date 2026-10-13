#pragma once
#include "scientific_detail.hpp"
#include "scientific_result_fields.hpp"
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
struct TemporalResultFields {
    scientific_detail::ResultFields &publication;
    std::vector<std::vector<std::uint64_t>> axes;
    std::vector<scientific_detail::ResultValue> declarations, coordinates;
    rfc::runtime::ExecutionLease lease;
};
Json temporal_field_history(const Json &query, const std::vector<TemporalReading> &readings,
                            const std::vector<TemporalMaps> &maps, bool validate);
scientific_detail::ResultValue temporal_history_fields(const Json &query,
    const std::vector<TemporalReading> &readings, const std::vector<TemporalMaps> &maps,
    TemporalResultFields &fields);
struct CensoredResultFields {
    scientific_detail::ResultFields &publication;
    std::uint64_t source, groups, contrasts;
};
scientific_detail::ResultValue censored_endpoint_bounds(const NativeStudy &study, const Json &query,
                              const std::map<std::string, Vector> &fields,
                              const std::vector<bool> &finite_upper_support,
                              const std::map<std::string, Json> &meanings, const Json &coordinates,
                              bool validate, CensoredResultFields *native = nullptr);
} // namespace gadget::research_detail
