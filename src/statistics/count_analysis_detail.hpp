#pragma once
#include "../scientific_detail.hpp"
#include "../scientific_result_fields.hpp"
#include <gadget/record_view.hpp>
#include <gadget/resources.hpp>
#include <gadget/statistics/count_analysis.hpp>
#include <gadget/statistics/counts.hpp>
#include <mutex>
#include <rfc/database/native_tensor_window.hpp>

namespace gadget::statistics::count_detail {
struct Sample {
    std::string id, unit, native_id;
    std::size_t coordinate;
};
struct Feature {
    std::string id, field;
};
struct Configuration {
    std::vector<Sample> samples;
    std::vector<Feature> features;
    std::vector<std::string> columns, normalization_fields, pool_fields;
    std::vector<rfc::numeric::RationalVector> rows;
    std::vector<double> contrast, factors, exposure;
    std::vector<std::size_t> reduced;
    std::string quantity, units, family, dispersion_method, normalization, inference, correction;
    bool center{}, fitted{};
    double alpha{}, confidence{};
    CountFitControl control;
    DispersionControl profile;
    std::optional<LogDispersionPrior> effect_prior;
    std::size_t workers{1}, work{}, bytes{}, request_bytes{};
};
Configuration admit_count_request(const NativeRecord &request);
double exact_count(const rfc::numeric::Rational &value);
struct Observation {
    std::vector<double> values;
    std::size_t missing{}, invalid{};
};
struct FeatureResult {
    CountFit fit;
    std::optional<CountFit> reduced_fit;
    std::optional<DispersionFit> profile, reduced_profile;
    std::optional<NormalEffectPosterior> posterior;
    std::string status;
    std::size_t missing{}, invalid{};
    double alpha{};
    std::optional<double> log_p, adjusted_log_p, statistic, z, low, high;
};
} // namespace gadget::statistics::count_detail
