#pragma once
#include <gadget/scientific.hpp>
#include <stop_token>

namespace gadget::statistics {
// Preparation retains the immutable study and exact observation accessions.
// Fits reuse one design and normalization; original count fields stay native.
class PreparedCountAnalysis {
    struct State;
    std::shared_ptr<State> state_;
    explicit PreparedCountAnalysis(std::shared_ptr<State> state) : state_(std::move(state)) {}

  public:
    static PreparedCountAnalysis prepare(const ScientificStudy &study, const NativeRecord &request,
                                         std::shared_ptr<CpuResources> resources = {},
                                         std::stop_token stop = {});
    NativeRecord execute(std::stop_token stop = {}) const;
};
Json count_result_view(const NativeRecord &result);
NativeRecord count_model_analysis(const ScientificStudy &study, const NativeRecord &request,
                                  std::shared_ptr<CpuResources> resources = {},
                                  std::stop_token stop = {});
} // namespace gadget::statistics
