#pragma once
#include <gadget/scientific.hpp>

namespace rfc::learning { class PreparedTransportPatternBank; }
namespace gadget {
// The canonical bank codec retains the original context stars and joint family.
// Source and evidence publications retain its declared training inputs.
NativeRecord publish_pattern_model(const rfc::learning::PreparedTransportPatternBank &bank,
    std::uint64_t bank_namespace, const NativeRecord &source, const NativeRecord &evidence,
    const NativeRecord &declaration);

class NativePatternModel {
    struct State;
    std::shared_ptr<const State> state_;
  public:
    explicit NativePatternModel(NativeRecord publication);
    const NativeRecord &native_publication() const;
    std::string digest() const;
    NativeRecord execute_record(std::shared_ptr<const ScientificStudy> source,
        const NativeRecord &request, std::shared_ptr<CpuResources> resources = {}) const;
};
Json pattern_record_view(const NativeRecord &publication);
NativeRecord replay_pattern_record(const NativeRecord &publication,
    std::shared_ptr<CpuResources> resources = {});
} // namespace gadget
