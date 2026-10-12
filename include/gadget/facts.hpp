#pragma once
#include <gadget/assay.hpp>

namespace gadget {
// Selected literal source facts and completed answer values. This owning view
// captures verified bytes; model comparisons cannot change the selected paths.
class ScientificFacts {
    struct State;
    std::unique_ptr<State> state_;

  public:
    ScientificFacts(const AssayWorkspace &workspace, const Json &request,
                    const fs::path &artifact_root = {});
    ScientificFacts(const AssayWorkspace &workspace, const NativeRecord &request,
                    const fs::path &artifact_root = {});
    // Select typed values from an actual completed native operation. Workflow
    // references keep their independent owners and original flat coordinates.
    ScientificFacts(const NativeRecord &operation, const NativeRecord &request);
    ~ScientificFacts();
    ScientificFacts(const ScientificFacts &) = delete;
    ScientificFacts &operator=(const ScientificFacts &) = delete;
    const NativeRecord &native_publication() const;
    Json catalog() const;
    // Complete selected values with shared owner references. Full queries and
    // original observations remain in catalog() and every checked receipt.
    Json model_view() const;
    // Checks literal comparisons and separately records omitted selections.
    // Agreement does not verify a biological interpretation of these values.
    Json check(const Json &proposal) const;
};
} // namespace gadget
