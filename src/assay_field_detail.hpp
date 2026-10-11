#pragma once
#include "assay_detail.hpp"

namespace gadget::assay_detail {
struct AssayFieldsSource;
// Holds the admitted native endpoint throughout coordinate inspection and use.
struct AssayField {
    std::shared_ptr<const ScientificStudy> owner;
    Vector values;
    std::vector<unsigned char> supported;
    Json meaning, ids, keys;
    std::uint16_t grade{};
};
std::vector<AssayField> load_assay_fields(std::shared_ptr<const ScientificStudy> owner,
                                          const std::vector<std::string> &fields);
AssayField load_assay_field(std::shared_ptr<const ScientificStudy> owner, const std::string &field);
Json assay_sections(const AssayWorkspace &workspace, const Json &request,
                    std::shared_ptr<CpuResources> resources,
                    std::optional<AssayFieldsSource> *retained = nullptr);
Json section_fields_view(const NativeRecord &record);
Json relational_gluing(const AssayWorkspace &workspace, const Json &request,
                       std::shared_ptr<CpuResources> resources);
} // namespace gadget::assay_detail
