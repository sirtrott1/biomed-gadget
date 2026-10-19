#pragma once
#include <gadget/biomedical.hpp>
#include <gadget/native_record.hpp>
#include <gadget/resources.hpp>
#include <stop_token>

namespace gadget {
// Thin dispatch to the native result producers. The source_record request
// pins one bounded native source archive; all calculations consume native owners.
NativeRecord read_analysis_record(const fs::path &path, const NativeRecord &request);
NativeRecord execute_result_record(const NativeRecord &analysis,
                                             const NativeRecord &request,
                                             const EvidenceStore *evidence = nullptr,
                                             std::shared_ptr<CpuResources> resources = {},
                                             std::stop_token stop = {});
} // namespace gadget
