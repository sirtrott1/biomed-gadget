#pragma once
#include <gadget/io.hpp>
#include <optional>

namespace gadget::pipeline_detail {
struct ResultOperation {
    std::string operation, analysis_binding;
    fs::path source_path;
    std::string source_sha256;
    std::size_t source_bytes{}, analysis_bytes{};
};
ResultOperation admit_result_operation(const Json &request);
void admit_result_step(const Json &node, const std::map<std::string, std::string> &declared);
// Read once under the declared byte bound and hash exactly the decoded image.
std::string verified_record_image(const fs::path &, const std::string &sha256, std::size_t limit);
} // namespace gadget::pipeline_detail
