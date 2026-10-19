#pragma once
#include <gadget/record_view.hpp>
namespace gadget::pipeline_detail {
// A presentation certificate of an actual native result and retained owners;
// no mathematical result is calculated or reconstructed by this adapter.
Json result_answer_view(const fs::path &directory, const Json &reference);
Json result_answer_reference(const fs::path &directory, const Json &result,
    const Json &request, const std::string &method, const std::string &analysis_sha256 = {});
} // namespace gadget::pipeline_detail
