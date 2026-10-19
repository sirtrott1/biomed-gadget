#pragma once
#include <gadget/io.hpp>
namespace gadget::statistics {
Json count_answer_view(const fs::path &directory, const Json &reference);
Json count_answer_reference(const fs::path &directory, const Json &result, const Json &request);
} // namespace gadget::statistics
