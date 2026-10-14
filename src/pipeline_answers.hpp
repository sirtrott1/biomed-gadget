#pragma once
#include <gadget/pipeline.hpp>

namespace gadget::pipeline_detail {
fs::path checked_answer_path(const fs::path &directory, const Json &record);
Json checked_answer_view(const fs::path &directory,const Json &record);
Json observable_catalog();
NativeRecord answer_admission_record(const Json &intent, const Json &catalog,
                                      const Json *included = nullptr,
                                      const std::string *question = nullptr);
Json answer_admission(const Json &intent, const Json &catalog, const Json *included = nullptr,
                      const std::string *question = nullptr);
Json answer_witnesses(const Json &planning, const Json &pipeline, const Json &results,
                      const fs::path &directory);
} // namespace gadget::pipeline_detail
