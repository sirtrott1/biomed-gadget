#pragma once
#include <gadget/io.hpp>

namespace gadget::claim_detail {
Json member_path(Json value, const Json &path);
Json literal_value(const Json &value, const std::string &type);
bool compare_claim(const Json &actual, const Json &expected, const std::string &type,
                   const std::string &comparison);
void admit_condition(const Json &condition);
Json check_literal_claim(const Json &value, const Json &claim);
void verify_seal(const Json &value, const std::string &name);
} // namespace gadget::claim_detail
