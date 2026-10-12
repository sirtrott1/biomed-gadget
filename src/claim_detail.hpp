#pragma once
#include <gadget/io.hpp>
#include <gadget/native_record.hpp>

namespace gadget::claim_detail {
Json member_path(Json value, const Json &path);
Json literal_value(const Json &value, const std::string &type);
bool compare_claim(const Json &actual, const Json &expected, const std::string &type,
                   const std::string &comparison);
void admit_condition(const Json &condition);
Json check_literal_claim(const Json &value, const Json &claim);
void verify_seal(const Json &value, const std::string &name);
struct NativeLiteralObservation {
    NativeRecordValue value, binding;
};
// Read the existing immutable accession, without constructing a new field or
// claiming an encoding component for this materialized view.
NativeLiteralObservation read_native_literal(const NativeRecord &value, const std::string &type,
                                             std::optional<std::size_t> coordinate = {});
bool compare_native_literal(const NativeRecord &actual, const Json &expected,
                            const std::string &type, const std::string &comparison);
} // namespace gadget::claim_detail
