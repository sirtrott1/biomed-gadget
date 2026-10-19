#pragma once
#include <gadget/io.hpp>
#include <rfc/database/native_interpretation.hpp>
#include <optional>
#include <stop_token>

namespace gadget::examples {
inline constexpr std::uint64_t pubchem_binding_namespace = 72141;
inline constexpr std::uint64_t pubchem_declaration_namespace = 72142;
struct PubChemEvidencePacket {
    Json request, manifest, receipt, records;
    std::map<std::string, Json> entities;
    rfc::database::InterpretationSchema schema;
    rfc::database::LogicalSnapshot logical;
    rfc::database::InterpretationLimits limits;
};
// Original source readers supply the actual typed native declarations.
PubChemEvidencePacket prepare_pubchem_evidence(const Json &request,
    const std::string &assay, const std::optional<std::string> &concise,
    std::uint64_t source_namespace, std::stop_token stop = {});
rfc::database::InterpretedSnapshot interpret_pubchem_evidence(const PubChemEvidencePacket &packet,
                                                            std::stop_token stop = {});
void validate_pubchem_query(const Json &query, const std::string &digest);
void validate_pubchem_context(const Json &request, const std::string &digest);
} // namespace gadget::examples
