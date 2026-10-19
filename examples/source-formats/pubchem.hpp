#pragma once
#include <gadget/io.hpp>
#include <stop_token>

namespace gadget::examples {
// Source format inspection only. This owner retains the original bytes; the
// JSON projection is an external view, not a native evidence publication.
struct PubChemSource {
    std::string bytes;
    Json original;
    Json projection;
};
struct PubChemSourceLimits {
    std::size_t source_bytes;
    std::size_t definitions;
    std::size_t substances;
    // Includes absent descriptor/substance coordinates, not only reported data.
    std::size_t readout_slots;
    std::size_t concise_rows;
    std::size_t links;
};
// Supports a full PC_AssaySubmit with a positive public AID/SID. External only
// identifiers, unknown value types and integers outside signed 64 bits refuse.
// Definitions, outcomes and revocations remain source declarations; this does
// not infer targets, modifier pairings, activity thresholds or value transforms.
PubChemSource inspect_pubchem_assay(std::string source, const std::string &expected_aid,
                                   const PubChemSourceLimits &limits,
                                   std::stop_token stop = {});
// Concise summaries contribute explicit identities only, never activity values.
// Repeated and conflicting links remain separate source occurrences.
PubChemSource inspect_pubchem_concise(std::string source, const std::string &expected_aid,
                                     const PubChemSourceLimits &limits,
                                     std::stop_token stop = {});
} // namespace gadget::examples
