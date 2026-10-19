#pragma once
#include "pubchem_model.hpp"
#include <gadget/biomedical.hpp>

namespace gadget::examples {
bool is_pubchem_evidence(const NativeRecord &publication);
void verify_pubchem_evidence(const NativeRecord &publication,
                            const rfc::database::InterpretedSnapshot &native);
NativeRecord query_pubchem_evidence(const NativeRecord &publication,
    const rfc::database::InterpretedSnapshot &native, const NativeRecord &request);
NativeRecord context_pubchem_evidence(const NativeRecord &publication,
    const rfc::database::InterpretedSnapshot &native, const NativeRecord &request);
} // namespace gadget::examples
