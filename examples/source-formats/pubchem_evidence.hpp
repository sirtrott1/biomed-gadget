#pragma once
#include <gadget/biomedical.hpp>
#include <stop_token>

namespace gadget::examples {
// Admit complete retained PubChem assay/optional concise byte owners. The
// request declares source identity, AID, original hashes and resource limits.
// This produces a typed native source, never a binding endpoint projection.
// The source/request owners (including their histories) remain archived;
// EvidenceStore queries retain this publication and the actual query request.
std::shared_ptr<const EvidenceAdapter> pubchem_evidence_adapter();
NativeRecord admit_pubchem_evidence(const NativeRecord &source, const NativeRecord &request,
                                   std::stop_token stop = {});
} // namespace gadget
