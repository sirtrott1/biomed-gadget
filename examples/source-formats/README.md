# Source format readers

These readers extract the scientific records in BindingDB bulk TSV, PubChem
assay JSON and concise CSV, and declared GTEx GCT and metadata tables. The
format declarations remain in the examples; working measurements and evidence
enter GADGET through the general native data interfaces.

Readers preserve original identities, occurrence order, numerical tokens,
missing values and source coordinates. BindingDB's reported comparisons and
PubChem's typed readouts remain distinct from ordinary point measurements.

| Reader | Preserved information |
| --- | --- |
| [BindingDB](bindingdb.hpp) | Original fields, row byte spans, reported affinity endpoints and comparisons, variable target chain columns |
| [PubChem](pubchem.hpp) | Typed readouts, original numeric spellings, missing readouts, outcomes, revocations and separate SID and CID links |
| [GTEx](gtex.hpp) | Versioned gene or transcript identities, sample coordinates, observed zero, missing support and donor metadata |

## Streaming input

Enable `GADGET_BUILD_EXAMPLES` and link the `gadget_source_formats` target from
this directory. A BindingDB visitor can consume a binary input stream without
retaining the complete file:

```cpp
#include "bindingdb.hpp"
#include <fstream>

using namespace gadget::examples;
std::ifstream input("binding.tsv", std::ios::binary);
BindingDBBulkLimits limits{16 * 1024 * 1024, 100000, 4096, 1024 * 1024};
const auto summary = visit_bindingdb_bulk(input, limits, [](const auto &row) {
    // Consume the row's fields, affinities and original byte offset.
});
```

The stream buffer is reused after each callback. Copy any row data needed
after that callback returns. For input already held in memory, the same visitor
accepts a string view; keep those source bytes alive while using its views.
Row callbacks are provisional until the complete visit succeeds. A later
invalid row, truncated source or cancellation can invalidate an unfinished
import. Set `expected_source_bytes` when an independent source length is known.

GTEx requires an explicit layout, feature kind and measurement kind. Its portal
layout is a declared profile; it does not accept every GCT 1.3 variant. Metadata
joins preserve repeated donors. Statistical independence needs its own check.

PubChem inspection retains a complete declared assay and its decoded views in
memory. Set source, readout and link budgets before inspection. The optional
native adapter publishes typed evidence through the common evidence store.

For working datasets, use the [native data API](../../tools/rcdb/DATA.md).
Retain original files once as provenance, then store admitted fields and shared
context in bounded native records.

## Source checks

With `BUILD_TESTING=ON`, `gadget_bindingdb_bulk_tests`,
`gadget_pubchem_assay_tests` and `gadget_gtex_source_tests` check source grammar,
meaning, resource refusal and cancellation. `gadget_source_formats_reference_contracts`
compares authored fixtures with independent Python readers, decimal values and
original byte spans. It does not qualify an external dataset population.

The reference scripts also accept original BindingDB or PubChem input and its
expected SHA256. Run `python tests/bindingdb_bulk_reference.py --help` or
`python tests/pubchem_assay_reference.py --help` for arguments. Reader JSON
output is an inspection view; numerical work requires native admission.

## Native workflows

[Measurement analysis](ANALYSIS.md) describes bounded native imports, exact count
analysis, finite pathway tests and source evidence correspondences.
