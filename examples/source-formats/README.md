# Source format readers

These C++ examples inspect BindingDB bulk TSV, PubChem assay JSON and concise
CSV, and declared GTEx GCT and metadata tables. They retain original identities,
occurrence order, numerical tokens, missing values and source coordinates.

| Reader | Preserved information |
| --- | --- |
| [BindingDB](bindingdb.hpp) | Original fields, row byte spans, reported affinity endpoints and comparisons, variable target chain columns |
| [PubChem](pubchem.hpp) | Typed readouts, original numeric spellings, missing readouts, outcomes, revocations and separate SID and CID links |
| [GTEx](gtex.hpp) | Versioned gene or transcript identities, sample coordinates, observed zero, missing support and donor metadata |

## Use a reader

Enable `GADGET_BUILD_EXAMPLES` and link the `gadget_source_formats` target from
this directory. For example, a BindingDB visitor receives one original row at
a time:

```cpp
#include "bindingdb.hpp"
#include <gadget/io.hpp>

using namespace gadget::examples;
BindingDBBulkLimits limits{16 * 1024 * 1024, 100000, 4096, 1024 * 1024};
const auto bytes = gadget::read_text("binding.tsv", limits.source_bytes);
const auto summary = visit_bindingdb_bulk(bytes, limits, [](const auto &row) {
    // Inspect row.fields, row.affinities and row.byte_offset here.
});
```

Keep the source bytes alive and unchanged while using views. Row callbacks are
provisional until the complete visit succeeds. A later invalid row or
cancellation can invalidate an unfinished import.

GTEx requires an explicit layout, feature kind and measurement kind. Its portal
layout is a declared profile; it does not accept every GCT 1.3 variant. Metadata
joins preserve repeated donors. Statistical independence needs its own check.

PubChem inspection retains a complete declared assay and its decoded views in
memory. Set source, readout and link budgets before inspection. These functions
inspect original input; they do not fit models or publish native evidence.

For working datasets, use the [native data API](../../tools/rcdb/DATA.md).
Retain original files once as provenance, then store admitted fields and shared
context in bounded native records.

## Verify original input

With `BUILD_TESTING=ON`, `gadget_bindingdb_bulk_tests`,
`gadget_pubchem_assay_tests` and `gadget_gtex_source_tests` check source grammar,
meaning, resource refusal and cancellation. `gadget_source_formats_reference_contracts`
compares authored fixtures with independent Python readers, decimal values and
original byte spans. It does not qualify an external dataset population.

The reference scripts also accept original BindingDB or PubChem input and its
expected SHA256. Run `python tests/bindingdb_bulk_reference.py --help` or
`python tests/pubchem_assay_reference.py --help` for arguments. Reader JSON
output is an inspection view; numerical work requires native admission.
