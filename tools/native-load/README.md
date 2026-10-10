# Native load profiling

Configure `GADGET_BUILD_PROFILING_TOOLS=ON` and build the
`gadget_native_load_profile` target. The executable links `gadget::native` and
uses the RexGraphCpp series decoder, publication pins, accession directory and
GMP storage accounting functions.

```sh
gadget_native_load_profile --input study.rfcxs --record --accessions \
  --rewrite verified.rfcxs > load.tsv
```

The input must be a complete canonical RFCXS series. `--record` admits its
GADGET application declarations and materializes retained original byte owners.
Omit this option for a series without GADGET record declarations. `--accessions`
prepares the selected publication's indexes, then measures reuse of that same
directory. `--cutoff N` selects an actual publication sequence; the default is
the latest endpoint. Publication sequence is separate from biological time.

`--rewrite FILE` publishes the complete input history to a new file and compares
every output byte with the original. It preserves all frames and checkpoints,
including when an earlier cutoff is selected for inspection. An existing output
file is refused. The input is never overwritten.

The output is tab separated. `phase` rows contain the phase name, elapsed seconds,
current resident bytes and the process peak resident bytes observed after that
phase. Memory observations come from `/proc/self/status` on Linux; other systems
report `unavailable`. The peak is cumulative for the process. Phase times exclude
the subsequent resource accounting unless the phase is named `account`.

`value` rows contain a measurement name and integer value. History measurements
count encoded component vectors once by actual owner. Source measurements count
the selected endpoint's cells, ordered occurrences, incidence arrays, reverse
indexes, identities, aspect planes, scale storage and lattice coefficients.
Field measurements separate logical coordinates, structural support, nonzero
entries, retained and materialized value fields, vector capacity and GMP limb capacity.
Accounting inspects existing allocations without materializing an encoded field.
Zero entry and limb capacities indicate that the coefficients remain in the
counted history component; they do not indicate missing data. Declaration tensors have their
own counts. Record byte owners are deduplicated by shared owner identity.

These measurements do not form a complete allocation total. Allocator overhead,
transient decoder work, accession index allocations, other endpoints retained
by the history, declaration registry storage and local block storage are not
included in the component capacity counts. Original bytes stored in Integer
declarations are separate from their materialized byte owners. Full support does
not imply a stored value at every coordinate; a supported zero remains distinct
from absence.

Inspect the scientific partition publication when measuring an assay's
mathematical complex. An administrative collection publication can contain
independent scientific archives in its original byte owners; its own source
counts do not describe those nested complexes.

## Raw import profiling

The same build option supplies `gadget_raw_import_profile`. It imports a complete
assay specification into a new directory and reports completed stages:

```sh
gadget_raw_import_profile import.rfcxs new-workspace > import.tsv
```

Stages separate matrix reading, original artifact retention, source construction,
native series publication, study declaration and admission, workspace publication,
and final application validation. Repeated assays and partitions produce repeated
stage names in execution order. A stage time starts after the preceding report;
`total` measures the complete call, including reporting. The process peak is
cumulative and cannot be attributed to one stage's retained allocations.
Resource-only workspaces have no matrix stages. Imports retain all original
records, axes, support and coefficients; profiling does not select a subset.

Specifications use the native record format; JSON specifications remain accepted
through the compatibility reader. Original specification bytes are retained in
either case.

The C++ `import_assay_workspace` function returns the admitted workspace and its
native publication. Subsequent operations can use that owner without reopening
the output file. Existing `import_assays` callers retain their void interface.
Both APIs accept an optional observer of completed stages. Observer exceptions
abort the import and remove its unpublished destination.
