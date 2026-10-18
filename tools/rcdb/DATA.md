# Native data I/O

`gadget_data` separates immutable original files, native working records and
selected mathematical owners. Core supplies typed arrays and relational
complexes. RCDB supplies storage, record envelopes and versions. RCQL supplies
query execution. GADGET supplies source declarations, coordinate maps and
materialization into its native publications.

Install the matching RexGraph packages described in [README.md](README.md),
then install this adapter in the same environment:

```sh
python -m pip install ./tools/rcdb
```

The `gadget-data` command is also available as `python -m gadget_data`.
HDF5 intake requires `h5py`. Export to GADGET requires the matching
`rexgraph_cpp` binding and compiled `gadget_native_field_import` executable.

## Intake and exploration

Read an original typed matrix using physical tiles:

```sh
gadget-data import-dense counts.npy counts.native \
  --reader npy --name counts --quantity reported_count --units count
gadget-data explore counts.native
gadget-data provenance counts.native \
  --field measurements --rows 10:20 --columns 4,1
gadget-data select counts.native selected.native \
  --field measurements --rows 10:20 --columns 4,1
```

NPY intake uses a read only mapping. HDF5 intake uses an explicitly named local
dataset and hyperslabs. Virtual datasets and external links require a separate
source admission. The `--missing nan` declaration maps NaN to missing support;
without that declaration, nonfinite coefficients refuse. Original IEEE values
retain their dtype and bits in native working records. Exact arithmetic uses
the declared binary interpretation.

Table intake uses Core's `DatasetDeclaration`, record readers and relation
construction. `import_records` retains complete native relation blocks.
`import_measurements` selects declared numeric columns and retains their source
row records, labels and context. Tables use two bounded passes, first to count
rows and then to publish. Record number rules, absence tokens, participant order
and relation declarations come from the supplied declaration. Long form data
uses an explicit contiguous group declaration and carries complete groups across
reader batches. A repeated later group refuses. Original typed input records
retain occurrence order even when the mathematical provider normalizes its
representation. A reader refuses
unsupported relations rather than splitting or coalescing their occurrences.

```python
from gadget_data import DataStore

with DataStore("counts.native") as data:
    mapping = data.provenance("measurements", rows=(19, 10), columns=(4, 1))
    selected = data.select("measurements", rows=(19, 10), columns=(4, 1))
    field = selected.tensor()
```

Catalogs and provenance queries open administrative index pages, without
decoding numerical field blocks. Selected reads open overlapping native blocks
and axis records. Original ordinals identify rows and columns. Repeated labels
remain separate coordinates. Field and dataset identities prevent equal
dimensions from identifying independent carriers. Every dependency records its
actual store, record, version and state digest.

`tensor()` returns a native exact `TensorField` with an actual source owner and
retained dependencies. Missing support requires an explicit scientific operation
before constructing a fully observed field. A supported zero remains observed.
`iter_select()` visits physical field windows without assembling a whole matrix.

`find_axis()` streams label and literal context matches on original coordinates.
Axis writers can explicitly index labels and selected context roles through
RCDB vocabulary. Search verifies actual pinned axis values, including after a
later record version. False, zero and missing context remain distinct. Broad
searches stream their axis records within the metadata budget.

## Sparse occurrences

`import_coordinate` streams general Matrix Market coordinate files. It preserves
entry order, duplicate coordinates, explicit zeroes, exact decimal values and
large integers. The declaration distinguishes implicit missing entries from
implicit observed zeroes. Symmetry expansion and duplicate aggregation require
separate declared operations.

`occurrences()` selects original entries. `iter_occurrences()` streams entry
blocks without expanding the logical matrix or its implicit zeroes. Selected
sparse materialization preserves original entry identities and axis maps.
Sparse occurrences do not silently become dense arithmetic inputs.

## Materialization and native analysis

Selected datasets are published as bounded native tiles. Each output axis keeps
its map to the parent coordinates. Each field block keeps the native versions
it consumed. The parent catalog is pinned by its immutable state identity.
Original files are not copied into every field or result.

Export an admitted selection to GADGET:

```sh
gadget-data export-native counts.native selected.rfcxs \
  --field measurements --rows 10:20 --columns 4,1 \
  --namespace 700100 \
  --importer /absolute/path/to/gadget_native_field_import
```

The bridge passes coefficients and support through the native C++ binding and
the SDK tensor codec into RFCXS. Its coordinate declaration contains source
references and maps. Structural sources can be transported in native components
with explicit vertex and occurrence maps. GADGET checks complete ordered
coverage and assembles them with the canonical builder. It does
not carry numerical field results through JSON. GADGET reads the actual native
complex and field, validates the input digest and namespace, and preserves
support while binding the declared matrix shape. The resulting `NativeRecord`
works with native field accessions and `gadget_native_field_query`. Biomedical
study declarations and metrics remain explicit application inputs.

Typed labels outside JSON scalar types retain their canonical Core value
payload. Coordinate entries explicitly reference that native identity payload.

`GADGET_DATA_FIELD` is an RCQL accession for a versioned field record. Install it
with `gadget_data.query.install()`. A query such as
`FROM RCDB_VERSION($db,$record,1) RETURN GADGET_DATA_FIELD($columns)` returns a
source bound exact `TensorField`. Columns are local to the selected physical
record. The field retains original global column identities and record version.

## Original sources and export

`SourceResolver` maps retained source identities to current paths and checks
their SHA256 and byte extent. Relative locators permit relocation under another
root. Explicit identity to path mappings support sources with repeated file
names. Opening an original pins and verifies the actual file descriptor.
Native calculations use native coefficients; original bytes remain provenance.

```sh
gadget-data verify counts.native --source-root /absolute/path/to/originals
gadget-data export-bundle selected.native complete.bundle \
  --dependency counts.native --source-root /absolute/path/to/originals
```

Export checks the complete linked dataset closure and preserves native RCDB
identities and history. Original files are included only when a source root is
supplied. The location index contains paths and identities, without scientific
coefficients. A bundle without originals still needs source resolution for
original byte verification. Missing or substituted parent versions refuse.

New datasets use RCDB's native SQLite provider. The existing local file provider
is available through `DataWriter(backend="local")` and remains readable. Some
local file reader profiles reserve the format byte limit during record loading;
native SQLite avoids that reservation. Both use the canonical RCDB record
engine and native mathematical payloads. Catalog indexes are paged, and the
reader retains administrative pages and axis labels within a declared cache
budget. Field values are not cached. Intake publishes index pages incrementally
instead of retaining one reference object for every field block.

`Limits` controls encoded block, decoded block and selected array admissions.
Selected axis labels, context and provenance references are charged against the
selection budget before retention.
CLI options `--record-mib`, `--decode-mib`, `--selection-mib` and `--catalog-mib`
set those budgets. Default dense tiles use short row bands and as many columns
as the block byte budget permits. Tiles, coordinate metadata and native export work have separate
admissions. The Python export API accepts `working_bytes` for exact transport.
These budgets govern declared work; process RSS also includes the runtime,
codecs and mathematical consumers. HDF5 compressed chunks must fit the decoded
budget. New dataset and export directories include an ignore file so generated
data is not added to a source repository accidentally.

The existing complete publication and assay interfaces remain available for
their declared formats. They retain complete archives. Use the dataset APIs for
bounded original intake, exploration and selected native materialization.

## Contracts

Run the native C++ contract through CTest. In an environment with the matching
Python providers and `h5py`, run the complete I/O contracts:

```sh
python tests/data_io.py /absolute/path/to/fresh-test-output \
  --importer build/presets/cpu/gadget_native_field_import \
  --consumer build/presets/cpu/gadget_data_field_tests \
  --selector build/presets/cpu/gadget_native_field_query
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python tests/data_io_scale.py /absolute/path/to/fresh-scale-output
```
