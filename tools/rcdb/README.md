# RCDB and RCQL

The [native data API](DATA.md) provides bounded original intake, paged catalogs,
source resolution, selected native fields, shared contexts, lookup across independent
axes and retained native calculation programs. Install it with
`python -m pip install ./tools/rcdb` after installing the matching providers.
It also installs `gadget_rcdb` for complete publication compatibility.

There are two storage interfaces. `gadget_data` publishes datasets in bounded
records so queries can retrieve selected fields. `gadget_rcdb` stores a complete
RFCXS publication as one native record for version selection and replay. Choose
that interface when the complete archive fits the declared record budget.

Both use Core, RCDB and RCQL from the active Python environment. The compiled
GADGET executable admits complete publications and runs GADGET analysis and
replay. RexGraph supplies the format, storage and query implementations.

The independent `rexgraph_cpp` bindings expose C++ SDK operations to Python.
Their source is in the top level `rexgraph_cpp/` directory of `rexgraph-python`;
`bindings/cpp/` builds the bridge against an external `rexgraphcpp` SDK.
Agent and System provide optional services and a browser application. The main
GADGET C++ application uses the SDK directly; these Python integrations have
their own installation and qualification profiles.

Follow the [three repository installation guide](../../DEPENDENCIES.md) to
build Core, RCDB, RCQL, the bridge and GADGET from their source checkouts using
one SDK prefix. The dataset adapter requires Core, RCDB and RCQL 1.3.2 or newer.
Install `./tools/rcdb[native]` after the bridge for native C++ Python operations.
Agent and System are not dataset dependencies.

The qualified CPython 3.14 Linux x86_64 profile is recorded in
[`requirements.lock`](requirements.lock). It selects six exact independent
wheels: RexGraph core, native bindings, RCDB, RCQL, Agent and System. Core,
RCDB, RCQL, Agent and System use version 1.3.2; the bridge uses version 0.1.0.
Use the [source installation](../../DEPENDENCIES.md#install-from-source) for
the public GitHub repositories. The lock file instead checks the exact recorded
wheel artifacts when those files are available; rebuilding the same source can
produce different wheel bytes. Install matching artifacts in a dedicated
environment:

```sh
python3.14 -m venv .venv-rcdb
.venv-rcdb/bin/python -m pip install \
  'numpy>=2' 'threadpoolctl>=3' 'scipy>=1.13' 'pandas>=2.2.2' \
  'safetensors>=0.3' 'fastapi>=0.110' 'uvicorn>=0.29'
.venv-rcdb/bin/python -m pip install --no-index --no-deps --require-hashes \
  --find-links /absolute/path/to/rexgraph-wheels -r tools/rcdb/requirements.lock
.venv-rcdb/bin/python -m pip check
```

The hash check covers the six project wheels; third party dependencies above
are installed separately. This profile includes the public Python packages and
independent native bridge. The native binding and GADGET use the same complete SDK header
inventory. The SDK remains an independent dependency. GADGET's separately built
C++ reader must pass the adapter contracts against the installed Python profile.
Rebuilding a native binding against another SDK requires its own qualification
and checkpoint compatibility checks.

Run the adapter contracts with this environment's interpreter. A noneditable
wheel installation allows imports to be checked from outside the source checkouts:

```sh
.venv-rcdb/bin/python -m pip show rexgraph rexgraph-cpp-bindings rexgraph-rcdb rexgraph-rcql rexgraph-agent rexgraph-system
.venv-rcdb/bin/python -I -c "import rexgraph, rexgraph_cpp, rcdb, rcql, agent, system; print(*(m.__file__ for m in (rexgraph, rexgraph_cpp, rcdb, rcql, agent, system)), sep='\n')"
```

To require the Python native binding and GADGET to use the same SDK headers,
check the selected environment against GADGET's recorded build inventory:

```sh
.venv-rcdb/bin/python -I tools/rcdb/check_native_stack.py \
  --build-info build/presets/formats/build-info.json --report native-stack.json
```

This checks every compiled `rexgraph_cpp` module's SDK receipt, including storage
and model providers. It refuses a different SDK snapshot even when package
versions match or archive interoperability tests pass. It reports the actual
project import paths and versions. It does not rebuild dependencies, change model
engine receipts or establish compatibility with an older checkpoint. A retained
Python runtime qualified against another SDK will refuse this same SDK check;
retain that runtime for its existing model records.

The default check requires Core, RCDB, RCQL and the bridge. Add `--services`
when Agent and System are part of the selected application environment.

Existing model checkpoints retain their original native engine identity. The
current graded transport provider refuses execution of a checkpoint requiring a
different engine. Its explicit `replay_graded_stalk_transport_model(state, source)`
API can produce a new state by replaying the original observations and checking
the complete retained transport family against the current provider. Preserve
the original checkpoint, source bytes and runtime; replay does not rewrite them
or automatically migrate stored records.

These native publication contracts require the Python core and storage/query
packages. Zstandard tests also require `zstandard`; zlib is built into Python.
Agent services, native Python bindings, model backends and physical device
execution have separate prerequisites and qualification.

The complete publication adapter stores the unchanged RFCXS archive. Its owning relational complex,
occurrences, exact fields, declaration channel, raw byte sources and publication
history remain in the C++ format. Admission invokes GADGET's native reader and
requires byte identical canonical re encoding. RCDB provides version selection,
integrity binding and durable storage through its existing record engine.

After installing the adapter and RCDB in the same environment, use the compiled
GADGET executable to admit complete publications:

```python
from rcdb import LocalStore
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from gadget_rcdb import install_codec, put_publication, read_publication

codec = install_codec(gadget="/absolute/path/to/gadget", max_bytes=256 * 1024 * 1024)
store = LocalStore("experiment.rcdb", compression=BlobCodecSpec(name="zstd"),
                   record_codecs=(*DEFAULT_RECORD_CODECS, codec))
try:
    put_publication(store, "coverage", "coverage.rfcxs", expected_version=0)
    snapshot = read_publication(store, "coverage", version=1)
    assert snapshot.value == open("coverage.rfcxs", "rb").read()
finally:
    store.close()
```

The example admits archives up to 256 MiB. `max_bytes` is supplied by the caller.
The installed RCDB `RecordCodec` applies
its own admission bounds. Archives outside that profile remain ordinary native
RFCXS files; this adapter does not partition a complex to satisfy a store limit.
The example selects RCDB's existing lossless Zstandard compression and requires
`zstandard`. `BlobCodecSpec(name="zlib")` selects zlib instead. Compression is
sealed in the store header; reopening never chooses a different codec implicitly.
Existing stores must already declare the publication codec in their immutable headers.
Every process reading these records must explicitly install the provider.

Repeated codec setup with the same resolved executable and byte budget reuses
the installed provider, including after reloading the adapter. A different
executable, budget, object type or admission capability refuses without changing
the occupied codec. `replace=True` explicitly replaces that process capability;
it does not change a store's immutable header. Initialize providers before
concurrent use.

The native record bridge uses the existing MemoryStore, LocalStore, SQLite
SQLStore and NativeObjectStore file and memory providers. Legacy FileStore
refuses generic native records. The backend contract checks complete bytes,
literal record identities, historical and valid time selection, tombstones,
conditional writers, the canonical copy boundary, typed RCQL results and
checkpoint/reopen behavior:

```sh
python -I tests/rcdb_codec_registration.py /path/to/gadget coverage.rfcxs
python -I tests/rcdb_backends.py /path/to/gadget coverage.rfcxs new-backend-tests
```

The backend test requires SQLAlchemy, fsspec and Zstandard. Select
`--compression none` or `--compression zlib` to exercise those declared profiles.
Repeated `--backend` arguments select providers. Memory providers retain state
within the process; file, local and SQLite histories are also reopened in a
fresh process. Remote object providers require their own publication capability.
`--result-max-bytes` declares a portable RCQL result budget when testing an
archive larger than the result codec's default profile.

For a whole native file, including one outside the RCDB record size profile,
copy it through the installed exclusive file writer:

```python
from gadget_rcdb import copy_publication_file

copy_publication_file(
    "source.rfcxs", "new-copy.rfcxs",
    gadget="/absolute/path/to/gadget", max_bytes=your_file_byte_budget,
)
```

The writer defaults to `gadget_publish_file` beside the resolved native
executable; `publisher=` can select its explicit installed path. Copying captures
one opened original file into a private fixed snapshot, admits that snapshot
with the selected native reader and compares every published byte against its
opened owner. Python does not allocate a complete archive buffer. Temporary disk
space must hold the capture and native admission output. A changed source during
capture or a substituted snapshot during admission refuses before final output.
Later changes to the original file locator do not change the captured owner.
Existing destinations, including links, remain unchanged. The underlying writer
requires a trusted parent and hard link publication support; it supplies no
fsync durability or automatic process stop recovery.

The file path retains the selected provider's admission contract. A native series
from a different application or SDK profile must pass the actual GADGET reader
before this API can copy it as a GADGET publication. It does not wrap unknown
artifacts, infer mathematical cell counts or weaken a socket admission budget.
`tests/rcdb_file_copy.py` checks fixed source ownership, fresh process replay,
exclusive writers, real output exhaustion and complete large files with measured
Python allocations.

`select_publication_field_file` offers the same file ownership route for original
field selection. Supply native source and request files, a fresh destination,
the existing `gadget_native_field_query` executable and separate source, request
and output byte budgets:

```python
from gadget_rcdb import select_publication_field_file

select_publication_field_file(
    "source.rfcxs", "request.rfcxs", "new-selection.rfcxs",
    selector="/absolute/path/to/gadget_native_field_query",
    gadget="/absolute/path/to/gadget",
    max_source_bytes=your_source_byte_budget,
    max_request_bytes=your_request_byte_budget,
    max_output_bytes=your_output_byte_budget,
)
```

The inputs are captured and admitted before invoking the native selector. The
result retains the complete original selected publication, including every
field, carrier, support flag and byte owner; a window changes the declared
selection and its counts, rather than partitioning that mathematical source.
The source namespace, native cutoff, tensor key, carriers, interpretation and
optional window execution limits stay in the existing native request schema.
Complete native output admission and byte comparison precede returning the new
file. Output budget exhaustion leaves the requested destination absent.
`tests/rcdb_field_files.py` checks complete source retention, independent window
counts, changed controls, exclusive publication and measured Python allocations.
Native reader and field computation memory remain separate from those Python
allocations. This API neither creates a bulk encoded field index nor raises the
RCDB record size limit.

The native biomedical query contract runs the existing gluing, field transfer
and allocation producers through exact RCQL workspace and request versions.
It selects every retained tensor and its final row window, compares the complete
original prefix and all tensor values, carriers and support in C++, and checks
foreign owners, future cutoffs, changed latest versions and store reopening:

```sh
python -I tests/rcdb_biomedical_queries.py \
  /path/to/gadget /path/to/gadget_native_field_query \
  /path/to/gadget_rcdb_field_request /path/to/gadget_native_field_tests \
  /path/to/gadget_native_biomedical_fields_tests /new/test/directory
```

Run this optional contract with the matching installed RCDB and RCQL packages.
Its independent generated inputs qualify native transport and typed execution;
complete biological inputs and socket upload admission have separate scopes.

An installed Agent can retrieve a publication that its local operator has already
admitted to the store. Use the existing record interface and select the literal
record ID and version:

```python
from agent.client import RexClient

client = RexClient(server_url, api_key=your_token, workspace=your_workspace,
                   frame_key=your_frame_key)
packet = client.rex_fetch_record(record_id, version=1)
archive = packet.snapshot().value
```

The receiving process must install the GADGET codec before opening the snapshot.
The packet retains exact metadata, store identity, record version and unchanged
RFCXS bytes. RCDB record versions and native publication cutoffs are separate.
Time selection uses `as_of` and `valid_at`; an exact version cannot be combined
with either time selector. The graph frame endpoint accepts Python graph objects
and does not interpret an RFCXS publication.

Uploading a GADGET packet through `rex_store_record` requires a declared native
cell admission capability. This adapter supplies byte bounds and canonical
archive validation, but currently supplies no cell counter. The Agent record
route refuses such uploads before native decoding or destination publication.
Local storage and retrieval do not qualify this missing upload capability.

The optional Agent contract uses its installed record routes, authentication,
workspace and compute budget middleware over an owned loopback listener. It
checks exact retrieval, upload refusal, complete bytes and restart. It requires
matching Agent and storage packages with the Agent server dependencies:

```sh
python -I tests/rcdb_agent.py /path/to/gadget /path/to/publication.rfcxs \
  /new/test/directory --compression zstd
```

The test creates its own configuration, identities, journals and stores. It
qualifies these record routes on the local host; other Agent services, remote
hosts and mathematical cell admission require their own checks.

Use RCDB's existing `copy_record` to transfer a selected native publication to
another store. Both stores must declare the GADGET codec. The destination can
use a different declared compression codec:

```python
from rcdb import copy_record

selected = store.read_record("coverage", version=1)
receipt = copy_record(
    store, destination_store, selected.record,
    expected_digest=selected.state_digest,
    expected_cursor=destination_store.change_cursor,
    preserve_signature=True, return_receipt=True,
)
```

The source version is selected before destination publication. The digest binds
that source object, and the destination cursor refuses a concurrent change.
The copy retains native bytes, metadata, tags and validity intervals while
allocating a destination record version. Its `CopyReceipt` names both store
identities, record addresses, versions and object digests; `to_bytes` and
`CopyReceipt.from_bytes` preserve that receipt. Native publication cutoffs stay
inside the unchanged archive. Copying a record does not transplant a store's
history or mutation commits. Generic native copies use the existing ungoverned
profile; use the destination cursor to make publication conditional.

`LocalStore.checkpoint()` creates a checked replay cache and retains the
authoritative journal and historical payloads. After any publication write error
or process interruption, close the affected handle and inspect the destination
before retrying: the journal can be published before a receipt is returned.
An ordinary `OSError` can also leave a torn tail; the exception type alone does
not establish whether a write committed. Ordinary opening refuses an incomplete
tail. With cooperating writers stopped, use the existing `LocalJournal.inspect()`
to identify its verified prefix and `repair_torn_tail()` to remove only the
incomplete final block. A damaged complete frame refuses repair. Keep the damaged
bytes for diagnosis before an explicit repair. Checkpoint corruption also refuses.
The local profile requires trusted storage directories and cooperating provider
writers.

Deleting a native record retains its earlier versions and creates a gap in its
transaction time selection. Recreating that literal record ID allocates the next
version; it does not reuse version one. An already returned snapshot keeps its
original bytes and metadata. The archive's native source namespace and cutoff
remain independent of the RCDB record version.

`tests/rcdb_recovery.py` checks these histories through fresh process checkpoint
replay and all three local compression profiles. Its POSIX failure controls use
actual file size write refusals at the journal append and checkpoint staging
boundaries, verify the committed prefix or previous checkpoint, then reopen and
retry. It also checks every incomplete length of a final journal frame and refuses
repair of damaged complete frames. Run it with a fresh output directory:

```sh
python tests/rcdb_recovery.py /path/to/gadget publication.rfcxs recovery-tests \
  --compression zstd
```

The `--original` profile checks the complete supplied publication through deletion,
recreation and fresh process replay, while omitting the small frame failure corpus.
These are local provider checks; they do not establish power loss behavior or
automatic recovery of a deployed service.

With an open store, RCQL retrieval returns a byte result:

```python
from rcql import BoundSource, SourcePolicy, call, query, run, source
from gadget_rcdb import install_rcql_operator

install_rcql_operator()
result = run(
    query(source("db"), call("GADGET_PUBLICATION_BYTES", "coverage", 1)),
    sources={"db": BoundSource(store, SourcePolicy.allow("identity", "records"))},
)
archive = result.values[0]
```

Omit the version argument to select the current record. RCQL records the selected
RCDB record identity, version and state digest. Its existing binary result codec
can retain the returned bytes and execution receipt. This operation retrieves a
publication; it does not claim that RFCXS has become a Python `RexGraph` or that
RCQL field operators accept an opaque archive.

The same provider can register RFCXS admission with RexGraph's file format
registry. It leaves the existing `.rcbd` and `.rcbf` handlers in place:

```python
from rexgraph import io
from gadget_rcdb import install_format

install_format(gadget="/absolute/path/to/gadget")
archive = io.load("coverage.rfcxs")
```

This reader returns validated native publication bytes. GADGET's native writer
provides binary publication; no Python tensor codec is added.

An application can submit the retrieved workspace and a native request to C++:

```python
from gadget_rcdb import execute_publication
result_bytes = execute_publication(workspace_bytes, request_bytes,
                                   gadget="/absolute/path/to/gadget")
```

The CLI equivalent is:

```sh
gadget assay-query --workspace workspace.rfcxs --request request.rfcxs \
  --out-record coverage.rfcxs
```

This command writes an exclusive native publication and produces no JSON display.
Use `gadget record-view --input coverage.rfcxs --member result` when a display is
needed. Controlled coverage publishes its physical segment and original interval
carriers, complete membership occurrences, and typed numerical fields directly.
Other recipes retain their existing execution implementations.

The `.rcbd` and `.rcbf` readers in RexGraph have different contracts. Renaming an
RFCXS archive would not make it compatible with those readers. This adapter uses
the existing C++ format without introducing another binary wire format.

For RCQL execution, store both the workspace and its request with the GADGET codec,
then explicitly install the executor:

```python
from gadget_rcdb import install_rcql_executor

install_rcql_executor(gadget="/absolute/path/to/gadget")
operation = run(
    query(source("db"), call("GADGET_ASSAY_QUERY", "workspace", 1, "request", 1)),
    sources={"db": BoundSource(store, SourcePolicy.allow(
        "identity", "records", "files", "gadget_execute"))},
)
result_bytes = operation.values[0]
```

A stored native scientific operation can also be replayed from its retained
original study and request. This calls GADGET's existing `replay-scientific`
provider; it returns a fresh native operation, preserving the complete source,
fields, supported zeroes and missing support. No numerical display is imported.

```python
from gadget_rcdb import install_rcql_scientific_replay

install_rcql_scientific_replay(gadget="/absolute/path/to/gadget")
replayed = run(
    query(source("db"), call("GADGET_SCIENTIFIC_REPLAY", "operation", 1)),
    sources={"db": BoundSource(store, SourcePolicy.allow(
        "identity", "records", "files", "gadget_replay"))},
)
replayed_bytes = replayed.values[0]
```

The record version is required. The receipt records the selected literal ID,
version, source state digest and replayed byte digest. `replay_publication(bytes,
gadget=...)` offers the same native delegation without an RCQL query. A native
archive without the admitted scientific source/request contract refuses replay.
Input archives and complete returned archives undergo native admission and exact
canonical byte checks. Replay uses a private temporary directory and leaves the
RCDB source unchanged; publish returned bytes explicitly when they should persist.

The same replay provider accepts retained scientific plans. Each ordered step
keeps its complete native operation and original study/request owners; the native
provider refuses missing steps, changed step identities and mixed study sources.
Fresh replayed child publications can have different operational namespaces.
Compare their actual fields, carriers, support and original source bytes rather
than requiring their enclosing archive buffers to match.

`tests/rcdb_scientific_plan.py` checks old and current plan versions under all
three local compression profiles, repeated replay, checked restart, malformed
plan refusal and canonical `copy_record` transfers. Its fixture uses public
`ScientificStudy`/`ScientificSession` APIs and original native query declarations.
The native comparison consumer checks every ordered child tensor and its exact
original study/request bytes.

Both input versions are required. The executor admits their native archives,
submits the retained request to GADGET, validates its native output, and records
the selected input versions and output byte digest. The operation requires an
explicit execution capability and declares temporary filesystem use. Numerical
execution and source ownership stay in GADGET and RexGraphCpp.

RCQL provider installation is process local. Repeating an installer with the
same declaration and resolved executable paths reuses the installed callable.
A different executable binding, conflicting declaration or incomplete
installation raises `ValueError` and preserves the existing registry entries.
Initialize these providers before concurrent query execution. If setup was
interrupted between declaration and callable registration, initialize a fresh
process; the adapter does not replace or repair an occupied operator name.

Native and compressed file format setup also reuses the installed reader when
its resolved executable, read only format declaration and declared byte budgets
match. A different reader, saver, executable, extension owner or budget raises
`ValueError` before changing the registry. Compression budgets must be positive
integers at setup; booleans are refused. The adapter inspects RexGraph's existing
format registry and keeps no separate registry. To change a registered file
reader, explicitly remove its name through `rexgraph.io.unregister_format()`
before installing the new profile. Initialize formats before concurrent use;
the underlying registry does not provide an atomic compare and register operation.

`tests/rcdb_format_registration.py` checks repeated setup, adapter reload, resolved
path aliases, occupied names and extensions, incomplete format declarations,
budget refusals and exact native/compressed reads. It requires a native archive,
the GADGET executable and its `gadget_publish_file` sibling, or an explicit
`--publisher` executable, with a fresh output directory.

Execution stages each complete input once. The selected native execution command
admits its original inputs through the canonical reader; the bridge avoids a
separate input validation and re encoding process. Output and file format
admission compare the canonical file with the captured complete byte buffer
without staging another archive copy. Every byte is checked in bounded transfer
buffers; native admission and complete result ownership remain unchanged.

Record operators remain observable reads. Repeated expressions select their
records again and retain each selected version and state digest in the trace.
They do not opt into query local memoization or persistent result reuse: current
versions and source policies can change, and native execution has filesystem
effects. RCQL's existing result cache refuses these store bound operations.
Explicit historical versions remain available after a newer publication.

To run the adapter contracts in an environment with RCDB and RCQL installed:

```sh
ctest --test-dir build/presets/cpu -R gadget_native_coverage_contracts
python3 tests/rcdb_adapter.py build/presets/cpu/gadget \
  build/presets/cpu/native-coverage-contracts/coverage.rfcxs \
  build/presets/cpu/rcdb-adapter-contracts
```

Use a fresh test output directory. The native execution contract accepts the
executable, a native assay workspace, a native assay request and a fresh output
directory as arguments to `tests/rcdb_native_execution.py`.

The registration contract checks repeated setup, executable bindings and partial
installation in fresh processes:

```sh
python3 tests/rcdb_registration.py build/presets/cpu/gadget \
  build/presets/cpu/gadget_native_field_query
```

The native transfer contract accepts a canonical publication and a fresh output
directory. It checks compressed destinations, checkpoints, conditional writers
and POSIX process interruption around local publication:

```sh
python3 tests/rcdb_transfer.py build/presets/cpu/gadget \
  coverage.rfcxs new-transfer-contracts
```

Native files outside the generic record profile can use RCDB's existing
compressed blob framing through the file adapter:

```python
from gadget_rcdb import compress_publication, read_compressed_publication

compress_publication("coverage.rfcxs", "coverage.rfcxs.rcblob",
                     gadget="/absolute/path/to/gadget", codec="zstd", level=3)
archive = read_compressed_publication(
    "coverage.rfcxs.rcblob", gadget="/absolute/path/to/gadget",
    max_encoded_bytes=your_encoded_byte_budget,
    max_output_bytes=your_decoded_byte_budget,
)
```

Compression passes a read only mapping of the fixed admitted source to RCDB's
existing compressor. It avoids a complete raw Python archive allocation and
produces the same codec frame as byte input. The compressor still materializes
its encoded output and may allocate its own output workspace; the mapping also
uses address space and resident file pages. Original file identity is checked
through compression, and final output uses the existing exclusive writer and
complete byte comparison. `tests/rcdb_compression_profile.py` compares both codec
frames and all decompressed original bytes against the prior byte input profile.

The adapter calls RCDB's `compress_blob` and `decompress_blob` functions. It
retains the complete RFCXS publication, including raw source bytes, and requires
canonical native admission. Compression changes storage bytes; decoding restores
the original archive byte for byte. The destination must be new. Encoded and
decoded admission budgets are supplied by the caller. This file transport does
not change the installed RCDB generic record bound.

Compressed output uses `gadget_publish_file`, normally installed beside the
resolved `gadget` executable. Set `publisher="/path/to/gadget_publish_file"` to
select another installed location. If that helper is unavailable, compression
refuses before creating output.

The helper passes the complete encoded byte stream to the SDK's staged
`CreateNew` operation. Short input, trailing input and ordinary write, flush or
close failures refuse without publishing a partial destination; retry uses the
same path. Concurrent exclusive writers produce one complete destination and
refuse the others. Existing files, directories and symbolic links are preserved.
After publication the adapter compares every output byte against the complete
encoded buffer, using bounded reads rather than another whole output copy.

This local profile requires trusted parent directories and filesystem hard link
support. It observes write, flush and close completion; it does not promise
fsync durability. An abruptly stopped process can retain an unpublished staging
directory. Its final destination stays absent or complete, and recovery of owned
staging residue is a separate operator action.

The generic helper also accepts a caller declared exact byte stream:

```sh
gadget_publish_file --out new-file.bin --bytes 123 < complete-input.bin
python -I tests/file_publication.py /path/to/gadget_publish_file \
  /path/to/file-publication-contracts /new/test/directory
```

The byte count must equal the complete input size. The helper stores opaque bytes
and does not infer a scientific format, build a complex or change its fields.

File admission captures one opened regular file within its measured size. RCDB
record admission checks the installed codec's byte budget on that descriptor
before reading; compressed admission checks the caller's encoded byte budget.
Growth and truncation during capture refuse before store publication or native
admission. Nonregular inputs refuse, including FIFOs without a writer on POSIX.
The captured bytes still undergo the existing complete native validation.
Record and compressed byte capture also check the opened file's identity, size
and modification metadata across the read, including changes of the same size.
These checks do not provide an atomic filesystem snapshot; inputs must remain
immutable during admission.

`install_compressed_format(gadget=..., max_encoded_bytes=...,
max_output_bytes=...)` registers `.rcblob` files with RexGraph's format registry.
It returns validated RFCXS bytes from `io.load`. GADGET's C++ executable accepts
the decoded `.rfcxs` file. This provider uses RCDB's whole buffer compression
API; compressed file size does not describe the memory required to admit the
complex or execute its field operations.

## Native field selection

`gadget_native_field_query` selects a tensor from an actual GADGET native
publication. The request names its exact source namespace, publication cutoff,
tensor key, nullable input/output carrier keys and interpretation:

```json
{
  "schema": "gadget.rcdb.field.query.v1",
  "source_namespace": "811",
  "cutoff": 4,
  "tensor": 888,
  "input_carrier": null,
  "output_carrier": null,
  "interpretation": "coefficient selection fixture"
}
```

Use the actual values from the selected publication. Equal dimensions do not
identify carriers. Publication sequence is separate from biological time.

```sh
gadget_native_field_query --input study.rfcxs --request selection.json \
  --out new-selection.rfcxs
```

The output is an administrative native operation record. Its
`result.selected_publication` retains the complete unchanged selected prefix,
including all original fields, carriers, support and byte sources. Reopen those
bytes with `NativeRecord::decode` to obtain the actual mathematical owner; the
receipt's own complex is not that owner. The receipt records the selected tensor
and original global coordinate bounds. It does not create a Python tensor or a
new source from a numerical view. Composite biomedical assertions require their
own declared selectors and do not become tensor fields implicitly.

An optional `window` has `lower` and `upper` arrays with zero based, half open
bounds on every original tensor axis. It requires an `execution` object naming
`name`, `records`, `rational_slots`, `sparse_entries`, `scalar_work` and
`integer_bits` limits for the existing SDK execution policy. These limits govern
the provider's named reference containers and work; they are not a process heap
quota. Supported zero and absence remain distinct. The current native window
provider materializes source entries during preparation, so a small coordinate
window does not imply low peak memory. Full original source ownership remains
in the output.

RCQL selection uses exact RCDB publication and native request versions:

```python
from gadget_rcdb import install_rcql_field_selector

install_rcql_field_selector(selector="/absolute/path/to/gadget_native_field_query",
                           gadget="/absolute/path/to/gadget")
selection = run(
    query(source("db"), call("GADGET_NATIVE_FIELD", "study", 1, "selection-request", 1)),
    sources={"db": BoundSource(store, SourcePolicy.allow(
        "identity", "records", "files", "gadget_field_select"))},
)
native_selection_bytes = selection.values[0]
```

The result type is `GadgetNativeFieldSelectionBytes`. RCQL records both selected
versions, their state digests and the native result digest. Native C++ retains
and validates the field owner; bytes remain distinct from Python mathematical
state. Existing RCDB and RCQL codec byte profiles still apply.

Native selection files use the SDK's `CreateNew` publisher. Competing writers
to one destination produce one complete winning file; later writers refuse
without changing it. Existing files, symbolic links and hard links refuse.
A failed staged write leaves no published output and can be retried at the same
destination. Parent directories must be trusted. This single file publication
does not provide fsync durability, a database transaction or crash recovery;
RCDB histories use their existing storage provider.
