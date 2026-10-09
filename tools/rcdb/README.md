# RCDB and RCQL

`gadget_rcdb.py` connects an installed RexGraph Python stack to GADGET's native
publication reader. It imports RCDB and RCQL from the active Python environment.
It does not install either package or implement their storage codecs.

The RCDB payload is the unchanged RFCXS archive. Its owning relational complex,
occurrences, exact fields, declaration channel, raw byte sources and publication
history remain in the C++ format. Admission invokes GADGET's native reader and
requires byte-identical canonical re-encoding. RCDB provides version selection,
integrity binding and durable storage through its existing record engine.

Add this directory to your application's Python module path. With GADGET and
RCDB installed:

```python
from rcdb import LocalStore
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from gadget_rcdb import install_codec, put_publication, read_publication

codec = install_codec(gadget="/absolute/path/to/gadget", max_bytes=your_byte_budget)
store = LocalStore("experiment.rcdb", compression=BlobCodecSpec(name="zstd"),
                   record_codecs=(*DEFAULT_RECORD_CODECS, codec))
try:
    put_publication(store, "coverage", "coverage.rfcxs", expected_version=0)
    snapshot = read_publication(store, "coverage", version=1)
    assert snapshot.value == open("coverage.rfcxs", "rb").read()
finally:
    store.close()
```

`max_bytes` is supplied by the caller. The installed RCDB `RecordCodec` applies
its own admission bounds. Archives outside that profile remain ordinary native
RFCXS files; this adapter does not partition a complex to satisfy a store limit.
The example selects RCDB's existing lossless Zstandard compression and requires
`zstandard`. `BlobCodecSpec(name="zlib")` selects zlib instead. Compression is
sealed in the store header; reopening never chooses a different codec implicitly.
Existing stores must already declare the publication codec in their immutable headers.
Every process reading these records must explicitly install the provider.

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

The same provider can register RFCXS admission with RexGraph's file-format
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

Both input versions are required. The executor admits their native archives,
submits the retained request to GADGET, validates its native output, and records
the selected input versions and output byte digest. The operation requires an
explicit execution capability and declares temporary filesystem use. Numerical
execution and source ownership stay in GADGET and RexGraphCpp.

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

Native files outside the generic-record profile can use RCDB's existing
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

The adapter calls RCDB's `compress_blob` and `decompress_blob` functions. It
retains the complete RFCXS publication, including raw source bytes, and requires
canonical native admission. Compression changes storage bytes; decoding restores
the original archive byte for byte. The destination must be new. Encoded and
decoded admission budgets are supplied by the caller. This file transport does
not change the installed RCDB generic-record bound.

`install_compressed_format(gadget=..., max_encoded_bytes=...,
max_output_bytes=...)` registers `.rcblob` files with RexGraph's format registry.
It returns validated RFCXS bytes from `io.load`. GADGET's C++ executable accepts
the decoded `.rfcxs` file. This provider uses RCDB's whole-buffer compression
API; compressed file size does not describe the memory required to admit the
complex or execute its field operations.
