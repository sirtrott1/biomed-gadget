# Complete publication storage

`gadget_rcdb` connects complete GADGET RFCXS publications to RCDB and RCQL.
GADGET's compiled reader validates the archive. RCDB stores its unchanged bytes
with record versions; RCQL retrieves or executes explicitly selected versions.

For datasets that should be read in bounded records, use [gadget_data](DATA.md).
This interface stores each complete archive as one record.

## Register the codec and store an archive

Run with the [matching providers](../../DEPENDENCIES.md) and a compiled GADGET
executable:

```python
from rcdb import LocalStore
from rcdb.header import BlobCodecSpec, DEFAULT_RECORD_CODECS
from gadget_rcdb import install_codec, put_publication, read_publication

codec = install_codec(gadget="/absolute/path/to/gadget", max_bytes=256 * 1024**2)
store = LocalStore("experiment.rcdb", compression=BlobCodecSpec(name="zlib"),
                   record_codecs=(*DEFAULT_RECORD_CODECS, codec))
try:
    put_publication(store, "coverage", "coverage.rfcxs", expected_version=0)
    snapshot = read_publication(store, "coverage", version=1)
    archive = snapshot.value
finally:
    store.close()
```

The byte budget applies alongside RCDB's own record bounds. Every reader process
must install the codec. Existing stores must already declare it in their headers.
Compression is part of the store declaration; reopening does not change it.
`zstd` additionally requires `zstandard`.

MemoryStore, LocalStore, SQLite SQLStore and NativeObjectStore support native
records. Legacy FileStore does not. See the providers' requirements for optional
SQL and object storage dependencies.

## Versions and copying

A record version and a native publication cutoff describe different things.
Record versions select stored objects; the cutoff remains inside the archive.
An exact version cannot be combined with `as_of` or `valid_at` time selectors.

Deleting a record preserves its history. Recreating that record ID allocates the
next version. Previously returned snapshots keep their original bytes.

Use RCDB's `copy_record` to transfer a selected version between stores that both
declare the codec:

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

The destination can use a different compression codec. The receipt retains both
store identities, versions and digests. Copying does not transplant store history.

## Query with RCQL

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

Omit the version argument to select the current record. The trace retains the
selected record identity, version and digest. Retrieval returns native archive
bytes; it does not turn an archive into a Python graph.

## Execute and replay

Supply workspace and request bytes to the compiled executor:

```python
from gadget_rcdb import execute_publication

result_bytes = execute_publication(
    workspace_bytes, request_bytes, gadget="/absolute/path/to/gadget",
)
```

The CLI equivalent writes a fresh native result:

```sh
gadget assay-query --workspace workspace.rfcxs --request request.rfcxs \
  --out-record coverage.rfcxs
gadget record-view --input coverage.rfcxs --member result
```

For versioned RCQL execution, store the workspace and request first:

```python
from gadget_rcdb import install_rcql_executor

install_rcql_executor(gadget="/absolute/path/to/gadget")
operation = run(
    query(source("db"), call("GADGET_ASSAY_QUERY", "workspace", 1, "request", 1)),
    sources={"db": BoundSource(store, SourcePolicy.allow(
        "identity", "records", "files", "gadget_execute"))},
)
```

Both versions are required. The executor admits their archives, executes the
request and validates the complete output. Publish returned bytes explicitly to
retain them in RCDB.

`install_rcql_scientific_replay(gadget=...)` registers
`GADGET_SCIENTIFIC_REPLAY(record_id, version)`. It requires `identity`, `records`,
`files` and `gadget_replay` capabilities. `replay_publication(archive, gadget=...)`
provides the same operation without RCQL.

Replay restores the original scientific study and request, including ordered
plan steps. Unsupported archives refuse. Fresh operational namespaces can
differ; compare fields, support and original sources rather than assuming equal
archive bytes.

## Files and compression

These APIs retain complete archives outside the generic RCDB record interface:

| Function | Use |
| --- | --- |
| `copy_publication_file` | Copy one admitted archive to a new file |
| `select_publication_field_file` | Apply a native field request and retain its source |
| `compress_publication` | Write an archive through RCDB blob compression |
| `read_compressed_publication` | Restore and admit the original archive bytes |
| `install_format` | Register RFCXS loading with `rexgraph.io` |

For a file copy:

```python
from gadget_rcdb import copy_publication_file

copy_publication_file(
    "source.rfcxs", "new-copy.rfcxs",
    gadget="/absolute/path/to/gadget", max_bytes=256 * 1024**2,
)
```

File selection requires explicit source, request and output byte budgets and the
`gadget_native_field_query` executable. Compression requires encoded and decoded
budgets when reading. These byte limits do not bound the native reader's entire
heap or raise RCDB's record size limit.

File copying and compression use `gadget_publish_file`, normally beside `gadget`.
They require fresh output paths, trusted parent directories and filesystem hard
link support. Publication checks complete bytes but supplies no fsync durability
or automatic recovery after a process is stopped.

The selected field file retains its complete original publication. A window
selects coefficients; it does not partition that owning mathematical source.
The [dataset API](DATA.md#materialization-and-native-analysis) publishes bounded
dataset selections through its separate interface.

## Registration and services

Install providers before concurrent use. Repeated setup with the same resolved
executables and declarations reuses the provider. Conflicting declarations,
executables or byte budgets refuse. File format changes require explicit removal
of the old format registration first.

Record operators perform observable reads. They select versions and retain a
trace each time; they do not use RCQL's persistent result cache. Execution also
uses temporary files.

An installed Agent can retrieve records already admitted to its store through
`rex_fetch_record(record_id, version=...)`. The receiving process must install
the GADGET codec before opening the snapshot.

The Agent upload route needs a mathematical cell admission capability that this
adapter does not supply. It refuses GADGET uploads. Local storage and retrieval
remain separate supported operations.

## Local recovery

`LocalStore.checkpoint()` caches replay while preserving the authoritative
journal. After a write error or interruption, close the handle and inspect the
destination before retrying; a write may have committed before its receipt was
returned.

With cooperating writers stopped, `LocalJournal.inspect()` identifies the
verified prefix. `repair_torn_tail()` removes an incomplete final block and
refuses a damaged complete frame. Preserve damaged bytes before explicit repair.
The local profile requires trusted directories and cooperating writers.

## Verification

[Testing instructions](../../tests/README.md) describe the shared environment.
The `tests/rcdb_*.py` contracts cover registration, backends, original byte
comparisons, native execution, replay, transfer, compression and recovery.
Optional Agent tests require its server dependencies and an owned listener.
