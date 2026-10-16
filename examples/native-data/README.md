# Native data and declared processing

The `gadget::native` CMake target provides `<gadget/native_data.hpp>`,
`<gadget/native_record.hpp>` and `<gadget/native_gluing.hpp>`. Its execution and linkage do not depend on JSON.
RexGraph owns the relational complex, tensor fields, declaration registry,
binary codecs, prepared mathematical actions and machine programs. GADGET
binds original data bytes and coordinates admitted work.

`NativeStudy::read(path)` detects RFCXS and RFCX from their canonical byte
signatures, regardless of file extension. RFCXS supplies its actual namespace;
an optional cutoff selects an existing endpoint. RFCX requires an explicit
namespace. Unknown bytes require an explicitly authored ingestion adapter and
are refused. Format detection supplies no biological interpretation.

`NativeStudy::from_publication(endpoint)` retains an actual SDK series pin and
prepares accession indexes on the first accession request. Reading its endpoint
or fields keeps that preparation deferred. `NativeStudy::read` uses this path
for RFCXS. Constructors accepting explicit accession limits validate and prepare
immediately. The directory is shared after preparation; its access operations
follow RexGraph's external synchronization contract.

Custom types belong in `TemporalEndpoint::declarations`: define the
`DeclarationRegistry`, attach `TemporalDeclaredValue` records to actual cells
or occurrences, bind source aspect factors, and retain authored
`TemporalSourceLayoutProgram` instructions. Persist the complete endpoint with
RexGraph's native RFCXS codec. GADGET accepts that endpoint or its canonical
binary directly, without recreating its topology from a semantic catalogue.

Original input bytes are `RawByteSource` owners, including arbitrary binary,
NUL bytes and invalid UTF8. A provenance record is an authored native Product
with exact Integer factors for source identity, byte offset and byte length.
`RawSpanDeclaration` names those stable factors and its actual type/version.
`NativeStudy::raw_span` resolves the native value against its explicitly bound
original bytes. The returned span owns both the bytes and the native accession;
it survives source file edits and destruction of the admitting wrapper. Raw
sources require distinct declared identities even if their contents are equal.

`decode_raw_span` applies an authored `CompiledLayout` directly to such a span
through RexGraph's `decode_declared_bits`. Integer, Rational, Product, Alternative,
Optional, Repeated and RationalTensor values use the canonical composite codec.
The resulting value retains the actual layout, native accession and original
bytes. This is declared format processing, rather than a GADGET parser that
guesses or normalizes the source. Arbitrary external formats still need an
explicit ingestion adapter or a matching declaration.

```cpp
#include <gadget/native_data.hpp>

auto study = gadget::NativeStudy::read(
    "experiment.rfcxs", {gadget::RawByteSource::read(83, "measurements.bin")});
gadget::RawSpanDeclaration binding{{71, 4, 2}, 500, 900, 700};
auto provenance = study.raw_span(33, binding);

auto program = study.program(27, 16, {{1}, {1000}});
rfc::composite::BinaryExecutionWorkspace workspace;
auto result = program.materialize(workspace);
```

Those identities refer to an authored example declaration. They are not global
GADGET keys. The caller supplies its actual declaration, data sources and
instruction word. `PreparedNativeSourceProgram` reads the original packed
aspect planes and executes the retained instructions with private scratch.
It preserves inactive occurrences, occurrence order and final word masks.
No adjacency construction or replacement algorithm is introduced.

`execute_native_gluing` accepts a prepared native gluing plan, authentic field
pins, explicit availability and a canonical CPU reservation. Results retain
the native preparation, source inputs and ordered boundary lineage. There is
no report generation in this operation. `native_gluing_endpoint` requires an
explicit output tensor/carrier identity for every chain and preserves supported
zeros. Persist the output through RFCXS and the original source square through
`write_native_atlas`; keep original raw files alongside those binaries.
Output publication namespaces remain explicit. This does not imply a signature
recoding or transport of arbitrary declared stalk payloads through the pushout.

The native contract executable exercises two different authored machine
programs over 130 original occurrences, custom factor ordering, raw opaque byte
provenance, immutable ownership, native byte detection/refusals and complete
binary gluing replay. It writes RFCXS and RFCXA examples in its output directory:

```sh
cmake --build build --target gadget_native_data_tests
build/gadget_native_data_tests "native example"
```

`NativeStudy::from_history(history, sources, cutoff)` selects an immutable
publication from an admitted series. Omitting `cutoff`, or supplying the final
sequence, reuses the actual retained tail. Earlier sequences use the SDK replay
service. The selected endpoint survives replacement or extension of the caller's
history. Accession indexes are prepared once when requested.

`NativeStudy::read(path, sources, namespace, cutoff, limits)` decodes the complete
archive and applies the same selection rules. RFCXS retains its declared
namespace; a supplied namespace must match. RFCX requires an explicit namespace.
The caller's byte and codec resource limits are enforced before admission.

The [native load profiler](../../tools/native-load/README.md) separates decoding
and accession costs and reports retained component capacity.

`NativeRecord` compiles application values into canonical composite declarations
on an actual endpoint. Objects use Products, alternatives retain their chosen
constructor, arrays declare their admitted cardinality, and exact numbers use
Integer or Rational. Cell and occurrence annotations bind actual native owners.
Text addresses retained bytes; opaque binary uses the native Integer codec and
an explicit byte length. Independent byte sources remain independent even when
their content is equal. No application wire codec is introduced.

`create(value, namespace, sources, endpoint)` retains the supplied complex,
fields, declarations and programs. `derive` preserves that publication and its
original byte owners while adding a new root. `retain_operation` records an
actual source publication, complete retained request and typed result. A record
created without a scientific endpoint is an administrative declaration and does
not assert biological topology.

`create_from_endpoint(value, namespace, sources, std::move(endpoint))` consumes
temporary endpoint storage before canonical freezing. Records retain their SDK
publication, so writing again reuses its encoded components. `RawByteSource::rekey`
assigns another declared identity to the same frozen byte owner. Transferred
vectors with possible mutable aliases are copied before that owner is shared.

For a large result collection, `NativeRecordValue::PackedArray` accepts a count
and a synchronous producer. The producer emits complete values in order:

```cpp
using Value = gadget::NativeRecordValue;
Value answers(Value::PackedArray{requests.size(), [&](const auto &emit) {
    for (const auto &request : requests)
        emit(compute_complete_value(request));
}});
auto publication = source.derive(std::move(answers));
publication.write("answers.rfcxs");
```

`source`, `requests` and `compute_complete_value` are supplied by the caller.
Each emitted value is encoded with RexGraphCpp's `CompiledLayout` and
`encode_declared_bits` before the next value is requested. The publication
retains the complete declaration registry, encoded bytes and ordered addresses.
It also shares equal generated text spans within the new byte pool. Explicit
original source identities and their bytes remain distinct.

Reading validates every packed value with `decode_declared_bits`, including
values retained by later publications. Selecting an element decodes that
element; member views retain its decoded owner after the parent view is
destroyed. Keep the selected `NativeRecord` while using references returned by
`value()` or `tensor()`; that view owns the decoded storage. Producer exceptions,
incorrect counts, invalid declarations, byte
addresses, bit lengths, padding and cyclic references are refused. The producer
must complete synchronously and must not retain the emitter.

Packed collections use GADGET application record revision 4 inside RFCXS.
Ordinary records continue to use revision 3. This reader accepts revisions 2,
3 and 4; older GADGET readers cannot admit revision 4. The RexGraphCpp RFCXS and
composite binary formats are unchanged. Packed collections reduce expanded
application value storage; the actual scientific complex, fields, source bytes
and complete results still occupy memory and storage.

New compact biomedical annotations pack their complete top level collections,
including fields, contexts and document bindings. The declared producer history
and its measured tensors are retained unchanged. Ordinary and earlier study
annotations remain readable. `packed_encoding(index)` exposes an array member's
actual native type, bit count and encoded byte span; retain the parent while
using that span. Declaration digests bind packed schemas and complete encoded
payloads together with their referenced original bytes. Digest traversal selects
one member at a time and preserves factor and collection order.

The consuming `derive` and `create_from_endpoint` overloads and
`retain_operation` accept an optional `NativeRecordObserver`. It reports
completed source declaration, composite value, declaration channel, canonical
freeze and record binding stages. An observer exception aborts construction.

Biomedical stores, scientific studies, assay workspaces, snapshots, workflow
histories and pipeline artifacts use these native owners. A snapshot directory
contains `study.rfcxs` by default. Assay workspaces publish `workspace.rfcxs`;
pipeline and partition receipts refer to `.rfcxs` artifacts. They can be replayed
without JSON sidecars or the original source filesystem. Hashes check integrity;
the actual retained bytes and declarations supply provenance.

The CLI can admit an existing external record and print a temporary view:

```sh
gadget native-record --input request.json --out request.rfcxs
gadget record-view --input request.rfcxs --format pretty
gadget record-view --input pipeline.rfcxs --member result
```

JSON remains an input adapter or exported view. Mathematical operations use the
original complex and its native fields, with the same authored carrier, support,
metric and arithmetic. A model receives admitted certificate data and proposes
queries; its proposal must pass ordinary source and method admission before any
mathematical operation runs. The model does not own or rewrite the source.
