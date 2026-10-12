# Scientific query sequences

`gadget_scientific_sequence` executes native scientific requests against one
admitted assay study. One `ScientificSession` retains the source field pins and
prepared context operations across requests. The output retains the source
publication, complete request bytes and every result certificate in RFCXS.

Build with `GADGET_BUILD_EXAMPLES=ON` or `BUILD_TESTING=ON`:

```sh
build/presets/cpu/gadget_scientific_sequence \
  workspace.rfcxs sequence.rfcxs answers.rfcxs
```

The output path must be new. A request uses `gadget.scientific.sequence.v1` and
names an assay. It supplies either an ordered `queries` array with unique `id`
and `query` members, or a `context_template`.

For example, author this request and convert it with `gadget native-record`:

```json
{
  "schema": "gadget.scientific.sequence.v1",
  "assay": "measurements",
  "context_template": {
    "contexts": ["region:1", "region:2"],
    "fields": ["all.expression:A", "all.expression:B"],
    "methods": ["local_sections", "hodge", "green_action"],
    "restriction": "scalar_identity",
    "metric": "identity",
    "green_lambda": "1"
  }
}
```

```sh
gadget native-record --input sequence.json --out sequence.rfcxs
```

These field and context names are examples; use the identities declared by the
study. Omit `contexts` to select all contexts that explicitly declare
`witness_scope: "block"`. Such a context retains exactly the cells and occurrences
of its native strict local block. Earlier context declarations without this
option retain their existing ambient witness behavior.

The template declares a scalar stalk for each selected cell and an identity
restriction for every retained occurrence. It supports grade zero sections,
Hodge actions, metric adjoints and shifted Green actions. For other stalks,
restrictions, grades or metrics, provide complete individual scientific queries
through `queries`. All requests are validated before execution.

The template keeps ordered context, field and method bindings, and constructs
each query when visited. It does not retain the expanded query collection.
After validation, each complete answer is encoded immediately into a native
packed collection. The output still contains every certificate, its original
request and the source publication. Both authored and template requests use the
same scientific session and native result reader.

Append `--profile` to report tab separated phase times and current and peak
resident bytes:

```sh
build/presets/cpu/gadget_scientific_sequence \
  workspace.rfcxs sequence.rfcxs answers.rfcxs --profile > phases.tsv
```

Phases separate workspace admission, query planning, validation, execution with
immediate encoding, declaration assembly, canonical freezing, record binding
validation and file writing. Execution and encoding occur within composite
value construction. `execution_and_encoding` reports that producer's completed
work; the following `composite_values` phase measures its remaining assembly.
Phase intervals are consecutive. `total` covers those same intervals and
precedes object destruction. Use an external process measurement for complete
elapsed time.
Linux memory observations use `/proc/self/status`; other systems report
`unavailable`. Peak memory is cumulative for the process.

Fields may restrict from their declared context to an explicitly selected block
only when the native carrier inclusion is valid. Required missing measurements
are rejected. Supported zeros remain measurements.

For `local_sections`, a section's `fields` array specifies the source bindings in
coordinate order. Its optional `axes` array names those coordinates in the stalk
frame. This permits separate fields to define independent sections on a shared
frame. Without `axes`, the field names are the axis names. Axis count, order,
uniqueness, quantity, units and complete field support are checked.

See the [spatial expression experiment](../bioinformatics/spatial-expression.md)
for a complete source, native chart construction and independent checks of
section projection, Hodge and Green identities.

`local_sections` publishes exact fields in `result.fields` when executed through
`ScientificSession::execute_record`. Each stalk keeps its actual cell identity
and a distinct native carrier. The bindings retain direct sum inclusions,
identity metrics and the source field realization maps. Occurrence maps retain
their original flag owner and direction. Compatible fields, residuals, defects
and complete cohomology frames use native tensor keys; the JSON view is an
optional rendering. Portable replay retains the original study and request.
