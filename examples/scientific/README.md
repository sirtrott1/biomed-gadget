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
