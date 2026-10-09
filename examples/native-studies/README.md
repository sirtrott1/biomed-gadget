# Native study declarations

`ScientificStudy::publish_annotations(producer, study, documents)` publishes
biomedical declarations after an existing `NativeRecord` history. Admission with
`ScientificStudy(publication)` validates the complete producer prefix, original
bytes, evidence bindings, field meanings and metrics. The annotation preserves
the source namespace, cells, ordered occurrences, carriers and numerical fields.

The `study` value uses `gadget.scientific.study.v3` and contains `id`, `design`,
`documents`, `cell_bindings`, `occurrence_evidence`, `contexts`, `fields` and
`field_sets`. The API supplies the actual `native` prefix binding; callers must
not supply a replacement. Construct these values with `NativeRecordValue`.
`design` declares the question, population, selection, covariates and limitations.
Fields declare their stored key, context, channel, quantity, units, method,
value type and byte evidence. Categorical fields cannot serve as linear amplitudes.

Each document declares its id, path label, SHA-256 and `utf8` or `binary` format.
The separate `documents` argument holds matching ids and complete `ByteSpan`
values from the producer's original raw sources. Admission compares every byte
against that source. Evidence columns select spans within these documents.

`cell_bindings` follows the complete native grade order. Each binding declares
`grade`, `id_prefix`, `label_prefix`, `kind` and an evidence rule. Cell ids append
their actual stable key to the prefix. Prefixes must not overlap and end with a colon. Pages resolve
these labels from the source; the publication does not store a label per cell.

| Evidence rule | Required binding |
|---|---|
| `spans` | Complete source-grade fields named by `document_index_field`, `start_field`, `end_field`, plus a `documents` dictionary |
| `documents` | One or more original document ids covering the complete bytes |
| `cofaces` | Adjacent `upper_grade`, whose directly evidenced original cofaces supply the spans |

Set `occurrence_evidence` to `parent_cell`. Each occurrence retains its own
identity, order, aspects and coefficient; its evidence resolves from the actual
parent relation. `source_index()` reports complete counts, and `read_source()`
returns caller-sized presentation pages without changing the mathematical carrier.

A context declares its id, strict-subcomplex block or ambient null block,
explicit time or null time, and one metric per grade. Metrics use `identity` or
`diagonal_field` with a stored field key. A diagonal field must have complete
support on that grade and strictly positive exact values. Its values are reused
when preparing field calculations. Temporal declarations retain their axis,
position, unit and population. `field_sets` binds declared stalk carriers, fields
and maps through the native field realization APIs. Optional `operators` binds stored maps or metric forms
to their actual input and output carriers.

`AssayWorkspace::publish_studies(request, id, studies)` creates a native
collection from admitted study owners. Aliases route requests; each study retains
its original history and namespace. The collection includes the complete
original request publication and binds its design through the actual retained
request byte owners. Controlled coverage uses these APIs through
`materialize_record()`. `replay_derivation()` restores the generating owners and
compares complete native source bytes, carrier declarations, fields, support and
original resources. Replay uses the generating request after collection aliasing.
Answer receipts select the generating result within the retained workspace.
`derivation_certificate()` returns its field bindings and scalar checks; it does
not enumerate the numerical source. `record-view --input WORKSPACE --assay ID`
renders this view. Fact readers validate the actual native publication and selection.
Compiler receipts and admitted source digests remain recorded
separately from that mathematical comparison.

Build and run `gadget_native_derivation_tests` for an executable example covering
annotation, evidence pages, field metrics, Hodge calculations, rejected invalid
bindings, native-only relocation, replay and collection aliases.
