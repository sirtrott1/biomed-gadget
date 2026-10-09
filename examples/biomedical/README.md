# Binding evidence and local workflows

This example imports declared ChEMBL and BindingDB responses, retains the
original assertions and evaluates queries on the resulting evidence store.
Build the database profile for bounded HTTPS acquisition:

```sh
cmake --preset databases
cmake --build --preset databases
```

Native import and query also work without network transport or a local model.
The example retains an explicit ChEMBL page and BindingDB affinity cutoff; it is
not a complete compound/target survey.

```bash
build/presets/databases/gadget databases
build/presets/databases/gadget fetch-databases \
  --request examples/biomedical/imatinib-binding-fetch.json \
  --out data/public/imatinib-binding-v1
build/presets/databases/gadget import-evidence \
  --snapshot data/public/imatinib-binding-v1 \
  --out data/evidence/imatinib-binding-v2
build/presets/databases/gadget evidence-query \
  --evidence data/evidence/imatinib-binding-v2 --query YOUR_QUERY.json
build/presets/databases/gadget query \
  --study data/evidence/imatinib-binding-v2/study/study.json --query YOUR_QUERY.json
```

Evidence queries: `inspect`, `coverage`, `affinity_evidence` (per record bound/point entailment for an explicit endpoint/unit maximum), `select` (exact endpoint/relation/unit/provider
category), `neighbors` (native subject key), `compare` (two assay identified equality
records with identical compound/target/endpoint/units). Inspect exposes the native
subject keys. Selection/traversal returns source pointers, receipts and original
records. Comparison describes reported tokens; shared assay IDs do not prove
experimental comparability or independence.

The scientific study has `all.Kd_nM`, `all.Ki_nM`, `all.IC50_nM`, `all.EC50_nM`
and declared source categories. Only reported equality points in nM inhabit these
numeric fields. `field_summary`, `support_mask` and `native_boundary` use the existing
scientific query schema. Missing, censored and other unit values remain in evidence.
The `record-map.json` correspondence links the mathematical record carrier to
native assertion keys. Entity identity is not inferred from coordinate dimensions.

Workflows may include `evidence` nodes alongside `query`/`model` nodes. Bind both
`source_digest` and `evidence_digest`, and pass `--evidence DIRECTORY`. A local role
can select `"task":"observer"` to read supplied query results or images
without receiving the entire numerical planning catalogue. Its statements remain
prospective. Planner roles retain strict mathematical proposal admission.
`--task observer` also works for local `ask` and `session` commands.

DrugBank needs licensed access; it is not included
in the open data connector set.

Prepare and run an entirely offline, source bound C++ workflow after import:

```bash
python tools/prepare_biomedical_workflow.py --binary /absolute/path/to/gadget \
  --evidence data/evidence/imatinib-binding-v2 \
  --out verification/binding-workflow --run
```

Add `--model /absolute/path/to/LOCAL.gguf` to register an observer and invoke it
conditionally after the exact evidence criterion. This helper writes declarations;
the C++ executable admits and executes every evidence and mathematical operation.
The example runs four worker canonical field actions and repeats the same action
through `cpu_cached_tiles` to verify exact result reuse. Local observers are
not qualified to reinterpret censored binding measurements; use the exact
certificate for numerical decisions and retain model statements as prospective.
Imports record their adapter version. Version two imports add standardized InChIKey
assertions; `--adapter-version 1` reproduces the earlier projection. An older store
is never silently extended when read.

For a numerical reading that GADGET checks, select `--reading reported_points`
with `--model`, `--context-tokens`, `--output-tokens` and `--threads`. The model
copies four declared field summary values; the next workflow node checks each
comparison against the original certificate. `--reading coverage` supplies the
complete coverage query result to a prospective observer. The default
`affinity_records` reading supplies every record of the threshold result and
needs a context large enough for that complete response. No reading silently
truncates a certificate.

The helper saves `operation.rfcxs` as well as the displayed receipt. When a
workflow uses an evidence store, the native operation retains the study and the
actual separate native evidence publication, including the original provider
response and header bytes. The source keys of those two publications do not
identify their entities or carriers with each other.

`EvidenceStore::source_context` accepts a native control request with schema
`gadget.biomedical.context.v1`, the evidence `source_digest`, the retained
`snapshot_sha256`, a nonempty `question`, and an ordered `assertion_keys` selection.
It retains every selected original assertion on its native owner. The catalog
binds the question and request with digests. Its availability profile is a frozen
snapshot: retrieval does not establish historical availability or biological time.

`EvidenceStore::admit_proposal` checks schema `gadget.biomedical.proposal.v1`
against that same request and question. A query is prospective and requires
separate execution. Comparisons and subject selections bind selected assertion
roles; categorical and affinity predicates require selected source evidence.
Category and affinity queries execute over the complete frozen snapshot when
subsequently passed to `execute`, rather than only the context selection.

Entity annotations name an original assertion key, role and exact source
identifier. Measurement annotations preserve its endpoint, units, relation and
original value token, including a censored bound or missing token. Every annotation
quotes exact UTF8 byte spans from the assertion's original response, using its
request ID as the document ID. Admission verifies these bindings and literal
quotes; the semantic interpretation of a quote remains unverified. Receipts retain
the context, proposal, original response and header byte owners for offline replay.
`proposals.cpp` demonstrates this API with native request and receipt files.
The request digest also binds its native declaration identity. Save the context
request as a native record and reuse that retained request when admitting a
proposal; importing the same text into a new native request creates a separate
request owner. The question digest separately binds the literal question.

`EvidenceStore::execute_record` takes a native query request and returns the
original source complex, selected owned assertions and exact native count or
rational results. Interval classification and point comparisons read the existing
native provider values directly. `execute(Json)` renders a compatibility
certificate after that calculation.

With the optional local model adapter, `LocalAgent` can bind an immutable
`EvidenceStore` and its retained native context request. `ask` accepts that
request's exact question and returns a prospective query, annotation or
clarification. `prospective_publication` retains the latest admitted native
proposal. The caller reviews the selected query and separately passes its native
`request.query` to `execute_record`; generation does not execute it. Retained
proposals can be replayed without another model call.

`EvidenceStore::bind_reported_field(study, request)` checks the existing reported
point projection before selecting a quantitative field. Its native request uses
schema `gadget.biomedical.field.request.v1`, the evidence `source_digest`, the
`study_digest`, a declared `field`, its explicit `endpoint` and `units: "nM"`.
The current projection includes the complete assay identified ChEMBL record
carrier. BindingDB compact assertions remain in the separate evidence owner.
The check verifies the retained snapshot and ordered record correspondence,
quantitative declaration, and every native value and support coordinate against
the original interpreted assertion payload. Censored bounds, missing tokens,
other units and other endpoints have absent support; reported zero has support.

The native binding retains the complete study archive and a separate archive at
the field's actual publication cutoff. Its `selector` uses the existing
`gadget.rcdb.field.query.v1` contract: exact source namespace, cutoff, tensor,
input/output carriers and interpretation. The evidence complex remains a
separate owner. An assertion's composite payload is not a quantitative tensor.
This verifies the source projection, without asserting assay comparability,
pooled affinity or a clinical interpretation.

The `fields.cpp` consumer writes `binding.rfcxs`, `field-source.rfcxs` and
`selector.rfcxs` into a new directory. With that example and the optional native
field selector built, run:

```bash
gadget_biomedical_fields EVIDENCE_DIRECTORY STUDY_DECLARATION REQUEST NEW_DIRECTORY
gadget_native_field_query --input NEW_DIRECTORY/field-source.rfcxs \
  --request NEW_DIRECTORY/selector.rfcxs --out NEW_DIRECTORY/selected.rfcxs
```

The selector receipt's `result.selected_publication` owns the unchanged native
field source. Reopen that archive for native calculations. The surrounding
selection receipt is an administrative publication; its namespace does not own
the selected field. The same selector request can be retained in RCDB and used
with the registered `GADGET_NATIVE_FIELD` RCQL operator; version selection and
the byte provider profile remain explicit. See `tools/rcdb/README.md`.

See [Biomedical response providers](providers.md) for the local original response
intake and provider extension contract, including source bytes, units, terms and
frozen replay requirements.
