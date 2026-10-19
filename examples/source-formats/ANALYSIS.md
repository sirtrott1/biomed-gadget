# Native measurement analysis

These examples connect original source files to GADGET's existing exact analyses.
Measurements use native fields. Original files remain separate provenance.

## Import

Build `gadget_source` with `GADGET_BUILD_EXAMPLES=ON` and install the
[Python data package](../../tools/rcdb/README.md). The reader request declares the
format and resource limits.

```sh
python examples/source-formats/import_native.py gtex counts.gct counts.native \
  --reader build/presets/cpu/examples/source-formats/gadget_source \
  --request intake.json --name counts
```

A GTEx request declares `layout`, `feature`, `measurement`, `sample_attributes`,
`subject_phenotypes` and `limits`. Layouts are `gct12` and `portal_tissue13`.
Measurements are `integer_counts`, `estimated_counts` and `tpm`. Reader limits
name `source_bytes`, `rows`, `columns`, `record_bytes`, `field_bytes` and `cells`.
Sources must be decompressed files with their original contents.

The `bindingdb` profile reads bulk TSV. Its limits name `source_bytes`, `rows`,
`columns` and `field_bytes`. `measurements` contains reported nM thresholds;
`comparison` retains equality, inequalities and approximation separately.
Compound identifiers, construct occurrences and all target chain fields remain
context. Context membership does not assign an interaction weight.

Both imports write bounded native records as they read. Late source errors
cancel publication. Original coordinate order, repeated occurrences, measured
zero and missing support remain distinct. Reader and record budgets are execution
admissions, not process memory guarantees.

## Analyze

[analyze_counts.py](analyze_counts.py) selects native rows and columns, binds an
observation chart and calls an existing scientific operation. Its binding names
`id`, `design`, `rows`, `columns`, `responses` and `covariates`.

A response is `{"id": "gene_a", "row": 0}`. The row selects the exported matrix;
its original coordinate and label remain available. Covariates declare their
`id`, `quantity`, `units`, `value_type` and `values` in the selected observation
order. Exact rational strings and null preserve values and missing support.

```sh
python examples/source-formats/analyze_counts.py counts.native analysis \
  --binding binding.json --query query.json --namespace 700100 \
  --binary build/presets/cpu/gadget \
  --importer build/presets/cpu/gadget_native_field_import
```

The chart has an explicit observation carrier. Its original source complex,
measurement tensor and coordinate maps remain retained. Repeated donors remain
available; statistical independence is an analysis assumption.

Use `linear_model_comparison` for declared product terms, complete coefficient
families and identifiable contrasts. Use `observed_response_model` for declared
count and exposure fields, or `constrained_response_model` for admissible response
constraints. Normalization, metrics and reference laws require declarations.
These are exact field analyses. They do not fit an NB2 likelihood or create Wald
intervals. Sampling confidence requires a justified reference law.

## Pathways and evidence

`gadget analyze-record` consumes a native result with a
[result operation](../../schemas/result-operation.schema.json). A pipeline
`result` step uses `result_input` to select an earlier native result. Source and
request hashes remain bound to the operation.

`pathway_overrepresentation` retains released GMT memberships, an explicit
correspondence table, original feature occurrences, background and tested family.
It constructs the exact finite intersection law for a uniform subset conditional
on selection size. Existing finite calibration computes inclusive tails.
Bonferroni, BH and BY policies use exact rational arithmetic. Their statistical
assumptions remain explicit; overlapping sets do not imply an independent joint
law. See the [request schema](../../schemas/pathway-overrepresentation.schema.json).

`feature_evidence_context` joins explicitly mapped features to actual assertion
participant occurrences. Source roles, repeated occurrences and unresolved or
ambiguous mappings remain available. A label match alone establishes no biological
correspondence. See the [request schema](../../schemas/feature-evidence-context.schema.json).

## PubChem

[pubchem_evidence.hpp](pubchem_evidence.hpp) admits typed readouts through the
canonical native interpretation. Definitions, outcomes, readouts, revocations,
SID and CID references remain separate. Definitions are referenced by identity.
Missing readouts do not copy the entire substance record.

Construct `EvidenceStore(publication, pubchem_evidence_adapter())` to use the
common query and context service. Admission and reopening check the original
source and typed interpretation. One declared assay is decoded in memory; its
limits must be set before admission.

The optional `gadget_source` executable also accepts `pubchem-query`,
`pubchem-result` and `pubchem-pipeline`. Each command selects the native PubChem
adapter and calls the common service. A query takes `EVIDENCE REQUEST OUTPUT`;
a result adds `ANALYSIS` before the output; a pipeline adds `WORKSPACE` instead.
The main `gadget pipeline` command accepts `--evidence` for its standard evidence
store.

## Tests

The standard C++ contracts include source readers, native evidence, exact pathways
and result pipelines. Enable `GADGET_TEST_NATIVE_PYTHON=ON` at configure time to
include measurement intake and count analysis through Core, RCDB, RCQL and the
bridge. The selected Python environment must contain those packages.
