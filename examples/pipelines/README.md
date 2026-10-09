# Assay pipelines

Choose a profile with `gadget pipelines`: `rna_seq`, `single_cell`, `atac_seq`,
`chip_seq`, `spatial`, `multi_omics`, `variants`, `sequence` or `custom`.
A profile lists operations appropriate to that assay family. It does not execute
all of them. The question, available measurements and declared parameters determine
the steps in a pipeline.

Defaults are selected recipes with explicit parameters. `gadget pipelines` lists
the recipes admitted by each profile and their `native_publication` capability.
Every selected query and derived recipe must have an admitted native producer.
All producer checks run before pipeline jobs start, including the scientific
operation inside a `native_query`. An unsupported producer refuses before a
destination is created. Each step retains its producer's native result owner;
the optional JSON certificate is a view of that result.

A goal can be the existing `quality` or
`coverage` string, or an object selecting a method and its arguments. The profile
restricts the available operations; only the selected goals execute.

The `indexed_sequence_alignment` recipe accepts `parameters.execution` with
positive `read_partition_records` and `reference_block_bases`. These widths
control normalized read batches and reference blocks. They preserve complete
placement results across block seams and do not change the edit budget or source
selection. The pipeline retains the supplied widths in its original native
request; the assay service admits and executes them.

| Profile | Entry data and recipe choices |
|---|---|
| RNA seq | FASTQ quality, alignment evidence, annotation and feature assignment bounds; admitted counts support matched expression contrasts and declared controls |
| Single cell | An admitted count matrix and cell dictionary; declared aggregation constructs pseudobulk, followed by matched contrasts or control analysis |
| ATAC seq and ChIP seq | Declared alignment/interval resources and counting units; interval coverage, matched scaled chromatin contrasts and motif occurrences |
| Spatial | Counts with an observation carrier, coordinates, frame and units; physical relations can feed native Green, Hodge or local section queries |
| Multi omics | Independently admitted assays with an explicit cross assay map; field transfer, aggregation and native field analysis |
| Variants | Variant records, reference and filter policy; genotype materialization and native field queries |
| Sequence | Read quality, motifs, exact pair/placement/edit/splice queries, observed k mer complexes and read overlap reconstruction; derived complexes support native structural analysis |

A recipe's `parameters` are the arguments of its existing operation, without
`schema` or `operation`. The operation contracts are
[assay queries](../../schemas/assays-query.schema.json) and
[scientific queries](../../schemas/scientific-query.schema.json). The
[default schema](../../schemas/pipeline-default.schema.json) references those
contracts. Native recipes select `assay` on the goal, or use the default's assay;
if the input contains exactly one assay, that identity is selected unambiguously.
Workspace recipes put their assay/resource identities inside `parameters`.

For pseudobulk, aggregation declares `count_policy: "disjoint_count_sum"`.
Admission requires unchanged feature identities with unit coefficients, disjoint
nonempty cell groups with unit coefficients, and nonnegative integer source
counts declared as `read_count` or `umi_count`, in their original channel and
units. The native map action then retains the count quantity, so its output can
feed matched expression analysis. General linear aggregation retains its separate
derived quantity; fractional weighting or overlapping groups fail the count sum
contract rather than acquiring count semantics.

A materializable goal may declare `output` with an assay identity, namespace,
design and channel. A later goal's `input` names that producing goal. Its native
assay selection then follows the actual derived assay identity. This connects
assignment bounds to field inspection, pseudobulk to contrasts, interval contrasts to
field analysis, and physical relations to complete native field actions.
Independent expression assignments, feature counts, controlled coverage,
spatial relations, aggregation and design aggregation materialize their native fields directly. Coupled
expression constraints use the existing constrained result path.
The native work result retains each admitted derived workspace. Dependent branches
and binding steps share that owner and verify its digest, so they do not decode the
same workspace again. A separate pipeline run admits its own inputs.
The result releases its workspace once every declared consumer has acquired it;
each active job retains its own reference until completion.

The supplied templates are [paired expression](paired-expression.json),
[pseudobulk with a paired contrast](pseudobulk.json),
[matched chromatin coverage](chromatin.json), and
[spatial Green response](spatial.json). Replace the declared identities,
experimental pairings and protocol parameters with those of the admitted data.
These templates enter at existing counts, alignments or intervals; protocol specific
barcode correction, genome alignment and other preprocessing require their own
qualified methods.

[Sequence structure](sequence-structure.json) connects observed k mer occurrences
to their native relational complex, a complete chain certificate and a Green
response. [Spatial structure](spatial-structure.json) retains a measured field on
declared physical incidence and reports its complete structural decomposition.
[Cross assay comparison](multi-omics.json) supplies a quantitative map between
independently measured carriers. Replace its example coefficients and quantities
with a justified measurement map; an identifier link alone does not supply one.

Selectable recipes include the existing `chain_certificate`,
`path_comparison`, `typed_moment`, `temporal_moment_change`, `structure_transition`
and `green_response` operations. `relational_persistence` consumes complete
source declared births and levels. Sequence and constrained allocation recipes
use their existing finite operation profiles. Each recipe refers to its original
operation schema, and materializable sequence/allocation goals support the same
derived input dependency mechanism as counts and spatial relations. These choices
also work in `plan-recipes`; the model selects declared names and their required
inputs without supplying coefficients, biological relations or filtration levels.
The [sequence intent](sequence-structure-intent.json) and
[cross assay intent](multi-omics-intent.json) supply those selectable bindings.

## Declared aggregation and cross assay maps

Aggregation applies the supplied feature and observation maps as `A X B^T`.
RexGraph owns both exact actions. The numerical intermediate includes only source
coordinates with nonzero map dependencies; original files, full axes, all declared
coefficients and output order remain in the producing certificate and lineage.
`parameters.workers` selects native CPU workers, with joined feature and observation
stages. Execution receipts report both stages and the required source dimensions.
One native CPU reservation covers source loading, map preparation and both joined
stages; the intermediate stages do not repeat hardware admission.
The pipeline infers that allocation for query and materialization goals unless
the goal supplies `workers`. An explicit allocation must accommodate the request.

Every output is supported only when all its nonzero dependencies are supported at
both stages. Cancellation does not turn an unknown input into an observed zero.
A declared zero map has no nonzero dependencies and returns a supported zero.
Disjoint pseudobulk count sums keep their separately checked count contract.
Group, map entry and intermediate dimensions follow the actual request instead
of the former fixed aggregation ceilings. Source formats retain their separate
admission profiles, and the complete output still costs space proportional to its
declared dimensions.

`field_transfer` retains complete numeric or indicator fields on independent
native SourceGrade carriers. Its explicitly supplied map acts from the right
carrier to the left; every mapped coordinate and residual remains available with
pairing, residual quadrance, rank and kernel dimension. RexGraph computes rank
once; rank nullity gives the dimension without constructing an unused kernel
basis. This method has no separate 256 coordinate ceiling. Complete support is
required on both supplied fields, including coordinates outside the map's image.

These observation maps do not automatically transfer a boundary or topology.
Accession, physical realization, reconstruction and verified chain maps have
separate declarations. Structural recipes consume the actual
source complex; aggregation cannot assert a quotient complex merely by grouping
its measurements.

For example, this default runs only measured support inspection:

```json
{
  "schema": "gadget.pipeline.default.v1",
  "profile": "rna_seq",
  "question": "Which measurements are reported in this expression assay?",
  "assay": "expression",
  "fields": ["all.sample_A"],
  "goals": ["coverage"]
}
```

Add goal `quality`, `reads` naming an admitted FASTQ resource and explicit
`quality_thresholds` to include a quality check. Prepare or execute the defaults:

```sh
gadget pipeline-prepare --workspace workspace/workspace.json --request defaults.json
gadget pipeline --workspace workspace/workspace.json --request defaults.json --out results/support
```

Custom designs use `gadget.pipeline.design.v1`. Every step names an `id`, `type`,
`input` and `request`. Types are `query`, `derive` and `workflow`. Input is `source`
or an earlier `derive` step. `depends_on` adds ordering dependencies. `workers`
declares a step's CPU demand; `concurrent_jobs` controls the pipeline's maximum
number of simultaneous jobs. RexGraph owns admission and joined completion.

For a partitioned `native_query`, the step's `workers` is its total CPU allocation.
The inner query's `workers` declares workers per partition, and the request's
`concurrent_jobs` requests parallel partitions. For example, 32 assigned workers,
four workers per query and eight concurrent partitions let the native work graph
schedule eight independent jobs within that allocation. Effective concurrency is
the minimum of the requested jobs, selected partitions and assigned capacity
divided by workers per partition. This replaces the former four job GADGET cap.
Whole carrier mathematics retains its existing partition checks: storage rows
do not define independent Green, homology or sheaf problems.

## Reusing field actions and selecting outputs

[Controlled expression](controlled-expression.json) shows a whole carrier
analysis with explicit controls, endpoint and four allocated CPU workers. Its
`parameters.workers` requests native field workers; the goal's `workers` supplies
that allocation. The fields must share one source carrier, with declared numeric
or indicator meanings. Categorical codes require separately declared indicators.

For an independent frame F of the declared controls (including the intercept),
GADGET composes RexGraph's existing maps, metric Gram and prepared exact solve to
evaluate `F solve(F*MF, F*Mx)`. The factorization is prepared once per request.
Each field retains its full input, represented part and residual. Reconstruction,
orthogonality and the solve residual are checked with rational arithmetic. Matched
contrast requests share the corresponding constant span preparation. Execution
receipts name the native owner, preparation count, applications and worker slots.

`gram_fields` selects an ordered subset of tested fields for the residual Gram
block. An empty array returns a 0 by 0 block; omitting it preserves the full tested
family. The `gram_fields` result names the block's columns. This selects an output,
without reducing or changing any field or the endpoint calculation. A full Gram
has quadratic output size in the number of fields, so requesting only the needed
block can avoid unnecessary output and arithmetic.

These two methods use actual carrier and field dimensions instead of the former
4,096 record and small field list ceilings. Import formats and the trajectory
method retain their separate admitted profiles. Exact bit growth, fill and full
output still contribute to cost. Certificate sealing accepts complete outputs through the JSON codec's integer
representation rather than imposing the former 16 MiB ceiling. Default CPU
admission uses observed affinity
and quota; an embedding application can supply its own native resource policy.
Automatic field calibration considers geometric worker counts through the actual
CPU and tile capacity, including that capacity.

Partition jobs write their complete certificates inside joined native jobs.
Only file receipts remain in the work graph, and the wave size follows effective
concurrency. The index records hashes, source identities and allocations. A failed
run removes its unpublished results. Whole carrier analyses still require one
carrier; writing separate storage partitions does not make their control, Green,
homology or sheaf problems independent.

## Planning from a research question

`plan-recipes` lets a local language model select and connect individual recipes
from a question. Supply an intent declaration with reusable bindings: each has
an `id`, a `purpose` describing the analysis it answers, and a fully parameterized
`goal` following the default schema. The binding supplies the goal identity.
Examples are [pseudobulk and expression](pseudobulk-intent.json) and
[spatial field response](spatial-intent.json).

```sh
gadget recipe-catalog --workspace workspace/workspace.json --intent intent.json
gadget plan-recipes --workspace workspace/workspace.json --intent intent.json \
  --model /path/to/local.gguf \
  --prompt 'Compare the declared matched units and retain their feature differences' \
  --threads 4 --context-tokens 8192 --output-tokens 1024 --out results/question
```

The model returns named goals. GADGET includes their declared input and ordering
prerequisites, then compiles them through `pipeline-prepare`'s existing default
compiler. A selected contrast can require pseudobulk; a selected Green action can
require physical relations. Unselected branches do not execute. Declare producers
before their consumers. The same declaration can answer several questions without
running every available analysis.

The catalog presents every binding's purpose, method, prerequisite identities and
goal hash. Parameters remain in the complete retained declaration. The model cannot
rewrite them. Missing pairings, metrics, coordinates or other experimental design
require an updated declaration; the planner returns clarification. An intent with
an empty `recipes` array supports this case without execution. No fixed binding or
graph depth ceiling is imposed; an oversized prompt fails the requested model
context rather than silently dropping bindings.

The proposal binds the workspace digest and intent digest. `recipe-proposal`
admits a proposal file through the same path without a model:

```sh
gadget recipe-proposal --workspace workspace/workspace.json --intent intent.json \
  --request proposal.json --prompt 'The original research question' --out results/replay
```

Successful execution retains `planning.json`, its hash in `pipeline.json`, the
original question, full intent declaration, proposal, included prerequisites and,
for local planning, model identity and timings. Mathematical results retain their
existing source bound requests and certificates. Planning decides which declared
analyses answer the question; numerical admission and calculation remain with the
assay services and RexGraph. `plan-pipeline` retains its separate query step
interface for direct inspection of existing sources.

```json
{
  "schema": "gadget.pipeline.design.v1",
  "profile": "rna_seq",
  "question": "Inspect read quality and the support of supplied count fields",
  "concurrent_jobs": 2,
  "steps": [
    {
      "id": "quality", "type": "query", "input": "source",
      "request": {
        "schema": "gadget.assays.query.v1", "operation": "read_quality",
        "resource": "reads", "quality_thresholds": [20, 30]
      }
    },
    {
      "id": "support", "type": "query", "input": "source",
      "request": {
        "schema": "gadget.assays.query.v1", "operation": "native_query",
        "assay": "expression", "query": {"operation": "coverage", "field": "all.sample_A"}
      }
    }
  ]
}
```

An existing materialization request can be a `derive` step. A following query can
then name the producing step as its input and select the derived assay. A
`workflow` step additionally selects `assay` and optional `partition`; its request
uses the scientific workflow schema. GADGET binds that workflow to the actual
selected study. The public C++ runner accepts existing model callbacks and an
optional evidence store for branching, recursive scientific workflows.

Each step retains its resolved `request.json`, result and hashes. Native
pointwise queries keep their per partition certificates. The final index contains
the source bound plan and receipts. Failure removes the staged pipeline; it does
not publish a completed result directory. To rerun, supply the retained input
workspace and saved plan with its matching source digest.

An expression recipe may supply normalization evidence directly, or explicitly
select `normalization_evidence_field`. The latter copies that field's admitted
source spans from the actual input assay into the resolved contrast request. This
allows a dependent contrast to reference a newly materialized assay. Those spans
identify the measurements used; the plan separately records the researcher's
normalization method and exact scales. The binding neither estimates scales nor
asserts that the source publication prescribed them. Supplied evidence and a
field binding are mutually exclusive.

Pipeline concurrency defaults to the selected step count and the available
native CPU allocation. Worker demand is admitted against that allocation.
Expression contrasts have no GADGET matched unit or selected field count ceiling;
all requested units must exist in the admitted carrier. Input format, source profile
and arithmetic budgets remain separately declared by the existing services.

## Mathematical definitions

The recipes consume RexGraphCpp field, metric, local section and solve APIs.

Matched expression retains each declared unit's scaled difference. With constant
frame `F` and the declared matched unit identity metric, RexGraph applies
`F solve(F* M F, F* M delta)` and returns both consensus and residual; their sum
reconstructs the complete difference vector and their metric pairing vanishes.
Control analysis uses the declared control span and the corresponding native
metric projector. Aggregation evaluates declared exact linear maps while
requiring support for each nonzero input dependency. Aggregated fields are derived
observations; they do not replace the source fields or define an inferred topology.

Physical neighborhood recipes append declared grade one incidences to retained
relations. Higher assertions require their own source declarations. Green response
uses the actual native boundary, adjoints and metric actions; its full rational
output satisfies `(I + lambda L_k)x = j` with a checked zero residual. Local section
analysis requires explicit stalks and restrictions. GADGET assembles these requests,
records experimental meaning and provenance, and schedules native jobs; it contains
no replacement projector, Green solver or sheaf implementation.

## Local question planning

```sh
gadget plan-pipeline --workspace workspace/workspace.json \
  --model /path/to/local.gguf \
  --profile rna_seq \
  --prompt 'Inspect the reported support of the expression assay' \
  --threads 4 --context-tokens 8192 --output-tokens 1024 \
  --out results/question
```

The model receives the admitted workspace catalog and profile choices. The
`--profile` option fixes the researcher's selected assay profile in the decoding
schema; without it, the model proposes a profile. Its structured proposal declares
selected requests; GADGET checks source binding,
profile, identities, dependencies and request contracts. An admitted proposal
runs when `--out` is supplied. A clarification produces no execution. Omitting
`--out` returns the proposal and prepared pipeline for inspection. This is a
planning interface, not an assertion of general biomedical interpretation quality.

Other model ecosystems can supply `gadget.pipeline.proposal.v1` through
`gadget pipeline-proposal --workspace ... --request proposal.json [--out ...]`.
Each proposal step supplies a structured `request`; encoded `request_json` is also admitted for existing external producers. The local decoder currently emits query steps on the source workspace. Derived workspace and recursive workflows use explicit designs or external structured producers. Programmatic designs
use ordinary JSON request objects. Assay specific schemas constrain the selected
profile; the operation's own admission remains authoritative.

## Native producer ownership

`AssayWorkspace::materialize_workspace` returns the admitted immutable workspace
after native publication. Coverage, spatial, independent assignment and aggregation producers use
this interface. `supports_native_materialization` reports the available producers;
unsupported operations refuse before creating a destination.

`AssayWorkspace::bind_workspaces` returns a bound workspace that retains the selected
studies. Aliases select those studies without changing their native namespaces,
source histories or field accessions. Dependent pipeline jobs consume these owners
through the existing producer leases. Native derivation results retain the actual
producer publication, and binding results retain the complete collection publication.
The existing materialization and binding receipt views remain available.
Aggregation retains the exact feature and observation maps, ordered group declarations,
output support and original source publications. Design aggregation also retains
the categorical row selections and original dictionary evidence. Derived feature
identifiers remain available to downstream queries. A later aggregation can consume
the retained fields directly; it does not read a generated matrix file. These maps
declare quantitative observations and do not transport the source topology.
The producer's `feature_map` and `observation_map` values identify native observation
tensors by key. Their input and output carriers declare the complete source and
group axes. The native tensor encoding retains full support and sparse coefficients.
In a pipeline, the native derivation `result.rfcxs` uses `gadget.operation.native.v1`;
the binding `result.rfcxs` uses `gadget.assay.workspace.native.v2`.

## Exact field actions

The native scientific query `green_action` evaluates `(I + lambda L_k)x = j`
through RexGraph's boundary and metric actions. It returns the complete response,
exact iteration count and a checked zero residual. No dense Hodge matrix or inverse
is constructed. `lambda` is a positive rational token.

```json
{"operation":"green_action","context":"all","field":"all.response","lambda":"1"}
```

The C++ `NativeStudy::green_action` interface has no application coordinate count
ceiling. Exact execution cost depends on retained incidences, the number of
distinct operator directions encountered and rational coefficient growth. The
JSON and assay source profiles retain their separate admission constraints.
The older GET `green` query now defaults to this native exact action without its
former 64 coordinate cap. An explicit `mode: "numerical"` retains the historical
floating point compatibility path.

Aggregation evaluates the feature map or observation map first according to the
smaller actual intermediate. Both orders use native exact matrix actions and
joined CPU tiles. Execution receipts report `first_map` and
`intermediate_coordinates`; source support and requested output order are preserved.
