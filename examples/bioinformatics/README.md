# Assay workspaces

An assay workspace contains original files, feature and observation dictionaries,
explicit correspondence records, canonical RexGraph tensor packets, and their
hashes. The same admission and query interfaces serve bulk or single-cell
expression, spatial measurements, chromatin regions and other quantitative
assays. Biological interpretation comes from the supplied assay declaration.

Use `gadget bioinformatics` to inspect the compiled formats and method limits.
The [sequencing and assay recipes](suite.md) cover source-backed feature counts,
splice junctions, variants, matched chromatin fields, physical spatial relations,
declared cell/feature aggregation and complete sequence-pair alignment.
H5AD support requires `GADGET_ENABLE_HDF5=ON` and HDF5 1.12 or later. Its independent
format tests additionally use Python packages `h5py` and `numpy`.

## Prepare and admit an assay

The feature and observation dictionaries are unquoted TSV files. Their first
column is `id`; later columns retain specimen, condition, batch, cell, time,
coordinate, feature, annotation or other source metadata. IDs are unique within
an axis. A dense matrix has `id` followed by observation IDs in its header and
feature IDs in the same order as the feature dictionary. Matrix Market uses
feature rows and observation columns. H5AD stores observation rows and feature
columns; GADGET verifies its internal `obs`/`var` IDs against both dictionaries.

```sh
python3 examples/bioinformatics/prepare.py \
  --matrix counts.mtx --features features.tsv --observations observations.tsv \
  --format matrix_market --id experiment --namespace 100000 \
  --kind expression --channel E --quantity read_count --units counts \
  --method 'Integer gene-assigned counts from the recorded counting pipeline' \
  --question 'Which matched expression changes are shared across specimens?' \
  --population 'The specimens in the supplied observation dictionary' \
  --carrier features --implicit-entries measured_zero --arithmetic exact_decimal \
  --out data/import
build/presets/cpu/gadget native-record \
  --input data/import/import.json --out data/import/import.rfcxs
build/presets/cpu/gadget import-assays \
  --spec data/import/import.rfcxs --out data/workspace
```

A sparse omitted value is a measured zero only when the caller declares that
policy. `unobserved` instead retains support only for reported coordinates,
including explicitly reported zeros. A missing-coordinate mask can override
otherwise implicit zeros in Matrix Market; the mask contains one-based row and
column indices separated by a tab. Dense matrices use declared missing tokens.
Duplicate coordinates reject admission rather than silently sum. The quantities `read_count` and `umi_count` require
nonnegative integers; UMI observations retain their original counting units. Other measured quantities can be exact decimal or rational
values. H5AD uses `ieee754_exact`, retaining the rational value of a stored float,
and an explicit `nan_policy`; this does not increase measurement precision.

Choose `carrier: observations` to make gene or protein measurements fields over
cells, spots or specimens. Declare `observation_fields` to expose numeric metadata columns as exact measured
covariates, or a specified category level as a dimensionless membership indicator.
Each declaration names its column, quantity, units, method and missing tokens.
Other category levels are supported zeros; missing metadata stays absent.
This exposes the existing `controlled_pairings`, Hodge,
Green and sheaf methods on those source identities. It retains the complete
measurement coordinates; it does not create a fitted embedding or select cell
labels. Field counts follow the complete source dictionaries and native keys.

Authored relation files admit ordered positive native grades. Participants name
an adjacent lower cell, with explicit existence, a `negative` orientation bit and
nonnegative `magnitude`. Relation arity does not choose its grade. Higher
attachments must satisfy the native chain/Bianchi checks and never arise merely
from overlap. Optional relation fields retain independent observed amplitudes on
positive grades, allowing higher attachments to enter Hodge and Green queries.
The [relation schema](../../schemas/assay-relations.schema.json) defines this input.

The original relation file is retained. A packet's canonical relation document
records precise byte spans for each relation, participant occurrence and measured
field. Occurrence evidence quotes its own participant record. Inactive and
cancelling occurrences retain their separate source records.

## Execute native methods

A workspace query wraps an existing native scientific query:

```json
{
  "schema": "gadget.assays.query.v1",
  "operation": "native_query",
  "assay": "experiment",
  "query": {
    "operation": "field_summary",
    "field": "all.sample-A",
    "support": "reported_only"
  }
}
```

```sh
build/gadget/gadget assay-query --workspace data/workspace/workspace.json \
  --request request.json
build/gadget/gadget assay-query --workspace data/workspace/workspace.json \
  --request request.json --out results/experiment
```

The latter writes each partition's full native certificate and a result index.
Add `concurrent_jobs: 2` to the workspace request for two independent jobs; the
native query's `workers` declares its own supported worker budget. RexGraph owns
both reservations and completion. No application thread pool is added.

Matrix dimensions use the canonical checked `IndexDomain`. Source dictionaries,
stored values, intervals and decoded genomic records have no additional
application count ceiling. Complete-file intake retains every original byte and
checks the actual string/stream representation; JSON retains json-c's signed
integer input representation and grows its parser stack from complete source
bytes. Local assay and native-query requests use the same complete-file intake.
This adapter loads files and native packets into
memory. It does not provide disk-backed atlas queries or genome-scale alignment.
The compiled capability catalog separates this admission from the remaining
finite sequence-search and allocation profiles.

Storage partitions are ordered, disjoint subsets of one declared axis.
Pointwise paired contrasts and reported-support summaries can run per partition.
Hodge, Green, subspace and sheaf analyses require the whole carrier in one packet;
GADGET rejects applying those methods independently to storage partitions.
Relations likewise require their indivisible declared carrier. An assay may set
`partition_size` to a positive native coordinate count. For independent storage
partitions the existing default is 1,024 coordinates; this partitions the source
without limiting its total size. A relation-bearing carrier remains whole, and a
requested partition size that would split it rejects admission.
Packets refer to the retained original-file hashes; preserve the complete workspace
for source reproduction. `reimport.json` points to the retained files inside the
workspace and can replay admission after relocation; `import-spec.json` preserves
the original declaration. A selected packet can be used independently for native
calculation, but its containing workspace holds the original assay files.

The [research method examples](../research/README.md) specify matched contrasts,
controls, finite changes, stalks, restrictions, support, metrics and admissible
normalization evidence. The packet's `assay_declaration` document supplies the
byte spans for caller-declared identity scaling. Estimated factors can instead
be admitted with their own producing source and method.

Measurement-map coordinates come from each selected field's native SourceGrade
carrier in its actual order. Grade-zero fields use the admitted feature or
observation IDs; positive-grade fields use their retained relation IDs, such as
`relation:e01`. Transport reports both carrier grades. The map's domain is the
right-hand field and its codomain is the left-hand field; it compares the left
field with the transported right field. Equal vector lengths do not substitute
for these coordinate identities.

## Physical files and assay-specific questions

`inspect-file --spec file.json` validates a hash-bound FASTA, four-line FASTQ,
BED3–6, GTF or unweighted Salmon `--dumpEq` file. A file declaration names schema
`gadget.bioinformatics.file.v1`, `format` and a `source` path/hash record. Genomic
resources also declare assembly, sequence hash, canonical coordinate convention
and contig lengths. GTF one-based closed coordinates are explicitly converted
to zero-based half-open intervals. Contig names remain unchanged. The currently
admitted text profiles require uncompressed inputs.

| Workspace method | Scientific object retained |
|---|---|
| `interval_coverage` | Constant-coverage intersection segments, every interval membership, exact weighted coverage, integrals, and regions across caller-declared thresholds |
| `sequence_motifs` | IUPAC-defined sequence occurrences with positions, orientation, and definite versus possible matches at ambiguous bases |
| `allocation_bounds` | Sharp target-count bounds over nonnegative assignments within retained Salmon candidate classes |
| `field_transfer` | Full source and transported fields, map rank/kernel, residual and exact pairing under a declared measurement map |

`allocation-bounds --request file.json` also accepts general candidate classes,
multiple samples, arbitrary candidate subsets as targets, explicit scales and
matched sample contrasts. Classes are primary multiway relations. A target gets
a fixed contribution from a class wholly inside it; a mixed class contributes
between zero and its full count. These marginal endpoints are attainable under
the independent-class constraint. Bounds for different targets are not generally
jointly attainable. This profile needs no allocation enumeration, EM fit or
pseudocount. Upstream candidate search and filtering remain part of its scope.

Each class has its own retained occurrence coordinates. RexGraph computes the
exact kernel of each class-sum constraint and the image under the declared
target map. GADGET joins those images in common target coordinates before
requesting the native ambiguity rank. This rank describes the homogeneous
linear constraint system; a zero-count class retains its structural kernel even
though nonnegativity fixes its observed allocation to zero. The count bounds
report that feasible restriction separately. Thus rank and nullity are available for
all admitted sources, without a small-matrix fallback. Bound-map construction
visits actual candidate memberships rather than every target/class pair.
Samples, contrasts and unweighted Salmon classes follow the declared source.

Set `"kernel_certificate": true` in a direct allocation or workspace query to
retain the sparse class kernel frames and their target images. The certificate
records candidate order, occurrence offsets, relation identities and kernel
column offsets. Its source order reconstructs the full retained constraint
kernel; target-image columns share the ordered targets in the request.
`tests/sequence_assignment_cli.py` independently checks these factors, exact
bounds and shared-image ranks, including sources above the former count limits.

[The allocation-factor example](allocation-factors.json) retains three transcript
candidates, two kinds of ambiguity and two count vectors:

```sh
gadget allocation-bounds --request examples/bioinformatics/allocation-factors.json
```

Its second-minus-first gene-A contrast lies in `[-12,-9]` for every admissible
assignment, whereas gene-B lies in `[-3,0]`. Transcript-A1 lies in `[-10,-2]`.
The two independent allocation directions have common target-image rank two.
These intervals describe the declared fragment counts and assignment constraints;
population comparisons require a separately declared sampling question.

Threshold regions describe observed coverage under a specified criterion; they
are not a population significance test. Motif occurrences do not imply enrichment
or regulation. Cross-assay transport requires ordered source IDs, domain and
codomain quantities/units, coefficients, method and evidence. Links such as
`same_specimen` do not authorize equating RNA and protein values.

## Local models and workflows

The existing bounded model and native work-graph interfaces can consume a selected
assay packet directly:

```sh
build/gadget/gadget workflow --workspace data/workspace/workspace.json \
  --assay experiment --partition 0 --request workflow.json --models local-roles.json
```

The workflow's source digest binds that actual native partition. Models propose
admitted scientific queries or interpret certificates; their observations remain
separate from assay measurements. Branching, dependency graphs, recursive calls,
multiple local model roles and literal claim checks use the existing workflow
contract. The full assay and its dictionaries stay in the portable workspace.
Model context includes the complete mathematical result when it fits the bounded
view, including section vectors and residuals. Larger views identify their omitted
coordinates and retain the complete result digest. The full certificate remains
in the workflow trace. A model's prose observation is not a verified numerical
claim; literal claims require the existing source-bound claim-checking interface.

## Qualification and remaining services

The [chromatin fragment walkthrough](chromatin-fragments.md) supplies pinned
public inputs, native requests and complete independent source verification
across chromosomes. The [spliced RNA walkthrough](spliced-rna.md) covers original
alignment blocks, splice occurrences and exact feature assignment bounds.
The [spatial expression walkthrough](spatial-expression.md) covers complete
measured expression, physical coordinates, native local charts, compatible
sections, Hodge actions and shifted Green responses.

`tests/assay_cli.py` independently enumerates varied allocation fixtures, checks
coverage at every reference base, verifies native paired fields and spatial sheaf
gluing, compares dense/CSR/CSC H5AD and Matrix Market intake, and exercises missing
support, source tampering, unit mismatches and malformed records. The tests use
unrelated authored nonclinical data. The broader existing native tests remain
required.

The toolkit integrates admitted assays, genomic codecs, indexed gapped
alignment, exact exon-chain discovery, observed assembly complexes, minimum read
covers, constrained finite allocations and rational persistence. Streaming
multigigabase alignment, gapped exon alignment, error-correcting consensus,
unbounded optimization and integral/zigzag persistence remain separate services.
The compiled capability catalog reports each profile and its current scope.
Pair alignment retains the complete edit grid with reusable native route
histories and both prefix/suffix dual certificates. Literal reference matching
uses complete disjoint-seed candidates; IUPAC and declared-junction matching
retain exhaustive verification. Both methods use source-sized requests without
separate dataset or result ceilings. [The sequence example](sequence-witnesses/README.md)
includes all source inputs and a default pipeline. `tests/sequence_witness_cli.py`
checks full original-coordinate fields against independent exact oracles, all
accepted placements and source-only replay.

Splice discovery, read covers and observed assembly use source-sized declarations
and intrinsic native coordinate domains. Read-cover subset keys retain every
source identity beyond a machine word; native reachability and suffix certificates
exclude only impossible starts and nonminimizing prefixes. All minimum layouts,
spellings and inherited Gram entries remain retained. The
[composite example](composite-sequences/README.md) includes complete source files
and a pipeline selecting the three methods. `tests/composite_sequence_cli.py`
checks full independent layouts, categorical restrictions, exact rank/counts and
larger source domains. Complete source and output remain resident, with actual
branching and elimination costs.

Coupled bounded allocations use native affine equality kernels and certified
independent parameter groups, retaining complete feasible and optimal bounds.
The [method guide](suite.md) explains explicit and factored witness contracts;
the [coupled example](coupled-allocation-factors.json) provides a portable request.
Alignment, assembly and constrained allocation publish their generating complexes
and numerical fields directly through `assay-query --out-record` and
`derive-assays`. Ordered route lifts, source occurrence maps and affine factors
remain available in the binary publication. `tests/native_sequence_fields.cpp`
checks native support, variable bindings, downstream ownership and relocated
replay for these operations and coupled expression.
Statistical significance requires an explicit sampling/null model; exact finite
calculations do not supply that model. Existing DESeq2 comparisons provide a
benchmark rather than a proof that an alternative is superior for every assay.

[RNA evidence and expression bounds](rna-evidence/README.md) provides a complete source bundle and a materializable RNA-seq recipe. It retains transcript ambiguity, explicit gene memberships, cross-sample contrasts, declared coupling constraints and portable replay.

[The paired-assay example](assay-sections/README.md) demonstrates declared measurement interfaces, complete weighted compatible and residual fields, missing support, native CPU execution and portable multi-source replay. Its selected protocol also supplies bindings for the local-model recipe planner.

[Source contexts and scientific annotations](source-context/README.md) shows record and native-variable references carried through RNA answers, with optional local-model proposals and detached replay.

The [scientific claim example](scientific-claims/README.md) checks selected source facts and completed analytical results, retaining contradictions and omitted facts.

The [relational field ensemble example](relational-ensembles/README.md) works on
complete fields of the retained complex. Shared gradient and curl preparations,
a checked harmonic complement, Dirac and Green actions, weighted spread and
sector coupling serve the same pipeline for a declared scientific question. Its authored primary
relations and real paired-expression fields retain all original coordinates
and declarations.
