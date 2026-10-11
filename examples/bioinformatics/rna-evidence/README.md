# From read evidence to expression bounds

This example runs the same expression method used for transcript sequences,
exon-compatible alignment records and unweighted equivalence classes. All input
files and their hashes are included. From the GADGET source root, with a built
`gadget` executable:

```sh
build/gadget/gadget import-assays \
  --spec examples/bioinformatics/rna-evidence/import.json \
  --out data/rna-evidence
build/gadget/gadget assay-query \
  --workspace data/rna-evidence/workspace.json \
  --request examples/bioinformatics/rna-evidence/query.json \
  > data/rna-evidence/expression.json
build/gadget/gadget pipeline-prepare \
  --workspace data/rna-evidence/workspace.json \
  --request examples/bioinformatics/rna-evidence/pipeline.json \
  > data/rna-evidence/plan.json
build/gadget/gadget pipeline \
  --workspace data/rna-evidence/workspace.json \
  --request data/rna-evidence/plan.json \
  --out results/rna-evidence
build/gadget/gadget replay-derivation \
  --workspace results/rna-evidence/expression/workspace/workspace.rfcxs \
  --assay expression --out results/rna-evidence-replay
```

The first sample contains five reads. One uniquely matches A1, two match both A1
and A2, one matches B, and one has no match. Each ambiguous read has three matching
windows: two in A1 and one in A2. It contributes one observation, retaining all
three placement witnesses. The second sample has three unique A1 reads, one
ambiguous A1/A2 read and one unmatched read.

| Observable | Baseline | Follow-up | Follow-up minus baseline |
|---|---:|---:|---:|
| A1 assigned counts | [1, 3] | [3, 4] | [0, 3] |
| A2 assigned counts | [0, 2] | [0, 1] | [-2, 1] |
| geneA, declared A1/A2 union | 3 | 4 | 1 |
| B assigned counts | 1 | 0 | -1 |
| Unassigned counts | 1 | 1 | 0 |
| All candidates | 5 | 5 | 0 |

GeneA changes by one count under every admissible allocation. Its individual
isoform changes are unresolved. This is an identification result for these
observations and the supplied transcript universe. It does not estimate a
population effect or supply a sampling significance test.

`intent.json` requests expression allocation ranges and inspection of the
baseline's stored quality scores. Both bindings have complete parameters. The
answer contract checks that a selection includes both outputs; it does not
choose a transcript allocation. Run the portable protocol with:

```sh
python3 examples/bioinformatics/replay-answers.py \
  --gadget build/gadget/gadget \
  --protocol examples/bioinformatics/rna-evidence \
  --out results/rna-answer-protocol
```

The run retains every read, candidate placement, class, target range and quality
histogram. `relocated-pipeline/answers.json` points to the full expression and
quality certificates. Reproduction moves the original source workspace and
regenerates the published expression assay from its retained input bytes and
generating request. No model is needed for this replay.

`expression_allocations` constructs shared classes from the complete candidate
sets. Every sample retains its local classes, original observations and mapping
to the common class coordinates. Targets are declared candidate sets; they may
be individual transcripts, gene unions or overlapping sets. RexGraph computes
sharp marginal bounds and the joint image rank of the class allocation kernels.
The kernel receipt concerns the signed linear ambiguity space. Nonnegativity and
observed counts determine the admissible allocation fibres; supported zero counts
remain distinct from missing observations.

Sample bounds in `expression` include the caller's positive rational scales.
The underlying `allocation_certificate` retains raw count bounds; its contrasts
already apply those scales. A scale must have a provenance declaration. There is
no estimated library-size factor, transcript-length correction or fitted
abundance. A positive lower contrast bound means the sign holds for every
admissible allocation; an interval containing zero retains the ambiguity.
Different targets' marginal endpoints need not occur in one joint allocation.

The optional `coupled-query.json` declares a between-sample A1 constraint and two
capacity bounds. Execute it with the same `assay-query` command. These authored
constraints fix A1 at two baseline and three follow-up counts. They are explicitly
additional assumptions, not findings inferred from the reads. Each coupling term
names a sample and a declared target and acts on **unscaled counts**. GADGET
translates those memberships and the original class conservation equations into
the existing bounded integer allocation service. Each output minimum and maximum
has a complete native solver request, certificate and original variable map.
The example explicitly requests factorized optimizer witnesses. A caller can
request the full explicit optimizer set instead. Infeasible constraints publish
an infeasibility certificate and cannot be materialized as measured expression.

Other evidence modes use the same `samples`, `targets`, `contrasts` and optional
`coupling` contract:

- `transcript_sequences` matches FASTA/FASTQ records against a hash-bound transcript
  FASTA under declared strands, IUPAC compatibility and a mismatch budget. Every
  contiguous placement is retained. Its coordinates are transcript coordinates;
  the current method does not align indels or infer new transcripts.
- `annotation_alignments` consumes SAM/BAM/CRAM with a declared annotation,
  reference frame, grouping, strand rule and filter. `junction_policy` is explicit.
  With `declared_exon_chain`, all aligned bases must lie in the feature's exon
  union, and every CIGAR N must exactly match a gap between adjacent declared
  exons. A deletion across that gap cannot stand in for a splice junction. A
  retained alternative mapping without a compatible feature remains an unassigned
  possibility. Paired ends are separate sequenced segments in this profile.
- `equivalence_classes` consumes unweighted `salmon_eq` fragment counts. Different
  source transcript orders and repeated equivalent class rows retain their source
  receipts and use an explicit exact aggregation map. Upstream candidate search
  and completeness remain declared evidence.

Annotation matching uses a complete contig/interval index, preserving original
feature order and all overlaps. Search sizes follow source geometry and checked
native coordinate representation. There is no feature-count, pair-comparison or
witness-count policy ceiling. Native kernels and coupled output can still require
substantial time and memory when the actual ambiguity is large.

Independent expression and feature assignment materialization retains the candidate
class complex and the exact class count, lower bound and upper bound fields.
Singleton targets use their original candidate coordinates. A grouped target has
an explicit target coordinate and a relation to its declared members. These
relations describe the assignment protocol. They do not assert tissue topology.
Class maps, bound maps and ambiguity images are retained as native tensors.
Supported zero counts remain distinct from coordinates outside a field's support.

The archive retains the original physical bytes and generating request. Pipeline
branches share the admitted source and fields. `replay-derivation` reconstructs
and compares the complete mathematical payload after the input workspace has
been removed. Coupled constraints consume typed class counts and exact solver optima. Every
solver retains its source complex, affine maps and complete optimizing factors
under a separate namespace. Explicit assignment maps connect those variables to
the original candidates and classes. Replay retains and checks each namespace.

`tests/expression_allocation_cli.py` independently enumerates finite assignments,
read placements, coupling constraints, extrema and optimizer residuals. It also
checks source-only replay, randomized candidate sets, large exact counts and
annotation dimensions exceeding the former application ceilings.
