# Complete sequence witnesses

This example includes the two FASTA inputs, their pinned hashes and a pipeline
request. Run it from the GADGET source directory with the compiled executable.

```sh
build/gadget/gadget import-assays \
  --spec examples/bioinformatics/sequence-witnesses/import.json \
  --out data/sequence-witnesses
build/gadget/gadget pipeline-prepare \
  --workspace data/sequence-witnesses/workspace.json \
  --request examples/bioinformatics/sequence-witnesses/pipeline.json \
  > data/sequence-witnesses/plan.json
build/gadget/gadget pipeline \
  --workspace data/sequence-witnesses/workspace.json \
  --request data/sequence-witnesses/plan.json \
  --out results/sequence-witnesses
```

The pair analysis retains the complete edit grid, an ordered minimum route,
every step on any minimum route, and prefix/suffix costs at every source
coordinate. Costs refer to retained native route histories. RexGraph checks both
dual fields and realizes the chosen route as its original relational field.

Reference matching retains all placements with at most one mismatch on either
strand, including both source spans, and lists the unmatched read. Disjoint
literal seeds identify a complete candidate set; the native mismatch operator
checks each candidate over the entire read. The result distinguishes the full
placement universe from the candidates requiring that operator action.

With `ambiguity: "literal"`, symbols match only equal tokens. On a reference
containing only A, C, G and T, an `N` in a query consumes one mismatch. Disjoint
seeds containing such a token cannot be an intact seed of an admitted placement;
the index omits those seeds and verifies each candidate over the entire read.
References containing ambiguous tokens retain the complete fallback for ambiguous
queries, so literal `N` against `N` can still match. The `possible` and `definite`
IUPAC rules are separate declared questions.

Save this request as `reference-matches.json` to select contiguous placements:

```json
{
  "schema": "gadget.assays.query.v1",
  "operation": "reference_sequence_matches",
  "query_resource": "queries",
  "reference_resource": "reference",
  "strands": ["+", "-"],
  "ambiguity": "literal",
  "maximum_mismatches": 1,
  "minimum_overhang": 1,
  "junctions": [],
  "evidence": "Every original query in the imported sequence example"
}
```

Retain its complete native result with:

```sh
mkdir -p results
build/gadget/gadget assay-query --workspace data/sequence-witnesses/workspace.rfcxs \
  --request reference-matches.json --out-record results/reference-matches.rfcxs
```

All accepted mismatch strata and repeated placements remain in the result;
the method also retains every unmatched query occurrence and its original bytes.

The native work graph executes these independent pipeline steps. The pipeline
index names both result files and their digests. The imported workspace retains
the original source files. The pair result can also be materialized with
`derive-assays` and reproduced with `replay-derivation` from its retained lineage.
See [the method guide](../suite.md) for the request and arithmetic contracts.

The separate `indexed_sequence_alignment` query accepts an optional `execution`
object with positive `read_partition_records` and `reference_block_bases` values.
Defaults are 256 records and 65536 bases. Every read and reference coordinate is
processed; these values set buffer widths. Seed scans continue across source
lines and block boundaries, reset at record boundaries and ambiguous bases, and
retain every repeated hit. Candidate windows address the original contig
coordinates. The result reports actual scans, batches and normalized buffer
sizes. Original byte sources, seed hits and complete results remain retained,
so these buffer sizes are not a bound on process memory.

Indexed edit search accepts `edit_histories: "one_minimum"` (the default) or
`edit_histories: "all_minimum"`. The second mode retains every distinct ordered
minimum edit history for each placement, including ties with equal costs or
aligned spellings. `minimum_edit_history_count` gives the exact family size,
and `minimum_edit_histories` contains the aligned query/reference words and
native route witnesses. Query and reference prefixes strictly advance along
these edit declarations, so zero cost matches cannot introduce cycles.

The native result also retains the prefix/suffix certificates, absent suffix
support and complete optimal transition support. Keep the complete parent
publication when selecting a child route: the parent owns the original query
and reference bytes and their coordinate maps. Complete families can require
large output storage. This indexed unit edit option preserves the existing
placement objective and does not define an error consensus.

`sequence_splice_discovery` accepts paired `maximum_edits` and
`scoring: "unit_edit_cost"` controls. This mode permits substitutions, insertions
and deletions within each exon. Each exon must consume at least
`minimum_exon_bases` query bases and reference bases. The existing intron bounds,
junction limit and transcript oriented motifs apply; an intron adds zero edits.
For each original query, strand, reference start, reference end and exact
junction count, the result retains every minimum ordered history within the
edit budget. Equal spelling or cost does not merge histories or source records.
Without these controls, the existing exact exon chain search applies.

Each returned history has complete aligned exon words, junctions and native
route witnesses. `minimum_history_groups` reports the exact family counts.
Tied histories reference one retained certificate group in the parent.
The parent publication retains original query, exon and intron byte spans,
including original framing and quality bytes. Native children retain the actual
biological state coordinates, operation fields, minimum certificates and absent
suffix support. Terminal observations have no biological prefix coordinate.
Exon span state coordinates are capped at the declared minimum; returned blocks
retain their complete query and reference intervals.
This complete search considers every reference start and may need substantial
state and output storage. It uses unit edits and declared sequence motifs;
annotation, affine penalties, quality weighting and biological consensus require
their own scientific declarations.

`read_edit_consensus` declares a complete read consensus profile. Required controls
are `resource`, per read allowed `strands`, total `maximum_edits`, positive
`minimum_paired_bases`, `objective: "minimum_length_then_total_unit_edits"`,
`ambiguity: "literal_ACGT"` and `evidence`. Every original read token is consumed
by a paired observation or a costed insertion, including the read's tails.
Each emitted position has at least one exact paired original token. The result
first minimizes output length, then total substitutions, insertions and
deletions, and retains every joint alignment and orientation tie.

The parent owns original reads, tokens, quality bytes and primary multiway
consensus evidence. A separate native complex retains the optimization routes
and certificates. Equal spellings retain distinct alignments; contained,
repeated and contradictory reads keep their source identities. Empty reads
and unresolved tokens are refused by this literal profile; impossible paired
requirements return `no_admissible_consensus`. This finite exact search can
require substantial state and output storage. Its declared result is a
consensus hypothesis under these constraints; quality, affine penalties and
other assembly objectives need separate declarations.
