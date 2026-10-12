# Exon chains, read covers and observed sequence structure

This source bundle supplies all three FASTA files, their pinned hashes and a
pipeline request selecting three related analyses. From the GADGET source root:

```sh
build/gadget/gadget import-assays \
  --spec examples/bioinformatics/composite-sequences/import.json \
  --out data/composite-sequences
build/gadget/gadget pipeline-prepare \
  --workspace data/composite-sequences/workspace.json \
  --request examples/bioinformatics/composite-sequences/pipeline.json \
  > data/composite-sequences/plan.json
build/gadget/gadget pipeline \
  --workspace data/composite-sequences/workspace.json \
  --request data/composite-sequences/plan.json \
  --out results/composite-sequences
```

The exon-chain result includes the complete spliced and contiguous query
placements, source spans, transcript-oriented junction motifs and ordered
native path fields. Exact query reconstruction checks each chain. The declared
exon, intron and junction parameters define the research question.

The read-cover result reconstructs `AACGTTGC` from the three distinct source
reads. Native reachability proves which reads can start a complete cover.
The search carrier retains every reachable prefix state from eligible starts,
with exact source subsets and last-read identities. Native routes certify its
minimum. Every minimum layout retains categorical gluing and restrictions that
recover all original reads. Alternative spellings, when present, retain their
complete shared-field Gram on the declared transition carrier.

The observed assembly result retains every literal 3-mer, every source
occurrence and each read's ordered walk. Native count actions and boundary rank
produce the occurrence field and cycle dimension. The resulting carrier can be
materialized with `derive-assays` for native field, Hodge, Green, local-section
and persistence queries. `replay-derivation` reproduces it from retained lineage.

The pipeline index names each complete result and its digest. It uses the native
work graph to join the independent analyses in declaration order. See
[the method guide](../suite.md) for their complete source and arithmetic contracts.

`intent.json` requests three distinct outputs: complete exact exon chains,
minimum read-cover reconstructions and observed k-mer structure. The source reads
have the unique minimum spelling `AACGTTGC` under the declared suffix-prefix
objective. That spelling is reconstructed by restricting it to every retained
read and checking complete overlap gluing. The k-mer carrier preserves branching
and repeated source occurrences separately from this reconstruction objective.

```sh
python3 examples/bioinformatics/replay-answers.py \
  --gadget build/gadget/gadget \
  --protocol examples/bioinformatics/composite-sequences \
  --out results/sequence-answer-protocol
```

`relocated-pipeline/answers.json` locates the full exon, reconstruction and k-mer
certificates. The replay includes only the declared import, intent and three
FASTA files. No previous workspace, development archive or model is required.

`read_protocol_evidence` extracts declared segments from every original FASTA or
four-line FASTQ record. Its query supplies `resource`, `resource_sha256` and a
`protocol` with `id`, `evidence`, `count_unit: "read_record_occurrence"` and
`segments`. Each segment declares `id`, `role` (`barcode`, `umi` or `sequence`),
zero-based half-open `start` and `end`, and `strand` (`+` or `-`). Reverse segments
use the DNA IUPAC reverse complement. Original source bytes retain case, headers,
line endings and qualities; extracted tokens use uppercase. A short record keeps
its occurrence and reports each unavailable segment separately.

The optional `protocol.barcode` declares the barcode `segment`, an original
FASTA `whitelist_resource` and `whitelist_sha256`,
`objective: "minimum_substitutions"`, a nonnegative `maximum_substitutions`, and
`symbol_rule: "literal_iupac_tokens"`. Every whitelist record must have the
declared segment length. Tokens are compared literally, including ambiguous
IUPAC symbols. All equally minimizing whitelist occurrences within the budget
remain candidates in original dictionary order. Distinct identifiers with the
same sequence remain distinct candidates. The reported minimum remains present
when it exceeds the budget; an unavailable barcode has no minimum.

UMI segments retain their observed tokens. This operation counts read record
occurrences and does not combine them into molecules. The supplied coordinates
and whitelist are protocol evidence; they require the caller's source hashes.
Every result retains its workspace and original native request owners. Write
the complete native result with:

```sh
build/gadget/gadget assay-query \
  --workspace data/composite-sequences/workspace.rfcxs \
  --request protocol-query.json \
  --out-record results/read-protocol.rfcxs
```

Omitting `--out-record` renders the optional complete view. Record indices and
original byte spans identify occurrences, while segment coordinates address
positions within the original read.
