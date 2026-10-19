# Controlled ChIP read coverage

This example compares two EGR1 ChIP libraries with two input control libraries
from the HudsonAlpha ENCODE K562 experiment. It retains every chromosome 21
read interval in the four complete published alignment archives and uses the
general `controlled_coverage` operation. The library identifiers and data
preparation belong to this example; the runtime method accepts caller declared
resources, samples, controls and contrasts.

## Data and design

The [UCSC archive](https://hgdownload.soe.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeHudsonalphaChipSeq/)
labels both EGR1 libraries as K562, PCR2x, treatment `None`, with control family
`K562/None/Input/PCR2x`. The two selected input libraries have the corresponding
cell, protocol and treatment. Replicate labels are retained as published;
the analysis does not infer additional pairing between the four libraries.
The adapter retains the original directory listing, file declarations,
publisher checksums and complete compressed alignment files.

The measured intervals are the original 25 base mapped reads. The separately
reported fragment length is not used to extend or shift them. Each original
row contributes one read, including duplicate rows. BED scores are retained
as original bytes and are not used as counts. The reference is the complete
[hg18 chromosome 21 sequence](https://hgdownload.soe.ucsc.edu/goldenPath/hg18/chromosomes/chr21.fa.gz).

| Library | Complete archive rows | Chromosome 21 rows | Read base integral |
|---|---:|---:|---:|
| EGR1 replicate 1 | 13,243,731 | 214,016 | 5,350,400 |
| EGR1 replicate 2 | 10,466,317 | 171,406 | 4,285,150 |
| Input replicate 1 | 16,377,339 | 266,884 | 6,672,100 |
| Input replicate 2 | 14,363,052 | 252,523 | 6,313,075 |

Chromosome 21 is a declared physical scope, not a random sample or a substitute
for genome wide qualification. The adapter scans every row of each complete
archive, checks the pinned SHA256 and publisher MD5, and writes selected rows
without changing their bytes or order. No read depth normalization, duplicate
removal or fitted significance model is applied. All sample divisors are one.

## Relational source and fields

GADGET builds the grade zero cells from the complete endpoint partition of the
observed interval union. Each original interval is a distinct grade one cell;
its ordered boundary occurrences name the physical cells it covers. Overlapping
reads and duplicate reads remain separate relations. An interval with many
members retains that arity directly. No pairwise clique or inferred higher
attachment is introduced.

The complete source, field carriers, measured support and original byte spans
are published through RexGraphCpp. A sample field stores the read coverage of
each physical cell. A zero is measured within the declared union; positions
outside that union are outside this carrier. Original strand values are
retained, while this query declares strand independent coverage.

For physical cell lengths `M`, the two input fields form the declared control
matrix `C`. RexGraphCpp obtains an independent image frame `F` from those actual
columns and computes the exact metric projection

```text
P = F (Fᵀ M F)⁻¹ Fᵀ M
represented = P y
residual = y − P y
Cᵀ M residual = 0
⟨y, y⟩M = ⟨represented, represented⟩M + ⟨residual, residual⟩M
```

All coefficients are rational. The runtime uses canonical image frame and
metric solve operations; it does not fit an embedding or compute an eigensystem.
The three declared fields `y` are EGR1 replicate 1, EGR1 replicate 2, and their
difference. Each result retains the complete input, represented and residual
fields, frame coordinates, quadrances and residual Gram matrix.

The represented component measures the part of a library's spatial coverage
in the span of the two observed input controls. The residual describes coverage
outside that span on this physical domain. Calling a residual a binding site
would require additional experimental evidence and a separate declared method.
The replicate difference exposes disagreement between the measured libraries
in their original units.

## Recorded result

The 904,829 original reads produce 1,279,572 physical cells and 5,454,865 ordered
boundary occurrences. The two original input fields have independent rank two.
The verifier checks every source row and every occurrence, including reads with
multiway membership in the physical partition.

| Contrast | Input quadrance | Represented quadrance | Residual quadrance | Represented share |
|---|---:|---:|---:|---:|
| EGR1 replicate 1 | 39,366,174 | 4,906,815.264 | 34,459,358.736 | 12.465% |
| EGR1 replicate 2 | 9,363,280 | 1,647,481.834 | 7,715,798.166 | 17.595% |
| Replicate 1 minus replicate 2 | 26,219,936 | 896,617.008 | 25,323,318.992 | 3.420% |

Quadrance uses squared read count units times reference bases. Displayed decimals
are rounded; the publication and verifier retain rational coefficients and
check the identities exactly. The two ChIP libraries have different spatial
coverage fields, and most of their difference remains outside the declared input
span. This is a comparison of these measured libraries at their original read
depths. It does not distinguish the experimental causes of their disagreement.

## Reproduce

Build GADGET against the independent RexGraphCpp project as described in the
[build instructions](../../cmake/README.md#build). Enable examples to build
the application and verifier. The following commands use repository relative
paths and download the complete public inputs into a new data directory.

```sh
cmake --preset cpu -DGADGET_BUILD_EXAMPLES=ON
cmake --build --preset cpu
python3 examples/bioinformatics/prepare_chip.py data/chip \
  --gadget build/presets/cpu/gadget --source-namespace 908105 --download
build/presets/cpu/gadget import-assays \
  --spec data/chip/import.rfcxs --out data/chip/imported
build/presets/cpu/gadget assay-query \
  --workspace data/chip/imported/workspace.rfcxs \
  --request data/chip/coverage.rfcxs --out-record data/chip/coverage-result.rfcxs
build/presets/cpu/gadget_verify_controlled_coverage \
  data/chip/coverage-result.rfcxs data/chip/imported/workspace.rfcxs \
  data/chip/coverage.rfcxs data/chip/input/egr1_1.bed data/chip/input/egr1_2.bed \
  data/chip/input/input_1.bed data/chip/input/input_2.bed
```

The source namespace is an explicit nonzero native identity. Select another
unused namespace when combining this publication with existing studies.
The preparation script uses a temporary compatibility declaration to author
the requests and removes it after native conversion. Requests and scientific
results use RFCXS; the original files remain the byte provenance.

Materialization makes the resulting source and fields available to general
scientific queries. Replay regenerates the publication from the retained inputs
and compares the complete source, carriers, fields, supports and evidence.

```sh
build/presets/cpu/gadget derive-assays \
  --workspace data/chip/imported/workspace.rfcxs \
  --request data/chip/materialize.rfcxs --out data/chip/derived
build/presets/cpu/gadget_verify_controlled_coverage --assay chip-coverage \
  data/chip/derived/workspace.rfcxs data/chip/imported/workspace.rfcxs \
  data/chip/materialize.rfcxs data/chip/input/egr1_1.bed data/chip/input/egr1_2.bed \
  data/chip/input/input_1.bed data/chip/input/input_2.bed
build/presets/cpu/gadget replay-derivation \
  --workspace data/chip/derived/workspace.rfcxs --assay chip-coverage \
  --out data/chip/replayed
build/presets/cpu/gadget_verify_controlled_coverage --assay chip-coverage \
  data/chip/replayed/workspace.rfcxs data/chip/imported/workspace.rfcxs \
  data/chip/materialize.rfcxs data/chip/input/egr1_1.bed data/chip/input/egr1_2.bed \
  data/chip/input/input_1.bed data/chip/input/input_2.bed
```

The independent verifier reads all original BED rows in declared sample order.
It reconstructs every physical cell by an event sweep, checks every interval
identity and byte span, and checks every ordered boundary occurrence. It then
verifies all sample fields, the complete declared control span, projection
coordinates, orthogonality, mass conservation, metric decomposition and
cross contrast residual Gram. It also checks the original request and workspace
bytes and exact canonical archive reencoding. The same executable accepts other
BED workspaces, rational divisors, dependent controls, an optional intercept,
separate strands and any declared number of samples or contrasts.

The native workspace shares one immutable payload for equal retained artifacts,
after verifying their complete byte equality. Each original path keeps its
separate declaration. This reduces repeated reference and provenance payloads
without merging scientific carriers. Portable operation publications still
retain their complete generating workspace; current admission is in memory.
The resulting native study remains queryable from `workspace.rfcxs` alone;
the recorded qualification also reopens it after removing generated side files.
