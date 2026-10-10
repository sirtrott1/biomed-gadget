# Spatial expression on declared local complexes

This experiment measures expression variation over the complete published
[10x Genomics breast cancer block A section 1](https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0)
count matrix and tissue coordinates. It asks whether each measured expression
field is compatible with declared restrictions on local spatial complexes, and
where its incompatible component occurs.

The input contains 36,601 genes and 3,798 tissue spots. Every released gene,
barcode and matrix record is retained. The selected 16 gene panel specifies the
questions asked; the other expression fields remain available in the native
study. Pixel coordinates refer to the original image frame. Counts are the
publisher's assigned UMI measurements. Source archives, original coordinate
rows and the image remain available through their complete original bytes.

## Reproduction

Build GADGET against the independent RexGraphCpp repository as described in the
[installation instructions](../../README.md). Enable examples and testing. Set
`GADGET_BUILD_DIR` to that build directory and `SPATIAL_DATA` to a new directory:

```sh
export GADGET_BUILD_DIR="$PWD/build/presets/cpu"
export SPATIAL_DATA="$PWD/data/spatial-expression"
python3 examples/bioinformatics/prepare_spatial.py \
  --data "$SPATIAL_DATA" --gadget "$GADGET_BUILD_DIR/gadget" --download \
  --source-namespace 983001 --radius-squared 75413 \
  --gene ESR1 --gene ERBB2 --gene EPCAM --gene KRT8 --gene KRT18 --gene KRT19 \
  --gene PTPRC --gene CD3D --gene CD3E --gene NKG7 --gene MS4A1 \
  --gene COL1A1 --gene COL1A2 --gene ACTA2 --gene VIM --gene MKI67
```

The preparation script verifies both complete archives by SHA256 before reading
them. Supply the archives under `originals/` and omit `--download` to run offline.
It writes four native request files. Any temporary JSON authoring files are
removed after conversion.

```sh
"$GADGET_BUILD_DIR/gadget" import-assays \
  --spec "$SPATIAL_DATA/import.rfcxs" --out "$SPATIAL_DATA/imported"
"$GADGET_BUILD_DIR/gadget" derive-assays \
  --workspace "$SPATIAL_DATA/imported/workspace.rfcxs" \
  --request "$SPATIAL_DATA/materialize.rfcxs" --out "$SPATIAL_DATA/derived"
"$GADGET_BUILD_DIR/gadget_verify_spatial_fields" --assay spatial-fields \
  "$SPATIAL_DATA/derived/workspace.rfcxs" "$SPATIAL_DATA/input/matrix.mtx" \
  "$SPATIAL_DATA/input/axis-observations.tsv" "$SPATIAL_DATA/input/axis-features.tsv" \
  "$SPATIAL_DATA/imported/workspace.rfcxs" "$SPATIAL_DATA/materialize.rfcxs"
"$GADGET_BUILD_DIR/gadget_scientific_sequence" \
  "$SPATIAL_DATA/derived/workspace.rfcxs" "$SPATIAL_DATA/sections.rfcxs" \
  "$SPATIAL_DATA/sections-result.rfcxs"
"$GADGET_BUILD_DIR/gadget_verify_spatial_sections" \
  "$SPATIAL_DATA/sections-result.rfcxs" "$SPATIAL_DATA/derived/workspace.rfcxs" \
  spatial-fields "$SPATIAL_DATA/section-summary.tsv"
```

To regenerate the derived study from its retained request and source:

```sh
"$GADGET_BUILD_DIR/gadget" replay-derivation \
  --workspace "$SPATIAL_DATA/derived/workspace.rfcxs" \
  --assay spatial-fields --out "$SPATIAL_DATA/replayed"
```

Run the field verifier on the replayed workspace with the same original inputs.
The verifier checks complete canonical encoding, every input measurement,
measurement support, original byte spans, ordered spatial occurrences and
declared chart membership. The section verifier independently checks every
reported section projection, coboundary defect, Hodge action and shifted Green
identity using exact rational arithmetic.

Generate the interactive figures from the independently verified table:

```sh
python3 examples/bioinformatics/render_spatial.py \
  --observations "$SPATIAL_DATA/input/axis-observations.tsv" \
  --features "$SPATIAL_DATA/input/axis-features.tsv" \
  --summary "$SPATIAL_DATA/section-summary.tsv" \
  --publication "$SPATIAL_DATA/sections-result.rfcxs" \
  --out "$SPATIAL_DATA/walkthrough.html"
```

The figures preserve the original coordinate frame and show measured counts,
anchor residuals, local residual quadrance and anchor Green responses for every
spot. Pointer inspection shows the original barcode and exact fractional values.

## Source and local topology

The completed run verifies 21,352,329 nonzero matrix entries and 82,852,292
assigned UMIs. All 139,010,598 gene and spot coordinates have declared
measurement support, including zeros. Checking every one of the 7,210,503 spot
pairs gives 11,032 physical neighborhood relations, 22,064 ordered occurrences
and 3,798 strict local charts. Direct publication, materialization and complete
replay pass the same independent source verification.

The spatial request declares the coordinate columns, frame, units, positive
diagonal metric, squared radius and relation type. A relation is added precisely
when the squared metric distance between two retained observations is at most
the declared squared radius. In this experiment, 75,413 squared pixels includes
immediate neighbors in the published spot layout. Each pair carries its two
ordered source occurrences with coefficients minus one and plus one. Proximity
does not create a higher cell.

Materialization preserves the original relational complex and expression
fields. New coordinate fields live on the observation carrier; squared distance
lives on the newly declared grade one relations. Existing grade one fields
retain their old support when the source carrier grows. New directions remain
absent from those fields. Original observation row spans are native tensor
fields pointing to the retained byte owner.

Each chart contains an anchor observation, all its incident grade one relations
and their complete boundaries. RexGraphCpp constructs its strict local closure.
The study declares `witness_scope: "block"` so local operators use exactly that
chart. Source cell identities remain primary keys. Compact display identifiers
such as `grade0:17` resolve to the actual native key and original barcode row.

The current materialization profile admits one complete observation carrier and
identity source metrics. Existing operators or other metrics require an explicit
transport when the carrier grows. A requested field must be fully measured on
the selected chart; unobserved entries are not filled with zero.

## Section, Hodge and Green calculations

For a gene field `x` on one chart, each cell has a scalar stalk and each retained
occurrence has a declared identity restriction. With the two oriented
occurrences, the section coboundary is the difference between the two measured
values. A compatible section is constant on each connected component.

RexGraphCpp computes the kernel of this coboundary and its metric orthogonal
projector. With the declared identity metric, a connected chart projects to its
exact mean `p`. The residual `r = x - p` retains every local observation. The
reported residual quadrance is the sum of squared residuals. It measures
departure from these declared compatibility conditions in squared UMI units.

The grade zero Hodge action is `L x = B B^T x`. It records the signed sum of
expression differences at each observation. The shifted Green field `y` solves
`(I + lambda L) y = x`; this experiment declares `lambda = 1`. The independent
verifier checks this identity coordinate by coordinate and checks preservation
of the constant section. Original counts, projected sections, residuals and
Green responses remain separate fields.

The summary table contains the actual anchor key, chart size, local mean,
residual quadrance, anchor residual, anchor Hodge action and anchor Green value
for each requested gene and chart. Comparisons use the same declared units,
metric and restrictions. A difference describes the observed source under
those declarations; its magnitude is not a population probability.
