# Matched RNA and antibody assignment fields

This workflow uses the complete captured 10x RNA/antibody matrix and its original
feature and barcode dictionaries. The retained capture has 713 barcodes, 33,538
RNA features and 17 antibody features. The modality matrices retain every count
and the original column order. Their independent native sources and unit
declarations remain separate.

Matching comes from each modality's projection of the same original matrix
column and barcode artifact. The preparation checks every declared artifact
hash, the complete partition by feature type, unique barcode identities and exact
column order. Matching labels from two unrelated captures would require a
separate supported correspondence. Repeated or nonbijective columns cannot use
this capture profile.

The selected question is the signed balance between a cell's total RNA
assignments and all its antibody assignments. Native aggregation sums every RNA
feature and preserves all 713 observations. This is an explicitly selected
library statistic, with its complete map and original source retained. Each
antibody field keeps all 17 original feature coordinates, including control
features and measured zeroes.

The comparison has an RNA stalk with one coordinate and an antibody
stalk with 17 coordinates. The declared coboundary is

```text
D = [1, -1, ..., -1]
```

Each stalk has the declared positive identity cost in its own measurement
domain. The common cost is a comparison choice, without biochemical calibration,
sampling variance estimation or cell type inference. The original RNA and
antibody UMI units remain attached to their fields.

For input `x`, let `d = x_RNA - sum(x_antibody)`. The residual is `d/18` on RNA
and `-d/18` on every antibody coordinate; its quadrance is `d^2/18`. Subtracting
that residual yields the compatible field. The complete kernel has 17 joint
directions: each direction changes RNA and one antibody together. Independent
coefficient intervals for each view would lose that correlation. The native result
retains the complete kernel, original coordinate maps, residual frame, metric,
restriction maps and all 713 observed decompositions.

From the repository, with a locally captured source directory containing
`paired-import.json` and its declared artifacts:

```sh
python3 tools/prepare_paired_assay_sections.py \
  --binary build/presets/cpu/gadget --source CAPTURE_DIRECTORY \
  --out NEW_OUTPUT_DIRECTORY --source-namespace 976300 --workers 2
```

The tool writes the complete explicit `correspondence.json`, pipeline design,
admitted plan and source receipt. The existing native pipeline imports the
original sources, derives the RNA library statistic, binds the retained owners
and calculates the declared sections. The complete operation is
`analysis/sections/result.rfcxs`; `result.json` is its optional display. Native
archives retain their source and request owners. The source receipt identifies
every supplied artifact hash and the selected scientific scope.

This is a library assignment comparison on one captured donor, preserving the
publisher's barcode selection and UMI assignment. It does not provide treatment,
clinical endpoint or modality calibration evidence. Capture attribution and
dataset distribution terms stay with the original source records and remain
separate from GADGET's Apache 2.0 license.
