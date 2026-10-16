# Compatibility of paired assay fields

This example separates the components of three quantitative assay fields that
satisfy declared measurement interfaces from the components that violate them.
It retains each assay's coordinates, units and evidence channel. The small RNA,
ATAC and protein matrices are constructed examples. Their restriction maps and
comparison costs illustrate the method; they are not experimental calibrations.

The five sections represent two complete observations, one incomplete
observation, a measured zero field and an entirely unobserved field. The interface
family includes two comparisons, a three assay balance, two cancelling occurrences,
an inactive occurrence and an occurrence with zero magnitude. All are retained.
The same operation accepts any number of declared views, observables and sections.

## Mathematics

Each view is a grade zero identity in a declared interface complex. Its stalk
contains all coordinates of the selected quantitative field from that view's
actual native source grade. These interface grades describe the comparison
construction; they do not relabel the biological source complex. Each interface
is a primary grade one relation of arbitrary arity. No pairwise expansion or
higher biological relation is inferred.

For occurrence `j`, the caller supplies a restriction matrix `rho_j` from a view
stalk to the interface's observable stalk. Its contribution to the cellular sheaf
coboundary is `e_j (1 - 2 o_j) q_j rho_j`, using its original existence,
orientation and magnitude. Repeated occurrences remain distinct even when their
contributions cancel. The three assay example declares a signed balance, rather
than pairwise equality of its inputs.

RexGraph's `LinearStalkBundle` and `LinearCellularSheaf` assemble this coboundary
`D`. The compatible section space is `ker(D)`. In the caller's positive diagonal
metric `M`, its orthogonal complement is `im(M^-1 D^T)`. An independent native
image frame `F` supplies the residual action

```text
r = F (F^T M F)^-1 F^T M x
z = x - r
```

The result checks `D z = 0`, `z^T M r = 0` and `x = z + r` exactly. It reports the
full `z` and `r` in every original view, their residual quadrance, the sparse
coboundary, independent frame, Gram system and coordinate maps. Local stalks, incidence restrictions and the declared metric determine the
compatible section and residual. Numerical assembly, rank, metric actions and
prepared solves are supplied by RexGraph C++.

With `complete_families: true`, each prepared support group also retains the complete `kernel_frame` and
`row_dependencies`. The kernel maps free coefficients to all observed stalk
coordinates and satisfies `D K = 0`; its columns number the compatible
dimension. They are joint directions, so coefficients cannot be chosen
independently for each view. Row dependencies satisfy `D^T L = 0` and retain
every relation among the assessed constraints. Their carriers name the original
selected coordinate and row maps. Omitting this output selection preserves the
earlier projection only contract and exact replay of existing native archives.

A comparison cost across different measurement units must be declared. Each view
metric names its measurement domain and the common quadrance codomain. The
example's diagonal coefficients are `1,2`, `3,4` and `5,6`. These are supplied
costs, not variance estimates or an implicit normalization of the assays.

## Missing measurements

An interface coordinate is assessed only when every nonzero, active occurrence
map dependency is observed. This rule applies before cancellation of repeated
occurrences. Missing measurements remain JSON `null` and remain absent in the
published native tensor. Supported zeros participate in the calculation.
Inactive and zero magnitude occurrences retain their maps but introduce no
numerical dependency.

For an incomplete section, the projection uses only observed coordinates and
assessed interface rows. Unassessed rows name the missing coordinates and have no
reported defect. Its compatible dimension and quadrance therefore refer to that
reported domain. The result does not assert compatibility at unseen coordinates.

## Results of the declared example

| Section | Observed coordinates | Assessed interface rows | Residual quadrance | Complete |
|---|---:|---:|---:|---|
| before | 6 | 6 | 1420/51 | yes |
| after | 6 | 6 | 0 | yes |
| partial | 5 | 4 | 0 | no |
| zero | 6 | 6 | 0 | yes |
| unobserved | 0 | 2 | 0 | no |

The complete `before` field violates the declared interfaces. For example, its
RNA and ATAC first coordinates are 6 and 2; their compatible values are both 3
under their different declared costs. The residuals are 3 and -1. The complete
`after` field satisfies all declared balances and remains unchanged.

The partial field has zero residual on the assessed domain, but its missing RNA
coordinate leaves two rows unassessed. It is not reported as a complete compatible
section. An unobserved field similarly has zero residual because no measured
coordinate participates. Only the inactive and zero magnitude rows are assessed.
This distinguishes an observed agreement from absence of comparable measurements.

## Run and reproduce

Build GADGET against the canonical RexGraph C++ checkout as described in the
[root README](../../../README.md). The mathematical example needs no model or
optional file format library. From the GADGET repository:

```sh
python3 examples/bioinformatics/assay-sections/replay.py \
  --binary build/gadget/gadget --out build/paired-assays
```

All required source matrices, dictionaries, hashes, queries and pipeline
declarations are in this directory. Outputs include the complete section
certificate, prepared pipeline, compatible RNA workspace and replay certificate.
The pipeline publishes the selected RNA compatible field and then inspects its
support. Publication retains every generating assay's original files and the
complete generating query and certificate. The replay rebuilds the calculation
from that archive after the generating workspace has moved.

`materialize.json` can instead publish this view directly with `derive-assays`.
Its `view` and `component` select one original view and `input`, `compatible` or
`residual`. The published quantity records the derivation, preserves the source
units and channel, and inherits the declared view metric. It preserves the
original assay complex, source grade, cell names and occurrence identities. The
complete interface complex, stalk restrictions and projection results remain in
a separate native publication inside the archive.

## Selecting the protocol

`assay_sections` is available in the `single_cell`, `multi_omics` and `custom`
pipeline profiles. Scientists supply the measurement maps and metrics and select
which sections to examine. Existing count preserving aggregation can first form
specified specimen or cell groups. Unrelated measurements need no comparison
unless the scientific question and declared maps call for one.

`intent.json` binds the example's method and support inspection for the existing
local model planner. A model may select those bindings; GADGET checks the source
and intent hashes and includes their required dependencies. With a build that
enables local llama.cpp inference and a supplied local model:

```sh
build/gadget/gadget plan-recipes \
  --workspace build/paired-assays/generating-workspace/workspace.json \
  --intent examples/bioinformatics/assay-sections/intent.json \
  --model /path/to/local-model.gguf \
  --prompt 'Inspect the compatibility of the declared assays and the support of the resulting RNA field' \
  --out build/model-selected-assays
```

Model selection does not define a new measurement map or populate missing values.
The declared numerical calculation can run without a model. Repeated sections
with identical support share a prepared exact solve. Native CPU tiles execute
independent section actions with joined completion and one caller inclusive
reservation. The factorized system has the constraint rank as its dimension;
complete source fields and result coordinates are retained.

## Read the retained sections

An assay query can write its complete owning native result with `--out-record`.
Use the generic section reading to inspect its original observations, exact
compatible and residual fields, missing coordinates, declared costs and complete
coefficient families:

```sh
build/presets/cpu/gadget assay-query --workspace workspace/workspace.rfcxs \
  --request examples/bioinformatics/assay-sections/query.json \
  --out-record sections.rfcxs
python3 tools/render_assay_sections.py --gadget build/presets/cpu/gadget \
  --operation sections.rfcxs --out section-reading
```

Open `section-reading/report.html`. Every declared observation appears once in
source order; its index is not biological time or an inferred spatial coordinate.
Charts convert values only for display. Tables preserve exact rational values
and absence, and the report links to the unchanged native owner and its byte hash.
Keep that native file with the reading when relocating results. No plotting or
model package is required by this reading tool.
