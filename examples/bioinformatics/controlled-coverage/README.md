# Coverage contrasts with declared controls

This example constructs a shared physical carrier from weighted BED observations,
compares two explicitly declared changes, and separates the changes represented by
an input control and a constant field. The `controlled_coverage` method is available
in the ChIP-seq, ATAC-seq and custom pipeline profiles. These are authored inputs
for verifying the method, rather than measurements from a biological experiment.

From the GADGET source root, with a built executable:

```sh
python3 examples/bioinformatics/replay-answers.py \
  --gadget /absolute/path/to/gadget \
  --protocol examples/bioinformatics/controlled-coverage \
  --out build/controlled-coverage-example
```

The reference FASTA, all four BED files, their hashes, import declaration and
answer-contract intent are included. The driver executes the declared coverage,
Hodge, Green and chain-structure methods, moves the source and result directories,
then replays the producing coverage method from its retained files and request.
Use a new output directory for each run.

Every interval endpoint divides the physical union into elementary segments.
Each sample field is the sum of its retained interval weights on each segment,
divided by its declared positive sample scale. A scale is an explicit input with
provenance; the method does not estimate it. `measured_zero_within_union` is an
explicit measurement declaration. Coordinates outside the supplied union remain
unobserved. Strands can be retained separately or explicitly ignored.

| Field | [0, 4) | [4, 8) | [8, 12) |
|---|---:|---:|---:|
| Input control | 1 | 1 | 2 |
| Baseline | 3 | 1 | 2 |
| Treated | 5 | 0 | 2 |
| Replicate | 4 | 1 | 2 |
| Treated minus baseline | 2 | -1 | 0 |
| Represented by controls | 1/2 | 1/2 | 0 |
| Treated residual | 3/2 | -3/2 | 0 |
| Replicate minus baseline | 1 | 0 | 0 |
| Replicate residual | 1/2 | -1/2 | 0 |

The metric is the positive diagonal of physical segment lengths. Let `F` be an
independent frame of the declared control fields, with a constant column only
when `intercept` is true. For each complete contrast `x`, RexGraph solves
`(Fᵀ M F)a = Fᵀ M x` and returns `Fa` and `x - Fa`. GADGET retains the
original redundant control matrix, independent frame, Gram matrix, coordinates,
reconstruction check and exact orthogonality check. An empty control span returns
zero represented field and the complete input as residual.

For this example the residual quadrances are 18 and 2, and their signed metric
pairing is 6. These quantities integrate the fields over reference length. The
first residual exceeds the declared threshold 1 only on [0, 4). They describe
observed fields under a specified control span; they are not population p-values
or proof that the input control causally explains the contrast.

Each original interval is a distinct grade-one relation over its complete list
of physical segments. Repeated intervals and the zero-weight treated interval
remain present. Relation arity follows the source geometry. Optional `source_kernel`
returns the full native boundary kernel in original interval-variable coordinates:
its signed combinations have zero physical coverage. This identifies coverage
redundancy, without asserting tumor feedback or biochemical cycles. No higher
relation is created from physical overlap.

The pipeline materializes input, represented and residual fields as separate
native tensors. Their column metadata records the contrast and representation.
The original intervals, occurrences, source spans, generating request and complete
certificate accompany them. The same physical-length metric is retained at grade
zero, with an explicitly declared identity metric for the interval-variable grade.
The signed contrast coefficients also become grade-one native fields on the original intervals; their native boundary exactly reproduces the input coverage contrast. The selected Hodge step applies the native Hodge operator to these original interval coefficients, and the Green step acts on the residual coverage field, both through RexGraph in the retained grade metrics.
The selected chain certificate also retains the gradient, curl and harmonic sectors, complete occurrence contributions and integral homology. On this interval source the harmonic sector consists of zero-coverage combinations of original interval coefficients; curl is zero because no grade-two relation has been declared. The researcher can select local sections, persistence and other native methods
through the existing scientific query contracts.

Interval processing retains all atoms and memberships and uses sparse native
sample actions. Spatial neighborhoods use a sorted coordinate axis to exclude
pairs only when that axis contribution already exceeds the declared squared
radius. All remaining distances are evaluated with RexGraph's exact metric. This
preserves complete pair order and occurrence identities, including coincident
points and dense neighborhoods. There are no application ceilings on samples,
contrasts, thresholds, atoms, memberships, coordinate dimensions or neighborhood
relations; the actual source and result still determine memory and runtime.

`tests/physical_geometry_cli.py` independently evaluates per-base coverage,
rational control projections, source kernels, pairings, thresholds and all-pairs
spatial neighborhoods. It also checks invalid declarations and detached-source
replay. `tests/geometry_scale.cpp` exercises more than two million complete
memberships and verifies every resulting native coverage coordinate.

`intent.json` requests the control decomposition and its Hodge, Green and chain
analyses as separate outputs. Selecting a downstream method includes its declared
coverage prerequisite; the complete answer contract still requires all four
outputs. Each answer points to the full certificate in the generating pipeline.
The control projection and the chain sectors refer to their declared carriers:
physical coverage segments and original interval variables, respectively.

```sh
python3 examples/bioinformatics/replay-answers.py \
  --gadget build/gadget/gadget \
  --protocol examples/bioinformatics/controlled-coverage \
  --out results/chromatin-answer-protocol
```

The replay stages only the declared import, intent and four BED files. It moves
the original workspace before regenerating the published fields and validates
every answer artifact from the relocated pipeline.
