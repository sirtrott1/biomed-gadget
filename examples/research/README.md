# Endpoint and temporal structure analysis

These queries operate on an admitted `ScientificStudy`, through the same C++
session used by the CLI, scientific plans and local model proposals. Their
implementation is in `gadget/research.hpp` and `src/research.cpp`. Dataset names,
source joins and study questions belong in adapters and query declarations.

`controlled_pairings` preserves the complete patient fields after projection
off a declared control span. For a control frame F and the declared identity
metric, RexGraph computes P = F(F*F)^−1F* and R = I − P. The result retains
the original, represented and residual fields; ⟨Rx, Ry⟩; both residual
quadrances; and the full residual Gram matrix. The coefficient
⟨Rx, Ry⟩/⟨Ry, Ry⟩ describes the measured field along the residual endpoint
direction. It is undefined when that endpoint direction is zero. With a binary
endpoint and an intercept alone, it is the observed mean difference between
the two groups, in the measurement's own units. Adding controls changes that
descriptive quantity. It supplies no sampling significance or causal effect.

```json
{
  "operation": "controlled_pairings",
  "fields": ["all.morphology", "all.expression"],
  "endpoint": "all.recorded_response",
  "controls": ["all.baseline_volume", "all.age"],
  "selection": ["all.cohort_membership"],
  "support": "joint_reported",
  "metric": "identity",
  "intercept": true
}
```

`trajectory_structure` takes two to six quantitative channels at two to eight
ordered visits. It constructs a common vocabulary of the observed signs of
their exact changes. Each joint state is a relation containing one signed
change category per channel. Four observed states can support an explicitly
declared coordinate exchange identity. These analytical identities supply the
higher boundary; they are not inferred biological interactions. Magnitudes
remain in the complete patient change fields and their state aggregates.

The result contains native integral homology and torsion, boundary matrices,
patient witnesses, gradient/curl/harmonic fields and shifted Green responses.
It also distinguishes classes supported by at least two patients from classes
carried by consecutive intervals' shared subcomplex. The common quotient
dimension is rank(H_left) + rank(H_right) − rank([H_left; H_right]), where the
two maps project the same shared cycles into the intervals' harmonic spaces.
Matching Betti numbers alone does not identify persistent classes.

Optional `controls` remove the same declared patient space span from each
endpoint and change field before aggregation. If E maps patient coordinates
to observed states, the controlled endpoint harmonic field is P_h E R y.
Both its patient evaluation E* P_h E R y and its projection onto the image of
the repeated support subcomplex are retained. Each change field has its own
controlled harmonic pairing with this endpoint. Aggregation retains patient
multiplicity; these are unnormalized pairings, with the stated field units.

```json
{
  "operation": "trajectory_structure",
  "series": [
    {"id": "shape", "fields": ["all.shape_0", "all.shape_1", "all.shape_2"]},
    {"id": "extent", "fields": ["all.extent_0", "all.extent_1", "all.extent_2"]}
  ],
  "selection": [],
  "controls": ["all.baseline_volume"],
  "support": "joint_reported",
  "time": {"axis": "acquisition", "unit": "visit", "positions": ["0", "1", "2"]},
  "faces": "coordinate_exchange_identities",
  "semantics": "observed_feature_change_relations",
  "lambda": "1",
  "minimum_witnesses": 2,
  "endpoint_fields": ["all.recorded_response"]
}
```

`structure_transition` instead consumes two existing declared contexts in the
same native complex. It uses their retained ambient identities, full fields,
metrics and time declarations. It reports the boundary and metric pullback
defects, all nine Hodge sector transports, field innovation, and the exact
Green response change decomposed into realization, innovation and operator
terms. A nonzero boundary defect prevents interpreting the correspondence as
an induced homology map. This query can work with source authored tissue,
pathway or interaction complexes; it does not require a feature sign recipe.

```json
{
  "operation": "structure_transition",
  "context": "baseline",
  "other_context": "followup",
  "field": "baseline.G",
  "new_field": "followup.G",
  "lambda": "1"
}
```

The generic table importer now accepts lexical `identity_column` values and
declared `contrasts` such as `{"id":"change","left":"visit_1","right":"visit_0"}`.
Native rational subtraction uses their joint observed support and matching
units. Categories retain a lexical codebook and cannot enter linear recipes.
The original table, dictionary and every row remain source evidence.

Run a query with `gadget query --study study/study.json --query query.json`.
Run a declared batch with:

```sh
python3 tools/scientific_batch.py --binary build/gadget/gadget \
  --study study/study.json --protocol protocol.json --out results
```

`--resume` appends queries to a completed batch after verifying every previous
receipt, the unchanged input study, the original query declarations and the
compiled GADGET/RexGraph fingerprints. Changed earlier queries or sources are
refused. The earlier manifest and protocol remain as provenance. A single leaf
per plan is the default to bound full field receipt size; `--plan-steps` permits
up to sixteen.

Each batch leaf has an `id`, `description` and scientific `query`. Batches use
native plans of at most sixteen leaves. A plan admits every leaf before any
executes. The batch manifest records completion or failure, source binding,
result hashes and elapsed time. The runtime and tests are bounded: 4096
patient coordinates, 128 observed trajectory states, and 128 coordinates per
context in a structure transition. This exact rational profile runs on CPU.

The I SPY2 source adapter and experiment declarations are
`tools/prepare_endpoint_study.py` and `tools/declare_endpoint_experiment.py`.
They exercise the general interface; the C++ research module contains no
I SPY2 identifiers. Contract tests use a separate nonclinical dataset and an
independent Fraction oracle.
