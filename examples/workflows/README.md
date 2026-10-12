# Source bound workflows

A workflow declares scientific queries, biomedical evidence queries, local model
roles, branches, dependencies and bounded recursive calls. Numerical leaves use
`ScientificSession` or `EvidenceStore`; graph jobs use the existing WorkGraph
scheduler and caller resource coordinator. The declaration's node, depth, model
call and parallel worker budgets control execution. Readiness and execution
clocks establish no biological observation time.

`execute_workflow_record` retains each completed leaf's actual native operation,
including its complete fields, carriers, support, free families and original
request. A leaf query's request owner is the original workflow publication.
Parallel results remain in declaration order. Independent publications remain
separate even when their dimensions or values agree. Native result references
use `gadget.workflow.publication.v1`, their namespace and cutoff, and a complete
RFCXS archive in `publication`. Repeated references to one actual owner reuse
its archive. Model JSON replies remain immutable `view_json` presentation text.
The workflow receipt preserves failed replies and completed source owners too.

The current strict scientific profile admits native model comparison, finite
reference calibration, Hodge decomposition, linear response observation,
declared response transition, metric relative novelty, observed/constrained
response models and relational field ensembles. Biomedical leaves use native
role/composite evidence operations. An unconverted scientific publisher is
refused during admission, including in an untaken branch or unused definition.
`execute_workflow` provides the existing compatibility display through the same
Runner; legacy publisher output remains an explicit migration boundary.
`workflow_record_view` displays a saved native workflow without executing it.
Generic record views expose the structural receipt and original archive refs.

A native branch or claim node reads actual typed values. A `condition` or claim
item declares `path`, `value_type`, `comparison` and expected `value`. Paths start
within the producing operation's scientific view or biomedical result. For a
nonscalar tensor add `coordinate`, its original flat coordinate, independently
of a compact display's selected-coordinate order. An unsupported coordinate
refuses a numerical comparison; reported zero remains comparable. Claims can
be contradicted and a later branch can select an explicit repair or base case.
The branch's retained `comparison` owns the literal check and its source pins.
Numerical agreement does not supply clinical or biological meaning.

A local role can return prospective JSON controls or use
`WorkflowModel::NativeInvoke` to return its actual native tool publication.
Strict workflows refuse JSON replies labelled as executed mathematical results.
The explicit workers/native callback constructor passes the assigned resource
scope to the callback. A role's context is an optional display with
`native_source_namespace` and `native_cutoff`; these identify the actual retained
producer, rather than a fresh record made from the display.

For `claims_from`, a prospective role supplies `claims` using
`gadget.claims.native.v1`: `source_namespace`, `cutoff` and `items`. The checker
requires the exact producer pin before comparing any proposed literal. Wrong
units, paths or interpretation cannot be repaired by inventing a measurement.
Scientific units and endpoint semantics must be selected from actual source
metadata; no scalar or field dimension establishes them automatically.

The retained native output can be moved or reopened independently of the source
folder. Reopen each `publication` archive to select its tensors or replay the
producing operation. The outer workflow is orchestration metadata, with its
original study and optional evidence publication; it does not merge independent
mathematical carriers. Local models supply prospective controls and explanations,
and native tools own the mathematical observations and results.
