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

`has_native_scientific_publisher` reports the implemented scientific inventory;
the selected source, carriers and query still require admission. Biomedical
leaves use native role/composite evidence operations. A query without a native
publisher refuses during admission, including in an untaken branch or unused
definition. `execute_workflow` provides the compatibility display through the
same Runner.
`workflow_record_view` displays a saved native workflow without executing it.
Generic record views expose the structural receipt and original archive refs.
The CLI uses that strict route when `--out` requests a native archive. Without
`--out`, the CLI displays the compatibility Runner result from the original
typed declaration. Mathematical leaves retain their actual native owners.

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
`LocalAgent::ask_record` and `ask_assigned_record` retain scientific query and
plan results directly, together with the original proposal and model receipts.
The CLI's local roles use this native interface. A workflow forwards the actual
query or complete ordered plan to dependent work and retains the original model
response in its trace. Every plan child must match its declared query and original
study owner. Replay of retained calculations requires no new model generation.
The explicit workers/native callback constructor passes the assigned resource
scope to the callback. A role's context is an optional display with
`native_source_namespace` and `native_cutoff`; these identify the actual retained
producer, rather than a fresh record made from the display.

For `claims_from`, a prospective role can supply `claims` using
`gadget.claims.native.v1`: `source_namespace`, `cutoff` and `items`. A local
`claim_reader` can instead supply its admitted prospective `claims` with the
original `source_digest`, `certificate_sha256` and `items`. The checker validates
that complete certificate binding and resolves the actual native producer pin
before comparing any proposed literal. Wrong
units, paths or interpretation cannot be repaired by inventing a measurement.
Scientific units and endpoint semantics must be selected from actual source
metadata; no scalar or field dimension establishes them automatically.

The retained native output can be moved or reopened independently of the source
folder. Reopen each `publication` archive to select its tensors or replay the
producing operation. The outer workflow is orchestration metadata, with its
original study and optional evidence publication; it does not merge independent
mathematical carriers. Local models supply prospective controls and explanations,
and native tools own the mathematical observations and results.

A model context reference given as a string displays the complete mathematical
result. To display only named native literals, use an explicit selection:

```json
{"from":"measured","selections":[{"path":["reported"],"value_type":"count"}]}
```

Use this object as an item in a model node's `context` array. The dependency must
own a completed scientific field publication. A selection supports `count`,
`integer`, `rational`, `text` or `boolean`. A tensor selection can declare an
original flat `coordinate`; it must be supported in that original field.
Supported zero remains zero, and a missing coordinate refuses the selection.
Exact rational values use rational strings. Paths, types, coordinates and actual
native bindings accompany every displayed value.

The context marks `result_complete` false and supplies `selections` in place of
`result`. Its certificate and result digests still identify the complete original
result. Its native namespace and cutoff still name the actual producer. The
workflow retains that complete publication independently of the display choice.
Ask a reader to use the displayed selection's original path for any proposed
claim. `check_claims` compares the claim against the complete native dependency.
A display selection cannot change source support, create a new measurement or
reduce the mathematical calculation.
