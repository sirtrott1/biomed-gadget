# Native context patterns

`gadget/native_pattern.hpp` publishes and executes the canonical
`rfc::learning::PreparedTransportPatternBank`. Training remains with the SDK.
The caller supplies the actual original source and complete evidence publications,
including the conditioning kernel, and a declaration of each input and output
frame's quantity and units. Publication checks those archives against the bank's
training provenance and retains the canonical bank codec bytes. Hash agreement
records the declared training inputs; it does not independently repeat training.

`NativePatternModel` admits the stored canonical bank once and owns it across
queries. Inputs name an actual scientific field and an explicit rational map from
its original coordinates into a bank input frame. Source digest, model digest,
quantities, units, map dimensions and required support must agree. A declared zero
coefficient requires no source observation. A nonzero coefficient cannot read
missing support. Output support includes known zero values.

The query schema is `schemas/pattern-query.schema.json`. The three operations are:

* `sparse_candidates`: separate complete candidate response families. An input
  outside an evidenced input domain has no response.
* `sparse_program`: ordered named stages connected by declared maps from original
  stable occurrence keys. Stages retain maps from the same bank parameters. Joint
  directions can cancel; unresolved input families remain open and block dependent
  stages.
* `sparse_gluing`: declared restrictions on original occurrence keys, with complete
  compatibility defects, a separating contradiction witness when inconsistent,
  and the conditioned parameter family, section and boundary response when
  consistent.

Response direction matrices use the SDK's reduced response bases. Programme
`parameter_directions` and `input_directions` use the original bank parameter
carrier. Gluing retains the map from conditioned parameters to that original
carrier. These coordinate systems have distinct identities even when their
sizes agree. Original bank input coordinates retain their actual basis carriers.

Run a query and replay its retained original owners:

```sh
gadget pattern-query --study study.rfcxs --model model.rfcxs \
  --request query.json --out result.rfcxs
gadget replay-pattern --record result.rfcxs --out replay.rfcxs
```

For a workflow, pass `--patterns roles.json`. A `gadget.pattern.roles.v1` file has
`roles` entries with `role`, `model` and `request`. Relative model and request paths
resolve against that file. A workflow `model` node selects the registered role;
its fixed admitted request executes through a native callback using the assigned
CPU resources. Local language model roles may be registered separately with
`--models`. Pattern results support native workflow conditions and claim checks.
The prompt never supplies numerical field values or executable source code.

`gadget_native_pattern_contracts` uses an explicitly synthetic training population
and the independent scientific source fixture. It checks original evidence,
observable versus conditioning uncertainty, missing support, frame and unit
refusals, shared unknown cancellation, exact contradictions, portable replay and
workflow ownership. Synthetic agreement is not biomedical model qualification.
