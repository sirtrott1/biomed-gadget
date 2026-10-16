# Gluing declared relational structures

Two assays may describe overlapping parts of a shared structure. A source overlap
requires a declared common complex and two complete embeddings. Equal labels,
values or matrix shapes do not establish that overlap.

The `relational_gluing` query consumes RexGraph's `PreparedNativeGluing`. Its
pushout retains the common cells and occurrence fibres once, preserves each
participant's remaining structure, and checks the complete chain family. Scalar
sections must have complete support, the same quantity, units and channel, and
identical values at corresponding interface cells. Disagreement refuses gluing.
Use `assay_sections` for an explicitly declared metric compatibility analysis.

This constructed example contains an interface edge, two triangular fillings,
two three way relations, repeated cancelling occurrences, an inactive occurrence
and a supported zero occurrence. The right participant deliberately uses a
different coordinate order and different names. Its explicit map supplies the
correspondence. Grade zero measurements, grade one signals and grade two signals
are separate sections. Boundary lineage continues to grade zero when requested.

From the GADGET source directory, with the executable on `PATH`:

```sh
python3 examples/bioinformatics/relational-gluing/make_example.py --out inputs
gadget import-assays --spec inputs/import.json --out workspace
gadget assay-query --workspace workspace/workspace.json --request inputs/query.json
gadget pipeline-prepare --workspace workspace/workspace.json --request inputs/pipeline.json
```

The glued complex has five grade zero cells, seven grade one cells and two
grade two cells. It retains 20 grade one occurrences and six grade two occurrences.
The shared edge's scalar value is 2. The two face signals are 7 and 11; their
boundary contributes 18 on the shared edge because a chain boundary adds both
incident contributions. Applying the next boundary gives supported zeros. This
distinguishes pointwise section gluing from chain pushforward or boundary action.

The compatibility result includes original publication namespaces and cutoffs,
original byte artifact bindings, inclusion maps for every original cell and occurrence, retained
attaching positions and aspects, exact coefficients, field meanings and ordered
boundary certificates. Generated keys address the prospective glued complex;
the inclusion maps address original publications. The result does not publish a
new native history. Heterogeneous stalk restrictions and signature recodings have
their own contracts and are outside this scalar gluing query.

The original documents and native binary remain in the workspace. Gluing resolves
primary identities and field pins directly in the admitted source; it does not
rebuild or embed a JSON catalogue. `assay-source-read` remains a separate optional
presentation interface. The independent verifier reads every original evidence
span, checks the complete byte artifacts and independently reconstructs the
source maps, occurrence fibres and exact field/boundary values:

```sh
gadget assay-query --workspace workspace/workspace.json --request inputs/query.json > result.json
python3 examples/bioinformatics/relational-gluing/verify.py --result result.json --reference inputs/reference.json --workspace workspace/workspace.json
```

Independent sections share one immutable native preparation and execute through
RexGraph CPU reservations and joined tiles. `workers` controls the requested
reservation. Preparation dimensions come from the actual sources. The generator's
`--vertices` and `--repetitions` options vary source size and section count.
`--modules` repeats the entire authored structural diagram, with separate original
identities and full interface maps. For example:

```sh
python3 examples/bioinformatics/relational-gluing/make_example.py --out incidence-inputs --modules 1000 --repetitions 8
```

This diagram has 4,000 vertices, 7,000 grade one relations, 2,000 grade two fillings
and 26,000 retained occurrences after gluing. It includes 2,000 three way relations.
The reference reconstructs the complete source maps, occurrence fibres and all
section/boundary coordinates from the authored records. These constructed diagrams
test software behavior; they are not biomedical datasets. Example defaults and
chosen benchmark sizes impose no application limits.

C++ consumers can use `gadget::execute_native_gluing` from
`<gadget/native_gluing.hpp>` to receive the owning native plan and complete native
chain lineage directly. `native_gluing_endpoint` exports every retained chain as
a source grade tensor over the constructed relational complex. Persist it through
RexGraph's RFCXS codec and persist the original input square through RFCXA;
neither requires a JSON result. See [native data admission](../../native-data/README.md).
