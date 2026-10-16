# Reading analytical results with local models

This protocol provides three explicit reading selections: 15 source and RNA
answer values, three RNA allocation contrast values, and 12 chromatin structural
values. It uses the same fact admission and comparison methods for each selection.
The selected values, types, paths and owners are fixed before a model reads them.

Run with a separately built GADGET executable:

```sh
python3 examples/bioinformatics/fact-readers/replay.py \
  --gadget /absolute/path/to/gadget --out build/fact-reader-example
```

The replay stages the used RNA and weighted BED inputs, executes the authored
pipelines, checks the selected values, replays their derivations and moves the
whole experiment. Every reading view must remain identical after that move.
The resulting `readers.json` records each workspace, fact request and result root
relative to the experiment directory. Optional local GGUF weights are supplied
separately.

The focused RNA selection reports geneA contrast bounds of one and the unresolved
A1 allocation direction. The chromatin selection reports a physical residual
`[3/2, -3/2, 0]`, residual quadrance 18, and eight independent zero coverage
directions among the original interval variables. It also retains the native
Green response and gradient, curl and harmonic quadrances. These are authored
method examples; their carriers and metrics are described in the
[coverage protocol](../controlled-coverage/README.md).

The Green values are `3/26`, `-3/26` and zero on the three physical segments.
The interval field has gradient quadrance `8/5`, curl quadrance zero and harmonic
quadrance `207/5`. These values concern different declared carriers: the physical
segment field and the original interval variable field. Their owners and complete
queries remain in the receipts. A model does not equate these carriers or infer
a biochemical feedback loop from the zero coverage sector.

Use `gadget read-facts` with the chosen entry's workspace, fact request and
artifact root, a local GGUF and the literal question in that fact request. The
[claim protocol](../scientific-claims/README.md) gives the command syntax.
The reader must report every selected value. Wrong comparisons, omissions and
false clarification remain in the original response and are checked separately.
A focused receipt covers its declared selection, rather than every possible
answer to the research question.
