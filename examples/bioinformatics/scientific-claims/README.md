# Checking source facts and analytical results

A source citation can be valid while the sentence attached to it is false.
GADGET checks literal claims against selected source descriptors and completed
analytical outputs. Researchers select each fact's path, type and scientific
meaning. Models can report values or comparisons; they cannot change those
selections.

This example reuses the RNA assignment and source-context protocols. It checks
15 facts: reported and missing RNA support, declared units, individual positive,
negative and inactive incidence coefficients, a transcript length, recorded read
quality, measured RNA and computed allocation contrast bounds. GeneA's assigned
read count increases by exactly one under every retained allocation. The A1
contrast remains unresolved. These are finite analytical statements about the
supplied demonstration inputs.

Run with a separately built GADGET executable:

```sh
python3 examples/bioinformatics/scientific-claims/replay.py \
  --gadget /absolute/path/to/gadget \
  --out build/scientific-claims-example
```

The driver stages the existing used inputs, executes the authored methods, checks
all selected facts and deliberately submits four incorrect readings. The result
receipts retain the incorrect proposals, actual observed values and contradictions.
It then moves the complete source and output directory and checks that the same
catalog and claims can be reconstructed. The RNA derivation also replays from its
retained source inputs.

`facts.json` uses `gadget.scientific.facts.v1`. A source selection names an admitted
context target and a descriptor path. An answer selection names a completed answer,
binding, producing step, result ordinal and path within that result. Answer file
hashes, result file hashes, certificate seals, source owners and effective queries
are checked before values enter the catalog. Paths resolve against the selected
objects; no claim-count, path-depth or array-index policy truncates the selection.
The admitted C++ object owns its captured values independently of subsequent file
edits. Opening a new object checks the files again.

The commands are `fact-catalog`, `claim-proposal` and, in a build with llama.cpp,
`read-facts`. For the latter, use the question recorded in `facts.json` as the
prompt, with a local GGUF and explicit context/output resources:

```sh
gadget read-facts --workspace WORKSPACE/workspace.json --facts facts.json \
  --artifacts RESULTS_ROOT --model LOCAL.gguf --prompt QUESTION \
  --threads 8 --context-tokens 16384 --output-tokens 2048
```

Inspect the same reading view without loading a model:

```sh
gadget fact-catalog --workspace WORKSPACE/workspace.json --request facts.json \
  --artifacts RESULTS_ROOT --view model
```

The version 2 view places the selected values and their types together and states
the reporting task separately from the research question. It preserves every
fact's meaning, source selector, path and seal. Repeated owners appear once; each
answer owner references its complete retained query by its hash. The full queries,
original observations and numerical certificates remain in the catalog and every
checked receipt. Owners are shared only when their complete captured bytes agree.
To give the model the full catalog instead, pass `--catalog-mode full`.

The reporting decoder provides one slot per selected fact, in selection order.
Its size follows the request; it has no fixed fact-count ceiling. Each slot fixes
the selected identity and equality comparison while leaving the reported value
to the model. Rational strings follow the native canonical integer/fraction
syntax. The backend independently rejects invalid values and checks every
comparison. The actual decoder schema and grammar hash are retained with the
model response.

Researchers can create a focused reading request by selecting fewer fact paths in
`facts.json`. Its catalog has its own seal, and the reader must report every fact
in that selection. The receipt makes that scope explicit; it does not claim that
a focused selection answers every part of the original research question.

Receipts distinguish `matched`, `contradicted`, `incomplete` and `selection_required`.
Each comparison uses the fact's declared type: exact rational ordering, integer
counts, signed integers, Boolean equality or literal text equality. Missing coordinates cannot be
read as numeric zero. Equal text does not identify cells in different native
sources or supply a chain map. Agreement verifies the stated literal comparison;
biological interpretation and the researcher's mapping from the question to these
facts remain separate. The complete captured evidence and original model proposal
remain available for review.

A model clarification retains its original reason and marks the selected facts as
unreported. The catalog already contains every selected value; a model
clarification does not establish that the data is unavailable. Rational claim
values must be literal strings, rather than objects carrying additional assertions.
