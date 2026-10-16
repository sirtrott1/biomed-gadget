# Scientific annotations and reproducible RNA answers

This example links a research question to retained source records, then runs the
existing RNA expression, read quality and measurement support methods. Source
references identify RNA features, quantitative fields, the whole assay context,
individual relation occurrences, transcript records and a FASTQ read. All inputs
are constructed method examples. The transcript to gene declaration belongs to
the read evidence protocol; the matrix retains its separate native namespace.

```sh
python3 examples/bioinformatics/source-context/replay.py \
  --gadget build/native/gadget --out build/source-context-example
```

The driver stages only the files used by this protocol from `rna-evidence/` and
`staged-assays/`, together with `relations.json`. Its output includes those seven
original inputs, their import declaration, a portable native workspace, the
source context request, annotation proposals and receipts, a typed pipeline
intent and the completed numerical results. It relocates the workspace and
results, checks identical source references, and regenerates the expression
analysis from its retained lineage. Reproduction needs the GADGET executable and
these public example files; it does not read private notes or earlier archives.

The native RNA fields preserve measured zeros and unobserved coordinates. The
reported field has two supplied coordinates; the unobserved field has none.
The authored relation retains two distinct occurrences of feature `f1` with
coefficients −1 and +1, and an inactive `f2` occurrence whose realized coefficient
is zero. Annotation references retain all three occurrences. Their cancellation
in the boundary leaves their source identities and roles available for inspection.

Read evidence analysis retains every admissible assignment to the supplied
transcripts. The declared geneA total increases by one assigned read under every
admissible allocation. The individual transcript changes remain unresolved.
The quality method reports the stored baseline quality symbols. The support
methods report the matrix's supplied coordinates. `answers.json` links each
requested output to its complete producing certificate and relevant source
mentions.

A source context request uses `gadget.assays.source.context.v1`. Each selected
native target names its assay, partition and study digest, together with its
cell, occurrence, field or context selector. A physical resource target names
the original artifact hash and, optionally, its decoded record index. Evidence
quotes use UTF8 byte offsets. Resource record quotes must lie within that
record's original span. The context catalog includes resolved owners, native
field coordinates and support, incidence roles, source hashes and verified spans.

```sh
gadget source-context --workspace WORKSPACE/workspace.json --request context.json
gadget annotation-proposal --workspace WORKSPACE/workspace.json \
  --context context.json --request annotations.json
```

A local model can propose the annotation step using the same context and exact
question. Replace `QUESTION` below with `context.json`'s question, unchanged.
The model selects declared target IDs, writes an interpretation and states its
polarity. Owners, selectors and evidence offsets come from the admitted context.

```sh
gadget annotate-source --workspace WORKSPACE/workspace.json \
  --context context.json --model /path/to/local-model.gguf --prompt 'QUESTION' \
  --threads 8 --context-tokens 16384 --output-tokens 1024 --gpu-layers 0
```

Annotation receipts mark interpretations as prospective and record
`interpretation_verified: false`. Byte verification establishes which record
was quoted. Evaluating the interpretation remains a scientific task. A
researcher can retain the proposal in a pipeline intent's optional
`source_context` declaration and explicitly map its mention IDs to typed answer
requirements. Planning and execution re admit that declaration; completed
answers retain the full mapping. The numerical methods use their original
parameters, fields, metrics and source declarations.

The original evidence and every selected target remain available. Model context
and output tokens are caller declared runtime resources. The annotation API has
no separate target, mention or evidence count ceiling.

Model responses retain the actual catalog shown to the model and its byte hash,
along with the system prompt hash. Compact annotation views keep every selected
target and quote while omitting the repeated request. Recipe views reference
sealed evidence and retain target descriptors, prospective mentions and the
researcher's mappings. Resource records share one decoded source during a request.
