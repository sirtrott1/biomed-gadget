# Chromatin fragment coverage

This example accounts for every published fragment on chromosomes 2 and 3 of
the 10x Genomics Next GEM v1.1 PBMC library. It preserves each interval, barcode
and reported read pair multiplicity. Chromosomes are complete analytical units;
rows and barcodes are not sampled or filtered. The complete publisher archive
is retained with the native workspace.

The inputs are the [published PBMC library](https://www.10xgenomics.com/datasets/500-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-1-1-0)
and matching [UCSC hg19 chromosome FASTAs](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/).
The preparation script checks pinned SHA256 hashes before using them. The
library is distributed under CC BY 4.0; UCSC supplies the reference sequences
for public use. Preserve the source attribution when sharing data or results.

## Run

Build GADGET with examples enabled and its separate RexGraphCpp dependency as
described in the root [README](../../README.md). This example needs the CPU
profile and Python 3. It does not require HTSlib, a model runtime or a GPU.
Run these commands from the GADGET repository:

```sh
cmake --preset cpu
cmake --build --preset cpu --parallel 1
python3 examples/bioinformatics/prepare_chromatin.py \
  data/chromatin --gadget build/presets/cpu/gadget --download
build/presets/cpu/gadget import-assays \
  --spec data/chromatin/import.rfcxs --out data/chromatin/imported
build/presets/cpu/gadget assay-query \
  --workspace data/chromatin/imported/workspace.rfcxs \
  --request data/chromatin/coverage.rfcxs \
  --out-record data/chromatin/coverage-result.rfcxs
build/presets/cpu/gadget_verify_coverage_fields \
  data/chromatin/coverage-result.rfcxs data/chromatin/input/fragments.bed \
  data/chromatin/imported/workspace.rfcxs data/chromatin/coverage.rfcxs
```

Use new destination paths for new publications. To reuse a downloaded complete
fragment archive, pass `--fragments PATH` to the preparation script. Its hash is
checked against the same source. The script writes native RFCXS declarations;
temporary JSON used by the authoring bridge is removed after conversion.

To admit the coverage source and its fields for subsequent scientific queries:

```sh
build/presets/cpu/gadget derive-assays \
  --workspace data/chromatin/imported/workspace.rfcxs \
  --request data/chromatin/materialize.rfcxs --out data/chromatin/derived
build/presets/cpu/gadget_verify_coverage_fields --assay fragment-coverage \
  data/chromatin/derived/workspace.rfcxs data/chromatin/input/fragments.bed \
  data/chromatin/imported/workspace.rfcxs data/chromatin/materialize.rfcxs
build/presets/cpu/gadget replay-derivation \
  --workspace data/chromatin/derived/workspace.rfcxs \
  --assay fragment-coverage --out data/chromatin/replayed
build/presets/cpu/gadget_verify_coverage_fields --assay fragment-coverage \
  data/chromatin/replayed/workspace.rfcxs data/chromatin/input/fragments.bed \
  data/chromatin/imported/workspace.rfcxs data/chromatin/materialize.rfcxs
```

The verifier selects the actual retained producer publication when reading an
annotated study. It does not substitute a later endpoint. For a partitioned
assay, verify each standalone publication with its corresponding original input.

## Source and field construction

Every original fragment is a grade one relation. The grade zero identities are
the physical atoms formed by interval endpoints within each declared contig and
strand frame. A relation retains one ordered occurrence for every atom it
contains. Identical interval coordinates remain distinct original relations.
Overlapping coordinates on different chromosomes never share a physical atom.

The boundary has coefficient one on those occurrences. Reported multiplicities
are a grade one field, rather than structural coefficients. RexGraphCpp's
`BoundaryOperator::apply` produces the exact grade zero coverage. No pair graph
or adjacency matrix is constructed. Arbitrary relation arity remains available
in the native source.

The metric is the exact physical length of each atom. For coverage `c` and
lengths `l`, the integral is `sum(l[i] * c[i])` and the quadrance is
`sum(l[i] * c[i] * c[i])`. These operations use RexGraphCpp rational arithmetic,
field actions and metric services. The request has an identity scale divisor
and an empty control span. Therefore its residual equals its measured field;
it is not a difference between treatment groups.

Native source fields retain interval and atom contigs, strands and coordinates,
original row indices and byte spans, barcode spans, multiplicities, support and
contrast fields. The original BED, reference FASTA, compressed publisher inputs,
generating workspace and request remain byte bound to the publication.

The observed union determines the field domain. Within that domain a measured
zero remains supported. Outside it, observations are absent. Neither state is
replaced by an implicit genome wide zero. This example does not assert higher
attaching maps from interval overlap.

## Complete verification

`verify_coverage_fields.cpp` is an independent qualifier for identity fragment
coverage with five column BED input. It checks every original row and barcode,
field coordinate, physical atom, relation and ordered occurrence. Its separate
event sweep checks exact coverage, metric, integral and quadrance on every
contig. It also checks complete measured field support, the empty control span,
threshold declarations, retained workspace/request bytes and complete canonical
archive reencoding. Its stdout tables are views of those checked results.

The qualifier accepts caller supplied contigs, intervals and thresholds. The
preparation adapter pins this particular published example; it does not restrict
the production operations to these chromosomes or this library. Controlled
comparisons with multiple samples, explicit rational scales and declared control
fields use the general `controlled_coverage` operation and its projection
certificates. The identity verifier rejects those different question contracts.

Reported read support is distinct from unique fragments and from cell counts.
All published barcodes are included, including those outside a filtered cell
set. An integral measures reported support across bases; a quadrance gives
greater weight to concentrated coverage. These quantities describe the supplied
mapped fragments. A treatment comparison requires its actual samples and design.

Portable operation files retain their complete generating context. Separate
publications can repeat that context. Their size is not a claim of minimum
storage or of streaming admission.
