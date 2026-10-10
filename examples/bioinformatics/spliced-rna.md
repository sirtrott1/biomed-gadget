# Alignment, splicing and assignment walkthrough

This walkthrough reads both complete chromosome 4 BAM files in
[pasillaBamSubset 0.48.0](https://bioconductor.posit.co/packages/3.22/data/experiment/html/pasillaBamSubset.html).
They contain untreated single end and paired end observations aligned with
TopHat 1.2.0. The package supplies chromosome 4 of the dm3 reference. The preparation
script checks it against the complete UCSC dm3 FASTA and retains every chromosome
4 record from the [UCSC ensGene annotation](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/genes/).
The published BAMs are chromosome subsets. No records within those files are
sampled. The original archive, compressed reference and annotation remain in
the data directory with their published terms and descriptions.

## Build and prepare

Install HTSlib 1.20 or later and the dependencies in the main README. This example
uses the formats preset with HDF5 disabled because its inputs are BAM, GTF and
FASTA. RexGraphCpp is imported from its separate checkout or installed package.

```sh
cmake --preset formats -DGADGET_ENABLE_HDF5=OFF
cmake --build --preset formats --target gadget gadget_alignment_fields_example
python3 examples/bioinformatics/prepare_spliced_rna.py \
  --gadget build/presets/formats/gadget --data ./rna-seq --download
build/presets/formats/gadget import-assays \
  --spec ./rna-seq/import.rfcxs --out ./rna-seq/workspace
```

Python 3.12 or later prepares inputs and performs independent verification.
Scientific operations run in C++. Requests and publications use RFCXS.
The script checks complete input hashes and generates native import and query
requests. Existing native requests are preserved; use a fresh data directory
for another run. All downloads have fixed byte hashes. A changed upstream file
fails preparation.

## Alignment evidence

For either `untreated1` or `untreated3`, run:

```sh
build/presets/formats/gadget assay-query \
  --workspace ./rna-seq/workspace/workspace.rfcxs \
  --request ./rna-seq/untreated1-alignment.rfcxs \
  --out-record ./rna-seq/untreated1-alignment-result.rfcxs
build/presets/formats/gadget_alignment_fields_example \
  ./rna-seq/untreated1-alignment-result.rfcxs ./rna-seq/untreated1-alignment
python3 examples/bioinformatics/verify_alignment_fields.py \
  --bam ./rna-seq/originals/pasillaBamSubset/inst/extdata/untreated1_chr4.bam \
  --prefix ./rna-seq/untreated1-alignment \
  --exclude-flags 3588 --minimum-mapq 0 --unknown-mapq include
```

The native complex has one grade zero identity per original alignment record.
Each aligned CIGAR block is a primary grade one relation with its record as its
boundary. Each distinct stranded CIGAR `N` junction is another primary relation;
its ordered occurrences point to every contributing record. All occurrences
remain present when a filter excludes a record. The selected count is the
canonical boundary transpose applied to the record selection field.

Typed fields retain flags, mapping quality, coordinates, query length and record
selection. Unmapped start and end coordinates are absent. A mapped start at zero
remains a supported zero. Original BAM bytes, record indices and decoded query
names accompany the fields. These are record and junction incidence relations;
the construction supplies no higher attaching maps.

The example writes complete record and relation tables. The verifier independently
decodes every BAM record and compares every coordinate, block, junction count and
ordered boundary occurrence. It performs no alignment or abundance fitting.

## Gene and transcript assignment

Replace `alignment` with `gene` in the query and output names, then verify:

```sh
build/presets/formats/gadget assay-query \
  --workspace ./rna-seq/workspace/workspace.rfcxs \
  --request ./rna-seq/untreated1-gene.rfcxs \
  --out-record ./rna-seq/untreated1-gene-result.rfcxs
build/presets/formats/gadget_alignment_fields_example \
  ./rna-seq/untreated1-gene-result.rfcxs ./rna-seq/untreated1-gene
python3 examples/bioinformatics/verify_feature_counts.py \
  --bam ./rna-seq/originals/pasillaBamSubset/inst/extdata/untreated1_chr4.bam \
  --annotation ./rna-seq/input/chr4.ensGene.gtf --prefix ./rna-seq/untreated1-gene \
  --attribute gene_id --overlap any --strand-policy ignore \
  --exclude-flags 3588 --minimum-mapq 0 --unknown-mapq include
```

For transcript queries, use `transcript` in the filenames and the verification
options `--attribute transcript_id --overlap all --exon-chain`. The transcript
rule requires all reported aligned bases inside the transcript exon union and
every crossed intron to match its declared exon chain. Deletions and skipped
introns contribute no aligned bases.

Each sequenced segment is identified by read group, query name and mate. Its
complete retained alternatives determine its candidate set. Alternative mappings
outside the annotation retain an unassigned outcome. The primary allocation
relations contain these candidate sets, including arities greater than two.
Class counts and lower and upper feature bounds are native fields on this complex.
The segment table maps each original observation to its class; overlap tensors
retain every base count and eligibility witness.

A lower bound counts segments assigned to a feature under every admissible
allocation. An upper bound counts segments that can be assigned to it. Equal
bounds identify the count under the specified evidence. Different bounds expose
unresolved assignment. The bounds are marginal; endpoints for different features
need not be attained together. Paired ends remain separate segments, and gene
and transcript queries answer different assignment questions.

The default requests include mapping quality 255, which means unavailable.
They exclude unmapped, failed quality, duplicate and supplementary records by
the explicit SAM flag mask 3588. Change these declarations for another research
question. The native paths retain the original inputs and generating workspace;
their portable archives include that complete context. Import catalogs contain
validated counts and headers. Use `inspect-file` for a complete compatibility
inspection.

## Native C++ access

`AssayWorkspace::execute_record` returns an operation publication. Its `result`
member contains the certificate; that certificate's `result` member contains the
field declaration. The operation endpoint owns the primary complex and fields.
Read them through the publication that owns their keys.

| Result schema | Field declaration | Source carrier |
|---|---|---|
| `gadget.alignment.fields.v1` | `record_fields` | Grade zero alignment records |
| `gadget.alignment.fields.v1` | `relation_fields` | Grade one blocks and stranded junctions |
| `gadget.feature.count.fields.v1` | `fields.class_count` | Grade one candidate classes |
| `gadget.feature.count.fields.v1` | `fields.lower`, `fields.upper` | Grade zero candidates |

Field declarations contain integer keys identifying endpoint tensors. Candidate
and class dictionaries preserve the ordering of their source coordinates.
Assignment `segments` and `overlaps` contain native tensors for observation
fibres and interval witnesses. Query name byte pools have offset and length
columns. The complete implementation of this reader is in
[alignment_fields.cpp](alignment_fields.cpp); it uses no alignment decoder.

The native record also retains original byte spans and the complete generating
workspace and request publications. `render_record_view` reconstructs the
established compatibility result when required. The alignment view needs HTSlib
to render original SAM text. Reading native field tables requires only the
ordinary CPU build.
