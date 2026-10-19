# Examples

Choose a guide for the input data or question. Each method declares the
measurements, source identities and comparison maps it needs. Shared labels
alone do not establish that two measurements can be compared.

## Data and assays

| Input or task | Guide |
| --- | --- |
| BindingDB, PubChem and GTEx source inspection | [Source format readers](source-formats/README.md) |
| Tables, matrices, intervals and assay workspaces | [Bioinformatics](bioinformatics/README.md) |
| Sequence alignment, variants and assembly | [Sequence and assay methods](bioinformatics/suite.md) |
| RNA evidence and expression assignment bounds | [RNA evidence](bioinformatics/rna-evidence/README.md) |
| Matched expression fields | [Paired contrasts](transcriptomics/README.md) |
| ChIP and ATAC coverage | [Controlled coverage](bioinformatics/controlled-coverage/README.md), [EGR1 example](bioinformatics/chip-coverage.md) |
| Spatial expression | [Spatial data](bioinformatics/spatial-expression.md) |
| Correspondence between assays | [Assay sections](bioinformatics/assay-sections/README.md) |
| Native count fields, exact pathways and evidence maps | [Native analysis](source-formats/ANALYSIS.md) |
| Biomedical databases and binding evidence | [Biomedical evidence](biomedical/README.md) |

## Analyses

| Question | Guide |
| --- | --- |
| Fit a declared linear model or constrained response | [Linear models](bioinformatics/linear-models/README.md), [constraints](bioinformatics/constrained-responses/README.md) |
| Inspect Hodge fields or Green responses | [Hodge fields](bioinformatics/hodge-frames/README.md) |
| Compare complete field realizations | [Field ensembles](bioinformatics/relational-ensembles/README.md) |
| Check consistency across local views | [Relational gluing](bioinformatics/relational-gluing/README.md) |
| Calibrate a statistic under an explicit finite law | [Calibration](bioinformatics/finite-calibration/README.md) |
| Compare declared structural omissions | [Ablations](bioinformatics/structural-ablations/README.md) |
| Analyze longitudinal measurements and censored endpoints | [Longitudinal analysis](bioinformatics/longitudinal-endpoints/README.md) |
| Follow changes in topology | [Homology diagrams](bioinformatics/homology-diagrams/README.md) |

## Execution and storage

- [Pipelines](pipelines/README.md) connect dependent assay operations.
- [Scientific sequences](scientific/README.md) reuse one source across ordered queries.
- [Native data](native-data/README.md) covers C++ records, original bytes and replay.
- [Python data I/O](../tools/rcdb/DATA.md) covers bounded intake and selected fields.
- [Source contexts](bioinformatics/source-context/README.md) connect queries and
  annotations to original records.

Build examples with `GADGET_BUILD_EXAMPLES=ON`; individual guides name their
required targets. Downloaded data and local models have separate requirements.
New experiments require new output directories.
