# biomed-gadget

biomed-gadget is a C++ library and application for biomedical analysis using
relational complexes. It brings measurements, biological relations and source
evidence into the same mathematical structure, so that an analysis can operate
on their relationships as well as their individual values.

A relation can involve several genes, specimens, molecules or observations.
Its participants, order, orientation and coefficients remain part of the
structure. Measurements are fields on that structure, with their original
coordinates, units and observation support. This allows a calculation to
distinguish a measured zero from an absent observation, or a shared label from
an established correspondence between assays.

The mathematics comes from [RexGraphCPP](https://github.com/sirtrott1/rexgraph-cpp).
GADGET supplies the biomedical declarations, data admission, analyses and
execution interfaces. The independent
[RexGraph Python library](https://github.com/sirtrott1/rexgraph-python) provides
Core, RCDB storage and RCQL queries for the optional Python data tools. Its
`rexgraph_cpp` bridge connects those tools to the C++ library.

## The library

The [scientific interface](include/gadget/scientific.hpp) supports linear
models, constrained and observed responses, Hodge fields, finite calibration,
model comparisons and temporal calculations. Where a model admits a family of
solutions, the native result can retain the complete coefficient family and
identified contrasts. Source complexes, metrics and maps remain available when
comparing results or refining a model with further evidence.

Conventional statistical methods live in [statistics](include/gadget/statistics/README.md).
These include Poisson and negative binomial count models, normalization,
dispersion estimation, Wald and likelihood ratio tests, effect shrinkage and
multiple testing correction. Each method has its own assumptions and numerical
diagnostics. [Pathway analysis](examples/source-formats/ANALYSIS.md) also supports
an exact finite reference law with an explicit feature universe and mapping.

[Native data records](tools/rcdb/DATA.md) store working measurements and shared
context separately from the original files retained as provenance. Selections
carry their source versions and coordinate maps into subsequent calculations.
[Pipelines](examples/pipelines/README.md) connect these operations and retain
their native results. Optional [local models](examples/workflows/README.md) can
propose queries and read checked results through the same execution interfaces.

The directory guides cover the available methods and their inputs:

| Area | Documentation |
| --- | --- |
| C++ library | [Interfaces and CMake targets](include/gadget/README.md) |
| Data import, selection and storage | [Python data tools](tools/rcdb/README.md), [C++ native records](examples/native-data/README.md) |
| Expression, chromatin, spatial and sequence data | [Assay workspaces](examples/bioinformatics/README.md), [sequence and assay methods](examples/bioinformatics/suite.md) |
| Source formats and biomedical evidence | [BindingDB, PubChem and GTEx](examples/source-formats/README.md), [binding evidence](examples/biomedical/README.md) |
| Regression and response models | [Linear models](examples/bioinformatics/linear-models/README.md), [observations](examples/bioinformatics/observed-responses/README.md), [constraints](examples/bioinformatics/constrained-responses/README.md) |
| Fields and relations | [Hodge fields](examples/bioinformatics/hodge-frames/README.md), [gluing](examples/bioinformatics/relational-gluing/README.md), [field ensembles](examples/bioinformatics/relational-ensembles/README.md) |
| Statistics | [Count models](include/gadget/statistics/COUNTS.md), [pathways and evidence maps](examples/source-formats/ANALYSIS.md) |
| Temporal data and learned patterns | [Field histories](examples/bioinformatics/versioned-history/README.md), [native patterns](examples/native_patterns/README.md) |
| Execution | [Scientific query sequences](examples/scientific/README.md), [pipelines](examples/pipelines/README.md), [workflows](examples/workflows/README.md) |

The [examples index](examples/README.md) links to the individual assay and
analysis guides.

## Installation

The Linux CPU build requires a C++20 compiler, CMake 3.24 or newer, Ninja,
`pkg-config`, GMP, `json-c` and OpenSSL development packages. Clone RexGraphCPP
beside GADGET, then build the application and the introductory example:

```sh
git clone https://github.com/sirtrott1/rexgraph-cpp.git rexgraphcpp
git clone https://github.com/sirtrott1/biomed-gadget.git biomed-gadget
cd biomed-gadget
cmake --preset cpu -DGADGET_BUILD_EXAMPLES=ON
cmake --build --preset cpu --target gadget gadget_scientific_example
export PATH="$PWD/build/presets/cpu:$PATH"
```

The [build guide](cmake/README.md) covers system packages, installed CMake
targets and optional formats. [DEPENDENCIES.md](DEPENDENCIES.md) records tested
revisions and installation of the Python packages and bridge. The C++
application can be built independently of those Python tools.

## First analysis

The introductory example creates a synthetic study with its source declarations.
Publish it as a native study, execute the supplied coverage query and inspect
the stored result. Each output path must be new.

```sh
gadget_scientific_example build/demo/study
gadget native-study --study build/demo/study/study.json --out build/demo/study.rfcxs
gadget query --study build/demo/study.rfcxs \
  --query examples/scientific/coverage.rfcxs --out build/demo/coverage.rfcxs
gadget record-view --input build/demo/coverage.rfcxs --member result
```

`gadget --help` lists the commands. For an existing dataset, the
[data guide](tools/rcdb/README.md) covers import, metadata exploration and selected
field retrieval. Some native calculations require the complete source carrier
in memory; the [resource controls](cmake/README.md#memory-and-scratch-space)
describe memory budgets and scratch storage. Build and test instructions are in
[cmake](cmake/README.md) and [tests](tests/README.md).

GADGET is licensed under [Apache 2.0](LICENSE).
Copyright 2026 biomed-gadget. See [NOTICE](NOTICE).
Dependencies and datasets retain their own terms.
