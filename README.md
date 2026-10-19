# biomed-gadget

biomed-gadget is a C++ library and command line application for biomedical data
analysis. It connects measurements, biological context and mathematical models
while keeping the original data and the steps that produced each result.

It is built for analyses that need to combine sources, inspect assumptions or
revisit results as evidence changes.

## Working with biomedical data

- Import count matrices, tables, sequences, genomic intervals and assay records.
- Analyze expression, spatial measurements, chromatin and declared relationships
  across assays.
- Run linear models, Hodge analysis, calibration, model comparisons and temporal
  queries with explicit inputs and assumptions.
- Store versioned data, retrieve selected measurements and replay analyses.
- Compose pipelines, with optional local models to propose queries and read
  checked results.

## What is different

**Relations keep all their participants.** Through RexGraph, a relation can join
several records without being replaced by pairwise edges. Repeated participants
and their order remain available for later queries.

**Results keep their scientific context.** Values retain original coordinates,
units and source mappings. Measured zero and missing observations stay distinct.
Comparing assays requires a declared correspondence, rather than matching array
shapes or labels alone.

**Supported exact models retain alternatives.** When evidence permits several
compatible solutions, an analysis can retain the complete solution family and
its constraints. Later evidence can refine that family through explicit maps.

**Stored results remain usable data.** Measurements and results can be queried
and passed into further calculations with their sources intact. Original files
remain provenance; working data use native records that can be read in parts.

[RexGraphCPP](https://github.com/sirtrott1/rexgraph-cpp) supplies the mathematics
and native formats. GADGET adds biomedical declarations, adapters and workflows.
The optional [Python stack](https://github.com/sirtrott1/rexgraph-python)
provides data access, RCDB storage and RCQL queries.

## Get started

The Linux CPU build needs a C++20 compiler, CMake 3.24 or newer, Ninja,
`pkg-config`, GMP, `json-c` and OpenSSL development packages.
[Build instructions](cmake/README.md) include system packages and optional
formats. The [dependency guide](DEPENDENCIES.md) gives tested revisions and the
complete C++ and Python installation.

```sh
git clone https://github.com/sirtrott1/rexgraph-cpp.git rexgraphcpp
git clone https://github.com/sirtrott1/biomed-gadget.git biomed-gadget
cd biomed-gadget
cmake --preset cpu -DGADGET_BUILD_EXAMPLES=ON
cmake --build --preset cpu --target gadget gadget_scientific_example
export PATH="$PWD/build/presets/cpu:$PATH"
```

### First analysis

Generate a small synthetic study, run a coverage query and inspect its stored
result. Output paths must be new.

```sh
gadget_scientific_example build/demo/study
gadget native-study --study build/demo/study/study.json --out build/demo/study.rfcxs
gadget query --study build/demo/study.rfcxs \
  --query examples/scientific/coverage.rfcxs --out build/demo/coverage.rfcxs
gadget record-view --input build/demo/coverage.rfcxs --member result
```

Use `gadget --help` to inspect commands. For real inputs, start with the
[examples](examples/README.md) or the [Python data guide](tools/rcdb/README.md).
Large native calculations can still require a complete source in memory;
[resource controls](cmake/README.md#memory-and-scratch-space) describe the limits.

## Guides

| Task | Start here |
| --- | --- |
| Load, explore and select data from Python | [Data I/O](tools/rcdb/README.md) |
| Choose an assay or analysis | [Examples](examples/README.md) |
| Build a pipeline | [Pipelines](examples/pipelines/README.md) |
| Use GADGET as a C++ library | [C++ interface](include/gadget/README.md) |
| Enable formats or local inference | [Build options](cmake/README.md#optional-features) |
| Run tests or contribute | [Testing](tests/README.md) |

## License

GADGET is licensed under [Apache 2.0](LICENSE).
Copyright 2026 biomed-gadget. See [NOTICE](NOTICE).
Dependencies and datasets retain their own terms.
