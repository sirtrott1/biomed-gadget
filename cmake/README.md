# Building GADGET

GADGET builds against the independent RexGraphCPP SDK, either from source or
from an installed CMake package. The C++ application and library can be built
without the Python data tools. [DEPENDENCIES.md](../DEPENDENCIES.md) covers
installation of the complete stack. Run the commands below from the GADGET root.

## Requirements

The Linux CPU build needs C++20, CMake 3.24 or newer, Ninja, `pkg-config`, GMP
with GMPXX, `json-c`, OpenSSL Crypto and threads. Use SDK 0.3.0 with the tested
revision in [RexGraph.lock.cmake](RexGraph.lock.cmake).

On Debian or Ubuntu:

```sh
sudo apt-get update
sudo apt-get install build-essential cmake ninja-build pkg-config git \
  libgmp-dev libjson-c-dev libssl-dev
```

## Build

With the SDK checked out at `../rexgraphcpp`:

```sh
cmake --preset cpu
cmake --build --preset cpu
build/presets/cpu/gadget --help
```

Use `-DREXGRAPH_SOURCE_DIR=/absolute/path/to/rexgraphcpp` for another location.
The nested layout `../rexgraph-cpp/rexgraphcpp` is also detected when the default
sibling is absent. Custom paths can be saved in `CMakeUserPresets.json`.

The CPU preset builds the application and data tools. Enable
`GADGET_BUILD_EXAMPLES=ON` for examples or `BUILD_TESTING=ON` for tests.
[tests/README.md](../tests/README.md) describes the test environment and suites.

## Installed SDK

After installing the SDK with its own CMake project:

```sh
cmake --preset installed -DCMAKE_PREFIX_PATH=/absolute/path/to/rexgraph-install
cmake --build --preset installed
cmake --install build/presets/installed --prefix /absolute/path/to/gadget-install
```

The [C++ guide](../include/gadget/README.md) covers downstream CMake targets.
Optional libraries installed outside system directories need their discovery
and runtime library paths too.

## Optional features

| Feature | CMake option | Dependency |
| --- | --- | --- |
| H5AD | `GADGET_ENABLE_HDF5=ON` | HDF5 1.12+ |
| SAM, BAM, CRAM, VCF and BCF | `GADGET_ENABLE_HTSLIB=ON` | HTSlib 1.20+ |
| HTTPS database retrieval | `GADGET_ENABLE_PUBLIC_DATABASES=ON` | libcurl |
| Local GGUF inference | `GADGET_ENABLE_LOCAL_MODEL=ON` | llama.cpp source |
| PNG input to models | `GADGET_ENABLE_MULTIMODAL=ON` | Local model support and libmtmd |
| Packed Vulkan operations | `GADGET_ENABLE_VULKAN=ON` | Vulkan and compatible SDK backend |

The `formats`, `databases`, `local-model` and `vulkan` presets extend `cpu`.
Combine features with explicit options when needed. Exact rational operations
use CPU execution; the Vulkan adapter has its own supported operations.

For HDF5 and HTSlib on Debian or Ubuntu:

```sh
sudo apt-get install libhdf5-dev libhts-dev
cmake --preset formats
cmake --build --preset formats
```

For HTSlib in a separate installation prefix:

```sh
GADGET_HTSLIB_PREFIX=/absolute/path/to/htslib-install
export PKG_CONFIG_PATH="$GADGET_HTSLIB_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
cmake --preset formats
```

The runtime loader must also find its shared libraries. On Linux, use a
configured loader path or set `LD_LIBRARY_PATH` to that prefix's `lib` directory.

For local inference:

```sh
cmake --preset local-model -DLLAMA_SOURCE_DIR=/absolute/path/to/llama.cpp
cmake --build --preset local-model
```

Model files and inference settings are supplied separately. See
[workflow examples](../examples/workflows/README.md) and
[model tests](../tests/README.md#optional-profiles).

## Memory and scratch space

Compilation and linking default to one process each. Increase
`GADGET_COMPILE_JOBS` only after measuring compiler memory use; a larger build
`--parallel` value alone does not raise the compile limit.

On Linux 4.7 or newer, the CLI limits additional data allocations to half the
currently available memory, accounting for cgroups. Use `--memory-mib N` or
`GADGET_MEMORY_MIB=N` to choose a budget within that reserve.

The allocation limit covers heap and private writable mappings. Read only file
mappings and shared memory are outside it. Dependency allocators can terminate
a process on allocation failure. Library callers manage their own limits;
`gadget::limit_process_memory` is available on Linux.

Native admission and some mathematical operations still need complete carriers
in memory. Dataset selections bound the working records they read, but do not
establish a memory bound for every subsequent calculation.

Native query batches use temporary disk space for complete result archives.
Set `TMPDIR` to a suitable scratch volume. A RAM backed temporary filesystem
uses host memory even when it is outside a process's RSS measurement.

## Dependency verification

Require the recorded complete SDK header inventory with
`GADGET_ENFORCE_DEPENDENCY_LOCK=ON`. Inspect a built application's inputs with:

```sh
python3 tools/check_dependency.py --build build/presets/cpu \
  --rexgraph /absolute/path/to/rexgraphcpp
```

`--consumed-only` checks every SDK header used by the compiler. A changed
consumed header requires rebuilding and relevant tests.

## Generated test data

Each preset reuses its build directory. Tests write to its `test-data/`
subdirectory. After reviewing results, remove generated fixtures with:

```sh
cmake --build --preset cpu --target gadget_clean_test_data
```

Keep research inputs and results outside build directories. The cleanup target
retains compiled executables.
