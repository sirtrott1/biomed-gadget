# Independent dependencies

Install GADGET with the two independent RexGraph repositories. These commands
create the directory layout used below:

```sh
git clone https://github.com/sirtrott1/rexgraph-python.git rexgraph-python
git clone https://github.com/sirtrott1/rexgraph-cpp.git rexgraphcpp
git clone https://github.com/sirtrott1/biomed-gadget.git biomed-gadget
cd biomed-gadget
```

The C++ repository is named `rexgraph-cpp` on GitHub. The local directory
`rexgraphcpp` matches the SDK's standalone examples and GADGET's default source
path. Explicit source paths and installed prefixes can use other names.

Repository directory names and Python imports serve different purposes:

| Source checkout | Installed interface | Role |
| --- | --- | --- |
| `rexgraph-python` | `rexgraph` | Python mathematics and compiled Cython kernels |
| `rexgraph-python/rcdb` | `rcdb` | Versioned native records and storage |
| `rexgraph-python/rcql` | `rcql` | Native queries and executable programs |
| `rexgraph-python/bindings/cpp` | `rexgraph_cpp` | Optional bridge to the independent C++ SDK |
| `rexgraphcpp` | `RFC`, `rfc::` CMake targets | C++ mathematics, formats and execution |
| `biomed-gadget/tools/rcdb` | `gadget_data`, `gadget_rcdb` | Dataset intake, selection, provenance and application integration |
| `biomed-gadget` | `GADGET`, `gadget::` CMake targets | Biomedical declarations and compiled application |

The bridge belongs to the Python repository. Its Python and Cython sources are
in the top level `rexgraph_cpp/` directory; its adapters and package build are in
`bindings/cpp/`. It compiles against an external SDK installation. GADGET uses
the SDK directly in C++ and the bridge for Python operations that require it.
Core, RCDB and RCQL also work without the bridge. Agent and System are optional
packages in the Python repository. The repositories retain separate histories,
builds and licenses. GADGET does not bundle their source or compiled libraries.
Core, the C++ SDK and the optional Python bridge publish Apache 2.0 licenses.
The bridge retains its own
[license notice](https://github.com/sirtrott1/rexgraph-python/blob/732f822aa367ad161d5131c747b91f69e05ae7a8/bindings/cpp/LICENSE).

## Install from source

The commands below start inside `biomed-gadget` on Linux. The source revisions
are pinned before compilation so the C++ application and Python bridge use the
same SDK. Existing checkouts can instead use these revisions in separate clones.

```sh
git -C ../rexgraph-python checkout --detach 732f822aa367ad161d5131c747b91f69e05ae7a8
git -C ../rexgraphcpp checkout --detach cc367801df8b2316be9d7930fa4cbcaabed655ad
```

On Debian or Ubuntu, install the system requirements:

```sh
sudo apt-get install build-essential cmake ninja-build pkg-config git \
  libgmp-dev libjson-c-dev libssl-dev libopenblas-dev python3-dev python3-venv
```

GADGET requires Python 3.11 or newer for the data adapter. Core compilation
requires Python development headers and BLAS/LAPACK; building the bridge also
requires the C++ toolchain and GMP. Pip obtains the declared build tools in
isolated build environments.

First install the independent C++ SDK:

```sh
export GADGET_RFC_PREFIX="$PWD/../rexgraph-install"
cmake -S ../rexgraphcpp -B ../rexgraphcpp/build/sdk -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DRFC_ENABLE_VULKAN=OFF \
  -DRFC_BUILD_APPS=OFF -DRFC_BUILD_LOCAL_TOOLS=OFF \
  -DRFC_BUILD_TESTS=OFF -DRFC_BUILD_BENCHMARKS=OFF \
  -DRFC_BUILD_CAPABILITY_DEMO=OFF -DRFC_BUILD_LEARNING=OFF
cmake --build ../rexgraphcpp/build/sdk --parallel 1
cmake --install ../rexgraphcpp/build/sdk --prefix "$GADGET_RFC_PREFIX"
```

Install Core, storage, queries and the bridge in one Python environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install ../rexgraph-python \
  --config-settings=compile-args=-j1
.venv/bin/python -m pip install '../rexgraph-python/rcdb[sql]' ../rexgraph-python/rcql
.venv/bin/python -m pip install ../rexgraph-python/bindings/cpp \
  --config-settings="cmake.define.CMAKE_PREFIX_PATH=$GADGET_RFC_PREFIX"
.venv/bin/python -m pip install './tools/rcdb[native]'
.venv/bin/python -m pip check
```

The repository names do not change the distribution names: Core installs as
`rexgraph` and the bridge as `rexgraph-cpp-bindings`. The Python repository also
provides `make install-cpp CPP_SDK_PREFIX=...` for its own selected environment.
Install Agent or System separately when an application uses those services.

Build GADGET against the same SDK prefix:

```sh
cmake --preset installed -DCMAKE_PREFIX_PATH="$GADGET_RFC_PREFIX" \
  -DGADGET_ENFORCE_DEPENDENCY_LOCK=ON
cmake --build --preset installed --parallel 1
cmake --install build/presets/installed --prefix "$PWD/../gadget-install"
.venv/bin/python -I tools/rcdb/check_native_stack.py \
  --build-info build/presets/installed/build-info.json
```

The final check reads every compiled bridge provider's SDK receipt. Equal
package versions alone do not establish equal SDK inputs. Add `--services` to
check Agent and System too. C++ consumers need both installation prefixes in
`CMAKE_PREFIX_PATH`. See the [C++ installation instructions](include/gadget/README.md#cmake).

## Check the installation

Add the installed executables and the selected Python environment to the shell
path. This also lets dataset export find the compiled native field importer.

```sh
export PATH="$PWD/.venv/bin:$PWD/../gadget-install/bin:$PATH"
gadget --help
gadget-data --help
python -I -c "import rexgraph, rexgraph_cpp, rcdb, rcql, gadget_data"
```

`export-native` requires an explicit `--importer`. With this shell path, pass
`--importer gadget_native_field_import`, or supply an absolute executable path.

The [first analysis](README.md#first-analysis) constructs a synthetic study and
executes a native query. The [data guide](tools/rcdb/DATA.md) covers bounded
intake, source lookup and selected fields. HDF5 intake requires `h5py`; install
it in the same Python environment when needed. Local models, remote database
transport and GPU adapters have separately selected dependencies.

## Dependency revisions

The source pins are Python
`5f9f4f9a7934c5c14e3f7344b15c0d8ea4147605`, version 1.3.2, and SDK
`4d52e37296a31bc6d7e6b43cdf091261e41518d6`, version 0.2.1.
The complete 470 header inventory has SHA256
`1e6adda1dec387c54bac9acc513a62806c35991cf593a163f416bb21bd2186a6`.
The SDK pin is in [cmake/RexGraph.lock.cmake](cmake/RexGraph.lock.cmake);
`GADGET_ENFORCE_DEPENDENCY_LOCK=ON` requires that inventory.

Retained model records keep their original engine identities. Preserve their
original runtimes. Installing a newer bridge does not certify an old model for
the new engine; explicit replay and complete field comparisons establish a new
qualified state. Dataset source identities, original coordinates, supported zero,
missing support and declared numerical enclosures remain part of the scientific
record across software updates.
