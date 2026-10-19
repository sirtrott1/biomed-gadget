# Testing

Run commands from the repository root. Tests need Python 3.11 or newer and
Git in addition to the [C++ build requirements](../cmake/README.md).

## CPU suites

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r tests/requirements.txt
cmake --preset contracts -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build --preset contracts
ctest --preset cpu
```

| Command | Scope |
| --- | --- |
| `ctest --preset cpu` | Ordinary CPU tests; excludes population scale, recovery and local inference |
| `ctest --preset contracts` | Complete suite registered in the CPU build |
| `ctest --test-dir build/presets/cpu -R NAME --output-on-failure` | Selected contracts |

The presets use one test process and one numerical library thread at a time.
Complete suites include workloads whose resource use depends on the input
population and model weights. Select them on a machine with sufficient memory
and scratch space.

## Optional profiles

For HDF5 and HTSlib, install the format dependencies and test requirements:

```sh
.venv/bin/python -m pip install -r tests/requirements-formats.txt
cmake --preset formats -DBUILD_TESTING=ON -DGADGET_BUILD_EXAMPLES=ON \
  -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build --preset formats
ctest --preset formats
```

Python integration contracts require the matching installed
[Core, RCDB, RCQL and optional bridge](../DEPENDENCIES.md).
[Dataset contracts](../tools/rcdb/DATA.md#contracts) exercise original values,
coordinate maps, missing observations and stored queries.

Local inference tests use `GADGET_TEST_MODEL=/path/to/model.gguf`.
`GADGET_TEST_NOVELTY_MODEL` can select a separate planner; it defaults to the
same model. Each model needs qualification for the roles it performs.

Physical GPU and multimodal tests require the selected device, libraries and
model. A software Vulkan run does not qualify a physical GPU.

## Contributing and reporting problems

Test the affected operations and compare results against independent expected
values. Describe input meanings, support, maps, arithmetic and failure conditions.
Use the existing RexGraph mathematical and storage providers.

A bug report should include the command, input scope, compiler, build profile,
`build-info.json`, dependency check and test output. Keep private datasets and
model weights out of reports.

To check compatibility with a retained scientific publication, configure
`GADGET_REPLAY_FIXTURE` with its file path and run
`gadget_native_boundary_replay_contracts`. The input may be RFCXS or gzip
compressed RFCXS. The contract verifies native replay and unchanged input bytes.
Keep publications containing private data outside the source tree.

Generated fixtures can be removed with the
[build cleanup target](../cmake/README.md#generated-test-data). CI configuration
is documented in [.github/CI.md](../.github/CI.md).
