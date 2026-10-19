# Python data tools

Use `gadget_data` to import scientific data, explore its metadata and retrieve
selected measurements. Original files remain provenance. Working records retain
values, coordinates, observation masks and the source versions they came from.

Core provides mathematical objects, RCDB stores versioned records and RCQL
executes queries. GADGET declares the data meanings and connects these operations
to the compiled application.

## Install

Follow the [dependency guide](../../DEPENDENCIES.md) to install the matching
RexGraph packages. From the GADGET root, use that Python environment to install:

```sh
python -m pip install ./tools/rcdb
```

Core, RCDB and RCQL 1.3.2 or newer are required. SQLite is the default working
store. HDF5 intake additionally needs `h5py`. Native C++ export needs the
`rexgraph_cpp` bridge and GADGET's compiled field importer. Agent and System
are optional services.

## Import and explore

Given an original `counts.npy` matrix, create a new dataset and inspect it:

```sh
gadget-data import-dense counts.npy counts.native \
  --reader npy --name counts --quantity reported_count --units count
gadget-data explore counts.native
gadget-data provenance counts.native \
  --field measurements --rows 10:20 --columns 4,1
```

Catalog and provenance queries read metadata before numerical fields. Select
only the rows and columns needed for a calculation:

```python
from gadget_data import DataStore

with DataStore("counts.native") as data:
    selected = data.select("measurements", rows=(19, 10), columns=(4, 1))
    values = selected.values
    observed = selected.observed
```

Repeated labels remain separate coordinates. `observed` distinguishes a measured
zero from a missing observation. Source mappings remain attached to the selection.

## Continue with an analysis

The [data API guide](DATA.md) covers declared tables, NPY, HDF5, Matrix Market,
shared context, lookup, RCQL programs and native field export. It includes
examples for retaining and restoring calculations with their original owners.

The [C++ native data guide](../../examples/native-data/README.md) covers the
compiled representation. Python export passes typed fields through the bridge;
JSON is an input or display interface.

## Store complete GADGET publications

`gadget_rcdb` stores a complete RFCXS archive as one RCDB record for version
selection and replay. Use it when the archive fits the declared record budget;
large datasets can instead use the bounded records in `gadget_data`.

The [publication API](PUBLICATIONS.md) covers codec setup, record retrieval,
execution, replay, copying and compression.

## Check the native stack

When using the bridge, check that every compiled provider agrees with the
GADGET build's SDK inventory:

```sh
python -I tools/rcdb/check_native_stack.py \
  --build-info build/presets/cpu/build-info.json
```

Add `--services` when Agent and System are installed. Package version numbers
alone do not establish equal SDK inputs. Existing models retain their recorded
engine requirements; [dependency compatibility](../../DEPENDENCIES.md#reproduction-inputs)
explains replay and runtime preservation.
