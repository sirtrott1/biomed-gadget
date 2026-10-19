# C++ interface

GADGET exposes biomedical analysis and native data interfaces. RexGraphCPP
provides the underlying mathematics, fields and formats.

## Use from CMake

Install both libraries as described in the [build guide](../../cmake/README.md).
Pass both installation prefixes in `CMAKE_PREFIX_PATH`, separated by a semicolon.

For biomedical declarations, queries and workflows:

```cmake
find_package(GADGET CONFIG REQUIRED COMPONENTS core)
target_link_libraries(my_analysis PRIVATE gadget::core)
```

For native records without the JSON and OpenSSL adapters:

```cmake
find_package(GADGET CONFIG REQUIRED COMPONENTS native)
target_link_libraries(my_analysis PRIVATE gadget::native)
```

Source consumers can use `add_subdirectory(path/to/gadget)` and the same targets.

## Interfaces

| Interface | Purpose |
| --- | --- |
| [ScientificStudy](scientific.hpp) | Admit a study and resolve its declared source identities |
| [ScientificSession](scientific.hpp) | Execute queries and plans on that study |
| [Workflows](workflow.hpp) and [pipelines](pipeline.hpp) | Coordinate dependent analyses and retain their results |
| [NativeStudy](native_data.hpp) | Read native publications, access fields and resolve original byte spans |
| [NativeRecord](native_record.hpp) | Retain declarations, results and byte owners in one publication |
| [Native gluing](native_gluing.hpp) | Combine declared local sources and fields through explicit interfaces |
| [Data fields](data_field.hpp) | Bind imported fields to original axes, units and record addresses |

`ScientificSession::execute_record` and `execute_plan_record` return native
query and plan records. The workflow and pipeline functions return their own
operation records. Results retain source publications and complete output fields;
session preparations are reused only for their study. See the
[scientific example](../../examples/scientific/README.md).

## Read selected coefficients

`NativeRecord::encoded_tensor(key)` reads a field at its recorded publication.
The reader keeps the original bytes and endpoint alive:

```cpp
const auto record = gadget::NativeRecord::read("analysis.rfcxs");
if (const auto field = record.encoded_tensor(field_key)) {
    const std::vector<std::size_t> coordinates{0, 2, 7};
    const auto selected = field->tensor().select(coordinates);
}
```

`field_key` and coordinates come from the actual publication. A sorted selection
keeps the original tensor domain and support; unrequested coordinates are absent.
A measured zero remains different from an absent coordinate.

Archive admission validates the complete history. A missing field key throws;
a field without a matching base or checkpoint encoding returns no reader.
Selecting coefficients does not imply that archive admission reads only those
coefficients.

## Ownership

Retain a selected `NativeRecord` while using its `value()` or `tensor()`
references. Likewise, keep a `RetainedRecord::native()` selection alive while
using a reference or byte span returned from it.

`NativeRecord::archive_source(key, path)` writes a new archive, maps it read only
and removes the temporary path. The returned owner keeps the complete bytes.
Use a new path on a suitable filesystem.

Original identities, carrier maps, field support and metrics belong to the
publication. Equal dimensions do not identify independent sources.
The [native data guide](../../examples/native-data/README.md) covers byte spans,
declarations, packed results and history selection.
