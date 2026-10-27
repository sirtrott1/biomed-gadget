# Continuous integration

The workflow runs the `cpu` test profile with GCC and Clang against the independent
public RexGraphCPP SDK. It enforces the header inventory and revision recorded
in [RexGraph.lock.cmake](../cmake/RexGraph.lock.cmake). Python dependency wheels
build alongside the compiler suites. Binding checks then use the GCC suite
binaries, with every native provider checked against the same SDK inventory.

Every push checks CI selection, README links, dependency verification and source
export. Eligible README changes can reuse a completed compiler qualification
with identical executable inputs. Reuse requires successful compiler and Python
binding checks. Missing qualification, other input changes and manual dispatch
run both compiler suites and the binding checks.

macOS also compiles all configured contracts, including the HTS and local model
adapters. It checks native CPU and memory observations before compiling the
remaining executables, then exercises count models, resource admission, native
records, file publication and source readers. macOS memory observations use
free plus inactive pages as advisory headroom; this is not a process heap quota.

## SDK source

The default is the public `sirtrott1/rexgraph-cpp` repository at the locked
revision. Repository variables can select another source or compiler runner:

| Variable | Purpose |
| --- | --- |
| `REXGRAPH_SOURCE_REPOSITORY` | GitHub repository containing the SDK |
| `REXGRAPH_SOURCE_REF` | Exact source revision |
| `GADGET_BUILD_JOBS` | Compiler concurrency |
| `GADGET_CI_RUNNER` | Runner label for the compiler jobs |

The complete header lock remains enforced when a source override is used.
Native dependency wheels are cached by source revision and toolchain. Compiler
caches are limited to 1 GB for each compiler and validate source, header and
compiler inputs before reusing an object. C++ consumers are transferred within
the current run so Python checks do not compile them again.

## Archive mode

Set `REXGRAPH_SDK_MODE=archive` and `REXGRAPH_ARCHIVE_SHA256` to use an SDK
archive. It must contain a source project or installed prefix at
`rexgraph-cpp/rexgraphcpp/`.

Choose one archive source:

- `REXGRAPH_ARCHIVE_URL`: HTTPS URL, supplied as a variable or secret.
- `REXGRAPH_ARCHIVE_REPOSITORY` and `REXGRAPH_ARCHIVE_REF`: repository and
  exact commit containing `rexgraph-sdk.zip`. A separate private repository can
  use `REXGRAPH_ARCHIVE_DEPLOY_KEY` for authentication.

Archive variables apply only in archive mode. Dispatch the workflow manually
after changing SDK settings to qualify the selected inputs.

## Resource use

The workflow checks numerical imports and figure rendering before compilation.
Linux compiler concurrency respects processor affinity, CPU quotas and available
host or container memory. The default request is three workers, admitted with
3 GiB per worker and 1 GiB of additional headroom.

Independent tests run two at a time. Numerical libraries use one thread per
process. A new push cancels the previous run for the same branch.

The ordinary `cpu` profile excludes tests labeled `scale`, `recovery` and
`local-inference`. Full source admission, population scale and durability
workloads remain in the `contracts` profile and can require substantially more
memory and scratch space. Select `contracts` when dispatching the workflow on a
suitable runner, or run `ctest --preset contracts` on a suitable host. The
qualification summary names the selected profile.
