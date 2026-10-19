# Continuous integration

The workflow runs CPU contracts with GCC and Clang against the independent
public RexGraphCPP SDK. It enforces the header inventory and revision recorded
in [RexGraph.lock.cmake](../cmake/RexGraph.lock.cmake).

Every push checks CI selection, README links, dependency verification and source
export. Eligible README changes can reuse a completed compiler qualification
with identical executable inputs. Missing qualification, other input changes
and manual dispatch run both compiler suites.

## SDK source

The default is the public `sirtrott1/rexgraph-cpp` repository at the locked
revision. These repository variables select another source:

| Variable | Purpose |
| --- | --- |
| `REXGRAPH_SOURCE_REPOSITORY` | GitHub repository containing the SDK |
| `REXGRAPH_SOURCE_REF` | Exact source revision |
| `GADGET_BUILD_JOBS` | Compiler concurrency |

The complete header lock remains enforced when a source override is used.
No compiled SDK objects are uploaded or cached by this workflow.

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
Independent tests run two at a time; source admission, trajectory scale and
durability fault tests run alone. Numerical libraries use one thread per process.
A new push cancels the previous run for the same branch.
