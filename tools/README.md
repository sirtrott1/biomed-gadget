# Tools

These tools support data preparation, dependency inspection and source export.
Dataset adapters describe the interpretation of particular external formats;
the production library uses general data and analysis operations.

| Task | Tool or guide |
| --- | --- |
| Import, explore and select native datasets | [Python data tools](rcdb/README.md) |
| Profile native loading, workspace admission or raw import | [Load profiling](native-load/README.md) |
| Check a build's SDK inputs | `check_dependency.py` |
| Check a local model build's llama.cpp inputs | `check_llama_dependency.py` |
| Export public source and available history | `export_source.py` |

## Source export

From the repository root, choose a new archive path:

```sh
python3 tools/export_source.py --out build/gadget-source.zip \
  --build build/presets/cpu --rexgraph /absolute/path/to/rexgraphcpp
```

The exporter checks the selected build and SDK, then includes public source,
licensing files, Git branches and tags and `SHA256SUMS`. It leaves the source
repository unchanged. Use `unzip` to restore executable file modes and inspect
history with `git log --all` in the extracted `gadget/` directory.

An export with history requires the source files and executable modes to match
the selected commit. Commit public changes first. Untracked source files and
missing tracked files are refused so the restored checkout matches its history.

`--without-history` creates a source snapshot. `--git-dir` selects Git metadata
managed separately from the worktree; that history must be readable. The export
preserves the selected branches, tags and their complete reachable history.
Local saved work, ref journals and unreachable objects are excluded.
Source and selected history must pass the publication check before export.

Research inputs, generated results, model weights, workspace notes and build
products are excluded. Git checkouts also exclude ignored files. Without Git,
keep generated files outside the public source directories. Existing tracked
fixtures remain included; check contents and distribution terms before adding
new binary fixtures.
