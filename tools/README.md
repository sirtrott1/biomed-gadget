# Tools

| Task | Tool or guide |
| --- | --- |
| Import, explore and select native datasets | [Python data tools](rcdb/README.md) |
| Profile native loading, workspace admission or raw import | [Load profiling](native-load/README.md) |
| Check a build's SDK inputs | `check_dependency.py` |
| Check a local model build's llama.cpp inputs | `check_llama_dependency.py` |
| Export public source and available history | `export_source.py` |

Dataset preparation scripts declare how particular external sources enter
GADGET. They belong to examples and tools; the production APIs use general
operations.

## Source export

From the repository root, choose a new archive path:

```sh
python3 tools/export_source.py --out build/gadget-source.zip \
  --build build/presets/cpu --rexgraph /absolute/path/to/rexgraphcpp
```

The exporter checks the selected build and SDK, then includes public source,
licensing files, available Git history and `SHA256SUMS`. It leaves the source
repository unchanged. Use `unzip` to restore executable file modes and inspect
history with `git log --all` in the extracted `gadget/` directory.

`--without-history` creates a source snapshot. `--git-dir` selects Git metadata
managed separately from the worktree; that history must be readable. The export
preserves every available object and ref while excluding local ref journals.

Research inputs, generated results, model weights, workspace notes and build
products are excluded. Git checkouts also exclude ignored files. Without Git,
keep generated files outside the public source directories. Existing tracked
fixtures remain included; check contents and distribution terms before adding
new binary fixtures.
