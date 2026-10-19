# Source intake

This runner exercises strict table and evidence admission with configurable
record, column, category and response counts. It queries the admitted fields,
coverage and source evidence after import. Original bytes remain available
for resolving the numerical values and categorical assertions.

Choose fixture sizes and a new output directory with:

```sh
python3 examples/bioinformatics/source-intake/replay.py --help
```

The generated inputs exercise admission behavior; their sizes are not limits
on supported scientific populations. For ordinary dataset intake, start with
[Python data tools](../../../tools/rcdb/README.md) or the
[assay workspace guide](../README.md).
