# Factored temporal maps

A temporal query can return its complete map as enumerated coefficients or as
a composition of native factors. The factored form preserves the map's domain,
codomain and action while avoiding expansion of the full coefficient array.

This runner executes both `map_output` choices for an existing case containing
`study.json` and a history query. Its verifier expands the factors and compares
them with the enumerated result:

```sh
python3 examples/bioinformatics/factored-temporal/replay.py --help
```

The default query filename is `history-query.json`; `--query` selects another
file in the case directory. See [field histories](../versioned-history/README.md)
for publication endpoints and correspondence declarations.
