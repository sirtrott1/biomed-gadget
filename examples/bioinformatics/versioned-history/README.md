# Versioned field histories

`temporal_field_history` connects publication endpoints through declared
correspondences and lineages. Each endpoint names an actual field on its source
complex. The query keeps those identities when applying an ordered sequence
of maps, even when the source coordinates change.

The runner builds an authored history or uses an explicitly supplied source
or fixture:

```sh
python3 examples/bioinformatics/versioned-history/replay.py --help
```

For a complete map in a more compact representation, see
[factored temporal maps](../factored-temporal/README.md). The factors preserve
the same map and original endpoint coordinates.
