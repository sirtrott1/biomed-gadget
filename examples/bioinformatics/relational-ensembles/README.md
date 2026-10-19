# Field ensembles

`relational_field_ensemble` analyzes a declared collection of measurement
fields with positive weights. Weighted means, covariance contractions and
effective dimensions use the declared metric. The source relations, field
support and grouping choices remain part of the calculation. The weights
describe the selected collection; a population probability model needs a
separate declaration.

The runner constructs multiway relations and an optional higher filling, or
reads a source with an explicit panel. It exposes the field weights and native
worker budget through its declarations:

```sh
python3 examples/bioinformatics/relational-ensembles/replay.py --help
```

`--unfilled` selects the authored case without its grade two filling. This
changes the declared complex rather than discarding a result component. See
[Hodge fields](../hodge-frames/README.md) for individual field decompositions.
