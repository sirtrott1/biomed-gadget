# Response transitions

`declared_response_transition` compares response fields at two source endpoints.
Each endpoint names its own publication, history cutoff and field. Explicit
maps relate their carriers, so changes in coordinates and changes in values
can be interpreted within the same calculation.

The runner constructs source publications, executes the transition and checks
its terms independently. It requires a GADGET executable and a native source
producer:

```sh
python3 examples/bioinformatics/declared-responses/replay.py --help
```

For the endpoint declarations, see the
[scientific query schema](../../../schemas/scientific-query.schema.json).
[Field histories](../versioned-history/README.md) extend these comparisons to
an ordered sequence of publications.
