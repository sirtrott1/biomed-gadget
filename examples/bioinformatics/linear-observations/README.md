# Linear response observations

`linear_response_observation` applies an observation map to an operator
response on a declared source carrier. The source operator, accession map,
observation map and metrics are supplied explicitly. The calculation can act
on the requested fields without constructing a full inverse.

The runner creates source publications and compares the reported observations
with independent expected values. Its producer and input options are:

```sh
python3 examples/bioinformatics/linear-observations/replay.py --help
```

Equal array dimensions do not identify two independently published carriers.
Use the source identities and maps in the
[query declaration](../../../schemas/scientific-query.schema.json).
