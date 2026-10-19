# Metric relative change

`metric_relative_novelty` compares fields from two publications through an
explicit carrier map and a metric on the new source. The original and mapped
fields remain available alongside the calculated change, so its coordinates
and scale can be inspected.

The runner checks the calculation independently with rational arithmetic and
can use an optional declared source. Its required executable and producer
arguments are:

```sh
python3 examples/bioinformatics/metric-novelty/replay.py --help
```

A shared label or matching array shape does not provide the comparison map.
Endpoint and metric declarations follow the
[scientific query schema](../../../schemas/scientific-query.schema.json).
