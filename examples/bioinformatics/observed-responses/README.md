# Observed response models

`observed_response_model` relates measurement fields to a declared design and
observation map. Count and exposure measurements, ordered categories and
measured dose responses use their respective readout meanings. The result
includes fitted fields, residuals and checks against declared readout ranges.

The example can construct an authored case or use an expression source with
an explicit panel. Its input and worker options are:

```sh
python3 examples/bioinformatics/observed-responses/replay.py --help
```

This is the native field response interface. Poisson and negative binomial
likelihood inference, Wald tests and dispersion estimation are available in
the separate [count statistics module](../../../include/gadget/statistics/COUNTS.md).
[Constrained responses](../constrained-responses/README.md) add physical shape
requirements to the field calculation.
