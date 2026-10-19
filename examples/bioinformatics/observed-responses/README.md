# Observed count, ordinal and dose response analysis

Fit declared response models with explicit measurements, exposures and observation support. The result retains complete fields and reports violations of declared readout ranges.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/observed-responses/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
