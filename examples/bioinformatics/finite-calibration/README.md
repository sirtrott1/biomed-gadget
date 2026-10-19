# Models under a declared finite reference law

Calibrate declared statistics against an explicit finite reference law. Preserve each reference state, its probability and the source ownership of the measured statistic.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/finite-calibration/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
