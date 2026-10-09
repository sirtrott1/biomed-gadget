# Longitudinal fields and censored endpoints

Analyze declared longitudinal fields and endpoint support. Event status, censoring and time coordinates remain distinct inputs to the calculation.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/longitudinal-endpoints/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](https://github.com/sirtrott1/biomed-gadget/blob/main/include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
