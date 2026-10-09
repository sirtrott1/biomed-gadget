# Structural reference calibration

Compare declared structural changes under a specified calibration law. Retain the original structures, observation maps and complete supported calibration fields.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/structural-calibration/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](https://github.com/sirtrott1/biomed-gadget/blob/main/include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
