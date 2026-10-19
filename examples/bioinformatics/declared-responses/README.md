# Declared response transitions

Evaluate response fields from declared designs and observation maps. Source identities and field support remain part of the native calculation.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/declared-responses/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
