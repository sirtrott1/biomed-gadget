# Complete source admission

Import a declared table with its field meanings, coordinate identities and observation support. Original bytes remain provenance for the working native records.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/source-intake/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
