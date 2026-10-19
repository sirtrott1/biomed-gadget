# Publication indexed relational fields

Read and replay versioned scientific records at explicit history cutoffs. Earlier results remain tied to their original source versions.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/versioned-history/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
