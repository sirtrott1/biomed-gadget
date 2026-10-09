# Physically constrained response models

Fit response fields under declared affine constraints. Preserve the constraint operators, source metrics and complete compatible coefficient family.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/constrained-responses/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](https://github.com/sirtrott1/biomed-gadget/blob/main/include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
