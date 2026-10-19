# Linear observations of operator responses

Declare observation maps on an existing source carrier and compare their resulting fields. Matching dimensions do not establish identity between independent sources.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/linear-observations/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
