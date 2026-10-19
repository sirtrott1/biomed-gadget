# Historical field execution and coupled event observations

Execute ordered scientific queries against declared source versions. The example retains complete fields and the dependencies between successive operations.

The example runner describes its required inputs and output directory:

```sh
python3 examples/bioinformatics/temporal-execution/replay.py --help
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
