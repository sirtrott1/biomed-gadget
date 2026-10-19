# Constrained response models

`constrained_response_model` fits response fields with the physical shape
constraints supplied in the request. Count, ordinal and measured dose inputs
keep their observation maps and units. Feasibility is checked on the declared
source, and the result includes the constraints and compatible coefficient
family.

The runner accepts an authored case or a source with a declared expression
panel. `--dose-response` adds a requested dose comparison. Its command options
are:

```sh
python3 examples/bioinformatics/constrained-responses/replay.py --help
```

The [observed response guide](../observed-responses/README.md) describes the
underlying field models. Request declarations are specified in the
[scientific query schema](../../../schemas/scientific-query.schema.json).
