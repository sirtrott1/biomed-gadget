# Hodge fields

`hodge_frame_decomposition` decomposes a field on a selected grade of a
relational complex using the supplied boundary operators and metrics. The
output includes the complete sector fields and reconstruction checks in the
original carrier coordinates.

The runner publishes authored sources and evaluates the decomposition with
unit or positive metrics. It uses a native source producer alongside the
GADGET executable:

```sh
python3 examples/bioinformatics/hodge-frames/replay.py --help
```

The query's publication and source sequence select the actual history being
analyzed. [Native study declarations](../../native-studies/README.md) describe
how field quantities, units and metrics are bound to that history.
