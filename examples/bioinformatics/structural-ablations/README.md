# Structural ablations

`source_structure_ablation` compares an original source with specified
omissions of relation roles and field spans. It retains the modified
structures, their maps and the resulting field comparisons. This makes each
omission an inspectable change to the declared source.

The runner can construct an authored case or read a source with a selected
panel. Its options are:

```sh
python3 examples/bioinformatics/structural-ablations/replay.py --help
```

An ablation describes the effect of a specified omission within the model;
it does not by itself identify a biological cause. A declared reference law
can be added through [structural calibration](../structural-calibration/README.md).
