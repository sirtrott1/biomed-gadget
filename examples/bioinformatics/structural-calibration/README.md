# Structural calibration

`structural_reference_calibration` evaluates structural statistics under a
specified exchange law. The example declares blocked permutations, including
which field axes move together and which coordinates can exchange within a
block. The null hypothesis and population remain part of that declaration.

The runner accepts authored measurements or a source and panel:

```sh
python3 examples/bioinformatics/structural-calibration/replay.py --help
```

The law's exchangeability assumptions must fit the experimental design.
[Finite calibration](../finite-calibration/README.md) describes the underlying
finite reference calculation, and [ablations](../structural-ablations/README.md)
describe the structural comparisons.
