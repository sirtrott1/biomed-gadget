# Finite reference calibration

`finite_reference_calibration` evaluates a model statistic under a supplied
finite reference law. The states, probabilities and observed statistic are
part of the request and result. Inclusive tails retain ties explicitly.

The example uses paired reference laws for a declared linear model. It can
construct an authored case or read an expression source and feature selection.
Run options are:

```sh
python3 examples/bioinformatics/finite-calibration/replay.py --help
```

Choosing the reference law requires a scientific justification for the
experiment. Exact arithmetic evaluates that law; it does not establish the
exchangeability or sampling assumptions. See
[structural calibration](../structural-calibration/README.md) for blocked
permutation declarations.
