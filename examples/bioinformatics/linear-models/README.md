# Linear models

`linear_model_comparison` represents measurement fields using declared
predictor terms, including products and category indicators. Multiple response
fields can share the same design. A deficient design retains its coefficient
family, while estimable contrasts can still have identified values.

This example constructs continuous and categorical comparisons, checks the
result through a pipeline and verifies that relocation preserves the source
owners. An optional NumPy benchmark compares numerical fits. Run it with:

```sh
python3 examples/bioinformatics/linear-models/replay.py --help
```

These field comparisons use the [scientific interface](../../../include/gadget/scientific.hpp).
Conventional count likelihoods and statistical tests have a separate
[statistics interface](../../../include/gadget/statistics/README.md).
