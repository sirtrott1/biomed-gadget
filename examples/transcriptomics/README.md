# Matched expression fields

`paired_field_contrasts` compares supported measurement fields on explicitly
matched units. Each pair declares its left and right fields, positive rational
scales and the evidence supporting those scales.

The result retains both scaled measurements, their difference, its shared
component and its orthogonal residual. A zero denominator leaves a mean ratio
undefined. Measured zero and missing support remain distinct.

Declare sample correspondence and normalization before calculation. The observed
contrast is not a fitted population effect or a sampling probability.

See the [assay query schema](../../schemas/assays-query.schema.json),
[pipeline guide](../pipelines/README.md) and [C++ interface](../../include/gadget/README.md).
