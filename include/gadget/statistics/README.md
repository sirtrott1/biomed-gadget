# Statistics and approximation

This module contains conventional statistical methods and numerical analysis
recipes. Its C++ interface uses `gadget::statistics`; implementations live in
[`src/statistics`](../../../src/statistics).

Methods consume native sources and retain their coordinates, support and
provenance. Each method declares its model assumptions. Approximate calculations
retain their numerical diagnostics. Exact arithmetic can evaluate a statistical
law exactly without establishing that the law describes the observations.

## Available interfaces

| Header | Purpose |
| --- | --- |
| [pathways.hpp](pathways.hpp) | Pathway overrepresentation under a declared uniform subset law, with exact tails and selected multiple testing correction |
| [execution.hpp](execution.hpp) | Contextual field analysis, prepared queries and numerical Green solves with retained residual diagnostics |

Link `gadget::core`. For a pathway request and its owning native records:

```cpp
#include <gadget/statistics/pathways.hpp>

auto result = gadget::statistics::pathway_overrepresentation(
    analysis, gene_sets, request);
```

The [measurement analysis guide](../../../examples/source-formats/ANALYSIS.md)
describes request declarations and native pipelines. Existing
`gadget/pathways.hpp` and `gadget/execution.hpp` includes remain available for
source compatibility. Rebuild C++ consumers after updating the library.

## Native mathematics

Exact field regression, coefficient families, constrained responses, Hodge
operations, metric moments and finite reference calculations retain their
[scientific interfaces](../scientific.hpp). Statistical methods call those
providers when applicable. Numerical kernels remain supplied by RexGraphCPP.

Poisson and negative binomial likelihood models, Wald inference, dispersion
pooling and effect shrinkage are not implemented in this module yet.
