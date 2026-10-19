# Statistics and approximation

`gadget::statistics` provides conventional statistical methods and numerical
analysis on GADGET's native fields. Count models use a declared design,
normalization and dispersion policy. Pathway tests use a declared feature
universe, mapping and reference law. Results carry the original coordinates,
support and source records alongside estimates, tests and numerical diagnostics.

The module is separate from the primary exact field and model interfaces. A
Poisson likelihood or normal reference tail introduces assumptions about the
observations; exact evaluation of a finite law establishes its arithmetic,
while its suitability for an experiment remains a scientific choice.

## Interfaces

| Header | Purpose |
| --- | --- |
| [count_analysis.hpp](count_analysis.hpp) | Native Poisson and NB2 models, normalization, dispersion, Wald or likelihood ratio tests and effect shrinkage |
| [counts.hpp](counts.hpp) | Prepared designs, count likelihoods, reference tails and exact or logarithmic multiple testing correction |
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

The [count model guide](COUNTS.md) covers request choices, assumptions and
resource controls. The [implementation guide](../../../src/statistics/README.md)
describes how these interfaces call the numerical and native providers.
