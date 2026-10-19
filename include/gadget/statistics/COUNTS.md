# Count models

`count_model_analysis` fits Poisson or NB2 models to declared independent
experimental units. The mean uses a log link. NB2 variance is
`mean + dispersion * mean²`.

This is a conventional statistical model over native fields. The exact field
models in `ScientificSession` retain their existing mathematical meaning.

## Choose the model

| Choice | Meaning |
| --- | --- |
| `poisson` | Variance equals the mean |
| `negative_binomial` | NB2 variance with positive dispersion |
| `fixed` dispersion | Use a declared dispersion |
| `profile` dispersion | Fit each feature over a declared interval |
| `pooled` dispersion | Fit one shared dispersion from declared reference fields |
| `log_prior` | Regularize the dispersion profile with a declared normal prior on log dispersion |
| `effect_prior` | Report a normal posterior approximation for the contrast using a declared normal prior |

Shared dispersion uses the sum of the declared reference marginal profile
likelihoods. It does not infer a joint law between reference features.

Dispersion search evaluates a grid, then refines around its best point. Results
retain the grid, convergence and boundary diagnostics. This search does not
certify a global maximum. Priors are supplied by the request; hyperparameters
are not estimated from the data.

## Declare the observations

Each sample names its original grade zero coordinate and an independent
experimental unit. Each feature names an existing quantitative scalar field.
Fields must share their actual source carrier and cutoff. Quantity and units
must match their source declarations.

Counts must be nonnegative integers exactly representable in binary64.
Supported zero remains an observation. Missing or invalid counts make that
feature unavailable; observations are never silently removed. Boundary fits, including ordinary all zero intercept models, retain a separate
status without an effect or reference tail.
Distinct unit labels are an explicit independence declaration, not evidence of
independence. Paired, repeated or correlated observations need a suitable model.

Normalization defaults to none. Supplied positive factors can optionally be
centered to geometric mean one. Median ratio normalization uses only declared
reference fields positive in every selected sample. It computes the ordinary
median of count to geometric mean ratios, then centers the factors. Missing or
invalid reference counts refuse the request. Exposure values, if present,
contribute an additional log offset.

## Fit and interpret

Fisher scoring uses RexGraphCPP matrix solves. Exact rank and design recession
constraints check whether finite coefficients exist before binary64 fitting.
Numerical score convergence alone cannot admit an infinite coefficient limit. Each feature retains convergence, score, likelihood and
Pearson diagnostics. Effect estimates and variances use the declared contrast
on the natural logarithm scale.

Wald tests use a normal reference. Likelihood ratio tests compare a reduced
design selected from the full design columns with a chi square reference.
Per feature dispersion is profiled separately in both models. Shared dispersion
is held fixed in that comparison. Likelihood ratio tests require an unpenalized
dispersion profile. Wald confidence intervals use the unshrunk contrast;
posterior credible intervals are reported separately.

These calibrations are asymptotic and conditional on the declared model.
Uncertainty from normalization and pooled dispersion is not propagated.
Shrinkage uses a normal approximation to the fitted contrast likelihood.
Numerical agreement with a reference implementation does not establish model
adequacy, biological significance or sampling calibration.

Probability correction supports none, Bonferroni, BH and BY. The declared family
includes every requested feature, treating unavailable tests as probability one
for adjustment while leaving their reported probabilities null. BH requires
valid null tails and independence or suitable positive dependence. BY permits
arbitrary dependence with valid null tails. Log probabilities preserve very
small tails; an underflow flag distinguishes a rounded zero from an exact zero.

## Run an analysis

The [request schema](../../../schemas/statistics-count-request.schema.json)
lists every required declaration and execution limit.

```cpp
#include <gadget/statistics/count_analysis.hpp>

auto prepared = gadget::statistics::PreparedCountAnalysis::prepare(study, request);
auto result = prepared.execute();
result.write("counts.rfcxs");
auto display = gadget::statistics::count_result_view(result);
```

```sh
gadget statistics --study study.rfcxs --request counts.json \
  --out-record counts.rfcxs
```

The command also accepts `--workspace workspace.rfcxs --assay ID` and an
optional `--partition`. A pipeline `statistics` step or the
`count_model_analysis` default recipe selects one assay partition. Pipelines
bind its actual study digest. A direct request must supply `source_digest`.
Result steps can select the resulting native feature array for pathway or
evidence operations. Planner answer contracts select `count_model_inference`.
Its answer retains checked native result, request and original study owners;
the conditional statistical assumptions remain part of the answer.

Numerical inputs come from native accessions. Results retain the request,
original study, field identities, support and the sample selection map.
Calculations reuse a prepared design and read one feature at a time. No dense
copy of the complete count matrix is created. `outputs.fitted_means` controls
whether every fitted mean is retained; coefficients and diagnostics remain
available either way.

Limits admit observations, features, design width, numerical work, iterations,
request bytes and estimated working plus output bytes. Shared CPU admission
controls workers and memory reservations. These are execution controls, not a
universal resident memory bound. A request is refused when its declared limits
cannot accommodate its complete scope.
