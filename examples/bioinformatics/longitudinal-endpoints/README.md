# Longitudinal fields and censored endpoints

`censored_endpoint_bounds` uses recorded lower and upper event times, endpoint
closure and finite support to bound event fractions and restricted durations.
The population, strata and time unit belong to the request. A censored interval
is retained as an interval rather than replaced with a single event time.

The runner also exercises native temporal histories. It accepts an authored
case or an explicit endpoint source:

```sh
python3 examples/bioinformatics/longitudinal-endpoints/replay.py --help
```

These are bounds for the declared observations and population. A fitted
survival distribution or population confidence statement requires additional
assumptions. See [field histories](../versioned-history/README.md) for the
separate temporal correspondence interface.
