# Predeclared Gaussian sequential stopping

This adds a reusable, dependency-free likelihood-ratio e-process and a frozen
norm-5 calibration example.  It is **methodology only**: it must not be applied
to a production run that started before the rule was committed.

For an aligned batch vector `X` and two frozen predictions with a shared frozen
covariance `Sigma`, the increment is

```text
log p_alt(X)/p_null(X)
  = (mu_alt-mu_null)' Sigma^-1 [X-(mu_null+mu_alt)/2].
```

Under the named Gaussian null, the exponential of the cumulative log ratio is
a nonnegative martingale.  Therefore Ville's inequality gives
`P_null(sup_t E_t >= 1/alpha) <= alpha`, including optional stopping.  Running
the reverse ratio gives the corresponding bound when the alternative is the
data-generating endpoint.  At the maximum batch count, failure to cross either
boundary is reported as `inconclusive`; it is not forced into a decision.

## Frozen example

The example uses 5M-replica batches, a minimum of 10 batches, a maximum of 100
(500M replicas), and `alpha=0.01` in each direction.  It includes H4, H12, and
zero as generating models and evaluates the three preordered pairwise tests.
The covariance is the independent child-only pilot scale.  Parent-amplitude
uncertainty is intentionally excluded, so this example cannot authorize a
production decision; a future preregistration must replace it with the complete
frozen covariance for its chosen primary statistic.

Reproduce the deterministic calibration with:

```bash
python scripts/sequential_gaussian_eprocess.py \
  predictions/norm5_sequential_stopping_20260829.yaml \
  --simulations 5000 --seed 126 \
  --output results/norm5-sequential-stopping/calibration.json
```

The JSON reports, for every ordered test under every declared generating
model, decision frequencies, expected batches, and compute fraction relative
to the fixed 500M design.  The mathematical error guarantee applies only at a
test's two declared endpoints.  Simulation under the third model is a behavior
check, not a new type-I guarantee.
