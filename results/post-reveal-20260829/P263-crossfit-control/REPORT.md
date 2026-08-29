# Issue #263 cross-fitted control audit

## Result

The revealed Phase-E sufficient statistics do **not** support a strict
variance-reduced score.  Cross-fitting the conditional mean of `J/2` on the
four event categories (`1234`, `12|34`, `14|23`, other) is an exact algebraic
decomposition of

```text
d_Q log P(14|23) = E[J/2 | 14|23] - E[J/2].
```

For any category function `m`, including one trained on the opposite batch
fold,

```text
E[J/2|H]-E[J/2]
= E[J/2-m(X)|H]-E[J/2-m(X)] + m(H)-E[m(X)].
```

The held-out residual and explained pieces therefore sum to the original
estimator in every fold.  On the revealed streams the maximum numerical
reconstruction errors are `4.55e-13` at level 1 and `1.82e-12` at level 2.

| revealed acquisition | conditional-score variance / primary variance | resulting frozen tangent-shape chi2 / df |
|---|---:|---:|
| level 1, 200k/geometry | 1.000000 | 1.097745746 / 3 |
| level 2, 500k/geometry | 1.000000 | 5.872171777 / 3 |

No new chi-square is claimed from these data.

## Why batch event counts are not strict controls

Every stored candidate auxiliary has an unknown finite-lattice expectation:

- `J` and `J^2` moments;
- the three event counts/probabilities;
- their event-weighted `J` and `J^2` moments.

A target-preserving control variate needs a known mean.  If a batch event
count is centered by the opposite fold, a common beta gives corrections
`beta*(X_A-X_B)` and `beta*(X_B-X_A)`, which cancel exactly.  Allowing different
fold betas leaves a beta-difference times mean-difference product; that is not
an exact zero-mean correction and can move the target after reveal.

Equivalently, the required training cross-matrix `Sigma_YW` cannot be formed
because the current raw files contain no exact-zero `W`.

## Minimal executable next-run control

The runner now retains one additional integer per synchronized
geometry/batch:

```text
sum_b = sum over samples of the number of open bonds.
```

For a geometry with `E_g` edges and `n_g` samples,

```text
W_g = (2*sum_b - n_g*E_g) / sqrt(n_g*E_g),
E[W_g] = 0 exactly at p=1/2.
```

Let `Y_b` be the three active, amplitude-anchored batch contributions to the
frozen tangent residual.  The declared even/odd crossfit is

```text
beta_even = Sigma_YW,even * Sigma_WW,even^+
beta_odd  = Sigma_YW,odd  * Sigma_WW,odd^+

Y_b^cv = Y_b - beta_odd * W_b   for even held-out batches,
Y_b^cv = Y_b - beta_even * W_b  for odd held-out batches.
```

Each beta is independent of its held-out bond control, whose expectation is
known to be zero.  Thus the correction preserves the target conditional on
the training fold.  The covariance is the sample covariance of the 100
cross-fitted contributions divided by 100, followed by the unchanged
three-dimensional GLS tangent score.

Only `sum_b` is required for this synchronized **batch-level** crossfit.  A
configuration-level analytic multi-geometry beta would additionally need
`sum J_g*b_h` and `sum I_p,g*J_g*b_h` for every geometry pair `(g,h)`.  The
control covariance itself is exact from the edge counts and declared
shared-edge RNG enumeration; `b_g*b_h` is only an optional audit moment.

## Implementation and validation

- `src/p263_boundary_qscore_pilot.cpp` emits `sum_b` without changing the
  lattice, connectivity classifier, `J`, or RNG function.
- `scripts/analyze_p263_crossfit_control.py` audits old streams, refuses to
  manufacture a strict score when `sum_b` is absent, and automatically runs
  the even/odd bond control on a future compatible stream.
- `experiments/p263_boundary_qscore_control_phaseF_20260829.yaml` freezes the
  schema and formula but deliberately leaves the next seed/counter domain
  undeployed.
- A synthetic known-control regression reduces covariance trace to
  `0.1174004` of the raw value while preserving the residual vector to floating
  precision.  This validates the executable algebra only; it is not a forecast
  of the physical P263 variance reduction.
- All 25 stacked `test_p263_*.py` tests pass.

## Scientific interpretation

The current variance problem cannot be solved by relabeling the already
stored event/J moments.  The information missing from Phase E is not another
ordinary connectivity probability; it is an auxiliary fluctuation with an
exact finite-lattice expectation.  Open-bond count is the smallest such
quantity already computed inside the runner.

This converts the observed "same rare-event count at doubled span" limitation
into a falsifiable next score: a small fresh pilot can directly report the
cross-fitted covariance-trace ratio before any larger acquisition is allowed.
If total bond count explains little of the `J=2k+b` noise, the next mechanism
candidate is a declared local cycle-rank proxy with a known Bernoulli mean,
not post-reveal regression on event counts.

## Claim boundary

- Exact revealed-data conclusion: conditional event-score regression is a
  variance-one algebraic no-op; primary scores are unchanged.
- Exact next-run mechanism: centered bond count has zero expectation at
  `p=1/2`, so opposite-fold beta application preserves the target.
- Exploratory only: the achievable physical variance reduction is unknown
  until a new seed/counter domain is frozen and acquired.
