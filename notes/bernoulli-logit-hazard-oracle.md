# Exact Bernoulli logit-hazard oracle

This Issue #100 method slice puts three finite-event derivative calculations in
one exact regression.  For a monotone event `A`, occupation count `K`, event
probability `F`, and natural parameter `eta=logit(p)`, direct differentiation
of the Bernoulli product measure gives

```text
dF/deta = Cov(1_A, K) = p(1-p) dF/dp.
```

Dividing by `F(1-F)` gives

```text
d logit(F)/deta = E[K | A] - E[K | not A].
```

The Fraction oracle enumerates the 16 configurations of the declared event
`K>=2` on four bits at `p=1/3,1/2,2/3`.  At every point it also recomputes
`dF/dp` as the sum of pivotal probabilities, so the Russo, score-covariance,
and conditional-occupation routes must agree without numerical tolerance.

Run:

```text
python3 scripts/bernoulli_logit_hazard_oracle.py
python3 -m unittest tests/test_bernoulli_logit_hazard_oracle.py -v
```

## Boundary

This is a finite threshold event, not a percolation scaling calculation.  It
introduces no Monte Carlo evidence, does not estimate a pivotal/four-arm
exponent or universal amplitude, and does not establish that a chosen
production observable has favorable logit conditioning.  Issue #100 remains
open.
