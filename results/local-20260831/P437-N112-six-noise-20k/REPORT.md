# N112 six-noise high-pass: nonzero topology, unresolved noisy readout

The C3 topology observer has **provably nonzero degree>=5 content**, but the
literal six-level common-noise Monte Carlo estimator is impractical at this
pilot scale. Stop at the frozen 20,000 replicas; do not solve this by blind
sample expansion.

## Estimand and exact new certificate

The measure is exactly the preceding `2402a33` N112 rho-child production:
224 independent square **bonds**, p=1/2, and three ordered period matrices
`[[8,8],[0,14]]`, `[[16,4],[0,7]]`, `[[16,12],[0,7]]`. This is not site
percolation at p~0.5927. On the same field, `E_i=1_rank0+1_rank2`, and
`F=(E_0+zeta^-1 E_1+zeta^-2 E_2)/3`.

Under stationary product resampling, `E[HF]=0` for every F. Therefore the
primary is the self-source spectral energy

`A = Re E[conj(F(X)) HF(X)] = sum_{j>=5} h_j ||F_j||^2`.

The imaginary part is zero. Subtracting any continuum constant from F does
not alter A. The comparator is the unbiased same-replica estimator
`|F(X)-F(Y)|^2/2` of `Var(F)`.

A 32-configuration exact certificate fixes bonds `140,168,196` open, all
other bonds closed except the five free bonds `0,28,56,84,112`. The fifth
mixed difference of F is exactly `-1/3`. Every degree<=4 polynomial has zero
such difference, so F is not low-degree. Because every configuration has
positive probability and h_j>0 for j>=5, A is strictly positive. This is an
algebraic existence result, not a lower bound of practically measurable size.

## Acquisition and results

X and Y are independent Bernoulli fields. Five shared fair retention masks
are cumulatively ANDed, giving nested retention probabilities
`1,1/2,1/4,1/8,1/16,1/32`; at each level retain X or replace by Y. All three
children share each resulting mask. The coefficients are exactly
`[1,-31,310,-1240,1984,-1024]`.

All 38 coordinates retain full 100-batch covariance. Errors below are one SE.

| Readout | Estimate | SE | Interpretation |
|---|---:|---:|---|
| High-pass self-energy A | -0.0122444 | 0.509505 | underpowered; one-sided p=.5096 |
| Unfiltered Var(F) | 0.0806472 | 0.000447155 | ordinary variance resolved |
| Imaginary self-energy | -0.902851 | 0.528138 | zero control z=-1.71 |
| K-centered Euler response | -0.214790 | 2.30619 | zero control z=-.093 |
| Degree-5 response | -8.05220 | 7.38096 | exact target .298004, deviation z=-1.13 |

The negative A estimate is not clipped. The population is positive; this
estimate simply has far too much acquisition noise. The Euler control is an
honest L3 site-Euler polynomial on the first nine noise bits, **not** a claim
to be the N112 physical Euler observer. Even the known degree-5 positive
control is not statistically resolved by this acquisition.

Wall time: **7.882 s**; summed worker CPU: **69.199 s**. Comparator
classification CPU: 19.517 s (excludes RNG/control/output overhead).
The high-pass response has **1,298,313 times** the comparator's sampling
variance. Squared-SNR/CPU-second is `8.35e-6` versus `1666.7`; their observed
ratio `5.01e-9` is noise-dominated and must not be interpreted as a precise
population efficiency estimate.

The exact bound `0<A<=Var(F)<=1/9`, combined with pilot variance, gives a
more useful optimistic envelope: even at the universal maximum A=1/9,
5-sigma readout projects to **10.51 million** replicas, and relative
efficiency is at most `4.12e-7` at these measured variances/costs. If A is all,
one tenth, or one hundredth of observed Var(F), the corresponding projections
are 19.96M, 1.996B, and 199.6B. These are fixed-estimator variance
extrapolations, not universal high-pass impossibility theorems. They do not
authorize new production.

## Lifecycle and provenance

- freeze: `7fd8aa0`; host amendment and exact execution code: `7a9ce54`;
- runner in raw metadata received the literal revision token `HEAD`; at
  execution this resolved to `7a9ce54`, with script SHA-256 recorded alongside
  raw files. Raw metadata is preserved rather than silently rewritten;
- parent baseline: main `eca7d4e`; prior N112 source `2402a33` not rerun;
- seed `43711231001`; counters `[43700000000,43700020000)`; 100 equal batches;
- host changed **before sampling** to local Mac after a separate seed0,
  counter0, 32-replica timing smoke took .0675 CPU seconds;
- HZ was started but rejected the registered SSH public key. No #437 remote
  job was launched and no key reset was attempted; existing cloud work is not
  counted as this dependency group;
- dependency group: `p437-N112-six-noise-fresh-20k-20260831`;
- status: fixed pilot complete; 4 focused tests pass; no PR created/merged.

Reproduce scoring:

```sh
python3 scripts/score_p437_high_pass_mc.py \
  results/local-20260831/P437-N112-six-noise-20k \
  --output results/local-20260831/P437-N112-six-noise-20k/score.json
```

## Scientific consequence

This closes the distinction between non-low-degree topology and its noisy
acquisition. It does not test a spatially independent JD-perp source, identify
a field, or count physical states. The next useful experiment should change
the estimator—positive/conditional mixed differences rather than million-fold
alternating noise cancellation—not increase this stream's sample count.
