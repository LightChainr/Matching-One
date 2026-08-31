# A fixed local fifth-order topology response is readable

The new positive endpoint is resolved in **20,000 fresh replicas**, without
scanning supports or trying to reconstruct the whole high-pass energy:

`B_S = (3.23893 +/- 0.21639) * 10^-6` (one batch SE).

Mean/SE is 14.97. There are **384 nonzero backgrounds (1.92%)**, spread across
all 100 batches. This is a measurable localized high-order topology response
under the original square-bond product measure, not merely an exceptional
configuration witness.

## Exact relation to the original high-pass energy

The support was fixed before acquisition to `S={0,28,56,84,112}`, the same five
bonds as the exact witness in `7fd8aa0`; it was not chosen from this sample.
On the three N112 rho children, retain

`F=(E0+zeta^-1 E1+zeta^-2 E2)/3`, `E=1_rank0+1_rank2`.

All other 219 bonds are independently Bernoulli(1/2). For each background,
enumerate all 32 support configurations and compute the normalized derivative

`D_SF=2^-5 sum_u (-1)^(5-|u|) F(u,X_outside)`.

In the orthonormal Rademacher expansion, `D_S` removes S from each character
that contains it. Consequently

`B_S=E|D_SF|^2=sum_{T contains S}|Fhat(T)|^2`.

For j>=5, `h_j>=h5=9765/32768`, so the **exact population inequality** is

`A_HP >= h5 B_S`.

There is no `C(224,5)` factor. The old witness has `D_SF=-1/96`, squared
magnitude `1/9216`, confirming normalization. B_S collects all Fourier orders
containing this support; it does not isolate a pure degree-five field.

## Measurements

| Quantity | Estimate | Batch SE |
|---|---:|---:|
| B_S | 3.2389323e-6 | 2.1639031e-7 |
| h5 B_S | 9.6521527e-7 | 6.4485211e-8 |
| Nonzero-background rate | .019200 | .00086082 |

The exact sample mean is `B_S=199/61440000`. The estimate of the population
lower-bound parameter is `h5 B_S=129549/134217728000`.
**Neither sample estimate is a statistically certain numerical lower bound.**
The population inequality is exact; its measured right side has Monte Carlo
uncertainty as reported above.

## Which fixed-support response is active?

Raw classes record the three integer signed sums `(n0,n1,n2)` before dividing
by 32. Across this block:

| Nonzero child | Backgrounds | Share of sampled C3 derivative energy |
|---|---:|---:|
| child0: 2omega | 360 | 514/597 = 86.10% |
| child1: omega/2 | 13 | 33/597 = 5.53% |
| child2: (omega+1)/2 | 11 | 50/597 = 8.38% |

No sample activates more than one child derivative. Thus the observed
nonzero complex derivatives lie on the signed character rays of the active
child. This is an observation in 20k, **not an exact prohibition of
coactivation**. Every integer value class and its batch count is preserved.

The chosen edge support is consecutive along the child0 horizontal cycle;
the same edge indices have different placement in the other children.
The 86% concentration is therefore a **fixed-support localized response**,
not a child-symmetric experiment and not evidence of a global C3 preference.
No rotated support was acquired or selected.

## Cost and relationship to the preceding pilot

- local wall time: **36.823 seconds**;
- summed worker CPU: **333.140 seconds**;
- 20k independent backgrounds x 32 exact support states x 3 child classifiers;
- all eight primary/child-quadratic coordinates retain full 100-batch covariance;
- raw responses use integer energy numerators divided by 36864, avoiding
  floating-point cancellation in the primary values.

This endpoint is a localized contribution, not an unbiased replacement for
the total A_HP estimator. Its resolution does not invalidate the old six-noise
pilot's stop decision. Nor should an efficiency ratio between the two
different estimands be advertised.

## Provenance and stop

Freeze/execution commit: `3c3fc574dc665a915e54f7897cf22848bb8b36ab`.
Parent derivation: `79988f8`; original witness: `7fd8aa0`.

Seed `43753111201`; counters `[43710000000,43710020000)`; 100 equal batches;
16 workers; local Mac only. A separate seed0/counter100000 timing smoke used
32 backgrounds and .299 CPU seconds, passing the frozen cost gate. No Huawei
job was launched. New dependency group:
`p437-N112-fixed-S5-lower-bound-fresh20k-20260831`.

Reproduce readout:

```sh
python3 scripts/score_p437_fixed_support.py \
  results/local-20260831/P437-N112-fixed-S5-20k \
  --output results/local-20260831/P437-N112-fixed-S5-20k/score.json
```

The two focused normalization/bound tests pass. Fixed 20k is complete; no
support rotation, sample extension, new PR, or merge was performed.
