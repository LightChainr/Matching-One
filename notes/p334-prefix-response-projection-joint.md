# P334: contact structure adds signed loading beyond source strength

The original8 readout now shares the complete original-20-batch covariance
factor with the prior hierarchy. This is the covariance-coordinator item of
`93ee4e98`, consuming the frozen projection `90226598` and prior joint
`ce20158a`. No projection is refitted here; all nonlinear uncertainty comes
from the producer's supplied delete-one-original-batch rows.

## Main result

Four exact contact descriptors capture substantially more **own-source
signed covariance loading** than source score energy alone. The paired
increment is 31–46 percentage points, with the uncertainty below computed
jointly, not by adding the two projection errors independently.

| N / physical receiver | Energy share | Four-contact share | Paired gain |
|---|---:|---:|---:|
| 325 / first | 51.68 ± 8.05% | 90.01 ± 9.29% | 38.32 ± 8.39 pp |
| 325 / second | 49.09 ± 9.25% | 80.28 ± 8.08% | 31.19 ± 9.01 pp |
| 425 / first | 51.86 ± 7.63% | 97.02 ± 6.97% | 45.16 ± 9.32 pp |
| 425 / second | 53.10 ± 7.60% | 99.00 ± 10.05% | 45.90 ± 6.66 pp |

Here the loading is
`2 Cov(mu_C,H_own C) - (1/2) Cov(mu_W,H_own W)`, pooled within rank cells
whose physical receiver is rank0, with the full 20,000-prefix denominator.
Energy-only residuals are `(4.331±1.209, 4.389±1.115, 3.978±0.813,
3.575±0.892) × 10^-8` in table order. Four-contact residuals become
`(0.896±0.920, 1.700±0.779, 0.246±0.591, 0.076±0.774) × 10^-8`.
Thus the strongest remaining loading discrepancy is N325/second, about
2.18 SE; a zero residual elsewhere is not an established closure.

## What remains after both baseline clocks

The clocks `mu_C,mu_W` reproduce their own clock-response loading by linear
projection identity. The non-tautological remaining readout is the partial
contact-response cross-moment

`Cov(Z,H | clocks)_linear = Cov(Z,H) - Cov(Z,clocks) Cov(clocks)^-1 Cov(clocks,H)`.

At N425/second receiver and physical second source, its C-response entries are:

- joint-safe mass: `(1.07475 ± 0.23633) × 10^-7`;
- source score energy: `(1.80780 ± 0.56397) × 10^-9`;
- own-safe degree: `(2.38216 ± 0.98852) × 10^-7`;
- own-safe loop: `(-0.22485 ± 3.19618) × 10^-8`.

The same own-C safe-mass entry is positive at roughly 2.4 SE in the other
three receiver/size cases. These descriptors are correlated views of the
same prefix ensemble. All 32 partial entries per size, including cross-source
and W responses, remain in the result rather than selecting only these rows.

The contact-after-clock projected response-variance increment is largest
relative to its SE for N425 second→second C:
`(4.80649 ± 1.74383) × 10^-11`. The other own-C increments are
`(3.43435±2.21289, 4.01394±2.67584, 2.67822±1.91481) × 10^-11`
for N325 first, N325 second and N425 first. All W and cross-source increments
are reported in the full table. These are in-sample plug-in projection
energies; their positivity is not a confidence lower bound on physical
response variance. No R-squared or additional omnibus test is assigned.

## Scientific interpretation and boundary

Within this fixed low-dimensional readout, contact composition carries
substantial signed loading beyond source strength. There is also a residual
own-C association after removing both baseline clocks, most clearly resolved
for N425/second. This is an exploratory relationship among prefix conditional
means and responses; it does not establish causality, a sufficient state,
a field identity, or out-of-sample prediction. Signed loading shares are not
variance fractions. No claim is made that loop contact alone vanishes.

The original8 stream is the sole input. The new64 reuse task has not been
read or pooled. Every receiver, physical source and C/W direction is kept
in one covariance factor with all preceding same-block coordinates.

## Reuse and scientific card

- Observer / sector: receiver-R0 within-rank-cell prefix covariance;
  physical-source C/W response, unrotated first/second coordinates.
- Source: original `e32a8593` forks, exact-score moments `375cd3a1`, exact
  descriptors `1cfa4ae8`; frozen projection `90226598`.
- Dependency group: the same original 20 paired batches per N as `ce20158a`;
  no new sampling, DP, raw replay, or independent evidence.
- Changed mechanism space: scalar source energy alone leaves substantial
  loading; four-contact structure accounts for much of this remaining loading.
- Not established: exact or predictive closure; positive confidence lower
  bounds for projection energy.
- Next observation: source-fixed conditional-stream reuse can distinguish
  persistent clock/contact cross-moments from finite-tail estimation noise;
  it must retain shared-prefix dependence and is handled separately.

Artifacts: [complete tables](../results/p334-prefix-response-projection-joint/REPORT.md),
[score and focused covariance](../results/p334-prefix-response-projection-joint/score.json).
The two compressed complete factors each contain 21,842 coordinates, including
all 358 supplied projection coordinates and 726 raw sufficient-statistic
coordinates, on the original 20 deleted-batch rows. Covariance is `F.T @ F`;
its rank is at most 19 and no inverse is taken.

Reproduce the thin join (no producer rerun):

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_prefix_response_projection_joint.py
```
