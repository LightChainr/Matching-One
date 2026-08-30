# C01--C03 same-N channels and frozen covariance evaluation

The general same-N engine now retains `cross`, configuration-level `both`,
`either`, `direction_0`, and `direction_1` for both the primal and
complementary matching graphs. Exhaustive `N=5,13` union-find/BFS comparisons
pass for every channel.

A 200,000-replica pilot (`seed=20260828`, counters `[0,200000)`) froze the
five-channel matching-odd GLS weights. Independent evaluation used counters
`[200000,2200000)` for each of `N=65` and `N=85`.

The most important C03 result is exact and negative: for this matching
construction, all five

```text
D_channel = (R_primal,channel - R_matching,channel) / 2
```

observables are identical configuration by configuration. The pilot therefore
froze symmetric weights `[0.2,0.2,0.2,0.2,0.2]`; optimized and equal-weight
evaluation variances are exactly equal. The measured variance-reduction ratio
is `1.0` at both sizes, so the C03 `>=2x` GPU gate does not pass.

For the matching-function orientation difference, the independent 2-million
evaluation gave

| N | difference | batch SE | z |
|---:|---:|---:|---:|
| 65 | `+1.3560e-3` | `6.20e-4` | `2.19` |
| 85 | `+4.3950e-4` | `6.46e-4` | `0.68` |

These lower-statistics values are compatible with the independent 30-million
run already stored under `../gaussian/` (`+1.0038e-3 +/- 1.68e-4` and
`+7.6033e-4 +/- 1.58e-4`). They do not add a stronger discovery claim.

The C02 matching-even dominance hypothesis is not supported. Even-sector
orientation differences depend strongly on the wrapping channel and may flip
sign between `cross` and `either`, while the matching-odd difference is
channel-identical. A held-out multi-angle confirmation and an exact angular
control are still required before calling the earlier matching-function signal
a confirmed spin-4 law.
