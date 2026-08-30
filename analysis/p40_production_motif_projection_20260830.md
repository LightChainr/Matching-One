# Issue #40: production motif projection

## Scientific card

- **Question:** can the exact fixed-`K` four-motif subspace remove the common-field orientation noise in `Delta q`?
- **Frozen design:** fresh 1,000,000 replicas at each of `N=65,85`; 100 batches; five-fold out-of-fold projection on `Delta E_mc`, `Delta(diagonal_mc)`, `Delta F0_mc`, and the exact single-family `Delta right_angle_mc`.
- **Result:** held-out variance reduction is `1.3630x` at `N=65` and `1.3314x` at `N=85`; the controls explain `26.63%` and `24.89%` of replica variance, respectively.
- **Mechanism:** edge and face enter with stable opposite coefficients; diagonal is small and the single translated right-angle family is nearly null. The four-control space is real but low-rank with respect to the target coupling.
- **Decision:** do not expand this phase to the remaining three declared pairs. The exact covariance certificate was necessary, but this production projection misses the `>=2x` information-value threshold.

The performance values above are exclusively omitted-fold scores. Training covariance is used only to fit each fold's coefficients and is not reported as achieved variance reduction.

## Frozen acquisition and gates

The prereveal manifest is `experiments/p40_production_motif_projection_20260830.yaml`, frozen in commit `f6b7411`. The runner was fixed first at `1b988f9`; binary SHA-256 is `639b234d0215a4e66fec9ed7e83d3b933d677547906ad60aa43ff9dcf118739c`.

Both runs used seed `4020260830`, counters `[40000000,41000000)`, and identical cyclic-label Bernoulli fields within each same-`N` pair. Each output contains 100 contiguous batches and 1,000,000 replicas. Euler identity and wrapping-channel equality both have total L1 failure `0` at both sizes. The local Clang binary records OpenMP as false, so `threads_requested=8` did not parallelize; this affects timing only, not the frozen RNG or score.

| N | pair | runner wall (s) | raw Var(Delta q) | held-out projected | held-out residual | unexplained | OOF VR |
|---:|---|---:|---:|---:|---:|---:|---:|
| 65 | `(8,1)/(7,4)` | 8.793 | 0.720707 | 0.191949 | 0.528758 | 0.733666 | 1.363018 |
| 85 | `(9,2)/(7,6)` | 11.555 | 0.761209 | 0.189467 | 0.571736 | 0.751096 | 1.331388 |

The two runner wall times sum to 20.348 seconds. The conditional OOF adjusted means are `0.0005424 +/- 0.0007272` at `N=65` and `0.0006205 +/- 0.0007561` at `N=85`; these are estimator diagnostics, not new finite-size physics claims.

## Full target/control covariance

Order: `Delta q`, `Delta E_mc`, `Delta diagonal_pair_mc`, `Delta face_mc`, `Delta right_angle_mc`.

`N=65`:

```text
 0.72070665   0.88974830  -0.93756198  -0.26281026   0.00830952
 0.88974830  15.13775055  -3.79385350   9.29933827   6.72723751
-0.93756198  -3.79385350  15.12471809   2.64766265   2.23444660
-0.26281026   9.29933827   2.64766265  10.47159047   6.59249156
 0.00830952   6.72723751   2.23444660   6.59249156   7.13955714
```

`N=85`:

```text
 0.76120940   1.07489598  -1.02161027  -0.23179843   0.04488503
 1.07489598  19.77131805  -4.91661918  12.17667875   8.79027131
-1.02161027  -4.91661918  19.77631568   3.47824664   2.94537837
-0.23179843  12.17667875   3.47824664  13.72628365   8.63299866
 0.04488503   8.79027131   2.94537837   8.63299866   9.34973763
```

All five fold coefficient vectors, raw sums, complete joint moments, metadata, and hashes are in `results/local-20260830/P40-production-motif-projection/`. Coefficients are exceptionally stable across folds; representative averages are approximately `(0.1753, 0.0141, -0.1870, 0.0043)` at `N=65` and `(0.1562, 0.0152, -0.1597, 0.0006)` at `N=85`.

## Interpretation boundary

This experiment establishes a reproducible partial projection, not a failed control-variate idea in general. It specifically says the four geometry-only, fixed-`K` motif contrasts do not span most of the target noise at these two sizes. A further experiment would need a new target-coupled observable class, not more sizes with the same four controls.
