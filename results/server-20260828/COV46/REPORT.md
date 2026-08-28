# COV46: hardened threshold-rank cross-size covariance audit

Status date: 2026-08-28 21:48 CST
Branch: `agent/cov46-harden` (not merged; original PR #46 research branch untouched)
Replay input: committed P33 histogram `results/server-20260828/P33/p33_all_sizes_arm_10m.hist.csv`
Original audit retained (not overwritten): `results/server-20260828/P33-cross-size-covariance/`

This is a hardening of the scoring implementation. Frozen predictions were not
retuned. Simulators and RNG were not changed. No N=185/265 or 1e8 campaign was
started.

## Four robustness gates

| Gate | Result | What was implemented |
|---|---|---|
| 1. Equal-batch-weight contract | **PASS** | Core `_orientation_batches()` now fails unless every aligned size, orientation, and batch shares one sample count. `covariance_of_mean` remains an unweighted equal-batch estimator and is no longer reachable on mixed weights. `--seed-label` is no longer declarative: CLI coupling requires metadata proving shared RNG schema, seed `2026093303`, and replica counters `[4000000000, 4010000000)` with per-batch ranges `first + batch * 100000`. |
| 2. Near-singular covariance | **PASS** | Scientific inverses are Cholesky-first with SVD/eigen pseudoinverse fallback. Eigenvalues, condition number, effective rank, and an eigenvalue-truncation sensitivity table are emitted on every score. Gauss-Jordan remains only as a cross-check. |
| 3. Finite-batch calibration | **PASS** | The quadratic form is labelled `plugin_asymptotic`. Each held-out score carries Hotelling/F (n=100, p=2) and an uncentered 2000-replicate batch bootstrap. These are diagnostics, not paper-facing p-values. |
| 4. Stronger synthetic regression | **PASS** | Exact analytic cases: identity covariance full=diagonal=6/5; equicorrelated 2-heldout full=2 vs diagonal=1; 3-size full=4/3 vs diagonal=2/3. Unequal cross-size samples fail closed. Rank-1 training covariance uses the SVD pseudoinverse. |

Covariance-related tests: **21 passed**. Full repository suite: **76 passed**.

## Empirical P33 replay (unchanged scientific conclusion)

Cross-size correlations remain real but modest (maximum absolute correlation about 0.220). Full vs diagonal held-out plugin chi-squares agree with the retained original archive to all printed digits, and the Cholesky score matches the Gauss-Jordan cross-check.

| metric | original diagonal | original full | hardened diagonal | hardened full | solver |
|---|---:|---:|---:|---:|---|
| A_M held-out plugin chi-square / 2 | 5.27335 | 5.53002 | 5.27335 | 5.53002 | cholesky |
| A_p held-out plugin chi-square / 2 | 5.29288 | 5.55159 | 5.29288 | 5.55159 | cholesky |
| root doubling joint residual / 2 | 3.44173 | 3.46246 | 3.44173 | 3.46246 | (from jackknife root-gap covariance) |

A_M training amplitude (full covariance): 0.529773 +/- 0.111119; held-out N=145,170 observations 1.08559, 1.53951.
A_p training amplitude (full covariance): 0.302567 +/- 0.0634742; held-out observations 0.621437, 0.881604.

Maximum |rho|: A_M 0.2204; A_p 0.2204.

## Eigenstructure

A_M covariance of means: eigenvalues [0.0150994, 0.0307078, 0.113252, 0.234373, 0.288523]; condition 19.11; effective rank 5/5.
A_p covariance of means: eigenvalues [0.00492313, 0.0100293, 0.0371128, 0.0767712, 0.0946202]; condition 19.22; effective rank 5/5.
A_M full held-out residual: eigenvalues [0.130467, 0.31069]; condition 2.381; effective rank 2/2.
A_p full held-out residual: eigenvalues [0.0427489, 0.101864]; condition 2.383; effective rank 2/2.

Neither 5x5 mean-covariance nor the 2x2 residual covariance is near-singular on this archive (condition ~19 and ~2.4). The SVD fallback is therefore unused on P33, but it is the path used for the rank-1 synthetic gate.

A_M full truncation sensitivity: cut=0 -> chi2=5.53002 (rank 2); cut=1e-14 -> chi2=5.53002 (rank 2); cut=1e-12 -> chi2=5.53002 (rank 2); cut=1e-10 -> chi2=5.53002 (rank 2); cut=1e-08 -> chi2=5.53002 (rank 2); cut=1e-06 -> chi2=5.53002 (rank 2).
A_p full truncation sensitivity: cut=0 -> chi2=5.55159 (rank 2); cut=1e-14 -> chi2=5.55159 (rank 2); cut=1e-12 -> chi2=5.55159 (rank 2); cut=1e-10 -> chi2=5.55159 (rank 2); cut=1e-08 -> chi2=5.55159 (rank 2); cut=1e-06 -> chi2=5.55159 (rank 2).

## Finite-batch calibration

The numbers above remain plug-in/asymptotic chi-squares with covariance estimated from 100 batches. They are not calibrated p-values.

A_M full: Hotelling F=2.737 on df=[2, 98]; F survival=0.06971; plugin chi2 survival=0.06298 (null calibration diagnostic, not paper-facing).
A_M full: uncentered batch bootstrap n=2000 seed=20260828: q50=6.918, q90=15.93, q95=18.7, q99=24.97; fraction >= observed 0.616 (sampling variability, not a null p-value).
A_p full: Hotelling F=2.748 on df=[2, 98]; F survival=0.06901; plugin chi2 survival=0.0623 (null calibration diagnostic, not paper-facing).
A_p full: uncentered batch bootstrap n=2000 seed=20260828: q50=6.937, q90=15.96, q95=18.73, q99=25.02; fraction >= observed 0.616 (sampling variability, not a null p-value).

Hotelling survival ~0.07 agrees with the asymptotic chi-square tail (~0.063) to the precision expected for n=100, p=2. The uncentered bootstrap places the observed statistic near the middle of the resampled distribution, so the modest full-vs-diagonal difference is not an artifact of a handful of batches. A centered/null parametric bootstrap is unnecessary given the Hotelling conversion; the empirical bootstrap here is a variability diagnostic.

## Coupling contract actually validated

Histogram sidecar metadata was required before treating equal batch ids as coupled:

- engine: same-N Gaussian threshold-rank Newman-Ziff
- rng: counter-derived SplitMix64 stream plus unbiased Fisher-Yates
- seed: 2026093303 (matches `--expected-seed`)
- replica counters: [4000000000, 4010000000)
- 100 batches x 100000 samples; per-batch range `4000000000 + batch * 100000`
- designs cover N=65,85,130,145,170
- all 1000 aligned (N, orientation, batch) cells share samples=100000

Without that metadata the CLI now exits rather than emitting a coupled score.

## Conclusion

The four merge-blocker gates pass. The P33 scientific takeaway stands: cross-size |ρ| ≲ 0.22 and full vs diagonal chi-squares are almost the same, so ignored coupling was not the source of the low-stat radial tension. The score is now fail-closed on batch weights, RNG coupling, rank, and finite-batch labelling.

Do not merge PR #46 from this agent branch; this is the hardening implementation for review.
