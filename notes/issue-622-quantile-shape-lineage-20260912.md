# #622 after #702: reflection residual and finite-size shape motion

Date: 2026-09-12. Base: `main` at `6edf775e` (includes #702–#704).
Machine: XPk2PZ (16 vCPU ARM, no GPU), job under `/workspace/mo-622/`.
Artifacts: `results/probe-invariant-shape/quantile-shape-lineage-622.json`,
`scripts/probe_invariant_shape/quantile_shape_lineage_622.py`,
`tests/test_quantile_shape_lineage_622.py`, `results/probe-invariant-shape/run-622.log`.

C2 reanalysis of already-committed blocks. No Monte Carlo, no GPU, no exponent
fit, no `docs/STATUS.md`, no issue closed.

## Object and order of operations

For each size the two orientations are inverted **separately** and only then
combined, `Q_spin0 = w_1 Q_1 + w_2 Q_2` with `w_1+w_2=1` and
`sum w_i cos4theta_i = 0` (checked to `6e-17`); CDFs are never mixed before
inversion. Then, on the nine frozen deciles with anchors `a=0.2`, `b=0.8`,

    W_N = Q_N(0.8) - Q_N(0.2) > 0
    Z_N(u) = [Q_N(u) - Q_N(0.2)] / W_N
    A_N(u) = Z_N(u) + Z_N(1-u) - 1

with `A` the corrected normalized-shape reflection residual (not
`M_N(p)+M_N(1-p)`). Every delete-one batch repeats the whole map; the reported
covariances are of one random object across the grid. Independent coordinates
`u = 0.1,0.3,0.4,0.5` are used because `A(0.2)=A(0.8)=0` and `A(1-u)=A(u)`
exactly; `Z(0.2)=0`, `Z(0.8)=1` are likewise fixed. Angles/covariance use the
`cos 4theta` convention of `scripts/threshold_quantile_lineage.py`
(blob `3b328f29`, the #655 lineage); the equal weighting is a sensitivity on the
same block, not a second experiment.

## Lineage, geometry, weights

Same channel (`rank-2 cross wrapping`) and same observable (`K_plus`, `K_minus`)
at all three sizes; the exact tiny site/bond labs are a different observable and
are not mixed in.

| N | reps (first, second) | shortest period `sqrt(N)` | seed | spin-0 weights |
|--:|:--|--:|--:|:--|
| 145 | (12,1), (9,8) | 12.0416 | 2026105003 | (0.50711806, 0.49288194) |
| 290 | (13,11), (17,1) | 17.0294 | 2026105004 | (0.50711806, 0.49288194) |
| 725 | (26,7), (23,14) | 26.9258 | 2026105011 | (0.53827771, 0.46172229) |

`site_count = |a+bi|^2`, so the shortest lifted period holds `sqrt(N)` sites.
All three spin-0 combinations are interpolations; distinct seeds mean the
cross-size covariance term is taken as zero on random-stream provenance.

## Reflection residual `A = 0` (spin0 primary)

| N | W (SE) | A(0.1) | A(0.3) | A(0.4) | A(0.5) | `\|A\|` (SE) | nominal chi2, 4 dof |
|--:|:--|--:|--:|--:|--:|:--|:--|
| 145 | 0.11921330 (3.7e-8) | -0.0204652 | 0.0093537 | 0.0138721 | 0.0152436 | 0.0305142 (2.6e-5) | 2.27e6 |
| 290 | 0.09220375 (2.8e-8) | -0.0159133 | 0.0072317 | 0.0107147 | 0.0117706 | 0.0236407 (2.9e-5) | 1.00e6 |
| 725 | 0.06553076 (2.4e-8) | -0.0113451 | 0.0051435 | 0.0076198 | 0.0083705 | 0.0168313 (3.5e-5) | 3.56e5 |

`p` is below double reference at every size. The tests are nominal Gaussian
references with estimated delete-one covariance. The 4x4 `A` covariance is
near-collinear (correlations up to 0.9999; condition number 2.4e8 / 6.3e7 /
1.7e7), so the full-inverse chi-square is quoted next to the correlation-free
diagonal chi-square; the full value is **smaller** than the diagonal value at
every size and every displacement test, so the rejection is carried by the
per-coordinate residuals (max `|t|` = 1354, 913, 525), not by inverting a
near-null direction. No pseudoinverse is used anywhere.

## Adjacent-size motion (pooled, spin0)

| pair | `\|dZ\|` (SE) | `\|dA\|` (SE) | chi2_Z, 7 dof | chi2_A, 4 dof |
|:--|:--|:--|:--|:--|
| 145→290 | 0.00537503 (2.1e-5) | 0.00687410 (3.9e-5) | 1.37e5 | 4.51e4 |
| 290→725 | 0.00499749 (2.7e-5) | 0.00680943 (4.5e-5) | 6.42e4 | 3.27e4 |

Both displacements are resolved from zero; `Cov(dZ)=Cov(Z_from)+Cov(Z_to)`.
The declared norm is the Euclidean `L2` on the independent coordinates.

Interval-to-interval change: `|dZ|` falls by `3.78e-4 +/- 3.45e-5` (10.9σ,
resolved); `|dA|` changes by `-6.47e-5 +/- 5.95e-5` (1.1σ, **not** resolved).

## Equal-weighting sensitivity

Equal weighting retains a spin-4 residue (`-1.37e-2`, `+1.37e-2`, `-4.11e-2`).
The difference from spin0 is at most `1.6e-5` in `A` and `8.7e-6` in `Z`. The
qualitative reading is unchanged; it is one analysis, not two.

## Measured profile

Load of the three blocks 1.0 s; N=725 pooled pass 1.05 s; one N=725 delete-one
batch 1.02 s. Full delivery (3 sizes × 2 weightings, 100 batches each,
single-threaded) 5 min 23 s wall, about 0.09 CPU-hours — the 4 CPU-hour cap is
not binding, so the complete delete-one path was delivered rather than a
substituted estimator. N=725 pooled `Q` reproduces the #655 decile values to
`4.1e-9` (that reference is published to 8 decimals).

## Outcome

* **Reflection residual: resolved.** `A != 0` at all three sizes at hundreds to
  >1000σ per coordinate; the normalized shape is not reflection-symmetric.
* **Size trend on this finite lineage: decreasing.** `|A| = 0.03051 -> 0.02364
  -> 0.01683`, each step resolved (178σ, 150σ). Three sizes cannot establish
  convergence or a limit; this is observed finite movement, and no exponent is
  fitted.
* **Adjacent-size shape change: resolved.** Non-affine displacement remains
  between adjacent sizes (`|dZ|`, `|dA|` both nonzero). Between the two
  available intervals `|dZ|` decreases slightly (10.9σ) while `|dA|` is
  consistent with constant (1.1σ).
* Weighting does not change the qualitative reading; no covariance-aware test
  is decided by a near-null direction.

No conclusion here identifies an exponent, a field, or a threshold; a larger-N
purchase would need a separate decision naming which surviving alternative it
separates, and none is authorized here. #622 stays open.
