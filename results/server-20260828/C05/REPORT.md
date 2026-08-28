# C05 / P33 — Thermal-coordinate tomography from threshold ranks

Engine: issue-9 Philox Fisher–Yates Newman–Ziff on the C00 `HomologyUnionFind` (exact `adj(P)/det(P)` windings). **Cross channel only**; either-wrap matching-function rows are not a second replication.

## Off-by-one convention

- `K_plus`: smallest black occupation `k` at which primal CROSS wrapping is true (`N+1` if never). Transition at rank `k` in `1..N` is the `k`-th uniform order statistic, `T|K=k ~ Beta(k, N+1-k)`.
- `K_minus = N - m* + 1` where `m*` is the first reverse-permutation matching occupation with CROSS wrapping (`0` if white never wraps).
- Matching function: `M(p)=P(K_plus <= m) - P(K_minus > m)` with `m ~ Binomial(N,p)`.
- Exact tiny tests (axis L=2,3 exhaustive perms; gaussian (2,1); axis L=4 subset vs published polynomial + MC): **PASS**. `K_minus <= K_plus` on every replica of those tests and of production.

## Frozen choices (N=65 only)

- u-grid: 0, 0.05, 0.1, 0.2, 0.4 (from N=65 `|Mbar|` coverage; p*_65=0.59275127).
- Scaling train: N=65,85,130. Held-out: N=145,170. N=145 is retained even if noisy.

## Per-size summary (2e6 CRN replicas each, 40 batches, 8 OpenMP threads)

| N | L | p*_Mbar | P4[D](p*) | se_batch | P4[S](p*) | se_batch | kappa3 | mean gap |
|---|---|---|---|---|---|---|---|---|
| 65 | 8.0622577 | 0.59275127 | 0.00082052952 | 9.9542946e-05 | -0.00096128541 | 9.2408798e-05 | -1.5407662 | 5.5484525 |
| 85 | 9.2195445 | 0.59276167 | 0.00077392017 | 8.8145792e-05 | -0.00032170315 | 7.4964498e-05 | -1.5627697 | 6.6089907 |
| 130 | 11.401754 | 0.59270641 | 0.00023543056 | 0.00010654116 | -0.00024108767 | 0.00010206238 | -1.590351 | 8.7002303 |
| 145 | 12.041595 | 0.59274231 | 0.00014242438 | 7.9579479e-05 | -0.00019385277 | 8.0086951e-05 | -1.5944238 | 9.335648 |
| 170 | 13.038405 | 0.59276929 | 0.0002591985 | 8.2756334e-05 | -3.5427991e-05 | 7.7209193e-05 | -1.6035977 | 10.335261 |

Signed effects are histogram-pooled P4 values; `se_batch` is the batch SD/sqrt(B) of the 40 CRN batches (covariance across orientations is already in the paired batches).

## Target tests

P4[D]_thermal_even at u=0 vs L^{-13/4} = -3.25:

- pairwise effective exponents: N=65→85 w_eff=-0.43599707, N=85→130 w_eff=-5.6017877, N=130→145 w_eff=-9.2052756, N=145→170 w_eff=7.5287978
- mean w_eff = -1.9285657 (target -3.25)
- looks-like-L^{-13/4} (mean w_eff within 1 of target): **no**
- train A for L^{-13/4}: 0.80754307
- held-out RMSE L^{-13/4}: 8.8763045e-05
- held-out RMSE L^{-13/4}(1+B log L), B=-0.23333906: 8.7698828e-05
- held-out RMSE free exponent (train e=-3.7840243): 8.7421905e-05
- log alternative improves held-out P4[D]: **yes**

P4[S]_thermal_even at u=0 vs L^{-2}:

- pairwise w_eff: N=65→85 w_eff=-8.1609325, N=85→130 w_eff=-1.3578723, N=130→145 w_eff=-3.9938296, N=145→170 w_eff=-21.36988
- held-out RMSE L^{-2}: 0.00015498765
- held-out RMSE L^{-13/4}: 8.3870669e-05

P4[D]_even values: N65=0.00082052952, N85=0.00077392017, N130=0.00023543056, N145=0.00014242438, N170=0.0002591985
P4[S]_even values: N65=-0.00096128541, N85=-0.00032170315, N130=-0.00024108767, N145=-0.00019385277, N170=-3.5427991e-05

## Remaining

- Axis L=8..32 production histograms were deprioritized (P33 same-N Gaussian first); exact axis L=2,3,4 tests did run and PASS.
- No 1e8-sample campaign. CPU pilot is 2e6/orientation-pair.
- N=1105 not started.
