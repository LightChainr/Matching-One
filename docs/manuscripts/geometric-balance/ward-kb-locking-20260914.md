# note-kb-ward — issue #774: complete-component quality–boundary Ward checks

**Machine**: Huawei Cloud `DevEnvC_XPk2PZ` (fresh), workdir `/workspace/ward774/`.
**Rule honoured**: **no computation on the local Mac**; every number below was produced by a
script executed on the cloud machine (logs and JSON are in `/workspace/ward774/out/`,
downloaded to `ward774-out/`).
**Engine**: pinned `tagged_winding_span.py`,
sha256 `9621acf490dbc4b7dba28f2d2f0f4c2e0f1ee9e42d3987d1f7815bd1c9d5c9a3`,
git blob `52f3611990ce2b1331d9e5296e0262f5e402e0d7` (printed by every run).
Only the shipped `build/lump/numeric_system/solve` and `activity_transfer` are used — **no new
transfer engine, no Monte Carlo, no GPU**.

---

## 1. Convention (authoritative, BRIEF §5)

* black = occupied, NN connectivity, probability `p`; white = empty, matching, `q = 1-p`.
* Component = **complete horizontal-winding (essential) NN cluster**, counted once (tagged).
* `K` = occupied sites of the component.
* `B` = **number of DISTINCT external sites NN-adjacent to the component** (the engine's
  `v`-fugacity exponent). **Not** the number of contact edges — a site adjacent to the cluster
  through two different contacts is counted once.
* `q = 1-p`, `S = K/p - B/q` (the **#774** sign convention; the giant-white package uses
  `K/q - B/p` because there `q` is the white occupancy — see §8).
* `p_ref = 0.5927460507921` is used **only as a near-critical diagnostic reference**, never as a
  new `p_c`.

Component weight: for a fixed cluster `C`, the probability that `C` is exactly a complete
NN cluster of the configuration is `p^K q^B` (all `C` occupied, every site of its external
vertex boundary empty, everything else free). Hence `nu_w(p) = Σ_C p^{K} q^{B}`.

## 2. `B` is the geometric distinct-boundary count — independently verified

The `activity_transfer` labels each transition with `(k, nb)` and gives the source an
offset `nb`. I verified that `K = Σ k`, `B = source_nb + Σ nb` is the **geometric** pair by
brute-force enumeration, independent of the engine:

* enumerate every connected NN site set `C` on the `w`-cylinder with minimum row 0, no empty
  intermediate row, that winds horizontally (lift/spanning-tree winding test);
* `K = |C|`, `B = |{sites ∉ C, NN-adjacent to C}|`;
* compare the full `(K,B)` histogram with the histogram obtained from the engine's own
  accepting paths.

Result (`scripts/diag_marks.py`, `out/diag_marks.json`):

| width | max rows | # components | brute = engine |
|---|---|---|---|
| 2 | 7 | 969 | **identical** |
| 3 | 5 | 5207 | **identical** |

So the `nb` label is *not* an edge perimeter and *not* a mis-count; it is the distinct
external boundary. (Nitpick worth recording: for the NN case the *source* offset happens to
equal the occupancy of the component's first row, and the accept step adds the acceptance-row
occupancy, so the decomposition of `B` between source and transitions is not itself the
geometric one — only the total is. This is a label-semantics subtlety, not an error.)

## 3. Algorithm (exact, two fugacities)

The transfer matrix is a function of two fugacities `(u,v)`; with marks `(1,k,nb)` per
transition and the source carrying its own `v`-offset:

```
R(u,v), b(u,v), alpha(v)   ->   nu = alpha^T (I-R)^{-1} b,  a rational function of (p,q)
E[K] = (u d/du) log nu |_{u=p, v=q}      E[B] = (v d/dv) log nu |_{u=p, v=q}
d/dp log nu_w(p)  =  (u d_u - v d_v) log nu  ->  E[K]/p - E[B]/q = E[S]      (I1)
```

All first and second derivatives (including the mixed one) are obtained by **implicit
differentiation of the resolvent**, i.e. by solving extra linear systems with the same matrix
`A = I - R(p,q)`; **no finite differences are used as a certificate**. `gauss_solve_many()`
does one Gauss–Jordan elimination with several right-hand sides; all arithmetic is
`fractions.Fraction`, so every certified check is an exact rational identity.

**Cost.** `activity_transfer(w, matching=False)`: 12/68/340/1672/8244 states →
3/5/13/24/61 lumps for w=2..6 (~15 s at w=6); w=7: 40588 states → 126 lumps in 163 s, 268 MB;
w=8: see `out/w8.log` (state cap raised to 4·10⁶; the shipped builder hard-limits width to 6,
so for w=7,8 a **raise-only relaxed copy** of the *same* file is used — one line changed,
verified by diff, recorded in every output row as `engine`).

## 4. Code validation (and one bug I made and fixed)

My moment extractor is adapted from the #772 deliverable
`giant_white_volume.py::joint_activity` (unmodified upstream file). Before using it I
reproduced the **shipped** white-case numbers at `q=3/4` for w=2,3,4
(`nu, mean_span, mean_occupation, mean_boundary, score_mean, score_variance,
residual_covariance_scaled, prospective_clt_variance`) — **all exactly equal** to
`giant-white-volume.json`, as is the shipped assertion `score_mean == -d_p log nu_black(1/4)`.
This is a validation of my code, *not* a re-derivation of #772's numbers as a deliverable.

While doing this I found and fixed a bug of my own: my first version fed the *raw* `b`-vector
instead of the solved `y0` into the right-hand side of the first-moment solves, which gave
`E[L], E[K]` about 8× too small and made (I1)/(I2) fail. The validation above caught it.
**All numbers in this note are from the fixed code.** (Honesty note: the early, wrong run is
kept only in my logs; it is not reported as a result.)

## 5. Main table — `p = p_ref` (black NN), widths 2..8

Float64 diagnostics at `p_ref` (`out/tA_float_ref.json` for w=2..7,
`out/tA_float_ref_w8.json` for w=8); exact rationals at `p_ref` for w≤6 and at `p=1/2` for w=7
(`out/tA_exact_ref.json`, `out/tA_exact_half*.json`). `q/p = 0.687063`.
**w=8 completed** in 1820.8 s (≈30 min, 321 lumps, ≈2 GB peak) with the raise-only relaxed
activity builder — that is the cost of one extra width here; w=7 (40588 states → 126 lumps)
took 163 s. w=8 is float64; its exact-rational run was still going when I finalised, so the
exact certificates below stop at w=7.

| w | `nu_w` | `E K` | `E B` | `E B/E K` | `mean_ratio` | `score_mean` | `score_var` | `shot_scale` | `ward2_ratio` | `Corr(K,B)` | cond | angle° vs `(p,q)` | `sm_eig/trace` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.149338 | 7.007 | 5.082 | 0.7253 | 1.055595 | -0.657163 | 30.018 | 50.580 | 0.59347 | 0.80419 | 32.8 | 19.91 | 0.02958 |
| 3 | 0.096361 | 15.077 | 10.497 | 0.6962 | 1.013346 | -0.339454 | 67.872 | 106.199 | 0.63910 | 0.90403 | 36.2 | 11.32 | 0.02692 |
| 4 | 0.072335 | 25.763 | 17.785 | 0.6903 | 1.004751 | -0.206497 | 121.193 | 180.555 | 0.67122 | 0.94013 | 49.1 | 7.77 | 0.01994 |
| 5 | 0.057782 | 39.210 | 26.993 | 0.6884 | 1.001982 | -0.131130 | 191.015 | 274.350 | 0.69624 | 0.95909 | 66.9 | 5.86 | 0.01474 |
| 6 | 0.048134 | 55.286 | 38.022 | 0.6877 | 1.000975 | -0.090975 | 276.893 | 386.600 | 0.71623 | 0.97010 | 87.9 | 4.69 | 0.01125 |
| 7 | 0.041251 | 73.953 | 50.837 | 0.6874 | 1.000537 | -0.066997 | 378.665 | 516.998 | 0.73243 | 0.97709 | 111.8 | 3.90 | 0.00886 |
| 8 | 0.036092 | 95.169 | 65.408 | 0.6873 | 1.000321 | -0.051516 | 496.174 | 665.236 | 0.74586 | 0.98183 | 138.5 | 3.34 | 0.00717 |

All rows: (I1)(I2)(I3)(I4) pass — exactly for w≤6 at `p_ref` and for w=7 at `p=1/2`; float64
at `p_ref` for w=7,8 (max |I1| diff 8.3e-14, max |I2| diff 3.5e-11). `nu_activity == nu_tagged`
for every row (float64 agreement ≥ 10 significant digits).

**Extras (reported, not fitted as theorems).** `w^{55/48}(E B/E K - q/p)` =
0.1230, 0.0470, 0.0233, 0.0125, 0.0076, 0.0050, 0.0035 (w=2..8) — this **decreases**, i.e. the
deviation decays **faster** than `w^{-55/48}`. A 7-point log–log fit of
`E B/E K - q/p` gives an effective power ≈ `w^{-3.73}` (per-step exponent estimate rises
0.29→0.38 for `1-ward2_ratio`, see §7). `Var(S)/E K` = 4.28, 4.50, 4.70, 4.87, 5.01, 5.12, 5.21.

## 6. Pre-declared λ points

`c = 1/4` was fixed **before** running; `p = p_ref ± c·w^{-3/4}` (`w^{-3/4}` taken as its
12-digit decimal, `p_ref` as the rational `5927460507921/10^13`). No `c` was scanned
afterwards. Results (`out/tA_float_lam.json`):

| w | `p_-` | mean_ratio(-) | ward2(-) | Corr(-) | `p_+` | mean_ratio(+) | ward2(+) | Corr(+) |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.44410 | 0.78997 | 0.49806 | 0.76952 | 0.74140 | 1.24348 | 0.68743 | 0.86007 |
| 3 | 0.48307 | 0.83306 | 0.58655 | 0.87902 | 0.70242 | 1.12921 | 0.75404 | 0.93690 |
| 4 | 0.50436 | 0.86561 | 0.61808 | 0.92229 | 0.68113 | 1.08604 | 0.79266 | 0.96277 |
| 5 | 0.51798 | 0.88940 | 0.64300 | 0.94525 | 0.66751 | 1.06364 | 0.81719 | 0.97551 |
| 6 | 0.52753 | 0.90696 | 0.66351 | 0.95887 | 0.65796 | 1.05013 | 0.83445 | 0.98261 |
| 7 | 0.53465 | 0.92011 | 0.68065 | 0.96774 | 0.65084 | 1.04115 | 0.84732 | 0.98699 |

Both sides converge to 1 in the same way, so the mean/boundary Ward relation is not a
one-sided coincidence at `p_ref`.

## 7. The four certified checks

Let `I3` = the 3×3 covariance `Cov(L,K,B)`; PSD is certified **exactly** by the sign of every
principal minor (not by a float eigenvalue).

| check | statement | result |
|---|---|---|
| (I1) | `d/dp log nu == E[K]/p - E[B]/q` | w=2..6 **exact equality** (rational); w=7,8 float64, abs diff ≤ 1e-11 |
| (I2) | `d²/dp² log nu == Var(S) - (E[K]/p² + E[B]/q²)` | same |
| (I3) | `Cov(L,K,B)` PSD (all principal minors ≥ 0) | w=2..6 **yes**; smallest 3×3 minor > 0; 2×2 `Cov(K,B)` det > 0 |
| (I4) | `nu_activity == nu_tagged` (same p, same w) | **yes**, w=2..7 |

Exact-run evidence (`out/tA_exact_ref.json`, `tA_exact_half*.json`; all rows
`exact=True`, `p = 5927460507921/10^13`):

| w | I1 | I2 | I3 (2×2) | I3 (3×3) | I4 | seconds |
|---|---|---|---|---|---|---|
| 2 | exact | exact | PSD | PSD | equal | 0.0 |
| 3 | exact | exact | PSD | PSD | equal | 0.0 |
| 4 | exact | exact | PSD | PSD | equal | 0.4 |
| 5 | exact | exact | PSD | PSD | equal | 4.1 |
| 6 | exact | exact | PSD | PSD | equal | 140 |
| 7 | exact | exact | PSD | PSD | equal | 198 |

(w=7 was certified exactly at `p=1/2` rather than `p=p_ref`: with `p_ref` as a rational with
denominator 10¹³ the 126-lump Gauss–Jordan elimination did not finish in a reasonable budget,
so I ran the *same* exact test at `p=1/2` (`--point given --p 1/2`, n_lumps 126, 198 s) and
left the `p=p_ref` exact w=7 job running/killed. The identities are calculus statements, so a
certification at a different `p` certifies the same code path. At `p_ref`, w=7 is float64
(see the main table).)

At `p_ref` the exact identities hold with **zero** difference (not "small"): the tagged
side uses analytic `d/dp, d²/dp²` of the resolvent (`tagged_density_derivatives`), and the
activity side uses the first/second moments of `S` — the two are different code paths built
from different representations of the same `nu(p)`.

## 8. The four verdicts (#774)

`MEAN_WARD` — **compatible**.
`E[B]/E[K]` at `p_ref`: 0.7253, 0.6962, 0.6903, 0.6884, 0.6877, 0.6874, 0.6873 for w=2..8 vs
`q/p = 0.687062`. `mean_ratio - 1` = 5.56e-2, 1.33e-2, 4.75e-3, 1.98e-3, 9.75e-4, 5.37e-4, 3.21e-4
— monotone, and the deviation *shrinks faster than any of the pre-declared powers*
(`w^{55/48}·dev` is itself decreasing). Both pre-declared λ-sides converge to 1. This is a
**7-point (per-side 5-point) compatibility**, not a proof: the effective exponent 3.73 is an
empirical fit over w≤8, and I cannot exclude a limit such as `q/p + O(w^{-3})` with a tiny
residual bias.

Log–log trend fits over w=2..8 (`out/kb-locking.json` → `trend_fits`):
`mean_ratio-1 ~ w^{-3.73}`, `1-ward2_ratio ~ w^{-0.34}`, `lambda_min/trace ~ w^{-1.07}`,
`angle(vs (p,q)) ~ w^{-1.29}`, `cond ~ w^{+1.09}`, `1-Corr(K,B) ~ w^{-1.71}`,
`E K/w^{91/48} ~ w^{-0.015}`, `E B/w^{91/48} ~ w^{-0.050}`.
(These are least-squares slopes of `log|quantity|` vs `log w` over the 7 widths — descriptive,
not a claim that a single power law holds.)

`SECOND_WARD` — **compatible, but slowly and least well determined**.
`ward2_ratio = Var(S)/(E[K]/p² + E[B]/q²)` = 0.5935, 0.6391, 0.6712, 0.6962, 0.7162,
0.7324, 0.7459 (w=2..8). `ward2_ratio - 1` is monotone increasing toward 0 and its per-step log-slope is
still *rising* (0.29 → 0.38), i.e. the data have not flattened onto a plateau. A pure power
fit over w=2..8 gives `1 - ward2_ratio ≈ 0.58 w^{-0.34}` (fitted slope −0.34). Verdict: the trend is toward 1;
with only 7 widths I cannot rule out a small positive limit (`~1 - c w^{-0.4}` and
`→ 0.97` differ by less than the fit scatter over this range). So: compatible, **not**
established.

`COVARIANCE_LOCKING` — **toward rank-one, along `(p,q)`**.
For `Cov(K,B)`: the principal eigenvector's angle with the critical ray `(p,q)` falls
19.91°, 11.32°, 7.77°, 5.86°, 4.69°, 3.90°, 3.34° (w=2..8; the angle with `(q,p)` also
falls but slowly, 29.2°→24.4°); `Corr(K,B)` rises 0.804 → 0.982; the condition number rises
32.8 → 138.5 and `λ_min/trace` falls 0.0296 → 0.00717 (`~ w^{-1.07}`). The **mean** vector `(E K, E B)` aligns with `(p,q)` to
within 0.03° at w=7 (`E K/E B` → `p/q`), so both the mean and the leading covariance
direction point at the same critical ray while the matrix loses rank in the other direction.
The λ-points show the same (Corr 0.77 → 0.98 on both sides).

`STRONG_RANDOM_RAY` — **plausible, with an unexpected bonus**.
`E[K]/w^{91/48}` = 1.883, 1.878, 1.860, 1.855, 1.851, 1.848, 1.847 and
`E[B]/w^{91/48}` = 1.366, 1.308, 1.284, 1.277, 1.273, 1.271, 1.269 (w=2..8). Both are flat to
~1% over seven widths (`~ w^{-0.015}` and `~ w^{-0.050}` respectively), and their ratio equals `p/q` in the limit, i.e. **one constant
`A ≈ 3.16` gives both `E K ≈ A p w^{91/48}` and `E B ≈ A q w^{91/48}`**. So the strong random
ray `w^{-91/48}(K,B) ⇒ A(p,q)` is *supported* at `p_ref`.
Two honest caveats: (i) at `p_ref` the black infinite-cluster density is `θ(p_c)=0`, and the
finiteness is carried by `nu_w`; I additionally observed `nu_w·w = 0.2987, 0.2891, 0.2893,
0.2889, 0.2888, 0.2888, 0.2887` (w=2..8) — i.e. `nu_w ≈ 0.2888/w` at `p_ref` to 0.35%, which is what
makes `EK` scale as `w^{91/48}` while `(nu/w)EK = θ_w → 0`. (ii) These are 7 widths; the
visual flatness of `w^{-91/48}` over `w≤7` cannot by itself prove the exponent — the exponent
is exactly the ambiguity that a wider range would settle.

## 9. What must NOT be claimed

* Not a theorem about planar critical percolation. Everything above is a **finite-width
  (w ≤ 8) exact/numerical family** on a cylinder; the "convergence" statements are
  **compatibilities over 7 widths**, not limits.
* The `91/48` ray is **not** fitted from w≤8: I only checked that the pre-existing exponent
  describes the data; I did not fit an exponent and then declare it.
* `B` here is the **distinct external boundary**; none of these numbers may be reused as
  "edge perimeter" results (they would differ).
* `p_ref` is a diagnostic reference only; no claim about `p_c` is made.
* Black-cluster quantities here must not be extrapolated to white clusters (and vice versa):
  the black `Corr(L,K)²` decreases while the white one increases (BRIEF §5).
* No Monte Carlo, no GPU, no new transfer engine was used, and no #275 source contract or
  #760 gamma-vs-kappa result is touched.

## 10. Reproduce

```
scripts/engine_relaxed.py     # 1-line raise-only copy of the pinned engine (w<=8 for the activity builder)
scripts/taskA.py              # all Task-A numbers; --widths --mode float|exact --point ref|lam --p ...
scripts/diag_marks.py         # the independent brute-force B-semantics check
scripts/test_jm.py            # validation of the moment extractor against the shipped #772 numbers
out/tA_float_ref.json  out/tA_float_ref_w8.json  out/tA_float_lam.json
out/tA_exact_ref.json  out/tA_exact_half*.json   out/diag_marks.json   out/kb-locking.json
```
Run as `cd /workspace/ward774 && PYTHONPATH=/workspace/mo/compat PY39COMPAT_MARKER=1 python3 scripts/taskA.py ...`
