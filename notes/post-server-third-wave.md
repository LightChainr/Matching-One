# Post-server third-wave interpretation

Status: updated after Huawei server checkpoint `28fd1d872f97a83c6ddd59bde996dae2343d3c77` on PR #21.

This note records what is now resolved, what was falsified as an execution gate, and the smallest set of follow-up experiments capable of distinguishing an angular coincidence from a genuine matching-odd spin-4 finite-size sector.

## 1. Current evidence state

### Small-Pell calibration is positive

For `(a,d)=(7,5)`, two independent 5,000,000-replica runs gave

```text
+4.6267e-4 +/- 5.06e-5
+4.0417e-4 +/- 5.10e-5
```

for the diamond-minus-axis matching-root gap. Their inverse-variance pool is

```text
+4.3365e-4 +/- 3.59e-5,
```

about twelve standard errors from zero. Halving the finite-difference step changed the result by less than `1e-6`, well below sampling error. Together with the exact `(3,2)` pair, the two-point effective root-gap exponent is about `4.23`; the `(7,5)` gap multiplied by the fourth power of the common physical scale is about `1.06`.

This is a successful power calibration, not an asymptotic exponent fit. Large-Pell fixed-`p` scans remain a poor discovery design.

### Same-N Gaussian orientation signal is resolved

For periods

\[
(a,b),\qquad(-b,a),\qquad N=a^2+b^2,
\]

the 30,000,000-replica production run measured the matching-function orientation difference at one frozen probability:

| `N` | representations | `Delta M` | SE | z | `N^(13/8) Delta M / Delta cos(4 theta)` |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)` vs `(7,4)` | `+1.00377e-3` | `1.68448e-4` | 5.96 | `0.65016 +/- 0.10911` |
| 85 | `(9,2)` vs `(7,6)` | `+7.60333e-4` | `1.58350e-4` | 4.80 | `0.65117 +/- 0.13562` |
| 145 | `(12,1)` vs `(9,8)` | `+3.28667e-4` | `1.95415e-4` | 1.68 | `0.55743 +/- 0.33143` |

All signs agree with `Delta cos(4 theta)`. The first two scaled amplitudes agree much more closely than their Monte Carlo errors require. After dividing out the angular factor, their two-point radial exponent is approximately `1.619`, close to `13/8=1.625`.

Checkpoint `28fd1d8` adds a disjoint-counter 2,000,000-replica evaluation at `N=65,85`:

```text
N=65: +1.3560e-3 +/- 6.20e-4
N=85: +4.3950e-4 +/- 6.46e-4
```

These are low-power but statistically compatible with the 30-million run. A descriptive inverse-variance pool gives approximately

```text
N=65: +1.0280e-3 +/- 1.626e-4
N=85: +7.4215e-4 +/- 1.538e-4
```

The two campaigns use disjoint counter ranges and different production provenance, but one was generated from source present in a working tree before the source commit. A clean-checkout, new-seed confirmation remains mandatory.

The scientifically defensible statement is now:

> A same-area, same-shape orientation signal is resolved in the square-site matching function at `N=65` and `N=85`, with the sign and first two normalized amplitudes predicted by a leading `cos(4 theta) N^-13/8` term.

This promotes the hypothesis to confirmation-stage P0. It does not yet determine the asymptotic exponent, exclude logarithmic mixing, or identify a CFT operator.

## 2. Sharpened working law

For fixed square-torus modulus,

\[
M_N(p_c;\theta)
=N^{-13/8}\left[A_0+A_4\cos(4\theta)+A_8\cos(8\theta)+\cdots\right]
+\text{subleading terms}.
\]

At identical `N`, the scalar part cancels in an orientation difference. The primary frozen model is

\[
\boxed{
\Delta M_N
=A_4\,\Delta\cos(4\theta)\,N^{-13/8}
\left[1+bN^{-\omega/2}+\cdots\right].
}
\]

The logarithmic alternative remains live:

\[
\Delta M_N
=\Delta\cos(4\theta)\,N^{-13/8}
\left[A_4+B_4\log N+\cdots\right].
\]

Because `M'_N(p_c)~N^(3/8)`, the associated root difference should scale as

\[
\Delta p_N^*\sim N^{-2}=L^{-4}.
\]

A high-value operator candidate is a matching-odd spin-4 field of total dimension `x=21/4`, possibly a level-four thermal-family descendant or quasiprimary in the `c=0` logarithmic theory. This is a conjecture, not an identification.

## 3. Gates changed by the checkpoint

### Matching-even dominance is rejected as a prerequisite

The matching-even orientation sector is channel-dependent and does not exhibit a stable larger harmonic. It may change sign between wrapping definitions. The earlier plan to require a large matching-even `L^-2` signal before accepting the matching-function result is therefore rejected as an execution gate.

The appropriate task is a joint decomposition by matching parity, rotation harmonic, thermal parity, and radial exponent, without assigning the dominant sector in advance.

### Wrapping-only GLS is structurally impossible

For the validated complementary configurations, the `cross`, `both`, `either`, and directional matching differences are identical configuration by configuration. The wrapping-only covariance matrix is rank one; optimized and equal-weight combinations have exactly the same variance.

Variance reduction is moved to issue #34: centered Euler/local-motif controls built around

\[
D_N=q+(V-Np)-(E-2Np^2)+(F_0-Np^4).
\]

### High-precision final-tail defect is closed

Checkpoint `28fd1d8` sets `mp.mp.dps` before parsing decimal CSV strings, adds a long-decimal regression, regenerates the summary, and records the corrected intercept

```text
0.59274605094603206266439366806726549
```

with held-out RMSE `1.3775861250986e-11`. Issue #32 is completed.

## 4. Threshold-rank representation is validated

C05 freezes `K_minus/K_plus` conventions and reconstructs `M`, `M'`, and the root from integer histograms. An axis `L=8`, 100,000-permutation pilot gives

```text
mean K_minus       35.59763
mean K_plus        41.10408
mean rank gap       5.50645
root                0.5925842499338915123
M'(p_ref)           8.33469658658750
```

with 83.85 seconds wall time in pure Python. The representation passes; the implementation does not pass the production-throughput gate. The next core engineering task is a deterministic C++/OpenMP port with bit-for-bit convention regressions.

Threshold ranks are now required for the orientation project. They allow the angular difference to be evaluated across the entire thermal window and distinguish a central root-shifting amplitude from an orientation-dependent slope times `p_ref-p_c`.

## 5. Kappa3 control status and modulus caveat

Square-bond and triangular-site exact-threshold score pipelines both pass exact regressions. The triangular production sequence on its natural 60-degree rhombic torus does not identify the correction exponent; the free-power scan chooses its lower search boundary.

The triangular torus has `tau=exp(i*pi/3)`, whereas the square-bond sequence has `tau=i`. Since `kappa3` is shape-dependent, a shared-intercept fit across those sequences is not a fixed-shape universality test. It is only a method diagnostic.

The next same-shape control should be union-jack/self-matching site percolation on a square torus. Triangular data should be extrapolated as their own fixed-modulus sequence, or compared to square-bond data implemented at the same modulus.

## 6. Finite-width extrapolation status

Leakage-safe Stage A selects degree-four polynomial `F(n^-2)` with `n_min=9`. Its held-out errors at widths 19--21 are

```text
+2.0916e-12, +4.8114e-12, +8.3517e-12
```

with RMSE `5.6943e-12`. This improves the earlier baseline but preserves the positive increasing drift. The best admissible Padé `[2/2]` is worse. Within the declared families, rational corrections do not cure the drift.

Stage B remains frozen until provenance-complete widths 22--24 are imported.

## 7. Highest-information next experiments

### A. Clean-checkout Gaussian confirmation

Use a new seed/counter range and committed source for

```text
N=65  (8,1)/(7,4)
N=85  (9,2)/(7,6)
N=130 (11,3)/(9,7)
N=145 (12,1)/(9,8)
N=170 (13,1)/(11,7)
```

Use at least 100 equal batches and power sample counts from the observed covariance. Do not select/drop sizes by significance. `N=130` and `N=170` provide better angular leverage than the older `N=205,425` confirmation order.

### B. C++ threshold-rank production

Port the frozen Python contract, validate exact histograms and counter concatenation, then run paired Gaussian orientations. Reconstruct `Delta M(p)`, derivatives, and root gaps after the run.

### C. Blind radial/thermal model test

Train on `N=65,85,130`; hold `N=145,170`. Compare:

1. fixed `13/8` constant amplitude;
2. fixed `13/8` plus one declared power correction;
3. fixed `13/8` plus logarithmic amplitude;
4. a free exponent selected on training sizes only.

Use threshold-rank curves to form a central reflection from levels of the orientation-averaged matching function rather than from a favored decimal `p_c`.

### D. Euler/motif controls

Implement exact microcanonical centering and pilot-frozen/cross-fitted weights. Duplicate wrapping indicators are rejected as controls.

### E. Same-modulus kappa3 control

Implement union-jack square-torus data before any cross-model universal/rational claim.

## 8. Quotient caveat for the magic orientation

The cyclic fast path requires `gcd(a,b)=1`. At `N=169`, `(13,0)` has quotient group `Z_13 x Z_13`, whereas `(12,5)` has `Z_169`. Both are exact square tori, but they do not possess a canonical shared cyclic labeling.

The pair remains useful as a held-out near-`pi/8` geometry test through the general period-matrix engine, provided Smith normal form and coupling conventions are reported explicitly.

## 9. Decision

The core hypothesis is now

\[
\boxed{
\Delta M_N\propto\Delta\cos(4\theta)N^{-13/8}
}
\]

with independent-seed, multi-size, full-thermal-curve, and eventually multi-angle confirmation required.

A successful confirmation would explain the empirical `L^-4` matching-root bias through a rotation-odd matching sector. The current evidence is strong enough to justify a production confirmation campaign, but not an exact exponent or operator claim.
