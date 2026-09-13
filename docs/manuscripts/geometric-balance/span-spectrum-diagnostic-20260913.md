# Span-spectrum diagnostic for the heat-kernel sewing conjecture (§6)

> **Partly superseded — see `span-spectrum-erratum-20260913.md`.** The solved width range is
> 2–7, not 2–8; the §2 tail-negligibility sentence is false for 4 of 30 configurations; the
> §4/§5 "additive span" mechanism reading is withdrawn; and two solver defects (an unnormalised
> `light` chain, an implicit certificate-versus-observable gap) are fixed there. The §3/§4
> tables themselves stand and are reproducible from the committed results JSON.

**Status: finite-width DIAGNOSTIC of the §6 research conjecture, delivered as the ticket asked.
The conjecture is neither proved nor refuted; its moment prediction is not supported at any
accessible fixed-p width, and its scaling exponent is near-diffusive but not yet settled.**

This note answers the testability demand of §6.3 directly: instead of fitting another exponent,
we measure the vertical-span distribution of the complete winding cluster — the object the
conjecture makes a sharp prediction about — exactly, at fixed subcritical p, for both
adjacencies, widths 2-8.

## 1. The measurable object

Let d_h(p) = nu_(w, span=h) be the anchored density of complete winding clusters whose vertical
span is exactly h (§3 of the sewing-with-memory note). The one-frontier retirement engine
(`scripts/cylinder_winding_intensity.py`, validated against all 18 published controls) is extended
with a per-component min-row offset: each retiring winding component reports its span, giving the
full spectrum d_1..d_D_MAX plus the exact tail mass (span > D_MAX).

Two exact facts hold at every width, and they are the engine's validation:

1. **Closure.** The depth-clamped chain is an exact lumping, so
   `sum_h d_h + tail = nu_w` EXACTLY, with nu_w the independently certified density of the
   winding-intensity engine. At NN w=4, p=1/4 the closure reproduces the #741 certified rational
   `1750262847.../(4^32)`-scale fraction to the last digit.
2. **Controls.** d_1 = p^w (1-p)^(2w) exactly (the full-row cluster); at NN w=4, p=1/2,
   `sum_{h<=3} d_h = 9087/1048576` — the §4.1 finite-height control of the sewing-with-memory
   note, reproduced digit for digit.

## 2. Data

Both adjacencies, p in {1/8, 1/4} (NN) and {1/16, 1/8} (matching), plus p=1/2 for w<=4;
widths w = 2..8 (D_MAX = 48 for w<=6, 20 for w=7, 16..20 for w=8). Tail fractions are
negligible everywhere (<= 2e-6, mostly < 1e-20). Machine-readable values are in
`results/geometric-consistency/span-spectrum-20260913.json`.

## 3. Finding 1 — the mean span is sublinear, near-diffusive but not settled

E[L] grows sublinearly at every fixed subcritical p:

| family | E[L], w=2..7 | E[L]/w, w=2 → 7 |
|---|---|---|
| NN p=1/4 | 2.08, 2.64, 3.14, 3.57, 3.94, (w7: 4.27) | 1.04 → 0.61 |
| NN p=1/8 | 1.52, 1.77, 2.01, 2.24, 2.45, (2.65) | 0.76 → 0.38 |
| matching p=1/8 | 2.05, 2.70, 3.24, 3.68, 4.06, 4.40 | 1.03 → 0.63 |
| matching p=1/16 | 1.85, 2.23, 2.63, 2.94, 3.22, 3.47 | 0.92 → 0.50 |
| NN p=1/2 (w<=4) | 3.62, 5.42, 7.09 | 1.81 → 1.77 (FLAT) |

Local log-log exponents at NN p=1/4: alpha ≈ 0.55-0.61 over w = 2..6 — near-diffusive,
drifting down, but the width range cannot yet separate alpha = 1/2 from alpha = 3/4.
At p = 1/2 the exponent is 1.0 (E[L] proportional to w): the dense regime is NOT diffusive.

## 4. Finding 2 — the moment prediction of §6 is not supported

The §6.3 conjecture predicts, once the sewing identifications hold,

    Var(L)/(E L)^2  ->  pi/3 - 1 = 0.0471975511966...

(Brownian-bridge range law). The measured ratios decline steadily and are 2-7x the target at
the largest widths:

| family | Var/(E L)^2, w=2 -> 7 |
|---|---|
| NN p=1/4 | 0.275, 0.250, 0.211, 0.180, 0.157, (0.140) |
| NN p=1/8 | 0.229, 0.232, 0.214, 0.192, 0.171, (0.152) |
| matching p=1/8 | 0.184, 0.200, 0.178, 0.157, 0.140, 0.127 |
| matching p=1/16 | 0.133, 0.122, 0.117, 0.105, 0.096, 0.089 |
| NN p=1/2 (w<=4) | 0.358, 0.333, 0.304 |

The decline is consistent with the CLUSTER-CHAIN picture: at fixed p the complete winding
cluster is a one-dimensional chain of irreducible pieces (CIV skeleton); its span is the sum of
~w/mu piece heights, so Var(L) grows like w while (E L)^2 grows faster, driving the ratio toward
zero — NOT toward the Brownian-bridge constant. No family shows the plateau at pi/3 - 1.

## 5. What this does to the conjecture

The §6.2 falsification list is exercised for the first time, and the data land on its own
named failure modes: the complete cluster carries MORE transverse fluctuation than a single
diffusing skeleton (bushes/branch Deaths add vertical extent piecewise), so the heat-kernel
trace over spans describes at most the SKELETON observable, not the complete-component span,
at fixed p. Concretely:

- The pure single-diffusion form `nu_(w,<=H)/nu_w -> F_range(H/sqrt(D_G w))` is not supported
  at any fixed p measured: neither the H/sqrt(w) collapse nor the moment constant appears.
- In the DILUTE joint limit (w p^2 -> 0) the proved Bessel regime already covers the
  near-straight rows; the conjecture's remaining content at fixed p would have to be a
  SKELETON-only statement (condition on no bushes / restrict to the irreducible chain), or
  absorb the piecewise-additive fluctuation into zeta_G(p).
- The measured slopes s(p) = d(E[L]-1)/d(w) — 0.24 (NN p=1/8), 0.50 (p=1/4), 1.75 (p=1/2);
  0.41 (matching p=1/8) — are themselves new exactly-measurable micro-objects feeding the
  amplitude question of #740.

## 6. Method and cost

- Engine: `scripts/span_spectrum_build.cpp` — the validated one-frontier advance() with an
  integer min-row offset per component (merge rule max, retirement span = offset), depth
  clamped at D_MAX; the clamped chain is an exact lumping (bins 1..D_MAX exact, tail exact).
- Solver: `scripts/span_spectrum_solve.py` — chunked streaming construction of the sparse
  chain (no 15 GB table ever resident), stationary solve by sparse LU (n <= 5e4) or power
  iteration, exact int64 residual certificate where the scale fits, exact-rational solve for
  n <= 300.
- Cost, measured: w=8 build 469 s (2 505 625 frontier states, 604M transitions, 15.4 GB table,
  ARM container); solves w<=6 minutes each on one 16-vCPU container; the full grid
  (38 (graph,width,p) configurations) ran across 10 containers in parallel.
- Every number above is a finite-width diagnostic; nothing here proves or refutes the
  conjecture asymptotically, and no new p_c, Monte Carlo, or GPU was used.
