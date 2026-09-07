# Probe MID — Does M(1/2) pin the shape on exact tori?

**Issue #625 · D1–D5 · 2026-09-07 · exact rational arithmetic throughout; no Monte Carlo; no L=4 bond; no arXiv; no N=725.**

Kill/promote probe for #622's **W4**: the odd moment `M(1/2)` either governs
the shape `Z_L` of the threshold law, or it does not. Census numbers
(`M(1/2) = -21/64` etc.) are **inputs** here, imported from #606/#627 (PR for
#619), not re-verified deliverables.

## Verdict (D5)

```text
W4 at L=3,4:  KILLED
evidence:     Game B on L=3: along the m-pinned line F_lambda(p)=1/2
              (m = E_lambda[r_b] - 1 = 0 EXACTLY for every lambda),
              Z_lambda(0.5; a=0.1, b=0.9) moves
              0.73817 -> 0.72795 -> 0.69151 -> 0.62685 -> 0.59329 ->
              0.54242 (lambda=1) -> 0.49120 -> 0.45693 -> 0.38844 ->
              0.34430 -> 0.32192 (lambda = 1/16 ... 16);
              spread over the 11-point lambda grid at interior u:
              0.416 (a=0.1/0.9) and 0.271 (a=0.2/0.8), both >> 1e-8.
              L=4 spot check (lambda = 1, 5/2, 2/5):
              Z(0.5; 0.1/0.9) = 0.53359 / 0.44405 / 0.62239,
              spread 0.178 >> 1e-8.
Z_3 vs Z_4:   ||Z_3 - Z_4||_inf on {0.1..0.9} = 0.0161 (anchors 0.2/0.8)
bond vs site: INCOMPARABLE at L=3 — bond X = r-1 law at p=1/2 is the
              3-atom measure (0.28786, 0.42429, 0.28786); inverse-CDF is a
              staircase; Z_bond undefined by the probe's own rule.
```

W4 is dead on exact L=3,4: `Z` is **not** a function of `m = E[r_b] - 1`.
What survives for #622 to absorb as a finite-L fact: the whole Game B
one-parameter family is a physical-family counterexample family, and the
two-size W1 datum `||Z_3-Z_4||_inf = 0.0161` sits between the probe's
thresholds: constancy is not killed at O(1e-3), not dead at O(1e-1) either.

## Programs

- **D1** config tables streamed once (`stream_joint`), joint
  `(n_black, r_b)` accumulated; Alexander `r_b+r_w=2` asserted on the stream
  (precondition, dual_fail = 0 at both sizes).
- **D2** physical `Z_L(u; 0.2/0.8)` exact at L=3,4:

| u | Z_3 | Z_4 |
|---|---|---|
| 0.1 | −0.270776433 | −0.275428409 |
| 0.2 | 0 (anchor) | 0 (anchor) |
| 0.3 | 0.201001621 | 0.199244567 |
| 0.4 | 0.372128340 | 0.367673381 |
| 0.5 | 0.528759509 | 0.522146446 |
| 0.6 | 0.679946684 | 0.672548493 |
| 0.7 | 0.833520795 | 0.827685264 |
| 0.8 | 1 (anchor) | 1 (anchor) |
| 0.9 | 1.203239386 | 1.219298955 |

  `||Z_3-Z_4||_inf = 0.01606` (worst u = 0.9).  p-choice reading (issue D2):
  "Z at p=p_L^H" = Z at u=1/2 (0.5288 vs 0.5221, |Δ| = 0.0066); "Z at p=1/2"
  = Z at u = F(1/2) (0.2649 vs 0.1813, |Δ| = 0.0837 — the same-p row is far
  apart because M(1/2) ≠ 0 moves the u where Q hits 1/2).  Self-symmetry
  rows: `Q(u)+Q(1-u)-1` max 0.1712 (L=3) / 0.1803 (L=4), plus the two
  p-choice rows `2p_L^H-1` = 0.1730 / 0.1813 and
  `Q(F(1/2))+Q(1-F(1/2))-1` = 0.1681 / 0.1762.
- **D3 Game A** (L=3 fake law, m pinned at −21/64 exactly): tilting
  INSIDE `{r_b=1}` by `s^{n-n_min}` keeps rank masses — hence m —
  invariant, yet the fake law `G = n_black/N` has quantiles whose shape
  swings Z(0.5; 0.1/0.9) from 0.333 (s=1) to 0.667 (s=2) to 1.0 (s=4);
  spread at fixed m = 0.667.  Even the fake-law version of "m pins Z" fails.
- **D3 Game B** (the game that matters): P(σ) ∝ p^n (1-p)^{N-n} λ^{r_b},
  λ = e^β > 0 rational (exact weights; β = ln λ reported as float only).
  On the line `F_λ(p) = 1/2` the odd moment is pinned at m = 0 exactly for
  every λ; the shape still moves (numbers above).  **W4 dies for the
  physical interpolation.**  Monotonicity of F_λ asserted on a 40-point
  p-grid for every λ before inversion.
- **D4** square-bond L=3, all 2^18 = 262144 configs, 3.3 s (time box 20
  CPU-min).  Rank counts [75460, 111224, 75460] match #627 C6 exactly (same
  `r` sanity, imported not re-derived).  `M_bond(1/2) = 0` EXACTLY and
  `P(r=0) = P(r=2)` EXACTLY: duality-oddness of X = r−1 is checked, not
  assumed.  `p_wrap(1/2-curve) = 0.419649` (solve P_p(r≥1) = 1/2 — bond
  curve, not a site claim).  The X-law at p=1/2 is 3-atom ⇒ `Z_bond`
  undefined ⇒ the bond/site comparison is INCOMPARABLE, with the jump
  structure recorded in JSON.

## Tolerance statement

All quantities are exact rationals; float displays are roundings.
W4's kill margin (ΔZ ≈ 0.4 at fixed m = 0) is ~4e7 times the 1e-8
tolerance of the issue; bisection runs 110 exact-Fraction iterations
(~2^-110 ≈ 1e-33), far below the required ≤1e-14.

## Inputs (cited, not re-verified)

`M_L(1/2) = -21/64` (L=3), `-13757/32768` (L=4); rank-pair census;
`p_L^H` = 0.586511455113 / 0.590672112331 — all from PR #606 /
PR #627 (issue #619 artifacts: M polynomials, Q(u) tables, joints).
D2's tables reproduce the #619 C3 quantile columns to 1e-9 at every
u ∈ {0.1..0.9} (that is the import-consistency check, not the job).

## Interface

#622 W4: **KILLED** by Game B (absorb the finite-L fact; the shape is not a
function of the odd moment on exact L=3,4).  #619 polynomials: imported as a
library (quantile columns match).  #608: `r` is the bond observable (rank
convention matched to #627 C6, counts identical).  #618/#620: not involved.
No productions scored.  Nothing enters `docs/STATUS.md`; no ticket closed.

## Scope

No Monte Carlo; no L=4 bond enumeration; no N=725; no exponent fits
(`Z_3, Z_4` not fitted to L^{-θ}); no STATUS edit.  Stop rules: none fired
(Game B coded without any new percolation engine; D1 precondition held;
D4 ran inside budget, so it was run even though Game B had already killed).
