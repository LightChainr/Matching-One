# #586 — Graded Q→1 tangent benchmark after #581

**Branch:** `theory/p586-graded-q-tangent-20260913`
**Issue:** #586 (graded tangent benchmark, B_even vs ambient-homology)
**Status:** Phase A (continuum internal check + exact lattice graded table) delivered and exact. Phase B bulk-logarithmic numeric declared a buy-back. Spin-4 reopening not yet earned.

## Method

Reuses the exact finite-volume gate of #581 (`scripts/qtangent_scale_decomposition`:
`enumerate_torus`, `covariance_split`, `exact_identity_gate`) for the **lattice
control**, and adds the **Phase A continuum internal check**: differentiate the
BPZ/(2,1)-degenerate differential equation in the coupling and verify the
Q-tangent satisfies the resulting inhomogeneous linear ODE. Deliverable is a
**graded tangent table**, not a Jordan yes/no verdict. The two measure-score
contributions (`B_even` duality-even Betti, `X = r−1` ambient homology) are
reported **separately** and only summed afterwards.

## Phase A — continuum internal check

Representative (2,1)-degenerate BPZ hypergeometric ODE
`L(h) G = λ(1−λ)G'' + (1−2λ)G' − h(1−h)G = 0`, with the Q/coupling entering
through `h`. Differentiating: `L(h) T = (1−2h) G` with `T = ∂_h G`. The finite-
differenced tangent `T` is solved on a 4001-point grid and the residual
`R = L(h0) T − (1−2h0) G(h0)` evaluated: **max|R| / ||G|| = 4.4e-6** (threshold
1e-3) → internal continuum check **passes**.

> Honesty note: this uses a representative (2,1)-degenerate BPZ ODE to demonstrate
> the *differentiate-the-equation* method and its inhomogeneous-ODE verification.
> The exact closed-form Cai arXiv:2603.28161 boundary four-point is **not**
> reproduced here; it is a declared buy-back. The check is the "internal continuum
> check independent of lattice sampling" the issue asks for, not a claim about the
> precise Cai amplitude.

## Phase A/B — exact lattice graded table (reused #581 gate, L=2,3)

`Cov(O, T) = ½ Cov(O, B_even) + ½ Cov(O, X)`, configuration by configuration,
at `p = 1/2` (every config equiprobable). Identity gate `T − T* = X` has **0
failures** at both L=2 (256 configs) and L=3 (262144 configs). The split is
**exact** (`split_is_exact = True`) for all six observables at both sizes.

Highlights (L=2; identical structure at L=3):

| observable | Cov(O,T) | B_even piece | X piece | topo frac |
|---|---|---|---|---|
| open_edges | 27/64 | **0** | 27/64 | 1.0 |
| wrap_either | 8405/65536 | −427/65536 | 69/512 | 1.051 |
| wrap_cross | 9259/65536 | **+427/65536** | 69/512 | 0.954 |
| components | 2911/65536 | 7903/65536 | −39/512 | −1.715 |
| cycle_rank | 30559/65536 | 7903/65536 | 177/512 | 0.741 |

Two structural facts the table establishes:

1. **`open_edges` loads purely on `X`** (B_even = 0): the ambient-homology source
   is a real, distinct channel, not a rename of B_even.
2. **`wrap_either` and `wrap_cross` carry opposite-sign B_even pieces**
   (∓427/65536) with identical X — the duality-even piece flips while ambient
   homology is shared. This is exactly #581's claim that B_even and X are **not**
   two algebraic names for the same direction: they separate under a concrete
   lattice observable.

The two measure-score contributions are therefore reported separately and only
combined after reporting, per the issue's hard rule.

## Phase B — known bulk logarithmic pair

The precise generic-Q bulk-logarithmic control (Vasseur–Jacobsen–Saleur /
Camia–Feng) is **not in the tree**. It is declared a **buy-back**; the graded
*method* and the exact L=2,3 lattice graded table are delivered. No bulk-log
signature is asserted without that data.

## Decision (per issue decision table)

- Boundary exact tangent (lattice control): **passes** (0 identity failures,
  exact split for all observables/sizes).
- Internal continuum ODE check: **passes** (residual 4.4e-6).
- Known bulk-log grading: **declared buy-back** (data absent).
- Therefore: promote the graded tangent machinery as an **exact method-level
  diagnostic** for #581's scale tomography. The narrow spin-4 tangent reopening
  of #263 is **not yet earned** — it remains conditional on the Phase B bulk-log
  positive control, which is a forward buy-back. No result is silently transported
  to square site; `B_even`/`X` are kept as typed lattice Q-lift objects.

## Claim boundary (preserved)

`B_even` is a duality-even measure tangent, not automatically the energy field;
`X` is an ambient-topology score, not automatically a defect insertion; a
successful tangent benchmark calibrates a method, it does not identify square-site
Matching One; the two graded terms are correlated parts of one derivative, not
independent evidence.
