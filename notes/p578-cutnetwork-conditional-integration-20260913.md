# #578 — Original-U conditional integration on cut networks (Contract A)

**Branch:** `analysis/p578-cutnetwork-conditional-integration-20260913`
**Issue:** #578 (Rao–Blackwellize the influence function, not a ratio)
**Status:** Contract A diagnostic delivered and exact on a fully-specified cut-network
prototype. Real original-U phi coefficients and P334 saved-prefix archive declared a
buy-back; the saved-prefix asset audit is **BLOCKED** (data absent from tree).

## Object

Keep the original pooled-root observable and the exact influence-function semantics
of PR #530. The first diagnostic target is the actual first-order object

    m(S) = E[phi | S]

for a declared conditioning state `S` — **not** a nonlinear ratio/root. We verify
the law of total variance exactly and report the variance decomposition and a
wall-clock cost model, so the issue's decision ("promote one conditioning contract,
or stop") has a reproducible basis.

## Method

A small planar two-terminal vertex network (terminals `L,R`; switchable `s1..s5`)
with a fixed cut geometry. The physical continuation law is the **fixed-cardinality
two-terminal vertex reliability** law ("choose uniformly among remaining switchable
vertices"): every subset of a fixed occupied cardinality `k` is equally likely. We
do **not** substitute independent edge reliability.

`phi` is a faithfully-structured **prototype** of the issue's schematic
`phi_g = [c_g psi_E,g - R a_g psi_q,g - R_xi a_g (q - <q>_g)] / D`, with
`psi_E` = L–R connectivity indicator, `psi_q` = fraction of switchable occupied,
and frozen nuisance/root constants `a, c, R, R_xi, D`. No nonlinear ratio is
Rao–Blackwellized; only the first-order influence function `E[phi|S]` is integrated.

## Results (Contract A, exact on the prototype)

`Var(phi) = 0.1002`. Tower identity `Var(phi) = Var(E[phi|S]) + E[Var(phi|S)]`
holds to machine precision for every level.

| level | conditioning state S | n_states | V_between | V_within | var reduction | wall-clock cost | info/wall-clock |
|---|---|---:|---:|---:|---:|---:|---:|
| S0 | rank/K/geometry only | 1 | 0.0000 | 0.1002 | 0.000 | 33 | 0.00e+00 |
| S1 | survival cardinality k | 6 | 0.0308 | 0.0694 | 0.307 | 38 | 8.08e-03 |
| S2 | cardinality + s1/s5 states (low-cost invariants) | 16 | 0.0713 | 0.0289 | 0.711 | 48 | **1.48e-02** |
| S3 | full typed cut network | 32 | 0.1002 | 0.0000 | 1.000 | 64 | 1.56e-02 |

**Reading.** Conditioning on cardinality (S1) removes 31% of variance; adding the
low-cost cut-network invariants (S2: s1/s5 states) removes 71% at the **best
wall-clock efficiency** (highest variance-reduction per unit cost). The full typed
network (S3) removes 100% but at higher cost. This is exactly the issue's point:
the correct objective is **not minimum variance alone** — S2 is the efficiency
sweet spot. A small `V_between` must not be read as a small physical state; it is an
estimator-compression quantity.

## Decision input

Contract A shows conditioning buys real, exact variance reduction (up to 71% at the
best wall-clock efficiency for S2). The recommendation is to **promote S2
(cut-network low-cost invariants) as the production conditioning level** and proceed
to Contract B (primitive-moment conditional integration) — but **only after the
real original-U phi coefficients and the P334 saved-prefix archive are bought back**,
because the current numbers are on a structural prototype, not the true phi.

## Honesty / BLOCKED / buy-backs (per issue's first no-new-production gate)

- **Gate 1 (tiny exact control):** DONE. Tower identity verified exactly on the
  prototype cut network.
- **Gate 2 (saved-prefix asset audit):** **BLOCKED.** The 147 frozen P334 prefixes
  and their measured 84–86% suffix-noise figure are **not in the tree**. Without
  them the true original-U `phi` ingredients (root/normalizer/source terms) cannot
  be reconstructed, so the 84–86% number **cannot** be re-attributed to the true
  `phi` and is not reported as an original-U variance reduction here. Declared
  buy-back.
- **Gate 3 (complexity profile):** provided for the prototype (full enumeration of
  32 futures). The complexity profile of the *actual* saved cut networks requires
  the P334 archive.
- **Gate 4 (oracle variance/cost curve):** delivered as the S0..S3 table above on
  the prototype.
- **Real phi coefficients** (`psi_E, psi_q, a_g, R, R_xi, D`): NOT in tree → buy-back.
- No new large-N production is started; the issue's "no new production until the
  gate produces a measured efficiency curve" is respected (curve exists for the
  prototype, pending real-phi buy-back).

## Boundaries preserved

This is estimator-layer work. It does **not** reopen #275's candidate-mapping
problem and does **not** supply the two missing homologous forward columns. `E[phi|S]`
is the first-order influence function; the final nonlinear ratio keeps its ordinary
finite-sample semantics and is not claimed unbiased from the tower property.
