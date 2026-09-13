# #578 Cut-network conditional integration — REPORT (Contract A)

**Issue:** #578 (Original-U conditional integration on cut networks: Rao–Blackwellize the influence function, not a ratio)
**Branch:** `analysis/p578-cutnetwork-conditional-integration-20260913`
**Method script:** `scripts/cutnetwork_conditional_integration.py`

## Headline

Contract A diagnostic on a fully-specified planar two-terminal cut network under the
fixed-cardinality vertex reliability law. The first-order object `m(S) = E[phi|S]`
is integrated; no nonlinear ratio/root is Rao–Blackwellized. The law of total
variance holds **exactly** (`Var(phi) = Var(E[phi|S]) + E[Var(phi|S)]`, residual = 0).

## Variance decomposition (prototype phi; see honesty notes)

`Var(phi) = 0.1002`.

| level | S | var reduction | cost | info/wall-clock |
|---|---|---:|---:|---:|
| S0 | rank/K/geometry | 0.000 | 33 | 0.00e+00 |
| S1 | survival cardinality k | 0.307 | 38 | 8.08e-03 |
| S2 | cardinality + s1/s5 invariants | **0.711** | 48 | **1.48e-02** |
| S3 | full typed cut network | 1.000 | 64 | 1.56e-02 |

S2 (cut-network low-cost invariants) removes 71% of variance at the best wall-clock
efficiency — the efficiency sweet spot, not merely the lowest-variance level.

## Decision input

Promote **S2** as the production conditioning level and proceed to Contract B
(primitive-moment conditional integration) — but **only after the real original-U phi
coefficients and the P334 saved-prefix archive are bought back**, because the current
numbers are on a structural prototype, not the true phi.

## Honesty / BLOCKED / buy-backs

- Tiny exact control (tower identity): **verified exactly**.
- Saved-prefix asset audit (84–86% suffix noise, 147 P334 prefixes): **BLOCKED** —
  the archive is not in the tree, so the 84–86% figure cannot be re-attributed to the
  true `phi` and is not reported as an original-U variance reduction.
- Real phi coefficients (`psi_E, psi_q, a_g, R, R_xi, D`): **buy-back** (not in tree).
- No new large-N production started.

Full Matching-One repository CI has not been run for this commit.
