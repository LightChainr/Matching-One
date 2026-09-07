# Literature officer — SKZ, Grimmett–Li, CLE one-arm, 2026-09-05

Three sources not previously opened in the literature-officer notes. Theory input; **does not enter** `docs/STATUS.md`. X this pass: only 2023–2025 posts recycling the Nolin–Sun backbone; nothing on algebraic square-site `p_c` or Pinson wrapping.

Do not re-read as new: Ziff 2006, Camia–Feng 2508.16047, Diamantis–Kleban 0905.1727, He 2411.18696, Nolin backbone, Akhunzhanov polynomials.

---

## Q2 — Simmons–Kleban–Ziff: the derivative of a primary *is* a primary

**Simmons, Kleban, Ziff,** [arXiv:0705.1933](https://arxiv.org/abs/0705.1933), *J. Phys. A* (2007). *Percolation crossing formulas and conformal field theory.* This is the paper whose three new crossings Diamantis–Kleban later proved are weight-0 second-order modular forms.

The load-bearing sentence is in the abstract:

> The main step in our approach implies the identification of the derivative of one primary operator with another. We present operator identities that support this idea and suggest the presence of additional symmetry in `c=0` conformal field theories.

Body: `L_{-1} ψ_1(x) = K ψ_3(x)` with `K = 3^{1/4}/(2√π)`. Integrating the new first crossing densities recovers Cardy `Π_h`, Watts `Π_{h\bar v}`, and Cardy’s mean cluster-crossing count `ℕ_h`, without a Potts limit or a higher null vector.

**Why this is the Q2 mechanism, not just another modular-form citation.** A logarithmic partner is, at the level of Virasoro, a state on which `L_0` is not diagonal; equivalently, a derivative that refuses to stay inside the primary’s conformal family. SKZ *identify* that derivative with a different primary. Diamantis–Kleban then show the resulting crossings live in the weight-0 second-order space on `Γ(2)`, not in holomorphic weight-4. The N=290 exclusion of `Ê4(2i)/Ê4(i)=11/4` is what this literature predicts for a Cardy-type leftover, and is a surprise only if the leftover is assumed to be an Eisenstein series.

This still does **not** identify matching-odd as `ψ_1` or `ψ_3`. It names the *operation* (derivative of a primary = another primary) that produces the function space the last freeze did not include. Next freeze: either a `Γ(2)` second-order competitor built from Cardy/Watts, or a written non-claim that matching-odd is not a crossing density.

---

## Matching / P2 — Grimmett–Li 2024, primary

**Grimmett and Li,** [arXiv:2205.02734](https://arxiv.org/abs/2205.02734) v3 (20 Feb 2024), *Random Struct. Algorithms* **65** (2024) 832–856. *Percolation critical probabilities of matching lattice-pairs.*

Site percolation, quasi-transitive plane graphs, Euclidean or hyperbolic.

- Theorem 1.2: `p_c^{site}(G_*) < p_c^{site}(G)` iff `G_*` contains a doubly-infinite non-self-touching path that uses some diagonal of `G`.
- Transitive case: strict inequality **iff `G` is not a triangulation**.
- Complementary identity (their (1.3), from the companion hyperbolic paper): `p_u^{site}(G) + p_c^{site}(G_*) = 1`.

The square lattice is not a triangulation, so `p_c(Sq8) < p_c(Sq)` is a theorem. In Euclidean 2-d, uniqueness and percolation thresholds coincide (`p_u = p_c`), and the identity collapses to Sykes–Essam `p_c(Sq) + p_c(Sq8) = 1`. That is the matching pair the project is named after. **Neither threshold is algebraised.** Square *bond* is `1/2`; square *site* is not. Do not cite this paper as giving `p_c(ℤ^2)=1/2` for site percolation — that is the bond value.

P2 consequence, unchanged but now with a theorem rather than a folklore matching relation: an algebraic form for square-site `p_c` would automatically give one for Sq8, and conversely. Exclusion of low-degree/height forms on one is exclusion on the other. Wierman height 6 still bounds the *published exact* class; it is not a bound on this pair, which has no published exact form.

---

## Q1 — Liu–Sun–Yu–Zhuang, July 2026

**Liu, Sun, Yu, Zhuang,** [arXiv:2410.12724](https://arxiv.org/abs/2410.12724) v2 (16 July 2026). *The bulk one-arm exponent for the CLE_{`κ'`} percolations.* Same SLE/LQG/welding machine as the backbone exponent. Solves the bulk one-arm of CLE_{`κ'`} (`κ' ∈ (4,8)`). Special case: bichromatic one-arm of critical 3-state Potts is **`4/135`**.

Not percolation Bernoulli, not square, not spin-4, not `x=21/4`. What it changes: the 2023–2025 backbone result is no longer a one-off. The method now produces a family of exact exponents, some rational (`4/135`), some transcendental (backbone). **Monochromatic `k>2` for Bernoulli percolation is still unsolved** — the original Nolin–Sun paper and Beffara–Nolin 2009 both left that open, and this 2026 paper does not close it.

Matching-odd should not be assumed to have a rational dimension just because polychromatic arms are rational. Backbone is the warning; this paper is the reminder that both kinds exist inside the same method.

---

## Opinions

1. **Q2 / #567 freeze.** Name a Cardy/Watts second-order competitor, or non-claim. SKZ is the operator reason that space exists; Diamantis–Kleban is the classification; Pinson `π({1,0})` (ratio 2.969 at `r=2`) is a third, different number. Three named things, not one `11/4`.
2. **Q1.** Still do not pay to re-derive energy–hull. Still colour-decompose before any `21/4` fit. Monochromatic `k>2` remaining open is the honest status, not a gap a lattice log-slope will fill.
3. **P2 / matching.** Cite Grimmett–Li for the strict inequality and the matching identity; cite Ziff 2006 for “square site is outside the cell/dual-cell method”; cite Wierman for height 6. Do not mix bond `1/2` into the site table.
4. **X.** Dead for this project this month. arXiv remains the source.

## Not established

- matching-odd = `ψ_1` or `ψ_3`;
- an exact algebraic `p_c` for Sq or Sq8;
- monochromatic `k>2` exponents;
- anything in the claim ledger.
