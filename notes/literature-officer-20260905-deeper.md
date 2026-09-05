# Literature officer — deeper pass, 2026-09-05

**Role:** theory input. Same epistemic level as `notes/literature-officer-brief-20260905.md` (#571) and an Astra answer. **Does not enter** `docs/STATUS.md`. **Does not score** any frozen block.

This pass is what the first brief listed as open holes and then stopped. Sources were read from arXiv HTML/PDF, not from indexes. Nothing below is in #571 except as a named hole.

`#566` stays open until #564 lands and the A-lattice JSON status flips; this PR does not touch that artifact.

## Ticket map

| Ticket | What this pass changes |
|---|---|
| #565 / Q4 / P2 §4.2 | The historical class is not `C(≤6, ≤4)`. Wierman bow-tie (1984) is degree 5 **height 6**. Ziff–Scullard 2006 generalised bow-ties go to degree 11, coefficient height 36. Q4's first branch ("the mechanisms bound degree and height") is dead on the published record, not just on A-lattice. |
| P2 Result F′ | `C(≤6, ≤4)` is a truncated snapshot. The census null on that class still stands as a null on that class. It is no longer a null against "every exact planar threshold we can source". Say so in §8. |
| #567 / Q2 | Name three literature objects *before* N=580: (i) Pinson wrapping `R(r)` via `η, Z_{m,n}`; (ii) Arguin homology type `{a,b} ⊂ H_1(T²)`; (iii) the Newman–Ziff wrapping-estimator exponent `L^{-11/4}`, which is `1/ν+θ` and is **not** the modular `11/4`. |
| matching-odd homology | Duncan–Kahle–Schweinhart (AIHP 2025) is the theorem language for giant cycles on `T^d`. Arguin 2002 is the percolation (`Q=1`) special case of subgroup probabilities. |
| engine plumbing | Akhunzhanov–Eserkepov–Tarasevich 2022: exact square-site wrapping polynomials to `L=12`. Zero-new-compute check of the torus channel against a published polynomial, not against Monte Carlo. |

X: Deng's group posted two 2026 long-range papers (`2608.20750`, `2608.15120`); Newman–Grassberger–Ziff `2607.24975` is directed/strongly-connected percolation. None of them discuss algebraic square-site `p_c` or matching-odd. The first brief's community read still holds.

---

## 1. Q4 / P2 — Wierman 1984 is height 6; the cell can be arbitrarily expensive

The first brief said secondary indexes "suggest height 6" and left a primary reading as a later hole. The primary papers:

**Wierman, J. Phys. A 17, L229 (1984)**, as reproduced by Ziff and Scullard, [cond-mat/0610813](https://arxiv.org/abs/cond-mat/0610813) = J. Phys. A 39, 15083 (2006), verbatim:

```text
1 - p - 6 p^2 + 6 p^3 - p^5 = 0
```

with solution in `[0,1]` `p_c = 0.404518…` (bow-tie **bond**). Degree 5. Coefficient list `(1, -1, -6, 6, 0, -1)`. **Height 6.**

A GitHub code search of this repository for `bow-tie`, `bowtie`, and `Wierman` on 2026-09-05 returned **zero** hits. The polynomial is not in the P2 table that motivated `C(≤6, ≤3)` and then `C(≤6, ≤4)`.

Same Ziff–Scullard paper, further exact triangle-triangle outputs:

| lattice | polynomial (as printed) | deg | height |
|---|---|---:|---:|
| Wierman bow-tie | `1-p-6p^2+6p^3-p^5` | 5 | **6** |
| bow-tie (d) | `1-2p^2-3p^3+4p^4-p^5` | 5 | 4 |
| martini-A (already in #564) | `p^5-4p^4+3p^3+2p^2-1` | 5 | 4 |
| bow-tie (b) | `1-p-2p^3-4p^4-4p^5+15p^6+13p^7-36p^8+19p^9+p^{10}+p^{11}` | **11** | **36** |
| bow-tie (c) | `1-2p^3-2p^4-2p^5-7p^6+18p^7+11p^8-35p^9+21p^{10}-4p^{11}` | **11** | **35** |

Ziff, Scullard, Wierman, Sedlock, [arXiv:1210.6609](https://arxiv.org/abs/1210.6609) = J. Phys. A 45, 494005 (2012), add asymmetric bow-ties. One is stated as rigorous (no negative-probability step):

```text
1 - p - p^2 - 4 p^3 - 2 p^4 + 15 p^5 - 10 p^6 + p^8 = 0
```

`p_c = 0.481216…`, degree 8, height 15. Two further degree-7/8 formulae in that paper are **conjectural** (the manifold argument is rigorous only on a region); do not put those in the exact class until the conjecture is closed.

**What this does to Q4.** The three mechanisms produce a polynomial whose degree and height grow with the cell. There is no published uniform bound. A-lattice broke height 3; Wierman 1984 already sat at height 6, twenty years earlier; generalised bow-ties sit at height 36. Q4's remaining live half is only: *is square site outside every reachable cell, not just the 2006 ones?* Ziff 2006 still answers 2006. Wierman does not answer square **site** either — these are bond thresholds on triangulated-cell graphs.

**What this does to P2.** Result F and F′ are nulls on `C(≤6, ≤3)` and `C(≤6, ≤4)`. Those nulls do not move. The **motivating sentence** "every exactly-known planar threshold we can source is in that class" is false, and was false in 1984. §8 should say:

- the census is exhaustive inside the class we chose;
- the class is not the historical record;
- square site is not known to lie in any larger cell of the same mechanisms, and is named by Ziff as outside the 2006 method.

Do not re-run the census at height 36 in this fire. That is a manuscript-scoping decision (#565), not a literature finding.

**Scullard, PRE 86, 041131 (2012),** "percolation critical polynomial as a graph invariant": on *unsolved* lattices the same method yields high-degree polynomials whose roots are **approximations** (kagome bond to `~10^{-6}`, not exact). Those polynomials do not belong in an "exact thresholds" table. Jacobsen–Scullard 2020 (Phys. Rev. Research 2, 012050) then treat the eigenvalue formulation as a numerical estimator, e.g. kagome `p_c = 0.52440499916744820(1)`. Distinct from an algebraic identity.

---

## 2. #567 / Q2 — three objects to name before N=580, not after

The frozen N=580 competitors are weight-4, bare aspect ratio, area, none. The first brief asked for Pinson as a fifth *class*. This pass supplies the actual language and a plumbing check.

### 2.1 Pinson wrapping, and Arguin's homology refinement

**Pinson, J. Stat. Phys. 75, 1167 (1994).** Wrapping / winding probabilities on the torus as combinations of Dedekind `η` and `Z_{m,n}(g; r)` at Coulomb-gas `g = 2/3`. Pruessner–Moloney, [cond-mat/0310361](https://arxiv.org/abs/cond-mat/0310361), confirm it numerically and record that at `r=2` their truncated series for `P((1,0), ≥1, r)` agrees with Pinson to relative `< 10^{-8}`.

**Arguin, J. Stat. Phys. 109, 301 (2002),** [hep-th/0111193](https://arxiv.org/abs/hep-th/0111193). The observable is not a single wrapping number. It is the subgroup of `H_1(T²) ≅ ℤ×ℤ` generated by the FK clusters: `{0}`, `ℤ×ℤ` (cross topology), or `{a,b}` with `a ∧ b = 1`. Pinson is the `Q=1` case of this. Modular identities:

```text
π_Q({(a,b)})|_τ     = π_Q({(a+b,b)})|_{τ+1}
π_Q({(a,b)})|_τ     = π_Q({(-b,a)})|_{-1/τ}
```

Monte Carlo on `τ = i`, `2i`, `i+1/2` agrees for `Q=1,2,3`. For `Q=4`, logarithmic corrections are suspected. Percolation is `Q=1`; matching-odd is not a Potts `Q` we have identified. The classification still applies: a wrapping-type lattice observable has a homology type, and the N=580 ladder is three values of `τ = ir` with `r ∈ {1,2,4}`.

**Morin-Duchesne–Saint-Aubin, [0812.2925](https://arxiv.org/abs/0812.2925).** Thin-torus (`τ_i → ∞`) decay of `π({1,0})` is governed by extended-Kac weights, including **half-integers**: for percolation, `γ_1 = 2 h_{1/2,0}` (Cardy's hull). This is the continuum prediction for how wrapping probabilities collapse as `r` grows, which is exactly the column the ladder measures. It is not `E_4`.

**Action on #567, before the run.** Either:

1. compute Pinson / Arguin `π({1,0})(ir)` at `r=1,2,4` from the published formula and add those three numbers as named competitors, or
2. write in the frozen file that the matching-odd slope is **not** claimed to be a Pinson wrapping, so those numbers are non-claims.

Scoring only against `E_4` and then reading `r` off the leftover is how N=290 acquired a post-hoc `bare_aspect_ratio`. The thin-torus expansion of `π({1,0})` is a function of `r` that can sit near `r` itself at these three points without being `E_4` and without being the identity map. That is the alternative the first brief asked for and this pass names.

No published Pinson table for Sq8 (NN+NNN) was found, again. Universality still says critical wrapping equals Pinson; the matching identity constrains the *difference*. Square wrapping vs Pinson, matching wrapping vs the same Pinson, difference vs Mertens–Ziff, remains the zero-new-compute check.

### 2.2 Do not identify two different `11/4`s

Newman–Ziff wrapping estimators for square-site `p_c` converge as `L^{-11/4}`. Martins–Plascak, [cond-mat/0304024](https://arxiv.org/abs/cond-mat/0304024), write it explicitly: `ν = 4/3` and correction-to-scaling `θ = 2` give `1/ν + θ = 11/4`.

The fingerprint's `11/4` is `r^2 / E_4(i)` at `r=2`, a modular-weight statement.

These are the same rational and different objects. A wrapping *threshold estimator* scaling as `L^{-11/4}` does not predict that a *spin-4 amplitude ratio* equals `11/4`. If the N=580 write-up conflates them, the competitor list becomes unreadable.

### 2.3 Exact wrapping polynomials on the square torus, `L ≤ 12`

**Akhunzhanov, Eserkepov, Tarasevich, J. Phys. A 55, 204004 (2022),** [2204.01517](https://arxiv.org/abs/2204.01517), Ziff-70th-birthday issue. Exact polynomials for square-**site** percolation probability:

- plane crossing, `L ≤ 17`;
- cylinder spanning, `L ≤ 16`;
- **torus wrapping along one direction, `L ≤ 12`.**

The polynomials live in the supplement (they are large). Divisibility identities are proved. Naive FSS on the wrapping polynomials gave `p_c = 0.59269`, worse than Jacobsen, as expected at these sizes.

This is the plumbing check the engine has not been run against: at `L ≤ 12` the torus wrapping channel has a published polynomial, not an MC error bar. If our wrapping observable, evaluated by exact enumeration or by the transfer engine at those sizes, disagrees with the polynomial, the channel is wrong. If it agrees, Pinson comparisons at production sizes are comparisons of a checked channel to a continuum formula.

They do not treat NN+NNN / Sq8.

---

## 3. Matching-odd homology has a theorem name

**Duncan, Kahle, Schweinhart,** *Ann. Inst. H. Poincaré Probab. Statist.* **61** (2025) 2235–2261, [arXiv:2011.11903](https://arxiv.org/abs/2011.11903). Homological percolation on `T^d`: giant `i`-cycles, sharp transition, and `p_c = 1/2` in middle dimension `i = d/2`. This is the high-d analogue of Kesten, in the language the project already uses (ambient homology, giant cycles, torus).

It does **not** identify matching-odd. It tells P3 / #275 what a published theorem about wrapping-as-homology looks like, and that the 2-d case is wrapping (`i=1` on `T^2`), which is Arguin/Pinson, not a new exponent.

A 2026 follow-up, [2601.00793](https://arxiv.org/abs/2601.00793) (Voronoi percolation, topological stability), cites them and repeats `p_c = 1/2` for `i`-dimensional giant cycles on the `2i`-torus. Orthogonal geometry, same definition.

---

## 4. Opinions for subsequent analysis

1. **#565 / Q4.** Rewrite the question so the first branch is closed: mechanisms do not bound height. The live question is square-site reachability. Do not freeze a new height cap; Wierman is 6 and generalised bow-tie is 36.
2. **P2 §8.** Drop "historical complexity range" as a law. Keep the census. Cite Wierman 1984 and Ziff–Scullard 2006 as the reason the class was a choice. Grimmett–Li 2024 and Nolin–Sun remain the right neighbours for matching and for "exact ≠ algebraic".
3. **#567, before any N=580 sample.** Add Pinson/Arguin `π({1,0})(ir)` at `r=1,2,4` as a competitor *or* an explicit non-claim. Name the Newman–Ziff `L^{-11/4}` so it cannot be read as the modular `11/4`.
4. **Engine.** Evaluate the committed wrapping channel on `L ≤ 12` against Akhunzhanov polynomials. That is cheaper than another 200 M-sample rung and is the only published exact torus polynomial for square site.
5. **Colour.** Still do the polychromatic vs monochromatic decomposition of matching-odd before fitting `x=21/4` (first brief). This pass adds nothing that licenses skipping it.
6. **#566.** Unchanged: primary A-lattice reading is on #571; JSON status flip waits on #564. This pass does not reopen it.
7. **Q1.** Arguin `Q=4` log corrections are a warning, not a coupling computation. He 2024 still owns `μ = -5/4`. Thin-torus `h_{1/2,0}` is a hull exponent, not the level-4 `x=21/4, s=4` pairing.
8. **Do not** put Scullard critical-polynomial roots for kagome / Archimedean lattices into the exact-threshold table; they are estimators.

## Not established, and not claimed here

- that matching-odd *is* Arguin type `{1,0}`;
- a numerical Pinson value at `r=2` or `r=4` (compute from the formula, do not guess);
- that the degree-11 bow-tie polynomials are all irreducible (height is read off the printed primitive polynomials; factoring would only *lower* degree);
- a new census at height 36;
- anything that belongs in the claim ledger.
