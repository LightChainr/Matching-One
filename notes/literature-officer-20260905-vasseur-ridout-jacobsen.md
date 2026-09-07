# Literature officer — Vasseur, Ridout t-field, Jacobsen 2015, 2026-09-05

Three primary sources not opened in earlier literature-officer notes. Theory input; **does not enter** `docs/STATUS.md`. X this pass: a 2024 joke equating the 3-d polymer exponent ~0.588 with square-site `p_c` ~0.593, and a 2022 undergrad talk on Wierman bounds. No researcher thread on algebraic square-site `p_c`, wrapping, or LCFT.

Do not re-read as new: Ziff 2006, Camia–Feng 2508.16047, He 2411.18696 `γ=−5/4`, Diamantis–Kleban, SKZ 0705.1933, Grimmett–Li, Nolin backbone, Akhunzhanov, Liu–Sun 2410.12724.

---

## Q1 — Vasseur–Jacobsen–Saleur: a lattice log that is *not* matching-odd

**Vasseur, Jacobsen, Saleur,** [arXiv:1206.2312](https://arxiv.org/abs/1206.2312), *Phys. Rev. Lett.* (2012). Bond percolation as `Q→1` Potts.

They construct a four-point cluster observable whose two-point function, after a subtraction that cancels the leading `r^{-5/2}` power, is a **pure logarithm**:

> `F(r) ≡ [ℙ₀(r) + ℙ₁(r) − ℙ≠²] / ℙ₂(r) ∼ θ + (2√3/π) log r`

The prefactor `2√3/π ≈ 1.1026` is universal. The log comes from mixing the energy operator with “the field that creates two propagating clusters.” Monte Carlo slope `1.15 ± 0.05`. This is the lattice ancestor of Camia–Feng’s energy–hull pair (same `x=5/4` channel).

**What it is not.** It is not matching-odd, not spin-4, not `x=21/4`, not square-site. A log-slope in a matching-odd readout is not automatically this object. To claim it is, the engine would have to be counting two-cluster / four-point FK events, which it is not.

He 2411.18696 (already noted for `γ^ᵖᵉʳᶜᵒ=−5/4`) adds one sentence that the earlier note underweighted: comparison with `c<1` Liouville “suggests the potential existence of **arbitrarily high rank Jordan blocks**.” The blocks He actually constructs are rank-2 (energy–hull) and rank-3 (`ε'` mixed with `T T̄`, coupling `a=−25/48`). Level-4 / spin-4 / `x=21/4` is still not constructed. The paid Q1 query is whether that potential is realised at spin-4, not whether energy–hull exists.

---

## Q2 — Ridout: the log partner of `T` is weight-4 *with a log*

**Ridout,** [arXiv:1303.0847](https://arxiv.org/abs/1303.0847). *Logarithmic Conformal Field Theory: Beyond an Introduction.* Boundary percolation CFT. The log partner `t` of the stress tensor satisfies

> `⟨t(z) t(w)⟩ = [A + (5/4) log(z−w)] / (z−w)^4`

(`A` depends on the choice of `|t⟩`; the `5/4` does not.) Anomaly numbers `b_{1,5}=−5/8`, `b_{3,1}=5/6`, later measured on the lattice.

This is the object whose two-point looks like a weight-4 density times a logarithm. Holomorphic `E₄` is the *diagonal* `c
eq0` cousin (`⟨T T⟩ ∝ c/z^4`), which **vanishes at `c=0`**. The thing that survives is `t`, and `t` is second-order (a log times `z^{-4}`), not a holomorphic Eisenstein series.

Together with SKZ (`∂ψ_1 = K ψ_3`) and Diamantis–Kleban (crossings are weight-0 second-order on `Γ(2)`), the modular-form menu for a `c=0` leftover is:

| object | weight | log? | `11/4`? |
|---|---|---|---|
| holomorphic `E₄` | 4 | no | yes, `E₄(2i)/E₄(i)=11/4` |
| Ridout `⟨t t⟩` | 4 × log | yes | no |
| Cardy/Watts crossings | 0, second-order, `Γ(2)` | via SKZ derivative | no |
| Pinson `π({1,0})` | wrapping, `r`-dependent | no | ratio 2.969 at `r=2` |

N=290 tested the first row and excluded it. The other three rows were not on the list. Next freeze should name at least one of them or non-claim all of them.

---

## P2 — Jacobsen 2015, in the author’s own split

**Jacobsen,** *J. Phys. A* **48** (2015) 454003, [arXiv:1507.03027](https://arxiv.org/abs/1507.03027). Critical polynomial `P_B(q,T)` on a basis of `n×m` cells.

> Moreover, in exactly solvable cases there is no finite-size dependence at all.

He then treats three cases that *do* have finite-size dependence, including **site percolation on the square lattice to `n_max=21`**, and reports the extrapolated estimator

> `p_c = 0.592 746 050 792 10(2)`

That number is the best published estimate (Akhunzhanov cites it as such). It is **not** a root of a finite-cell polynomial that is independent of cell size. By Jacobsen’s own criterion, square site is outside the exact class. Keep it out of the exact-`p_c` table. Ziff 2006 already said square site is outside the cell/dual-cell method; this is the same split, from the author of the method that later produced the 14-digit estimate.

Kagome bond, by contrast, is also in the extrapolated class (`0.524 404 999 167 439(4)`). The exact Archimedean cases remain the ones with no finite-size drift.

---

## Opinions

1. **Q1.** Vasseur’s `F(r)` is the lattice log to compare against, and it is energy–two-cluster, prefactor `2√3/π`. Matching-odd is a different counting. Colour first; do not fit `21/4`; do not pay to re-derive energy–hull. A Q1 query is only the high-rank/spin-4 question He left open.
2. **Q2 / #567.** Four named competitors, one tested. Put Ridout-`t`, Cardy/Watts, and `pinson_pi10_ratio` on the next list, or write the non-claim.
3. **P2 / #566.** Jacobsen 2015 is the primary source for “critical polynomial of square site is an estimator, not an exact form.” Quote the no-finite-size-dependence sentence. Wierman height 6 still bounds the published exact class.
4. **X.** Still empty. Stop expecting it this month.

## Not established

- matching-odd = Vasseur `F(r)` or Ridout `t`;
- a square-site algebraic `p_c`;
- rank-4 Jordan at `x=21/4`;
- anything in the claim ledger.
