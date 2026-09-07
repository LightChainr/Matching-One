# Literature officer synthesis — 2026-09-05

**Role:** theory input. Same epistemic level as an Astra answer. **Does not enter** `docs/STATUS.md`. **Does not score** a frozen block.

This file is the decision map. Primary quotes live in the per-source notes listed at the end. First-pass names that this packet only filled: `notes/literature-officer-brief-20260905.md` (#571).

X / community this month: Komargodski 2026-09-03 on `θ(p_c)=0` is 3-d continuity (2-d is Harris–Kesten; 3–10 still open in [2511.01851](https://arxiv.org/abs/2511.01851) v2). Nobody is discussing algebraic square-site `p_c`, matching-odd, or Pinson wrapping. arXiv is the channel.

---

## Ticket map

| Ticket | What the literature now forces | What it does **not** |
|---|---|---|
| **Q1** | The only constructed log pair is energy–hull (`x=5/4`, `γ=−5/4`). Connectivity 4-points have **no even spin**. Backbone CFT has **complex** subleading dimensions. Colour-decompose before any `21/4` fit. | matching-odd = spin-4 / rank-4 / `x=21/4` |
| **Q2 / #567** | N=290 tested holomorphic `E₄` (`11/4`) and excluded it. The leftover, if Cardy/wrapping-type, lives in a different function space. Freeze a named competitor or non-claim **before** N=580. | leftover = `E₄`; leftover = Pinson without naming which `R` |
| **Q3** | Discrete holomorphy / C4 (Cardy) is a different lattice field than matching-odd. He’s rank-3 is `ε'`–`T T̄`, not level-4. | a lattice C4 fit as Q1 |
| **Q4 / #565** | Degree/height of exact thresholds grow with the cell. Wierman bow-tie is deg 5 **height 6**. Square site is outside Ziff 2006 *and* outside Jacobsen’s “no finite-size dependence” class. | a theorem that all exact `p_c` live in `C(≤6,≤4)` |
| **P2 / #566** | A-lattice quintic and height 4 stand (Ziff PRE 73, 016134, quotes in #574 note). Square site has **no** published exact form. Jacobsen 2015 `0.59274605079210(2)` is an estimator. | square-site algebraic `p_c`; height 4 as a law |
| **#576 / wrapping** | Four published numbers at `r=1`. Engine must name which. Aspect `r=2,4` changes them. | unnamed wrapping vs `11/4` |
| **matching-odd homology** | Duncan–Kahle giant 1-cycle **is** wrapping, `p_c=1/2` only for square **bond** / triangular **site**. Matching-odd is the Sq/Sq8 involution, a different object. | matching-odd = giant cycle; square-site wrapping threshold = `1/2` |

---

## Q1 — do not pay to re-derive energy–hull; do not fit `21/4` blind

Constructed log pair, now a theorem on the triangular lattice:

- He [2411.18696](https://arxiv.org/abs/2411.18696): `γ^perco = −5/4` at energy–hull (`h=5/8`). Rank-3 of `ε'` with `T T̄`, coupling `a=−25/48`. Liouville comparison “suggests arbitrarily high rank” — **not constructed** at spin-4.
- Camia–Feng [2508.16047](https://arxiv.org/abs/2508.16047) v2 (1 Jun 2026): lattice energy field + four-arm partner; two-point `|z|^{-5/2}`. Triangular site. No `21/4`, no spin-4, no square.
- Vasseur–Jacobsen–Saleur [1206.2312](https://arxiv.org/abs/1206.2312): lattice pure log `F(r) ∼ θ + (2√3/π) log r`, energy mixed with two-cluster field. Matching-odd is not this counting.

Colour / spin:

- Tan 2019: `x=21/4` is scalar `P_{4s}` (spin 0), not spin-4.
- Picco–Ribault–Santachiara [1607.07224](https://arxiv.org/abs/1607.07224): bootstrapped connectivity 4-points live in `ℳ_{2ℤ,ℤ+1/2}`, **no even spin**. Leading `(5/96, 5/96)`.
- Radhakrishnan–Tassion [2410.23250](https://arxiv.org/abs/2410.23250): mono `>` poly, strictly.
- Nolin backbone: transcendental, `0.356666…`. Sun–Xu–Zhuang [2410.04767](https://arxiv.org/abs/2410.04767): the same elementary equation has **countably many complex roots** as subleading annulus exponents.
- Liu–Sun–Yu–Zhuang [2410.12724](https://arxiv.org/abs/2410.12724) v2 (16 Jul 2026): CLE_κ' one-arm; 3-state Potts `4/135`. Monochromatic Bernoulli `k>2` still open.

**Paid Q1 query is only:** is the level-4 spin-4 pair fixed by `μ=−5/4`? Everything else is in the notes.

---

## Q2 / #567 — four named competitors, one tested

N=290 excluded `E₄(2i)/E₄(i)=11/4` at 4.9σ (`1.880 ± 0.177`). That is the holomorphic `c
eq0` cousin of `⟨T T⟩ ∝ c/z^4`, which **vanishes at `c=0`**.

| object | weight / type | log? | `11/4`? | source |
|---|---|---|---|---|
| holomorphic `E₄` | 4, `SL(2,ℤ)` | no | yes | tested, excluded |
| Ridout `⟨t t⟩` | 4 × log | yes | no | [1303.0847](https://arxiv.org/abs/1303.0847): `[A+(5/4)log]/(z−w)^4` |
| Cardy/Watts crossings | 0, second-order, `Γ(2)` | via SKZ `∂ψ_1=K ψ_3` | no | [0705.1933](https://arxiv.org/abs/0705.1933), [0905.1727](https://arxiv.org/abs/0905.1727) |
| Pinson wrapping | `η`, `Z_{m,n}`, `g=2/3` | no | ratio **2.969** at `r=2` | 1994; #576 note |

SKZ: at `c=0` the derivative of a primary *is* a primary. That is why the crossings are second-order, not Eisenstein.

**Freeze before N=580:** at least one row other than `E₄`, or a written non-claim that matching-odd is not a crossing/wrapping density. Leaving `E₄` as the only modular shape is how N=290 excluded a list that did not contain the actual functions.

---

## #576 — name the wrapping event

Newman–Ziff (PRE 64, 016706) quoting Pinson at **square aspect, criticality**:

| symbol | event | value | Duncan–Kahle analogue |
|---|---|---|---|
| `R_h` | specified direction | **0.521 058 290** | A, specified |
| `R_e` | either or both | **0.690 473 725** | A |
| `R_b` | both directions | **0.351 642 855** | S (homology basis) |
| `R_1` | one specified axis, not the other | **0.169 415 435** | A minus S |

Identities: `R_b = R_e − 2 R_1`, `R_h = R_e − R_1`. The #576 computation of `π({1,0})(r=1)` **is** `R_h` (matched to 12 digits). At `r=2`, `π({1,0})` ratio vs `r=1` is **2.969**, not `11/4=2.75`.

Akhunzhanov [2204.01517](https://arxiv.org/abs/2204.01517): exact square-site wrapping polynomials, torus `L≤12`. Engine comparison is remaining compute, not literature.

Duncan–Kahle [2011.11903](https://arxiv.org/abs/2011.11903) / Schweinhart–Shuman [2601.00793](https://arxiv.org/abs/2601.00793): giant 1-cycle = wrapping; `p_c=1/2` for square **bond** and triangular **site**, not square site. Event A ≠ event S.

**#576 freeze:** `R_h` / `R_e` / `R_b` / `R_1`, and whether the ladder is at `r=1` only. Matching-odd is not claimed to be any of them.

---

## Q4 / P2 / #566 — square site is outside the exact class

Exact published class (do not freeze as a law):

- Ziff 2006 A-lattice: `p^5-4p^4+3p^3+2p^2-1=0`, height 4. Square site named **outside** the cell/dual-cell method. Quotes verified in the #574 note.
- Wierman bow-tie 1984: `1-p-6p^2+6p^3-p^5=0`, deg 5 **height 6**. Generalised bow-ties deg 11 height 36. One asymmetric deg 8 height 15 (2012).
- Jacobsen, *J. Phys. A* **48** (2015) 454003: “in exactly solvable cases there is **no finite-size dependence** at all.” Square-site critical polynomial is extrapolated to `n=21` → `p_c = 0.59274605079210(2)`. Author’s own split: not exact.
- Grimmett–Li RSA 2024 [2205.02734](https://arxiv.org/abs/2205.02734): `p_c(G_*)<p_c(G)` iff not a triangulation. Square is not, so `p_c(Sq8)<p_c(Sq)` is a theorem. With `p_u=p_c` this is Sykes–Essam `p_c(Sq)+p_c(Sq8)=1`. **Neither is algebraised.** Do not cite as site `p_c(ℤ^2)=1/2` (that is bond).
- Mertens–Ziff [1603.07289](https://arxiv.org/abs/1603.07289): finite-size matching identity; unique root of the matching polynomial ∼ `L^{-4}`.

P2 novelty is unchanged: certified exclusion of low-degree/height forms against the published intervals, not a new exact `p_c`.

---

## Subsequent analysis (ordered)

1. **#576 / #567 freeze (before N=580).** Write which of `{R_h, R_e, R_b, R_1, Cardy/Watts second-order, Ridout-t, non-claim}` is on the competitor list. Include `pinson_pi10_ratio` at `r=2` (2.969) if the ladder is wrapping-flavoured. Do not retest only `11/4`.
2. **Colour decomposition, before any `21/4` fit.** Tan scalar vs Tassion mono>poly vs Picco no-even-spin. If the leftover is even-spin it is not the bootstrapped connectivity 4-point.
3. **Do not pay Q1** to re-derive energy–hull. Query is rank-4 / spin-4 at `μ=−5/4`, or skip.
4. **P2 / #566.** Cite Ziff 2006 (outside the method) + Jacobsen 2015 (estimator, not exact) + Wierman height 6 (published exact class is already taller than 4). Keep `0.59274605079210(2)` out of the exact table.
5. **Engine (compute, not literature).** Akhunzhanov `L≤12` wrapping polynomials vs the wrapping channel, once the channel is named.
6. **Matching-odd homology.** Treat as Sq/Sq8 involution. Do not import Duncan `p_c=1/2`.
7. **X.** Watch Komargodski / Deng / Sun. This month is empty for our tickets.

---

## Source notes (primary quotes)

| file | fills |
|---|---|
| `literature-officer-brief-20260905.md` (#571) | first map; Ziff 2006 |
| `ziff-2006-a-lattice-primary-reading.md` | #566 scientific half |
| `literature-officer-20260905-deeper.md` | Wierman height 6; Arguin |
| `literature-officer-20260905-issue574-quotes.md` | Ziff / Suding / Scullard quotes |
| `literature-officer-20260905-issue576-wrapping.md` | Pinson `π({1,0})` at `r=1,2,4`; Akhunzhanov |
| `literature-officer-20260905-q1-he-q3.md` | He `γ=−5/4` |
| `literature-officer-20260905-q1-noise.md` | Tan `21/4` scalar |
| `literature-officer-20260905-q1-q2-followup.md` | Camia–Feng 2026; Diamantis–Kleban |
| `literature-officer-20260905-skz-grimmett-cle.md` | SKZ derivative=primary; Grimmett–Li; Liu–Sun 2026 |
| `literature-officer-20260905-vasseur-ridout-jacobsen.md` | Vasseur `F(r)`; Ridout `t`; Jacobsen 2015 split |
| `literature-officer-20260905-picco-x.md` | Picco no even spin; Komargodski `θ(p_c)` is 3-d |
| `literature-officer-20260905-homology-giant.md` | Duncan giant = wrapping; A ≠ S |
| `literature-officer-20260905-newman-annulus.md` | four `R_*`; complex annulus roots |
| `literature-officer-20260905-mertens-matching.md` | finite-size matching identity |

---

## Not established, and not claimed here

- matching-odd = Vasseur `F(r)`, Ridout `t`, Picco `R_σ`, Pinson `R_h`, or Duncan giant cycle;
- a square-site algebraic `p_c`;
- rank-4 Jordan at `x=21/4`;
- `θ(p_c)=0` in 3-d as a refereed theorem;
- anything in the claim ledger.
