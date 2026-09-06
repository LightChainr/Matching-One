# Literature officer synthesis — 2026-09-05

**Role:** theory input. Same epistemic level as an Astra answer. **Does not enter** `docs/STATUS.md`. **Does not score** a frozen block.

This file is the decision map. Primary quotes live in the per-source notes listed at the end. First-pass names that this packet only filled: `notes/literature-officer-brief-20260905.md` (#571).

**Erratum (same packet):** `notes/literature-officer-20260905-erratum-pi-vs-Rh.md`. `π({1,0})(i)=0.169415435321` is **not** Newman–Ziff `R_h=0.521058290`. The #576 wrapping note was right; an earlier draft of this synthesis was not.

X / community: Komargodski 2026-09-03 on `θ(p_c)=0` is 3-d continuity (2-d is Harris–Kesten; 3–10 still open in [2511.01851](https://arxiv.org/abs/2511.01851) v2). Nobody is discussing algebraic square-site `p_c` or matching-odd. arXiv is the channel.

Overnight context (not this PR): N=580 has been run and replayed ([#577](https://github.com/LightChainr/Matching-One/pull/577)); the next modular freeze is [#585](https://github.com/LightChainr/Matching-One/issues/585) (map-resolved torus, solution spaces not isolated rays), not a re-freeze of N=580 as if it had not happened.

---

## Ticket map

| Ticket | What the literature now forces | What it does **not** |
|---|---|---|
| **Q1** | The only constructed log pair is energy–hull (`x=5/4`, `γ=−5/4`). Connectivity 4-points have **no even spin**. Backbone CFT has **complex** subleading dimensions. Colour-decompose before any `21/4` fit. | matching-odd = spin-4 / rank-4 / `x=21/4` |
| **Q2 / #567 / #585** | N=290 tested holomorphic `E₄` (`11/4`) and excluded it. Leftover, if Cardy/wrapping-type, lives in a different function space. Fit solution spaces, not one ray. | leftover = `E₄`; leftover = Pinson without naming `π` vs `R_h` |
| **Q3** | Discrete holomorphy / C4 (Cardy) is a different lattice field than matching-odd. He’s rank-3 is `ε'`–`T T̄`, not level-4. | a lattice C4 fit as Q1 |
| **Q4 / #565** | Degree/height of exact thresholds grow with the cell. Wierman bow-tie is deg 5 **height 6**. Square site is outside Ziff 2006 *and* outside Jacobsen’s “no finite-size dependence” class. | a theorem that all exact `p_c` live in `C(≤6,≤4)` |
| **P2 / #566** | A-lattice quintic and height 4 stand (Ziff PRE 73, 016134). Square site has **no** published exact form. Jacobsen 2015 `0.59274605079210(2)` is an estimator. | square-site algebraic `p_c`; height 4 as a law |
| **#576 / wrapping** | `π({1,0})` → `0.169…` (`R_1`-class). `R_h` → `0.521…` (Akhunzhanov). They are different. Name which. `pinson_pi10_ratio` at `r=2` is **2.969**. | unnamed wrapping vs `11/4`; `π = R_h` |
| **matching-odd homology** | Duncan–Kahle giant 1-cycle **is** wrapping, `p_c=1/2` only for square **bond** / triangular **site**. Matching-odd is the Sq/Sq8 involution. | matching-odd = giant cycle; square-site wrapping threshold = `1/2` |

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

## Q2 / #585 — four named competitors, one tested; fit spaces not rays

N=290 excluded `E₄(2i)/E₄(i)=11/4` at 4.9σ (`1.880 ± 0.177`). That is the holomorphic `c\neq0` cousin of `⟨T T⟩ \propto c/z^4`, which **vanishes at `c=0`**.

| object | weight / type | log? | `11/4`? | source |
|---|---|---|---|---|
| holomorphic `E₄` | 4, `SL(2,\mathbb{Z})` | no | yes | tested, excluded at N=290 |
| Ridout `⟨t t⟩` | 4 × log | yes | no | [1303.0847](https://arxiv.org/abs/1303.0847): `[A+(5/4)log]/(z−w)^4` |
| Cardy/Watts crossings | 0, second-order, `Γ(2)` | via SKZ `∂\psi_1=K \psi_3` | no | [0705.1933](https://arxiv.org/abs/0705.1933), [0905.1727](https://arxiv.org/abs/0905.1727) |
| Pinson `π({1,0})` | `{1,0}` class | no | ratio **2.969** at `r=2` | 1994; #576 note |

SKZ: at `c=0` the derivative of a primary *is* a primary. That is why the crossings are second-order, not Eisenstein. #585 is the ticket that already asks to fit modular-covariance **solution spaces**, which is this table rather than another isolated `11/4` ray.

N=580 (#577) was underpowered: `bare_aspect_ratio` survives; weight-4 is not excluded at 3σ on `r=4`. That does not license skipping the named wrapping competitor on the next design.

---

## #576 — `π({1,0})` and `R_h` are different numbers

Newman–Ziff (PRE 64, 016706) quoting Pinson at **square aspect, criticality**:

| symbol | event (their words) | value |
|---|---|---:|
| `R_h` | wraps horizontally (specified direction, **including** the other) | **0.521 058 290** |
| `R_e` | wraps either horizontally or vertically, or both | **0.690 473 725** |
| `R_b` | wraps in both directions | **0.351 642 855** |
| `R_1` | wraps around one specified axis **but not** the other | **0.169 415 435** |

Identities (digits): `R_b = R_e − 2 R_1`, `R_h = R_e − R_1`.

The #576 computation of `π({1,0})(r=1)` is **0.169415435321**. That matches `R_1` to 12 digits. It does **not** match `R_h`. Akhunzhanov’s torus polynomials are `R_h` (specified-direction wrap, including both). See the erratum.

`pinson_pi10_ratio` is `π({1,0})(ir)/π({1,0})(i)`:

| `r` | `π({1,0})` | ratio |
|---:|---:|---:|
| 1 | 0.169415435321 | 1 |
| 2 | 0.503035897695 | **2.969** |
| 4 | 0.855969321054 | **5.052** |

At `r=2` this sits 8% from `11/4=2.75`. N=290’s `1.880 ± 0.177` is itself ~6σ from `2.969`, so Pinson does not explain that run; it still has to be **named** on the next list.

Duncan–Kahle [2011.11903](https://arxiv.org/abs/2011.11903) / Schweinhart–Shuman [2601.00793](https://arxiv.org/abs/2601.00793): giant 1-cycle = wrapping; `p_c=1/2` for square **bond** and triangular **site**, not square site. Event A (at least one giant cycle) is closer to `R_e`; event S (homology basis) is closer to `R_b`. Analogues, not identities.

**#576 freeze, already written in the wrapping note:** matching-odd is a **non-claim** vs Pinson wrapping. `pinson_pi10_ratio` is a named competitor. Keep `R_h` (→ `0.521`) and `π({1,0})` (→ `0.169`) apart. Engine vs Akhunzhanov `L=3,4` is remaining compute (`L=10` ancillary file is corrupt).

---

## Q4 / P2 / #566 — square site is outside the exact class

- Ziff 2006 A-lattice: `p^5-4p^4+3p^3+2p^2-1=0`, height 4. Square site named **outside** the cell/dual-cell method. Quotes in the #574 note. #566 is still the P2 blocking primary-reading status flip (JSON), not this packet.
- Wierman bow-tie 1984: `1-p-6p^2+6p^3-p^5=0`, deg 5 **height 6**. Generalised bow-ties deg 11 height 36.
- Jacobsen, *J. Phys. A* **48** (2015) 454003: “in exactly solvable cases there is **no finite-size dependence** at all.” Square-site critical polynomial is extrapolated to `n=21` → `p_c = 0.59274605079210(2)`. Author’s own split: not exact.
- Grimmett–Li RSA 2024 [2205.02734](https://arxiv.org/abs/2205.02734): `p_c(G_*)<p_c(G)` iff not a triangulation. Square is not. With `p_u=p_c` this is Sykes–Essam. **Neither threshold is algebraised.** Do not cite as site `p_c(\mathbb{Z}^2)=1/2` (that is bond).
- Mertens–Ziff [1603.07289](https://arxiv.org/abs/1603.07289): finite-size matching identity; unique root ∼ `L^{-4}`.

P2 novelty is unchanged: certified exclusion of low-degree/height forms against the published intervals, not a new exact `p_c`.

---

## Subsequent analysis (ordered)

1. **#585 / next ladder freeze.** Name `{pinson_pi10_ratio, R_h, Cardy/Watts second-order, Ridout-t, non-claim}` on the competitor list. Do not retest only `11/4`. N=580 already ran (#577); this is the *next* design, not a rewrite of that one.
2. **Keep `π({1,0})` and `R_h` apart.** 0.169 vs 0.521. Akhunzhanov checks `R_h`. The ratio competitor is `π`.
3. **Colour decomposition, before any `21/4` fit.** Tan scalar vs Tassion mono>poly vs Picco no-even-spin.
4. **Do not pay Q1** to re-derive energy–hull. Query is rank-4 / spin-4 at `μ=−5/4`, or skip.
5. **P2 / #566.** Cite Ziff 2006 + Jacobsen 2015 + Wierman height 6. Keep the 14-digit estimator out of the exact table. This PR does not flip the JSON status.
6. **Engine (compute).** Akhunzhanov `L=3,4` vs the wrapping channel, once the channel is named as `R_h`.
7. **Matching-odd homology.** Sq/Sq8 involution. Do not import Duncan `p_c=1/2`.

---

## Source notes (primary quotes)

| file | fills |
|---|---|
| `literature-officer-brief-20260905.md` (#571) | first map; Ziff 2006 |
| `ziff-2006-a-lattice-primary-reading.md` | #566 scientific half |
| `literature-officer-20260905-deeper.md` | Wierman height 6; Arguin |
| `literature-officer-20260905-issue574-quotes.md` | Ziff / Suding / Scullard quotes |
| `literature-officer-20260905-issue576-wrapping.md` | **source of truth** for `π` vs `R_h`; Akhunzhanov |
| `literature-officer-20260905-erratum-pi-vs-Rh.md` | this check |
| `literature-officer-20260905-q1-he-q3.md` | He `γ=−5/4` |
| `literature-officer-20260905-q1-noise.md` | Tan `21/4` scalar |
| `literature-officer-20260905-q1-q2-followup.md` | Camia–Feng 2026; Diamantis–Kleban |
| `literature-officer-20260905-skz-grimmett-cle.md` | SKZ; Grimmett–Li; Liu–Sun 2026 |
| `literature-officer-20260905-vasseur-ridout-jacobsen.md` | Vasseur `F(r)`; Ridout `t`; Jacobsen 2015 |
| `literature-officer-20260905-picco-x.md` | Picco no even spin |
| `literature-officer-20260905-homology-giant.md` | Duncan giant = wrapping; A ≠ S |
| `literature-officer-20260905-newman-annulus.md` | four `R_*`; complex annulus roots |
| `literature-officer-20260905-mertens-matching.md` | finite-size matching identity |

---

## Not established, and not claimed here

- matching-odd = Vasseur `F(r)`, Ridout `t`, Picco `R_σ`, Pinson `π({1,0})`, `R_h`, or Duncan giant cycle;
- `π({1,0}) = R_h`;
- a square-site algebraic `p_c`;
- rank-4 Jordan at `x=21/4`;
- `θ(p_c)=0` in 3-d as a refereed theorem;
- anything in the claim ledger.
