# Literature officer brief — 2026-09-05

**Role:** theory input. Same epistemic level as an Astra answer or a
derivation we have not checked. **Does not enter** `docs/STATUS.md`.
**Does not score** any frozen block.
**Companion primary reading:** `notes/ziff-2006-a-lattice-primary-reading.md`
(closes the scientific half of #566).

This is the first pass that went outside the Ziff/Scullard/Jacobsen citation
island. Everything below is here because it changes a ticket, a competitor
list, or a sentence in the P2 draft. Papers that are merely "related to
percolation" are omitted.

## Ticket map

| Ticket / file | What the literature does |
|---|---|
| #566, P2 §1.1 / §4.2 | A-lattice quintic confirmed from Ziff PRE 73, 016134. Height 4 stands. |
| #565, Astra Q4 | Degree/height regularity is **not** a theorem of the three mechanisms. Cell size is unbounded; A-lattice already broke height 3. Square site is named by Ziff as outside the 2006 method. |
| #567, ROADMAP item 2 | Do not score the ladder against `E₄`. Continuum wrapping is Pinson 1994 (`η`, `Z_{m,n}`), and crossing densities are **second-order** modular forms on `Γ(2)`, not weight-4 holomorphic forms. Name that as a competitor or a non-claim before the run, not after. |
| Q2 | Additive-shape ambiguity has a literature name: Kleban–Zagier / Diamantis–Kleban. The `11/4` fingerprint is a conjunction with a normalization the literature does not force. |
| Q1, #275 | `μ = -5/4` at `h = 5/8` is He 2024. Watts' crossing needs extra Jordan modules (Ridout 2009). Lattice log observables exist (Vasseur–Jacobsen–Saleur; Hu–Deng `C_n ~ L^{2y_t-d} log L`). None of this identifies the lattice overlap. |
| matching-odd `x = 21/4` | Polychromatic 8-arm / 4-cluster watermelon, `(ℓ²-1)/12`. Monochromatic 8-arm has **no closed form**. Backbone (mono 2-arm) is transcendental. Colour-decompose before fitting 21/4. |
| matching pair itself | Grimmett–Li 2024: `p_c(G_*) < p_c(G)` iff `G` is not a triangulation. Square is not. |
| P2 novelty | No published PSLQ/Sturm census of square-site `p_c`. Short-range connectivities *can* be rational (`g_n = 3/4`, `g_{nn} = 11/16`) while `p_c` stays unknown. |

X / community: Deng's group is posting 2026 long-range percolation; nobody
is discussing algebraic square-site `p_c` or matching-odd. The Quanta
flooding proof is orthogonal. `@ZiffRziff` is not Robert M. Ziff.

---

## 1. Q4 / P2 — the height bound is a small-sample accident

**Ziff, PRE 73, 016134 (2006)** = arXiv:cond-mat/0510245. Primary reading in
the companion note. The A-lattice bond threshold is the unique root in
`(0,1)` of `p⁵ − 4p⁴ + 3p³ + 2p² − 1 = 0` (height 4). Same paper:

> the method does not appear to work for … site percolation on the square
> and honeycomb lattices, and bond percolation on the kagomé lattice.

**Scullard, PRE 73, 016107 (2006):** martini site `p⁴ − 3p³ + 1 = 0`
(height 3). Confirmed.

**What Q4 should not assume.** The three mechanisms (self-dual, self-matching,
star-triangle / cell-dual-cell) produce a polynomial whose degree and height
grow with the cell. There is no uniform bound independent of `n`. A-lattice
was the seventh sample and broke height 3; a larger cell will break 4. The
census class `C(≤6, ≤4)` is a **historical snapshot**, not a theorem about
the mechanisms. Q4's first branch ("the mechanisms bound degree and height")
is the one the record already weakens. The useful remaining half is: *is
square site provably outside every reachable cell of those mechanisms, not
just the 2006 ones?* Ziff 2006 answers 2006. It does not answer "every".

**Nolin–Qian–Sun–Zhuang, arXiv:2309.05050.** The backbone exponent (monochromatic
two-arm) is the unique root in `(1/4, 2/3)` of an equation involving `sin`,
numerically `0.356666…`, and is transcendental. So "exact" 2d percolation
quantities need not be algebraic. P2 discussion can say this without claiming
transcendence of `p_c`.

**Grimmett–Li, Random Structures & Algorithms 65 (2024) 832,**
arXiv:2205.02734. Matching lattice-pairs: for transitive plane `G`,
`p_c(G_*) < p_c(G)` iff `G` is not a triangulation. Square is not. This is
the theorem behind `p_c(Sq8) = 1 − p_c(Sq)`, not just Sykes–Essam 1964.
van den Berg 1981 already had a mosaic counterexample; Grimmett–Li is the
modern if-and-only-if.

P2 §8 can cite Grimmett–Li and Nolin–Sun. It should not claim the census
class is forced by the mechanisms.

---

## 2. Q2 / #567 — `E₄` is the wrong modular object

The N=290 ratio `A(2i)/A(i) = 1.880 ± 0.177` excluded area-normalized
weight-4 `E₄` at 4.9σ. The literature on what torus/crossing observables
actually are:

| Paper | What it is | Bearing |
|---|---|---|
| Pinson, J. Stat. Phys. **75**, 1167 (1994) | Exact wrapping / winding probabilities on the torus, combinations of Dedekind `η` and `Z_{m,n}(g; r)` at Coulomb-gas `g = 2/3` | This is the continuum prediction for wrapping as a function of aspect ratio `r`. Our fingerprint is a finite-size version of a Pinson observable, not of `E₄(τ)` |
| Pruessner–Moloney, cond-mat/0310361 | Numerical confirmation of Pinson windings | Exists; we do not need to redo the continuum curve |
| Kleban, cond-mat/9911070; Kleban–Zagier | Derivatives of `Π_h`, `Π_{hv}` form a **vector modular form**: one ordinary component, one "form × integral of a form" | Second-order, not holomorphic weight 4 |
| Diamantis–Kleban, arXiv:0905.1727, CNTP **3** (2009) 677 | Simmons–Kleban–Ziff crossing densities are weakly holomorphic **second-order modular forms of weight 0 on `Γ(2)`**, type involving `η⁴` | `E₄` is weight 4 on `SL(2,ℤ)`. These are not that |
| Mertens–Ziff, PRE **94**, 062152 (2016), arXiv:1603.07289 | Exact finite-torus identity relating cluster numbers to wrapping via the matching polynomial `χ_□ = p − 2p² + p⁴` | Already in the repo. Gives the *difference* `R − R̂`, not the Pinson curve of Sq8 itself |

**For #567.** The frozen competitors are weight-4, bare aspect ratio, area,
none. Literature adds a fifth *class*, not a sixth number: Pinson wrapping
`R(r)` is a function of `r` through `η` and theta, and is not proportional
to `E₄`. If the matching-odd slope is a wrapping-type observable, the
ladder should be compared to Pinson (or declared not to be). Either
sentence is a result. Scoring only against `E₄` and then reading `r` off
the leftover is how the post-hoc `bare_aspect_ratio` happened at N=290.

No published Pinson table for the matching lattice (NN+NNN / Sq8) was
found. Universality says the *critical wrapping probabilities* are the
same; the matching identity constrains the difference. That is the
zero-new-compute check: square wrapping vs Pinson, matching wrapping vs
the same Pinson, difference vs Mertens–Ziff.

---

## 3. Q1 / `x = 21/4` — two different 5.25s, and colour

Coulomb-gas watermelon / polychromatic arm exponents
(Saleur–Duplantier; Aizenman–Duplantier–Aharony, PRL **83**, 1359 (1999)):

```text
x_ℓ^P = (ℓ² - 1) / 12
```

| ℓ | x | geometry |
|---:|---|---|
| 4 | 5/4 | four-arm / energy / pivotal |
| 8 | **21/4** | 8-leg watermelon = 4 disjoint clusters |

The repository `x = 21/4` candidate is the 8-leg (even, four-cluster)
exponent, not a new number. Four-arm `5/4` is the energy/Jordan layer
He 2024 pins with `γ^perco = μ = -5/4`.

**Beffara–Nolin, Ann. Probab. 39 (2011):** monochromatic `j`-arm exponents
satisfy `α_j < α'_j < α_{j+1}` and have **no closed form** for `j ≥ 2`
except `j = 2` (backbone), now known to be transcendental. If
matching-odd is same-colour wrapping, `21/4` is the wrong formula.

**Reeves–Sosoe, arXiv:2009.07029.** Polychromatic arm probabilities on the
**square** lattice are equivalent up to constants (colour-switching works
despite the dual lattice). Square site may use the polychromatic formula
without hiding on the triangle. That does not license skipping the colour
test.

**He, arXiv:2411.18696,** SciPost Phys. 19, 008: `μ = -5/4` at `h = 5/8`.
This is Q1's primary-level input. Descent to the level-4 `x = 21/4, s = 4`
state is still the question; the literature does not compute that coupling.

**Ridout, Nucl. Phys. B 810 (2009) 503.** Watts' both-direction crossing
uses a primary that does not sit in the usual extended Kac table; extra
rank-2 modules appear. Q1's Jordan story is richer in the Watts channel
than in Cardy's.

**Hu–Blöte–Ziff–Deng, PRE 90, 042106 (2014),** arXiv:1406.0130. Square
**bond** critical connectivities: `g_n = 3/4` exact, `g_{nn} = 11/16`
(Mitra–Nienhuis conjecture, MC to `0.6875000(2)`), surface `g_n = 5/8`.
Next-nearest neighbour is a matching-lattice edge. Algebraic numbers in
this model live in short-range connectivities; they need not live in
`p_c`. Fluctuations `C_n ~ L^{2 y_t - d} log(L/L_0)` are the Vasseur
logarithmic observable on the lattice. If matching-odd amplitudes spit
out a simple fraction, check `3/4`, `11/16`, `5/8` before inventing one.

**Feng–Deng–Blöte, PRE 78, 031136 (2008).** Square NN+NNN (our matching
lattice) has Coulomb-gas `X_{t2} = 4` corrections **and** a logarithmic
factor; triangular site does not. Matching-lattice FSS is dirtier than
square site. Do not fit matching-odd as a pure power.

---

## 4. What P2 can say that the literature does not already say

Nobody has published a certified exhaustive exclusion of low-degree
integer polynomials against the square-site intervals. PSLQ in the
experimental-mathematics literature is used for `π`, MZVs, QFT periods,
not for this constant. That is the paper's actual novelty, together with
the interval-disjointness and the sensitivity-certified nulls.

Do not sell "the historical complexity range" as a law. Sell: every
exactly-known planar threshold we can source is in `C(≤6, ≤4)`; square
site is not in that class on any of the four published intervals; the
2006 method that produced those thresholds does not apply to square
site, by the method's author.

---

## 5. Open literature holes (worth a later pass, not this PR)

- Pinson wrapping curve evaluated on Sq8 (NN+NNN). Not found.
- Independent 8-leg watermelon MC on square site. Not found; triangle and
  Coulomb gas only.
- A primary reading of Wierman bow-tie / larger-cell polynomials that
  would push height past 4. Secondary indexes already suggest height 6;
  that would move F′ the way A-lattice moved F, and is the reason Q4
  should not freeze a bound.
- Jacobsen 2015 TL-eigenvalue paper is already in the interval table;
  no new algebraic form for square site there either
  (`p_c = 0.59274605079210(2)`).

## Not established, and not claimed here

- identification of matching-odd with the 8-leg watermelon;
- that the N=290 ratio is a Pinson number;
- a bound on degree or height for all future exact thresholds;
- transcendence of square-site `p_c`;
- anything that belongs in the claim ledger.
