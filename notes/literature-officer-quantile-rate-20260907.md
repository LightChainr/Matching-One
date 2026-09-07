# Literature officer — location-to-rate: inverse-CDF windows, DKS descendants, F1 connecting lemma

**Date:** 2026-09-07
**Ticket:** [#620](https://github.com/LightChainr/Matching-One/issues/620). Brief: `notes/probes/probe-LIT-retrieval-quantile-rate-20260907.md` on `docs/probes-full-law-20260907` (PR #616).
**Role:** theory input. Same epistemic level as #601 / #592 / the #613 note. **Does not enter** `docs/STATUS.md`. **Does not close** #276, #321, #613, #566, #592, #601, #606, #615. **Does not mix** with PR #602, #604, #606, #607, #615, or with the three computational probes (#617 / #618 / #619).
**Access:** every theorem below is a primary quote with PDF-stable coordinates, or is tagged **NOT FOUND**. Paywalled PDFs that were not retrieved are flagged as such; the published abstract is then quoted and is *not* a substitute for a theorem-body quote.
**Lattice discipline:** triangular-site theorems are not square-site theorems. Gate 1 “five widths” and #582 “five transitions” are not literature.

The rate probe (MID) is forbidden from searching. This note is the search.

---

## Decision table

| Q | Tag | One-line |
|---|---|---|
| **Q1** | **FOUND-BUT-NOT-INVERSE-CDF** (window existence) / **NOT FOUND** (percolation window, or inverse-CDF of `r_G`) | Boolean-analysis and DKS give *some* window. DKS’s `λ□(N)` is the single `u=1/2` quantile of `{ϕ_* ≠ 0}`. Nothing names nine interior quantiles of `r_G`, and nothing gives `N^{-3/4}`. |
| **Q2** | **FOUND** (location + exponential for `i=1`) / **NOT FOUND** (power window; square-site+matching in print) | DKS Thm 1–3: sharp = `λ□(N)±ε` of *fixed* `ε`, not `L^{-θ}`. Exponential decay of `P(A□)` below threshold for `i=1` and `i=d−1`. `d=2,i=1` plaquette = square **bond**. Site model = permutohedral = triangular in `d=2`. |
| **Q3** | **FOUND** (FK / Bourgain / Hatami / DKS) / **NOT FOUND** (sharp-threshold ⇒ inverse-CDF window lemma; never `L^{-1/ν}` without RSW/arms) | FK96 window is `c log(1/2ε)/log n`. Bourgain (1.2) is `τ/p → 0`. Hatami generalizes Friedgut; does not shrink the window. |
| **Q4** | see table in Q4 | Square-site **RSW is FOUND** (Zeng Thm 1.1 unpublished; KST 2023 Comment 1). Square-site **conformal invariance remains OPEN**. These two are not the same statement; this is **not** in tension with the #606 ledger. Inverse-CDF of crossings: **NOT FOUND** on every lattice. |
| **Q5** | **FOUND** triangular (LSW / SW) / **NOT FOUND** square bond/site four-arm; **NOT FOUND** connecting lemma | `α_4 = −5/4`, `ν = 4/3` on triangular site. Connecting lemma “four-arm ⇒ `M' ≍ L^{3/4}` on a torus”: **NOT FOUND**. |
| **Q6** | **FOUND** Cardy *limit* / **NOT FOUND** finite-`L` expansion as a theorem | Simulations are not theorems. |
| **Q7** | **FOUND** | (1.1) is motivation, not a theorem of GL. (1.3) is the theorem, in the hyperbolic companion. Amenability + uniqueness ⇒ (H3). |
| **Q8** | **FOUND** (DCT site paragraph) / AB87 PDF not retrieved; published abstract covers the claim | DCT Thm 1.1(3) + §1.2 “Site percolation. As in [AB87]…”. AB87 abstract: all translation-invariant independent percolation on homogeneous lattices. |
| **Q9** | **FOUND** location-type descendants / **NOT FOUND** rates, matching, quantiles | Voronoi homological 2601.00793 is still `p=1/2` on a `2i`-torus. Bobrowski–Skraba thermodynamic-limit location. |
| **Q10** | **NOT FOUND** as a named lemma | Elementary tautology if the window is already uniform on `[ε,1−ε]`. Single-`u` sharpness does not give nine quantiles. |
| **Q11** | **FOUND** Cardy limit / **NOT FOUND** finite-`L` expansion of `Q_L(u)−1/2` | Self-duality locates `p_L^H = 1/2`; it does not give a shape rate. |
| **Q12** | **FOUND** GPS on the triangular *plane* / **NOT FOUND** torus wrapping conversion | `E[# pivotals in Q] ≍_Q η^{-2} α_4^η(η,1) = η^{-3/4+o(1)}` on `ηT`. Extra for a torus: periodic images, wrapping definition of pivotals. |
| **Q13** | **NOT FOUND** | No paper claims “because `{r_G>0}` is sharp, `Q_N(u)−p_c = O(N^{-α})` with percolation `α` from boolean influences”. |
| **Q14** | unused | Primary papers were not thin. |

---

## Reuse table (already quoted in `notes/p613-quantile-convergence-20260907.md`)

Copy once, then stop re-summarising. New retrieval starts after this table.

| Source | Coordinates | Quote already used for #613 | What it does *not* contain |
|---|---|---|---|
| DCT, CMP **343** (2016) 725–745; arXiv:1502.03050v3 | Thm 1.1 item 3, p. 2; §1.2 “Site percolation”, p. 3 | “If `(J_{x,y})` is finite-range, then for any `β<β_c`, there exists `c=c(β)>0` such that `P_β[0↔Λ_n^c] ≤ e^{-cn}` for all `n≥0`.” / “As in [AB87], the proof may be adapted to site percolation on transitive graphs.” | No rate at `β_c`. No inverse-CDF. No torus. No four-arm. |
| Grimmett–Li, RSA **65** (2024); arXiv:2205.02734v3 | p. 2–3, displayed (1.1)–(1.3) | (1.1) `p_c^{site}(G)+p_c^{site}(G*)=1` is Sykes–Essam *motivation*. (1.3) `p_u^{site}(G)+p_c^{site}(G*)=1` is the theorem (proved in the companion). “When `G` is amenable, we have `p_c^{site}(G)=p_u^{site}(G)`, in agreement with (1.1).” | No finite-size window. No quantiles. |
| Grimmett–Li hyperbolic, arXiv:2203.00981v3 | Thm 1.1(a); Remark 1.4 | `p_u^{site}(G1)+p_c^{site}(G2)=1` for a matching pair. Amenable ⇒ `p_c=p_u` by uniqueness of the infinite cluster. | Same. |
| van den Berg, JMP **22** (1981) 152–157, DOI 10.1063/1.524747 | published abstract (PDF paywalled; see Q7) | Sykes–Essam “suggested that the above relation holds for all mosaics … **we have constructed a counterexample**”; derivation for “a more restricted class of graphs … based on the usual assumption that below the critical probability the mean cluster size is finite”. | No window, no torus, no inverse-CDF. |
| DKS, arXiv:2011.11903v4; AIHP 2025 | abstract + Thm 1–3 | Giant cycles (`im H_i → H_i(T^d) ≠ 0`) have a sharp threshold. `d=2,i=1` plaquette = square **bond** on the torus. Site model is permutohedral (triangular in `d=2`). | No `L^{-θ}` window. No square-site+matching. No multi-`u` inverse-CDF of rank. |
| Burton–Keane, CMP **121** (1989) 501–505 | cited by GL Remark 1.4 as uniqueness; PDF not retrieved this pass | Uniqueness of the infinite cluster under stationarity + finite energy (standard input to `p_c=p_u` on amenable lattices). | No rate. |
| Smirnov / LSW as named in the #606 ledger | pinned in Q5 | Triangular site only. | Square site / square bond four-arm **open**. |

Theorem L itself (`Q_N(u)→p_c` uniformly on compact subsets of `(0,1)`, no rate) is a **repository lemma**, not literature. It is not upgraded by quotation density.

---

## Method (queries)

Queries actually run (arXiv / Crossref / AMS / AIP / Springer). Negative results are recorded at the question they belong to.

1. `homological percolation torus sharp threshold` → DKS 2011.11903, Bobrowski–Skraba 2005.14011, Schweinhart–Shuman 2601.00793.
2. `inverse CDF quantile finite-size scaling percolation` / `p_c(L) - p_c four-arm` → Kesten scaling (via SW), BCKS (box, `|C_max|`), GPS pivotals. No inverse-CDF of `r_G`.
3. `Friedgut Kalai sharp threshold` → Proc. AMS 124 (1996) 2993–3002, DOI 10.1090/S0002-9939-96-03732-X (PDF retrieved). AMS “S0002-9939-96-03515-2” is the wrong identifier.
4. `Bourgain sharp threshold appendix Friedgut JAMS 1999` → Friedgut JAMS 12 (1999) 1017–1054 + Bourgain appendix (PDF retrieved).
5. `Hatami small total influences` → arXiv:1008.1021v3.
6. `RSW site percolation Z^2` → Zeng arXiv:1309.2273; Köhler-Schindler–Tassion arXiv:2011.04618 / Duke Math. J. 172 (2023).
7. `four arm exponent percolation triangular` → LSW math/0108211, Smirnov–Werner math/0109120.
8. `Cardy formula finite size correction theorem` → Smirnov 0909.4499 (limit only).
9. `Grimmett Li matching 2205.02734` / `hyperbolic site 2203.00981` → PDFs retrieved.
10. `van den Berg matching lattices 1981` → published abstract (AIP 403 on PDF).
11. `Aizenman Barsky CMP 108` → published abstract (Springer PDF paywalled).
12. `Garban Pete Schramm pivotal 1008.1378` → PDF retrieved.
13. `Borgs Chayes Kesten Spencer birth infinite cluster` → published abstract CMP 224 (2001); arXiv IDs `math/9811085`, `math/9912001`, `math-ph/0103012` are **wrong papers**. Full BCKS PDF not retrieved.
14. `sharp threshold inverse CDF` / `quantile window boolean function percolation` → **NOT FOUND**.

---

## Q1 — Inverse-CDF / quantile finite-size scaling

**Tag:** **FOUND-BUT-NOT-INVERSE-CDF** for the existence of *some* window; **NOT FOUND** for a theorem that an inverse of a sharp-threshold family

```text
Q_N(u) := inf{ p : P_N(event_p) ≥ u }
```

inherits a *percolation* window (`N^{-3/4}` / `L^{-1/ν}`) for crossing events or for `{r_G ≥ 1}` / `{r_G = 2}` on a torus; **NOT FOUND** for any theorem about nine interior quantiles of `r_G`.

### What exists (single-quantile / some window)

**DKS, arXiv:2011.11903v4, p. 5–6, displayed (1) and Theorem 3.** Observable: the event `A□` that `ϕ_*: H_i(P) → H_i(T_N^d)` is nontrivial. They *define* a single quantile

> For each `N ∈ ℕ`, let `λ□(N,i,d)` satisfy `P_{λ□(N,i,d)}(A□) = 1/2`.

and immediately warn:

> We should mention that this choice of `λ□(N,i,d)` is somewhat arbitrary. We could replace `1/2` in Equation 1 with any constant strictly between 0 and 1, for example, and the sharp threshold results we use would apply just as well.

Theorem 3 then says: for every `ε>0`, `P_{λ□(N)−ε}(A□) → 0` and `P_{λ□(N)+ε}(S□) → 1`. That is a *fixed-ε* window about a *single* (or, by the remark, any one) quantile. It is not a power of `N`, and it is not a simultaneous statement about a vector of interior `u`.

**Friedgut–Kalai, Proc. AMS 124 (1996) 2993–3002, p. 2993, Theorem (and Thm 2.1 p. 2994).** Observable: a symmetric monotone subset `A ⊂ {0,1}^n`.

> **Theorem.** For every symmetric monotone `A`, if `μ_p(A)>ε` then `μ_q(A)>1−ε` for `q = p + c_1 log(1/2ε)/log n`. (`c_1` is an absolute constant.)

Window width `O(log(1/ε)/log n)`, never `L^{-3/4}`. Because of the `log(1/ε)`, the *same* argument is uniform on compact subsets of `(0,1)`, at the price of a larger implied constant. This is a window for the *event*, not a named theorem about the inverse-CDF of homological rank.

**Bourgain appendix to Friedgut, JAMS 12 (1999), p. 1017–1018 (definitions) and appendix (1.2).** Sharp vs coarse:

> one observes in many cases a “sharp” threshold phenomenon, in the sense that `μ_p(A)` jumps from near 0 to near 1 in an interval `τ = τ(n)` which is small with respect to the threshold `p = p(n)` when `n→∞`, thus `τ(n)/p(n) → 0`. (1.2)

The random-graph giant-component example they give is `p = 1/m`, `τ ∼ 1/m^{4/3}` — a power, but for *Erdős–Rényi*, not for torus wrapping, and not derived from four-arm inputs.

**Hatami, arXiv:1008.1021v3, abstract.** Structure of Boolean functions with small total influence; “generalizes the core of Friedgut’s seminal work … and improves the result of Bourgain in his appendix.” Does not produce a smaller window, and does not mention inverse-CDFs.

### What does not exist

- A theorem whose *stated observable* is `Q_N(u) := inf{p: P_N(r_G ≥ 1) ≥ u}` or `F_N^{-1}` with `F_N = E[r_G]/2`, for more than one `u`.
- A theorem that a boolean-analysis window of width `w_N` upgrades, without RSW/arm inputs, to the percolation window `L^{-1/ν}`.
- BCKS (published abstract, CMP **224** (2001) 153–204; **PDF not retrieved**; arXiv IDs tried were the wrong papers) studies **bond** percolation in a **box** of side `n`, observable the size of the largest clusters, and “the scaling window in which the system behaves critically”. That is finite-size scaling for `|C_max|`, not the inverse-CDF of homological rank on a torus.

**Rate probe may write.** **T:** “No named theorem gives a percolation-scale window for the inverse-CDF of `r_G`.” **C:** “DKS `λ□(N)` is a `u=1/2` quantile of `{ϕ_* ≠ 0}` on the *plaquette* (square-bond) torus, with a fixed-`ε` rather than `L^{-θ}` window; FK96 gives a `1/log n` window for symmetric monotone events.” **Blocked:** any `O(L^{-3/4})` or “nine quantiles share the four-arm window”.

---

## Q2 — DKS and descendants: rates, not just location

**Tag:** **FOUND** for location-sharpness and for exponential decay away from `p_c` when `i=1` or `i=d−1`. **NOT FOUND** for a `L^{-θ}` rate for `P(im H_i ≠ 0)`. **FOUND-BUT-WRONG-LATTICE** for square-site + matching (not in print).

### What “sharp” means in DKS

**DKS arXiv:2011.11903v4, abstract, p. 1:**

> We show that for every `i` and `d` there is a sharp transition from nonexistence of giant cycles to giant cycles spanning the homology of the torus. … we prove that `p_c = 1/2` in the case of middle dimension `i = d/2` for both models. This gives finite-volume high-dimensional analogues of Kesten’s theorems that `p_c = 1/2` for bond percolation on a square lattice and site percolation on a triangular lattice.

**Theorem 1 (p. 5).** `d = 2i`, `char F ≠ 2`: `P_p(A□) → 0` for `p < 1/2` and `P_p(S□) → 1` for `p > 1/2`.

**Theorem 2 (p. 5).** `i = 1`: the same with threshold `p̂_c` = bond-percolation critical point of `ℤ^d`. `i = d−1`: threshold `1 − p̂_c`. And (same page):

> In the above, we also show that the decay of `P_p(A□)` below the threshold and `P_p(S□)` above the threshold is exponentially fast for both `i = 1` and `i = d−1`.

**Theorem 3 (p. 6).** For every `ε > 0`, `P_{λ□(N)−ε}(A□) → 0` and `P_{λ□(N)+ε}(S□) → 1`. The window is an *additive constant `ε` independent of `N`*, not a power of `N` (or of `L`).

Lattice identification, from the introduction (p. 3–4): the plaquette model on `T_N^d = ℤ^d/(Nℤ)^d` with `i=1` “is” bond percolation on the cubical torus. In `d=2` that is square-**bond**. The site model is permutohedral tiling; in `d=2` that is the hexagonal/triangular duality, i.e. triangular **site**, not square site.

### Square-site + matching in print

DKS do not write the NN / NN+NNN matching pair on the square torus. Theorem L’s square-site lift is the repository’s (H3) plus DKS’s union-bound mechanism. **Do not claim novelty for that lift; do record that nobody has written the site+matching version in print.** Tag: **NOT FOUND** in print.

### Descendants with a *power* for homological percolation

**Schweinhart–Shuman, arXiv:2601.00793v2, abstract / p. 1.** Voronoi percolation on a torus:

> As a consequence, we prove a sharp phase transition for the emergence of `i`-dimensional giant cycles in Voronoi percolation on the `2i`-dimensional torus.

Location-type, threshold `1/2` on a self-dual `2i`-torus. No `L^{-θ}` window, no quantiles of rank.

**Bobrowski–Skraba, arXiv:2005.14011v1, abstract.** Continuum percolation on the flat torus, thermodynamic limit `n r^d = λ`: giant `k`-cycles appear (and, for `k=1`, the trivial-to-surjective transition is identified). Exponential decay of probabilities *outside* the thermodynamic window. Observable is still `{Im H_k ≠ 0}`, not an inverse-CDF, and the scaling is continuum-intensity not `L^{-1/ν}`.

**Rate probe may write.** **T:** “DKS prove location-sharpness of giant cycles; for `i=1` the subcritical decay of `P(A□)` is exponential in `N`. They do not prove a power window.” **C:** “`d=2,i=1` plaquette is square bond; permutohedral site is triangular. Square-site+matching is not in DKS.” **Blocked:** “DKS give `L^{-3/4}`”; “DKS already did square site”.

---

## Q3 — Window-width theorems vs percolation window

**Tag:** **FOUND** for each named boolean-analysis theorem (with the window they actually give). **NOT FOUND** for “sharp threshold ⇒ inverse-CDF window” as a lemma. None of these windows can be as small as `L^{-1/ν}` without RSW/arm inputs.

| theorem | window actually given | hypothesis class | can it be `L^{-1/ν}` without RSW/arms? |
|---|---|---|---|
| Friedgut–Kalai, Proc. AMS 124 (1996), Thm p. 2993 / Thm 2.1 p. 2994 | `q = p + c log(1/2ε)/log n` | symmetric monotone subsets of `{0,1}^n` (transitive group) | **No.** `1/log n`. For percolation `n = Θ(L^2)` this is `1/log L`. |
| Friedgut, JAMS 12 (1999) | sharp ⇔ `τ(n)/p(n) → 0`; coarse otherwise | monotone graph properties; Bourgain appendix for general monotone | **No.** The *definition* of sharp is `o(p)`, not a percolation power. Their ER-giant example `τ ∼ m^{-4/3}` uses graph-specific inputs. |
| Bourgain appendix (Friedgut JAMS; displayed (1.2) p. 1018) | `τ(n)/p(n) → 0` | monotone Boolean | **No.** |
| Hatami, arXiv:1008.1021 | structure of small-total-influence Boolean functions; does not claim a smaller window | product spaces, not necessarily monotone | **No.** |
| DKS Thm 3 | `λ□(N) ± ε` for any fixed `ε>0` | homological events on the torus, using Friedgut–Kalai-type sharpness (their Thm 6) | **No.** Fixed `ε`, not a power. |
| “sharp threshold ⇒ inverse-CDF window” lemma | — | — | **NOT FOUND.** |

Without RSW (to get a uniformly positive crossing probability at `p_c`) and four-arm / Kesten scaling (to convert Russo’s formula into `ν(2−ξ_4)=1`), boolean analysis cannot see `L^{-3/4}`. Duminil-Copin, *Sixty years of percolation*, arXiv:1712.04651, records the classical RSW ⇒ Kesten `p_c=1/2` route, not a boolean-analysis route to `ν`.

**Rate probe may write.** **T:** “FK / Bourgain / Hatami / DKS-Thm-3 windows are `1/log n` or fixed-`ε`. There is no lemma that upgrades sharpness to an inverse-CDF of width `L^{-1/ν}`.” **Blocked:** any citation of FK96 as F1.

---

## Q4 — RSW / box-crossing: who has what lattice

**Tag:** mixed; see table. Square-site RSW is **FOUND**. Square-site conformal invariance is still **OPEN**. Inverse-CDF of crossings is **NOT FOUND** on every row.

| lattice / model | RSW? | source (Thm, page) | implies what for crossing *probabilities* | implies what for *inverse-CDF* |
|---|---|---|---|---|
| **triangular site** | **FOUND** | Smirnov, arXiv:0909.4499, abstract: conformal invariance of crossing probabilities and Cardy’s formula, critical site percolation on the triangular lattice. RSW is classical for this lattice (self-matching); used as input by LSW math/0108211 p. 2. | Crossing probabilities of conformal rectangles have a conformally invariant *scaling limit*. At finite `L`, RSW gives `c(λ) ≤ P_{p_c}(horizontal crossing of n × λn) ≤ 1−c(λ)`. | **NOT FOUND** as a theorem for `Q_L(u)−1/2`. Cardy is a limit, not a finite-`L` expansion (Q6, Q11). |
| **square bond** | **FOUND** | Russo, *Z. Wahrsch.* 43 (1978) 39–48; Seymour–Welsh, *Z. Wahrsch.* 47 (1978) 215–221 — as recorded by Duminil-Copin, arXiv:1712.04651 p. 2: “an important result obtained simultaneously by Russo [74] and Seymour-Welsh [77]”; Kesten 1980 uses it for `p_c=1/2`. KST arXiv:2011.04618 p. 3: “In 1978, the first result of this kind was proven in the case of Bernoulli percolation by Russo [Rus78, Rus81] and independently by Seymour and Welsh [SW78].” **Russo / SW PDFs not retrieved this pass; the existence of the theorem is recorded by the PDFs we do have.** | Box-crossing property at `p=1/2`: probabilities bounded in `(c,1−c)` depending only on aspect ratio. | **NOT FOUND** for inverse-CDF. Self-duality gives the *median* crossing quantile exactly `1/2` (see Q11). |
| **square site** | **FOUND** (RSW) | **Zeng, arXiv:1309.2273v1, Theorem 1.1 p. 1** (unpublished; the PDF calls itself “the current dissertation”): “For any `λ>0`, there exists `c=c(λ)>0`, such that for all `n≥1`, `c ≤ P_{p_c}[there exists a horizontal crossing of any n by λn box] ≤ 1−c`, where `P_{p_c}` stands for critical site percolation on `ℤ^2`.” Uses Kesten box-crossing and the matching pair `(ℤ^2, ℤ^{2,*})`; Remark 2.5 cites `p_c(ℤ^2)+p_c(ℤ^{2,*})=1` from Kesten’s book. **Köhler-Schindler–Tassion, arXiv:2011.04618v1, Theorem 1 p. 2**, published *Duke Math. J.* **172** (2023) 809–864: for every `ρ≥1` there is a homeomorphism `ψ_ρ` such that for every invariant positively associated measure `P` on bond configurations of `ℤ^2` and all `n≥1`, `P[C(ρn,n)] ≥ ψ_ρ(P[C(n,ρn)])`. **Comment 1, p. 2:** “While we show the theorem and its proof in the framework of bond percolation on `ℤ^2` for presentational purposes, we would like to point out that **the same proof applies to more general lattices … and more general processes (including site percolation, continuum models, or level lines of random fields).**” | Critical square-site box-crossing probabilities stay in `(c,1−c)`. This is RSW, **not** conformal invariance, **not** Cardy’s formula, **not** SLE. | **NOT FOUND.** |
| **FK-Ising / critical Ising** | **FOUND** (RSW + conformal covariance of fermionic observables) | Chelkak–Smirnov, arXiv:0910.2045v2, abstract: discrete holomorphic fermions for 2D Ising / FK-Ising at criticality on isoradial graphs; “universal and conformally invariant scaling limits”. Duminil-Copin–Tassion, arXiv:1901.08294, abstract: renormalization of crossing probabilities for planar random-cluster, including a “critical continuous” alternative in which crossing probabilities stay bounded away from 0 and 1. DHN arXiv:0912.4253 is cited by KST as the FK RSW; **DHN PDF text extraction was garbled (Type-3 fonts); do not quote the body.** | Uniform (in scale) bounds on crossing probabilities at criticality, plus conformal invariance of *observables*, not of square-site percolation. | **NOT FOUND.** |
| **homological torus events** | **FOUND** location, **NOT FOUND** RSW-window | DKS Thm 1–3 (Q2). Tassion, arXiv:1410.6773v2, Theorem (abstract): “the standard Russo-Seymour-Welsh theory is valid for Voronoi percolation. This implies that at criticality the crossing probabilities for rectangles are bounded by constants depending only on their aspect ratio.” Voronoi is not the square-site torus. | DKS: `P(A□)→0` or `1` off `p_c`, exponentially for `i=1`. No `(c,1−c)` bound *at* `p_c` for wrapping events is stated as RSW. | **NOT FOUND.** |

### Square-site RSW vs the #606 ledger

The #606 ledger’s “square-site conformal invariance open” is about the **scaling limit / Cardy / SLE_6**. Zeng Thm 1.1 and KST Comment 1 are about **box-crossing probabilities bounded in `(c,1−c)`**. RSW is a strictly weaker statement than conformal invariance; Smirnov’s theorem is the one that is open on square site. **There is no tension.** Do not close #276 or #321. A pointer comment quoting Zeng Thm 1.1 is the brief’s only permitted issue-comment.

**Rate probe may write.** **T:** “Square-site RSW is in the literature (Zeng Thm 1.1; KST Comment 1). Square-site conformal invariance is not.” **C:** “RSW gives uniform-in-scale crossing bounds at `p_c`; it does not give `Q_N(u)−p_c = O(L^{-θ})`.” **Blocked:** treating Zeng as Cardy; treating triangular Cardy as square-site Cardy.

---

## Q5 — Four-arm `α_4 = −5/4` and `ν = 4/3`

**Tag:** **FOUND** on triangular site. **NOT FOUND** (open) on square bond and square site. **NOT FOUND** for the connecting lemma four-arm ⇒ `M' ≍ L^{3/4}` on a torus.

### Triangular site (proved)

**Lawler–Schramm–Werner, arXiv:math/0108211v1, Theorem 1.1 p. 1.** One-arm, critical site percolation on the triangular grid:

> `P[0 ↔ C_R] = R^{-5/48+o(1)}, R → ∞`.

And p. 2:

> It is conjectured that the theorem holds for any planar lattice. However, Smirnov’s results mentioned above have been established only for site percolation on the triangular lattice, and hence we can only prove our result in this case.

**Smirnov–Werner, arXiv:math/0109120v2, Theorem 1 p. 3**, triangular site, `p → 1/2`:

> (i) `θ(p) = (p−1/2)^{5/36+o(1)}` as `p→1/2+`.
> (iii) `ξ(p) = (p−1/2)^{-4/3+o(1)}`.
> (iv) the same for `ξ*`.

provided, at `p=1/2`, `P[A^1_R] = R^{-5/48+o(1)}` and `P[A^2_R] = R^{-5/4+o(1)}`. The four-arm (polychromatic two-arm in their numbering of alternating blue/yellow) exponent `5/4` is the `α_4` of the ledger; `ν = 4/3` is (iii). The paper attributes the scaling relations to **Kesten 1987** (their [13]); **Kesten CMP 109 PDF was not retrieved this pass.** The SW theorem is the combination of Kesten’s relations + LSW + Smirnov.

### Square bond / square site

LSW p. 2, quoted above: **only triangular site**. DMT, arXiv:2011.15090, records Kesten’s relation `ν(2−ξ_4)=1` for random-cluster `q=1` (Bernoulli) as an *input* (their (R7) discussion p. around 468), and does not compute `ξ_4` on the square lattice. **NOT FOUND** as a theorem for square bond or square site.

### Connecting lemma (four-arm ⇒ `M' ≍ L^{3/4}` on a torus)

**NOT FOUND.** Russo’s formula says `d/dp P(crossing) = E[# pivotals]`. GPS convert four-arm to that expectation **in a planar quad** (Q12). Passing to `E[r_G]` on a *torus*, with wrapping pivotals and the matching term `pivotal_primal(p)+pivotal_matching(1-p)`, is not a named theorem. That is the hole between LSW and F1-as-the-ledger-writes-it.

**Rate probe may write.** **T:** “`α_4=−5/4` and `ν=4/3` are theorems on triangular site (LSW Thm 1.1, SW Thm 1). They are open on square site and square bond.” **C:** “GPS give four-arm ⇒ expected planar-quad pivotals; the torus/wrapping conversion is not literature.” **Blocked:** F1 on square site; F1 on the repository torus.

---

## Q6 — Finite-size *corrections* on lattices that *do* have conformal invariance

**Tag:** **FOUND** for the *limit* of crossing probabilities (Cardy / Smirnov). **NOT FOUND** as a theorem for `p_c(L)−p_c`, for crossing probability at `p_c` minus its limit, for any quantile of a homological event, or for any named `L^{-θ}`.

**Smirnov, arXiv:0909.4499v1, abstract:**

> We study scaling limits and conformal invariance of critical site percolation on triangular lattice. … As a particular case we obtain conformal invariance of the crossing probabilities and Cardy’s formula. Then we prove existence, uniqueness, and conformal invariance of the continuum scaling limit.

Cardy’s formula is a statement about `η → 0`, not about a finite-`L` expansion.

BCKS (published abstract only; PDF not retrieved) give a scaling *window* for `|C_max|` in a box for **bond** percolation, not a named `θ` for an inverse-CDF.

Simulations (Ziff, Newman–Ziff, Jacobsen, …) are **not** theorems. If they appear in #566 they remain empirical. Do not launder a fit into F1.

**Rate probe may write.** **T:** “Even on triangular site, there is no theorem expanding `Q_L(u)−p_c` in powers of `L`.” **Blocked:** any finite-`L` Cardy expansion; any Ziff exponent as T.

---

## Q7 — Grimmett–Li (1.1) vs (1.3), primary PDF

**Tag:** **FOUND**. Verification against arXiv:2205.02734v3 (RSA 65 (2024)). The #613 note does not over-claim.

**Grimmett–Li, arXiv:2205.02734v3, p. 2–3, verbatim:**

> Sykes and Essam presented motivation for the exact relationship
> (1.1) `p_c^{site}(G) + p_c^{site}(G*) = 1`,
> and this has been verified in a number of cases when `G` is amenable (see [6, 16]).
> Note that, since `G` is a subgraph of `G*`, it is trivial that
> (1.2) `p_c^{site}(G*) ≤ p_c^{site}(G)`.
>
> … it is proved, amongst other things, that
> (1.3) `p_u^{site}(G) + p_c^{site}(G*) = 1`,
> where `p_u^{site}` is the critical probability for the existence of a unique infinite open cluster.
> When `G` is amenable, we have `p_c^{site}(G) = p_u^{site}(G)`, in agreement with (1.1) (see [18, Chap. 7] …).

**Companion, Grimmett–Li, arXiv:2203.00981v3, Theorem 1.1(a) p. 3:**

> Let `(G1, G2)` be a matching pair derived from the mosaic `M ∈ Q`. We have that (1.1) `p_u^{site}(G1) + p_c^{site}(G2) = 1`.

**Remark 1.4 (Amenability), p. 4:**

> If `G ∈ Q` is one-ended and in addition amenable, by the uniqueness of the infinite cluster [1, 13], we have `p_c^{site}(G) = p_u^{site}(G)`; in this case, `p_c^{site}(G) ≥ 1/2` by (1.2). If `G` is transitive, we have `p_c^{site}(G) = 1/2` if and only if `G` is the usual amenable, triangular lattice.

The chain the #613 note uses is therefore exactly the published chain: (1.3) + amenability/uniqueness ⇒ (1.1). (1.1) is **not** a theorem of 2205.02734.

**van den Berg, JMP 22 (1981) 152–157, DOI 10.1063/1.524747.** **PDF not retrieved** (AIP 403). Published abstract, quoted from the AIP page:

> In 1964 Sykes and Essam obtained the relation `P_c^{(s)}(L) + P_c^{(s)}(L*) = 1`, where `L` and `L*` are a pair of matching lattices and `P_c^{(s)}` denotes the critical probability (site-case). The proof was not complete, but based on certain assumptions about the mean number of clusters. Though Sykes and Essam suggested that the above relation holds for all mosaics (i.e., multiply-connected planar graphs) and decorated mosaics, **we have constructed a counterexample**. Subsequently, **for a more restricted class of graphs**, an alternative derivation of the Sykes–Essam relation is given, this time based on the usual assumption that **below the critical probability the mean cluster size is finite**.

The class is “more restricted than all mosaics”; the extra assumption is finite mean cluster size below `p_c` (sharpness). That is the input Menshikov / AB87 / DCT later supply. **H3 has known counterexamples outside that class**; it holds here because `ℤ^2` is amenable and in the restricted class. Do not carry (H3) to an arbitrary planar matching pair.

**Rate probe may write.** **T:** “(H3) is (1.3) plus `p_c=p_u` on amenable `ℤ^2`, not (1.1) as a theorem of GL.” **Blocked:** citing (1.1) as GL’s theorem.

---

## Q8 — DCT site adaptation

**Tag:** **FOUND**. The #613 note does not over-claim. AB87 PDF not retrieved; the published abstract covers the note’s use.

**DCT arXiv:1502.03050v3, Theorem 1.1 item 3, p. 2:**

> If `(J_{x,y})_{x,y∈V}` is finite-range, then for any `β < β_c`, there exists `c = c(β) > 0` such that `P_β[0 ↔ Λ_n^c] ≤ e^{-cn}` for all `n ≥ 0`.

**§1.2, paragraph “Site percolation”, p. 3:**

> As in [AB87], the proof may be adapted to site percolation on transitive graphs.

Both `G` (NN on `ℤ^2`) and `Ĝ` (NN+NNN on `ℤ^2`) are locally finite, transitive, finite-range.

**Aizenman–Barsky, CMP 108 (1987) 489–526, DOI 10.1007/BF01212322.** **PDF not retrieved** (Springer/Euclid stubs). Published abstract:

> The equality of two critical points — the percolation threshold `p_H` and the point `p_T` where the cluster size distribution ceases to decay exponentially — is proven for **all translation invariant independent percolation models on homogeneous `d`-dimensional lattices** (`d ≥ 1`).

“All translation invariant independent percolation models on homogeneous lattices” includes finite-range **site** percolation on `ℤ^2` with NN or NN+NNN adjacency. The #613 note’s claim that AB87 covers the site adaptation is supported by the abstract; a theorem-body quote was not obtained.

Menshikov 1986 is the classical companion (Dokl. Akad. Nauk SSSR **288**); PDF not retrieved this pass; DCT’s bibliographical comment p. 3 names it jointly with AB87 for hypercubic Bernoulli percolation.

**Rate probe may write.** **T:** “(H2) is DCT Thm 1.1(3) plus the §1.2 site sentence, for both graphs.” **C:** “AB87 abstract covers independent percolation on homogeneous lattices; the journal PDF was not in hand.” **Blocked:** claiming a rate as `p ↑ p_c` (`c(β)→0` is not quantified).

---

## Q9 — Homological percolation after DKS

**Tag:** **FOUND** for location-type descendants. **NOT FOUND** for rates / window widths of power type, for square-site vs plaquette vs bond as a rate statement, for matching/dual *quantiles*, and for quantiles of a homological height / rank.

Forward search from arXiv:2011.11903:

- **Schweinhart–Shuman, arXiv:2601.00793v2.** Voronoi, `2i`-torus, threshold `1/2`. Location. No power window, no inverse-CDF.
- **Bobrowski–Skraba, arXiv:2005.14011** (and the earlier experimental 1803.06637 was a *wrong arXiv id* in one download pass — 1803.06637 is an unrelated PDE paper; the homological papers are 2005.14011 and the experimental companion cited by DKS as [BS20]). Thermodynamic-limit location of giant `k`-cycles; exponential tails outside that regime. Continuum, not square-site matching.
- DKS themselves already record Langlands–Pouliot–Saint-Aubin (1994), Pinson (1994), MDSA (2009) as discussion, not as rate theorems.

**NOT FOUND:** a paper that studies `Q_N(u)` for `u ∈ (0,1)` of `rank im H_1`, or a power `N^{-θ}` for `P(im H_i ≠ 0)` inside a shrinking window.

**Rate probe may write.** **T:** “Post-DKS homological percolation is still location-type.” **Blocked:** citing 2601.00793 or BS as a rate.

---

## Q10 — Does monotonicity + sharpness give a common `ω` for all interior `u`?

**Tag:** **NOT FOUND** as a named lemma in percolation or in probability textbooks under a special name (lead-lag of quantiles, Bahadur representation, monotone rearrangement were searched; none is this statement).

What *is* true, and is already the quantile step of Theorem L in the #613 note: if `F_N` is continuous and strictly increasing `[0,1]→[0,1]` and there is `w_N` with `F_N(p_c − w_N) → 0` and `F_N(p_c + w_N) → 1`, then for every compact `K ⊂ (0,1)` one has eventually `Q_N(K) ⊂ [p_c−w_N, p_c+w_N]`. That is the definition of the inverse of a monotone function, not a theorem that needs a name.

What is **not** automatic:

- A window proved only at a *single* `u` (DKS’s `λ□(N)` at `1/2`) does not, without a modulus of continuity of `F_N` or a uniform-in-`u` sharpness statement, control other quantiles at the *same* width. FK96 *does* give uniformity on `[ε,1−ε]` because of the `log(1/ε)` factor; that uniformity is still the FK window `1/log n`, not `L^{-1/ν}`.
- A `u_N → 0` (or `→1`) sequence is excluded by Theorem L’s compact-subset hypothesis; no literature lemma saves it.

The rate probe’s R5 (one exponent for the 9-vector is not implied by location) stands. Location + monotonicity give a *common qualitative limit*, not a common *exponent*.

**Rate probe may write.** **T:** “No named lemma upgrades single-`u` sharpness to a common percolation exponent for nine interior quantiles.” **C:** “If a window is already uniform on `[ε,1−ε]`, sharing it among inverse-CDFs is the definition of `F^{-1}`. That is the #613 quantile step, not a citation.” **Blocked:** “boolean sharpness ⇒ one `ν` for the 9-vector”.

---

## Q11 — Self-dual / self-matching finite-size laws

**Tag:** **FOUND** that Cardy is a *limit* at `p_c`. **NOT FOUND** for a finite-`L` expansion of `Q_L(u)−1/2` or of crossing-probability profiles at self-duality.

On triangular site and square bond, self-duality / self-matching locates the critical point at `1/2` (Kesten; Harris). The repository’s `p_L^H = 1/2` exactly is that fact plus finite-volume self-duality, not a shape theorem.

Smirnov 0909.4499 (Q6) is the scaling-limit Cardy formula. No PDF retrieved in this pass states a theorem of the form `P_{p_c}(crossing of L×L) = Cardy + C L^{-θ} + o(L^{-θ})` or `Q_L(u) − 1/2 = L^{-1/ν} F^{-1}(u) + …`.

**NOT FOUND** is useful: F2 ≡ 0 (vanishing of an operator-level correction) does not, by itself, give a rate for the shape of `Q_L`.

**Rate probe may write.** **T:** “Self-duality locates the median, not the profile.” **Blocked:** finite-`L` Cardy; F2⇒shape rate.

---

## Q12 — Connecting Russo pivotals to four-arm on a torus

**Tag:** **FOUND** on the triangular *plane*, for a bounded quad. **NOT FOUND** for torus-wrapping pivotals.

**Garban–Pete–Schramm, arXiv:1008.1378v5, p. 5 (PDF pp. 5–6):**

> a site `x` (i.e., a hexagon) is pivotal in a given configuration for left-right crossing if and only if there are open arms from it to `ab` and `cd` and closed arms to `bc` and `da`, i.e., if it satisfies the alternating four-arm event to `∂Q` … For most points `x ∈ Q` (in the “bulk”), the probability of this event is comparable … to `α_4^η(η, 1)`, the alternating four-arm probability from distance `η` to `1` on `ηT`. … it is not hard to show that the **expected number of pivotal points in `Q` is `≍_Q η^{-2} α_4^η(η, 1)`**, which is **known to be `η^{-3/4+o(1)}` on `ηT`**.

Conversion on `ℤ^2` (actually on `ηT`, triangular): `#` sites in a unit quad is `≍ η^{-2}`; each is four-arm with probability `α_4(η,1)`; product `η^{-2} · η^{5/4+o(1)} = η^{-3/4+o(1)}`. Russo then says the `p`-derivative of a crossing probability is this expectation.

**Extra to pass to a torus** (not in GPS, not found elsewhere as a theorem):

1. The event is wrapping / `r_G ≥ 1`, not a quad-crossing with four marked boundary arcs.
2. Periodic images: a site can be pivotal via an arm that wraps, so the relevant four-arm event is on the torus (or in the cover, to distance `∼ L`), not to the boundary of a simply-connected `Q`.
3. The repository identity `M'(p) = pivotal_primal(p) + pivotal_matching(1−p)` has no literature counterpart; GPS treat one colouring, not a matching pair.
4. “Not hard to show” in GPS is for the bulk of a planar quad with RSW/separation of interfaces (their Appendix), which is triangular-site technology.

**Rate probe may write.** **T:** “GPS convert four-arm → expected planar-quad pivotals on triangular site. The torus wrapping conversion is not a theorem.” **Blocked:** F1 on the repository torus; `M' ≍ L^{3/4}` as a citation.

---

## Q13 — What Theorem L cannot inherit from sharp-threshold boolean analysis

**Tag:** **NOT FOUND.**

No paper was found that claims: because `{r_G > 0}` (or `{ϕ_* ≠ 0}`) is a sharp monotone event, `Q_N(u) − p_c = O(N^{-α})` for an `α` taken from boolean influences, with `α` equal to the percolation `1/ν = 3/4`.

What boolean analysis *does* give, and DKS *use*, is a window of width `1/log N` (FK) or a fixed `ε` (DKS Thm 3). Theorem L already knows the qualitative limit; boolean analysis does not upgrade it to a percolation rate.

**Rate probe may write.** **T:** “Boolean sharpness is not a rate for `Q_N(u)−p_c`.” **Blocked:** any `O(N^{-α})` with `α` from influences.

---

## Q14 — X / Twitter

Unused. Q1–Q13 are not all NOT FOUND for every nearby object (windows exist; inverse-CDFs of `r_G` do not). No tweet is a quote.

---

## Square-site RSW: the #276 / #321 pointer (do not close)

Zeng, arXiv:1309.2273v1, Theorem 1.1, p. 1, unpublished:

> For any `λ > 0`, there exists `c = c(λ) > 0`, such that for all `n ≥ 1`
> `c ≤ P_{p_c}[ there exists a horizontal crossing of any n by λn box ] ≤ 1−c`,
> where `P_{p_c}` stands for critical site percolation on `ℤ^2`.

KST, arXiv:2011.04618, Comment 1: the published RSW homeomorphism for invariant positively associated *bond* percolation on `ℤ^2` “applies to … site percolation”.

This is RSW, not conformal invariance. It does not close #276 or #321. It does not prove F1.

---

## One-page sheet: what the rate probe may now write

### Allowed **T** sentences

1. Theorem L (`Q_N(u) → p_c` uniformly on compact subsets of `(0,1)`, no rate) is a scoped repository corollary of DCT + GL (1.3) + amenability + digital Alexander duality. It is not a new theorem and it is not a rate.
2. There is **no** named theorem that the inverse-CDF of `{r_G ≥ 1}` / `{r_G = 2}` / `E[r_G]/2` on a torus inherits a window of width `L^{-1/ν}` or `N^{-3/4}`.
3. DKS prove location-sharpness of giant cycles. Their “sharp” is `λ□(N) ± ε` (any fixed `ε`), plus exponential decay of `P(A□)` below threshold for `i=1` and `i=d−1`. `λ□(N)` is a single (`u=1/2`) quantile. `d=2,i=1` plaquette is square **bond**. The site model is permutohedral / triangular. Square-site + matching is not in print.
4. Friedgut–Kalai / Bourgain / Hatami give windows of size `p/log n` or `τ/p → 0`. They cannot produce `L^{-1/ν}` without RSW/arm inputs. There is **no** “sharp threshold ⇒ inverse-CDF window” lemma.
5. `α_4 = −5/4` and `ν = 4/3` are theorems on **triangular site** (LSW Thm 1.1, Smirnov–Werner Thm 1). They are **open** on square bond and square site.
6. The connecting lemma “four-arm ⇒ `M' ≍ L^{3/4}` on a torus” is **NOT FOUND**. GPS give the planar-quad conversion on triangular site.
7. Cardy’s formula is a **scaling limit**. There is **no** theorem giving a finite-`L` expansion of `Q_L(u)−p_c` or of `P_{p_c}(crossing) − Cardy`, on any lattice.
8. (H3) is GL (1.3) plus `p_c = p_u` on amenable graphs, not (1.1) as a theorem of the RSA paper. van den Berg has counterexamples to matching-sum-to-one outside a restricted class.
9. (H2) is DCT Thm 1.1(3) plus §1.2 site adaptation. AB87’s published abstract covers independent percolation on homogeneous lattices.
10. Post-DKS homological work (Bobrowski–Skraba, Schweinhart–Shuman) is location-type. **NOT FOUND** for quantiles of homological rank.
11. Monotonicity + single-`u` sharpness do **not** give a common percolation exponent for nine interior quantiles. Sharing a window that is already uniform on `[ε,1−ε]` is the definition of `F^{-1}` (already in #613).
12. Boolean sharpness does **not** imply `Q_N(u)−p_c = O(N^{-α})` with percolation `α`.
13. Square-site **RSW is FOUND** (Zeng Thm 1.1; KST Comment 1). Square-site **conformal invariance remains OPEN**. These are not in tension. Do not treat RSW as F1.

### Allowed **C** sentences (extra named)

1. **C, extra = RSW (Zeng / KST) + four-arm (open on this lattice):** crossing probabilities at `p_c` stay in `(c,1−c)` on square site; this still does not give a rate for `Q_N(u)`.
2. **C, extra = GPS conversion + torus wrapping (not a theorem):** *if* one assumes a torus analogue of “expected wrapping pivotals `≍ L^2 α_4(1,L)`” *and* `α_4 = −5/4`, then Russo would give `M' ≍ L^{3/4}`. The first ‘if’ is C, not T. On square site the second ‘if’ is also C.
3. **C, extra = Kesten scaling relations (Kesten 1987; PDF not in hand; quoted via SW):** on a lattice where polychromatic four-arm and one-arm exponents exist, they determine `ν`. Existence of those exponents is T only on triangular site.
4. **C, extra = BCKS finite-size scaling for `|C_max|` in a box (abstract only):** a scaling *window* for a different observable on bond percolation in a box. Not a citation for `Q_N(u)`.

### Still blocked

- F1 on square site, or F1 on the repository torus.
- Any `O(L^{-3/4})` / `O(N^{-3/4})` for `Q_N(u)−p_c` as a theorem.
- Fitting `ν` or `α_4` to #582 / N=725.
- Treating triangular-site theorems as square-site theorems.
- Treating Gate 1 “five widths” or #582 “five transitions” as literature.
- Citing (1.1) as Grimmett–Li’s theorem.
- Closing #276, #321, #613. Location is already a scoped corollary; quotation density does not upgrade it to a rate.
- Mixing this PR with #602 / #615 / #617 / #618 / #619.
- `docs/STATUS.md`.
- Affine charts, 1.55, the 4% residual, weighting (HIGH probe).
- Re-enumerating L=3,4.

---

## PDF-not-retrieved register

| item | what was obtained | what was not |
|---|---|---|
| van den Berg JMP 22 (1981) | published abstract (AIP), DOI 10.1063/1.524747 | journal PDF (HTTP 403) |
| Aizenman–Barsky CMP 108 (1987) | published abstract, DOI 10.1007/BF01212322 | journal PDF (Springer/Euclid stubs) |
| Burton–Keane CMP 121 (1989) | cited by GL Remark 1.4 and by Duminil-Copin 1712.04651; uniqueness under stationarity + finite energy | journal PDF |
| Borgs–Chayes–Kesten–Spencer CMP 224 (2001) | published abstract, DOI 10.1007/s002200100521 | full PDF (Microsoft TR 403; arXiv `math/9811085`, `math/9912001`, `math-ph/0103012` are different papers) |
| Kesten, *Scaling relations for 2D-percolation*, CMP 109 (1987) | used by Smirnov–Werner as [13]; scaling relations `θ,χ,ξ` from arm events | journal PDF |
| Russo 1978 / Seymour–Welsh 1978 | existence recorded by Duminil-Copin 1712.04651 and by KST 2011.04618 | journal PDFs |
| DHN, arXiv:0912.4253 | PDF downloaded; Type-3 fonts, body not quoteable | usable text extract |
| Menshikov 1986 | named by DCT §1.2 | Doklady PDF |

None of these holes changes a T sentence above. The AB87 abstract is enough to confirm the site-adaptation pointer; the vdB abstract is the counterexample sentence the brief asked for; BCKS is the wrong observable anyway.

