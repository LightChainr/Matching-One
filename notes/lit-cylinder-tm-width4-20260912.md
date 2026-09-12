# #711 retrieval: cylinder/TM literature vs the width-4 all-p matching identity

2026-09-12. Retrieval only, independent of #708/#710. Parent #650. No census, no
transfer engine, no Huawei. This note is the deliverable
`notes/lit-cylinder-tm-width4-20260912.md` on a draft PR against `main`.

Subject of the comparison (from #711, not re-derived here): the draft all-p identity
for square-site matching on circumference 4,
`M_{4,m}(p) = (tr B_15^m - tr B_5^m)/(1+t)^{4m}`, generic scalar order 16, cylinder
root equal to Jacobsen 2015 Table 2 n=4, finite roots strictly below that root.

**Tripwire honoured.** `B_5`, `B_15`, `B_16` are *our* objects. They do not appear in
any source read here, and this note nowhere attributes them to a published paper. No
quotation below is reconstructed; each is copied from the fetched text or PDF.

## 0. Sources fetched, and their marking

| # | Source | Access route | Mark |
|---|---|---|---|
| S1 | Jacobsen 2015, *J. Phys. A* **48** 454003, arXiv:1507.03027v1 | arXiv abs + arXiv HTML v1 + arXiv PDF v1 | **PRIMARY_TEXT_READ** (§§2,4,6.1–6.2,7,8,9; Table 1; Table 2) |
| S2 | Scullard & Jacobsen 2012, *J. Phys. A* **45** 494004, arXiv:1209.1451v1 | arXiv abs + arXiv PDF v1 | **PRIMARY_TEXT_READ** (defs, Eq. (5); §3 transfer matrix) |
| S3 | Jacobsen & Scullard 2013, *J. Phys. A* **46** 075001 | IOP landing abstract | **ABSTRACT_ONLY** |
| S4 | Scullard & Jacobsen 2015, arXiv:1511.04374v1 (*Potts-model critical manifolds revisited*) | ar5iv HTML + arXiv PDF v1 | **PRIMARY_TEXT_READ** (§§2–3; Eq. (2)) |
| S5 | Mertens & Ziff 2016, *Phys. Rev. E* **94** 062152, arXiv:1603.07289v2 | ar5iv HTML + arXiv PDF v2 | **PRIMARY_TEXT_READ** (Eqs. (11)–(12),(20)–(24),(31),(39)) |
| S6 | Akhunzhanov, Eserkepov & Tarasevich 2022, *J. Phys. A* **55** 204004, arXiv:2204.01517v1 | arXiv abs + arXiv HTML v1 | **PRIMARY_TEXT_READ** (defs, Eqs. (4)–(9)) |
| S7 | Yang & Zhou 2024 Comment, *J. Phys. A* **57** 258001, DOI 10.1088/1751-8121/ad4d2c | Crossref/OpenAlex abstract; IOP landing | **ABSTRACT_ONLY** |
| S8 | Jacobsen 2024 Reply, *J. Phys. A* **57** 258002, DOI 10.1088/1751-8121/ad4d33 | Crossref/Semantic Scholar abstract; IOP landing | **ABSTRACT_ONLY** (body not reachable — see §5) |

Six primary texts were read; three of them (S1, S5, S6) were read from the full text,
not the abstract. S3, S7, S8 are marked honestly as abstract-only.

**Equation-numbering caveat for S1.** The arXiv HTML (LaTeXML) and the arXiv PDF v1 of
1507.03027 number the same displayed equations differently from roughly §5 onward. This
ticket (and `notes/p681-cylinder-sector-bridge-20260912.md`) uses the **HTML/LaTeXML**
numbers, so those are the ones printed here. For orientation, the PDF numbers of the
two equations the ticket names are: HTML (13) = PDF (13); HTML (50) = PDF (44). The
square-site R-matrix is HTML (32) = PDF (26); the free-energy definitions are HTML (48)
= PDF (42); the `o(n^-2)` statement is HTML (49) = PDF (43). Content, not number, is
what is quoted.

---

## 1. Jacobsen 2015 (S1): what is actually proved for the n=4 square cylinder

### 1.1 The object is a graph polynomial, defined as a *difference* at finite basis

Eq. (4) of S1 defines the critical polynomial as a difference of two partition-function
channels on a finite `n x m` basis:

> `P_B(q,v) = Z_2D - q Z_0D .`  (4)

So in S1 itself the "matching/graph polynomial is a difference" statement is **finite-m,
not a limit and not a root statement**. The three event weights `Z_2D`, `Z_1D`, `Z_0D`
are defined by the connectivity of a configuration when `B` is tiled into the infinite
lattice (S1, §2, Figure 1). The phrase "difference of two traces" is *not* used; the
difference is of two **partition functions / event weights**.

### 1.2 The eigenvalue identity (Eq. 13) and the block decomposition (Eq. 9)

S1 §4 takes `m -> infinity` first, so `B` is a semi-infinite cylinder of circumference
`n`, and replaces the transfer matrix by a direct sum indexed by the number of strings
`s`:

> `T~ = ⊕_{k=1}^{n} T^{(s=2k)} ⊕ T_open ⊕ T_closed .`  (9)

> "This is so precisely because contributions to `Z_1D` are excluded from (4), implying
> that loops winding around the cylinder carry the weight `n_wind = 0`."

Each of the `s = 0` sectors `T_open`, `T_closed` acts on reduced states whose count is
Eq. (8):

> `(1/2) C(2n,n) ~ 4^n .`  (8)

The main result is then an **eigenvalue identity**, proved by a monotonicity /
intermediate-value argument:

> `P_B(q,v) = 0  ⇔  Λ_open = Λ_closed ,`  (13)
>
> "valid for a basis B of size `n x m`, with `n` finite and `m -> infinity`."

The proof text immediately before (13):

> "For `v >> 1` the dominant contribution to (3) will be `A = E`, and hence
> `Λ_open > Λ_closed` by direct computation. Conversely, for `v << 1` the dominant
> contribution is `A = ∅`, whence `Λ_open < Λ_closed`. Since both terms in (4) grow
> exponentially in `m`, the factor of `q` is unimportant, and the intermediate value
> theorem implies our main result"

**What is proved vs assumed.** The *existence and uniqueness of the crossing*
`Λ_open = Λ_closed` is what (13) proves (for fixed `n`, `m -> ∞`); the *identification of
that crossing with the true `p_c`* is a separate conjecture about the method, as S2/S4
state explicitly (see §2.3). S1's own (10) is an ordering statement:

> `Λ_open , Λ_closed > Λ^{(2)} > Λ^{(4)} > ... > Λ^{(2n)} .`  (10)

Neither a 5+15+16 split nor a `tr B_15^m - tr B_5^m` identity appears anywhere in S1.
S1's decomposition is by **string number `s`**, with two extra `s = 0` sectors.

### 1.3 Convergence, and what is observed vs proved

S1 §4.3, on fixed-`n` convergence in `m`:

> "As expected, the results converge rapidly to the `m = infinity` limit, the rate of
> convergence being exponential in `m`."

No closed-form `c ρ^m/m` displacement is given for this regime.

On the `n`-dependence of the pseudo-critical point, S1 defines `f_open`, `f_closed`
(Eq. 48) and states:

> `f_open(n) - f_closed(n) = o(n^{-2})`  (49)
>
> "vanishes fast as `n -> infinity`, right at the critical point `p = p_c`."
>
> "This is a suggestive argument, but it does not quite explain the convergence
> properties of the eigenvalue method."

and then

> `p_c(n) - p_c = O(n^{-4}),`  (50)
>
> "and moreover the corrections appear to be `O(n^{-6})`, `O(n^{-8})`, and so on."
>
> "It is clear that more work would be required to establish whether (49) can be shown
> — obviously using more ingredients — to actually imply (50)."

So **(50) is an observed convergence law, explicitly not deduced from (49)**. The
positive-`n` counterpart is the finite-size-scaling ansatz Eq. (34)
`p_c(n) = p_c + Σ_k A_k / n^{Δ_k}` and its refined form Eq. (40)
`p_c(n) = p_c + Σ_k A_k / n^{2(k+1)}`, with `Δ_1 = 4.0001(2)`, `Δ_2 = 6.00(1)`.

### 1.4 The n=4 datum and the R-matrix

S1 Table 2 (site percolation, square lattice, `n x ∞` bases) contains, for `n = 4`:

> `4    0.5914171708531384817988341017359231779642`

The maximum size reached is `n_max = 21` ("Using the eigenvalue method we have obtained
the thresholds on `n x ∞` bases up to `n_max = 21`."). S1 §6.1 gives the square-site
R-matrix:

> `Ř_i = E_{i+2} E_i + v E_{i+1} ,`  (32)
>
> "This `Ř`-matrix contains only two out of fourteen possible terms ..."

where `E_i` are periodic-Temperley–Lieb generators.

### 1.5 Answer to Q1

* For site percolation on the square cylinder `n = 4`, S1 proves the **eigenvalue
  identity** (13) (`P_B = 0 ⇔ Λ_open = Λ_closed`, finite `n`, `m -> ∞`) and reports the
  numerical **Table 2** datum `0.5914171708531384817988341017359231779642`.
* It is **not** phrased as a wrapping-probability result, and it is **not** a closed
  finite-`m` trace identity. The graph-polynomial definition (4) is a finite-basis
  difference of `Z_2D` and `q Z_0D`; the eigenvalue criterion is its `m -> ∞` form.
* The `n`-convergence statement `O(n^{-4})` (50) is **observed**, and S1 says so.

---

## 2. Scullard–Jacobsen critical polynomials, 2012–2014 (Q2)

### 2.1 SJ 2012 (S2): the defining equality at finite basis

S2 defines the three event probabilities with normalization

> `P(0D;B) + P(1D;B) + P(2D;B) = 1`  (4)

and the criticality condition

> `P(2D;B) = P(0D;B) .`  (5)
>
> "Despite its apparent simplicity, eq. (5) is the main result of this paper."

with the interpretation

> "the unique root of `P_B(p)` in `[0,1]` either gives the exact percolation threshold
> for the lattice, or provides an approximation that becomes more accurate with
> appropriately increasing size of `B`."

The `1D` event is defined and used to classify configurations but does not enter (5).
S2 states the equivalence with contraction–deletion as an **open problem**, and the
exactness of the root as a property of the **infinite-`B` limit** (or exactly solvable
cases), not of finite `B`.

### 2.2 SJ 2015 (S4): the finite-`m` difference with the `q` coefficient

S4 §2 defines, for the `q`-state Potts model on a finite basis `B`,

> `P_B(q,v) = P_2D(q,v) - q P_0D(q,v) .`  (2)

S4 §3 is explicit that the transfer matrix computes the two weights separately and that
the root is the point where they balance:

> "Our earlier transfer matrix computations of critical polynomials would compute the
> weights of the 0D, 2D and 1D configurations and we could then set `P_2D = P_0D`.
> However, this is wasteful because `P_1D` is never used for anything."

### 2.3 Answer to Q2

* The difference `P_2D - q P_0D` is the **definition of `P_B` at finite basis**, in
  both S2 (as the equality `P(2D)=P(0D)`, no `q`) and S4/S1 (as `P_2D - q P_0D`,
  Eq. (2)). It is **not** introduced only at `m -> ∞` and **not** only at the root.
* The root is where the difference **vanishes** (`P_2D = q P_0D`; for percolation
  `q = 1`), and the interpretation of that root as `p_c` is the infinite-`B` limit (S2)
  or exact solvability.
* **It is not our object.** The two channels are `2D` vs `0D` global-connectivity event
  weights, not the digital-Alexander matching channels `P_2`,`P_0` of #705/#710, and
  not traces of two matrices. The structure "critical polynomial = difference of two
  channel weights, finite basis" is in print; the identification with
  `tr B_15^m - tr B_5^m` is not.

S3 (JS 2013) adds only the probabilistic-definition/transfer-matrix framing; its IOP
abstract reads, verbatim: "we give a probabilistic definition of `P_B(q,v)`, which
facilitates its computation, using the transfer matrix, on much larger `B` than was
previously possible."

---

## 3. Published finite-length displacement of a cylinder estimator (Q3)

The ticket asks who has `p_{n,m} - q_n ~ c ρ^m/m` for percolation wrapping/matching.
The retrieved components are:

**Fixed `n`, varying `m` (S1).** Convergence is "exponential in `m`" (§4.3, quoted in
§1.3). S1 gives no closed form and no `ρ^m/m` coefficient.

**Varying torus size `L` (S5, Mertens–Ziff).** The matching function is a finite-`L`
object:

> `M_L(p) = R^x_L(p) - R̂^x_L(1-p),  x ∈ {c,b,e,h}`  (20)
>
> "This is the main result of this paper."

and its root obeys

> "The matching function `M_L(p)` has a unique root `p*_L ∈ (0,1)` which converges to
> the critical density `p_c` as `L -> infinity`. Empirically, the rate of convergence is
> `p*_L - p_c ~ L^{-w}` with `w ≈ 4` [9,10]. This is significantly faster than the
> convergence of estimators derived from wrapping probabilities in the primary lattice
> alone, which converge like `p - p_c ~ L^{-2.75}` [6]."

S5's scaling analysis predicts (Eq. 39) `p*_L - p_c ~ L^{2-x-1/ν}`, numerically
`w = -4.17`, and Eq. (40) `p*_L - p_c ~ L^{2-y-3/ν} ≈ L^{-1.55}`; the integral estimator
(Eq. 41) converges as `L^{-1.65}`. **All of these are power laws in the linear size `L`,
not exponential in a length `m`.**

**Varying cylinder/torus size `L` (S6).** The exact cylinder spanning and torus wrapping
polynomials are computed up to `L = 16` (cylinder) and `L = 12` (torus); thresholds are
extracted by the universal-value, inflection, and crossing estimators, (4)–(6), and
extrapolated with the power-law series

> `p_c(L) = p_c(∞) + Σ_k A_k L^{-Δ_k}`  (9)

There is **no** exponential-in-`L` correction formula in S6.

**Finding for Q3.** No published formula of the form `p_{n,m} - q_n ~ c ρ^m/m` was found
for percolation wrapping/matching. What is in print is (i) exponential convergence in
the cylinder length `m` at fixed `n`, stated without a rate (S1 §4.3), and (ii) power-law
finite-size scaling in the linear size (S5, S6). The `ρ^m/m` form is the standard
root displacement of a difference of two exponential modes with equal leading
coefficients (as in the width-2 note's Eq. (3)); the searched literature applies the
exponential free-energy-correction machinery to the **Ising free energy**, which the
ticket explicitly excludes. This is a **negative retrieval result**, not a no-go
theorem: an absent formula is not evidence that no one has written it.

---

## 4. Transfer-matrix block decompositions and the "difference of two traces" (Q4)

**S1 (Jacobsen 2015) does decompose the transfer matrix — by string number.** Eq. (9)
(quoted in §1.2) gives `T~ = ⊕_{k=1}^{n} T^{(s=2k)} ⊕ T_open ⊕ T_closed`. For `n = 4`
this is **six** sectors: `T^{(2)}, T^{(4)}, T^{(6)}, T^{(8)}, T_open, T_closed`. Each of
`T_open`, `T_closed` acts on Eq. (8)'s `(1/2) C(2n,n)` reduced states, i.e. **35** at
`n = 4`; the full "complete" state count grows like `16^n` (S1 §3, citing Eq. (14) of
[13]). So a block decomposition of the cylinder transfer matrix at `n = 4` **is** in
print — but it is an `s`-indexed decomposition plus two `s = 0` sectors, **not a
5+15+16 split**, and the two distinguished sectors are pairings *within* `s = 0`.

**A matching function as a difference of two topological-sector quantities is in print**
(S5). `M_L(p)` is a difference of two wrapping probabilities (Eq. 20), and S5 identifies
it with the Scullard–Jacobsen criterion

> `R^c_L(p) - R^0_L(p) = 0`  (21)
>
> "This condition says that the probability of wrapping both ways is equal to the
> probability of wrapping neither way ..."
>
> "But the probability of no wrapping on the lattice is equal to the probability of
> cross-wrapping on the dual lattice `R^0_L(p) = R̂^c_L(1-p)`, and thus we see that (21)
> is identical to the right-hand side of (12) being equal to 0."

and S5 proves `M_L(p_c) = 0` exactly at finite `L` for self-matching lattices
(Eq. 22, "for all values of `L`") and for self-dual bond lattices (Eq. 24). This is the
closest published structure to "a matching function as a difference of two topological
sector quantities at finite size".

**Finding for Q4.** No 5+15+16 split was found; no identity of the form
`tr B_15^m - tr B_5^m` was found; no "difference of two traces" with those or comparable
trace dimensions was found. What is published is (a) S1's `s`-graded decomposition with
`T_open ⊕ T_closed`, and (b) S5's finite-`L` wrapping-probability difference. Whether
(a) and (b) can be composed into the width-4 trace difference of #710 is **not**
answered by any source read here.

---

## 5. Jacobsen 2024 Reply (Q5) — ABSTRACT_ONLY

**S8, Jacobsen 2024 Reply, *J. Phys. A* 57 258002, DOI 10.1088/1751-8121/ad4d33.**
The published abstract (Crossref, verbatim) is a single sentence:

> "The authors replies to the comment made by Yang and Zhou (2024 *J. Phys. A: Math.
> Theor.*) on his 2015 paper entitled 'Critical points of Potts and O(*N*) models from
> eigenvalue identities in periodic Temperley-Lieb algebras' (Jacobsen 2015 *J. Phys. A:
> Math. Theor.* **48** 454003)."

The **body is not reachable in this session**: the IOP PDF endpoint returns a Radware bot
page, the IOP landing page renders only page chrome, the paper is not indexed on arXiv
(an arXiv API title query returned 0 results), and the HAL record is behind an Anubis
challenge. **Marked ABSTRACT_ONLY.** A publisher search-index snippet (opening line only,
`[LIT]`, not primary-read) reads: "In their comment, Yang and Zhou [3] have extended the
series `p_c(n)` to `n = 24`. They also studied the same model with helical boundary ...".
This snippet is reported as a search-index artifact, not as a verified quotation, and no
argument in this note depends on it.

For context, the Reply answers the **Comment S7** (Yang & Zhou 2024, *J. Phys. A* 57
258001, DOI 10.1088/1751-8121/ad4d2c, ABSTRACT_ONLY; Crossref abstract read, verbatim):

> "We present an algorithm to compute the exact critical probability `h(n)` for an
> `n x ∞` helical square lattice with random and independent site occupancy. The
> algorithm has time complexity `O(n^2 c^n)` and space complexity `O(c^n)` with
> `c = 2.7459...` and allows us to compute `h(n)` up to `n = 24`. Since the extrapolation
> result of `h(n)` is inconsistent with the current best estimation of `p_c`, we also
> compute and extend the exact critical probability `p_c(n)` for an `n x ∞` cylindrical
> square lattice to `n = 24`. Our calculation shows that the current best result of
> `p_c = 0.592 746 050 792 10(2)` by Jacobsen (2015 *J. Phys. A: Math. Theor.* **48**
> 454003) is incorrect and the corrected value should be `0.592 746 050 7896(1)`."

Neither S7 nor S8 supplies a finite-length displacement formula of the Q3 type; their
subject is the disputed last digits of the `n -> ∞` extrapolation and the helical
boundary condition.

---

## 6. Compact answers

| Q | Answer | Sources |
|---|---|---|
| 1 | Square-cylinder site percolation `n=4`: proved **eigenvalue identity** (13) and a numerical Table 2 value; graph-polynomial defined as finite-basis difference (4); `O(n^{-4})` (50) **observed**, not deduced from (49). `B_5/B_15/B_16` do **not** appear. | S1 |
| 2 | `P_B = P_2D - q P_0D` is a **finite-basis definition** of the critical polynomial (S2 Eq. 5; S4/S1 Eq. 2), not only a limit or a root statement. Root is where it vanishes. It is *not* a difference of the digital-Alexander channels and *not* a difference of two traces. | S1, S2, S4 |
| 3 | No published `p_{n,m}-q_n ~ c ρ^m/m` found for percolation wrapping/matching. In print: exponential-in-`m` convergence at fixed `n` (no rate), and power-law-in-`L` scaling (`w≈4`, primary-lattice `L^{-2.75}`). Negative retrieval result. | S1, S5, S6 |
| 4 | A **string-graded** block decomposition with `T_open ⊕ T_closed` *is* in print; for `n=4` it is 6 sectors, not 5+15+16. A **matching function as a difference of two wrapping-probability sectors** is in print at finite `L`. No 5+15+16 split, no `tr B_15^m - tr B_5^m`. | S1, S5 |
| 5 | Reply body unreachable: **ABSTRACT_ONLY** (one-sentence published abstract). Comment also ABSTRACT_ONLY. | S7, S8 |

## 7. Explicitly not in print (checked, negative)

* No appearance of `B_5`, `B_15`, `B_16`, or `tr B_15^m - tr B_5^m`.
* No `5 + 15 + 16` sector split of a cylinder transfer matrix.
* No "matching polynomial = difference of two matrix traces" identity with those trace
  dimensions.
* No `c ρ^m/m` finite-length displacement for a percolation wrapping/matching cylinder
  estimator.
* No claim that the width-4 finite roots being strictly below the `n=4` root is novel —
  S1 itself exhibits `p_c(n) < p_c(n+1)` rising to the limit, and S5 exhibits a unique
  finite-`L` root converging to `p_c`.

## 8. Checks and caveats

* Every quotation above was copied from the fetched text (HTML) or from `pdftotext`
  output of the arXiv PDF; none was reconstructed from memory.
* S1's Table 2 `n=4` value `0.5914171708531384817988341017359231779642` was read from
  the arXiv PDF v1, line printed in §1.4, and agrees with the arXiv HTML rendering.
  This matches the `q4` quoted in #711 to the printed digits. That is a **citation of a
  published constant**, not a novelty claim.
* The equation-number discrepancy between the arXiv HTML and PDF renderings of S1 is
  recorded in §0 so that a later reader can reconcile this note with the PDF.
* "Matching" in S5 is the **matching lattice** (square lattice plus face diagonals); it
  is not the digital-Alexander matching observable of this repository. The two uses of
  the word are kept distinct throughout.
