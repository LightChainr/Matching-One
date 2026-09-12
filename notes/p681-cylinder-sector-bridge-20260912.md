# #681 — restricted finite-torus to cylinder-sector map: dictionary, primary passages, one verdict

Date: 2026-09-12. Base read: `main@6edf775e` (includes #702–#704). Retrieval/proof
only. No transfer engine, no census, no production run. This note does **not**
enter `docs/STATUS.md` or `docs/ROADMAP.md`.

Scope reminder used throughout: `D = r_b − 1` and `M_{n,m} = P_2 − P_0` are already
theorems on honest tori (`notes/digital-alexander-duality-proof.md`,
`notes/exact-foundations-correction-20260912.md`). They are cited, not re-proved.
The unrestricted #675 operator no-go is withdrawn and is **not** cited as a theorem.

---

## 1. Typed dictionary

| object | definition used here | source |
|---|---|---|
| torus size | `(n,m)`: honest periodic square-cell torus, both directions identified; committed exact data at `n = m = L ∈ {3,4}`; the cylinder estimators use circumference `n` and `m → ∞` | repo; Jacobsen 2015 §3 |
| site weights | Bernoulli, product measure on `{0,1}^{nm}`: a configuration of size `k` has weight `p^k (1−p)^{nm−k}`; `v = p/(1−p)` (Potts/FK variable) | Mertens–Ziff (3); Jacobsen 2015 (3) |
| state space | all `2^{nm}` black/white colourings; `B` = black vertex set | repo enum |
| local black law | NN 4-connected induced graph `G_B`, ambient image `A ⊆ H_1(T^2;Q)` | digital-Alexander note |
| local white law | NN **and** NNN (8-connected) induced matching graph, complementary convention; ambient image `C = A^⊥` | digital-Alexander note |
| rank | `r_b = rank A ∈ {0,1,2}`, `r_w = rank C = 2 − r_b` (configuration-by-configuration theorem) | digital-Alexander note |
| rank-2 event | `A_2 = {r_b = 2}` (black spans both directions), mass `P_2` | repo |
| rank-0 event | `A_0 = {r_b = 0}` (black contractible, white spans both directions), mass `P_0` | repo |
| rank-1 event | `A_1 = {r_b = 1}` (one-direction wrapping), mass `P_1`; **does not appear in `M`** | repo |
| observable | `q = 1[r_b>0] − 1[r_w>0] = r_b − 1`; finite matching function `M_{n,m}(p) = E[q] = P_2 − P_0`; normalized `F = (1+M)/2` | repo; Mertens–Ziff |
| closure functional | Jacobsen's `n_wind = 0`: loops of non-trivial homotopy are removed with weight `n_wind`; this excludes `Z_1D` and forces `s = 0` states to be either **open** or **closed**, never mixed | Jacobsen 2015 §3, (9) |
| transfer blocks | `T̃ = ⊕_{k≥1} T^{(s=2k)} ⊕ T_open ⊕ T_closed`; ordering `Λ_open, Λ_closed > Λ^{(2)} > Λ^{(4)} > …` | Jacobsen 2015 (9), (10) |
| open sector | `Z_2D ~ (Λ_open)^m`: FK cluster spanning (both ways on the cylinder) — the black rank-2 event | Jacobsen 2015 (11) |
| closed sector | `Z_0D ~ (Λ_closed)^m`: dual FK cluster spanning — the white rank-2 / black rank-0 event | Jacobsen 2015 (12) |
| square-site R-matrix | `Ř_i = E_{i+2} E_i + v E_{i+1}`, with `q = 1`, `p = v/(1+v)` | Jacobsen 2015 (32) |
| spin sectors | `μ = 0` untwisted: dominant `Λ_0`, first subleading `Λ_0^{(1)}`; `μ = 1` twisted: `Λ_1`; ordered `Λ_0 > Λ_1 > …`, `Λ_μ = Λ_{q−μ}` | Jacobsen 2015 §9, (56) |
| abstract encoding (contrast) | configuration-indexed diagonal `W_p` with projectors `P_j` onto `r_b = j` gives `M = tr((P_2 − P_0) W_p)` — trivially possible, exponentially large, **not** a row transfer, **not** pTL, **not** a leading-eigenvalue identity | #702 note §2 |

**Terminology disambiguation.** The phrase "one-pTL sector criterion" is **not**
Jacobsen's and does **not** occur in arXiv:1507.03027 (checked against the full
HTML body). The criterion intended by the ticket is identified here as Jacobsen's
in-algebra statement `P_B = 0 ⇔ Λ_open = Λ_closed` (Eq. 13) inside the periodic
Temperley–Lieb (pTL) algebra, equivalently the spin-representation statement
`P_B = 0 ⇔ Λ_0^{(1)} = Λ_1` (Eq. 57). Both are stated for **finite `n`**, in the
**`m → ∞`** limit. The abstract states the same thing qualitatively: "`T_c(n)` is
determined by equating the largest eigenvalues of two topologically distinct
sectors of the transfer matrix."

---

## 2. Primary-source passages

### 2.1 Jacobsen 2015 — arXiv:1507.03027, J. Phys. A 48 (2015) 454003 — `PRIMARY_TEXT_READ`

Read from the arXiv HTML rendering of v1 (full body; 33 pp.). Verbatim:

- **(4)** `P_B(q,v) = Z_2D − q Z_0D.` Preceding definition: `Z_0D` sums edge
  subsets whose clusters all have trivial homotopy on toroidal boundary
  conditions; `Z_2D` "corresponds to clusters that wrap both periodic
  directions"; `Z_1D` ("wrapping one but not the other periodic direction") does
  not appear.
- **(9)** `T̃ = ⊕_{k=1}^{n} T^{(s=2k)} ⊕ T_open ⊕ T_closed.`
- **(10)** `Λ_open, Λ_closed > Λ^{(2)} > Λ^{(4)} > ⋯ > Λ^{(2n)}.`
- **(11)/(12)** `Z_2D ∼ (Λ_open)^m`, `Z_0D ∼ (Λ_closed)^m`.
- **(13)** `P_B(q,v) = 0 ⇔ Λ_open = Λ_closed`, "valid for a basis `B` of size
  `n × m`, with `n` finite and `m → ∞`." Preceding sentence: "If (4) is to have a
  (positive, unique) zero as a function of `v` … there must exist some
  `v_c(n) > 0` so that `Λ_open = Λ_closed`."
- **(22)–(23)** `n = 1` square case: lower block-triangular `T` with blocks
  `s = 2` and the two `s = 0` blocks, each 1-dimensional; eigenvalues
  `Λ^{(2)} = 4x^2`, `Λ_closed = (n_loop + 2x)^2`, `Λ_open = x^2(2 + n_loop x)^2`.
- **(24)** `Λ_open − Λ_closed = n_loop(x^2 − 1)(n_loop x^2 + 4x + n_loop)`, "is
  proportional to the graph polynomial `P_B(q,v) = (v^2 − q)(q + 4v + v^2)` for
  the `n × m = 1 × 1` basis. In this case, where all blocks are one-dimensional,
  it is obvious … that the relevant roots of `P_B(q,v)` are independent of `m`."
  **This is the only place in the source where proportionality is stated, and it
  is explicitly restricted to the 1-dimensional-block case.**
- **(32)** square-site R-matrix `Ř_i = E_{i+2} E_i + v E_{i+1}`.
- **(49)** `f_open(n) − f_closed(n) = o(n^{-2})`; **(50)**
  `p_c(n) − p_c = O(n^{-4})`, "and moreover the corrections appear to be
  `O(n^{-6})`, `O(n^{-8})`, and so on." NOTE: this is the **cylinder** estimator
  `p_c(n)` from the open/closed free-energy difference, not the torus matching
  root.
- **(51)–(55)** spin representation; **(56)**
  `Z_{μ,ν} = c_{μ,ν}(Λ_μ)^m + c_{μ,ν}^{(1)}(Λ_μ^{(1)})^m + …`, ordered
  `Λ_μ > Λ_μ^{(1)} > ⋯`, `Λ_0 > Λ_1 > …`.
- **(57)** `P_B(q,v) = 0 ⇔ Λ_0^{(1)} = Λ_1`, "valid for finite `n`, in the limit
  `m → ∞`."

### 2.2 Mertens–Ziff 2016 — arXiv:1603.07289, Phys. Rev. E 94, 062152 — `PRIMARY_TEXT_READ`

Read from the arXiv HTML rendering of v2 (full body, 8 pp.). Verbatim:

- **(2)** square matching polynomial `χ_□(p) = p − 2p^2 + p^4`.
- **(4)** `N_L(p) − N̂_L(1−p) − L^2 χ(p) = R_L^x(p) − R̂_L^x(1−p)`, where `x`
  runs over wrapping types. **(12)** the cross-wrapping form; **(17)**
  `M_L(p) = R^b_L(p) − R̂^b_L(1−p)`; **(20)** "`M_L(p) = R_L^x(p) − R̂_L^x(1−p)`,
  `x ∈ {c,b,e,h}` … This is the main result of this paper."
- **(11)** `N_L − N̂_L − (V − E + F_0) = +1 / −1 / 0` according as the black
  cross-wraps / the matching cross-wraps / neither.
- **(13)** `χ(p) = L^{-2}(⟨V⟩ − ⟨E⟩ + ⟨F_0⟩)` "is the matching polynomial of
  Sykes and Essam."
- **(15)** `M_L(p) = N_L(p) − N̂_L(1−p) − L^2 χ(p)` "the matching function."
- **(21)** Scullard–Jacobsen condition `R^c_L(p) − R^0_L(p) = 0`; the text proves
  it is "identical to the right-hand side of (12) being equal to 0", because
  `R^0_L(p) = R̂^c_L(1−p)` (no wrapping on the primary = cross-wrapping on the
  dual). Abstract: "The criterion that follows is related to the criterion
  Scullard and Jacobsen use to find precise approximate thresholds."
- **root.** "The matching function `M_L(p)` has a unique root `p*_L ∈ (0,1)`
  which converges to the critical density `p_c` as `L → ∞`. Empirically, the rate
  of convergence is `p*_L − p_c ∼ L^{−w}` with `w ≈ 4` [16,17]." They then give
  **(39)** `p*_L − p_c ∼ L^{2−x−1/ν}` and state the numerical value
  `w = 2 − x − 1/ν = −3.42 − 3/4 = −4.17`, "somewhat larger than the value 4
  suggested by Jacobsen [17]"; the measured slope is `−4.07`.
- **transfer/eigenvalue.** No "difference of leading eigenvalues" statement
  occurs in this paper; the only transfer-matrix mentions are methodological.

### 2.3 Jacobsen 2024 reply — J. Phys. A 57 258002, DOI 10.1088/1751-8121/ad4d33 — `ABSTRACT_ONLY`

Honest retrieval status:

- The reply is **not** on arXiv (checked the arXiv author listing and API; only
  v1 of 1507.03027 exists). The HAL record `hal-04980206` carries metadata and a
  one-sentence abstract only, with no deposited file: *"The authors replies to the
  comment made by Yang and Zhou (2024 J. Phys. A: Math. Theor.) on his 2015 paper
  entitled 'Critical points of Potts and O(N) models from eigenvalue identities in
  periodic Temperley–Lieb algebras' (Jacobsen 2015 J. Phys. A: Math. Theor. 48
  454003)."*
- The IOP full text is behind a Radware bot gate; the Wayback copy is rate-limited.
  Three short **body fragments** were recovered only from the publisher's own
  search index and are quoted here as fragments, not as primary body:
  "We illustrated the method in three different cases, providing, in particular,
  precise numerical values for these estimators `p_c(n)` for site …"; "In their
  comment, Yang and Zhou [3] have extended the series to `n = 24`. They also
  studied the same model with helical boundary …"; "They also studied the same
  model with helical boundary conditions, again giving estimators up to order
  `n = 24`. Based on these two …".

Consequence: the reply's **verdict cannot be inherited**. Nothing below depends
on it. Its fragments are consistent with Jacobsen treating `p_c(n)` as a
finite-size estimator, which is already established from §2.1.

Supporting (not one of the three named sources): Yang–Zhou 2024 Comment
J. Phys. A 57 258001, DOI 10.1088/1751-8121/ad4d2c — `ABSTRACT_ONLY`. Abstract:
they compute exact critical probabilities `h(n)` for an `n × ∞` helical square
lattice and extend the cylindrical `p_c(n)` to `n = 24`, concluding that
Jacobsen's `p_c = 0.59274605079210(2)` "is incorrect and the corrected value
should be `0.5927460507896(1)`." The body was not retrievable. This is a
numerical disagreement about the 13th decimal of an extrapolation; it does not by
itself bear on the operator bridge, and is recorded only as context.

---

## 3. Strongest equality actually justified at finite `n,m`

The exact finite statement available is Mertens–Ziff **(20)**:
`M_L(p) = R^x_L(p) − R̂^x_L(1−p)` at every `p` and every periodic size `L`, with
`x` any of cross/both/either/horizontal. Combined with repo theorem
`M_{n,m} = P_2 − P_0`, this identifies the finite torus matching function with a
**signed wrapping-amplitude difference**. It is a probability-root balance, not a
difference of leading eigenvalues.

Mertens–Ziff **(21)** + `R^0_L = R̂^c_L(1−p)` then identifies the same quantity
with the Scullard–Jacobsen critical polynomial, whose zero is Jacobsen's
`P_B = 0`. So at the level of the **polynomial and its zero**, the chain

```
root of M_{n,m} = 0   ⟺   R^c − R̂^c = 0   ⟺   Scullard–Jacobsen P_B = 0
```

is a theorem at finite size (Mertens–Ziff), and Jacobsen 2015 (13) identifies the
`m → ∞` limit of that same `P_B` with `Λ_open = Λ_closed`. What is **not** a
theorem at finite `m` is the replacement of `Z_2D − Z_0D` by `Λ_open − Λ_closed`;
Jacobsen's proportionality (24) is proved only for the `n = 1` case where every
block is 1-dimensional.

**Coefficient identity against already-committed data (a check, not a proof for
all `n,m`).** The committed exact honest-torus polynomials are
(`notes/probe-exact-controls-new-map-20260907.md`, reproduced by
`scripts/homological_balance/exact_torus_enum.py`):

```
L=3: M(p) = −1 + 6p^3 − 18p^7 + 18p^8 − 4p^9
L=4: M(p) = −1 + 8p^4 + 32p^6 − 64p^7 + 172p^8 − 704p^9 + 1104p^10
             − 608p^11 − 56p^12 + 128p^13 + 16p^14 − 32p^15 + 6p^16
```

with the committed exact values `M(0) = −1`, `M(1) = +1`, `M(1/2) = −21/64`
(L=3) and `−13757/32768` (L=4), `M'(1/2) = 225/64` (L=3), `4209/1024` (L=4), and
unique roots `p_L^H = 0.586511455…` (L=3), `0.590672112…` (L=4). Re-running the
committed enumerator reproduces `dual_fail = 0`, rank-pair counts
`(0,2)/(1,1)/(2,0) = 259/162/91` at L=3, and the same root. These are the
coefficient identities available to check a proposed formula; they do not
establish any statement for all `n,m`.

---

## 4. The `m → ∞` step

- **Which amplitudes overlap.** In Jacobsen's FK construction the two blocks that
  survive are `Z_2D ∼ c_open Λ_open^m` and `Z_0D ∼ c_closed Λ_closed^m` (11)–(12),
  and `Z_1D` is removed by `n_wind = 0`. The dominant `Λ_open^m` does **not**
  cancel between `Z_2D` and `Z_0D` (unlike the spin form (55), where `Λ_0^m`
  cancels and the survivors are `Λ_0^{(1)}` and `Λ_1`). The crossing compared in
  (13) is therefore between the two `s = 0` blocks directly.
- **Nonzero overlap.** The intermediate-value argument behind (13) requires the
  coefficients `c_open`, `c_closed` to be nonzero as functions of `v` near the
  crossing (nonzero top/bottom-slice overlap of the spanning and dual-spanning
  sectors). The sources assert this only implicitly; for the square-site
  `Ř_i` (32) it is not displayed.
- **Leading growth rates.** Eq. (10) orders `Λ_open, Λ_closed > Λ^{(2)} > …`;
  this strict dominance at `q = 1`, `v > 0` is what lets the two `s = 0`
  amplitudes govern the limit. Degeneracies (§9: "some of the inequalities might
  not be sharp when `q` takes particular values") are hypotheses to be stated,
  not conveniences.
- **Passing from probability-root balance to a crossing.** At finite `m`,
  `Z_2D − Z_0D = c_open Λ_open^m − c_closed Λ_closed^m + …` is **not**
  proportional to `Λ_open − Λ_closed` unless `c_open = c_closed` and the
  subleading `s ≥ 2` terms drop; only as `m → ∞` does the zero of the finite-`m`
  polynomial converge to the crossing (13). This is precisely the #702 caution,
  and it is why the finite-torus `M` must not be swapped for an eigenvalue
  difference.
- **Limit order.** Jacobsen is explicit: `n` finite, `m → ∞` first. The repo's
  `M_{n,m} = P_2 − P_0` lives at finite `m`. The two limits (`m → ∞` versus
  `n → ∞`) do not commute by assumption.
- A failure of one closure (e.g. wrong boundary word) is **not** a no-go for all
  operators: #702's abstract diagonal representation already exhibits *an*
  operator encoding `M`, and #675's blanket no-go is withdrawn.

---

## 5. Three claims kept separate

1. **`p_n^TL → p_c` (Jacobsen cylinder estimator).** Committed CSV
   `data/jacobsen_2015_square_site_cylinder.csv` (n = 1…21) is the eigenvalue
   identity sequence; n=1 is exactly `1/2` (consistent with (24)), n=21 is
   `0.592744551481…`, trending to the paper's `0.59274605079210(2)`. Jacobsen
   claims `p_c(n) − p_c = O(n^{-4})` from `f_open − f_closed = o(n^{-2})`
   ((49)–(50)); the committed first differences are consistent with a leading
   `~n^{-4}`.
2. **`p_L^H → p_c` (honest-torus matching root).** Committed exact L=3,4 roots
   `0.586511455`, `0.590672112`; the root of `M_L` is unique in `(0,1)` (Sturm).
   Mertens–Ziff report the rate only **empirically** (`w ≈ 4`), with their own
   derived `2 − x − 1/ν = −4.17` and measured `−4.07`, and note Jacobsen's `4`
   is a suggestion. **No matching-root exponent-4 theorem is found in any of the
   three named sources.** The rate is left open.
3. **Equality of the two finite estimators or of their leading correction laws.**
   They are **not equal at finite size**: at n=3 the cylinder value is
   `0.5888806999` and the torus matching root is `0.5865114551` (gap ≈ 2.37e-3);
   at n=4, `0.5914171709` vs `0.590672112…` (gap ≈ 7.45e-4). Mertens 2022's
   committed `p_med`/`p_cell` columns (`n^{-7/4}` sequences, `ν = 4/3`) are yet
   other finite estimators and also differ from both. The arithmetic
   `13/4 + 3/4 = 4` is conditional on **two separate** asymptotic inputs
   (`M_L(p_c) ∼ L^{-13/4}` and `M'_L(p_c) ∼ L^{3/4}`); it is not proof of either,
   and #613 supplies qualitative location only.

---

## 6. One scoped verdict: **conditional cylinder bridge**

There is a genuine, explicitly definable link, but it is conditional at the
operator level.

- **Proved (primary text).** Mertens–Ziff (20)/(21) identify the finite
  matching-function root with the Scullard–Jacobsen critical-polynomial zero;
  Jacobsen 2015 (13)/(57) identify that polynomial's zero — in the `m → ∞`,
  fixed-`n` limit — with the crossing of two topologically distinct `s = 0`
  sectors of the periodic Temperley–Lieb transfer matrix.
- **Conditional.** The step from the finite honest-torus object
  `M_{n,m} = P_2 − P_0` to the eigenvalue crossing requires, and is not supplied
  as a theorem by the sources: (i) `q = 1`; (ii) the `n_wind = 0` closure, which
  is exactly the projection onto the rank-2/rank-0 pair and the dropping of the
  `Z_1D` (rank-1) channel that `M = P_2 − P_0` already omits; (iii) the limit
  order `n` finite, `m → ∞`, with the strict dominance (10) and no `q = 1`
  degeneracy; (iv) coefficient nondegeneracy `c_open, c_closed ≠ 0` and, for any
  *finite-`m`* identification, `c_open = c_closed`, which is proved only for
  `n = 1` (24).

So: **a conditional cylinder bridge**. The polynomial-level chain is explicit;
the finite-operator identification is not, and the matching-root rate is not
established.

---

## 7. One smallest remaining obligation

**O1 (small exact algebra, no engine, no census).** At `n = 2` — the smallest
circumference whose `s = 0` open/closed blocks are **not** 1-dimensional —
compute exactly, from the square-site `Ř_i` of Jacobsen (32):

1. the two dominant `s = 0` eigenvalues `Λ_open(v)`, `Λ_closed(v)`;
2. the `m`-coefficients `c_open(v)`, `c_closed(v)` in `Z_2D`, `Z_0D` (56); and
3. test the finite-`m` identity at `q = 1`
   `Z_2D − Z_0D = c_2(v) (Λ_open^m − Λ_closed^m)` for all `m ≥ 1`,
   equivalently `c_open(v) = c_closed(v)` (or the weaker
   `Λ_open − Λ_closed ∝ P_B(q=1,v)` with a nonzero proportionality generalizing
   (24)).

Outcome split:

- **Holds** → the finite-torus `M_{n,m} = P_2 − P_0` connects to the sector
  crossing *explicitly at finite `n,m`*; the bridge upgrades from conditional to
  explicit.
- **Fails** (`c_open ≠ c_closed`) → the finite-torus root differs from the
  crossing by the coefficient ratio, the bridge is irreducibly a limit
  statement, and the limit order must be carried as a stated hypothesis.

This is the single step none of the three sources supplies, and it is decidable
by hand-sized exact computation. It does **not** decide the matching-root
exponent, which remains open independently.

---

## Provenance

| source | status | what was read |
|---|---|---|
| Jacobsen 2015, arXiv:1507.03027, J. Phys. A 48 454003 | `PRIMARY_TEXT_READ` | full HTML v1 body; Eqs. (3),(4),(9),(10),(11),(12),(13),(22),(23),(24),(32),(49),(50),(51)–(57) |
| Mertens–Ziff 2016, arXiv:1603.07289, Phys. Rev. E 94 062152 | `PRIMARY_TEXT_READ` | full HTML v2 body; Eqs. (2),(3),(4),(11)–(15),(17),(20),(21),(23),(31),(39)–(41) |
| Jacobsen 2024 Reply, J. Phys. A 57 258002, DOI 10.1088/1751-8121/ad4d33 | `ABSTRACT_ONLY` | HAL metadata + one-sentence abstract; three publisher-index body fragments; full text behind bot gate, not on arXiv |
| Yang–Zhou 2024 Comment, J. Phys. A 57 258001, DOI 10.1088/1751-8121/ad4d2c (supporting) | `ABSTRACT_ONLY` | abstract only |
| Mertens 2022, arXiv:2109.12102, J. Phys. A 55 334002 (committed data only) | not re-read as primary here | definitions of `p_med` (14a), `p_cell` (14b), `p_pol` (26); committed CSV columns |

No STATUS/ROADMAP edits. No merge. No production. The unique science lives in
this PR.
