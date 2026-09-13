# Thermal `m λ^m` vs Jordan/LCFT diagnostics: literature retrieval

**Date:** 2026-09-12
**Ticket:** #714 (parent #650), retrieval-only
**Claim level:** literature status only; no repository claim is upgraded, downgraded or moved
**This note does not enter** `docs/STATUS.md`
**Scope guard:** this note does **not** adjudicate #275, does not power on any
machine, does not enumerate, and does not build a transfer matrix.

#710 records that `∂_p M_{4,m}(q4)` contains `m λ_*^m` produced by two simple
Perron branches crossing, and that the block is semisimple (toy:
`diag(λ+ε,λ-ε)` gives `d/dε[(λ+ε)^m-(λ-ε)^m]|_{ε=0} = 2mλ^{m-1}` with no Jordan
block anywhere). #275's thermal-jet line has treated `m λ^m` as possible
LCFT/Jordan evidence. This note records the **published status** of that
inference, from primary texts, with the actual diagnostics quoted.

## What was read

| # | Source | Where read | Tag |
|---|---|---|---|
| 1 | Vasseur, Jacobsen, Saleur, *Logarithmic observables in critical percolation*, [arXiv:1206.2312](https://arxiv.org/abs/1206.2312), J. Stat. Mech. **L07001** (2012) | full ar5iv HTML | PRIMARY_TEXT_READ |
| 2 | Cardy, *Logarithmic correlations in quenched random magnets and polymers*, [arXiv:cond-mat/9911024](https://arxiv.org/abs/cond-mat/9911024) | full ar5iv HTML | PRIMARY_TEXT_READ |
| 3 | Gurarie, *Logarithmic operators in conformal field theory*, [arXiv:hep-th/9303160](https://arxiv.org/abs/hep-th/9303160), Nucl. Phys. B **410**, 535 (1993) | full ar5iv HTML | PRIMARY_TEXT_READ |
| 4 | Creutzig, Ridout, *Logarithmic conformal field theory: beyond an introduction*, [arXiv:1303.0847](https://arxiv.org/abs/1303.0847), J. Phys. A **46**, 494006 (2013) | full ar5iv HTML, diagnostic sections | PRIMARY_TEXT_READ |
| 5 | Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley–Lieb algebras*, [arXiv:1507.03027](https://arxiv.org/abs/1507.03027), J. Phys. A **48**, 454003 (2015) | full ar5iv HTML | PRIMARY_TEXT_READ |
| 6 | Mertens, Ziff, *Percolation in finite matching lattices*, [arXiv:1603.07289](https://arxiv.org/abs/1603.07289), Phys. Rev. E **94**, 062152 (2016) | full arXiv HTML | PRIMARY_TEXT_READ |
| 7 | Bamieh, *A tutorial on matrix perturbation theory*, [arXiv:2002.05001](https://arxiv.org/abs/2002.05001) | full arXiv HTML | PRIMARY_TEXT_READ |
| 8 | Kato, *Perturbation Theory for Linear Operators*, Grundlehren **132**, Springer (reprint of the 1980 edition) | front matter + complete TOC from the Springer Classics scan; **body pages not accessed** | [LIT: primary text not verified] |
| 9 | J. O. Smith (Stanford CCRMA), *Eigenvalue sensitivity example* (Wilkinson-type), [ccrma.stanford.edu](https://ccrma.stanford.edu/~jos/matdoc/Eigenvalue_sensitivity_example.html) | full page | PRIMARY_TEXT_READ (secondary educational source) |
| 10 | Betcke, *Perturbation results for eigenvalue problems*, MATH0058 lecture notes (UCL), [tbetcke.github.io](https://tbetcke.github.io/math0058_lecture_notes/eigenvalues_perturbation_theory.html) | full page | PRIMARY_TEXT_READ (secondary educational source) |
| 11 | Pinson, *Critical percolation on the torus*, J. Stat. Phys. **75**, 1167 (1994), [doi:10.1007/BF02186762](https://link.springer.com/article/10.1007/BF02186762) | Springer, ADS and Scilit all bot-blocked; abstract not verified | [LIT: primary text not verified] |

## 1. The published LCFT/Jordan diagnostic (Q1)

What the primary texts actually require before calling something a Jordan
block / logarithmic operator:

**Gurarie 1993** — the defining object is a Jordan cell of `L0`, not a
factor in a derivative. Verbatim, his eqs. (13)/(19):

```text
L0|C,n>   = (h_C+n)|C,n>
L0|C_1,n> = |C,n> + (h_C+n)|C_1,n>
```

> "Ordinary primary operators are known to be the eigen vectors of the `L0`
> operators, and their eigen values are the dimensions of these operators.
> It will be shown that those "new" operators, which I will call
> pseudo-operators, are the basis of the Jordan cell for `L0`."

with the resulting two-point functions (his eq. (22))

```text
<C_1(z)C_1(w)> = -2/(z-w)^{2h_C} [log(z-w) + λ']
<C(z)C_1(w)>   =  1/(z-w)^{2h_C}
```

and the criterion

> "…we must include logarithmic operators in the theory if it possesses at
> least two operators the product of which when expanded according to the
> fusion rules … contains the contribution of at least two operators with
> the same dimension."

**Creutzig–Ridout 2013** — same statement in modern language:

> "… the corresponding field-theoretic models require, in addition, certain
> reducible, but indecomposable, representations. Such models have come to
> be known as logarithmic conformal field theories because the type of
> indecomposability required leads to logarithmic singularities in
> correlation functions."

> "… what happens if the primary field ϕ(z) corresponds to a state |ϕ⟩ which
> has a Jordan partner |Φ⟩ under the `L0`-action: `L0|Φ⟩ = h|Φ⟩ + |ϕ⟩`."

giving two-point functions of the form `⟨ΦΦ⟩ = (C − 2B log(z−w))/(z−w)^{2h}`
(their eq. (1.10)), with `B` physical and `C` basis-dependent.

**Cardy 1999** — the replica-limit mechanism, and what is *forbidden*:

> "… if two scaling dimensions `x_i` and `x_j` become degenerate in such a
> way that `A_ii ~ −A_jj → ∞` with `A_ii(x_i − x_j)` remaining finite, the
> leading terms will cancel leaving a logarithmic term proportional to
> `r^{−2x_i} ln r`."

> "Such operators should not occur in unitary conformal field theories, such
> as correspond to pure critical systems with positive Boltzmann weights…"

**Vasseur–Jacobsen–Saleur 2012** — the sharpest published diagnostic, at
`Q = 1`: two operators (energy `ε` and 4-leg/2-hull `ψ̂`) whose scaling
dimensions **degenerate** (`Δ_ε = Δ_ψ̂ = 5/4`); the log coefficient is fixed
by a **derivative of the dimension difference at the degeneracy**,

```text
lim_{Q→1} (Δ_ψ̂ − Δ_ε)/(Q−1) = √3/π        (their eq. (9))
⟨ψ̃(r)ψ̃(0)⟩ = 2A(1) r^{-5/2} [ … + (4√3/π) log r ]   (their eq. (8))
```

and the Jordan claim is then proven from the **scale-transformation mixing**

```text
ψ̃(Λr) = Λ^{-5/4} ( ψ̃(r) + (2√3/π) log Λ ε(r) )   (their eq. (11))
```

> "In other words, the scale transformation generator (or Hamiltonian) is
> non-diagonalizable, with a rank-2 Jordan cell mixing the two fields `ψ̃_ab`
> and `ε`."

Every published identification quoted above runs through **(i)** a degeneracy
of two scaling weights, **(ii)** a mixing of the two fields, and **(iii)** a
logarithm in a *correlation function* — never through a factor `m` in a
parameter derivative of a transfer-matrix trace.

## 2. Published warnings that crossings/derivatives are not Jordan evidence (Q2)

**The LCFT literature itself draws exactly the boundary #710 asks about.**
Vasseur–Jacobsen–Saleur 2012, verbatim:

> "In conclusion, it is important to stress that logarithmic terms such as
> those we have identified would not be present for generic `Q`, and occur
> solely because of the special degeneracies present at `Q=1`. This is of
> course quite different from logarithmic dependencies in other non-local
> quantities — see e.g. [27, 28] — **which are obtained as derivatives of
> correlation functions with respect to the Boltzmann weights (such as
> `Q`)**."

And Cardy 1999 says the same about the percolation connectivities, which are
precisely `q`-derivatives at `q = 1`:

> "… The connectivities of the percolation problem are given by the
> derivatives with respect to `q` at `q=1`, and are finite. But in this case
> **there are no logarithmic terms of the above form** [17]."

So the primary LCFT literature classifies parameter-derivative structures as
a **different** phenomenon from the Jordan-cell logarithms, not as evidence
for them.

**The algebra side.** Standard analytic perturbation theory treats a split
pair entirely inside the semisimple case. Bamieh 2020 (full text read),
verbatim hypothesis:

> "Throughout this note, we will assume the semi-simple case, i.e. that
> `A_ε` has a full set of eigenvectors (i.e. diagonalizable) for each `ε` in
> some neighborhood of zero."

Under this hypothesis a degenerate eigenvalue splits into **analytic
branches** `λ̄ + ε μ_j` (`μ_j` = eigenvalues of the perturbation restricted to
the eigenspace; his §3.2, Example 1: `A0 + εA1 = I + εM` gives
`λ_ε1 = 1 + εα, λ_ε2 = 1 + εβ`), and the first-order shift of each simple
branch is the familiar `λ_1i = w_0i* A_1 v_0i` (his eq. (25)). No Jordan
block enters; the word "Jordan" does not appear in the document. Applying
`x ↦ x^m` to analytic branches and differentiating is then ordinary calculus:

```text
d/dε [ (λ+ε)^m − (λ−ε)^m ]|_{ε=0} = 2mλ^{m-1}
```

which is exactly the #710 toy, with `tr(B^m)` eigenvalues `λ_i^m` of the
semisimple block. Nothing in this chain requires or produces a nontrivial
Jordan form.

**Kato** (textbook; [LIT] — body not verified, structure verified from the
book's own TOC in the Springer Classics scan): Chapter Two,
*"Perturbation theory in a finite-dimensional space"*, §1 *"Analytic
perturbation of eigenvalues"* with subsections *"Singularities of the
eigenvalues"* (§1.2) and *"Remarks and examples"* (§1.6) is where (i)
holomorphic dependence of simple eigenvalues and (ii) the singular
(fractional-power) behavior attached to genuinely defective eigenvalues are
treated. The operative contrast — **linear split ⇒ semisimple; square-root
split ⇒ genuine Jordan cell** — is the textbook content, but the exact
theorem text could not be fetched and is not quoted here.

**The numerical warning.** The Stanford (J. O. Smith, Wilkinson-type)
example read in full: along a generic perturbation direction a pair of
simple eigenvalues of the diagonalizable test matrix **coalesces and leaves
the real axis** (`eig(A + 0.5E) = 2.4067 ± 0.1753i`), with eigenvector
condition `cond(X) = 3.4×10^9` at the near-coalescence point. This is the
avoided-crossing/collision behavior of a *generic* (non-commuting) split —
the opposite limit from the commuting `diag(λ+ε,λ−ε)` toy, and in neither
direction does a Jordan block appear. Betcke's notes (read in full) make the
qualitative point directly:

> "If we have a polynomial with a multiple root then a small arbitrary
> perturbation in the coefficients will turn a multiple root into several
> simple roots."

No source read claims that a factor `m` from a crossing *is* Jordan evidence;
the published boundary runs the other way.

## 3. Do the finite-size matching/wrapping papers read `1/m` as Jordan? (Q3)

**No.** Both papers were read in full.

**Jacobsen 2015** ([arXiv:1507.03027]) is, ironically, built on an eigenvalue
crossing of exactly the #710 shape: the critical point is determined by
equating the largest **Perron–Frobenius** eigenvalues of two transfer-matrix
sectors,

> "`T_c(n)` is determined by equating the largest eigenvalues of two
> topologically distinct sectors of the transfer matrix."

i.e. `P_B(q,v) = 0 ⟺ Λ_open = Λ_closed` (his eq. (13)), unique and positive
by Perron–Frobenius away from criticality. The paper **never** uses the
words crossing, degeneracy, Jordan, LCFT or indecomposable for this; the only
multiplicity statement is combinatorial direct-sum multiplicity of identical
blocks. "Logarithmic corrections" appears once, for the `q = 4`
marginally-irrelevant operator — an unrelated sense of "log". Its Table 2,
row `n = 4`, reads

```text
0.5914171708531384817988341017359231779642
```

tabulated as the **square-lattice site percolation threshold estimate
`p_c(n)`** from an `n × ∞` basis — confirming #710's annotation that this
number is a finite-size estimate, "not a new `pc`".

**Mertens–Ziff 2016** ([arXiv:1603.07289]) relates average cluster numbers to
wrapping probabilities (`M_L(p) = R^x_L(p) − R̂^x_L(1−p)`, "This is the main
result of this paper."); the text contains **no** `ρ^m/m` sums, no
eigenvalue analysis at all, and no Jordan/log/indecomposable language
(logarithms appear only as log-log plot coordinates in figure captions).

Neither paper interprets any `1/m` displacement `ρ^m/m` as Jordan structure —
the question does not arise in their texts.

## 4. Kato: simple vs non-simple eigenvalues (Q4)

Summary of the published position (see §2 above for tags):

- **Simple eigenvalue** of a holomorphic family: analytic branch, first-order
  shift `λ_1 = w*A_1 v` (Bamieh eq. (25), read verbatim; Kato Ch. II §1,
  [LIT]).
- **Semisimple degenerate eigenvalue**: several analytic branches
  `λ̄ + ε μ_j`, obtained by diagonalizing the perturbation on the eigenspace
  (Bamieh §3.2/A.4, read verbatim). A split pair of simple branches crossing
  linearly is this case: the `m` in `2mλ^{m-1}` comes from the outer power,
  not from defectiveness.
- **Non-semisimple (defective) eigenvalue**: singular behavior; fractional
  powers of the perturbation (Kato Ch. II §1.2 *Singularities of the
  eigenvalues*, §1.6 *Remarks and examples*, [LIT]; not quoted because the
  body text was not accessible).

The diagnostic contrast is therefore published, but its `m λ^{m-1}` side is
the **semisimple** signature: linear split, analytic branches, no Jordan
cell. A nontrivial Jordan cell announces itself by *fractional* powers, not
by `m`.

## 5. Consequence for #714 / #275

| Proposition | Published status after this reading |
|---|---|
| Percolation `c=0` contains genuine LCFT Jordan cells | yes, established (Gurarie; Cardy; Creutzig–Ridout; Vasseur–Jacobsen–Saleur) — diagnosed by weight degeneracy + mixing + log in **correlations** |
| A log/factor produced by **differentiating with respect to a Boltzmann weight** is an LCFT/Jordan diagnostic | **no** — VJS: "quite different"; Cardy: percolation connectivities (`q`-derivatives) have "no logarithmic terms of the above form" |
| A split pair of simple eigenvalue branches gives `m λ^{m-1}` in the parameter derivative of the `m`-th power | yes — elementary algebra on analytic (semisimple) branches; Bamieh read in full; Kato Ch. II §1 [LIT] |
| Linear split ⇒ semisimple; fractional-power split ⇒ Jordan | published contrast (Kato Ch. II §1.2/§1.6 [LIT]; textbook-level) |
| Mertens–Ziff / Jacobsen read `1/m` in `ρ^m/m` as Jordan | **no such interpretation exists in either paper** (both read in full) |
| #275's `m λ^m` thermal jet is thereby resolved | **not claimed** |

## Not established

- that the published status settles the *empirical* question of what
  #275's thermal jet is; this note only records what the literature's
  diagnostics do and do not license;
- any statement about `∂_p M_{4,m}` beyond what #710 already derives
  internally (no transfer matrix was built here);
- Kato's theorem wording (body text not accessed; [LIT]);
- Pinson 1994 content beyond the citation (all full-text routes
  bot-blocked on 2026-09-12; [LIT]);
- anything about the repository's own thermal spin-4 LCFT candidate
  (`thermal-jordan-spin4-descendant.md`) — its Virasoro-side Jordan pair is a
  repository construction and is neither confirmed nor refuted by the sources
  above.

## Sources

1. arXiv:1206.2312 — Vasseur, Jacobsen, Saleur (2012), J. Stat. Mech. L07001.
2. arXiv:cond-mat/9911024 — Cardy (1999).
3. arXiv:hep-th/9303160 — Gurarie (1993), Nucl. Phys. B 410, 535.
4. arXiv:1303.0847 — Creutzig, Ridout (2013), J. Phys. A 46, 494006.
5. arXiv:1507.03027 — Jacobsen (2015), J. Phys. A 48, 454003.
6. arXiv:1603.07289 — Mertens, Ziff (2016), Phys. Rev. E 94, 062152.
7. arXiv:2002.05001 — Bamieh (2020/2022).
8. T. Kato, *Perturbation Theory for Linear Operators*, Grundlehren 132,
   Springer; Classics in Mathematics reprint of the 1980 edition (TOC
   verified from the Springer scan).
9. J. O. Smith, *Eigenvalue sensitivity example*, Stanford CCRMA online
   matdoc (accessed 2026-09-12).
10. T. Betcke, MATH0058 lecture notes, *Perturbation results for eigenvalue
    problems* (accessed 2026-09-12).
11. H. T. Pinson, J. Stat. Phys. 75, 1167 (1994), doi:10.1007/BF02186762
    ([LIT]).
