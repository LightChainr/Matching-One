# #613 — independent scope review of the #276 critical-point bridge

2026-09-08. WorkBuddy CLI (local machine). Bounded proof/literature task:
discharge the five verification items of #613 against the proof in
`docs/astra/ANSWER-610-20260907.md` (Q4, commit `65f5db1f`) and the
repository's digital Alexander duality
(`notes/digital-alexander-duality-proof.md`, workspace commit `8b5bc840`).
The note `notes/p613-quantile-convergence-20260907.md` (commit `f0981a98`,
merged via PR #614) already delivers Theorem L; this review re-derives and
checks its five verification items independently from the primary sources,
and records two scope corrections.

Verdict in one line: **the bridge survives — all five items discharge, with
two scope clarifications that must travel with the theorem.**

---

## Item 1 — the local lift/embedded-ball argument, with explicit radius

**Claim to check.** A nonzero ambient cycle (a loop representing a nonzero
class of `H_1(T_N; Q)`) has a lift exiting an embedded ball; inside such a
ball the occupation law coincides with the infinite-graph product law; a
union bound over `N` starting vertices gives
`P_p(r_G > 0) <= C(p) N exp[-c(p) ell_N] -> 0` for `p < p_c(G)`.

**Discharged.** Let `T_N = R^2/P_N Z^2` with shortest nonzero Euclidean
period `ell_N`. Because `pi_1(T^2)` is abelian, a loop is null-homologous
iff null-homotopic, so a nonzero class has a lift to the universal cover
with distinct endpoints `x, x+v`, `v in P_N Z^2 \ {0}`.

*Graph vs Euclidean distance (NN).* An NN path from `x` to `x+v` in the
cover has `l^1`-length at least `|v|_2 >= ell_N`, so it reaches
NN-graph distance `>= ell_N` from its start.

*Graph vs Euclidean distance (matching, NN+NNN).* A diagonal step covers
Euclidean distance `sqrt(2)`, so a matching path reaches Ghat-graph
distance `>= |v|_2/sqrt(2) >= ell_N/sqrt(2)`.

*Conservative embedded radius.* Take `R_N = floor(ell_N/4)`. An
NN-ball of graph radius `R` lies in the Euclidean ball of radius `R`
(since `|y|_1 <= R => |y|_2 <= R`), diameter `<= 2R_N <= ell_N/2 < ell_N`,
so the covering projection is injective on it. A Ghat-ball of graph
radius `R` lies in the Euclidean ball of radius `sqrt(2) R`, diameter
`<= 2 sqrt(2) R_N <= 0.708 ell_N < ell_N` — also injective. The escaped
distances (`ell_N` and `ell_N/sqrt(2)`) both exceed `2 R_N`, so any lifted
nonzero cycle exits the ball. One radius serves both graphs; any
`R_N <= 0.35 ell_N` would do.

*The bound.* On the embedded ball the restricted measure is the product
measure of the infinite graph (independence is local; embeddedness makes
the identification legitimate). For `p < p_c(G)`, exponential decay
(Item 2) gives `P^{Z^2}_p(0 <-> dist n) <= A(p) e^{-c(p) n}`; union bound
over the `N` starting vertices:

```text
P_p(r_G > 0) <= N A(p) e^{-c(p) R_N}
             <= C(p) N exp[-(c(p)/4) ell_N],
```

which tends to zero exactly when `ell_N/log N -> infinity`. The rate is
per-fixed-`p` and degenerates as `p ↑ p_c` — no uniformity up to `p_c` is
claimed, and none is needed for fixed quantiles.

**Scope correction 1 (sufficient, not necessary).** `ell_N/log N -> infinity`
is a *sufficient* condition for this fixed-off-critical union-bound
argument, not a proved *necessary* condition for quantile convergence. A
fixed-width torus genuinely fails (quasi-one-dimensional behaviour, no
transition at `p_c(Z^2)`), but a failure at one width `O(log N)` does not
establish failure for every width `O(log N)`; closing that would need a
separate argument and is out of scope. Keep the hypothesis as stated, but
label it sufficient.

*Honest quotients.* The cellwise digital Alexander proof requires each
periodic unit square to be an embedded cell with four distinct corners;
`ell_N > sqrt(2)` guarantees this. Tiny quotients with local
identifications are covered only by the separate finite oracle
(`scripts/digital_alexander_*_oracle.py`), not by the cellwise argument.
Short-period degeneracies are correctly outside the theorem.

## Item 2 — exponential decay for site percolation on both graphs

**Claim to check.** Subcritical site connection probabilities decay
exponentially on both the NN graph `G = Z^2` and the NN+NNN matching graph
`Ghat`, with no unqualified bond theorem substituted.

**Discharged.** Primary source:

> H. Duminil-Copin, V. Tassion, *A new proof of the sharpness of the phase
> transition for Bernoulli percolation and the Ising model*,
> Commun. Math. Phys. **343** (2016) 725–745; arXiv:1502.03050;
> DOI 10.1007/s00220-015-2480-z.

* **Theorem 1.1, item 3** (finite-range case, verbatim): "If `(J_{x,y})` is
  finite-range, then for any `β < β_c`, there exists `c = c(β) > 0` such
  that `P_β[0 <-> Λ_n^c] <= e^{-cn}` for all `n >= 0`." The main theorem is
  stated for **bond** percolation with general transitive couplings.
* **§1.2, remark "Site percolation"** (verbatim): "As in [AB87], the proof
  may be adapted to site percolation on transitive graphs", with the
  modified crossing functional `φ_p(S) = Σ_{x∈S} Σ_{y∉S, {x,y}∈E}
  P_p[0 ↔ x]`. This is the site adaptation #613 asked to be discharged,
  and it is a remark, not a separately numbered theorem — the citation
  should say so.
* The remark's `[AB87]` is M. Aizenman & D. J. Barsky, *Sharpness of the
  phase transition in percolation models*, Commun. Math. Phys. **108**
  (1987) 489–526, which treats general finite-range independent percolation
  (site included) directly; M. V. Menshikov, Dokl. Akad. Nauk SSSR **288**
  (1986) 1308–1311 is the classical companion. Either route covers both
  `G` and `Ghat`: both are locally finite, transitive, finite-range, so
  the site adaptation applies to each at every `p` below its own site
  critical point.

The DCT statement is a **graph**-distance bound; Item 1's embedded-ball
conversion turns it into the Euclidean `ell_N` form. No unqualified bond
theorem is used anywhere.

## Item 3 — the amenable matching critical-point relation

**Claim to check.** `p_c^site(G) + p_c^site(Ghat) = 1` for `G = Z^2`,
sourced from a probability theorem, not from Mertens–Ziff's analyticity
heuristic.

**Discharged, with a citation precision the issue's pointer does not
capture.** The issue points at Grimmett–Li, DOI 10.1002/rsa.21226
(G. R. Grimmett, Z. Li, *Percolation critical probabilities of matching
lattice-pairs*, Random Structures & Algorithms **65** (2024), arXiv:2205.02734;
Crossref-verified: received 6 May 2022, accepted 15 Feb 2024). That paper's
Theorem 3.1 proves the **strict inequality** criterion
(`p_c(G*) < p_c(G)` iff `G` is not a triangulation, for transitive `G`) and
states (1.1)–(1.3) in its introduction, where (1.1) is
`p_c^site(G) + p_c^site(G*) = 1` — but (1.1) is presented there as the
Sykes–Essam relation "verified in a number of cases when `G` is amenable",
not as that paper's theorem.

The theorem actually proving the equality is the companion:

> G. R. Grimmett, Z. Li, *Hyperbolic site percolation*,
> Random Structures & Algorithms (2024, accepted v3); arXiv:2203.00981.

Its abstract: "If `(G1,G2)` is a matching pair derived from some
quasi-transitive mosaic `M`, then `p_u(G1) + p_c(G2) = 1`", and "when `G`
is amenable we have `p_c(G) = p_u(G)`". The discharge chain is:

```text
p_u^site(G) + p_c^site(G*) = 1        [Grimmett–Li, arXiv:2203.00981]
p_c^site(G) = p_u^site(G)  for amenable G
                              [uniqueness of the infinite cluster;
                               Burton–Keane, Comm. Math. Phys. 121
                               (1989) 501–505]
=> p_c^site(Z^2) + p_c^site(Ghat) = 1,
   with G* = Z^2 + face diagonals = Ghat.
```

Classical alternative: J. van den Berg, *Percolation theory on pairs of
matching lattices*, J. Math. Phys. **22**(1) (1981) 152–157,
DOI 10.1063/1.524747, derives the relation for a restricted class under the
(finite mean cluster size) assumption later proved by Menshikov/Aizenman–
Barsky — so van den Berg plus sharpness also yields the equality for `Z^2`.

**Warning that must stay in the record.** van den Berg's abstract records
that Sykes–Essam's suggested general relation has a constructed
counterexample outside the restricted (amenable) class. So (H3) is a
genuine hypothesis with known failures in greater generality; it holds here
specifically because `Z^2` is amenable, transitive and one-ended. Do not
carry (H3) to an arbitrary planar matching pair. Mertens–Ziff
(arXiv:1603.07289) supply the *finite* matching identity used in the
repository's duality note; their introduction's analyticity motivation is,
as the issue says, not the probability theorem and is not used as one.

## Item 4 — scope against the existing literature; lemma, not new theorem

The matching-root / homological-percolation comparison:

> P. Duncan, M. Kahle, B. Schweinhart, *Homological percolation on a torus:
> plaquettes and permutohedra*, arXiv:2011.11903 (v4 2023); published,
> Ann. Inst. H. Poincaré Probab. Statist. (2025).

Their giant cycles are exactly the event `im[H_i(S) -> H_i(T^d)] ≠ 0` that
the repository calls `r_G > 0`. They prove a sharp transition from
nonexistence of giant cycles to giant cycles spanning the ambient homology,
and — directly relevant — "we also prove convergence of the threshold
function to a constant in certain cases", with `p_c = 1/2` at middle
dimension. In `d = 2, i = 1`, their plaquette model **is** bond percolation
on the square torus, so the union-bound mechanism reconstructed in Item 1
is already in the literature for square-bond tori. Their site model is the
**triangular** lattice, not the square; the square-site NN/matching-pair
instance is not literally in DKS.

**Classification.** Theorem L is a **scoped repository lemma**: the
sharp-threshold mechanism is Duncan–Kahle–Schweinhart / standard finite-size
percolation theory (Margulis–Russo for strict monotonicity; the quantile
squeezing is elementary); the genuinely repository-specific additions are
(i) the square-site matching instance, which additionally consumes (H3),
and (ii) the digital Alexander duality (H1) that converts the vacant-set
matching rank into the complementary rank — the 16-pattern regular-
neighbourhood certificate in `notes/digital-alexander-duality-proof.md` is
the piece with no literature counterpart. No novelty claim is supportable,
and none should be made. Cite as: a corollary of Duminil-Copin–Tassion
(site sharpness/exponential decay), Grimmett–Li (matching critical-point
equality), Duncan–Kahle–Schweinhart (giant-cycle threshold mechanism), and
the repository's digital Alexander duality.

## Item 5 — what fails, and stays outside

1. **Thin tori (`ell_N = O(log N)` or bounded width).** The union bound no
   longer tends to zero; at bounded width the model is quasi-one-dimensional
   and has no transition at `p_c(Z^2)` at all. Also `ell_N <= sqrt(2)`
   breaks the honest-cell hypothesis, so (H1) itself is not stated there.
   Both the growth hypothesis and `ell_N > sqrt(2)` are load-bearing.
   (Cf. Scope correction 1: sufficiency, not necessity.)
2. **Unbounded extrapolation weights.** A weighted quantile combination
   `Σ_k w_k Q_N^{(k)}` with `Σ_k w_k = 1` converges iff
   `sup_N max_k |w_k| < ∞` (each `Q_N^{(k)}` is covered individually;
   mixed signs are fine because quantile functions, not CDFs, are
   combined). With unbounded weights the `O(1)` finite-size error can be
   amplified without bound — boundedness is not cosmetic.
3. **A u-grid approaching 0 or 1 with `N`.** The squeeze needs a *fixed*
   `eps` with `F_N(p_c - eps) < u_N < F_N(p_c + eps)`. For `u_N -> 0` no
   such `eps` exists; `Q_N(u_N)` is then only constrained to the correct
   side of `p_c` and can converge to anything in `[0, p_c]`. A quantified
   extension would need a lower bound on `u_N` against the
   super-exponentially small tail `N C e^{-c' ell_N}` — a separate estimate,
   not claimed.

Non-claims restated: no rate at or near criticality (`c'(p) -> 0` as
`p ↑ p_c`; `F_N(p_c)` unconstrained in `[0,1]`); no numerical value or
certified interval for `p_c^site(Z^2)` (that is #112's separate finite-event
route); no correction-to-scaling shape, CFT operator content, or the value
3/4 for square-site thermal scaling.

---

## Verdict

All five items discharge:

| item | status | carrier |
|---|---|---|
| 1. lift/embedded-ball union bound | **proved**, explicit `R_N = floor(ell_N/4)` | this note §Item 1; consistent with Theorem L §1 |
| 2. site exponential decay, both graphs | **cited**: DCT Thm 1.1(3) + §1.2 site remark (Aizenman–Barsky/Menshikov route) | this note §Item 2 |
| 3. amenable matching relation | **cited**: Grimmett–Li *Hyperbolic site percolation* (equality) + Burton–Keane amenability; RSA 2024 paper states but does not prove (1.1) | this note §Item 3 |
| 4. novelty scope | **repository lemma, not a new theorem**; DKS holds the giant-cycle threshold (bond case) | this note §Item 4 |
| 5. failure modes | thin tori, unbounded weights, `u_N -> {0,1}` all correctly excluded | this note §Item 5 |

**Expected result confirmed, not falsified.** The #276 Phase A/B bridge now
rests on stated, sourced inputs: qualitative `Q_N(u) -> p_c^site(Z^2)` for
honest square-cell tori with `ell_N/log N -> infinity`, no CFT input, no
fitted exponent. Two scope corrections travel with it: (a)
`ell_N/log N -> infinity` is sufficient, not proved necessary, for the
quantile conclusion; (b) the RSA-paper pointer for item 3 should be cited
through the *Hyperbolic site percolation* theorem, with the van den Berg
counterexample caveat attached to any generalization beyond amenable pairs.

No Monte Carlo was run; no simulation artifact is touched. STATUS.md
untouched; nothing merged; no ticket closed.

## Sources

| ref | item |
|---|---|
| [DCT16] | Duminil-Copin & Tassion, CMP **343** (2016) 725–745; arXiv:1502.03050; DOI 10.1007/s00220-015-2480-z. Thm 1.1(3); §1.2 "Site percolation". |
| [AB87] | Aizenman & Barsky, CMP **108** (1987) 489–526. |
| [Men86] | Menshikov, Dokl. Akad. Nauk SSSR **288** (1986) 1308–1311. |
| [GL24a] | Grimmett & Li, *Matching lattice-pairs*, Random Struct. Alg. **65** (2024); DOI 10.1002/rsa.21226; arXiv:2205.02734. States (1.1)–(1.3); proves strict-inequality criterion (Thm 3.1). |
| [GL24b] | Grimmett & Li, *Hyperbolic site percolation*, Random Struct. Alg. (2024); arXiv:2203.00981. Proves `p_u(G1)+p_c(G2)=1` and `p_c=p_u` for amenable `G`. |
| [vdB81] | van den Berg, J. Math. Phys. **22**(1) (1981) 152–157; DOI 10.1063/1.524747. Restricted-class derivation; records the Sykes–Essam counterexample. |
| [SE64] | Sykes & Essam, J. Math. Phys. **5** (1964) 1117–1127. |
| [BK89] | Burton & Keane, CMP **121** (1989) 501–505. |
| [DKS25] | Duncan, Kahle & Schweinhart, arXiv:2011.11903; Ann. Inst. H. Poincaré Probab. Statist. (2025). |
| [MZ17] | Mertens & Ziff, arXiv:1603.07289 (finite matching identity). |
| repo | `notes/digital-alexander-duality-proof.md` (H1, 16-pattern certificate). |
| repo | `notes/p613-quantile-convergence-20260907.md` (Theorem L, commit `f0981a98`). |
| repo | `docs/astra/ANSWER-610-20260907.md` Q4 (the proof under review, commit `65f5db1f`). |
