# #613 — verifying #276's critical-point bridge

2026-09-07. Verdict: **the bridge holds as stated, under three imported
hypotheses, one of which is not automatic.** The proof supplied in
`docs/astra/ANSWER-610-20260907.md` Q4 is correct in structure; this note
supplies (i) an explicit conservative radius with the graph/Euclidean
conversion made precise, (ii) exact source links for the two external inputs,
(iii) the scope comparison that decides it is a **repository lemma, not a new
theorem**, and (iv) three failure modes that must stay outside the conclusion.

Verdict in one line: **not falsified, not novel — a scoped corollary.**

---

## 0. Statement to be proved

**Theorem L (repository lemma).** Let `T_N = R^2 / P_N Z^2` be honest periodic
square-cell tori with `N` vertices and shortest nonzero Euclidean period
`ell_N = min{|v| : v in P_N Z^2, v != 0}`, with

```text
(H0)  ell_N > sqrt(2)   and   ell_N / log N -> infinity.
```

Let `r_G` be the rank of `im[H_1(occupied NN graph) -> H_1(T_N; Q)]` and let
`r_hat` be the corresponding rank for the vacant set in the NN+NNN matching
graph. Assume

```text
(H1)  r_G + r_hat = 2   configuration-wise          (digital Alexander duality)
(H2)  for each of the two transitive finite-range graphs G (NN) and Ghat
      (NN+NNN), and every p below its own site critical point,
      P_p(0 <-> distance n) decays exponentially in n
(H3)  p_c^site(G) + p_c^site(Ghat) = 1
```

Put `F_N(p) = E_p[r_G]/2`. Then `F_N` is a polynomial, strictly increasing on
`(0,1)`, with `F_N(0)=0`, `F_N(1)=1`; and for every fixed `u in (0,1)`

```text
Q_N(u) := F_N^{-1}(u)  ->  p_c^site(G)   as N -> infinity,
```

uniformly on every compact subinterval of `(0,1)`. Any bounded family of
orientation weights `w_k` with `sum_k w_k = 1` gives
`sum_k w_k Q_N^{(k)}(u) -> p_c` as well.

For primitive Gaussian square tori `P_N = [[a,-b],[b,a]]`, `a^2+b^2 = N`,
`gcd(a,b)=1`, every nonzero lattice vector has squared length at least `N`, so
`ell_N = sqrt(N)` and (H0) holds.

---

## 1. The local lift and the union bound, made precise

### 1.1 What `r_G > 0` means

`r_G >= 1` iff some cycle in the occupied NN graph represents a nonzero class
in `H_1(T_N; Q)`. Because `pi_1(T^2)` is abelian, a loop on `T^2` is
null-homologous iff it is contractible, so a loop with nonzero class has a lift
to `R^2` running from `x` to `x + v` for some **nonzero** period `v in P_N Z^2`.
Conversely such a lifted path projects to a noncontractible loop. Hence

```text
{r_G > 0}  =  {exists a vertex x and a nonzero period v with x occupied-connected to x+v by an occupied NN path in the lift}
```

(the lift is along the covering `R^2 -> T_N`).

### 1.2 Graph distance versus Euclidean distance

This is the step the supplied proof left implicit. Both graphs embed in `R^2`
with explicit metric distortion:

| graph | edge set | `d(0,y)` bounds in terms of `|y|_2` |
|---|---|---|
| `G` (NN) | `(±1,0),(0,±1)` | `|y|_2 <= d_G(0,y) = |y|_1 <= sqrt(2) |y|_2` |
| `Ghat` (NN+NNN) | adds `(±1,±1)` | `|y|_2/sqrt(2) <= d_Ghat(0,y) <= sqrt(2)|y|_2` |

So a lifted NN path from `x` to `x+v` reaches `G`-graph distance at least
`|v|_1 >= |v|_2 >= ell_N`, and a lifted matching path reaches `Ghat`-graph
distance at least `|v|_2/sqrt(2) >= ell_N/sqrt(2)`.

### 1.3 The conservative radius

Take

```text
R_N = floor( ell_N / 4 ).
```

Two things must hold, and both do:

* **Embedded.** `B_G(x,R) ⊆ B_E(x,R)` because `|y|_1 <= R` implies `|y|_2 <= R`;
  so `B_G(x,R_N)` has Euclidean diameter `<= 2 R_N <= ell_N/2 < ell_N`, and the
  covering map is injective on it. For the matching graph
  `B_Ghat(x,R) ⊆ B_E(x, sqrt(2) R)`, diameter `<= 2 sqrt(2) R_N <= 0.708 ell_N
  < ell_N`, also injective. Injectivity is what licenses "the occupation law
  inside the ball is the infinite-graph product law".
* **Escaped.** By 1.2 the lifted path reaches graph distance `>= ell_N` (NN) or
  `>= ell_N/sqrt(2) = 0.707 ell_N` (matching), both `>= 2.8 R_N`, so it leaves
  the ball in either metric.

`R_N = ell_N/4` is therefore a single conservative radius for **both** graphs,
and any `R_N <= 0.35 ell_N` would do. The constant `1/4` is not optimal; it is
chosen so that the same number works for the NN and the NN+NNN adjacency.

### 1.4 The bound

Write `Λ_n^c` for the set of vertices at graph distance `> n`. By the union
bound over the `N` starting vertices, and using embeddedness to identify the
restricted law with the infinite-graph product law,

```text
P_p(r_G > 0)  <=  N * P_p^{Z^2}( 0 <-> Λ_{R_N}^c )  <=  N * A(p) * exp(-c(p) R_N)
               <=  N * A(p) * e^{c(p)} * exp(-c(p) ell_N / 4).
```

Setting `C(p) = A(p) e^{c(p)}` gives the form claimed in Q4,
`C(p) N exp[-c'(p) ell_N]` with `c' = c/4`; the factor `1/4` is absorbed into
the rate, which is why the two forms are the same statement. Finally

```text
N exp(-c'(p) ell_N) = exp( log N - c'(p) ell_N ) -> 0
    <=>   ell_N / log N -> infinity,
```

which is exactly (H0). Note the rate `c'(p)` depends on `p` and degenerates as
`p ↑ p_c`, so the conclusion is per-fixed-`p`, never uniform up to `p_c`.

### 1.5 The supercritical half

Under `P_p` the vacant set has the law of an independent configuration of
density `1-p`, so `P_p(r_hat > 0) = P_{1-p}^{Ghat}(r_hat > 0)`. Apply 1.4 to
`Ghat` at density `q = 1-p`: if `q < p_c(Ghat)` it tends to zero. By (H3),
`q < p_c(Ghat) <=> p > 1 - p_c(Ghat) = p_c(G)`. Since `r_G + r_hat = 2` by (H1),
`{r_G < 2} = {r_hat > 0}`, so

```text
P_p(r_G < 2) -> 0     for every fixed p > p_c(G).
```

Together with `F_N(p) = E[r_G]/2 = P(r_G=1)/... ` more directly:
`F_N(p) -> 0` for `p < p_c` and `F_N(p) -> 1` for `p > p_c`.

---

## 2. The quantile step

`r_G` is monotone under adding occupied sites (the occupied graph grows, so the
image of `H_1` grows), takes values in `{0,1,2}`, and satisfies `r(∅)=0`,
`r(V)=2`. Hence `F_N` is a polynomial in `p` with `F_N(0)=0`, `F_N(1)=1`.

**Strict increase.** Margulis–Russo: `d/dp E_p[f] = Σ_v E_p[δ_v f]` with
`δ_v f(ω) = f(ω ∪ {v}) - f(ω \ {v})`. Monotonicity gives `δ_v r >= 0`
everywhere. Take any maximal chain `∅ = ω_0 ⊂ … ⊂ ω_N = V`; since `r` goes from
`0` to `2`, at some step `r(ω_{k+1}) > r(ω_k)`, so with `η = ω_k`,
`v = v_{k+1}` we get `δ_v r(η) > 0`. That configuration has Bernoulli weight
`p^{|η|}(1-p)^{N-|η|} > 0` for every `0<p<1`, so `F_N'(p) > 0` on `(0,1)`.

Therefore `F_N` is a continuous strictly increasing bijection `[0,1] -> [0,1]`
and `Q_N = F_N^{-1}` is well defined and continuous on `(0,1)`.

**Convergence.** Fix `u in (0,1)` and `eps > 0` with `p_c ± eps in (0,1)`. By
§1, `F_N(p_c - eps) -> 0 < u` and `F_N(p_c + eps) -> 1 > u`, so for large `N`

```text
F_N(p_c - eps) < u < F_N(p_c + eps)  =>  p_c - eps < Q_N(u) < p_c + eps.
```

**Uniformity.** If `u` ranges over a compact `K ⊂ (0,1)`, the same two
inequalities hold simultaneously for all `u in K` once `F_N(p_c-eps) < min K`
and `F_N(p_c+eps) > max K`. The same argument gives `Q_N(u_N) -> p_c` for any
`u_N` that stays in a compact subset of `(0,1)`; in particular the homological
balance point `Q_N(1/2) -> p_c`, which is #276's `p_N`.

**Bounded orientation weights.** The repository combines **quantile functions**,
not distribution functions: `Q_N = Σ_k w_k Q_N^{(k)}` with `Σ_k w_k = 1`. Each
orientation's `Q_N^{(k)}` is covered by the theorem (its own `F_N^{(k)}` obeys
every step above), so `Q_N(u) -> p_c` provided `sup_N max_k |w_k| < ∞`. No
positivity of the weights is needed — which matters, because the `spin0`
combination has mixed signs. (Had one combined the *CDFs* instead, mixed-sign
weights would break monotonicity and the inverse would not exist; combining
quantile functions avoids that.)

---

## 3. The two external inputs, with sources

### (H2) Exponential decay, site percolation, both graphs

**Primary.** Duminil-Copin & Tassion, *A new proof of the sharpness of the
phase transition for Bernoulli percolation and the Ising model*,
Commun. Math. Phys. **343** (2016) 725–745; arXiv:1502.03050;
DOI 10.1007/s00220-015-2480-z.

* Abstract: "The proof applies to infinite range models on arbitrary locally
  finite **transitive** infinite graphs. … For **finite-range** models, we also
  prove that for any `β < β_c`, the probability of an open path from the origin
  to distance `n` decays exponentially fast in `n`."
* **Theorem 1.1, item 3** (verbatim): "If `(J_{x,y})_{x,y∈V}` is finite-range,
  then for any `β < β_c`, there exists `c = c(β) > 0` such that
  `P_β[0 <-> Λ_n^c] <= e^{-cn}` for all `n >= 0`."
* **§1.2, paragraph "Site percolation"** (verbatim): "As in [AB87], the proof
  may be adapted to site percolation on transitive graphs."
  This is exactly the site adaptation #613 asks to be discharged. Both `G`
  (NN on `Z^2`) and `Ghat` (NN+NNN on `Z^2`) are locally finite, transitive and
  finite-range, so (H2) holds for each.
* The reference `[AB87]` inside that sentence is M. Aizenman & D. J. Barsky,
  *Sharpness of the phase transition in percolation models*, Commun. Math.
  Phys. **108** (1987) 489–526, which treats general finite-range independent
  percolation models and therefore covers site percolation on both adjacencies
  directly. The classical companion is M. V. Menshikov, *Coincidence of
  critical points in percolation problems*, Dokl. Akad. Nauk SSSR **288**
  (1986) 1308–1311.

`P(0 <-> Λ_n^c)` in DCT is a **graph**-distance statement; §1.2 above converts
it to the Euclidean `ell_N` used in §1.4.

### (H3) The matching critical-point relation — and a warning

The relation (1.1) of

**G. R. Grimmett & Z. Li, *Percolation critical probabilities of matching
lattice-pairs*, Random Structures & Algorithms **65** (2024) 832–…;
DOI 10.1002/rsa.21226; arXiv:2205.02734** (verified: received 6 May 2022,
accepted 15 Feb 2024)

reads (verbatim from the published PDF):

```text
(1.1)   p_c^site(G) + p_c^site(G*) = 1
(1.2)   p_c^site(G*) <= p_c^site(G)                      (trivial, G ⊆ G*)
(1.3)   p_u^site(G) + p_c^site(G*) = 1
```

and the paper says of (1.1): "Sykes and Essam **presented motivation** for the
exact relationship (1.1), and this has been verified in a number of cases when
`G` is amenable (see [6, 14])." It then states: "When `G` is amenable, we have
`p_c^site(G) = p_u^site(G)`, in agreement with (1.1)."

So (1.1) is **not** a theorem of that paper; the theorem is (1.3), proved in the
companion **Grimmett & Li, *Hyperbolic site percolation*, arXiv:2203.00981**
(their ref. [12]), namely `p_u(G) + p_c(G*) = 1`. The chain we actually use is:

```text
(1.3)  p_u^site(G) + p_c^site(G*) = 1                 [Grimmett–Li, arXiv:2203.00981]
   +   p_c^site(G) = p_u^site(G)  for amenable G      [stated in the RSA paper; rests on
                                                       uniqueness of the infinite cluster,
                                                       Burton–Keane (1989), Comm. Math.
                                                       Phys. 121, 501–505]
   =>  (1.1)  p_c^site(G) + p_c^site(G*) = 1
```

with `G = Z^2` (amenable, transitive, one-ended, not a triangulation) and
`G* = G` plus all face diagonals `= Z^2` with NN+NNN `= Ghat`. This is exactly
the pair in the theorem, so (H3) is discharged.

Classical route, for readers who prefer it: Sykes & Essam (1964), J. Math.
Phys. **5**, 1117–1127, motivated (1.1); **J. van den Berg, *Percolation theory
on pairs of matching lattices*, J. Math. Phys. **22**(1) (1981) 152–157,
DOI 10.1063/1.524747** (= ref. [6] of Grimmett–Li) gave an alternative
derivation "for a more restricted class of graphs … based on the usual
assumption that below the critical probability the mean cluster size is
finite". That assumption is precisely sharpness, proved later by Menshikov
(1986) / Aizenman–Barsky (1987), so van den Berg's derivation plus sharpness
also yields (1.1) for `Z^2`.

**⚠ Warning that must stay in the record.** van den Berg's abstract states that
Sykes and Essam "suggested that the above relation holds for all mosaics …
**we have constructed a counterexample**", and that the derivation holds only
for a restricted class. So (H3) is a genuine hypothesis with **known
counterexamples in greater generality**; it is satisfied here specifically
because `Z^2` is amenable and lies in van den Berg's restricted class. Do not
carry (H3) to an arbitrary planar matching pair.

### (H1) is not imported

`r_G + r_hat = 2` is the repository's own digital Alexander duality
(`notes/digital-alexander-duality-proof.md`), proved there for honest periodic
square-cell tori. It is **assumed** here, not sourced from the literature. If
that note's cellwise argument is ever revised, Theorem L falls with it.

---

## 4. Scope: a repository lemma, not a new theorem

This is the item where the answer is "no novelty", and it should be said
plainly.

**Paul Duncan, Matthew Kahle & Benjamin Schweinhart, *Homological percolation
on a torus: plaquettes and permutohedra*, arXiv:2011.11903, published in
Ann. Inst. H. Poincaré Probab. Statist. (2025)** already prove that giant cycles
— their term for exactly the event `im[H_i(S) -> H_i(T^d)] != 0` that #276
calls `r_G > 0` — have a sharp threshold. Their abstract (v1 and v4) states:
"1-dimensional and `(d-1)`-dimensional plaquette percolation on the torus have
similar sharp thresholds at `p̂_c` and `1 - p̂_c` respectively, where `p̂_c` is
the critical threshold for bond percolation on `Z^d`", and "we prove that
`p_c = 1/2` in the case of middle dimension `i = d/2`". In `d = 2`, `i = 1`,
1-dimensional plaquette percolation **is** bond percolation on the square
torus, so

> the statement "the ambient-homology-rank event has a sharp threshold at
> `p_c`" is already in the literature for square-**bond** percolation on the
> torus, proved by essentially the union-bound argument reconstructed in §1.

Their site-percolation model is the **triangular** lattice, not the square
lattice, so the square-**site** case with the NN / matching pair used here is
not literally in DKS. It is obtained by the same argument plus (H3), which is
the one input DKS do not need. Related prior discussions they cite are
Langlands–Pouliot–Saint-Aubin (1994), Pinson (1994), MDSA (2009), and
Bobrowski–Skraba [BS20, BS22].

**Conclusion on novelty.** Theorem L is a **scoped repository lemma**: the
sharp-threshold mechanism is DKS's / standard finite-size percolation; the
square-site matching instance and the passage from a sharp threshold to
*convergence of every fixed quantile* are the repository's own additions. It
should be cited as "a corollary of Duminil-Copin–Tassion, Grimmett–Li, and the
digital Alexander duality", not submitted as a new theorem. Nothing here
supports a novelty claim.

---

## 5. What fails — three modes to keep out

1. **Thin tori / `ell_N = O(log N)`.** The bound `N e^{-c' ell_N}` no longer
   tends to zero, and the conclusion is genuinely false: on tori of bounded
   width the model is quasi-one-dimensional and has no phase transition at
   `p_c(Z^2)` at all. Also, for `ell_N <= sqrt(2)` the quotient is not honest:
   a unit cell can have two corners identified and an edge of the NN+NNN graph
   can close into a loop, so (H1) itself is not stated. **Keep `ell_N/log N ->
   infinity` and `ell_N > sqrt(2)` in the hypothesis.**
2. **Unbounded extrapolation weights.** With `sup_N max_k |w_k| = ∞` the
   weighted combination can amplify the `O(1)` finite-size error without bound
   and the limit is not controlled. Boundedness is load-bearing.
3. **A `u`-grid that approaches 0 or 1 with `N`.** The proof needs
   `F_N(p_c - eps) < u_N` and `F_N(p_c + eps) > u_N` for a *fixed* `eps`. For
   `u_N -> 0` (or `-> 1`) no such `eps` is available: `Q_N(u_N)` is only
   constrained to lie on the correct side of `p_c`, and can converge to
   anything in `[0, p_c]`. A quantified extension would require a lower bound
   on `u_N` against the (super-exponentially small) tail
   `F_N(p_c - eps) <= N C e^{-c' ell_N}`; that is a separate estimate and is
   **not** claimed here.

Two further non-claims, restated because they are easy to smuggle in:

* No **rate**. At `p = p_c` the argument says nothing; `F_N(p_c)` may converge
  to any value in `[0,1]`, and `c'(p) -> 0` as `p ↑ p_c`.
* No **value**. Theorem L locates the limit at `p_c^site(Z^2)`. It does not
  compute it, does not certify a numerical interval (that is #112's separate
  finite-event route), and says nothing about correction-to-scaling shape —
  which is what #610 is free to study now without pretending to locate `p_c`.

---

## 6. Verdict

**Not falsified.** No stated infinite-graph input fails for this model, and the
ambient-cycle event *is* bounded by embedded-ball exits, with the explicit
conservative radius `R_N = floor(ell_N/4)` and the graph/Euclidean conversion
of §1.2. `(H2)` is Duminil-Copin–Tassion Theorem 1.1(3) plus their §1.2 site
adaptation (with Aizenman–Barsky / Menshikov as the classical route). `(H3)` is
Grimmett–Li (1.3) plus `p_c = p_u` for amenable `G`, equivalently van den Berg
(1981) plus sharpness — with the recorded caveat that the relation has known
counterexamples outside the amenable/restricted setting. `(H1)` is the
repository's own digital Alexander duality and remains an assumption.

**Not novel.** The sharp threshold for giant cycles on a torus is
Duncan–Kahle–Schweinhart (arXiv:2011.11903). What is added here is the
square-site matching instance and the quantile formulation.

**What it buys #276.** Phase A/B of #276 now rests on a stated, sourced
qualitative theorem `Q_N(u) -> p_c` for the declared geometries, with no CFT
input and no fitted correction exponent. The remaining #276 content — the
*shape* of `Q_N(u) - p_c` — is untouched by this note and is exactly where
#610's finite-threshold-law work belongs.

## 7. Sources

| ref | item |
|---|---|
| [DCT16] | H. Duminil-Copin, V. Tassion, *A new proof of the sharpness of the phase transition for Bernoulli percolation and the Ising model*, Commun. Math. Phys. **343** (2016) 725–745. arXiv:1502.03050, DOI 10.1007/s00220-015-2480-z. Thm 1.1(3); §1.2 "Site percolation". |
| [AB87] | M. Aizenman, D. J. Barsky, *Sharpness of the phase transition in percolation models*, Commun. Math. Phys. **108** (1987) 489–526. |
| [Men86] | M. V. Menshikov, *Coincidence of critical points in percolation problems*, Dokl. Akad. Nauk SSSR **288** (1986) 1308–1311. |
| [GL24] | G. R. Grimmett, Z. Li, *Percolation critical probabilities of matching lattice-pairs*, Random Structures & Algorithms **65** (2024) 832–; DOI 10.1002/rsa.21226, arXiv:2205.02734. Eq. (1.1)–(1.3). |
| [GL22] | G. R. Grimmett, Z. Li, *Hyperbolic site percolation*, arXiv:2203.00981 (proves (1.3)). |
| [vdB81] | J. van den Berg, *Percolation theory on pairs of matching lattices*, J. Math. Phys. **22**(1) (1981) 152–157, DOI 10.1063/1.524747 (ref. [6] of [GL24]; contains the counterexample). |
| [SE64] | M. F. Sykes, J. W. Essam, *Exact critical percolation probabilities for site and bond problems*, J. Math. Phys. **5** (1964) 1117–1127 (ref. [16] of [GL24]). |
| [BK89] | R. M. Burton, M. Keane, *Density and uniqueness in percolation*, Commun. Math. Phys. **121** (1989) 501–505 (`p_c = p_u` on amenable transitive graphs). |
| [DKS25] | P. Duncan, M. Kahle, B. Schweinhart, *Homological percolation on a torus: plaquettes and permutohedra*, arXiv:2011.11903; Ann. Inst. H. Poincaré Probab. Statist. (2025). |
| repo | `notes/digital-alexander-duality-proof.md` — (H1). |
| repo | `docs/astra/ANSWER-610-20260907.md`, Q4 — the proof under review. |
