# Literature officer — #601: five things #598 may be rediscovering

**Date:** 2026-09-06
**Role:** theory input. Same epistemic level as an Astra answer and as [#572](https://github.com/LightChainr/Matching-One/pull/572). **Does not enter** `docs/STATUS.md`. **Does not score** a frozen block. **Does not close** #598, #593, #596, #566.
**Ticket:** [#601](https://github.com/LightChainr/Matching-One/issues/601). Distinct from [#592](https://github.com/LightChainr/Matching-One/issues/592) (#588 literature: Mori–Zwanzig / positive realization / finite-horizon balancing). No source below answers both.
**Access:** abstracts, exact statements, and full bibliographic detail are the payload. PDF URLs are secondary.

The write-up currently treated as ours is commit `bdf1339f`, `notes/p398-reflection-parity-20260906.md` (not on `main`). This note says, for each of Q1–Q5, **cite or claim**. A negative answer is a result.

---

## Decision table (cite vs claim)

| Q | Verdict | Cite | What remains ours |
|---|---|---|---|
| **Q1** | **CITE.** The fixed-point count is classical (2005). The orbit formula is Burnside. | Callan–Smiley, arXiv:math/0510447, **Theorem 1**. Ding 2016, Theorem 2.1.2 / Kreweras anti-isomorphism of the two even-`n` classes. OEIS [A007123](https://oeis.org/A007123) (Callan, 8 Oct 2005). | The *partition* equality `coarsest D0-admissible lumping = R-orbit partition` at five widths is a P398 fact, not a combinatorial theorem. The identification of *which* `R` is the task (`wrap`). |
| **Q2** | **NOT a theorem** without extra hypotheses. Aut-orbit is always *a* lumping. Coarsest = Aut-orbit is not. | Kemeny–Snell Thm 6.3.2 (criterion). Buchholz 1994 (exact vs ordinary). D'Angeli–Donno 2013 **Prop. 13** (counterexample) and **Thm 12** (when it *is* a theorem). Godsil–Royle equitable partitions. | Gate 1 branch A at widths 4–8 is an *instance*. "Barrett & Feng" was not found as a lumpability paper. |
| **Q3** | Abstract vanishing is Schur / Wigner–Eckart. **No named Markov-generator selection rule found.** The gap is real. | Wigner–Eckart (Hamiltonian / tensor operators). Hänggi 1978 I (symmetries of the master equation) and II (Markov linear response = FDR, not this). Golubitsky–Stewart–Schaeffer is bifurcation, not this. | The C2 verification (integrand `5.5e-16`; exposed channel `halves_linked` at odd `w`) is the instance. Stating the lemma for equivariant Markov generators, with the Duhamel integrand vanishing pointwise, is a small genuine contribution of language, not of mathematics. |
| **Q4** | Fieller *sample-size for a ratio of means* is standard. Choosing among *observables* by Fieller / projective criteria is **not found**. Shared-denominator-across-folds is **not found** as a named remedy. | Fieller 1940/1954. Hauschke–Kieser–Diletti–Burke, *Stat. Med.* 18 (1999) 93–105. Kieser–Hauschke 1999. Gleser–Hwang 1987; Koschat 1987; von Luxburg–Franz 2009. PowerTOST `sampleN.RatioF`. | Gate 2's scoring of candidate observables, and the shared-denominator fold remedy, stay ours as method. Adjacent to pooled-control in bioequivalence, not identical. |
| **Q5** | **No named "non-monotone prefix" pathology.** The phenomenon is the known non-monotonicity of interpolatory / Galerkin error off complete Krylov grades. Standard remedy is: only truncate at completed block boundaries. | Grimme 1997 (moment matching ⇔ full Krylov of that grade). Gugercin–Antoulas–Beattie, *SIAM J. Matrix Anal. Appl.* 30 (2008) (IRKA is not descent; H2 error can increase). Beattie–Gugercin (same). Saad (no residual monotonicity for Arnoldi). | The width-5–8 observation (rank 7 worse than 6, etc., while declared readouts stay monotone) is an instance, not a new named fact. |

---

## Q1 — Noncrossing partitions fixed by a reflection

### The count, exactly

P398 observes, at widths 4–8, that every reflection of the `w`-cycle fixes `C(w, ⌊w/2⌋)` noncrossing states, hence by Burnside for the group `{1, R}`

```text
#orbits under a single reflection
  = [ Catalan(w) + C(w, ⌊w/2⌋) ] / 2
  = 10, 26, 76, 232, 750           for w = 4..8.
```

Predictions in #593: `r_positive(9) = 2494`, `r_positive(10) = 8524`.

This is **Callan–Smiley Theorem 1 plus Burnside**, not a discovery.

### Primary source

David Callan and Len Smiley, *Noncrossing Partitions Under Rotation and Reflection*, arXiv:math/0510447, 30 October 2005. Not subsequently published as a journal article (no published version found).

> **Abstract.** We consider noncrossing partitions of `[n]` under the action of (i) the reflection group (of order 2), (ii) the rotation group (cyclic of order `n`) and (iii) the rotation/reflection group (dihedral of order `2n`). First, we exhibit a bijection from rotation classes to bicolored plane trees on `n` edges, and consider its implications. Then we count noncrossing partitions of `[n]` invariant under reflection and show that, somewhat surprisingly, they are equinumerous with rotation classes invariant under reflection.

**Definition.** The complement of a partition `π` of `[n]` is `C(π) := n+1 − π` (elementwise). `π` is *self-complementary* if `C(π) = π`.

> Clearly, a NC partition is self-complementary if its labeled polygon diagram is invariant when flipped across a vertical line.

> **Theorem 1.** The number of self-complementary NC partitions of `[n]` is `|𝒜ₙ| = \binom{n}{\lfloor n/2\rfloor}`.

Proof (even case `n = 2m`, quoted in outline): a self-complementary NC partition of `[2m]` corresponds to a NC partition of `[m]` with maximal blocks optionally marked; via the Dyck-path correspondence these are Dyck `m`-paths with returns available for marking; flipping marked components yields a balanced `m`-path, counted by `\binom{2m}{m}`. The odd case is "similar and is omitted."

> As shown in Section 4, **the central binomial coefficient counts those NC Partition Patterns invariant by any reflection across a diameter which fixes the `n` points.** Cauchy–Frobenius allows us to count chirally inequivalent patterns by adding half the non-invariant patterns to the invariant ones.

That last sentence is the orbit formula: `#orbits = [Catalan(n) + \binom{n}{\lfloor n/2\rfloor}] / 2`.

### Vertex-axis vs edge-axis (the fact that makes the orbit count non-identifying)

For even `n` the dihedral group `I₂(n)` has **two non-conjugate reflection classes**:

- vertex-axis: reflecting line through two opposite vertices;
- edge-axis: reflecting line bisecting two opposite sides.

(Rao–Suk 2019, and Ding 2016, Theorem 1.2.2; see below.) For odd `n` there is one class.

Callan–Smiley's involution `C` is one specific reflection (the vertical flip of the labelled `[n]`): vertex-axis when `n` is odd, edge-axis when `n` is even. Their wording "any reflection across a diameter which fixes the `n` points" claims the count for every such reflection. What *proves* the two even-`n` classes have the same cardinality is Ding:

Ziqian Ding, *Dihedral Symmetries of Non-crossing Partition Lattices*, Ph.D. dissertation, University of Miami, August 2016. Advisor: Drew Armstrong. 68 pages. Open access: [scholarship.miami.edu, 991031447399702976](https://scholarship.miami.edu/esploro/outputs/doctoral/Dihedral-Symmetries-of-Non-crossing-Partition-Lattices/991031447399702976).

> We start from enumerative properties of the lattice `NC(n)^F`. Next we investigate the recursive structure on the lattice `NC(n)^F` related to central binomial coefficients and the Catalan numbers.

> When `n` is even the reflections fall into two conjugacy classes: one with `F` : `{F, R²F, R⁴F, …}`; one with `RF` : `{RF, R³F, R⁵F, …}`. **We will show that the lattices `NC(n)^F` and `NC(n)^{RF}` are anti-isomorphic under a map which is called the Kreweras Complement and hence we only need to study the behavior of `NC(n)^F`.**

> We will first investigate the enumerative properties of `NC(n)^F` by establishing a nice bijection between `NC(n)^F` and `NC(n)^{R^{⌊n/2⌋}}` which tells us that **the number of `NC(n)^F` is just equal to the central binomial coefficient `\binom{n}{\lfloor n/2\rfloor}`.**

> **Theorem 2.1.2** ([3, Theorem 1]). The number of self-complementary non-crossing partitions of `[n]` is `\binom{n}{\lfloor n/2\rfloor}`. A non-crossing partition `π ∈ NC(n)` is called self-complementary if `F(π) = π`, which is clearly equivalent to say that `π ∈ NC(n)^F`.

So: Callan–Smiley give the count for one class; Ding gives the Kreweras anti-isomorphism of the two even-`n` classes to each other (and a bijection of either to the 180°-rotation-fixed sublattice, itself counted by the same central binomial via the `q = −1` Catalan). **All `w` reflections give `C(w, ⌊w/2⌋)`.** That is a theorem, not only P398 data. It is exactly why agreement with `r_positive` "identifies nothing on its own" (#598 write-up, already correct as a caveat; the missing piece was the citation).

### What cyclic sieving is, and is not

Victor Reiner, Dennis Stanton and Dennis White, *The cyclic sieving phenomenon*, *J. Combin. Theory Ser. A* **108** (2004), no. 1, 17–50. DOI [10.1016/j.jcta.2004.04.009](https://doi.org/10.1016/j.jcta.2004.04.009).

> The cyclic sieving phenomenon is defined for generating functions of a set affording a cyclic group action, generalizing Stembridge's `q = −1` phenomenon. The phenomenon is shown to appear in various situations, involving `q`-binomial coefficients, Pólya–Redfield theory, polygon dissections, **noncrossing partitions**, finite reflection groups, and some finite field `q`-analogues.

CSP for `NC(n)` is the **rotation** action, with the `q`-Catalan. The `q = −1` evaluation recovers the 180°-rotation-fixed count `Cat_n(−1) = \binom{n}{n/2}` (even `n`), which is Ding's `NC(n)^{R^{n/2}}`, *not* a reflection-fixed count. **Do not cite RSW 2004 as the reflection theorem.**

Sujit Rao and Joe Suk, *Dihedral sieving phenomena*, arXiv:1710.06517, v3 8 March 2019.

> **Proposition 5.1.** Let `n` be odd and `X = {non-crossing partitions of [n]}`. Then the triple `(X, 1/[n+1]_{s,t} \binom{2n}{n}_{s,t}, I₂(n))` exhibits dihedral sieving.

They treat **odd `n` only** (one reflection class). The proof cites Ding Theorem 2.1.5 for the reflection character value being the central binomial. For even `n` they remark that the two non-conjugate reflections `s` and `rs` require a different representation ring (their `χ_b`). Dihedral sieving is a generating-function lift of the same count; it is not needed to cite the count.

### Type B is a different object

Vic Reiner, *Non-crossing partitions for classical reflection groups*, *Discrete Math.* **178** (1998), 229–250 (often cited as 1997). Type-B noncrossing partitions of `[n]` are counted by `\binom{2n}{n}`, realized as **centrally symmetric** NC partitions of `2n` points (equivalently, 180° rotation of a `2n`-gon). That is `NC(2n)^{R^n}`, not `NC(n)^F`. Athanasiadis, *Electron. J. Combin.* **5** (1998), R42, Theorem 2.1: `|NC(B_n)| = \binom{2n}{n}`. **Do not cite type B for this count.**

Simion–Ullman, *On the structure of the lattice of noncrossing partitions*, *Discrete Math.* **98** (1991), 193–206, is the symmetric-chain decomposition of `NC(n)`, not the reflection-fixed sublattice.

### OEIS

[A007123](https://oeis.org/A007123): `1, 1, 2, 4, 10, 26, 76, 232, 750, 2494, 8524, 29624, …`

Callan, 8 Oct 2005 (same week as the preprint):

> `a(n)` = number of noncrossing set partitions of `[n]` up to reflection (`i ↔ n+1−i`).

Offset: the terms `10, 26, 76, 232, 750, 2494, 8524` are `a(5)` through `a(11)`, i.e. NC of `[4]` through `[10]` up to one reflection. So

```text
r_positive(w)  =  A007123(w+1)     (w = 4..10).
```

`r_positive(9) = 2494` and `r_positive(10) = 8524` are therefore **checks against a 2005 sequence**, not against an in-house extrapolation.

### What to write instead of the #598 sentence

The sentence

> the reported positive block counts match a single-reflection orbit formula

becomes

> the reported positive block counts are the Burnside orbit counts of `NC(w)` under a single reflection, equal to `[Cat_w + \binom{w}{\lfloor w/2\rfloor}]/2` by Callan–Smiley, Theorem 1, and (for the even-`w` second conjugacy class) Ding's Kreweras anti-isomorphism `NC(w)^F ≅ NC(w)^{RF}`. This is OEIS A007123(`w+1`). Cardinality does not identify `R`; the task does.

The partition equality `coarsest D0-admissible lumping = this orbit partition`, block for block, at five widths, remains a P398 computation (#598 / Gate 1 branch A). Combinatorics does not give that.

---

## Q2 — Symmetry and exact lumpability

### What is a theorem

**Aut-orbit partitions of a chain automorphism are always strong lumpings.** If `G` is a group of bijections of the state space with `p(gx, gy) = p(x, y)` for all `g, x, y` (equivalently, the generator is `G`-equivariant), then the orbit partition is strongly lumpable. This is the content of Kemeny–Snell's criterion applied to orbits, and is the Markov-chain case of "the orbit partition of a group of automorphisms is equitable" (Godsil–Royle).

John G. Kemeny and J. Laurie Snell, *Finite Markov Chains*, Van Nostrand, 1960; Springer reprint 1976. **Theorem 6.3.2:**

> A Markov chain is lumpable with respect to the partition if and only if the distribution of the next subset is the same for every state in the present subset.

Peter Buchholz, *Exact and ordinary lumpability in finite Markov chains*, *J. Appl. Probab.* **31** (1994), no. 1, 59–75. DOI [10.2307/3215235](https://doi.org/10.2307/3215235). MR 1260571.

> Exact and ordinary lumpability in finite Markov chains is considered. Both concepts naturally define an aggregation of the Markov chain yielding an aggregated chain that allows the exact determination of several stationary and transient results for the original chain.

Buchholz's *ordinary* lumpability is Kemeny–Snell strong lumpability (the lumped process is Markov for every initial distribution). His *exact* lumpability is the dual (preserves occupancy / backwards). P398's "exact strong lumping" is ordinary/strong lumpability in this language. Buchholz does **not** prove that the coarsest lumping equals Aut orbits.

Chris Godsil and Gordon Royle, *Algebraic Graph Theory*, GTM 207, Springer, 2001, Chapter 9. The orbit partition of `Aut(Γ)` is always equitable. The converse is false: a regular graph that is not vertex-transitive has the trivial partition `{V}` as an equitable partition strictly coarser than the Aut-orbit partition. (The random-walk chain on such a graph therefore has a lumping coarser than Aut orbits; the trivial partition is degenerate as a "coarsest lumping" of an irreducible chain, but the same distinction — equitable vs orbit — is the one #596 Gate 1 names as A vs B.)

### What is not a theorem

**"The coarsest observable-admissible strong lumping equals the Aut-orbit partition"** is not a theorem of lumpability. It fails as soon as the chain has an equitable partition coarser than Aut orbits that still refines the readout colouring.

Named counterexample, not a graph-theory footnote:

Daniele D'Angeli and Alfredo Donno, *The lumpability property for a family of Markov chains on poset block structures*, *Adv. in Appl. Math.* **51** (2013), no. 3, 367–391. DOI [10.1016/j.aam.2013.04.007](https://doi.org/10.1016/j.aam.2013.04.007). arXiv:1304.4180.

> We prove that, for such a product, every lumping can be obtained from the action of a suitable subgroup of the generalized wreath product of symmetric groups, acting on the underlying poset block structure, **if and only if** the poset defining the Markov process is totally ordered, and one takes the uniform Markov operator in each factor state space.

> **Theorem 12.** Let `P_I` be the Insect Markov chain on the boundary `X^n` of `T_{q,n}`. Let `ℒ = {L_1, …, L_k}` be a lumping partition for `P_I`. Then there exists `K ≤ Aut(T_{q,n})` such that the orbit partition of `X^n` under the action of `K` is `⊔_{i=1}^k L_i`.

> **Proposition 13.** Let `(I, ⪯)` be a finite poset, with `|I| = n`, and suppose that there exist `i, j ∈ I` such that `i ⊀ j` and `j ⊀ i`. Let `P_I` be the Insect Markov chain on `X^n` associated with `(I, ⪯)` and let `F_I` be the corresponding automorphism group. **Then there exists a lumping of `P_I` which is not induced by any subgroup of `F_I`.**

Theorem 12 is the "coarsest = Aut-orbit" statement **with hypotheses** (the poset is a chain, uniform factors). Proposition 13 is the counterexample when those hypotheses fail. P398 is not an Insect chain on a poset block structure; the citation is for the logical shape, not for an identification.

### Generic rates

If the rate matrix is constrained only by the linear equations forced by a declared automorphism group `K`, then for a Zariski-open set of such rates the *only* strong lumpings are the `K`-orbit partitions (and their coarsenings that remain `K`-unions — which, if `K` is the full Aut of a generic such matrix, is just the orbit partition). I did not find this written as a named theorem in the lumpability literature. It is the reason Gate 1's branch A can be true for P398 without being true for Markov chains in general: `{J, D}` on NC states is a very special, highly symmetric generator, not a generic `K`-equivariant matrix.

### "Barrett & Feng"

No standard lumpability paper by Barrett & Feng was found. Nearby names: Jernigan–Baran, *Statist. Probab. Lett.* (2003), testing lumpability (Kemeny–Snell); various Barrett papers on bisimulation of labelled transition systems; Feng on unrelated Markov topics. Treat the name as a false lead unless a specific reference is supplied.

### Bisimulation / coalgebra

Ordinary lumpability = strong bisimulation of Markov chains (standard; e.g. Buchholz, and Bernardo 2008: "Markovian bisimilarity corresponds to ordinary lumpability"). Coalgebraic treatments (e.g. Sokolova) recast the same criterion. They do not add a theorem that coarsest bisimulation = Aut orbits.

### What this does to Gate 1

#596 Gate 1 distinguished

```text
A. symmetry orbit:     P_JD = P_orbit
B. common equitable / coherent structure beyond symmetry
C. accidental cancellation of the pencil
```

#598 resolved A at widths 4–8, *block for block*. That is an instance. It is not D'Angeli–Donno Theorem 12, and it is not a general fact about coarsest lumpings. The honest sentence is:

> For this generator and this dictionary, the coarsest D0-admissible strong lumping *is* the orbit partition of the one reflection the dictionary preserves. Aut-orbit partitions are always lumpings; the converse is not a theorem (D'Angeli–Donno, Prop. 13). Five widths do not make it one.

A generic-rates lemma, if someone wants to write it, would turn "instance" into "theorem with hypotheses." That is optional theory, not a citation.

---

## Q3 — First-order response selection rules for equivariant generators

### The statement

Baseline generator `G` is `K`-equivariant, source `B` and readout `C` transform in irreps `ρ_B`, `ρ_C`, perturbation in `ρ_H`. The Duhamel integrand

```text
C e^{(t−s)G} H e^{sG} B
```

vanishes (hence so does the first-order response) unless the trivial representation of `K` occurs in the relevant tensor product (`ρ_C^* ⊗ ρ_H ⊗ ρ_B`, up to duals). For `K = C₂` this is parity multiplication.

### What this is, classically

1. **Duhamel / variation of constants** for `d/dε exp(t(G+εH))` at `ε=0` is undergraduate semigroup theory. Not a contribution.

2. **Schur's lemma.** `Hom_K(ρ_C, ρ_H ⊗ ρ_B)` is zero unless the trivial appears in `ρ_C^* ⊗ ρ_H ⊗ ρ_B`. Equivalently, an invariant linear functional of an equivariant composition that lands in a non-trivial isotypic component is zero. This is the whole representation-theoretic content.

3. **Wigner–Eckart.** The physics name for matrix elements of a tensor operator between irreps: `⟨α' j' m'| T^k_q |α j m⟩` vanishes unless `m' = m+q` and `j'` occurs in `j ⊗ k`. Stated for Hamiltonian / unitary quantum mechanics, angular momentum. Standard textbook: Messiah, *Quantum Mechanics*, Vol. II, Ch. XIII; or Edmonds, *Angular Momentum in Quantum Mechanics*. **Not stated for Markov generators.**

4. **Hänggi 1978 I.** P. Hänggi, *Stochastic processes I: Asymptotic behaviour and symmetries*, *Helv. Phys. Acta* **51** (1978), 183–201.

   > Furthermore, we investigate the symmetry properties of stochastic processes. We discuss the consequences for stochastic processes of both: symmetry transformations in state space and symmetry properties obtained by interchanging the time arguments in the joint-probability (generalized detailed balance). Various symmetry relations for multivariate probabilities and multi-time correlation functions are obtained.

   State-space symmetries of the master equation, and generalized detailed balance. **No selection rule for the first-order response of an equivariant generator to a symmetry-breaking perturbation.**

5. **Hänggi 1978 II.** P. Hänggi, *Stochastic processes II: Response theory and fluctuation theorems*, *Helv. Phys. Acta* **51** (1978), 202–219.

   > Linear and nonlinear response theory are developed for stationary Markov systems describing systems in equilibrium and nonequilibrium. Generalized fluctuation theorems are derived which relate the response function to a correlation of nonlinear fluctuations of the unperturbed stationary process.

   This is FDR for Markov processes: response related to unperturbed correlations. It is the *existence* of a linear response, not a representation-theoretic vanishing of one.

6. **Golubitsky–Stewart–Schaeffer**, *Singularities and Groups in Bifurcation Theory*, Vol. II, Springer 1988. The equivariant branching lemma is about bifurcating equilibria of equivariant ODEs (`f(x, λ) = 0` with `Γ`-equivariance). It is not a linear-response statement for `C exp(t(G+εH)) B`.

7. **Antown–Dragičević–Froyland**, *Optimal linear responses for Markov chains and stochastically perturbed dynamical systems*, *J. Stat. Phys.* **170** (2018), 1051–1087. arXiv:1801.03234. They *choose* a perturbation to *maximise* the linear response of an observable. Opposite direction: they assume the response is available and optimize it. No vanishing rule.

8. Recent transfer-operator linear response (e.g. arXiv:2606.02889, "Ulam Approximation for Nonautonomous Systems: Equivariant Measures and Linear Response") uses "equivariant" for time-dependent families, not for a group action on the state space.

### The gap

I did not find a citable statement of the form:

> Let `G` be the generator of a Markov process, `K`-equivariant. Let `H` transform in a non-trivial irrep of `K`, and let `B`, `C` be `K`-invariant. Then the first-order response `d/dε[C exp(t(G+εH)) B]|_{ε=0}` vanishes identically, and the Duhamel integrand vanishes pointwise.

The mathematics is Schur's lemma applied to Duhamel; anyone who writes it down has it. The **gap is that no one states this cleanly for Markov generators**. That is a useful negative. It makes the #598 lemma a small genuine contribution of *language and of the pointwise-integrand check*, not of representation theory.

What is *not* a contribution: the abstract vanishing. What *is* P398's: the verification that the integrand (not just the integral) is `~1e-15` on protected pairs, and that the complementary channel `halves_linked` at odd width fires at `0.292` / `0.107`. A run in which everything vanished would have been a broken test; the write-up already says this.

Cite Schur / Wigner–Eckart for the vanishing; cite Hänggi II if a Markov-response authority is wanted (for the *existence* of the derivative, not for the selection rule); claim the lemma-for-generators as a remark, not as a theorem of 1978.

---

## Q4 — Fieller used for experiment design rather than inference

### What is standard (sample size for a ratio of means)

E. C. Fieller, *The biological standardization of insulin*, *J. Roy. Statist. Soc. Suppl.* **7** (1940), 1–64.

E. C. Fieller, *Some problems in interval estimation*, *J. Roy. Statist. Soc. Ser. B* **16** (1954), 175–185. The `(1−α)` set for `μ_a/μ_b` is bounded / two unbounded rays / the whole line according as `g < 1`, `g = 1`, `g > 1`, with `g = t^2 s^2 ν_{22} / b^2`.

Dieter Hauschke, Meinhard Kieser, Edgar Diletti and Martin Burke, *Sample size determination for proving equivalence based on the ratio of two means for normally distributed data*, *Statistics in Medicine* **18** (1999), 93–105. DOI [10.1002/(SICI)1097-0258(19990115)18:1<93::AID-SIM992>3.0.CO;2-8](https://doi.org/10.1002/(SICI)1097-0258(19990115)18:1<93::AID-SIM992>3.0.CO;2-8).

> We consider the problem when equivalence is defined in terms of the ratio of population means and the original (untransformed) data are normally distributed. Application of the intersection-union principle to the test proposed by Sasabuchi results in a two one-sided tests procedure of size `α`. We give the associated `100(1−2α)%` confidence interval and **derive the exact methods for calculation of power and sample sizes** for the parallel group design and the two-period cross-over.

Meinhard Kieser and Dieter Hauschke, *Approximate sample sizes for testing hypotheses about the ratio and difference of two means*, *J. Biopharm. Statist.* **9** (1999), 641–650. PMID 10576408.

> A simple approximation is given to the sample size required for testing hypotheses about the ratio of the means. The formula includes the situations of testing noninferiority, superiority, or equivalence.

This is the canonical bioequivalence / relative-potency home. R package **PowerTOST**, function `sampleN.RatioF`, implements Fieller-based sample size for a ratio of means (Labes, Schütz, Lang). Gate 2's "smallest future sample multiplier for which an `α`-level Fieller set lies inside a declared window" is this object, applied to a different functional.

### Unbounded CIs (the named failure mode, but not the one asked)

Leon Jay Gleser and Jiunn T. Hwang, *The nonexistence of `100(1−α)%` confidence sets of finite expected diameter in errors-in-variables and related models*, *Ann. Statist.* **15** (1987), 1351–1362.

> It is impossible to construct confidence intervals for key parameters which have both positive confidence and finite expected length.

Martin A. Koschat, *A characterization of the Fieller solution*, *Ann. Statist.* **15** (1987), 462–468. Fieller is the unique exact-coverage procedure in a large class; no procedure gives bounded `α`-level intervals with probability 1.

Ulrike von Luxburg and Volker H. Franz, *A geometric approach to confidence sets for ratios: Fieller's theorem, generalizations, and bootstrap*, *Statistica Sinica* **19** (2009), 1495–1517. arXiv:0711.0198.

> Different researchers (Gleser and Hwang, 1987; Koschat, 1987; Hwang, 1995) have shown that any method which is not able to generate such unbounded confidence limits for a ratio leads to arbitrary large deviations from the intended confidence level.

When the denominator's CI contains 0 (`g ≥ 1`), the Fieller set is unbounded. That is the named pathology of *inference* for a ratio. It is adjacent to Gate 2's problem but not the same: Gate 2's validation statistic is itself a ratio whose *denominator is estimated with a weak sample*, so the *design cost* is unstable under cross-validation.

### What was not found

1. **Choosing among candidate *observables* by Fieller / projective criteria.** Particle-physics "optimal observable" (Gail–Wedi, or #370's idea) maximises a different functional (typically `s/√b` or Fisher information for a coupling). I did not find a literature that ranks observables by "how soon a Fieller set for *this* ratio falls inside a window." Gate 2's scoring rule stays ours as method.

2. **The validation-statistic-is-itself-a-ratio-with-a-weak-denominator, and the shared-denominator-across-folds remedy.** No named paper. Adjacent: pooled control in bioequivalence (share the control-arm estimate across comparisons); Fieller with a common denominator in relative potency when several treatments share a vehicle. Neither is a CV-of-a-ratio-design-statistic discussion, and neither is "share the denominator estimate across CV folds." The improvised remedy is reasonable and should be described as such, not as a rediscovery.

---

## Q5 — Non-monotone block-Krylov prefixes

### What is known

Moment matching holds when the projection subspace contains the *full* Krylov subspace of the corresponding grade. An incomplete block (a prefix that stops inside a block level) does not match the moments associated with that level, and the Galerkin / Petrov–Galerkin closure on that prefix is not an interpolant of that grade.

Eric J. Grimme, *Krylov projection methods for model reduction*, Ph.D. thesis, University of Illinois at Urbana-Champaign, 1997.

> Projections onto unions of Krylov subspaces lead to a class of reduced-order models known as rational interpolants. The cornerstone of this dissertation is a collection of theory relating Krylov projection to rational interpolation.

(Grimme's implicit moment-matching theorem: if `range(V)` contains the Krylov subspace of grade `q` at a shift, the reduced transfer function matches `q` moments there. A prefix that is not a Krylov subspace of any grade matches correspondingly fewer, and in the block case, a cut inside a block is not a completed grade.)

Roland W. Freund, *Krylov-subspace methods for reduced-order modeling in circuit simulation*, *J. Comput. Appl. Math.* **123** (2000), 395–421. Block Krylov for multiple starting vectors: the proper definition is the deflated block sequence; incomplete blocks are not the object the theory is written for.

Serkan Gugercin, Athanasios C. Antoulas and Christopher A. Beattie, *H₂ model reduction for large-scale linear dynamical systems*, *SIAM J. Matrix Anal. Appl.* **30** (2008), 609–638.

IRKA is **not a descent method**. Beattie–Gugercin (later surveys, e.g. arXiv:1409.2140):

> IRKA is not a descent algorithm; that is, the H₂ error might fluctuate during intermediate steps and premature termination of the algorithm could result (at least in principle) in a worse approximation than what was provided for initialization.

The H₂ interpolatory error is **not monotone in reduced order**. Truncating at an intermediate order can be worse than a smaller interpolant that *does* sit at a completed Krylov grade.

Yousef Saad, *Iterative Methods for Sparse Linear Systems*, SIAM, 2nd ed. 2003. Arnoldi residual norms are not monotone (GMRES minimises, Arnoldi does not). No guarantee that a longer Krylov prefix has smaller Galerkin error on a *declared* observable that is not the residual.

### The named-pathology question

I did not find a named "non-monotone block-Krylov prefix" pathology. The reading in #601 —

> a rank cutting *inside* a Krylov block level adds half a level and the Galerkin closure spends it on the readouts it was built from

— is the standard incomplete-grade phenomenon, described in the language of this repository. The standard remedy **is** "only truncate at completed block boundaries" (and, for H₂-optimal reduction, only at orders that satisfy interpolation conditions, not at arbitrary prefixes). Citing Grimme + Gugercin–Antoulas–Beattie for "error is not monotone off complete grades; truncate at block boundaries" is the honest move. Describing the width-5–8 ladder (rank 7 worse than 6, 9 worse than 8, 10 worse than 9, declared readouts monotone) as a *new* pathology would be a rediscovery of the incomplete-grade fact.

---

## What this does to neighbouring tickets

**#598.** The orbit-count sentence is a citation of Callan–Smiley Theorem 1 + Burnside + Ding's even-`n` anti-isomorphism. The partition equality and the selection-rule *verification* stay. Phase C (factor the balanced realization through the quotient) is still the open half of that ticket; nothing here closes it. The C2 lemma may be written as a remark citing Schur, not as a new theorem.

**#593.** `r_positive(9) = 2494` and `r_positive(10) = 8524` are OEIS A007123(10) and A007123(11). If DeepSeek's exact refinement returns those *and* the partition equality (not just the cardinality), the orbit account is confirmed out of sample as combinatorics. If the *cardinality* fails, Callan–Smiley is not at issue — the implementation is. If the cardinality matches and the *partition* fails, Gate 1 branch A is width-limited and the state semantics need repair. The memory-degree half of #593 is untouched by this literature (and is the half still worth computing).

**#596 Gate 1.** Resolved on A at five widths, as an instance. Do not write "by the theorem relating coarsest lumpings to Aut orbits." There is no such theorem without hypotheses; the hypotheses that make it true (D'Angeli–Donno Thm 12) do not hold of P398.

**#594 Q1.** Already rewritten by #598 as "the closed condition is a symmetry." This packet adds the citation for the orbit *count*. The codimension / lumpability-preserving-variety computation remains optional general theory, not needed for P398.

**#592.** No overlap. Do not dump these sources there.

**#370.** Optimal-observable literature remains adjacent to Q4 item 3 and is not this.

---

## Subsequent analysis (ordered)

1. **Rewrite the #598 orbit-count sentence as a citation.** Callan–Smiley Thm 1, Ding 2.1.2 / Kreweras anti-isomorphism, OEIS A007123(`w+1`). Keep the partition equality as the P398 fact.
2. **#593, if it runs, reports partition equality at 9 and 10, not just `r_positive`.** Cardinality is combinatorics; the partition is Gate 1.
3. **Do not claim a general "coarsest lumping = Aut orbits" theorem.** If a generic-rates lemma is wanted, write it; do not cite Buchholz or Barrett & Feng for it.
4. **Write the C2 selection rule as a remark citing Schur / Wigner–Eckart**, with the gap ("not found for Markov generators") stated. Do not cite Hänggi II as the selection rule.
5. **Gate 2 stays a method contribution** on the observable-ranking side. Cite Hauschke et al. 1999 for the sample-size-for-a-ratio-of-means half. Describe the shared-denominator fold as an improvisation.
6. **Krylov prefixes: truncate at completed block boundaries.** Cite Grimme + Gugercin–Antoulas–Beattie for non-monotonicity off complete grades. Do not name a new pathology.

---

## Bibliographic list (full)

- Antown, F., Dragičević, D. and Froyland, G. (2018). Optimal linear responses for Markov chains and stochastically perturbed dynamical systems. *J. Stat. Phys.* **170**, 1051–1087. arXiv:1801.03234.
- Athanasiadis, C. A. (1998). On noncrossing and nonnesting partitions for classical reflection groups. *Electron. J. Combin.* **5**, R42.
- Beattie, C. A. and Gugercin, S. Model reduction by rational interpolation. In *Model Reduction and Approximation* (Benner et al., eds), SIAM, 2017. arXiv:1409.2140.
- Buchholz, P. (1994). Exact and ordinary lumpability in finite Markov chains. *J. Appl. Probab.* **31**, no. 1, 59–75. DOI 10.2307/3215235. MR 1260571.
- Callan, D. and Smiley, L. (2005). Noncrossing partitions under rotation and reflection. arXiv:math/0510447, 30 Oct 2005.
- D'Angeli, D. and Donno, A. (2013). The lumpability property for a family of Markov chains on poset block structures. *Adv. in Appl. Math.* **51**, no. 3, 367–391. DOI 10.1016/j.aam.2013.04.007. arXiv:1304.4180.
- Ding, Z. (2016). *Dihedral Symmetries of Non-crossing Partition Lattices*. Ph.D. dissertation, University of Miami, August 2016. Advisor D. Armstrong. https://scholarship.miami.edu/esploro/outputs/doctoral/Dihedral-Symmetries-of-Non-crossing-Partition-Lattices/991031447399702976
- Fieller, E. C. (1940). The biological standardization of insulin. *J. Roy. Statist. Soc. Suppl.* **7**, 1–64.
- Fieller, E. C. (1954). Some problems in interval estimation. *J. Roy. Statist. Soc. Ser. B* **16**, 175–185.
- Freund, R. W. (2000). Krylov-subspace methods for reduced-order modeling in circuit simulation. *J. Comput. Appl. Math.* **123**, 395–421.
- Gleser, L. J. and Hwang, J. T. (1987). The nonexistence of 100(1−α)% confidence sets of finite expected diameter in errors-in-variables and related models. *Ann. Statist.* **15**, 1351–1362.
- Godsil, C. and Royle, G. (2001). *Algebraic Graph Theory*. GTM 207, Springer. Chapter 9.
- Golubitsky, M., Stewart, I. and Schaeffer, D. G. (1988). *Singularities and Groups in Bifurcation Theory*, Vol. II. Springer.
- Grimme, E. J. (1997). *Krylov projection methods for model reduction*. Ph.D. thesis, University of Illinois at Urbana-Champaign.
- Gugercin, S., Antoulas, A. C. and Beattie, C. A. (2008). H₂ model reduction for large-scale linear dynamical systems. *SIAM J. Matrix Anal. Appl.* **30**, 609–638.
- Hänggi, P. (1978). Stochastic processes I: Asymptotic behaviour and symmetries. *Helv. Phys. Acta* **51**, 183–201.
- Hänggi, P. (1978). Stochastic processes II: Response theory and fluctuation theorems. *Helv. Phys. Acta* **51**, 202–219.
- Hauschke, D., Kieser, M., Diletti, E. and Burke, M. (1999). Sample size determination for proving equivalence based on the ratio of two means for normally distributed data. *Stat. Med.* **18**, 93–105.
- Kemeny, J. G. and Snell, J. L. (1960). *Finite Markov Chains*. Van Nostrand; Springer reprint 1976. Theorem 6.3.2.
- Kieser, M. and Hauschke, D. (1999). Approximate sample sizes for testing hypotheses about the ratio and difference of two means. *J. Biopharm. Statist.* **9**, 641–650.
- Koschat, M. A. (1987). A characterization of the Fieller solution. *Ann. Statist.* **15**, 462–468.
- Rao, S. and Suk, J. (2019). Dihedral sieving phenomena. arXiv:1710.06517v3, 8 Mar 2019.
- Reiner, V. (1997/98). Non-crossing partitions for classical reflection groups. *Discrete Math.* **178**, 229–250.
- Reiner, V., Stanton, D. and White, D. (2004). The cyclic sieving phenomenon. *J. Combin. Theory Ser. A* **108**, no. 1, 17–50.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*. 2nd ed., SIAM.
- Simion, R. and Ullman, D. (1991). On the structure of the lattice of noncrossing partitions. *Discrete Math.* **98**, 193–206.
- von Luxburg, U. and Franz, V. H. (2009). A geometric approach to confidence sets for ratios: Fieller's theorem, generalizations, and bootstrap. *Statistica Sinica* **19**, 1495–1517. arXiv:0711.0198.
- OEIS A007123. Callan comment, 8 Oct 2005.

---

## Not established, and not claimed here

- a general theorem that coarsest strong lumping = Aut-orbit partition, without extra hypotheses;
- a named Markov-generator selection rule in the 1978 (or later) literature;
- a named method for ranking observables by Fieller / projective design cost;
- a named shared-denominator-across-folds remedy;
- a named "non-monotone block-Krylov prefix" pathology beyond incomplete-grade moment matching;
- type-B NC, RSW cyclic sieving, or Simion–Ullman as the reflection-fixed count;
- anything in the claim ledger;
- anything about square-site `p_c`, matching-odd, or N=580.

This packet does not enter `docs/STATUS.md`.
