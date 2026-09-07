# Literature officer — #601 round 2: neighbourhoods that are not theorems

**Date:** 2026-09-06
**Role:** theory input. Same epistemic level as an Astra answer and as [#572](https://github.com/LightChainr/Matching-One/pull/572). **Does not enter** `docs/STATUS.md`. **Does not score** a frozen block. **Does not close** #598, #593, #596, #566, #601.
**Ticket:** follow-up to [#601](https://github.com/LightChainr/Matching-One/issues/601) / [#602](https://github.com/LightChainr/Matching-One/pull/602). Distinct from [#592](https://github.com/LightChainr/Matching-One/issues/592). No source below answers both.
**Access:** abstracts, exact statements, and full bibliographic detail are the payload. PDF URLs are secondary.

This is **not a recap** of round 1. Round 1 (`notes/literature-officer-20260906-issue601.md`, commit `0c6e6b14`) already cited Callan–Smiley Theorem 1, Ding 2016, D'Angeli–Donno Prop. 13 / Thm 12, Hauschke et al. 1999, Grimme 1997, Gugercin–Antoulas–Beattie 2008. Those sources are not re-quoted here. The brief was: retrieve further related research. Negative answers remain results.

A neighbourhood of a gap is not a theorem filling the gap.

---

## What moved, and what did not

| Q | Round 1 | Round 2 | Still ours |
|---|---|---|---|
| **Q1** | **CITE.** Callan–Smiley Thm 1 + Ding even-`n` anti-isomorphism. OEIS A007123. | Two 2024–2026 NC papers (Dougherty–Root; Ebrahimi-Fard et al.) are **not** competitors for the orbit count. Callan–Smiley remains unpublished as a journal article. | Partition equality at five widths. Identification of *which* `R`. |
| **Q2** | **NOT a theorem.** Aut-orbit always a lumping; coarsest = Aut fails (D'Angeli–Donno Prop. 13). | The classical name of Gate 1 A vs B is **Schurian vs non-Schurian** (Higman / Cameron / Babai) plus Kudose 2009 Prop. 2.4 (McKay graph). Watanabe–Wolfer arXiv:2412.08400v3 (accepted *Information Geometry* 2026) is adjacent to generic-rates and **does not** contain Aut / orbit / Zariski language. Best "Barrett & Feng" candidate is **Barr & Thomas 1977** (Barr, not Barrett; Thomas, not Feng). | Gate 1 A at widths 4–8 is still an instance. Generic-rates lemma still not found named. |
| **Q3** | Abstract vanishing is Schur. **No named Markov-generator selection rule.** | Diaconis 1988 Ch. 3E is the *group-walk* case of the same Schur fact. Santos Gutiérrez–Zagli–Carigi 2025 still *optimise* a Markov linear response (same direction as Antown–Dragičević–Froyland 2018). The Duhamel-integrand form for a general equivariant Markov generator is still unnamed. | The C2 pointwise check. The language of the lemma. |
| **Q4** | Sample-size for a ratio of means is Hauschke et al. 1999. Ranking observables / shared-denominator-across-folds **not found**. | **Sitter–Wu 1993 F-optimality** is the named "Fieller as a design criterion" paper. It chooses **dose-support points of a binary-response experiment**, not among observables. Minkin–Kundhal 1999 is the same island. Shared-denominator-across-folds still not found. | Gate 2's scoring of candidate *observables*, and the shared-denominator fold remedy. |
| **Q5** | No named "non-monotone prefix" pathology. Remedy: truncate at completed block boundaries. | Freund's **exact deflation** of a block-Krylov sequence: a prefix that is not a deflation boundary is not a Krylov subspace of any grade. Restatement of Grimme, not a name. | The width-5–8 observation. |

**#600 Phase C** (not one of Q1–Q5, but the neighbouring compute ticket opened the same hour): the expected result — nonzero task Hankel data identical between the full space and the `R`-orbit quotient; `R`-odd sector exactly uncontrollable and/or unobservable — is the **Kalman decomposition** applied to a `C2`-equivariant realization with even `B` and even `C`. Cite Wonham. Not a new theorem. The odd-width `halves_linked` trap is what makes the test non-vacuous.

---

## Q4 — F-optimality is named, and is the wrong object

Round 1 asked whether Fieller is used for *design* rather than inference, and whether anyone chooses among *observables* by a Fieller / projective criterion. The first half now has a name. The second half still does not.

### The named criterion

Randy R. Sitter and C. F. J. Wu, *Optimal designs for binary response experiments: Fieller, D, and A criteria*, *Scandinavian Journal of Statistics* **20** (1993), no. 4, 329–341. JSTOR [4616288](https://www.jstor.org/stable/4616288). ISSN 0303-6898.

> **Abstract.** For binary response experiments, we consider the use of a confidence interval based on Fieller's theorem as a design criterion (**F-optimality**). For symmetric distributions and under mild conditions, we show that F-optimal, D-optimal, and A-optimal designs have two or three support points. A complete characterization of these designs is given. The possibility of having 4-point symmetric designs is discussed.

What is being designed: the **dose-support points** of a binary-response (quantal) experiment, typically for ED50 / relative potency. The decision variables are locations on the dose axis, not a menu of observables.

Companion paper, same authors, same year, *inference* rather than design:

R. R. Sitter and C. F. J. Wu, *On the accuracy of Fieller intervals for binary response data*, *J. Amer. Statist. Assoc.* **88** (1993), no. 423, 1021–1025. DOI [10.1080/01621459.1993.10476370](https://doi.org/10.1080/01621459.1993.10476370). JSTOR 2290794.

> Finney proposed the use of a fiducial interval for the median response dose based on Fieller's theorem. An alternative is to use the asymptotic confidence interval. … Our results show that Fieller intervals are generally superior.

Same island, still ED50 of a binary-response curve.

Same island, later:

Salomon Minkin and Kiran Kundhal, *Likelihood-based experimental design for estimation of ED50*, *Biometrics* **55** (1999), no. 4, 1030–1037. DOI [10.1111/j.0006-341X.1999.01030.x](https://doi.org/10.1111/j.0006-341X.1999.01030.x).

> In selecting the best dosage choice for the estimation of ED50, it is natural to try to minimize the length of the confidence intervals. In this presentation, the dose allocation that minimizes the length of the likelihood-based confidence intervals is presented and compared with alternative allocations that have been proposed based on the length of different types of confidence intervals, such as those based on the asymptotic variance or on **Fieller's Theorem**.

Again: dose allocation, not observable selection.

V. Guiard, *Some remarks on the estimation of the ratio of the expectation values of a two-dimensional normal random variable (correction of the theorem of Milliken)*, *Biometrical Journal* **31** (1989), no. 6, 681–697. DOI [10.1002/bimj.4710310605](https://doi.org/10.1002/bimj.4710310605). (Fieller region is exact LR; Milliken's "conservative" claim is false.) Cited from the Sitter–Wu trail; still inference for a ratio, not ranking of observables.

### What this is not

Gate 2 (#595 / #596) scores candidate *histogram channels / p / frozen linear combinations* by the future sample multiplier that puts an `α`-level Fieller set for `(A4, delta)` inside a declared contamination window. That is:

```text
choose among observables,   not among dose-support points;
the functional is a ratio of estimated amplitudes, not an ED50.
```

Sitter–Wu F-optimality must be cited as the named "Fieller-as-design-criterion" paper, and immediately distinguished. Citing it as the Gate 2 method would be a misidentification.

### Shared-denominator-across-folds

Still not found as a named remedy. Adjacent objects that are *not* this:

- pooled-control / pooled-variance in bioequivalence (standard, different problem);
- Rubin's pooling of multiply-imputed point estimates across completed datasets (MICE; a CV-folds p-value question on CrossValidated recites this and is not a paper);
- Gleser–Hwang / Koschat / von Luxburg–Franz unbounded CIs when the denominator CI contains 0 (round 1; *inference* pathology, not a CV-design remedy).

#600's second job ("the archived block has 100 batches, so a 33-batch validation fold estimates `A4²` — the denominator of the cost — very poorly. … a proper block bootstrap over batches would replace the workaround with an interval") is the same defect, stated as compute rather than as literature. A block bootstrap of the design cost is a valid replacement for the improvisation; it is not a citation of a named shared-denominator theorem.

**Verdict, upgraded:** F-optimality is named (Sitter–Wu 1993) and is the wrong object. Ranking observables by Fieller, and shared-denominator-across-folds, remain ours as method.

---

## Q2 — Schurian vs non-Schurian is the classical name of Gate 1 A vs B

Round 1: Aut-orbit partitions are always lumpings; coarsest = Aut-orbit is not a theorem (D'Angeli–Donno Prop. 13); it *is* a theorem on a chain/tree (their Thm 12). Gate 1 A at five widths is an instance.

The same distinction is older than lumpability, and has a name.

### Equitable ≠ orbit (the graph-theory sentence)

Satoru Kudose, *Equitable partitions and orbit partitions*, VIGRE REU paper, University of Chicago, 2009. Open: [math.uchicago.edu/~may/VIGRE/VIGRE2009/REUPapers/Kudose.pdf](https://www.math.uchicago.edu/~may/VIGRE/VIGRE2009/REUPapers/Kudose.pdf).

> **Abstract.** We consider two kinds of partition of a graph, namely orbit partitions and equitable partitions. Although an orbit partition is always an equitable partition, the converse is not true in general.

> **Proposition 2.4.** An orbit partition is an equitable partition.

> The graphs in Fig. 2 show that the converse is false in general. The first graph in Fig. 2 is called **McKay's graph** and is probably the most popular graph of this type in the literature. Note that the black cell is not an orbit: no automorphism takes an outer black vertex to an inner black vertex, since automorphisms must preserve cycles. The second graph is 3-regular, so the trivial partition is an equitable partition. However, since it is not vertex-transitive, the trivial partition is not an orbit.

This is Godsil–Royle Chapter 9 in undergraduate language, with a named counterexample. The random-walk chain on McKay's graph therefore has a strong lumping (the colouring) strictly coarser than Aut-orbits. That is Gate 1 B as a *graph*, not as P398.

### Schurian vs non-Schurian (the configuration-theory sentence)

A coherent configuration is **Schurian** if it is the orbital configuration of a permutation group — i.e. if its relations are Aut-orbits on ordered pairs. Otherwise it is non-Schurian. The distinction is Higman's, named after Schur; standard expositions:

- László Babai, unpublished notes on graph isomorphism / quasipolynomial algorithms, §3.1.3. Observation 3.1.12: orbitals of `G ≤ S(Ω)` form a coherent configuration `X(G)`. "We say that a coherent configuration is Schurian if it is the orbital configuration of some permutation group." Remark 3.1.14: "Not every coherent configuration is Schurian; in fact, there are large families of strongly regular graphs with no non-identity automorphisms (line graphs of Steiner triple systems, point graphs of Latin squares)."
- Peter J. Cameron, standard surveys of coherent configurations / association schemes. Same definition.
- D. G. Higman, the original combinatorial relaxation of permutation groups.

Gate 1 of #596 distinguished

```text
A. symmetry orbit:     P_JD = P_orbit
B. common equitable / coherent structure beyond symmetry
C. accidental cancellation of the pencil
```

**A is the Schurian case. B is the non-Schurian case.** Kudose Prop. 2.4 is A ⇒ equitable; McKay's graph is a B that is not an A. D'Angeli–Donno Prop. 13 remains the Markov-chain counterexample (lumpable, not Aut-induced). None of these is a theorem that P398's `P_JD` equals `P_orbit`. Five widths remain five widths.

### Watanabe–Wolfer: adjacent to generic-rates, not the lemma

Shun Watanabe and Geoffrey Wolfer, *Characterization of exponential families of lumpable stochastic matrices*, arXiv:2412.08400v3, 26 November 2025. Tokyo University of Agriculture and Technology. Accepted *Information Geometry*, 22 February 2026; published 6 April 2026. DOI [10.1007/s41884-026-00194-7](https://doi.org/10.1007/s41884-026-00194-7). HTML: [arxiv.org/html/2412.08400v3](https://arxiv.org/html/2412.08400v3).

> **Abstract.** It is known that the set of lumpable Markov chains over a finite state space, with respect to a fixed lumping function, generally does not form an exponential family of stochastic matrices. In this work, we explore efficiently verifiable necessary and sufficient conditions for families of lumpable transition matrices to form exponential families. To this end, we develop a broadly applicable dimension-based method for determining whether a given family of stochastic matrices forms an exponential family.

From the v3 text:

> **Corollary 4.1** (No multi-row merging block criterion). If `(𝒴, ℰ)` has no multi-row merging block with respect to `κ`, then `𝒲_κ(𝒴, ℰ)` forms an e-family.

> **Theorem 4.1** (Redundant merging block criterion). If `(𝒴, ℰ)` has a multi-row merging block with respect to `κ` which is redundant, then `𝒲_κ(𝒴, ℰ)` does not form an e-family.

> **Definition 3.1** (Merging block). Let `(x, x′) ∈ 𝒟` be such that for some `y ∈ 𝒮_x` it holds that `|{ (y, y′) ∈ ℰ : y′ ∈ 𝒮_{x′} }| > 1`. Then `(x, x′)` is a merging block of `(𝒴, ℰ)` with respect to `κ`. Furthermore, when `|𝒮_x| ≥ 2`, the merging block is multi-row.

A redundant merging block is one for which there exists `𝒯 ⊂ 𝒳` containing the block with the subgraph on `𝒯` remaining strongly connected after the block is deleted.

Author announcement (G. Wolfer, 14 March 2026): accepted at *Information Geometry*; the paper studies when a family of Markov chains compatible with a given lumping map forms an exponential family.

**What this paper does not contain.** Searched: Aut, orbit, Zariski, generic rates, coarsest lumping. None appear. The object is: *fix* a lumping map `κ`, vary the rates inside the `κ`-lumpable family, and ask whether that family is an e-family. The object Gate 1 / Q2 wants is the converse direction: *fix* a symmetry group `K` (so the rates are `K`-equivariant), vary the rates inside that linear family, and ask whether the *coarsest* strong lumping of a generic such matrix equals the `K`-orbit partition. Watanabe–Wolfer is a paper one would cite *next to* a generic-rates lemma, not *as* one.

Round 1's sentence stands: a Zariski-open generic-rates lemma, if someone wants to write it, would turn the five-width instance into a theorem with hypotheses. That remains optional theory, not a citation.

### "Barrett & Feng"

Round 1: no lumpability paper by Barrett & Feng. Best candidate, found this round:

- Marlin U. Thomas and Donald R. Barr, *An approximate test of Markov chain lumpability*, *J. Amer. Statist. Assoc.* **72** (1977), no. 357, 175–179. (Chi-squared test of a proposed lumping when `P` is unknown.)
- D. R. Barr and M. U. Thomas, *Technical Note — An eigenvector condition for Markov chain lumpability*, *Operations Research* **25** (1977), no. 6, 1028–1031. DOI [10.1287/opre.25.6.1028](https://doi.org/10.1287/opre.25.6.1028).

> Under certain conditions the state space of a discrete parameter Markov chain may be partitioned to form a smaller "lumped" chain that retains the Markov property. Existing statements of conditions on the transition matrix `P` of the original chain characterizing such lumpability are difficult to verify in practice, especially when the dimension of `P` is large. An alternate approach, based on the eigenvectors of `P`, is presented and illustrated with examples.

This is a computational test for a *declared* partition, not a theorem that coarsest = Aut-orbit. The names are Barr and Thomas, not Barrett and Feng. Nearby: Jernigan–Baran, *Statist. Probab. Lett.* (2003), chi-squared test of lumpability, correcting Thomas–Barr's degrees of freedom. Treat "Barrett & Feng" as a false lead unless a specific reference is supplied. Do not cite Barr–Thomas as that lead.

**Verdict, unchanged as a logical claim, upgraded as a dictionary:** coarsest = Aut-orbit is still not a theorem without extra hypotheses. The extra hypotheses that *do* make it a theorem are D'Angeli–Donno Thm 12 (chain/tree). The classical name of the distinction Gate 1 drew is Schurian vs non-Schurian / orbit vs equitable (Kudose Prop. 2.4, McKay graph). Watanabe–Wolfer is a 2026 neighbourhood of generic-rates, not the lemma.

---

## Q3 — Diaconis group-walk Fourier is still Schur

Round 1: vanishing of `C e^{(t−s)G} H e^{sG} B` is Schur's lemma on `Hom_K(ρ_C, ρ_H ⊗ ρ_B)`. Wigner–Eckart is the physics name, for Hamiltonians. Hänggi 1978 I–II is master-equation symmetries and Markov FDR, not this. No named Markov-generator selection rule. The gap is real.

### The group-walk case, which is not the general case

Persi Diaconis, *Group Representations in Probability and Statistics*, IMS Lecture Notes–Monograph Series Vol. 11, 1988. Chapter 3, *Random Walks on Groups*, §E *The Markov chain connection* (p. 48). Open: [projecteuclid.org, lnms/1215467408](https://projecteuclid.org/ebooks/institute-of-mathematical-statistics-lecture-notes-monograph-series/Group-representations-in-probability-and-statistics/toc/lnms/1215467399).

The Fourier transform of a probability `Q` on a finite group `G`, evaluated at an irrep `ρ`, is `Q̂(ρ) = ∑_{s} Q(s) ρ(s)`. Convolution of the walk becomes multiplication of these matrices. In a basis adapted to the irreps of `G`, the regular representation (and therefore the transition operator of a *Cayley walk*) is **block-diagonal by irreps**. Class functions (conjugation-invariant driving measures) make each block a scalar. Invariant observables live in the trivial isotypic component.

This is Schur's lemma applied to a group-walk. It is the reason an invariant function of a Cayley walk has vanishing first-order response to a perturbation that lives in a nontrivial isotypic. It is **not** a statement about a general `K`-equivariant Markov generator on a set that is not `G` itself, and it is not the Duhamel integrand `C e^{(t−s)G} H e^{sG} B` vanishing pointwise.

A later Diaconis paper in the same neighbourhood, still group-walk / lumping-by-cosets, not the generator form:

Persi Diaconis, Arun Ram and Mackenzie Simper, *Double coset Markov chains*, arXiv:2208.10699, 23 August 2022.

> Let `G` be a finite group. Let `H, K` be subgroups of `G` and `H \ G / K` the double coset space. Let `Q` be a probability on `G` which is constant on conjugacy classes. The random walk driven by `Q` on `G` projects to a Markov chain on `H \ G / K`. This allows analysis of the lumped chain using the representation theory of `G`.

Useful for Q2's "Aut-orbit is always a lumping" in the group-walk case. Not Q3.

### 2025 still optimises the response, still the opposite direction

Manuel Santos Gutiérrez, Niccolò Zagli and Giulia Carigi, *Markov matrix perturbations to optimize dynamical and entropy functionals*, arXiv:2507.14040, 21 July 2025.

> An important problem in applied dynamical systems is to compute the external forcing that provokes the largest response of a desired observable quantity. For this, we investigate the perturbation theory of Markov matrices in connection with linear response theory in statistical physics. We use perturbative expansions to derive linear algorithms to optimize physically relevant quantities …

This is the 2025 continuation of Antown–Dragičević–Froyland 2018 (round 1): *choose* `H` to maximise the response. Q3 asks for a vanishing theorem when `H` is symmetry-breaking and `C`, `B` are invariant. Opposite direction.

**Verdict, unchanged:** the gap is real. Diaconis 1988 is the citation for the group-walk special case. Stating the lemma for equivariant Markov generators, with the Duhamel integrand vanishing pointwise, remains a contribution of language, not of representation theory. The C2 verification remains the instance.

---

## Q1 leftovers — two papers that do not compete

Round 1 already closed the count. Two 2024–2026 papers sit next to NC combinatorics and do not touch the orbit formula.

Michael Dougherty and Gina Root, *Noncrossing partitions from hull configurations*, arXiv:2604.14458, 15 April 2026. 14 pages. MSC 05A18, 06A07, 52A10, 52C35.

> Each finite configuration of points in the plane determines a corresponding lattice of noncrossing partitions. When these points form the vertex set of a convex polygon, the associated lattice is the classical noncrossing partition lattice (introduced by Kreweras in 1972) … If, on the other hand, all points of the configuration lie on a common line segment, the result is a Boolean lattice. In this article, we examine the more general class of hull configurations … We prove that the corresponding lattices of noncrossing partitions are unions of maximal Boolean subposets and, under certain circumstances, have symmetric chain decompositions.

> **Theorem A** (Theorem 2.7). If `P` is a hull configuration with at least one blank side, then `NC(P)` has a symmetric chain decomposition.

> **Theorem B** (Theorems 3.6 and 3.9). If `P` is a hull configuration, then `NC(P)` is a union of maximal Boolean subposets, which are in one-to-one correspondence with the noncrossing trees on `P` with convex geodesics.

No reflection-fixed count, no Burnside orbit formula, no A007123. Simion–Ullman is the SCD of classical `NC(n)`; Dougherty–Root extend SCD to hull configurations. Different object.

Kurusch Ebrahimi-Fard, Loïc Foissy, Joachim Kock and Frédéric Patras, *Noncrossing arithmetic*, arXiv:2407.17660, v2. 20 pp. Final author version, to appear in *Category Theory in Computational Mathematics and Theoretical Informatics*.

> Higher-order notions of Kreweras complementation have appeared in the literature … The present article aims at offering a simple account of various aspects of higher-order Kreweras complementation … we exhibit the lattice of noncrossing partitions as the decalage of a partial monoid structure on noncrossing partitions encoding higher-order Kreweras complements.

Ding 2016 used the ordinary Kreweras complement as an anti-isomorphism `NC(n)^F ≅ NC(n)^{RF}` of the two even-`n` reflection classes. Higher-order Kreweras is a different (free-probability) lift of that involution. It does not re-prove, and does not compete with, Callan–Smiley Theorem 1.

Callan–Smiley, arXiv:math/0510447, remains unpublished as a journal article. Armstrong, *Generalized noncrossing partitions and combinatorics of Coxeter groups*, arXiv:math/0611106, cites it as [42]. That is still the publication status.

**Verdict, unchanged:** cite Callan–Smiley Thm 1 + Ding. Do not cite Dougherty–Root or Ebrahimi-Fard for this count.

---

## Q5 — Freund exact deflation is Grimme in block language

Round 1: no named "non-monotone prefix" pathology; moment matching holds iff the projection contains the *full* Krylov of that grade (Grimme 1997); H2 interpolatory error is not monotone (Gugercin–Antoulas–Beattie 2008); standard remedy is truncate at completed block boundaries.

Roland W. Freund, *Krylov subspaces associated with higher-order linear dynamical systems*, arXiv:math/0501484, 27 January 2005; cf. also *Model reduction methods based on Krylov subspaces*, *Acta Numerica* **12** (2003), 267–320.

The block-Krylov matrix with `m` starting vectors is scanned left to right; each column linearly dependent on columns to its left is deleted (**exact deflation**). The `n`-th block Krylov subspace `𝒦_n(M, R)` is the span of the first `n` columns of the *deflated* matrix, not of the undeflated block matrix. A prefix that cuts inside a block before a deflation boundary is therefore not a Krylov subspace of any grade.

This is Grimme's complete-grade requirement, written for multiple starting vectors. It is not a named "non-monotone prefix" pathology. Implicit restart (Sorensen / ARPACK) discards nonessential modes at restart boundaries for the same structural reason: a truncated Krylov basis that is not a Krylov subspace does not inherit the interpolatory guarantees.

**Verdict, unchanged as a claim, restated:** cite Grimme + Freund exact deflation + Gugercin–Antoulas–Beattie. Only truncate at completed (deflated) block boundaries. Do not name a new pathology.

---

## Neighbouring ticket #600 Phase C — Kalman, not a new theorem

#600 asks to build the finite-horizon balanced realization two ways (full Catalan(`w`) space vs `R`-orbit quotient) and compare Hankel singular data at `w = 4..8`. Expected result, stated in advance:

```text
all nonzero task Hankel singular data identical between A and B;
the anti-invariant (R-odd) sector exactly uncontrollable and/or unobservable.
```

This is the Kalman decomposition of a linear system that is equivariant under a linear involution, with even sources and even readouts. Literature pointer, not a P398 theorem.

### Unobservable / unreachable subspaces

W. M. Wonham, *Linear Multivariable Control: A Geometric Approach*, Springer, 1979 (2nd ed.; 1st ed. 1974). Standard definitions, as compiled e.g. in lecture notes of the KTH geometric-control course:

- The **unobservable subspace** of the pair `(C, A)` is `ker Ω`, where `Ω` is the observability matrix stacked as `C, CA, …, CA^{n−1}`. It is `A`-invariant.
- The **reachable subspace** of the pair `(A, B)` is `⟨A | im B⟩ = im [B  AB  …  A^{n−1}B]`. It is `A`-invariant.
- A subspace `V` is **`(A, B)`-invariant** (controlled invariant) iff `A V ⊂ V + im B`, iff there exists feedback `F` with `(A+BF)V ⊂ V` (Wonham–Morse 1970; Basile–Marro 1969).

The Kalman decomposition splits the state space into the four combinations of (un)reachable × (un)observable. The transfer function, and therefore every input–output invariant (Hankel singular values, balanced order, response matrices), depends only on the reachable-and-observable summand.

For a `C2`-equivariant realization `R A R = A`, `R B = B` (even sources), `C R = C` (even readouts), the odd eigenspace of `R` is unreachable from even `B` and/or unobservable from even `C`. Hence the task Hankel data of the full space and of the even (orbit-quotient) space agree on all nonzero singular values, and disagree exactly when a readout that is not even is included.

That last clause is #600's trap, already correctly stated on the ticket: at odd `w`, `halves_linked` is not `R`-even, so A and B **must** disagree on that readout (first-order integrands 0.292 and 0.107 at `w = 5, 7`). A run where B reproduces A on every readout at odd width is a bug.

A closer geometric-control neighbour, still not a P398 theorem:

Xuefeng Shen and Melvin Leok, *Geometric symmetry reduction of the unobservable subspace for Kalman filtering*, arXiv:1901.03474, 11 January 2019.

> In this article, we consider the implications of unobservable subspaces in the construction of a Kalman filter. In particular, we consider dynamical systems which are invariant with respect to a group action, and which are therefore unobservable in the group direction.

Continuous symmetry, filtering, not discrete `C2` on a Markov generator. Cite as the sentence "group-invariant directions are unobservable", not as Phase C.

**What #600 should write if the numbers come out as expected:**

> The nonzero finite-horizon Hankel singular data of the declared even sources and even readouts are identical between the Catalan space and the `R`-orbit quotient, as required by the Kalman decomposition of a `C2`-equivariant realization (Wonham). The odd-width `halves_linked` disagreement is the non-vacuous check. This is not a new theorem.

If the numbers disagree on even readouts, the deliverable is which convention broke it (source lifting, readout restriction, contrast subspace, Gramian metric), as the ticket already says. Do not patch the task to rescue the theorem.

---

## Map to open tickets

| Ticket | What round 2 adds | What it does not do |
|---|---|---|
| **#601** | Neighbourhoods of Q2–Q5; Kalman pointer for the sibling compute ticket. | Does not close #601. Owner closes against #602 if the notes are accepted. |
| **#598** | Still: orbit-*count* sentence is a citation (round 1). Partition equality stays a P398 fact. Selection-rule *verification* stays a P398 fact. | Does not rewrite the #598 note. Does not close. |
| **#593** | Unchanged: `r_positive(9)=2494`, `(10)=8524` are OEIS A007123(`w+1`). If that ticket runs, report **partition equality** at 9 and 10, not just cardinality. | Does not run widths 9/10. |
| **#596 Gate 1** | A vs B now has the classical names Schurian / orbit vs non-Schurian / equitable (Kudose Prop. 2.4, McKay graph). Five widths still do not make A a theorem. | Does not flip the gate. |
| **#596 Gate 2 / #595** | Cite Sitter–Wu 1993 F-optimality as the named Fieller-as-design paper, and distinguish it: dose-support of a binary-response experiment, not ranking of histogram channels. Shared-denominator-across-folds still unnamed; #600's block-bootstrap of the design cost is a valid replacement for the improvisation. | Does not score a channel. Does not buy N=650. |
| **#600 Phase C** | Expected result = Kalman decomposition / Wonham unobservable subspace, applied to even `B`, even `C`, involution `R`. Odd-width `halves_linked` is the teeth. | Does not run the SVDs. Compute ticket stays compute. |
| **#594 / Q5** | Freund exact deflation restates Grimme. Still no named "non-monotone prefix" pathology. | Does not close #594. |
| **#592** | Still no overlapping source. Watanabe–Wolfer is lumpability/e-families, not Mori–Zwanzig / positive realization / finite-horizon balancing. Wonham is geometric control, adjacent to balancing only through Kalman. | Does not answer #592. |
| **#370 / #579** | Particle-physics optimal observables and #579's denominator-free projective logic remain the adjacent *method* citations for Gate 2. Sitter–Wu is not a substitute. | Unchanged. |
| **#566** | Untouched. | Untouched. |

---

## X / arXiv sweep

- **X**, this round, keyword search on lumpability / NC / Fieller / Krylov / orbit / equitable: no hits in the indexed window. The one related post found via the web, not via X search: G. Wolfer (@geowolfer), 14 March 2026, announcing acceptance of Watanabe–Wolfer at *Information Geometry* ([arxiv.org/abs/2412.08400](https://arxiv.org/abs/2412.08400)).
- **arXiv 2025–2026**, lumpability / equitable / Aut-orbit / generic rates: Watanabe–Wolfer 2412.08400v3 is the live paper. No Zariski-open coarsest=`K`-orbits lemma found. Santos Gutiérrez–Zagli–Carigi 2507.14040 continues the *optimise-the-response* line. Dougherty–Root 2604.14458 and Ebrahimi-Fard 2407.17660 are NC combinatorics, not the orbit count. Diaconis–Ram–Simper 2208.10699 is double-coset lumping of group-walks.

A negative sweep is a result: the gaps named in round 1 (Markov-generator selection rule; generic-rates lemma; ranking observables by Fieller; shared-denominator-across-folds; named non-monotone-prefix pathology) were not filled by 2025–2026 literature either.

---

## What this does not do

- Not a STATUS edit, not a claim-ledger edit, not a close of #598 / #593 / #596 / #566 / #601 / #600 / #595 / #594 / #592.
- Not a recap of Callan–Smiley, Ding, D'Angeli–Donno, Hauschke, Grimme, Gugercin–Antoulas–Beattie. Those stay in the round-1 note.
- Not a theorem that coarsest D0-admissible lumping = `R`-orbit partition. Five widths remain an instance. Schurian / Kudose / Watanabe–Wolfer / Barr–Thomas do not change that.
- Not an identification of Sitter–Wu F-optimality with Gate 2's observable ranking.
- Not a named Markov-generator selection rule. Diaconis 1988 is the group-walk case.
- Not a Phase C computation, and not a new theorem if Phase C comes out as expected — that would be Kalman.
- Not an #592 packet.

Owner may close #601 against PR #602 once both notes are accepted. I am not auto-closing.
