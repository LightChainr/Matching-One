# Literature officer — #601 round 3: the reflection count was published in 2004, and other corrections

**Date:** 2026-09-08
**Role:** theory input. Same epistemic level as [#602](https://github.com/LightChainr/Matching-One/pull/602). **Does not enter** `docs/STATUS.md`. **Does not score** a frozen block. **Does not close** #598, #593, #596, #600, #601.
**Ticket:** [#601](https://github.com/LightChainr/Matching-One/issues/601), third retrieval round. Distinct from [#592](https://github.com/LightChainr/Matching-One/issues/592). No source below answers both.
**Access:** web retrieval in this session was **not** arXiv-blocked; primary sources (arXiv abs pages, Crossref, Semantic Scholar API, OEIS internal, publisher pages) were fetched directly. Quotes below are from fetched text.
**Method note:** one round-3 lead arrived without verifiable primary text and was **dropped** after checking; flagged inline. Everything quoted below was verified against a fetched primary source or a publisher/Crossref record.

This is **not a recap** of rounds 1–2 (`notes/literature-officer-20260906-issue601.md`, `notes/literature-officer-20260906-issue601-round2.md`). Round 1/2 claims are restated only where round 3 **corrects or strengthens** them.

---

## Headline: Q1's count is published, peer-reviewed, 2004 — and prior rounds missed it

Round 1's verdict for Q1 ("cite Callan–Smiley, arXiv:math/0510447, unpublished as a journal article") is **incomplete**. The primary published source is:

Shu-Chiuan Chang, Jesper Lykke Jacobsen, Jesús Salas and Robert Shrock, *Exact Potts Model Partition Functions for Strips of the Triangular Lattice*, **J. Statist. Phys. 114 (2004) 763–823**. DOI [10.1023/B:JOSS.0000012508.58718.83](https://doi.org/10.1023/B:JOSS.0000012508.58718.83). arXiv:cond-mat/0211623 (v2 "Version published in journal").

Callan–Smiley themselves say so. From the fetched HTML of arXiv:math/0510447v3, **Note added:**

> Note added: A result equivalent to Theorem 3 and its consequences appeared in [6]. In particular, the formula we give for NC Dihedral Classes is **Corollary 2.1 in that paper**, which is, we believe, its first occurrence in the literature. We thank the authors of [6] for notifying us of these facts.

with reference [6] being exactly the Chang–Jacobsen–Salas–Shrock (hereafter CJS) paper above. CJS derive the count through Potts-model transfer-matrix connectivity sums on an `L`-cycle with reflection, in the process **splitting the even-`L` reflections into the two types** (axis through vertices vs axis through edge midpoints) and deriving their generating functions separately — which is precisely the "vertex-axis vs edge-axis" question #601 Q1 asked about, answered in a published paper three years before Callan–Smiley and twelve years before Ding's thesis.

**Consequences for the citation chain of Q1:**

| Claim | Round 1 cite | Round 3 correction |
|---|---|---|
| Reflection-fixed NC count = `C(w, ⌊w/2⌋)` | Callan–Smiley Thm 1 (unpublished) | **CJS 2004, Corollary 2.1** (published); Callan–Smiley as the NC-language restatement |
| Both even-`w` reflection classes give the same count | Ding 2016 thesis (Kreweras anti-isomorphism) | CJS 2004 derives **both classes' generating functions separately** in the published text; Ding's lattice anti-isomorphism remains the structural explanation |
| Orbit formula `[Cat_w + C(w,⌊w/2⌋)]/2` | "Cauchy–Frobenius" in Callan–Smiley | unchanged; also **OEIS A007123 formula line**, verified this round in the OEIS internal page: `a(n+1) = (Catalan(n) + binomial(n, floor(n/2)))/2` |

The OEIS verification also confirms the sequence's standing interpretation as "number of noncrossing set partitions of [n] up to reflection (i<->n+1-i)" (Callan's 8 Oct 2005 comment is still on the entry), so `r_positive(w) = A007123(w+1)` stands, now with a published primary source behind it.

What stays ours is unchanged: the **partition equality** (coarsest D0-admissible lumping = `R`-orbit partition, block for block, widths 4–8) is P398 computation, not combinatorics, and cardinality still identifies no particular `R`.

---

## Q2 — two nearest-neighbour theorems, still no conjunction

Round 1/2 verdict stands: Aut-orbit partitions are always lumpings; "coarsest strong lumping = Aut-orbit partition" is not a theorem without hypotheses (D'Angeli–Donno Prop. 13 / Thm 12); Gate 1 branch A at five widths is an instance. Round 3 adds the two closest works in either direction, neither of which was in rounds 1–2:

1. **Orbits-are-lumpable direction (Markov-chain, graph-automorphism form).** P. L. Simon, M. Taylor, I. Z. Kiss, *Exact epidemic models on graphs using graph-automorphism driven lumping*, **J. Math. Biol. 62 (2011) 479–508**. DOI [10.1007/s00285-010-0344-x](https://doi.org/10.1007/s00285-010-0344-x). Their theorem (orbit classes of the graph automorphism group acting on the `2^N` CTMC state space "yield a lumping of the system") is the epidemic-modelling home of "Aut-orbit ⇒ lumpable" — one direction only. They do **not** prove the orbit partition is the *coarsest* lumpable partition (coarseness is only illustrated on the complete graph, where the orbit partition happens to be the population-count partition). Cite this as the Markov-chain prior art for the direction that *is* a theorem, next to Godsil–Royle Ch. 9.

2. **Generic-rates direction (rigidity).** J. Eilertsen, V. G. Romanovski, S. Schnell, S. Walcher, *Lumping of reaction networks: Generic and critical parameters*, **arXiv:2606.28895** (v1 27 Jun 2026, v2 2 Jul 2026; math.DS). For parameter-dependent mass-action networks: at **generic** parameters (an open set in parameter space), exact linear lumping yields only "obvious" reductions; nontrivial lumpings exist only on **critical algebraic subvarieties**, identified algorithmically via finitely many polynomial equations. This is the closest published object to the Zariski-open generic-rates lemma round 1 asked for — same logical shape (generic ⇒ no lumping beyond the forced ones) — but phrased for chemical reaction networks and linear lumpings of ODEs, **not** for Markov generators and **not** identifying the coarsest lumping with symmetry orbits. It is a neighbourhood, not the lemma. Consistent with round 2's verdict on Watanabe–Wolfer: cite *next to* a generic-rates lemma, not *as* one.

Clean negatives this round, checked directly: Rubino–Sericola weak-lumpability papers (no orbit/symmetry content); Cardelli–Tribastone–Tschaikowski–Vandin CRN bisimulation (coarsest bisimulation by partition refinement, never identified with Aut-orbits); Stochastic Well-formed Nets (Chiola–Dutheillet–Franceschinis–Haddad 1993 and successors — symmetry exploited for aggregation, no coarsest=orbits theorem); Derisavi–Hermanns–Sanders 2003 (algorithmic coarsest lumping, no symmetry content). An arXiv full-text query for `"orbit partition" AND "Markov chain"` returned one unrelated hit (Choi–Lim–Wang, arXiv:2512.13067, group-averaged MCMC kernels — projection chains of group-averaged kernels, not lumpability of equivariant generators).

---

## Q3 — the gap survives, and the strongest review has both halves but never joins them

Round 1/2 verdict stands: abstract vanishing is Schur; no named selection rule for first-order response of an equivariant **Markov generator** was found. Round 3 adds:

1. **The review that contains both halves, and never connects them.** P. Hänggi and H. Thomas, *Stochastic processes: Time evolution, symmetries and linear response*, **Phys. Rep. 88(4) (1982) 207–319**. Verified bibliographic record (ScienceDirect / ADS). This is the flagship review covering *both* master-equation symmetries (a symmetry group of the process, ergodic classes, orbit structure — via group actions on state space, **no irreducible representations, no Schur/Wigner–Eckart argument, no selection rules**) *and* Markov linear response (Dyson equation, susceptibility, spectral representations — developed for a **completely arbitrary** perturbation, no symmetry input). That the single most complete treatment of "symmetries of stochastic processes" and "linear response of stochastic processes" never states the vanishing rule is direct, checkable evidence that the gap is real in the classical literature. **This supersedes round 1's pair of Hänggi 1978 I/II citations as the canonical "review has both halves" reference** (1978 I/II can remain for the FDR point specifically).

2. **Nearest quantum neighbours (classification, not response rules).** The Lindblad-literature neighbours of the statement are classification works: symmetry classifications of Liouvillian dynamics and tenfold-way-type schemes for open systems (e.g. Kawabata–Kulkarni–Li–Numasawa–Ryu, *Symmetry of open quantum systems: classification of dissipative quantum chaos*, PRX Quantum 4, 030328 (2023); arXiv:2212.00605; and the many-body Lindbladian tenfold-way line). These classify Liouvillian symmetries and spectral statistics; none states a response-vanishing selection rule by irrep content. Similarly the Lindblad linear-response/Kubo literature (e.g. Albert–Bradlyn–Fraes–Jiang-type work on geometry and response of Lindbladians, PRX 6, 041031 (2016); arXiv:1512.08079) develops response formulas without a symmetry-based vanishing theorem. One quantum-adjacent name from the search round ("symmetries and steady states of open quantum systems", Buyskikh et al.) **could not be verified to exist** on arXiv and is **dropped** from this packet; the two verified classification/PRX cites above carry the point.

3. **2025–2026 FRR line still opposite-direction.** The nonequilibrium fluctuation–response relation literature for Markov jump processes (Marconi–Puglisi–Rondoni–Vulpiani 2008 lineage; 2025–2026 continuations) gives exact response–correlation relations for **generic** perturbations — the existence direction, never a symmetry-forced vanishing. Stochastic-resonance literature treats symmetry breaking as *producing* response. No selection rule anywhere.

**Verdict, unchanged:** the gap is real. The lemma-for-generators (Duhamel integrand vanishing pointwise, `C₂`-parity case verified at 1e-15 in #598) remains a small genuine contribution of language. Cite Schur/Wigner–Eckart for the vanishing, Hänggi–Thomas 1982 as the review that demonstrably has both halves without joining them, and Diaconis 1988 Ch. 3E for the group-walk special case (round 2).

---

## Q4 — simultaneous Fieller intervals with a shared denominator exist as a named literature; the two design gaps stay open

Round 1/2 verdict stands on both gaps. Round 3 adds a literature that prior rounds called "adjacent, not found":

1. **Shared-denominator Fieller is a developed method — as simultaneous inference.** G. Dilba, F. Bretz, V. Guiard, *Simultaneous confidence sets and confidence intervals for multiple ratios*, **J. Statist. Plann. Inference 136(8) (2006) 2640–2658**. DOI [10.1016/j.jspi.2004.11.009](https://doi.org/10.1016/j.jspi.2004.11.009) (verified via Crossref). Generalizes Fieller from one ratio to **simultaneous Fieller-type intervals for multiple ratios sharing a control denominator** (Dunnett-type many-to-one contrasts; multivariate-t joint distribution whose correlations depend on the shared control mean and residual variance). Antecedent: C. W. Dunnett, *A multiple comparison procedure for comparing several treatments with a control*, **JASA 50 (1955) 1096–1121**. DOI [10.1080/01621459.1955.10501294](https://doi.org/10.1080/01621459.1955.10501294) (verified via Crossref). — joint fiducial-type limits for potencies of several drugs **relative to a common standard**, i.e. the shared-denominator structure. Software embodiment: R package **mratios** (`sci.ratio`, `n.ratio` for sample size on relative margins; cites Dilba–Bretz–Hothorn–Guiard, Statist. Med. 25 (2006) 1131–1147).

   **What this does and does not do for the improvised remedy.** Gate 2's cross-validation instability is: the *design-cost statistic* is itself a ratio with a weakly estimated denominator, so per-fold costs are unstable; the improvisation shares one denominator estimate across folds. Dilba–Bretz–Guiard / Dunnett share the denominator across **comparisons** to get simultaneous coverage — not across **CV folds** to stabilise an estimated design criterion. So the shared-denominator structure is named and classical, but the design/CV use remains an improvisation, as round 2 concluded. Cite Dilba–Bretz–Guiard + Dunnett as "sharing the denominator is classical for simultaneous ratio inference", and keep the fold-stabilisation as ours.

2. **Fieller-as-design lineage has one more entry.** I. Fainaru, *F-optimal designs for binary response experiments*, M.Sc. thesis, Carleton University, 1994 (DOI 10.22215/etd/1994-02810) — extends Sitter–Wu 1993 within the same dose-support setting. No observable-selection variant.

3. **Clean negatives, rechecked:** no metrology/assay/sensor literature ranks **observables** by Fieller-length or interval-width design criteria; no named "cross-validation of a design criterion" theory for ratios; the weak-denominator → unbounded Fieller-set pathology is documented in the Mendelian-randomization weak-instrument literature (e.g. Nat. Commun. 2022, DOI 10.1038/s41467-022-28553-9) — the closest statement of the *failure mode*, still with no shared-denominator-across-folds remedy. Finney's *Statistical Method in Biological Assay* remains unverified at primary level this round; Dunnett 1955 is the verifiable anchor.

---

## Q5 — clean negative confirmed with one practitioners' artifact

Round 1/2 verdict stands: no named "non-monotone prefix" pathology; the phenomenon is incomplete-grade moment matching (Grimme) plus non-monotone interpolatory/H₂ error (Gugercin–Antoulas–Beattie); remedy is truncate at completed (deflated) block boundaries (Freund). Round 3 adds:

1. **Practitioners' folklore, unnamed.** Siyang Hu, Nick Wulbusch, Alexey Chernov, Tamara Bechtold, *Error Estimation and Stopping Criteria for Krylov-Based Model Order Reduction in Acoustics*, **arXiv:2412.10559** (13 Dec 2024; math.NA). They analyse the industrial stopping estimator ‖G_{r+1} − G_r‖/‖G_r‖ (consecutive-Krylov-order ROM difference; heuristic from Bechtold–Rudnyi–Korvink, J. Micromech. Microeng. 15 (2004) 430, implemented in Ansys MOR) and report that the true error **converges to a minimum (a floor), not monotonically to zero**, and recommend — quoting — "due the oscillation of the estimator, one may consider some kind of low pass filtering, e.g., moving average." So error oscillation along the Krylov-order ladder is acknowledged, smoothed away, never named or explained, and no earlier non-monotonicity citation is given. This is the strongest "practitioners know this" artifact found; cite it as such.

2. **The structural contrast is implicit, not named.** Balanced truncation carries the monotone a priori bound ‖G − G_r‖_∞ ≤ 2·Σ_{j>r} σ_j (Enns 1980; Glover, Automatica 20 (1984) 285–300; Antoulas 2005 Ch. 7) — monotone in `r` by construction of the sorted Hankel tail. Krylov/interpolatory projection carries **no** a priori monotone bound (Gugercin–Antoulas–Beattie 2008: IRKA is not descent). No source states "BT monotone, Krylov not" as a theorem; it follows by combining the bound with GAB. The honest framing of the width-5–8 ladder (rank 7 worse than 6, etc.) is: **an instance of the known absence of any monotone guarantee off complete Krylov grades**, not a new pathology.

3. **One more "more order can hurt" precedent.** Z. Bai and R. W. Freund, *A partial Padé-via-Lanczos method for reduced-order modeling*, **Linear Algebra Appl. 332–334 (2001) 139–164**. DOI [10.1016/S0024-3795(00)00291-3](https://doi.org/10.1016/S0024-3795(00)00291-3) (verified via Crossref). Full-PVL models of growing order can **acquire unstable poles**; the remedy is a controlled partial truncation of the Krylov chain. About stability, not held-out error, but the earliest published "higher Krylov order can degrade the model; stop at a controlled point" pattern found.

---

## Corrections ledger (what round 3 changes)

1. **Q1 cite list re-anchored to published literature.** CJS 2004 (J. Statist. Phys. 114, 763–823) is the primary published source for the reflection-fixed count including the two even-`w` reflection classes; Callan–Smiley is the NC-language restatement and credits CJS Corollary 2.1 in its own Note added; Ding 2016 remains the lattice-structural explanation. Any write-up citing "unpublished Callan–Smiley" should be corrected.
2. **Rao–Suk has a journal version** missed by round 1: *Dihedral sieving phenomena*, **Discrete Math. 343(6) (2020) 111849**, DOI 10.1016/j.disc.2020.111849 (verified via Crossref; its reference list includes Ding 2016). Round 1 cited only arXiv:1710.06517v3.
3. **Q3's "review with both halves" citation upgraded** from Hänggi 1978 I/II to Hänggi–Thomas 1982 Phys. Rep. 88(4), with the explicit observation (verified against the review's structure) that its symmetry section and its linear-response section never connect — the checkable evidence for the gap.
4. **Q4's shared-denominator story upgraded** from "not found, adjacent to pooled control" to "classical as simultaneous Fieller inference (Dunnett 1955; Dilba–Bretz–Guiard 2006; mratios), still not a CV/design-stability remedy."
5. **Q2 gains two nearest neighbours** (Simon–Taylor–Kiss 2011 Thm: orbits lumpable, one direction; Eilertsen et al. 2026 arXiv:2606.28895: generic-parameters rigidity in CRN lumping) — neither fills the coarsest=orbits gap; the gap verdict is unchanged.
6. **Q5 gains one practitioners' artifact** (Hu et al. arXiv:2412.10559 — oscillation acknowledged, smoothed, unnamed) and one precedent (Bai–Freund 2001 partial PVL). Verdict unchanged.

## What this does not do

- Not a STATUS edit, not a claim-ledger edit, not a close of #598 / #593 / #596 / #600 / #601.
- Not a recap of round 1/2 sources, which stand unless named above.
- Not a theorem that coarsest D0-admissible lumping = `R`-orbit partition; five widths remain an instance.
- Not a Markov-generator selection rule; the Q3 gap stands.
- Not an identification of Sitter–Wu / Fainaru F-optimality or Dilba–Bretz–Guiard with Gate 2's observable ranking.
- Not a named non-monotone-prefix pathology; Q5's verdict stands.

Owner may close #601 against the round-1/2 PR once the notes are accepted; this packet amends, not reopens. I am not auto-closing.
