# Venue retrieval for P4 (matching-odd scalar closure) — wrapping/crossing literature

Ticket: [#696](https://github.com/LightChainr/Matching-One/issues/696). Retrieval and notes only.
The P4 draft read for this note is `docs/manuscripts/p4-matching-odd-scalar-closure/README.md`
on branch `pr688` (= `pull/688/head`), unmerged, read with `git show pr688:...`.
No P4 result was recomputed, no Monte Carlo was run, no `p_c` claim is made, and
P3's projective statistic is not imported. `docs/STATUS.md` and `docs/ROADMAP.md` untouched.

Literature marking follows the repo convention (p2 README / p3 manuscript):
`PRIMARY_TEXT_READ` = the full text was fetched and read in this session;
`ABSTRACT_ONLY` = bibliographic record plus (part of) the abstract only;
`[LIT]` = known citation, not fetched.

## 0. What was actually fetched

| # | paper | venue | how read | marker |
|---|---|---|---|---|
| S1 | Langlands, Pouliot, Saint-Aubin, *Conformal invariance in two-dimensional percolation* | Bull. Amer. Math. Soc. **30**(1), 1–61 (1994); arXiv:math/9401222 | full text (arXiv PDF, 60 pp.) | PRIMARY_TEXT_READ |
| S2 | H. T. Pinson, *Critical percolation on the torus* | J. Stat. Phys. **75**, 1167–1177 (1994), DOI 10.1007/BF02186762 | bibliographic record + opening abstract sentence; Springer/ADS full text bot-blocked | ABSTRACT_ONLY |
| S3 | M. E. J. Newman, R. M. Ziff, *Fast Monte Carlo algorithm for site or bond percolation* | Phys. Rev. E **64**, 016706 (2001); arXiv:cond-mat/0101295 | full text (APS harvest PDF, 16 pp.) | PRIMARY_TEXT_READ |
| S4 | S. Mertens, R. M. Ziff, *Percolation in finite matching lattices* | Phys. Rev. E **94**, 062152 (2016); arXiv:1603.07289 | full text (ar5iv HTML) | PRIMARY_TEXT_READ |
| S5 | G. Pruessner, N. R. Moloney, *Winding clusters in percolation on the torus and the Möbius strip* | J. Stat. Phys. **115**, 839–853 (2004); arXiv:cond-mat/0310361 | full text (ar5iv HTML) | PRIMARY_TEXT_READ |
| S6 | J. J. H. Simmons, P. Kleban, R. M. Ziff, *Percolation crossing formulae and conformal field theory* | J. Phys. A **40**, F771–F784 (2007); arXiv:0705.1933 | full text (ar5iv HTML) | PRIMARY_TEXT_READ |
| S7 | J. Cardy, *Critical percolation in finite geometries* | J. Phys. A **25**, L201 (1992) | citation only | [LIT] |
| S8 | G. M. T. Watts, *A crossing probability for critical percolation in two dimensions* | J. Phys. A **29**, L363–L368 (1996), DOI 10.1088/0305-4470/29/14/002 (Crossref) | citation only | [LIT] |
| S9 | R. K. Akhunzhanov, A. V. Eserkepov, Yu. Yu. Tarasevich, *Exact percolation probabilities for a square lattice: site percolation on a plane, cylinder, and torus* | J. Phys. A **55** (2022); arXiv:2204.01517 | abstract + bibliographic record | ABSTRACT_ONLY |

All quotes below were taken from the fetched texts of S1, S3–S6 (or, for S2, from
the abstract fragment visible on the Springer/ADS record).

## Q1. Has a matching-odd / orientation-odd wrapping contrast been published?

**Answer: not as a cos 4θ orientation contrast of a matching function, in anything I
could read. The two halves of the observable each have a published home; their
oriented difference does not.** Caveat marked below.

The published "odd" objects I found are odd in **p**, not in orientation:

- Mertens–Ziff (S4) define the matching function on the torus and prove the exact
  wrapping identity (their Eq. 20, quoted verbatim):

  > "M_L(p) = R^x_L(p) − R̂^x_L(1−p),  x ∈ {c, b, e, h} (20) … **This is the main
  > result of this paper.**"

  This is the published ancestor of P4's matching/primal pairing `S = (R_G + R̂)/2`,
  `D = (R_G − R̂)/2`, and of the finite Russo identity the P4 draft uses as a spine
  (MZ derive `N_L(p) − N̂_L(1−p) = L²χ(p) + O(1)` from it). It is odd under
  `p → 1−p` combined with lattice duality — it says nothing about lattice
  orientation θ.
- Newman–Ziff (S3) define the four torus wrapping probabilities `R^(h), R^(v), R^(b), R^(e)`
  and the channel difference (their Eqs. 11–12):

  > "R^(1) = R^(h) − R^(b) = R^(e) − R^(h) = ½(R^(e) − R^(b)), (12)"

  a difference of wrapping channels on the **same** lattice — the closest published
  thing to P4's "first-minus-second pairing" logic, but the difference is over
  topological channel, not over orientation, and it is not matching-odd.

Orientation enters the published torus literature in two other ways, neither of
which is a spin-4 contrast of the matching function:

- **Per homology class, not harmonic.** Langlands–Pouliot–Saint-Aubin (S1, §3.7)
  ask, in their words:

  > "We can ask for the probability, always at criticality, that a given subgroup Z
  > of H1(S) is contained in the image, and expect that the response depends only on
  > the conformal class of S."

  and simulate the square torus `S1 = C/(Z + ZiZ)` at mesh 1/500 plus a branched
  double cover, reporting `π̂(1,0), π̂(0,1), π̂(1,1), π̂(1,−1), π̂(0), π̂(H)`
  (their Table 3.7). Diagonal classes (1,±1) are measured, i.e. 45°-oriented
  wrapping exists as data — but no first-vs-second orientation pairing is formed
  and no cos 4θ coefficient is extracted. Their Table 3.7 values for the two
  isomorphic tori (e.g. `π̂(1,0)` 0.1693 vs 0.1700) are presented as a
  conformal-invariance check, with differences at the level of their statistical
  error.
- **Aspect-ratio (modulus) dependence, not orientation harmonic.** Pinson (S2,
  ABSTRACT_ONLY) computes "the various crossing probabilities defined by
  R. Langlands, P. Pouliot, and Y. Saint-Aubin" on the torus with CFT; her
  winding-cluster formulas were tested numerically by Pruessner–Moloney (S5),
  who confirm the topological prediction `𝒫̂(X, r) = 𝒫̂(0, r)` ("This is in
  perfect agreement with our numerical results") and check winding numbers
  `(1,0), (1,±1), (1,±2), (1,±3)` across 14 aspect ratios 30/30 … 900/1. Again:
  an aspect-ratio (modulus) axis, no orientation-odd contrast.

A targeted search for a `cos 4θ` / four-fold harmonic of a percolation wrapping
observable returned nothing (closest hits were exact wrapping-polynomial
enumerations, S9, which are again orientation-unresolved).

**Caveat to keep:** S2 (Pinson 1994) could not be read in the primary; her paper is
the one published text where a harmonic decomposition of torus crossing data could
still be hiding. The P4 bib pass (P4 README §10 item 2) must read Pinson 1994 in the
primary before the paper claims channel novelty. Until then the correct statement is:
**new among the texts read here, with S2 unverified.**

Bonus cross-check (supports P4's channel semantics): P4's exact map
`DeltaS_cross = −DeltaS_either` (even) and `D_either = D_cross` (odd) is the
discrete counterpart of identities already in the published record —
MZ Eq. (19) `R^1_L(p) = R̂^1_L(1−p)` and NZ Eq. (11) `R^(e) = 2R^(h) − R^(b)`.
The P4 erratum (#108) reclassification is consistent with these; nothing in the
fetched texts contradicts the draft's channel algebra.

## Q2. Finite-size scalar closure: what published papers test, and what they do when it fails

Three of the read texts are genuine analogues for "signal exists; one scalar does
not close it":

1. **Mertens–Ziff 2016 (S4)** is the closest structural analogue, and it *succeeds*
   with one scalar + corrections. The scaling-limit odd matching function collapses:

   > "M_L(p) = f(z) − f(−z) (36) in the scaling limit. A scaling plot of M_L(p) is
   > shown in the inset of Fig. 5. This curve is universal for systems of this shape
   > (square torus), except for a scale factor on z."

   But the pre-asymptotic data require **two extra correction exponents** fitted
   beyond the known `1/ν = 3/4` (their Eq. 38 analysis): "These plots give
   2−x = −3.42 and 2−y = 0.705", and their summary of the residuals is worth
   quoting at the P4 draft directly:

   > "Finite-size effects explain why many of our observations and analyses agree
   > only approximately."

   This is the published version of "the amplitude-like mode closes, the shape mode
   does not": MZ need a two-exponent correction tower to make their scalar story
   consistent, and they say so in exactly the register P4 uses for its §6 failures.
2. **Pruessner–Moloney 2004 (S5)** is a *published scalar-closure failure that stayed
   a failure* — the honest precedent for P4 §6.3. Testing Pinson's CFT amplitudes
   for multiple winding clusters:

   > "one might be tempted to find a systematic dependence of C((1,0),n) and
   > α((1,0),n) on n, such as an exponential and a second order polynomial,
   > respectively. However, we were unable to identify these functions."

   and, critically for P4's methodology, their error budget statement:

   > "the main source of error is not statistical, but systematic, namely in the
   > choice of the specific function"

   A frozen-prediction paper can cite this line verbatim: at 10⁶–10⁹ samples the
   residual is the model, not the noise — which is precisely P4's situation at
   500M replicas.
3. **Newman–Ziff 2001 (S3)** shows the standard fallback when a single clean law is
   not available: use a statistic whose limit is known without a fit. Their
   threshold estimator is the maximum of the channel-difference `R^(1)`, and in 3D
   they fit only after the fact: "we can estimate p_c by varying the scaling
   exponent until an approximately straight line is produced." Published practice
   tolerates this; P4's frozen-prediction discipline is *stricter* than what NZ did
   in the same journal.

**Refinement of the ticket's question:** none of the read texts treats a scalar
failure as evidence for a non-scalar (Jordan-type) mechanism — MZ absorb it into
correction exponents, P–M leave it unidentified, NZ avoid it by statistic choice.
P4's low-rank mixing claim is therefore *stronger* than anything published in this
niche, and the paper should present the rank-2/Jordan ordering (its §7) as a
hypothesis ranking diagnostic, not as an identification — the draft already does
this (its §6.2 "compatible, not identified" clause). The MZ corrections tower
(`L^(−x)`, `L^(−y)` with x ≈ 5.42 unknown) is the conservative alternative
explanation a PRE referee will reach for first, and the draft's §7 evidence
(one projected derivative channel closing at −0.009σ while its companion misses at
+2.70σ) is exactly the observation that distinguishes the two. That contrast should
be made explicit in the manuscript; it currently lives only in P50 REPORT.md prose.

## Q3. Lattice sizes in the published wrapping/crossing FSS record

Actual linear sizes from the fetched texts:

| paper | engine | linear sizes | samples |
|---|---|---|---|
| NZ 2001 (S3) | MC (Newman–Ziff algorithm) | 128×128 (pc estimate), **1024×1024** (stretched-exponential), 512×512 (Fig. 10 context) | "more than 7.3×10⁹ separate samples, about half of which were for systems of size 128×128" |
| MZ 2016 (S4) | exact enumeration + MC | exact **L = 3–11**; MC "L = 16, 24, 32, 48", "up to L = 128" | exact, plus standard MC runs |
| LPS 1994 (S1) | MC at criticality | square torus mesh 1/**500** (500×500); branched cover square side **282**; rectangles/parallelograms to 1000×1000 (Table 3.2) | "The sample size for all our experiments was at least 10⁵, and very often" more (10⁶) |
| P–M 2004 (S5) | MC, patch-parallel | patches L = 10, 100, 1000; total lattices **300² to 30000²** (9×10⁸ sites) | "at least 10⁶ realizations"; "of the order of 10⁹ samples" for r=1 |
| SKZ 2007 (S6) | hull-walk MC | **512×512 bonds** primary; FSS at L = 64, 128, 256, 1024 | "we were able to generate 3.3×10¹¹ hulls on a lattice of 512×512 bonds" |

So the published MC record for wrapping observables sits at **L ≈ 512–1024 for the
headline finite-size statements**, with one paper (P–M) pushing to 3×10⁴, and the
exact-enumeration branch deliberately working at L ≤ 128 where exactness or frozen
predictions substitute for size.

P4's N = 145–425 sites with 500M replicas per point is therefore:

- **above** the exact-enumeration branch (MZ: L = 3–11 exact; Akhunzhanov et al. S9
  exact wrapping polynomials at small L, ABSTRACT_ONLY) where small L is
  compensated by exactness;
- **below** the MC branch's comfort zone. 500M replicas is comparable to NZ's
  7.3×10⁹ samples and P–M's 10⁹, so the *statistical* budget is normal for this
  literature; what is not normal is spending it entirely below L ≈ 512.

**Verdict on the ticket's question ("normal computational paper, or will referees
demand L ~ 10³–10⁴?"):** referees will demand *more L*, specifically at least one
lineage at L ≥ 512, because every MC FSS claim in the table above that a referee
would recognize made it at 512+. They will not demand 10⁴ — that was needed only
for P–M's winding-number tail probabilities (~10⁻⁹), not for amplitude/shape
tests. P4's own §8 item 8 ("a second independent doubling curve (new geometry, not
a rerun) is the natural referee ask") is correct, and the literature adds a second
axis: the new geometry should also be **larger**, not just different.

## Q4. Venue fit, independently

- **Physical Review E — fit.** Scope: statistical, multiphase, and soft-matter
  physics, including computational statistical mechanics. The entire
  wrapping-probability MC literature this paper would cite lives here: NZ (PRE 64,
  016706) and MZ (PRE 94, 062152). A methods-plus-percolation paper with frozen
  predictions and a negative scalar-closure result is a natural PRE
  Statistical Physics article. Analogue: **Mertens–Ziff 2016** — same journal,
  same observables, likewise without a new exponent claim at its core.
- **JSTAT / J. Stat. Phys. — fit.** Scope: statistical physics theory and
  computation (JSTAT), and the older JSP carries exactly the torus lineage: Pinson
  (JSP 75) and Pruessner–Moloney (JSP 115). Analogue: **Pruessner–Moloney 2004** —
  torus winding observables, CFT-null testing, honest systematic-error statement.
  If the paper leans on the CFT/Pinson connection, JSTAT/JSP is arguably the more
  topical home; if it leans on the frozen-prediction methodology, PRE.
- **J. Phys. A — fit (narrower).** Scope: mathematical and theoretical physics,
  including statistical mechanics and its computational study. Precedent: SKZ
  crossing formulae (J. Phys. A **40**, F771, Fast-Track) and Watts (J. Phys. A
  **29**, L363). Best if the final paper foregrounds the exact channel algebra and
  the conformal-candidate reading (`x = 21/4`, STATUS line 68 level only);
  weaker if the paper's center of mass is the frozen-prediction governance.
- **Physical Review Letters — stretch (recommend against on current evidence).**
  The ticket's tripwire applies: I fetched no PRL percolation analogue that
  carries a finite-size signal with no new exponent, no operator identification,
  and no `p_c`. The PRL-adjacent percolation results in the fetched set all had a
  sharp universal number at stake (NZ's stretched-exponential exponent 4/3 vs
  Gaussian; SKZ's exact crossing densities at 3.3×10¹¹ hulls). P4's headline is a
  rejection + three falsifications on N ≤ 425 — the PRL editors' "beginning a new
  line of important follow-up work" bar is not met by a paper whose own §8 lists a
  mandatory second lineage as missing. **Do not recommend PRL.**
- **SciPost Physics — plausible secondary, not first choice.** Scope: core
  scientific physics with open peer review; percolation FSS papers appear there.
  The open-review format would actually suit P4's heavy artifact-citation style
  (referees can check JSON fields). But it adds no reach beyond PRE, and the
  frozen-prediction governance story will read as unusual to a general
  statistical-physics readership. Agree with the prior assessment that it is not
  the home.

**Summary:** PRE first, JSTAT/JSP second, J. Phys. A if the exact-algebra angle
grows, PRL no, SciPost acceptable fallback. This independently agrees with the
prior assessment's point 1.

## Q5. Readiness: is the second independent doubling curve actually required?

What published practice shows:

- **MZ 2016** published with a single lattice family (square torus) — but with
  exact data L = 3–11 *plus* MC to L = 128, i.e. size coverage, not just replica
  coverage.
- **P–M 2004** ran the same analysis for site *and* bond percolation as an
  internal cross-check, and across 14 aspect ratios; their published confidence
  rests on that breadth.
- **LPS 1994** checked the same torus probability against an isomorphic branched
  cover — a deliberate second realization of the same conformal class (their
  Table 3.7, S1 vs S2 columns).
- **NZ 2001** validated the algorithm across site/bond, several lattices, and 2D/3D.

So the discipline in this niche is: a single-lineage claim passes only when
something else is exact or independently varied. P4 has 500M-replica precision and
frozen predictions (its equivalent of MZ's exactness), but Block III — where two of
the three scalar failures live — is one lineage (N = 65→130, 85→170 trained;
145→290 held out), and the draft says so (§8.8).

**Verdict:** a clearly labelled single-lineage limitation *can* pass at PRE (MZ's
"agree only approximately" is precedent for publishing with acknowledged
finite-size systematics), but the referee will ask for the second doubling curve
and the paper is stronger pre-empting that with one new geometry at N ≥ 512 (which
is simultaneously the Q3 answer). The draft's §10 item on this is correctly
listed as referee-ask, not as optional polish. Also worth keeping from Q2: the
paper should surface the MZ-style correction-exponent alternative explanation and
show where its data discriminate against it — that pre-empts the most likely
theory-side referee objection and costs no new computation (it is analysis of
committed artifacts).

## Stress-test of the prior Grok assessment (agree / disagree / refine)

1. **"PRE or JSTAT after a real writing pass; not PRL / SciPost main."** — **AGREE.**
   Independently derived above from fetched analogues (MZ → PRE; P–M/Pinson →
   JSP/JSTAT; SKZ/Watts → J. Phys. A as the theory-leaning alternative). The PRL
   tripwire was checked and holds: no fetched PRL analogue lacks a new
   exponent/`p_c`.
2. **"Most 'new physical fact' of the four objects *if* the signal survives larger N."**
   — **AGREE, with a refinement.** "Larger N" in this literature means L ≥ 512
   (every MC headline result in the table above was made at 512+), which is a
   different simulation regime from N ≤ 425. If the signal does *not* survive, the
   salvageable paper is the methodology-plus-falsification structure (P–M 2004 is
   precedent that an honest "we could not identify the scalar law" paper is
   publishable), not the signal claim.
3. **"N ~ 10² is the obvious referee objection; 500M replicas do not buy systematic
   finite-size control."** — **AGREE, and it can now be sharpened with numbers:**
   500M replicas is normal for this literature (NZ 7.3×10⁹ samples; P–M ~10⁹);
   the linear-size deficit (425 vs 512–1024) is the objection, and the
   correction-exponent alternative (MZ's `2−x = −3.42`, `2−y = 0.705` tower) is
   the form the objection will take. Referees will not require 10⁴ sites — that
   demand belongs to tail-probability studies (P–M's ~10⁻⁹ winding probabilities).
4. **"Mechanism claim must stop at low-rank non-scalar mixing; no LCFT module."**
   — **AGREE.** Nothing in the fetched texts supports an operator identification;
   the closest published treatments (MZ's correction exponents; P–M's refusal to
   identify `C(n)`, `α(n)`) both stop before any operator claim. P4's draft §7
   already stops exactly there; keep it.
5. **"Do not use the N=580 aspect ladder as physics evidence."** — **AGREE.** The
   draft's own §8.2 records the ladder verdict (`underpowered`) and keeps it a
   limitation; nothing in this retrieval touches that.

Net: five agreements, two refinements (the L ≥ 512 number; the correction-exponent
alternative explanation that should be named and tested in the manuscript's
analysis, not just anticipated).

## Claim boundary restated

- No `p_c` inference anywhere above. No STATUS/ROADMAP promotion.
- All numbers quoted from fetched literature are identified by paper and, where
  relevant, by equation/table; nothing was recomputed from P4 artifacts.
- P3's projective statistic appears nowhere above as a P4 result.
- S2 (Pinson 1994) is the one unresolved primary; it is flagged rather than guessed
  at. If the P4 bib pass cannot obtain the text, the channel-novelty claim must be
  softened to "new among surveyed primary sources."
