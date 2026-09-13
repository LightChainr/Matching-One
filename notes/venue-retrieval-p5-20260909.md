# Venue retrieval — P5 / AP-2026-0311 after the JAP/AAP desk reject (issue #697)

Date: 2026-09-09. Ticket: #697 (parent #650). Manuscript object: PR #688 draft
`docs/manuscripts/p5-exact-continuation-bounded-state/README.md` (read from local branch
`pr688`, **not merged**, not checked out as this branch) and the submitted proof PDF
AP-2026-0311, *"Homological Occupation Growth: Branching Tests and Cut-Network
Representations"* (Ying/Yu, submitted 2026-09-03, 17 pp., desk-rejected 2026-09-09 by
JAP/AAP with *"not in a suitable scope … insufficient methodological novelty from an
applied probability view"*).

**Retrieval and notes only.** No STATUS/ROADMAP edits, no `references.bib` edits, no
gadget/rank computation, no Huawei. The AP proof PDF was read locally via `pdftotext`
and is **not** committed here.

## 0. What the two objects are

The submitted 17-page paper and the #688 draft are **one paper object**, but they are not
identical in coverage (this matters for the split question, §3):

| Content | AP-2026-0311 (submitted) | P5 #688 draft |
|---|---|---|
| Cut-network representation + converse realization (torus, rank-one) | Thms 4.1, 4.4, 5.1, 5.2, 6.1 (submitted core) | Part I, Thms A–D |
| Depth-as-resource / no-bounded-radius / positive-vs-signed rank | not in submission | Thms B, C, D |
| Kalman finite-horizon no-go (`p_c = 1/2` vs `1/3`, identical responses) | not in submission | Part II, Thm E |
| P398 reflection quotient (w=4..8 instance) | not in submission | Part III, Thms F–H |
| 4×4 torus witness `(1,7/8,9/14,5/14,4/35,…)` + fork gap 1/98 | §§3, 6 (Thms 3.1, 3.3, 6.1) | Part I opening fact |

Everything below treats them as one object; §3's split recommendation never implies two
overlapping submissions.

## 1. Sources fetched, with marking

| Source | URL / location | Marking |
|---|---|---|
| JAP about page (Cambridge, Applied Probability Trust) | cambridge.org/core/journals/journal-of-applied-probability | **PRIMARY_TEXT_READ** |
| JAP Vol 63 Issue 3 (Sept 2026) table of contents | cambridge.org/core/journals/journal-of-applied-probability/issue/E52E939B7D9921109FE29020814EB70D | **ABSTRACT_ONLY** (titles + one-line summaries; full texts not read) |
| AAP about page (Cambridge) | cambridge.org/core/journals/advances-in-applied-probability | **ABSTRACT_ONLY** (page fetched; contains no scope text — AAP scope is stated on the JAP page as the joint-submission companion) |
| AAP "Phase Transitions" collection preface (June 2025) | cambridge.org/core/journals/advances-in-applied-probability/collections/june-2025-phase-transitions | **PRIMARY_TEXT_READ** (preface only) |
| AAP "Branching Processes" collection (Aug 2025) article list | …/collections/august-2025-branching-processes | **ABSTRACT_ONLY** |
| CPC about page (Cambridge) | cambridge.org/core/journals/combinatorics-probability-and-computing | **PRIMARY_TEXT_READ** (scope quoted verbatim below) |
| EJP journal page (IMS) | imstat.org/journals-and-publications/electronic-journal-of-probability/ | **PRIMARY_TEXT_READ** (scope statement quoted below) |
| ALEA homepage + Vol XXIII TOC | alea.impa.br | **ABSTRACT_ONLY** (no scope statement on page) |
| EJC homepage + Vol 33 Issue 3 TOC | combinatorics.org | **ABSTRACT_ONLY** (issue TOC read; scope statement not on the page) |
| *Linear Algebra and its Applications* aims & scope (Elsevier) | sciencedirect.com/journal/linear-algebra-and-its-applications/about/aims-and-scope | **PRIMARY_TEXT_READ** |
| Random Structures & Algorithms (Wiley) | onlinelibrary.wiley.com/journal/10982418 | **NOT_FETCHED** (two attempts returned navigation chrome only; no scope text). Tripwire is satisfied by CPC + EJP. |
| J. Phys. A scope page (IOP) | iopscience.iop.org/journal/1751-8121/page/Scope | **NOT_FETCHED** (404) |
| Callan–Smiley, *Noncrossing partitions under rotation and reflection*, arXiv:math/0510447 (2005) | arxiv.org/abs/math/0510447 | **ABSTRACT_ONLY** (abstract read; confirms reflection-invariant NC partition count via the Kreweras involution; Thm 1 text not read in this pass) |
| D'Angeli–Donno, *The lumpability property for a family of Markov chains on poset block structures*, Adv. Appl. Math. 51(3):367–391 (2013), arXiv:1304.4180 | arxiv.org/abs/1304.4180 | **ABSTRACT_ONLY** |
| AP-2026-0311 proof PDF (17 pp.) | local, `/Users/lc/Downloads/AP-2026-0311_Proof_hi (1).pdf`, read-only | **PRIMARY_TEXT_READ** (full text extracted; not committed) |
| P5 README on `pr688` | `git show pr688:docs/manuscripts/p5-exact-continuation-bounded-state/README.md` | **PRIMARY_TEXT_READ** |
| Kemeny–Snell *Finite Markov Chains*; Wonham (Kalman decomposition); Schur/Wigner–Eckart; Ding 2016 (Miami PhD); Yannakakis / Fiorini et al. (slack-matrix nonnegative rank); Balle–Panangaden–Precup LICS 2015; Larsen–Skou 1991; Deng–Feng 2017 | — | **[LIT]** — cited from the AP PDF's own reference list and the #688 draft; primaries not fetched in this pass |

Tripwire check: JAP scope fetched; CPC and EJP scope fetched. All three quoted below.

## 2. Q1 — Was the desk reject predictable from the journals' own scope?

**Fetched scope text (JAP, verbatim, Cambridge Core):**

> "With a publication record spanning more than five decades, the *Journal of Applied
> Probability* is the oldest journal devoted to the publication of research in the field
> of applied probability. It is an international journal published by the Applied
> Probability Trust, and it serves as a companion publication to the *Advances in
> Applied Probability*. Its wide audience includes leading researchers across the entire
> spectrum of applied probability, including biosciences applications, operations
> research, telecommunications, computer science, engineering, epidemiology, financial
> mathematics, the physical and social sciences, and any field where stochastic modeling
> is used.
> A submission to Applied Probability represents a submission that may, at the
> Editor-in-Chief's discretion, appear in either the *Journal of Applied Probability* or
> the *Advances in Applied Probability*. Typically, shorter papers appear in the
> *Journal*, with longer contributions appearing in the *Advances*."

Note the scope statement is *application-domain* broad, not *method* broad. Nothing in it
excludes reliability, graphs on tori, or finite Markov chains — and indeed the AP PDF
itself cites two JAP signature papers (Boland, JAP **38** (2001) 597–603;
Navarro–Samaniego–Balakrishnan, JAP **47** (2010) 235–253) — so "reliability is foreign
here" would be the wrong reading of the reject.

**The journal's actual diet.** From the fetched JAP 63(3) (Sept 2026) TOC — 20 articles —
the percolation / branching / Markovian-network items are (all **ABSTRACT_ONLY**):

1. M. Marivain, *"Brochette first-passage percolation"* — a new FPC model; time constant
   and shape theorem.
2. D. Bertacchi, F. Zucca, *"Strong survival and extinction for branching random walks
   via new order for generating functions"* — classification theorem for multitype
   branching random walks.
3. B. Lodewijks, *"A star is born: Explosive Crump–Mode–Jagers branching processes"* —
   explosion criterion, applied to preferential attachment.
4. D. Goldsztajn, K. Avrachenkov, *"Asymptotically optimal policies for weakly coupled
   Markov decision processes"* — fluid limits, asymptotic optimality.
5. Ü. Işlak, B. Yeşiloğlu, *"Analysis of clustering and degree index in random graphs
   and complex networks"* — limit theorems for network indices in ER/BA/WS models.

The AAP side shows the same shape: the fetched June-2025 *Phase Transitions* collection
preface frames the journal's interest as "bulk properties … undergo an abrupt change as a
parameter of the system moves across a critical value" in "probabilistic models of
statistical physics … percolation and spin systems … stochastic models of networks,
branching and replication, epidemics, algorithms".

**Comparison verdict.** Every sampled analogue is an *asymptotic* paper: a shape theorem,
a survival/extinction classification, an explosion criterion, a limit theorem, a CLT.
AP-2026-0311 contains no asymptotic statement of any kind — its theorems are exact finite
representation results (cut along an essential cycle ↔ two-terminal connectivity; every
cofacial plane network realizes on a torus; k+1 predictive classes on 8k vertices) plus
finite witnesses. It also *concedes in its own §8.1* that the testing instrument is
prior: "The abstract ability of copying tests to refine linear trace information is
therefore prior theory." So the desk-reject letter's two clauses — scope, and
"insufficient methodological novelty from an applied probability view" — are both
consistent with what the journals' own pages and TOCs show. **The reject was predictable
as an audience call, and it is an audience call, not a finding that any construction is
false.**

Per the tripwire: among the fetched items above, **none is a construction/gadget-based
paper like AP-2026-0311** — the closest (Marivain) introduces a new model but immediately
proves asymptotics. Therefore this note does **not** recommend resubmitting the same text
to JAP/AAP, in any form. (Allowed answer: **no**.)

## 3. Q3 — Split vs keep-one (taken before Q4 so the venue rows have a target)

Independent of Grok, and keeping the one-paper-object rule:

- **Part I content (= essentially the submitted AP paper + Thms B/C/D).** This is a
  self-contained representation theorem with exact witnesses. It is the submission-ready
  unit. It should go to **one** venue as one paper (see §4; best ticket-list fit: EJC).
- **Part II (Kalman no-go).** Theorem E is genuinely disjoint from Part I: different
  object (declared `(G,B,C)` tasks), different method (realization theory), different
  literature (control, identifiability). As a standalone counterexample paper it does
  not overlap the Part I submission — no shared text, no shared theorems. A separate
  submission of Part II alone is therefore **not** a violation of the one-paper rule,
  but it is only worth doing if the honesty ladder (§4c of the #688 draft) survives
  review as a framing, since Theorem E alone is a short paper.
- **Part III (P398 quotient).** Not ready: the quotient identification is proved at five
  widths (`w = 4..8`, blocks `10, 26, 76, 232, 750`) and is stated as an **instance**,
  not a theorem; the generic-rates lemma is unwritten (#688 draft §8.1). No venue
  recommendation for Part III in its current state. It should **not** be split off now,
  and should not ride along as a third of a submission either — if Part I goes out
  alone, Part III stays in the repo until the lemma exists.

Net: **keep-one for what was submitted (Part I, possibly enriched by Thms B/C/D); one
optional separate short paper for Part II; nothing for Part III yet.** This never puts
the same theorems in front of two editors.

## 4. Q4 — Venue fit, with fetched scopes

- ***Combinatorics, Probability & Computing.*** Scope (fetched, verbatim): "devoted to
  the three areas of combinatorics, probability theory and theoretical computer science.
  Topics covered include classical and algebraic graph theory, extremal set theory,
  matroid theory, probabilistic methods and random combinatorial structures;
  combinatorial probability and limit theorems for random combinatorial structures; the
  theory of algorithms …". Plausible on paper, but the center of gravity is random
  structures and asymptotics. No analogue of an exact finite representation theorem with
  witnesses was fetched from CPC. **Second-best ticket-list fit for Part I** (the
  graph-theory clause is real), but it would be a framing stretch.
- ***Random Structures & Algorithms.*** **NOT_FETCHED** (see §1) — no recommendation
  beyond Grok's own, which this retrieval can neither confirm nor contradict. Flagged
  for a future pass.
- ***Electronic Journal of Probability / ALEA.*** Scope (fetched, EJP, verbatim):
  "publishes full-length research articles in probability theory. Short papers, those
  less than 12 pages, should be submitted first to its sister journal, the Electronic
  Communications in Probability (ECP)." The board (fetched) is strongly
  interacting-particles/percolation weighted (Toninelli ed.; Arguin, Manolescu, Camia,
  Teixeira, …). The fetched ALEA Vol XXIII TOC is likewise limit-theorem-heavy. Fit:
  **Part II**, reframed as an exact identifiability/counterexample result for finite
  Markov chains — and only via ECP if kept under 12 pages. No fetched EJP/ALEA analogue
  of either part; EJP for Part I would be a poor fit (no asymptotics).
- ***Electronic Journal of Combinatorics.*** No scope statement on the fetched page
  (ABSTRACT_ONLY), but the current issue (33(3), 2026, fetched TOC) contains exact
  finite combinatorics of exactly AP-2026-0311's flavor: vertex-connectivity and
  reliability-adjacent graph papers, enumeration with small witnesses, Cayley-graph
  random-walk papers. **Best ticket-list fit for Part I** — the cut-network theorems and
  the seven-vertex witness are finite combinatorics of networks, and the paper's
  reliability citations ([6] Brown et al., *Networks* 76 (2020)) show where the
  neighboring literature lives. Risk to state honestly: EJC editors will ask for the
  combinatorial payoff to be front-and-center, i.e. the "homological" wrapper becomes
  context, not headline.
- ***Linear Algebra and its Applications*** (and *Systems & Control Letters*, not
  fetched). LAA scope (fetched, verbatim): "publishes articles that contribute new
  information or new insights to matrix theory and finite dimensional linear algebra …
  and … significant applications of matrix theory or linear algebra to other branches of
  mathematics and to other sciences provided they contain ideas and/or statements that
  are interesting from linear algebra point of view." Theorem E is literally a statement
  about controllable/observable ranks of `(G,B,C)` triples — **LAA is a credible Part II
  home**, arguably better than EJP because the audience already owns Kalman machinery.
  SCL (**NOT_FETCHED**) would be the control-native alternative.
- ***J. Phys. A / JSTAT.*** **NOT_FETCHED** (404 on the IOP scope page). No
  recommendation recorded. Note the draft's own boundary — P398 is not percolation —
  makes these venues a framing trap for Part III anyway.
- ***JAP/AAP as a resubmit.*** **No.** Per §2: no fetched analogue as
  construction-based as this paper, so the tripwire forbids recommending it.

## 5. Q2 — Classical ingredients: used or reproved?

| Ingredient (primary) | Status in P5/AP | Verdict |
|---|---|---|
| Kalman controllable/observable rank, finite-horizon realization (Wonham; standard control texts) **[LIT]** | Part II *uses* it: exact task order **defined** as the Kalman CO dimension; Phase H factorization is "the Kalman decomposition applied to a C₂-equivariant realization — cite Wonham" (#688 draft §5) | **Uses** |
| Exact lumpability of Markov chains (Kemeny–Snell **[LIT]**) | AP PDF Thm 3.3 *applies* the Kemeny–Snell criterion ("Direct exact backward refinement … splits four survival classes … Hence the survival partition fails the Kemeny–Snell criterion"); P5 Part III *uses* coarsest-lumping refinement as a tool | **Uses** |
| Schur / Wigner–Eckart vanishing of odd matrix elements **[LIT]** | P5 Part III Thm G: the vanishing is the classical statement; the paper's claim is **language + pointwise verification** on a declared grid, "claimed as a remark citing Schur, not a theorem" (#688 draft §5) | **Uses, with a new pointwise-formulation remark** |
| Callan–Smiley reflection-fixed noncrossing counts (arXiv:math/0510447 — **ABSTRACT_ONLY**, abstract confirms the reflection-invariance count via the Kreweras involution); Ding 2016 **[LIT]**; Burnside | P5 Part III: counts `[Catalan(w)+C(w,⌊w/2⌋)]/2 = 10,26,76,232,750` are **cited, not claimed** — "The paper cites; it claims nothing here" (#688 draft §5). Draft boundary §7.1 is explicit | **Uses** |
| D'Angeli–Donno 2013 (AAM 51(3):367–391; arXiv:1304.4180 — **ABSTRACT_ONLY**) | Fetched abstract confirms the load-bearing reading: "every lumping can be obtained from the action of a suitable subgroup … **if and only if** the poset … is totally ordered" — i.e. lumpings not arising from group actions exist, so "coarsest lumping = Aut-orbit" needs hypotheses P398 does not satisfy. P5 cites Prop. 13 as the *negative* boundary on its own w=4..8 instance | **Uses (as the reason the instance is not a theorem)** |
| Nonnegative vs ordinary rank of slack matrices (Yannakakis; Fiorini et al., regular n-gon **[LIT]**) | P5 Part I Thm D **reproduces the standard regular-n-gon separation** as a self-contained family (`rank ≤ 3`, nonnegative rank `n`) inside the response-matrix story; the draft itself says "Background the paper must cite but not claim" for the neighboring reduction | **Reproves a known instance as background** — must stay cited-not-claimed in prose |
| Probabilistic testing / copying tests (Larsen–Skou; van Breugel et al.; Deng–Feng **[LIT]**) | AP PDF §8.1 states the fork "is the same algebraic pattern used in probabilistic testing" and "the abstract ability of copying tests to refine linear trace information is therefore prior theory" | **Uses** |

**Published cousin of the no-go "bounded task rank does not determine a threshold"?**
Within what was fetched: **none found**. The closest fetched adjacency is D'Angeli–Donno
(lumping structure ≠ group structure) and the weighted-automaton canonical-form line
(Balle–Panangaden–Precup **[LIT]**, cited by the AP PDF itself as prior for
predictive/trace representations). No fetched source states "identical finite-horizon
responses, different closing-gap thresholds (`1/2` vs `1/3`)". This negative is bounded
to this pass's searches — hidden-Markov identifiability and interacting-particle-systems
literature were **not** swept and remain the open place a cousin could live. The #688
draft's own framing (Theorem E as a two-point separation with the honesty ladder) is the
correct hedge: even if a cousin exists, the torus/occupation connection keeps Part I
distinct, and Part II would still be an *explicit, verified* counterexample.

## 6. Q5 — Title/framing

The submitted title, *"Homological Occupation Growth: Branching Tests and Cut-Network
Representations"*, leads with a coinage ("homological occupation growth") that no
probability or combinatorics editor can map onto a known theorem class. What the paper
actually proves: (a) complete survival laws (lifetime signatures) do not determine
branching continuations, with gap 1/98; (b) the missing state is a cut-network
two-terminal connectivity state, update-closed; (c) every cofacial plane two-terminal
network realizes on a torus; (d) k+1 predictive classes on 8k vertices inside one
survival class. Suggested titles, theorem-naming, no marketing:

1. **"Complete survival laws do not determine branching continuations: cut-network
   states and a converse realization theorem for rank-one occupation growth on the
   torus."**
2. **"Lifetime signatures versus copying tests: an exact cut-network representation and
   unbounded branching-predictive classes on 8k vertices."**
3. **"Two-terminal network states for torus occupation growth: representation,
   realization, and the limits of reliability summaries."**

Each names either the separation (survival law vs branching) or the representation
theorem, in vocabulary an editor already parses. For the full P5 object (if/when Part II
is bundled), keep "no-go" language out of the title until the critical-sector-completeness
hedge is stated inside the abstract, not the introduction.

## 7. Stress-testing the prior Grok assessment

Point by point, with the fetched evidence:

1. **"The JAP desk reject is a correct audience call, not a finding that the
   constructions are false."** — **Agree**, and the fetch sharpens it: JAP's scope
   statement is application-broad but its sampled TOC (§2) is uniformly
   asymptotic-theorem-shaped; AP-2026-0311 has no asymptotic statement. The reject
   letter's own second clause matches the paper's §8.1 concession that the testing
   instrument is prior.
2. **"Methods are classical; packaging is new."** — **Agree with a refinement.** The
   *instrument* (copying tests, lumpability criterion, Kalman ranks) is classical and the
   paper says so. But Grok's phrasing undersells Thms 4.1/5.2: the cut-network
   representation + converse realization is a genuinely new two-sided bridge between
   rank-one homological continuation and planar two-terminal reliability. It is new
   *structural* methodology — it simply is not *applied-probability* methodology, which
   is what the EiC's clause was measuring. The right venue is one whose audience owns
   structure, not limits (§4: EJC first).
3. **"Split: cut-network rank + no bounded radius → CPC or RSA; Kalman no-go → EJP or
   ALEA; P398 quotient waits."** — **Partially agree; refine two ways.** (i) The
   submitted paper is *already* the Part-I unit; splitting "cut-network rank" from "no
   bounded radius" would thin an already-complete object — enrich, don't shatter. (ii)
   CPC/RSA: CPC's fetched scope lists graph theory, so it is *arguable*, but its
   probabilistic clauses are all random-structures/asymptotic; RSA could not be fetched,
   so no RSA recommendation is recorded here. EJC (fetched TOC, §4) is the better
   ticket-list target for Part I. (iii) EJP/ALEA for Part II: plausible, but the
   fetched EJP scope routes <12-page papers to ECP, and LAA's fetched scope makes LAA at
   least as good a fit because Theorem E is a Kalman-rank statement. (iv) P398 waits:
   **agree exactly** — five widths are an instance, the generic-rates lemma is unwritten,
   and D'Angeli–Donno (fetched) is precisely why it cannot be written as a general
   theorem.
4. **"The submitted title hid the theorem."** — **Agree**; §6 proposes replacements that
   name the separation and the representation theorem.
5. **"Not PTRF / Ann. Probab. / CMP."** — **Agree**, trivially: no asymptotics anywhere
   in the object; nothing fetched suggests otherwise. Not re-litigated.

## 8. Claim boundaries kept (restated for the record)

- P5 and AP-2026-0311 are one paper object; §3's split never double-submits theorems.
- The w=4..8 lumping–orbit equality is an **instance**, not a theorem (D'Angeli–Donno
  fetched and quoted as the reason).
- No wrapping identities (`D=r-1`, `F=(1+M)/2`) appear here or are recommended anywhere.
- No JAP/AAP appeal; no JAP/AAP resubmission of the same text (no fetched analogue;
  §2's analogue list is explicitly not construction-based).
- The AP proof PDF is not committed; no `references.bib` entries were added.
- `docs/STATUS.md` and `docs/ROADMAP.md` untouched; #644, #650, #688 untouched.

## 9. Not fetched / open for a follow-up pass

- **RSA scope** (Wiley page returned no text) — needed before any RSA call.
- **SCL, J. Phys. A, JSTAT scope** — needed before a Part II control-physics call.
- Primaries for Kemeny–Snell, Wonham, Ding 2016, Yannakakis / Fiorini et al.,
  Larsen–Skou, Deng–Feng, Balle et al. — marked **[LIT]** above; a venue-prep pass
  should fetch and transcribe them into `references.bib` per the #601 packet convention.
- Hidden-Markov / identifiability literature sweep for a published cousin of Theorem E.
