# P2 venue retrieval — independent primary-source pass (issue #694)

**Date:** 2026-09-09. **Machine:** local. **Scope:** retrieval and notes only. No census re-run, no
PSLQ re-run, no Huawei, no edits to `docs/STATUS.md` / `docs/ROADMAP.md` / the P2 draft.

**Input read in the clone:** `docs/manuscripts/p2-algebraic-exclusion/manuscript.md` (abstract, §1.3, §1.4,
§2.1, §4.2–4.3, §6.5–6.6, §8, Target-venue paragraph), `docs/manuscripts/p2-algebraic-exclusion/README.md`
(claims A–F′, open items), `data/literature_threshold_sources.json` (provenance manifest).

The "Prior Grok assessment" in the ticket was stress-tested, not copied. Every claim below carries a source
status tag: `PRIMARY_TEXT_READ` (full or substantive text fetched and read), `ABSTRACT_ONLY`
(publisher/abstract page fetched, full text not read), `[LIT]` (identified, not fetched), `[SECONDARY]`
(compilation page, not a primary source). URLs were fetched this session unless tagged otherwise.

---

## Q2. The four source papers and this journal

All four confirmed on fetched IOPscience article pages, plus the journal scope page.

| Source | Verified citation | Status |
|---|---|---|
| Jacobsen 2015 | Jesper Lykke Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley–Lieb algebras*, **J. Phys. A: Math. Theor. 48 (2015) 454003**, published 20 Oct 2015. DOI 10.1088/1751-8113/48/45/454003 | `PRIMARY_TEXT_READ` (abstract page; abstract states `p_c = 0.592 746 050 792 10(2)` for square-site percolation, `n_max = 21`) — https://iopscience.iop.org/article/10.1088/1751-8113/48/45/454003 |
| Mertens 2022 | Stephan Mertens, *Exact site-percolation probability on the square lattice*, **J. Phys. A: Math. Theor. 55 (2022) 334002**. DOI 10.1088/1751-8121/ac4195 | `PRIMARY_TEXT_READ` (arXiv 2109.12102 full text) — IOP page: https://iopscience.iop.org/article/10.1088/1751-8121/ac4195 |
| Yang–Zhou 2024 (Comment) | Yi Yang and Shuigeng Zhou, *Comment on 'Critical points of Potts and O(N) models…'*, **J. Phys. A: Math. Theor. 57 (2024) 258001**, published 6 Jun 2024. DOI 10.1088/1751-8121/ad4d2c | `ABSTRACT_ONLY` — https://iopscience.iop.org/article/10.1088/1751-8121/ad4d2c |
| Jacobsen 2024 (Reply) | Jesper Lykke Jacobsen, *Reply to Comment on 'Critical points of Potts and O(N) models…'*, **J. Phys. A: Math. Theor. 57 (2024) 258002**, published 6 Jun 2024. DOI 10.1088/1751-8121/ad4d33 | `ABSTRACT_ONLY` (IOP page fetched; abstract is a one-line cross-reference, also confirmed via Semantic Scholar) — https://iopscience.iop.org/article/10.1088/1751-8121/ad4d33 |

Mertens 2022 final estimates, verified from the arXiv full text (eqs. 25a/25b):

> `p_c = 0.592 746 050 786(3)  (p_med),   p_c = 0.592 746 050 60(8)  (p_cell).`

These match `data/literature_threshold_sources.json` exactly.

Yang–Zhou 2024 abstract (fetched from the IOP page), quoted verbatim in the decisive clause:

> "Our calculation shows that the current best result of pc=0.59274605079210(2) by Jacobsen (2015 J. Phys.
> A: Math. Theor. 48 454003) is incorrect and the corrected value should be 0.5927460507896(1)."

### J. Phys. A scope, from the publisher about page (fetched 2026-09-09)

From https://publishingsupport.iopscience.iop.org/journals/journal-of-physics-a-mathematical-and-theoretical/about-journal-physics-mathematical-theoretical/ :

> "To be acceptable for publication in the journal, papers must make significant, original and correct
> contributions to one or more of the topics within these sections. **Mathematical papers should be clearly
> motivated by actual or potential application to physical phenomena. Algorithmic papers are encouraged and
> should contain clear physical background and theoretical motivations.**"

Relevant section listings on the same page:

- Statistical physics — "exactly solvable models in statistical mechanics"; "statistical mechanics, phase
  transitions and critical phenomena"; "**numerical and computational methods, analysis of algorithms**";
  "lattice models, random walks and combinatorics".
- Mathematical physics — "**algebraic structures and number theory**"; "**numerical approximation and
  analysis**"; "**computational methods**".

Note: the plain IOPscience `/journal/1751-8121/about` URL 404s; the scope page above is the publisher's
current canonical "About" page (fetched successfully). The `/journal/1751-8121` home page resolves but
carries no scope text.

### Correction for the record: Scullard 2006 venue

`PRIMARY_TEXT_READ` — Christian R. Scullard, *Exact site percolation thresholds using a site-to-bond
transformation and the star-triangle transformation*, **Phys. Rev. E 73, 016107 (2006)**, published
10 Jan 2006. Not J. Phys. A. Abstract, quoted from the fetched APS full text:

> "The primary lattice studied here, the 'martini lattice,' is a hexagonal lattice with every second site
> transformed into a triangle. The site threshold of this lattice is found to be 0.764826…, i.e., the
> solution to p4 − 3p3 + 1 = 0, while the others have (√5 − 1)/2 (the inverse of the golden ratio) and 1/√2."

This is the primary the draft cites as `[LIT: primary text not verified]` in §1.1 for the martini-lattice
thresholds. The fetch confirms the threshold claims; the draft's `[LIT]` marker can be cleared by a future
draft-revision ticket (not done here — this ticket does not edit the draft's claims).

---

## Q1. Comparable published work (certified vs heuristic, decimal vs interval)

Three primary papers fetched and read; none of them certifies a null against an *interval*. That absence is
itself the finding: the manuscript's §1.2 diagnosis matches the primary texts, and its "as far as we know,
absent from the integer-relation literature" claim for sensitivity-certified nulls (§8, item 3) survived
this search.

**1. Bailey–Borwein–Girgensohn 1994 — `PRIMARY_TEXT_READ`.** David H. Bailey, Jonathan M. Borwein, Roland
Girgensohn, *Experimental evaluation of Euler sums*, **Experimental Mathematics 3(1) (1994) 17–30**.
Full text fetched from https://www.davidhbailey.com/dhbpapers/eulsum-em.pdf (journal copy:
https://www.tandfonline.com/doi/abs/10.1080/10586458.1994.10504573). Content: PSLQ-based identification of
Euler sums as rational linear combinations of zeta values and log 2. Certification status, quoted:

> "It should be emphasized that the results in Table 3 are not established in any rigorous mathematical
> sense by these calculations. However, in each case the 'confidence level' (see Section 3) of these
> detections is better than 10^-50, and in most cases is in the neighborhood of 10^-100."

Target form: a single high-precision decimal of an exact constant — not an interval, not a measured
physical quantity. Heuristic, positive-only.

**2. Bailey–Broadhurst 2001 — `PRIMARY_TEXT_READ`.** David H. Bailey, David J. Broadhurst, *Parallel
Integer Relation Detection: Techniques and Applications*, **Mathematics of Computation 70 (2001)
1719–1736**. Full text fetched from https://www.davidhbailey.com/dhbpapers/ppslq.pdf (AMS page:
https://www.ams.org/journals/mcom/2001-70-236/S0025-5718-00-01278-3/). Content: parallel PSLQ
("PSLQ++"), applications to QFT sums and a chaos-theory constant. Certification status, quoted:

> "The ratio between the smallest and the largest y entry when a relation is detected can be taken as a
> 'confidence level' that the relation is a true relation and not an artifact of insufficient numeric
> precision. Very small ratios at detection, such as 10^-100, almost certainly denote a true relation,
> although of course such results do not constitute a rigorous proof."

Note the paper's framing of what an integer-relation algorithm *can* do in the null direction: "or can
produce bounds within which no integer relation exists" — but the applications in the paper are all
detections, and no certified non-existence theorem over an enumerated class is reported. Target form:
high-precision decimals of exact constants; no interval.

**3. Scullard 2006 — `PRIMARY_TEXT_READ`** (see Q2). Comparable in *spirit* rather than method: it
certifies exact thresholds by construction (site-to-bond + star–triangle), i.e. the kind of certified
result about percolation thresholds that physics venues publish. It performs no search against a measured
constant.

**4. Jacobsen 2015 — `PRIMARY_TEXT_READ` (abstract)** is the nearest in-field analogue: a certified
transfer-matrix/eigenvalue *computation* feeding a *non-certified extrapolation* to produce the target
itself. Useful for the venue argument; not a relation search.

**[LIT] identified, not fetched (recorded so nobody has to re-derive them):**

- W. P. Orrick, B. Nickel, A. Guttmann, J. H. H. Perk, *The susceptibility of the square lattice Ising
  model: new developments*, **J. Stat. Phys. 102 (2001) 1113** — identification of exact Ising
  susceptibility structure from high-precision numerics; candidate J. Stat. Phys. analogue, unfetched.
- Heuristic decimal-target identifications of physical constants (fine-structure-constant style claims)
  exist widely in the literature but none certified; none fetched, so none is quoted here.

**Summary answer to Q1.** The published tradition is (a) PSLQ/integer-relation *detection* at a rounded or
high-precision decimal of an exact constant, explicitly non-rigorous (both Bailey-school quotes above);
(b) certified *derivations* of exact thresholds (Scullard-type), which search nothing. A **complete
certified census of an enumerated polynomial class against a measured constant at its quoted precision**
has no fetched precedent, and none of the comparables searches an interval rather than a decimal. This
confirms the manuscript's §1.2 framing and supports its novelty claim for the boundary-degree criterion
(§4.3) and the planted-root null control (§6.4).

---

## Q3. Is the disjointness remark known?

**Refine, in both directions.**

*What the literature already states (pairwise, informally):*

- Mertens 2022 (`PRIMARY_TEXT_READ`, arXiv full text, after eq. 25):
  > "The first value deviates from the reference value (19) by 2 errorbars, the second by 2.5 errorbars."
  where (19) is Jacobsen's `0.592 746 050 792 10(2)`. Mertens also reports his `p_pol` estimator,
  `p_c = 0.592 746 050 792 0(4)` (eq. 26), which "agrees with Jacobsen's value (19) within the errorbars".
- Yang–Zhou 2024 (`ABSTRACT_ONLY`) flatly asserts Jacobsen 2015 "is incorrect", with corrected value
  `0.5927460507896(1)`; Jacobsen 2024's Reply (`ABSTRACT_ONLY`) exists — the pair is a published,
  unresolved dispute inside one journal issue.

So "some of these estimates disagree" is **not** new: two of the four source papers say so about specific
pairs, in print.

*What this search did not find:* any source that assembles the intervals at their quoted precisions and
states the **four-way pairwise disjointness**, or states the consequence that **no pooled interval is
defensible and every exact-search claim must be reported per interval**. The community compilation
practice, checked on the Wikipedia *Percolation threshold* page (`[SECONDARY]`,
https://en.wikipedia.org/wiki/Percolation_threshold, fetched 2026-09-09), is to list fifteen square-site
determinations side by side in one table row — including `0.59274605079210(2)` [Jacobsen15],
`0.592746050786(3)` [ref 28], `0.5927460507896(1)` [ref 29] and an additional `0.59274605079016(1)`
[ref 30, identity not confirmed in this pass] — without remarking that several of these mutually exclude
each other at their own quoted uncertainties. Listing-not-pooling is the implicit acknowledgment of the
problem; no fetched source makes the methodological statement.

**Refinement for the draft (goes to the PR discussion, not the draft):** §2.1's novelty should be worded
as *systematizing a disagreement each author noted only pairwise, at one estimator's own precision*, not
as a first observation that the estimates differ. Mertens's own `p_pol` value agreeing with Jacobsen is
worth a sentence in that discussion: the disjointness is between *estimator families*, and one of
Mertens's three estimators sits with Jacobsen — which strengthens, rather than weakens, the
"targets must be plural" conclusion of §8 item 1, since even one paper's internal estimators span
incompatible intervals.

No consensus `p_c` value pooling these four was found in this pass.

---

## Q4. Venue fit, independently assessed

**Journal of Physics A: Mathematical and Theoretical — fit (primary).**
Independent grounds, beyond the Grok argument:
- The publisher scope page *explicitly* invites exactly this hybrid: "Algorithmic papers are encouraged
  and should contain clear physical background and theoretical motivations", and both relevant sections
  list "numerical and computational methods, analysis of algorithms" (Statistical physics) and
  "algebraic structures and number theory" / "computational methods" (Mathematical physics). A certified
  census against a percolation threshold sits in the intersection of those listings.
- Closest analogue: **Mertens 2022** — a single-author exact-algorithm paper about square-site percolation
  whose estimates section is the paper's physics payoff — same journal, same constant, same structure
  (algorithm + the constant). `PRIMARY_TEXT_READ`.
- All five documents the paper corrects or consumes (2015 paper, 2022 paper, 2024 Comment, 2024 Reply)
  are in this journal, verified above.

**Experimental Mathematics — fit (fallback, genuine).** Scope, quoted from the publisher aims page
(`ABSTRACT_ONLY`, https://amstat.tandfonline.com/action/journalInformation?show=aimsScope&journalCode=uexm20):
> "Experimental Mathematics publishes original papers featuring formal results inspired by
> experimentation, conjectures suggested by experiments, and data supporting significant hypotheses."
> … "New theorems proved with the help of experimental results are highly acceptable" … "Computer
> experiments should be reported in such a way that they can be repeated by other researchers."

The census is a formal result *of* a computer experiment with a repeatability supplement (§9) — the scope
fits almost verbatim. Closest analogue: **Bailey–Borwein–Girgensohn 1994, Exp. Math. 3(1) 17–30**
(`PRIMARY_TEXT_READ`) — but note it is heuristic and positive; **neither venue has a published precedent
for a certified *null* integer-relation census**, which is exactly the gap this paper occupies wherever
it lands.

**Mathematics of Computation — stretch.** Analogue: Bailey–Broadhurst 2001, Math. Comp. 70 (2001)
1719–1736 (`PRIMARY_TEXT_READ`) — but that paper is an *algorithm* paper (parallel PSLQ) with
applications; MoC would want the screening method itself as the contribution, and Theorem 2 is a bound
whose proof is elementary (MVT + rounding). The percolation target would read as application. Possible,
not recommended first.

**Journal of Statistical Physics / JSTAT / PRE — wrong as written, with one nuance.** Nuance: PRE *does*
publish exact-threshold derivations (Scullard 2006 is PRE 73, 016107, verified above), so the *field*
publishes there. But this manuscript derives nothing physical: no mechanism, no new estimate, no critical
phenomenon — it certifies the absence of a form. Closest J. Stat. Phys. analogue: Orrick–Nickel–
Guttmann–Perk 2001 (`[LIT]`, unfetched — identification work, and 25 years old; no recent analogue
found). JSTAT: no analogue found.

**Journal of Algebra / number-theory journals — wrong.** No analogue found in this pass: the exclusion
theorem's content is entirely about the complexity class vs a physics constant; the number-theoretic
skeleton (Prop. 1, Thm 2) is elementary, and nothing in the paper is a statement about polynomials as
such. Agree with Grok point 5, now with the scope-page asymmetry as evidence: JPA's own scope lists
"algebraic structures and number theory" *motivated by physical phenomena*; a pure venue has no such
bridge.

---

## Q5. Desk-reject risks — what an editor would actually say

- **"Insufficient physics" (JPA).** Real but mitigable: the disjointness finding (§2.1) is a
  physics-literature correction, and JPA published the pure-algorithm Mertens 2022. The framing note in
  the draft's Target-venue paragraph (lead with the census, not the methodology) is the right defense.
  Moderate risk, not a desk-reject.
- **"The class is arbitrary."** Real at every venue, and the draft already carries the two-part defense
  (§4.2.1 concedes the historical range is a selected sample; §4.3 supplies a literature-free stop rule).
  An editor who reads only the abstract can still bounce it; the abstract should keep the
  "degree 4 = boundary degree" sentence, since it is the answer to the objection.
- **"Negative computational result."** Real at Math. Comp. and pure venues; explicitly *not* real at
  Experimental Mathematics (scope quote above welcomes formal results from experiments); at JPA it
  converts into the disjointness correction. Venue-dependent, not intrinsic.
- **"The height bound is already broken by bow-tie thresholds."** True as a fact (§4.2: Wierman 1984
  height 6; Ziff–Scullard 2006 degree 11 at height 36) and already conceded in the draft — §4.2, §6.6
  and §8.1 all restate it. A skimming referee could still raise it; a footnote-level pointer in the
  abstract ("no form at the complexity of any *known exact planar threshold*…" is already gone from the
  draft — good) is the only remaining exposure. The draft's withdrawn phrasing ("we withdraw that
  phrasing", §6.6) means this objection no longer has a target in the text.

---

## Stress-test of the prior Grok assessment, claim by claim

1. **"First venue should be J. Phys. A because three of the four intervals come from that journal."**
   **Agree, with a stronger and a weaker finding.** Stronger: all four *plus* the Reply are verified JPA
   (Q2 table), and the publisher scope page explicitly solicits algorithmic papers with physical
   background — Grok's argument rests on reader proximity alone and misses that the journal's stated
   scope covers the paper's genre. Weaker: the "correction does not reach the people it is about"
   premise needs care — the *disagreement* has already reached them (Mertens's 2/2.5-errorbar remark;
   the Comment/Reply exchange). What has not reached anyone is the *systematic four-way statement and
   its no-pooling consequence*. JPA remains right; the novelty pitch should be the systematization
   (see Q3 refinement).
2. **"Fallback Experimental Mathematics."** **Agree**, now grounded in the fetched aims-and-scope page:
   "formal results inspired by experimentation", repeatability expectations, and the precedent of BBG
   1994. One refinement: Grok says the boundary-degree criterion is "closer to ExpMath's centre of
   gravity" — supported, but the deeper fit is that ExpMath explicitly publishes experiments *without*
   complete formal closure, so the manuscript's sensitivity-certified nulls (§6.4) and frozen-contract
   method (§5) are themselves genre-expanding there.
3. **"A respectable specialized paper, not CMP/PTRF/Annals; not a transcendence result."** **Agree.**
   No fetched source contradicts; the theorems are elementary (Prop. 1, Thm 2 verified in §3 of the
   draft), the claim boundary (§1.4, §8) explicitly forswears transcendence, and nothing fetched elevates
   the census above specialized.
4. **"Methodological pieces are the durable novelty; the physics-facing news is the pairwise disjoint
   intervals."** **Refine.** The methodological pieces (plural targets, cross-interval resolution,
   planted-root sensitivity, boundary degree) are durable, and the retrieval found no precedent for any
   of them (Q1). But the disjointness is *partially* pre-known: Mertens 2022 states the pairwise
   deviations in print, and the Comment/Reply dispute is public. The news is the four-way exact statement
   + "no pooled interval is defensible", and the retrieval strengthens the methodological side
   accordingly: the transferable toolkit is the part with no literature at all.
5. **"Not to PRE as a physics-mechanism paper; not to pure number-theory journals."** **Agree**, with the
   Scullard nuance: PRE publishes exact-threshold *derivations* (Scullard 2006 is PRE), so the objection
   to PRE is not "PRE never publishes threshold papers" but "this paper contains no mechanism or
   derivation, which is what PRE threshold papers have". J. Algebra: no analogue found; wrong.

---

## Source log (everything fetched this session)

| # | Source | URL | Status |
|---|---|---|---|
| 1 | Jacobsen 2015, JPA 48 454003 | https://iopscience.iop.org/article/10.1088/1751-8113/48/45/454003 | PRIMARY_TEXT_READ (abstract page) |
| 2 | Mertens 2022, JPA 55 334002 (IOP page) | https://iopscience.iop.org/article/10.1088/1751-8121/ac4195 | ABSTRACT_ONLY (page) |
| 3 | Mertens 2022 (arXiv 2109.12102 full text) | https://arxiv.org/pdf/2109.12102 | PRIMARY_TEXT_READ |
| 4 | Yang–Zhou 2024 Comment, JPA 57 258001 | https://iopscience.iop.org/article/10.1088/1751-8121/ad4d2c | ABSTRACT_ONLY |
| 5 | Jacobsen 2024 Reply, JPA 57 258002 | https://iopscience.iop.org/article/10.1088/1751-8121/ad4d33 | ABSTRACT_ONLY (page + Semantic Scholar record) |
| 6 | J. Phys. A about/scope page | https://publishingsupport.iopscience.iop.org/journals/journal-of-physics-a-mathematical-and-theoretical/about-journal-physics-mathematical-theoretical/ | PRIMARY (publisher page) |
| 7 | Scullard 2006, PRE 73 016107 | https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevE.73.016107/fulltext | PRIMARY_TEXT_READ |
| 8 | Bailey–Borwein–Girgensohn 1994, Exp. Math. 3(1) 17–30 | https://www.davidhbailey.com/dhbpapers/eulsum-em.pdf | PRIMARY_TEXT_READ (author copy) |
| 9 | Bailey–Broadhurst 2001, Math. Comp. 70 1719–1736 | https://www.davidhbailey.com/dhbpapers/ppslq.pdf | PRIMARY_TEXT_READ (author copy) |
| 10 | Experimental Mathematics aims & scope | https://amstat.tandfonline.com/action/journalInformation?show=aimsScope&journalCode=uexm20 | ABSTRACT_ONLY (publisher page) |
| 11 | Wikipedia, Percolation threshold (compilation practice) | https://en.wikipedia.org/wiki/Percolation_threshold | [SECONDARY] |
| 12 | Orrick–Nickel–Guttmann–Perk 2001, J. Stat. Phys. 102 1113 | not fetched (Springer/OpenAlex/S2 attempts blocked or rate-limited) | [LIT] |

Blocked/not-fetched notes: `iopscience.iop.org/journal/1751-8121/about` 404s (scope obtained from the
publisher support page instead, #6); `link.springer.com` blocked the Orrick fetch (IP block);
Semantic Scholar/OpenAlex search endpoints rate-limited during the Orrick query. Nothing else was
blocked; no NEED_HUAWEI was required.

**Bottom line.** The tripwire is met: one fetched publisher scope page (#6, plus #10), and four fetched
comparable papers with title/year/venue/content (#7, #8, #9, and #3 as the in-field algorithmic
analogue). The draft's target-venue decision (JPA primary, Exp. Math. fallback) is **confirmed
independently**, with one substantive refinement: §2.1's disjointness finding should be positioned as a
systematization of an in-literature disagreement, not a first observation of it.
