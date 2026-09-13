# Venue and literature retrieval for the P3 draft ("Test the ray, not the ratio")

**Date:** 2026-09-09
**Ticket:** #695 (parent #650; manuscript ticket #579)
**Kind:** retrieval and notes only. No Monte Carlo, no rescoring, no Huawei, no edits to the
manuscript's scientific claims, no `docs/STATUS.md`.
**Method:** primary sources were fetched over the web on 2026-09-09. Each source is marked
`PRIMARY_TEXT_READ` (full text or page fetched and the quoted sentence read in it),
`ABSTRACT_ONLY` (only the abstract fetched), `[LIT]` (not fetched; known from a secondary
source's reference list), or `API_QUERY` (search-API hit counts, not full text). No quotation
below is from memory: every quoted string was read in the fetched document named beside it.

---

## 1. What the draft claims (the target of this retrieval)

§8 of `docs/manuscripts/p3-projective-inference/manuscript.md` states:

> "The practice of testing amplitude ratios against predicted ratios is ubiquitous — it is how
> conformal-invariance predictions for aspect-ratio and boundary-condition dependence are compared
> with lattice data [6] **[LIT]**"

and reference [6] is the draft's only marked hole: Cardy 1984–1986, "specific citation to be
fixed against the primary." The ticket asks whether primary literature actually **ratio-tests**
finite-size amplitudes (forms `Y/X` from measured amplitudes and tests it against a predicted
ratio), and whether the word "ubiquitous" survives contact with that literature.

## 2. Q1 — papers that form a ratio of finite-size amplitudes and test it (the [6] hole)

**Headline result: the practice is real, but it is not "ubiquitous", and the draft's sentence
must be weakened** (§6 below). The dominant framing in the crossing-probability literature is a
**full-curve comparison** of the measured probability against the conformal formula as a function
of aspect ratio; explicit ratio-of-amplitudes tests exist but are the exception, and most sit in
the universal-amplitude-ratio literature rather than the aspect-ratio literature.

### 2.1 Ratio tests, primary text read

**R1. Simmons, Kleban & Ziff (2007), "Percolation crossing formulae and conformal field theory",
*J. Phys. A: Math. Theor.* **40**, F771. `PRIMARY_TEXT_READ`**
(ar5iv full text of arXiv:0705.1933; venue confirmed on the IOPscience page
https://iopscience.iop.org/article/10.1088/1751-8113/40/31/F03).

The cleanest genuine ratio test found: a measured fraction of Monte Carlo events is compared
against a predicted **ratio of two integrated conformal densities**. Verbatim, from the
numerical-verification section:

> "The above fraction of multiple crossing events, 0.07018, is somewhat above the predicted
> value 0.069189, which is found by integrating the formulas for ν_h and π̄ᵇ_h and taking the
> ratio of the integral of the former to the sum of the integrals of the former and latter."

and, on the finite-size extrapolation of that ratio:

> "The data fit very well to a straight line when plotted as a function of 1/L, with an intercept
> of 0.06928, quite close to the predicted value."

Measured `Y/(X+Y)` vs predicted ratio, extrapolated to the infinite-lattice limit — a ratio test
in exactly the sense the draft describes. It is also honest evidence for the counter-claim below,
because the same paper's main test is a full-curve comparison ("The relative difference between
the two curves is on the order of 2%…").

**R2. Izmailian & Yeh (2009), "Ising model with mixed boundary conditions: universal amplitude
ratios", *Nucl. Phys. B* **814**, 573. `PRIMARY_TEXT_READ`**
(ar5iv full text of arXiv:1005.1712; journal reference read on the arXiv abs page).

Boundary-condition amplitude ratios, formed and tested against a conformal prediction. Verbatim:

> "In this paper we present the new set of the universal amplitude ratios for the mixed boundary
> universality class. The results are in perfect agreement with a perturbated conformal field
> theory scenario proposed by Cardy [cardy86]."

The ratios are of the form `r_s(k) = b_k/a_k` (correlation-length over free-energy first-order
amplitudes): a ratio of two finite-size amplitudes with different boundary conditions, compared
with the ratio predicted by perturbed CFT. Caveat stated honestly: these amplitudes come from
the **exact solution** (Euler–Maclaurin of exact eigenvalues), not from Monte Carlo measurement.

**R3. Kamieniarz & Blöte (1993), "Universal ratio of magnetization moments in two-dimensional
Ising models", *J. Phys. A* **26**, 201. `ABSTRACT_ONLY`**
(IOPscience abstract page
https://iopscience.iop.org/article/10.1088/0305-4470/26/2/009/meta; full text paywalled).

Verbatim from the abstract:

> "The authors calculate the universal critical-point ratios of the square of the second and the
> fourth moment of the magnetization for ferromagnetic Ising models on the square and on the
> triangular lattices. … This analysis is also applied to rectangular systems with arbitrary
> aspect ratios."

A ratio of two finite-size amplitudes from transfer-matrix finite-size data, including
aspect-ratio-dependent systems, against a universal value. Because only the abstract was read,
this citation cannot supply a verbatim `Y/X` sentence; it is included as ABSTRACT_ONLY support,
not as one of the three load-bearing quotations.

### 2.2 Primary papers that compare with lattice data but do **not** ratio-test

These carry the "conformal-invariance predictions … are compared with lattice data" half of the
draft's sentence — and they show that the comparison is normally a curve, not a ratio.

**C1. Langlands, Pouliot & Saint-Aubin (1994), "Conformal invariance in two-dimensional
percolation", *Bull. Amer. Math. Soc.* **30**, 1–61. `PRIMARY_TEXT_READ`**
(PDF fetched from https://publications.ias.edu/sites/default/files/conf-invar-ps.pdf).

Their §3.2 is literally titled "Experimental verification of Cardy's formula." Verbatim:

> "The goal of the first experiment is twofold: to verify again Cardy's prediction for the
> function π_h(r) on M0 and to obtain values of π_hv(r, M0) suitable for comparison in other
> experiments."

> "Figure 3.2. Comparison of 81 measured values of ln π_h/(1 − π_h) (dots) with Cardy's
> prediction (curve)."

> "Despite the systematic error, the agreement is remarkable and we shall compare the results of
> the following experiments with π^cft instead of π̂ when the former is applicable."

The comparison is of the measured function of aspect ratio against the predicted function; the
closest thing to a ratio is a log-odds transform of a single probability, not a ratio of two
amplitudes.

**C2. Watts (1996), "A crossing probability for critical percolation in two dimensions",
*J. Phys. A* **29**, L363. `PRIMARY_TEXT_READ`** (arXiv HTML of cond-mat/9603167).

Watts compares his formula for π_hv against the Langlands *et al.* Monte Carlo data. Verbatim:

> "Comparing this function with the numerical data for π_hv presented in [2], our agreement is as
> good as that of Cardy's formula for π_h."

> "In figure 1 we plot log(π_hv) against log(r) for both the numerical results obtained in [2]
> and for the function 𝓕(x) above. The agreement is excellent."

Full-curve comparison again. Watts also records the historical fact that Cardy's fit to LPS's
data was so good that the numericians adopted it in place of their own data — which cuts against
any claim that ratio-testing is the routine quantitative step:

> "In fact, this fit was so good, that Langlands et al. decided to adopt Cardy's result as π_h,
> rather than their own numerical data, in their later paper on more complicated calculations [2]."

**C3. Ziff (1995), "On Cardy's formula for the critical crossing probability in 2D percolation",
*J. Phys. A* **28**, 6479. `ABSTRACT_ONLY`** (IOPscience page
https://iopscience.iop.org/article/10.1088/0305-4470/28/5/013; DeepBlue full text blocked by
Cloudflare, IOP full text paywalled). The abstract establishes that the paper expands Cardy's
formula as a function of aspect ratio; the full text could not be read, so it is used here
neither for nor against.

**[LIT] Langlands, Pichet, Pouliot & Saint-Aubin (1992), "On the universality of crossing
probabilities in two-dimensional percolation", *J. Stat. Phys.* **67**, 553** — known here from
Watts (1996)'s reference list [2]; not fetched. The multi-lattice-family universality test and
the natural next primary to read if [6] is rebuilt.

**[LIT] Pinson (1994), "Critical percolation on the torus", *J. Stat. Phys.* **75**, 1167** —
the torus-amplitudes companion to LPS. Springer's CDN blocked both PDF and HTML fetches from
this machine on 2026-09-09 (HTTP 406/403, repeated), and Semantic Scholar records
`openAccessPdf: CLOSED`. Not read; listed so the reader knows it was attempted and is the next
thing to fetch for the boundary-condition half of [6].

### 2.3 Verdict on the [6] hole

The ticket asked for three primary papers that *actually ratio-test*. Honest scoring:

| candidate | ratio of measured amplitudes vs predicted ratio? | status |
|---|---|---|
| Simmons–Kleban–Ziff 2007 | **yes**, explicitly (measured event fraction vs predicted ratio of integrated densities, extrapolated in 1/L) | PRIMARY_TEXT_READ |
| Izmailian–Yeh 2009 | **yes**, but the amplitudes are exact-solution, not simulated | PRIMARY_TEXT_READ |
| Kamieniarz–Blöte 1993 | **yes** (universal ratio from finite-size data incl. aspect ratios), but abstract only | ABSTRACT_ONLY |
| LPS 1994; Watts 1996; Ziff 1995 | no — full-curve comparisons | PRIMARY_TEXT_READ ×2, ABSTRACT_ONLY ×1 |

Three papers can be cited — but the third only on its abstract, and two of the three are not
aspect-ratio ladder measurements in the draft's sense. What the primary reading actually shows
is that the **routine** comparison in this literature is the full curve of a probability or
amplitude against its conformal prediction as a function of aspect ratio (LPS 1994; Watts 1996),
with ratio-of-amplitudes tests appearing as secondary or special-purpose checks. The draft's
"ubiquitous" therefore does not survive, and §8 must be weakened as in §6.

## 3. Q2 — is the projective statistic new?

**Bottom line: the exact statistic as a *test* appears new to the FSS/percolation literature,
but narrowly, and the draft must engage four pieces of nearby prior art it does not currently
cite.**

**Negative evidence (novelty survives):**

- Inspire HEP API, query `Fieller confidence interval`: **0 hits**. `API_QUERY`
- arXiv metadata search, query "Fieller": 4 hits, all statistics — Owen (2025,
  arXiv:2510.00389), Tang (2021, arXiv:2110.12636), von Luxburg & Franz (2009, arXiv:0711.0198),
  Franz (2007, arXiv:0710.2024). None in FSS/percolation/stat-mech. `API_QUERY` /
  `ABSTRACT_ONLY`
- "Fieller" + "scaling amplitude": 0 hits; "Hotelling" + "finite-size scaling": 0 hits;
  "finite-size scaling" + covariance + amplitude fits: 0 hits. `API_QUERY` (arXiv search is
  metadata-only, so absence is weak-but-indicative)
- Gleser–Hwang (1987): only statistics sources found; no physics citation evidence. Negative
  result.
- `pyfssa` documentation, https://pyfssa.readthedocs.io/en/stable/reference/fssa.html.
  `PRIMARY_TEXT_READ`. The standard FSS data-collapse tool's quality statistic is
  covariance-free: "This is the reduced chi-square statistic for a data fit except that the
  master curve is fitted from the data itself" — its parameters are `x, y, dy` per-point
  standard errors, with no covariance input (quality function from Houdayer & Hartmann,
  *Phys. Rev. B* **70**, 014418 (2004)).

**Prior art that must be cited and distinguished:**

**P1. Parisen Toldin (2011), "Improvement of Monte Carlo estimates with covariance-optimized
finite-size scaling at fixed phenomenological coupling", *Phys. Rev. E* **84**, 025703(R).
`PRIMARY_TEXT_READ`** (arXiv abs + HTML of arXiv:1104.2500; quotes independently re-verified
directly in this retrieval). Verbatim (abstract):

> "Within this scheme of finite-size scaling, we exploit the statistical covariance between the
> observables in a Monte Carlo simulation in order to reduce the statistical errors of the
> quantities involved in the computation of the critical exponents."

Verbatim (method section; the exact phrase "orthogonal projection" does not appear, but the
orthogonal/parallel decomposition does):

> "A positive-definite scalar product <,> between two vectors δÂ and δB̂ can be defined by the
> covariance <δÂ, δB̂> ≡ COV[δÂ, δB̂] = E[δÂ δB̂] − E[δÂ] E[δB̂]"

with δO decomposed into components "orthogonal δO⊥ and parallel δO∥ to the space 𝒱". This is the
covariance-metric geometry of the draft's statistic, already inside FSS — but it *minimizes the
variance of an error-propagated observable* (norm minimization under the covariance metric,
projection onto an affine subspace of dimension N−1, not a ray): no χ², no hypothesis test, no
null inversion. The draft's claim survives as "the projective statistic as a *test statistic*
is new", and §1.4/§8 should cite Parisen Toldin (2011) as the nearest neighbour. The abstract's
"the natural statistic is then…" should acknowledge that a covariance scalar product is already
in FSS use for error minimization.

**P2. von Luxburg & Franz (2009), "A Geometric Approach to Confidence Sets for Ratios: Fieller's
Theorem, Generalizations, and Bootstrap", *Statistica Sinica* **19**, 1095; arXiv:0711.0198.
`ABSTRACT_ONLY`.** A geometric generalization of Fieller in the statistics literature; the
draft's §2.4 (confidence regions for directions) and §1.4 ("it is not new mathematics") should
name it, because it is the statistics-side generalization closest to what the draft builds.

**P3. Bruno & Sommer (2023), "On fits to correlated and auto-correlated data", *Computer Physics
Communications* **285**, 108643. `PRIMARY_TEXT_READ`** (abstract + introduction + conclusions
fetched; PDF paywalled). Verbatim (abstract):

> "Standard χ2 tests require a reliable determination of the covariance matrix and its inverse
> from correlated and auto-correlated data, a challenging task often leading to close-to-singular
> estimates. These motivate modifications of the definition of χ2 such as uncorrelated fits."

Correlated-χ² fitting is standard practice in lattice QCD; this supports the draft's framing of
the pseudo-inverse/rank-reporting discipline (§2.2) rather than competing with it.

**P4. Wang, Iyer & Hannig (2003), "Uncertainty Calculation for the Ratio of Dependent
Measurements", *Metrologia*. `PRIMARY_TEXT_READ`** (NIST publication page). Verbatim (abstract):

> "Although an exact confidence interval procedure, known as the Fieller interval, is available
> for this problem, practitioners often use various non-exact methods."

Fieller is alive in metrology but absent from physics FSS — which is exactly the draft's premise,
now documented rather than assumed.

**Standard-practice control.** Mertens, Jensen & Ziff (2017), "Universal features of cluster
numbers in percolation", *Phys. Rev. E* **96**, 052119; arXiv:1602.00644. `PRIMARY_TEXT_READ`.
Verbatim: "The metric factor b cancels out in the dimensionless ratio R = A₁C₁/B₁² =
f(0)f″(0)/2f′(0)² which is predicted to be universal for systems of a given shape." The
dimensionless-ratio practice the draft criticizes, in the percolation literature, formed by
cancelling the denominator — no covariance, no ray test anywhere in the fetched text.

## 4. Q3 — annihilating functionals / divided differences as class-level tests

**Bottom line: honest null.** Searches run: OpenAlex fulltext "divided differences" +
"finite-size scaling" (4 hits, all quantum Monte Carlo papers with no FSS-amplitude content);
arXiv metadata "divided differences" + "scaling amplitudes" (0 hits); web searches for
divided-difference or linear-combination class tests of critical/scaling amplitudes (only
generic scaling reviews and Kenna-type hatted-exponent amplitude **relations**, which are built
from ratios, not from a functional annihilating a class); "divided differences" lattice operator
(numerical-analysis and lattice-quasi-PDF usages only). `API_QUERY` and search-result level —
the only `[LIT]`-grade items here are the Kenna-type amplitude-relation papers
(arXiv:2402.00427 / *Entropy* **26**, 221 (2024)), which are ratio-based, not annihilator-based.

The closest construct in the fetched literature is Izmailian & Yeh (2009)'s universal linear
combinations `r_s(k) = b_k/a_k` — but those test individual predicted ratios, not a functional
that annihilates an entire class whatever the amplitudes. The draft's §8 sentence "we have not
seen it deployed as a class-level tool for finite-size ladders" stands as phrased ("we have not
seen"), and should not be strengthened to a proven negative: Google Scholar and the paywalled
PRB/PRE archives were not fully searchable.

## 5. Q4/Q5 — venue fit and the one-dataset reanalysis precedent

### 5.1 Q4: journals (scope sentences verbatim from the publisher pages named)

**Physical Review E — FIT.** Scope (journals.aps.org/pre/about):

> "Physical Review E is a trusted and interdisciplinary journal for high-quality, significant
> developments in the interrelated areas of statistical, nonlinear, complex systems, biological,
> soft matter, fluids, plasma and computational physics, including machine learning and
> artificial intelligence."

Acceptance criterion (same page): "Add to the existing knowledge related to statistical,
nonlinear, biological, and soft matter physics." A covariance-weighted χ² test of FSS amplitude
vectors in lattice percolation sits in PRE's statistical/computational core, and PRE explicitly
publishes reanalyses (see §5.2). **Analogue:** Jensen & Ziff, "Universal amplitude ratio Γ−/Γ+
for two-dimensional percolation", *Phys. Rev. E* **74**, 020101(R) (2006) — an amplitude-ratio
determination paper in PRE (`ABSTRACT_ONLY`, arXiv:cond-mat/0607146).

**JSTAT — FIT.** Scope (IOP publishing-support page, the direct /about page being blocked):

> "JSTAT is targeted to a broad community interested in different aspects of statistical
> physics, which are roughly defined by the fields represented in the conferences called
> 'Statistical Physics'."

Topically exact. Closest verified analogue: Xun *et al.*, *J. Stat. Mech.* (2025) 123301
(extended-range percolation, extensive MC + FSS); **no true amplitude-ratio-testing JSTAT
analogue was located — noted honestly.**

**J. Stat. Phys. — MILD STRETCH.** Springer's aims-and-scope page was bot-blocked; the only
description obtained is Wikipedia's ("All areas of statistical physics as well as related
fields concerned with collective phenomena in physical systems are covered"), which is **not a
publisher quote**. Topically a fit, but JSP's culture favours formal theory over a
Monte-Carlo-reanalysis methodology paper.

**Computer Physics Communications — STRETCH.** Scope (sciencedirect.com
aims-and-scope): "Computer Physics Communications publishes research papers and application
software in the broad field of computational physics" … "The focus of CPC is on contemporary
computational methods and techniques and their implementation, the effectiveness of which will
normally be evidenced by the author(s) within the context of a substantive problem in physics."
Acceptable only if the method ships as software with the paper. Analogue: Hu, Chen, Izmailian &
Kleban, *CPC* **126**, 77–81 (2000) (Monte Carlo methods for percolation).

**The American Statistician — STRETCH.** Scope (journal landing page): "articles that address
practical problems at the interface between statistical methodology and areas of application."
The generalized-Fieller statistic fits that brief, but the physics framing must lead or the
audience match fails; and the paper has no simulation study, so a statistics journal would
treat it as a note. Analogue: Hirschberg & Lye, "A Geometric Comparison of the Delta and
Fieller Confidence Intervals", *The American Statistician* **64**(3), 234–241 (2010).

**J. Phys. A — FIT.** Scope (IOP publishing-support page): acceptable papers must contribute to
listed sections, which include verbatim "statistical mechanics, phase transitions and critical
phenomena" and "numerical and computational methods, analysis of algorithms". Analogue: Delfino,
Viti & Cardy, "Universal amplitude ratios of two-dimensional percolation from field theory",
*J. Phys. A* **43**, 152001 (2010).

### 5.2 Q5: does PRE accept a reanalysis of one dataset that inverts a verdict?

**Yes — precedent exists.** Three fetched examples:

1. **Wang, Zhou, Zhang, Garoni & Deng (2013), "Bond and site percolation in three dimensions",
   *Phys. Rev. E* **87**, 052107. `PRIMARY_TEXT_READ`** (full text via harvest.aps.org).
   Revises the long-standing simple-cubic bond threshold 0.248 812 6(5) down to
   0.248 811 82(10). Verbatim: "Figure 1 shows that our estimate of pc lies slightly below the
   central value 0.248 812 6 reported in Ref. [21]."
2. **Mertens (2022), "Exact Percolation Probability on the Square Lattice", *J. Phys. A*
   **55**. `PRIMARY_TEXT_READ`** (ar5iv of arXiv:2109.12102). Weakens every prior Monte Carlo
   estimate of the square-lattice site threshold. Verbatim (abstract): "We use the data to
   compute estimates for the percolation threshold pc that are several orders of magnitude more
   precise than estimates based on Monte-Carlo simulations." Verbatim (text): "The first value
   deviates from the reference value [19] by 2 errorbars, the second by 2.5 errorbars." And:
   "as Bob Ziff put it, 'Monte-Carlo is dead' [40] when it comes to computing pc for 2d
   systems."
3. **Jensen & Ziff (2006), *Phys. Rev. E* **74**, 020101(R). `ABSTRACT_ONLY`.** Fixes a
   previously scattered verdict on a percolation amplitude ratio. Verbatim (abstract):
   "Putting all methods together we find a consistent value Γ−/Γ+=162.5±2, a significant
   improvement over previous results that placed the value of this ratio variously in the range
   of 14 to 220."

So the draft's Q5 worry — that PRE would not take a methods paper whose worked example cannot
name a surviving law — is answerable from these: PRE has published reanalyses whose whole point
was that a previously accepted number or verdict was wrong, including one (Jensen–Ziff) whose
result is precisely the fixing of a scattered amplitude-ratio verdict. The draft's undetermined
eighth verdict, honestly reported as undetermined, is within that precedent.

## 6. Stress-test of the prior Grok assessment

1. **"PRE is the right first venue; the amplitude-ladder literature lives there" — AGREE.**
   PRE's scope sentence (§5.1) plus the fetched PRE record: Jensen–Ziff 2006 (amplitude ratio),
   Mertens–Jensen–Ziff 2017 (universal ratio of cluster-number amplitudes), Wang *et al.* 2013
   (threshold reanalysis). J. Phys. A is a defensible second; JSTAT third.
2. **"The math is not new; novelty is diagnosing the design pathology" — AGREE, with one
   refinement.** The math is not new, but the retrieval shows the novelty claim must be narrowed
   and cited, not just asserted: Parisen Toldin (2011) already uses the covariance scalar
   product and an orthogonal decomposition inside FSS (for error optimization), and von
   Luxburg–Franz (2009) already generalize Fieller geometrically (in statistics). The novel
   content is the *test statistic with a χ² tail on rank(S) − dim(V) d.o.f.*, the diagnosis of
   the self-sealing denominator choice, and the worked inversion — not the geometry.
3. **"Without 3–4 primary ratio-testing citations the literature-practice claim dies" —
   PARTIALLY CONFIRM, and go further.** Three primaries exist (§2.1), so the practice is not
   undocumented; but the primary reading shows the *dominant* practice is full-curve comparison,
   so "ubiquitous" dies even with the citations in hand. Weaken the sentence (§6 of the
   manuscript changes below) rather than propping it up with the three papers.
4. **"Statistics journals will call this a note" — AGREE.** TAS's own scope ("practical
   problems at the interface…") plus the absence of any simulation study makes a statistics
   venue a stretch; Hirschberg & Lye (2010) is the analogue shape there.
5. **"Do not claim a percolation-threshold or modular-weight result" — AGREE, consistent with
   the draft's own §5.3.** Nothing in this retrieval touches that boundary.

## 7. What this retrieval supports changing in the draft

**The one sentence that must change** is §8's

> "The practice of testing amplitude ratios against predicted ratios is ubiquitous — it is how
> conformal-invariance predictions for aspect-ratio and boundary-condition dependence are compared
> with lattice data [6]"

Proposed replacement, consistent with what was actually read:

> "Conformal-invariance predictions for aspect-ratio and boundary-condition dependence are
> routinely compared with lattice data, usually by full-curve comparison of the measured
> probability against the predicted one (Langlands, Pouliot & Saint-Aubin 1994; Watts 1996);
> explicit tests of a ratio of measured amplitudes against a predicted ratio exist (Simmons,
> Kleban & Ziff 2007; Izmailian & Yeh 2009; Kamieniarz & Blöte 1993) but are the exception rather
> than the routine form of the comparison."

Proposed bibliographic replacement for [6] (quotes in §2):

- [6a] Cardy, J. L. (1992), *J. Phys. A* **25**, L201 — the theory source the sentence wants. `[LIT]`
- [6b] Langlands, Pouliot & Saint-Aubin (1994), *Bull. Amer. Math. Soc.* **30**, 1–61. `PRIMARY_TEXT_READ`
- [6c] Watts, G. M. T. (1996), *J. Phys. A* **29**, L363. `PRIMARY_TEXT_READ`
- [6d] Simmons, Kleban & Ziff (2007), *J. Phys. A* **40**, F771. `PRIMARY_TEXT_READ`
- [6e] Izmailian & Yeh (2009), *Nucl. Phys. B* **814**, 573. `PRIMARY_TEXT_READ`
- [6f] Kamieniarz & Blöte (1993), *J. Phys. A* **26**, 201. `ABSTRACT_ONLY`

Additional citation debts found (for the authors' consideration, not patchable here):

- **Parisen Toldin (2011)**, PRE **84**, 025703(R), in §1.4/§8 — nearest prior use of the
  covariance-metric geometry in FSS (error optimization, not testing). `PRIMARY_TEXT_READ`
- **von Luxburg & Franz (2009)**, *Statistica Sinica* **19**, 1095, in §2.4 — the statistics-side
  geometric generalization of Fieller. `ABSTRACT_ONLY`

Per the ticket's tripwire: this note does **not** repeat the "ubiquitous" claim, and the "we
correct a literature practice" framing should be softened to "we recommend a statistic whose
two-entry case is the known-correct special case of the ratio practice" unless the full texts of
Pinson (1994) and Langlands *et al.* (1992) turn up further ratio-testing on reading.

## 8. Sources attempted and blocked

- Springer CDN (Pinson 1994 PDF/HTML; JSP aims-and-scope): HTTP 406/403 IP-block and bot
  challenge, 2026-09-09.
- University of Michigan DeepBlue (Ziff 1995 full text): Cloudflare interstitial, twice.
- ADS (Pinson abstract): human-verification interstitial.
- journals.aps.org direct article pages: blocked/404; PRE full text obtained via
  harvest.aps.org.
- IOPscience /about pages for JSTAT and J. Phys. A: direct fetch failed; scope quotes taken
  from IOP's official publishingsupport.iopscience.iop.org pages instead.
- APS (a guessed Kamieniarz–Blöte DOI page): 404; the IOP abstract page was used instead.
- Reachable and used: arXiv (abs/HTML), ar5iv, IOPscience abstract pages, IAS publications,
  KCL Pure, harvest.aps.org, Inspire HEP API, OpenAlex API, Semantic Scholar API, NIST
  publication pages, pyfssa docs, publisher aims-and-scope pages for Elsevier and TAS.
