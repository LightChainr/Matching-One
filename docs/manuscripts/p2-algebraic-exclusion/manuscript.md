# Certified exclusion of low-degree algebraic relations for a critical threshold

**Manuscript draft — P2 of the publication portfolio. Issue [#551](https://github.com/LightChainr/Matching-One/issues/551).**

**Status:** assembly draft. Every numerical statement below is generated from committed result artifacts by
`scripts/p2_manuscript_evidence_table.py`; the tables live in [`tables.md`](tables.md) and are regenerated, not
hand-edited. Sections still needing a literature pass are marked **[LIT]**. See [`README.md`](README.md) for the
section-by-section readiness ledger.

---

## Abstract (draft)

The site-percolation threshold of the square lattice, `p_c^site(Z^2) ≈ 0.5927460508`, has no known closed form,
and informal integer-relation searches are periodically reported for it. We replace heuristic search with a
*complete certified census*. For every primitive integer polynomial of degree 1 through 4 with coefficient height
at most 100 — 158,062,321,920 polynomials per target — we decide exactly, in rational and Sturm arithmetic,
whether it has a root in each of four method-specific intervals constructed from the published estimates of
Jacobsen (2015), Mertens (2022, two estimators) and Yang–Zhou (2024). Degrees 1–3 are excluded on all four
intervals. At degree 4 the two narrowest intervals are excluded, while the two widest retain 1 and 15 surviving
quartics respectively.

We show that the surviving quartics are artifacts of interval width rather than candidate formulas: each of the
16 survivors has a root in exactly one of the four intervals and is certified outside the other three. We also
give a certified measure of how closely the search class can approach the targets at each degree — 1.5·10⁻⁴,
7.0·10⁻⁸, 7.1·10⁻¹⁰, and below 10⁻¹² at degrees 1, 2, 3, 4 — which locates degree 4 as the boundary degree for
this height: the first degree at which the class's approach resolution reaches the width of the published
intervals, and therefore the last degree at which a negative result carries information.

The conclusion is a bounded one. Within the declared finite complexity class, no low-degree algebraic relation
holds. This is not evidence of transcendence, not a claim that `p_c` is non-algebraic, and not a claim that no
exact representation of another type exists.

---

## 1. Introduction

### 1.1 Exact thresholds in two dimensions **[LIT]**

*Needs a literature pass.* The section should survey the tradition of exactly known planar percolation
thresholds — bond percolation on the square, triangular and hexagonal lattices; site percolation on the
triangular lattice; the star–triangle/duality mechanisms that produce them; and the kagome-site root of
`1 - 3p² + p³` — and then state why square-lattice *site* percolation has resisted the same treatment.
The repository's own exact anchor for the matching pair,

```text
p_c^site(Z^2) + p_c^site(NN+NNN) = 1,
```

constrains the pair without determining either member, and should be cited here as the reason a closed form is
natural to look for and has not been found.

Required citations at minimum: Jacobsen 2015 (J. Phys. A 48, 454003), Mertens 2022 (J. Phys. A 55, 334002),
Yang–Zhou 2024 (J. Phys. A 57, 258001), Jacobsen 2024 Reply (J. Phys. A 57, 258002). See `references.bib`.

### 1.2 Why heuristic integer-relation searches are not enough

Integer-relation searches (PSLQ and relatives) on a decimal estimate of a physical constant are widely reported
and weakly informative. Three defects recur:

1. **The target is not a number.** A published threshold estimate is a method-bound quantity with a quoted
   uncertainty, not a definition of the infinite-lattice constant. Selecting one rounded decimal and searching
   against it silently promotes an estimator to a definition.
2. **The search is not complete.** A heuristic scan reports what it finds; it does not certify what it did not
   find, so a null result has no stated scope.
3. **The arithmetic is not certified.** Binary floating point can neither establish that a residual is nonzero
   nor that a root lies outside an interval.

This paper's contribution is to remove all three defects simultaneously for a bounded but explicitly stated
class, and to report the negative result as a theorem about that class.

### 1.3 Scope statement

We prove statements of exactly one form:

> No primitive integer polynomial of degree `d ≤ D` and coefficient height `≤ H` has a real root in the interval
> `I`, where `I` is constructed from one named published estimate at its own quoted precision.

We do **not** claim that `p_c^site(Z^2)` is transcendental, that it is non-algebraic, that it has no exact
representation of another type (a root of a higher-degree or larger-height polynomial, a closed form in
transcendental constants, an infinite product or series), or that the four source intervals can be replaced by a
single pooled interval. Section 8 restates these boundaries.

---

## 2. Canonical provenance of the target intervals

We do not search against "the" value of `p_c`. We search against four intervals, each constructed from one
published estimate at that source's own stated precision, with no confidence-level homogenization and no
pooling. [Table 1](tables.md#table-1--canonical-provenance-of-the-frozen-method-intervals) is the canonical
provenance table; the machine-readable manifest is `data/literature_threshold_sources.json`, and every interval
endpoint is frozen in `analysis/pslq_search_contract.json`.

The estimators are not interchangeable. Mertens's `p_med` and `p_cell` are defined by
`R_{n,n}(p_med) = 1/2` and `R_{n,n}(p_cell) = R_{n-1,n-1}(p_cell)` on the finite `n × n` square; Jacobsen's
estimate comes from an eigenvalue identity between topological sectors of the periodic Temperley–Lieb algebra on
a semi-infinite cylinder; Yang–Zhou's corrected estimate uses exact helical `h(n)` and an extended cylindrical
estimator. These are four different finite-size observables with four different extrapolations.

### 2.1 The four intervals are pairwise disjoint

Assembling the intervals exactly gives a fact that should be stated plainly rather than smoothed over
([Table 2](tables.md#table-2--interval-ordering-and-separation)):

```text
mertens-2022-p-cell  <  mertens-2022-p-med  <  yang-zhou-2024-corrected  <  jacobsen-2015-eigenvalue
```

with adjacent gaps `1.03·10⁻¹⁰`, `5·10⁻¹³` and `2.38·10⁻¹²`. **No two of the four intervals intersect.** At their
own quoted precisions the four published estimates are mutually inconsistent, and at least three of the four
intervals must fail to contain `p_c`.

This is not a criticism of any source — quoted uncertainties from different extrapolation schemes are not
comparable statistical confidence intervals, which is exactly why the manifest refuses to homogenize them. But it
does have a hard methodological consequence for this paper, and for any integer-relation search on this constant:

> There is no defensible pooled interval. A search reported against a single interval is a search against one
> estimator's extrapolation, and its negative result transfers to the others only by explicit re-run.

Every result in Section 6 is therefore reported per interval and never pooled.

### 2.2 Provenance limitations carried forward

Two provenance items remain open and are stated in the manuscript rather than resolved silently:

- the Yang–Zhou `h(n)` and cylindrical `p_c(n)` tables through `n = 24` are not transcribed from a primary
  full-text source; only the corrected estimate in the publisher-verified abstract is used, and the interval is
  marked `primary_abstract_verified_table_pending`;
- the numerical estimate attributed by secondary compilations to the Jacobsen 2024 Reply is deliberately not
  promoted to a search target, and no fifth interval is constructed from it.

---

## 3. Completeness theorem for the search class

### 3.1 The class

Fix a degree `d ≥ 1` and a height bound `H ≥ 1`. Let

```text
C(d, H) = { a = (a_0, ..., a_d) in Z^(d+1) : gcd(a_0, ..., a_d) = 1, a_d >= 1, |a_i| <= H }
```

be the primitive, sign-normalized integer polynomials of degree exactly `d` and height at most `H`. Sign
normalization is harmless: `P` and `-P` have the same roots. Because the minimal polynomial of an algebraic
number is primitive up to sign, every real algebraic number of degree `≤ D` and height `≤ H` is a root of some
member of `C(1, H) ∪ ... ∪ C(D, H)`.

**Proposition 1 (counting).** `|C(d, H)| = Σ_{g=1}^{H} μ(g) ⌊H/g⌋ (2⌊H/g⌋ + 1)^d`.

*Proof.* Tuples with `a_d ∈ [1, H]` and `a_i ∈ [-H, H]` that are divisible by `g` number
`⌊H/g⌋ (2⌊H/g⌋+1)^d`; Möbius inversion over `g | gcd(a)` gives the claim. ∎

At `H = 100` this yields `12,175`, `3,355,121`, `749,507,743` and `157,309,446,881` for `d = 1, 2, 3, 4`, i.e.
`158,062,321,920` polynomials per interval and `632,329,518,400` declared interval comparisons over the four
targets. These counts are committed independently in `results/pslq-look-elsewhere-ledger/latest.json` and are
cross-checked against every census artifact.

### 3.2 The screening bound

The quartic census is a certified meet-in-the-middle screen followed by exact decisions. Let `[l, u] ⊂ (0, 1)`
be a method interval with midpoint `m`, let `S = 10^15`, and set `w_k = round(S · m^k)` for `k = 1, ..., 4`. For
`a ∈ C(4, H)` define the integer

```text
T(a) = S·a_0 + sum_{k=1..4} a_k · w_k.
```

Let `ρ = Σ_{k=1..4} |S·m^k − w_k| ≤ 2` be the total weight rounding error and `D = H·(1+2+3+4) = 1000` the global
derivative bound, valid because `|P_a'(x)| ≤ Σ_k k|a_k| ≤ D` for `|x| ≤ 1`.

**Theorem 2 (screening completeness).** Put `B = ⌈ S·D·(u−l)/2 + H·ρ ⌉`. If `P_a` has a root in `[l, u]`, then
`|T(a)| ≤ B`.

*Proof.* First, `|T(a) − S·P_a(m)| ≤ Σ_k |a_k| · |S·m^k − w_k| ≤ H·ρ`. Second, if `P_a(ξ) = 0` for some
`ξ ∈ [l, u]`, the mean value theorem gives `|P_a(m)| = |P_a(m) − P_a(ξ)| ≤ D·|m − ξ| ≤ D·(u−l)/2`. Combining,
`|T(a)| ≤ S·D·(u−l)/2 + H·ρ ≤ B`. ∎

Theorem 2 is the completeness statement the census rests on: the enumeration may discard any `a` with
`|T(a)| > B` without an exact test, because no such polynomial can have a root in the interval. Every `a` with
`|T(a)| ≤ B` is retained and receives an exact decision. The screen additionally retains a wider *near* set at
`B + 10^9` (that is, `B + 10⁻⁶` at scale `S`), which certifies the globally closest polynomial whenever its
minimum residual on the interval is below `10⁻⁶` — a hypothesis verified at run time for every interval.

### 3.3 The exact decisions

Retained candidates are decided in exact rational arithmetic, never in binary floating point:

- **root decisions** by Sturm sequences on the square-free part, with root isolation to 120 bits;
- **monotonicity** by isolating the roots of `P_a'` on the interval; a candidate with an internal stationary
  point would require the exact algebraic range path, and the census asserts at run time that none of the
  retained quartics has one (`near_candidates_with_stationary_point = 0` on all four intervals);
- **residuals** as exact rationals at the interval endpoints, so a nonzero residual interval is a proof of
  exclusion rather than a numerical impression.

Degrees 1–3 use exact integer-scaled endpoint and derivative-monotonicity arguments over the same class, with
independent Sturm certificates recorded for the closest witness on each interval.

### 3.4 The four levels of completeness, separated

The manuscript should keep these apart, since reviewers conflate them:

| Level | Statement | Where established |
|---|---|---|
| Counting completeness | `C(d, H)` has exactly the stated cardinality | Proposition 1; look-elsewhere ledger |
| Screening completeness | no root-containing polynomial is discarded by the screen | Theorem 2 |
| Exact certification | every retained candidate gets a proof-level root decision | Sturm/rational decisions, §3.3 |
| Implementation verification | the code realizes the above | §7, regression tests, artifact digests |

---

## 4. Motivating the search bounds

### 4.1 What the bounds must not be

A negative result over a class chosen for convenience is uninteresting. The manuscript needs a reason that
`degree ≤ 4, height ≤ 100` is the scientifically right class, and — equally important — a reason to stop there.

### 4.2 The historically proposed forms **[LIT]**

*Needs a literature pass.* Survey the algebraic forms that have actually been proposed or conjectured for planar
percolation thresholds, and record their degrees and heights. The exactly-known planar thresholds are of very low
complexity: `1/2`, `2 sin(π/18)`, and the kagome-site root of `1 - 3p² + p³` (degree 3, height 3). The claim to
support is that any closed form in the same tradition as the known exact thresholds lies inside `C(≤4, 100)`,
so its absence is a substantive negative result about that tradition.

### 4.3 The approach-resolution argument, and why degree 4 is the stopping point

Independently of the literature survey, the census itself supplies a quantitative reason to stop at degree 4 at
this height. For the closest polynomial of each degree on each interval we compute a certified lower bound on the
distance from the interval to that polynomial's nearest root: if `|P| ≥ r` throughout the interval and
`|P'| ≤ D` on `[0,1]`, then every root of `P` in `[0,1]` is at distance at least `r/D` from the interval. Call
this the class's *approach resolution* at that degree
([Table 4](tables.md#table-4--certified-approach-resolution-of-the-search-class)):

```text
degree 1:  1.53 x 10^-4
degree 2:  6.98 x 10^-8
degree 3:  7.10 x 10^-10
degree 4:  0  (roots inside two intervals; floors 1.90 x 10^-12 and 9.24 x 10^-13 on the other two)
```

The four interval widths span `4·10⁻¹⁴` to `1.6·10⁻¹⁰`. Comparing the two scales gives the structure of the whole
result:

- At degrees 1, 2 and 3 the closest polynomial of the class stays *further than one full interval width* from
  every target — the minimum ratio over all twelve degree-and-interval pairs is 4.78. Exclusion at these degrees
  is therefore not a near miss; the class simply cannot reach the resolution of the published estimates.
- At degree 4 the approach resolution first enters the band of the interval widths. Exactly there, survivors
  appear — and they appear on the two *widest* intervals (15 on `p_cell`, width `1.6·10⁻¹⁰`; 1 on `p_med`, width
  `6·10⁻¹²`) while the two narrowest intervals remain excluded, with margins of only 47.5 and 4.62 interval
  widths.

Degree 4 at height 100 is thus the **boundary degree** for this problem: the first degree whose closest member no
longer clears the interval widths, and the last degree at which a negative result is informative. A degree-5
census at the same height would have an approach resolution well below even the narrowest interval width, so it
would produce survivors on every interval *by counting alone*, whatever the true nature of `p_c`. Extending the
search would enlarge the negative counts and reduce their evidential content.

This is the manuscript's own stop rule, and it is quantitative rather than budgetary.

> **Boundary of this argument.** The approach resolution is a certified statement about the committed closest
> polynomials, not an equidistribution theorem for roots of bounded-height polynomials. It gives no null
> distribution and no p-value, and the survivor counts in §6.2 are reported as counts, not as significance.

---

## 5. Method

1. **Freeze the contract before searching.** `analysis/pslq_search_contract.json` fixes the intervals, the search
   class, the arithmetic requirements, the false-positive controls, and the result policy, and records the SHA-256
   of the provenance manifest. It contains no results. Every census artifact re-verifies the provenance digest at
   run time and fails on drift.
2. **Construct intervals from quoted precision only.** `central_value ± quoted_uncertainty`, no confidence-level
   homogenization, no pooling, each source searched separately.
3. **Enumerate the class.** Sign-normalized primitive tuples, degree exactly `d`, height `≤ 100`.
4. **Screen (degree 4).** Certified fixed-point meet-in-the-middle per Theorem 2, implemented in C++
   (`scripts/degree4_fixed_point_screen.cpp`) and driven from Python, with candidate and root-filter counts
   cross-checked between the two layers.
5. **Decide exactly.** Rational endpoint evaluation, derivative monotonicity, and Sturm isolation at 120 bits for
   every retained candidate; exclusion only when the certified residual interval excludes zero.
6. **Emit an immutable artifact.** One JSON per degree and interval, carrying the search parameters, the exact
   counts, the closest witness with its exact residual, and every root witness with its isolating bracket.

---

## 6. Results

### 6.1 The exclusion table

[Table 3](tables.md#table-3--exclusion-results-by-degree-and-method-interval) is the complete result. Summarized:

| Degree | Class size | Excluded on |
|---:|---:|---|
| 1 | 12,175 | all four intervals |
| 2 | 3,355,121 | all four intervals |
| 3 | 749,507,743 | all four intervals |
| 4 | 157,309,446,881 | `jacobsen-2015-eigenvalue`, `yang-zhou-2024-corrected` |

**Result A.** No primitive integer polynomial of degree 1, 2 or 3 with height at most 100 has a root in any of
the four method intervals.

**Result B.** No primitive integer quartic of height at most 100 has a root in the Jacobsen 2015 interval or in
the Yang–Zhou 2024 corrected interval.

**Result C.** Exactly 1 primitive quartic of height at most 100 has a root in the Mertens `p_med` interval, and
exactly 15 have a root in the Mertens `p_cell` interval.

Results A and B are exclusions. Result C is a census, not a set of candidate formulas — §6.2 is the reason.

### 6.2 The surviving quartics are width artifacts, not candidates

The 16 survivors of Result C are listed with certified root brackets in
[Table 5](tables.md#table-5--cross-interval-resolution-of-the-16-surviving-quartics). Re-deciding each survivor
against all four intervals in exact Sturm arithmetic gives:

**Result D.** Each of the 16 surviving quartics has a root in exactly one of the four method intervals, and is
certified to have no root in the other three. No quartic survives two intervals.

The separations are not marginal. The single `p_med` survivor, `58x⁴ + 99x³ − 7x² + 99x − 84`, has its root at
`0.5927460507870630607956...`, which is `2.4·10⁻¹²` below the Yang–Zhou interval and further still from
Jacobsen's. The 15 `p_cell` survivors sit `1.2·10⁻¹⁰` to `2.6·10⁻¹⁰` below every other interval.

Two consequences follow, and the manuscript should state both.

First, **no survivor can be promoted**. Promoting one would require asserting that its own interval is the
correct one, which by §2.1 contradicts at least two of the other three published estimates.

Second, the survivor counts track interval width and not physics
(`survivor_count_is_monotone_in_width = true`): 15, 1, 0, 0 against widths `1.6·10⁻¹⁰`, `6·10⁻¹²`, `2·10⁻¹³`,
`4·10⁻¹⁴`. The implied densities, `9.4·10¹⁰` and `1.7·10¹¹` survivors per unit length, agree to within a factor
of two across a 27-fold change in width — exactly what a width artifact looks like, and not what a genuine
algebraic relation would look like.

A pointed illustration: the *same* quartic `58x⁴ + 99x³ − 7x² + 99x − 84` is the closest height-100 quartic on
three of the four intervals. It is the best algebraic approximant the class contains near this constant. On
`p_med` it counts as a survivor; on Yang–Zhou and Jacobsen it is a certified miss. Its status is decided entirely
by which published extrapolation one adopts.

### 6.3 The exact-relation controls

- **Positive control.** The kagome-site threshold `1 − 2 sin(π/18)` is recovered by the same pipeline as the
  unique physical root of `1 − 3p² + p³` (degree 3, height 3) with a 120-bit isolating interval, confirming that
  the machinery finds a genuine low-height relation when one exists
  (`results/pslq-kagome-exact-control/latest.json`).
- **Lattice-native candidates.** Four frozen lattice-native values — kagome-site, the (3,12) site value, and two
  martini-descendant constants — are certified disjoint from all four method intervals, with exact separation
  lower bounds (`results/pslq-lattice-native-candidates/latest.json`).
- **Standard-constant relations.** The frozen pairwise bases `a + b·p + c·C` for
  `C ∈ {π, e, log 2, √2, √3, √5}` at height 100 — 20,057,676 relations per interval — contain no relation on any
  interval, with the closest witnesses re-verified at the declared confirmation precision and at all three
  interval stability points (`results/pslq-standard-constant-pairwise/latest.json`,
  `results/pslq-standard-constant-stability/latest.json`).
- **Synthetic false-positive calibration.** A seeded set of 100 random decimals on `[0.55, 0.65]` was run through
  the degree-1 filter; none produced a relation within the declared resolution floor
  (`results/pslq-synthetic-false-positive-calibration/latest.json`).
- **Look-elsewhere ledger.** The exact cardinalities of every frozen search family are committed separately, so
  the size of the hypothesis space is a recorded number rather than an afterthought
  (`results/pslq-look-elsewhere-ledger/latest.json`).

---

## 7. Calibration: what to add, and what not to

Per the portfolio's stop rule, calibration is extended only where it answers a likely reviewer objection.

**Worth adding:**

- an independent second implementation of the final exact filter (a different Sturm/root-isolation code path, or
  a computer-algebra cross-check on the retained candidate sets), since the census's credibility rests on the
  exact decision layer;
- deterministic synthetic constants constructed to *have* a known degree-4 relation at comparable height, run
  through the full quartic pipeline as a positive control at the boundary degree — the existing synthetic
  calibration is degree 1 only;
- interval-perturbation sensitivity for the degree-4 near hits: how the survivor sets move as the interval
  endpoints are varied within the sources' quoted precision.

**Not worth adding, and explicitly out of scope:**

- expanding the constant library to generate more negative counts;
- raising the degree or the height by default (see §4.3 — this reduces evidential content);
- constructing a fifth, pooled or "consensus" interval.

---

## 8. Discussion and scope

The result is a bounded certified exclusion. Restated as the manuscript must leave it:

```text
excluded within the declared finite complexity class
  != transcendental
  != non-algebraic
  != no exact representation of another type
```

Specifically, this work does **not** establish that `p_c^site(Z^2)` is transcendental or irrational; that it is
not algebraic of degree ≤ 4 with height > 100; that it is not algebraic of degree > 4; that it has no closed form
involving transcendental constants beyond the six in the frozen library; or that any one of the four method
intervals contains it.

What it does establish is stronger than what an informal PSLQ report establishes, in a way that is worth stating
directly. A heuristic search reports what it found. This census reports, with proof, what cannot be there: a
complete class, exactly enumerated, exactly decided, against named targets at their own stated precision, with
the negative result's scope fixed in advance by a frozen contract.

Three further contributions are methodological and transfer beyond this constant:

1. **Targets must be plural.** The four published estimates are pairwise disjoint at their quoted precisions
   (§2.1). Any exact-search claim about this constant that reports a single interval has silently chosen an
   estimator.
2. **Survivors must be cross-resolved.** A polynomial surviving one interval and excluded by three others is a
   width artifact (§6.2). Cross-interval resolution is a cheap, decisive test that heuristic searches do not
   perform, and it converts an awkward positive count into a clean negative statement.
3. **Search classes have a boundary degree.** The approach resolution of a bounded class can be certified and
   compared to the target's width (§4.3). Beyond the degree where these meet, negative results are guaranteed to
   fail and positive ones are guaranteed to appear, independently of the physics. This gives an objective place
   to stop, and it applies to any integer-relation search against a measured constant.

### 8.1 Future work that would change the conclusion

A rigorous narrowing of the threshold interval by a *proved* bound rather than an extrapolation would sharpen
every statement here, since exclusion strength is governed by interval width. That is the appropriate next step,
and it is a probability-theory problem, not a search problem.

---

## 9. Reproducibility supplement

The full machine-readable specification, artifact list, and SHA-256 digests are in
`results/p2-algebraic-exclusion-manuscript/latest.json`.

**Frozen inputs**

| Artifact | Role |
|---|---|
| `data/literature_threshold_sources.json` | canonical provenance manifest (digest pinned in the contract) |
| `data/jacobsen_2015_square_site_cylinder.csv` | Jacobsen Table 2, `n = 1..21`, decimals as printed |
| `data/mertens_2022_square_site_estimators.csv` | Mertens Tables 4 and 5, `n = 1..24` |
| `analysis/pslq_search_contract.json` | frozen intervals, class, arithmetic and controls; contains no results |

**Census code**

| Path | Role |
|---|---|
| `scripts/degree1_rational_exclusion.py` | degree-1 exact exclusion |
| `scripts/degree2_polynomial_exclusion.py` | degree-2 exact exclusion |
| `scripts/degree3_interval_exclusion.py` | degree-3 exact exclusion, per interval |
| `scripts/degree4_interval_exclusion.py` | degree-4 census driver, per interval |
| `scripts/degree4_fixed_point_screen.cpp` | certified fixed-point meet-in-the-middle screen |
| `scripts/exact_polynomial_root_certificate.py` | Sturm sequences, root isolation, stationary classification |
| `scripts/p2_manuscript_evidence_table.py` | manuscript assembly; renders `tables.md`, no census computation |

**Results** — ten census artifacts under `results/pslq-degree{1,2,3,4}-*/latest.json` and six control artifacts
under `results/pslq-*/latest.json`, all digested in the manuscript artifact.

**Regeneration**

```bash
python3 scripts/p2_manuscript_evidence_table.py --output results/p2-algebraic-exclusion-manuscript/latest.json
python3 scripts/p2_manuscript_evidence_table.py --markdown --output docs/manuscripts/p2-algebraic-exclusion/tables.md
python3 -m unittest tests.test_p2_manuscript_evidence_table
```

The census itself is not re-run by the manuscript pipeline. To reproduce a census artifact from scratch, run the
corresponding degree script with `--output` and compare digests; the degree-4 intervals are the expensive ones.

---

## Target venue **[LIT]**

*Needs a literature pass, and cross-checking with whoever owns literature search.* The natural candidates are the
journal that published all four source papers (J. Phys. A), for direct comparability, or an
applied/computational number theory venue. The trade-off to resolve: the physics venue gives the result its
audience, while the number-theory venue is the better home for Theorem 2 and the boundary-degree argument.
