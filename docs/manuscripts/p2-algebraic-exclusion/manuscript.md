# Certified exclusion of low-degree algebraic relations for a critical threshold

**Manuscript draft — P2 of the publication portfolio. Issue [#551](https://github.com/LightChainr/Matching-One/issues/551).**

**Status:** complete draft. Every numerical statement below is generated from committed result artifacts by
`scripts/p2_manuscript_evidence_table.py`; the tables live in [`tables.md`](tables.md) and are regenerated, not
hand-edited. Remaining literature items are marked **[LIT]**. See [`README.md`](README.md) for the
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

Three further results make the census interpretable. First, the four published intervals are **pairwise
disjoint** at their own quoted precisions, so at least three of them do not contain `p_c` and no pooled interval
is defensible. Second, each of the 16 surviving quartics has a root in **exactly one** of the four intervals and
is certified outside the other three, which identifies the survivors as artifacts of interval width rather than
candidate formulas. Third, a planted-root control shows the census path detects a height-100 quartic root that is
present, and rejects one that is absent, at **both** widths where the census reported a null — so those
zero-survivor results are sensitivity-certified rather than blind spots.

A companion exhaustion closes the one gap the height-100 class leaves. The exactly-known planar percolation
thresholds are algebraic of degree at most 6 and coefficient height at most 4: the `(3,12²)` site value
`√(1 − 2 sin(π/18))`, a root of `x⁶ − 3x⁴ + 1`, supplies the largest degree, and the Ziff "A lattice" bond
threshold, the root in `[0,1]` of the irreducible quintic `p⁵ − 4p⁴ + 3p³ + 2p² − 1`, supplies the largest
height. Both fall outside `C(≤4, ≤100)`. Exhausting all 2,351,328 primitive integer polynomials of degree ≤ 6
and height ≤ 4 excludes every one of them on all four intervals, and the class's certified approach resolution
still clears every interval width by a factor of at least 54, so the null carries information. No form at the
complexity of any known exact planar threshold has a root in any published interval.

We also give a certified measure of how closely the search class can approach the targets at each degree —
1.5·10⁻⁴, 7.0·10⁻⁸, 7.1·10⁻¹⁰, and below 10⁻¹² at degrees 1 to 4 — which locates degree 4 as the *boundary
degree* for this height: the first degree at which the class's approach resolution reaches the width of the
published intervals, and therefore the last degree at which a negative result carries information.

The conclusion is a bounded one. Within the declared finite complexity class, no low-degree algebraic relation
holds. This is not evidence of transcendence, not a claim that `p_c` is non-algebraic, and not a claim that no
exact representation of another type exists. We also note explicitly that the class is a choice rather than an
inherited bound: the three mechanisms that produce exact planar thresholds bound neither the degree nor the
height of what they produce — published bow-tie thresholds reach degree 11 at height 36 — and square-site is a
lattice where all three mechanisms fail.

---

## 1. Introduction

### 1.1 Exact thresholds in two dimensions

Two-dimensional percolation is unusual among lattice models in that a handful of critical thresholds are known
exactly. Sykes and Essam established the matching relation

```text
p_c^site(L) + p_c^site(L*) = 1
```

for a matching pair of planar lattices, together with the exact bond thresholds `2 sin(π/18)`, `1/2` and
`1 − 2 sin(π/18)` for the triangular, square and honeycomb lattices [Sykes–Essam 1964]. The mechanism behind
these values is duality combined with the star–triangle transformation, and later work extended the same
machinery: site-to-bond and star–triangle constructions give the exact martini-lattice thresholds
[Scullard 2006] **[LIT: primary text not verified]**, and a generalized cell/dual-cell transformation generates
in principle an unbounded family of exactly solvable lattices [Ziff 2006] **[LIT: primary text not verified]**.
That last reference carries more weight here than the others: §4.2 takes one of its exact thresholds — the "A
lattice" bond value — as the polynomial that sets the height bound of the class exhausted in §6.6.
Site thresholds for the Archimedean lattices are catalogued in [Suding–Ziff 1999] **[LIT: primary text not
verified]**, where most entries are numerical and only a few are exact.

Square-lattice *site* percolation is the conspicuous gap. It is the most-studied planar percolation problem, it
has been computed to thirteen digits by several independent methods, and it is not known in closed form. The
matching identity constrains it — the square site lattice and its NN+NNN matching partner satisfy the Sykes–Essam
relation — but the relation determines only the *sum* of the pair, not either member, so it does not by itself
produce a formula. That is precisely why a closed form is natural to look for, and why heuristic searches for one
keep being reported.

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

This paper removes all three defects simultaneously for a bounded but explicitly stated class, and reports the
negative result as a theorem about that class.

### 1.3 Contributions

1. A **complete certified census** of primitive integer polynomials of degree ≤ 4 and height ≤ 100 against four
   independently sourced threshold intervals, with a screening-completeness theorem (§3) and exact Sturm
   decisions throughout. Degrees 1–3 are excluded on all four intervals; degree 4 is excluded on the two
   narrowest (§6.1). A second exhaustion over degree ≤ 6 at height ≤ 4 is excluded everywhere (§6.6). That
   second class contains every entry of our exact-threshold table together with the A-lattice quintic, but it
   is a class we chose and defend by the resolution criterion of §4.3, not a bound the literature supplies:
   §4.2 gives published thresholds at heights 6 to 36, one of them at degree 11. The height was itself
   corrected twice — the census ran at height ≤ 3 first (§6.5) — and we report that sequence rather than the
   final choice alone.
2. The observation that the four published intervals are **pairwise disjoint** (§2.1), which makes per-interval
   reporting a necessity rather than a conservatism, and which is what makes the degree-4 survivors interpretable.
3. A **cross-interval resolution** of every surviving quartic showing each survives exactly one interval (§6.2),
   converting an awkward positive count into a clean negative statement.
4. A **sensitivity control** at the widths where the census reports a null, establishing that those nulls are
   not blind spots (§6.4).
5. The **boundary-degree** criterion (§4.3): a certified comparison between a search class's approach resolution
   and the target's width, which gives an objective place to stop that applies to any integer-relation search
   against a measured constant.

Contributions 2, 3, 4 and 5 are methodological and transfer beyond this constant. In our view they are the more
durable part of the paper; the census is the evidence that forced them.

### 1.4 Scope statement

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
| Implementation verification | the code realizes the above | §6.4 sensitivity control; regression tests; artifact digests |

---

## 4. Motivating the search bounds

### 4.1 What the bounds must not be

A negative result over a class chosen for convenience is uninteresting. The manuscript needs a reason that
`degree ≤ 4, height ≤ 100` is the scientifically right class, and — equally important — a reason to stop there.

### 4.2 The complexity of the known exact thresholds

Every planar percolation threshold that is known exactly is an algebraic number of very low complexity. Taking
the minimal polynomials from the repository's certified lattice-native candidate artifact
([Table 6](tables.md#table-6--algebraic-complexity-of-the-exactly-known-planar-thresholds)):

| Lattice | Closed form | Minimal polynomial | Degree | Height |
|---|---|---|---:|---:|
| square bond; triangular site | `1/2` | `2x − 1` | 1 | 2 |
| triangular bond | `2 sin(π/18)` | `x³ − 3x + 1` | 3 | 3 |
| honeycomb bond; kagome site | `1 − 2 sin(π/18)` | `x³ − 3x² + 1` | 3 | 3 |
| (3,12²) site | `√(1 − 2 sin(π/18))` | `x⁶ − 3x⁴ + 1` | 6 | 3 |
| martini bond | `1/√2` | `2x² − 1` | 2 | 2 |
| martini descendant | `(√5 − 1)/2` | `x² + x − 1` | 2 | 1 |

The triangular-bond row is derived from the certified kagome-site row by the matching substitution `p → 1 − p`,
the same Sykes–Essam relation that anchors this problem — a small internal consistency check on the table.

One row belongs in this table and is not generated with it. The generalized cell/dual-cell construction of
[Ziff 2006] yields an exact bond threshold for the lattice it calls **A**, namely the root in `[0,1]` of

```text
p⁵ − 4p⁴ + 3p³ + 2p² − 1
```

which we certify here to be irreducible over `Q`, of degree 5 and **coefficient height 4**, with its unique root
in `(0,1)` isolated at 120 bits and agreeing with the separately quoted decimal `0.625457`
(`results/ziff-a-lattice-complexity/latest.json`). It is not in the repository's frozen lattice-native candidate
library, which is why it is absent from the generated Table 6; the library is pinned by a contract digest and we
have not reopened it for a value that is not a square-site candidate.

That row refutes the height-3 bound we had assumed. It is important to say what it does **not** do, because we
made this error once and the natural repair is to make it again one unit higher. It does not install height 4 as
the record's ceiling. The bow-tie bond threshold of [Wierman 1984], reproduced in [Ziff–Scullard 2006], is the
root in `[0,1]` of

```text
1 − p − 6p² + 6p³ − p⁵
```

with `p_c = 0.404518…`: degree 5 and **coefficient height 6**, twenty-two years before the A-lattice. The
generalized bow-ties of the same 2006 paper reach degree 11 at heights 35 and 36, and the asymmetric bow-tie of
[Ziff–Scullard–Wierman–Sedlock 2012] is degree 8 at height 15. The three mechanisms of §4.2.1 bound neither the
degree nor the height of what they produce; the cell can be made larger, and the polynomial grows with it. There
is no published uniform bound, and we are aware of no argument that one exists.

The consequence for this paper is a change of standing, not of result. **A search class is a choice we make and
must defend, not a range the literature hands us.** Every census below remains an exhaustive certified null on
the class it names. None of them is a null on "the complexity of exactly-known thresholds", because that phrase
does not denote a bounded set.

We keep the sequence visible rather than presenting the final class as if it had been chosen first: we assumed
height ≤ 3 from a table we assembled ourselves, found the A-lattice quintic at height 4 and widened, and then
found — from the primary texts, on the second look — that the record had reached height 6 in 1984 and height 36
by 2006. Twice the error was the same one: reading a bound off an incomplete table. What survives that is the
exhaustion, not the framing.

The degree bound is the one place where the height-100 class does not dominate the tradition. The `(3,12²)` site
value has degree 6, because it is the square root of a lower-degree threshold. A "`(3,12²)`-like" closed form for
square-site percolation — a square root, or another low-degree radical, of a simple algebraic number — would have
degree 5 or 6 and is **not** covered by `C(≤4, ≤100)`.

Rather than raise the height without limit, which §4.3 shows would destroy the negative result's content, we
close that gap at **degree ≤ 6 at height ≤ 4**. That is the largest class in this family whose certified approach
resolution still clears every published interval width — the criterion of §4.3, and the only defensible reason to
stop anywhere — and it contains every entry of Table 6 together with the A-lattice quintic. It does not contain
the bow-tie thresholds above, and we do not claim it does. That class holds 2,351,328 polynomials per
interval, `1.5·10⁻⁵` of the class already searched, and is reported in §6.6. Its approach resolution is checked
there and still clears every interval width, so this particular widening does not cost the null its content.

### 4.2.1 The historical range is a selected sample, and we say so

Table 6 is not a random sample of planar percolation thresholds. It is the list of
thresholds that are **known exactly**, and a planar threshold becomes known exactly
when one of three mechanisms delivers it: self-duality, self-matching, or a
star-triangle (Yang–Baxter) reduction, including the generalized self-dual-cell
constructions. Every row of Table 6 is an output of one of those three.

Square-site percolation is precisely a case where all three fail. Its matching lattice
is the square lattice with both diagonals, which is not the square lattice, so the
Sykes–Essam self-matching argument that fixes the triangular-site value at `1/2` does
not close; the site problem on `Z²` is not self-dual; and no star-triangle reduction
for it is known.

So "degree ≤ 6, height ≤ 4" is the complexity range **reachable by those mechanisms**,
observed on the lattices where they work. Whether the mechanisms bound the degree and
height of their own output — whether the regularity is a theorem in disguise rather
than a coincidence among six numbers — is, as far as we are aware, open. We do not
assume it, and a reader should not read Table 6 as a prior over the algebraic
complexity of an arbitrary threshold.

The A-lattice row makes the point concretely. It is a mechanism output like the others —
a generalized self-dual-cell construction — and it sits a full unit of height above every
row of Table 6. Six numbers were enough to suggest "height ≤ 3"; the seventh broke it. We
have no reason to believe the seventh is the last one.

Two things follow, and they pull in opposite directions.

The negative result of §6.6 should be read for what it is: *no form at the complexity
of a mechanism-produced threshold fits any published interval.* That is a meaningful
statement about the square-site value — it says the constant does not look like the
ones we can derive — and it is weaker than "no simple algebraic number fits", which we
do not claim.

And it is why §4.3 exists. The approach-resolution criterion given there is derived
from the search class itself and makes no reference to the literature at all. It
supplies a stopping degree that does not rest on a selected sample, and it is the
bound we would keep if Table 6 turned out to carry no information.

### 4.3 The approach-resolution argument, and why degree 4 is the stopping point

Independently of the literature, the census itself supplies a quantitative reason to stop at degree 4 at this
height. For the closest polynomial of each degree on each interval we compute a certified lower bound on the
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
search in height would enlarge the negative counts and reduce their evidential content.

This is the manuscript's own stop rule, and it is quantitative rather than budgetary. Note that it constrains
*height* extension specifically; a degree extension at *low* height is a different proposition, and §8.1 treats
it separately.

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
- **Lattice-native candidates.** Four frozen lattice-native values — kagome-site, the (3,12²) site value, and two
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

### 6.4 The nulls are sensitivity-certified, not blind spots

Results A and B are null results at widths as small as `4·10⁻¹⁴`. A referee is entitled to ask whether the
pipeline could have found a quartic root at that width had one been there. We answer it directly
([Table 7](tables.md#table-7--quartic-census-sensitivity-where-the-census-returned-a-null)).

Take a quartic root that is *known* to exist — a root witness the census itself already committed — and plant it
inside a synthetic interval of the width in question, placed off-centre. Then run the **unmodified** census path
on that interval. The negative twin shifts the same interval one full width away from the planted root, so the
root is certified outside.

The control covers exactly the intervals whose census result was a null, and it reads that set from the census
artifacts rather than naming it. On `p_med` and `p_cell` the census itself found 1 and 15 roots, so its
sensitivity at those widths is already demonstrated and a planted root would add nothing; the question has force
only where the answer was zero.

**Result E.** Over two planted quartics × the two null-result widths × two polarities, all 8 trials behave as
required: every positive trial reports the planted quartic among its root witnesses, and no negative trial
reports it (`results/pslq-degree4-synthetic-boundary-control/latest.json`).

In every positive trial the interval returns exactly one root — the planted one — and every negative trial
returns none. So at `4·10⁻¹⁴` and `2·10⁻¹³` these regions are genuinely sparse and the pipeline resolves them:
the zero-survivor results of Result B are certified nulls, not blind spots.

This control also serves as the implementation-verification level of §3.4: it exercises the screen, the
monotonicity check and the Sturm decisions end to end, on inputs whose correct answer is known in advance.

### 6.5 The degree ≤ 6, height ≤ 3 class is closed

Every entry of Table 6 has degree ≤ 6, and exactly one — the `(3,12²)` site value, root of `x⁶ − 3x⁴ + 1` —
lies outside `C(≤4, ≤100)`. We exhaust that degree range directly.

This section reports the exhaustion at **height ≤ 3**, which was the bound we believed the published record
respected when the census was designed and run. It does not: §4.2 gives the A-lattice quintic at height 4 and
the bow-tie thresholds at heights 6, 15, 35 and 36, the last of them at degree 11. §6.6 reports the census at
the widened height. We keep this section as it stands, rather than rewriting it, because the height-3 result is
what the committed artifacts, the sensitivity control and the independent replication of §7 actually cover, and
because the sequence — a class chosen from an incomplete reading of the literature, widened when the reading
was corrected, then found still not to bound anything — is the sort of thing a reader is entitled to see rather
than have tidied away.

**Class.** `C(d, 3)` for `d = 1..6`: `15 + 129 + 975 + 7,041 + 49,935 + 351,489 = 409,584` primitive
sign-normalized polynomials per interval.

**Method.** The class is small enough to evaluate every member exactly at both interval endpoints, in integer
arithmetic on the common denominator. Two consequences of `|P'| ≤ D` on `[0,1]`, with `D = 3d(d+1)/2`, do the
rest. If `|P(l)| > D·(u−l)` then `P` cannot vanish on `[l, u]`, since `|P(x)| ≥ |P(l)| − D·(x−l)`; only
polynomials failing that test receive an exact Sturm decision. And for a root `ξ` of `P` outside `[l, u]`, either
`l − ξ ≥ |P(l)|/D` or `ξ − u ≥ |P(u)|/D`, so `min(|P(l)|, |P(u)|)/D` is a certified lower bound on the distance
from the interval to the nearest root — no monotonicity assumption required.

**Result F.** No primitive integer polynomial of degree ≤ 6 and height ≤ 3 has a root in any of the four method
intervals. The exclusion is not marginal: the certified screen retains **zero** candidates at every degree on
every interval, and the closest member of the whole class stays `9.23·10⁻⁸` away — between `5.8·10²` and
`2.3·10⁶` interval widths, depending on the target.

See [Table 8](tables.md#table-8--exhaustion-of-the-historical-complexity-range-degree--6-height--3) for the
per-degree class sizes, closest polynomials, certified distance floors and floor-to-width ratios.

The degree-6 minimiser is `x` times the degree-5 minimiser — at this height the class simply contains no
irreducible sextic that comes closer, which is itself a statement about how coarse the historical complexity
range is near this constant.

**Sensitivity (Result G).** Because the screen retains nothing, the Sturm decision path never executes during the
exclusion, and a null produced that way needs the same treatment as §6.4. We plant the `(3,12²)` polynomial
itself — the exact form that motivated this section, taken from the certified lattice-native artifact — inside a
synthetic interval of each of the four method widths, at both polarities, and run the **unmodified** scan. All 8
trials pass: every positive trial retains exactly one candidate, decides it, and reports `x⁶ − 3x⁴ + 1`; every
negative trial retains the same candidate and correctly returns no root
(`results/pslq-degree6-low-height-control/latest.json`).

### 6.6 The degree ≤ 6, height ≤ 4 class is closed as well

Result F is stated at height ≤ 3, and §4.2 shows nothing in the literature licenses that bound: the A-lattice
quintic `p⁵ − 4p⁴ + 3p³ + 2p² − 1` has height 4. We therefore widen to `C(d, 4)` for `d = 1..6` and exhaust it.
This class is chosen, not inherited — §4.2 gives published thresholds well outside it — and §4.3 is why we stop
here rather than higher: it is the largest class in the family whose certified approach resolution still clears
every published interval width.

**Class.** `23 + 265 + 2,639 + 24,913 + 229,703 + 2,093,785 = 2,351,328` primitive sign-normalized polynomials
per interval — 5.7 times the height-3 class, and still `1.5·10⁻⁵` of the quartic census of §6.1. The derivative
bound becomes `D = 4d(d+1)/2 = 2d(d+1)`; nothing else in the method of §6.5 changes, and the same script runs it
with the height threaded through as a parameter, so the two censuses share one implementation by construction.

**Result F′.** No primitive integer polynomial of degree ≤ 6 and height ≤ 4 has a root in any of the four method
intervals. As at height 3, the certified screen retains **zero** candidates at every degree on every interval,
so every one of the 9,405,312 interval-polynomial pairs is decided by the exact endpoint bound alone. The
closest member of the whole class stays `8.70·10⁻⁹` away from the nearest interval — a margin of **54 interval
widths** at its narrowest, against `5.8·10²` at height 3.

See [Table 10](tables.md#table-10--exhaustion-of-the-corrected-historical-range-degree--6-height--4) for the
per-degree class sizes, closest polynomials, certified distance floors and floor-to-width ratios.

That factor of about eleven in the margin is the price of one extra unit of height, and it is worth naming
plainly rather than burying in a table. One unit of height cost an order of magnitude of approach resolution;
the null is not indefinitely robust to widening the class, and §4.3 is where that limit is made quantitative.
At height ≤ 4 the margin still clears every published width by more than fifty times, so this widening does not
cost the result its content — but a reader should not extrapolate the null to height 10.

**Sensitivity.** The planted-root control of Result G is not re-run at height 4. It plants a height-3 sextic,
and it exercises the scan path — the same screen and the same exact Sturm decision, over a strictly larger
class — rather than the enumeration bound. What it certifies (that a root inside a frozen-width interval is
retained, decided and reported) carries here unchanged. What it does not certify is the height-4 enumeration
itself; that is checked instead against independently counted class sizes, degree by degree, in the assembly
script. We state this asymmetry rather than claim a control we did not run.

**Replication.** The independent second implementation of §7 covers height ≤ 3 only. Result F′ therefore rests
on a single implementation of its enumeration, and the two implementations agreeing at height 3 is evidence
about the shared scan path rather than about this class. Extending the replication is cheap — it is the same
height parameter — and is listed in §7 as recommended rather than done.

Taken with Results A–C, this gives the paper's cleanest statement:

> **No integer polynomial of degree ≤ 6 and coefficient height ≤ 4, and none of degree ≤ 4 and height ≤ 100,
> has a root in any of the four published intervals for `p_c^site(Z²)`.**

An earlier draft of this paper stated that sentence as "no algebraic form at the complexity of any exactly-known
planar percolation threshold", and we withdraw that phrasing. It is not supportable: §4.2 gives published
thresholds at heights 6, 15, 35 and 36, one of them at degree 11, so "the complexity of an exactly-known
threshold" is not a bounded set and no finite census can exhaust it. The statement above says exactly what was
searched, which is both weaker and true. Everything the census actually establishes survives the change of
wording; what does not survive is the claim to have covered the literature.

The A-lattice row itself is no longer the weak link it was. Its verification status is `PRIMARY_TEXT_READ`
(§4.2): the primary text was read and confirms both the quintic and the value. What that row establishes is the
falsity of our original height-3 bound — not a new bound of its own.

## 7. Calibration: what was added, and what was not

Per the portfolio's stop rule, calibration is extended only where it answers a likely reviewer objection.

**Added (§6.4):** a degree-4 planted-root sensitivity and specificity control at the two widths where the census
reports a null. The pre-existing synthetic calibration was degree 1 only, while the paper's interpretation is
load-bearing at degree 4 — the boundary degree — and the strongest claims (Results A and B) are nulls at the
narrowest widths. This was the gap most likely to be challenged.

**Added (§6.5, partial):** an independent second implementation, covering the degree ≤ 6 height ≤ 3 census but
not the quartic census. `scripts/degree6_independent_replication.py` was written separately from the same frozen
protocol and enumerates `C(d,3)` in its own code, screening at the interval *midpoint* where
`scripts/degree6_low_height_exclusion.py` screens at both *endpoints*. Both screens are certified consequences
of `|P'| ≤ D` on `[0,1]`, but they evaluate different points against different bounds
([Table 9](tables.md#table-9--agreement-of-two-independent-implementations-degree--6-height--3)). Every
interval-by-degree cell agrees on class size, screen survivors, root-containing polynomials and distinct roots;
the exclusion verdict agrees on all four intervals; and both implementations single out the same closest member
of the class, coefficient for coefficient.

Equal counts alone would be weak evidence, since two implementations can agree by both being empty. The check
that carries the weight is the residual. The two evaluate that minimiser at *different points*, so their
residuals must differ — and the mean value theorem caps the difference at `D(u-l)/2` for that polynomial's own
`D = Σ k|a_k|`. Table 9 gives the observed gap and the allowance on each interval; the gap is non-zero, tracks
the interval width, and stays inside the allowance everywhere. The two implementations are computing the same
quantity, not merely reporting the same zero.

The replication is partial in a way worth stating precisely. Both import
`scripts/exact_polynomial_root_certificate.py` unchanged, so the Sturm path is shared rather than replicated.
For this census that shared code contributes nothing: both screens retain zero candidates at every degree on
every interval, so root isolation never runs, and the null is produced end to end by two independently written
certified screens. It does run in the planted-root controls of §6.5, where the replication does not reach.

**Recommended, not done:**

- an independent second implementation of the *quartic* census (§6.1–6.3), whose C++ fixed-point
  meet-in-the-middle screen and whose Sturm decisions on the 16 retained candidates are both single-implementation
  — this is where Results A–D live, and it is the more valuable of the two replications;
- extending the existing second implementation from height ≤ 3 to height ≤ 4, so that Result F′ — which now
  carries the paper's closing statement — has the replication that Result F has. This is the same height
  parameter and about 70 s of arithmetic. We did not thread it ourselves: the replication's value is that it was
  written by another party against the frozen protocol and committed unchanged, and an edit by us to make it
  cover our own later class would have to be disclosed as such, which costs more than it buys. The extension is
  worth asking its author for;
- interval-perturbation sensitivity for the degree-4 near hits: how the survivor sets move as the interval
  endpoints are varied within the sources' quoted precision.

**Explicitly out of scope:**

- expanding the constant library to generate more negative counts;
- raising the *height* by default (see §4.3 — this reduces evidential content);
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
not algebraic of degree ≤ 4 with height > 100; that it is not algebraic of degree 5 or 6 at height > 3, nor of
degree > 6 at any height; that it has no closed form involving transcendental constants beyond the six in the
frozen library; or that any one of the four method intervals contains it.

What it does establish is stronger than what an informal PSLQ report establishes, in a way that is worth stating
directly. A heuristic search reports what it found. This census reports, with proof, what cannot be there: a
complete class, exactly enumerated, exactly decided, against named targets at their own stated precision, with
the negative result's scope fixed in advance by a frozen contract, and with the sensitivity of the null
independently certified.

Four contributions are methodological and transfer beyond this constant:

1. **Targets must be plural.** The four published estimates are pairwise disjoint at their quoted precisions
   (§2.1). Any exact-search claim about this constant that reports a single interval has silently chosen an
   estimator.
2. **Survivors must be cross-resolved.** A polynomial surviving one interval and excluded by three others is a
   width artifact (§6.2). Cross-interval resolution is a cheap, decisive test that heuristic searches do not
   perform, and it converts an awkward positive count into a clean negative statement.
3. **Nulls must be sensitivity-certified.** Planting a known root at the width where a null was reported and
   re-running the unmodified pipeline (§6.4) distinguishes "nothing is there" from "the method cannot see it".
   This costs seconds and is, as far as we know, absent from the integer-relation literature.
4. **Search classes have a boundary degree.** The approach resolution of a bounded class can be certified and
   compared to the target's width (§4.3). Beyond the degree where these meet, negative results are guaranteed to
   fail and positive ones are guaranteed to appear, independently of the physics. This gives an objective place
   to stop, and it applies to any integer-relation search against a measured constant.

### 8.1 Future work

**Beyond the class we closed.** §6.6 closes degree ≤ 6 at height ≤ 4, which covers every entry of Table 6 and
the A-lattice quintic. What remains untested starts immediately above it, and includes published thresholds:
the bow-tie bonds of §4.2 sit at heights 6 to 36, one of them at degree 11. Degree 5–6 above height 4, and
degree > 6 at any height, are both open here. Neither is a free extension. By the
boundary-degree criterion of §4.3, raising the height on a fixed degree drives the class's approach resolution
below the interval widths, at which point survivors appear by counting alone and a null carries no information;
degree ≤ 6 at height ≤ 10 is already `890,350,944` polynomials per interval, with degree 6 alone at
`848,419,937`, comparable to the height-100 cubic census. Any such extension should therefore begin by computing
the class's approach resolution and confirming it still clears the widths — a cheap calculation that decides
whether the census is worth running at all.

**Is there a mechanism bound at all?** Earlier drafts treated the low complexity of
Table 6 as a regularity that might turn out to be a theorem. On the published record it
is not one. The bow-tie thresholds of §4.2 are outputs of the same star-triangle
mechanism and reach degree 11 at height 36, and the cell can be enlarged further; the
degree and height grow with it, with no published uniform bound. So the "coincidence"
branch is the one the evidence supports, and its consequence is the one we have taken:
our class is a choice defended by §4.3's resolution criterion, the census stands
unchanged as a fact about that class, and its interpretation is what §4.3 supports on
its own. The question that remains live, and that we regard as the most valuable one
adjacent to this work, is narrower and better posed: **can square-site percolation be
shown to lie outside every cell reachable by these mechanisms, not merely outside the
cells tried so far?** [Ziff 2006] states that his method does not reach it, in those
words, but that is a statement about one method in one year. A negative answer would
give the square-site case a reason rather than an observation. Our own history here is
a caution on how such a question should be approached: we read a bound off an
incomplete table twice, at height 3 and then at height 4, and both times the next paper
we read broke it. A bound on these mechanisms will have to come from an argument about
the cell, not from a longer table.

**A rigorous interval.** Exclusion strength is governed by interval width, and every interval used here comes
from an extrapolation rather than a proof. A rigorous narrowing of the threshold interval by a proved bound would
sharpen every statement in this paper. That is a probability-theory problem, not a search problem, and it is the
step that would matter most.

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

**Census and control code**

| Path | Role |
|---|---|
| `scripts/degree1_rational_exclusion.py` | degree-1 exact exclusion |
| `scripts/degree2_polynomial_exclusion.py` | degree-2 exact exclusion |
| `scripts/degree3_interval_exclusion.py` | degree-3 exact exclusion, per interval |
| `scripts/degree4_interval_exclusion.py` | degree-4 census driver, per interval |
| `scripts/degree4_fixed_point_screen.cpp` | certified fixed-point meet-in-the-middle screen |
| `scripts/degree4_synthetic_boundary_control.py` | planted-root sensitivity control (§6.4) |
| `scripts/degree6_low_height_exclusion.py` | degree-1..6 exhaustion at height 3 (§6.5) and height 4 (§6.6) |
| `scripts/ziff_a_lattice_complexity.py` | certifies the A-lattice quintic's degree, height and root (§4.2) |
| `scripts/degree6_low_height_control.py` | planted `(3,12²)` sensitivity control (§6.5) |
| `scripts/degree6_independent_replication.py` | second, independently written implementation of that census (§7) |
| `scripts/degree6_implementation_agreement.py` | cell-by-cell comparison of the two implementations (§7) |
| `scripts/exact_polynomial_root_certificate.py` | Sturm sequences, root isolation, stationary classification |
| `scripts/p2_manuscript_evidence_table.py` | manuscript assembly; renders `tables.md`, no census computation |

**Results** — ten census artifacts under `results/pslq-degree{1,2,3,4}-*/latest.json`, four historical-range
artifacts under `results/pslq-degree6-low-height-*/latest.json` and four more at the corrected height under
`results/pslq-degree6-height4-*/latest.json`, the A-lattice certification in
`results/ziff-a-lattice-complexity/latest.json`, four replication artifacts under
`results/pslq-degree6-low-height-replication-*/latest.json`, their comparison in
`results/pslq-degree6-implementation-agreement/latest.json`, and the control artifacts under
`results/pslq-*/latest.json` including both sensitivity controls, all digested in the manuscript artifact.

**Regeneration**

```bash
python3 scripts/degree4_synthetic_boundary_control.py \
    --output results/pslq-degree4-synthetic-boundary-control/latest.json
python3 scripts/degree6_low_height_exclusion.py --all
python3 scripts/degree6_low_height_exclusion.py --all --height 4
python3 scripts/ziff_a_lattice_complexity.py \
    --output results/ziff-a-lattice-complexity/latest.json
python3 scripts/degree6_low_height_control.py \
    --output results/pslq-degree6-low-height-control/latest.json
for interval in jacobsen-2015-eigenvalue mertens-2022-p-med \
                mertens-2022-p-cell yang-zhou-2024-corrected; do
    python3 scripts/degree6_independent_replication.py "$interval" \
        --output "results/pslq-degree6-low-height-replication-$interval/latest.json"
done
python3 scripts/degree6_implementation_agreement.py
python3 scripts/p2_manuscript_evidence_table.py \
    --output results/p2-algebraic-exclusion-manuscript/latest.json
python3 scripts/p2_manuscript_evidence_table.py --markdown \
    --output docs/manuscripts/p2-algebraic-exclusion/tables.md
python3 -m unittest tests.test_p2_manuscript_evidence_table tests.test_degree4_synthetic_boundary_control
```

The census itself is not re-run by the manuscript pipeline. To reproduce a census artifact from scratch, run the
corresponding degree script with `--output` and compare digests.

---

## Target venue

**Decision: *Journal of Physics A: Mathematical and Theoretical*, with *Experimental Mathematics* as the
fallback.**

The methodological contributions of §8 are the more durable part of this paper, and they would be at home in a
computational number theory or methods venue. That argues for the fallback, not the primary. The deciding
consideration runs the other way:

The paper's central claim is about **one physical constant**, and its most consequential readers are precisely
the people who produce the estimates it tests. Three of the four intervals come from papers in this journal —
Jacobsen 2015, Mertens 2022, Yang–Zhou 2024 — as does the 2024 Reply. §2.1 shows those four published intervals
are pairwise disjoint at their own quoted precisions, which is a correction to that literature's reporting
practice and not merely an input to our census. Published anywhere else, that correction does not reach the
people it is about. The methodology travels by citation regardless; the correction does not.

Secondary considerations pointing the same way: the venue has an established remit for exact and
transfer-matrix computation in statistical mechanics, so Theorem 2 and the boundary-degree criterion need no
special pleading; and direct comparability with the four source papers is what lets a reader check the interval
table against its origins in one sitting.

Take *Experimental Mathematics* if the physics venue judges the work insufficiently physical — its remit covers
certified computational searches and negative results explicitly, and the boundary-degree criterion is closer to
its centre of gravity than to J. Phys. A's.

**Framing note for submission.** Lead with the census as the work and §8 as what the work forced. A submission
that leads with the methodology invites the reading that the constant is a worked example, which understates
both the completeness of the search and the significance of the disjointness finding for this literature.
