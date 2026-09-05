# Q4 — Is the degree/height regularity a theorem about the solvable mechanisms?

**Blocks:** the interpretation of this project's main exclusion result, and the
scope paragraph of the P2 manuscript.
**Bears on:** the original target of the whole project — the square-site
percolation threshold.

---

## The observation

Every planar percolation threshold that is known exactly is an algebraic number of
strikingly low complexity. We have certified the minimal polynomials by exact Sturm
isolation at 100–120 bits:

| lattice | minimal polynomial | degree | height |
|---|---|---:|---:|
| square bond, triangular site | `2x - 1` | 1 | 2 |
| triangular bond | `x³ - 3x + 1` (as `2 sin(π/18)`) | 3 | 3 |
| honeycomb bond, kagome site | `x³ - 3x² + 1` (as `1 - 2 sin(π/18)`) | 3 | 3 |
| (3,12²) site | `x⁶ - 3x⁴ + 1` | 6 | 3 |
| martini bond | `2x² - 1` | 2 | 2 |
| martini descendant | `x² + x - 1` | 2 | 1 |
| Ziff 2006 "A lattice", bond | `p⁵ - 4p⁴ + 3p³ + 2p² - 1` | 5 | **4** |

Each of these comes from one of exactly three mechanisms: **self-duality**,
**self-matching**, or a **star-triangle / Yang–Baxter reduction** (including the
generalized Scullard–Ziff self-dual-cell constructions).

We have used that regularity as a *heuristic* to choose what to exclude: our
exhaustive census covers degree ≤ 6 at height ≤ 4, and separately degree ≤ 4 at
height ≤ 100, and finds no root in the narrowest published interval for the
square-site threshold.

**A datum that bears directly on the question, and that we found the hard way.**
The last row is a late addition. For some time we worked with "degree ≤ 6 **and
height ≤ 3**", which held across the first six rows, and we chose a census class
from it. It is false. The A-lattice threshold is an output of the generalized
self-dual-cell mechanism like the others, and it sits a full unit of height above
every other row. Six mechanism outputs supported a bound that the seventh broke.

That is what an accident of a small sample looks like. It is not proof that no
bound exists — a correct bound might simply be looser than the one we guessed —
but it does mean the regularity should not be assumed, and it is the reason this
question is worth asking rather than answering by inspection. (Our A-lattice row
is corroborated from two independent indexes rather than read from the primary
text; if you know it to be wrong, saying so is itself a useful answer.)

## The question

> **Is the degree-and-height regularity a consequence of the mechanisms, rather
> than an accident of the small sample?**
>
> Concretely: for a threshold produced by a self-duality, self-matching, or
> star-triangle reduction on a planar lattice with a fundamental cell of `n` bonds
> or sites, is there a bound on the degree and on the coefficient height of the
> resulting minimal polynomial, in terms of `n` and the cell's combinatorics?
>
> And then: **does the square site lattice provably fail to admit any threshold in
> the reachable set of those mechanisms?**

The second half is the one that matters. The square-site matching lattice is the
square lattice with both diagonals, which is *not* the square lattice, so
self-matching does not apply; the square site problem is not self-dual; and no
star-triangle reduction is known for it. Those are three separate negative facts.
We want to know whether they are three instances of one obstruction, and whether
that obstruction is provable.

## Why the answer decides something, in both directions

| Answer | What changes |
|---|---|
| **The mechanisms do bound degree and height, and square-site is outside the reachable set** | Our exclusion stops being a numerical census of an arbitrary class and becomes the numerical half of a real statement: *the square-site threshold is not reachable by any known exact mechanism, and is not a low-complexity algebraic number either.* That is a different and much stronger paper. |
| **The mechanisms bound nothing; the regularity is small-sample coincidence** | Our choice of class was arbitrary and we should say so in print, rather than presenting "the historical complexity range" as though it meant something. The census stands as a fact; its interpretation shrinks. |
| **The mechanisms bound degree/height and square-site is *inside* the reachable set** | The best outcome: it says where to look, and the exclusion tells us the answer must have degree ≥ 7, or height > 4 at degrees 5–6, or height > 100 at degrees ≤ 4, on the Jacobsen interval. That is a search space, not a mystery. |

There is no answer to this question that leaves our position unchanged, which is why
it is worth an expensive query.

## If you have an exact value

If the work of answering this produces a candidate exact value or closed form for the
square-site threshold, state it as a **minimal polynomial with integer coefficients**,
or as a closed form precise enough to evaluate to 50 digits. We have a filter that
places any such claim against four published intervals and two exhaustive censuses in
one command (`scripts/threshold_claim_intake.py`), and it will tell us within seconds
whether the claim is already refuted by a committed certificate. Please do not round.

Note that the four published intervals are **pairwise disjoint at their own quoted
precisions**, so at least three of them do not contain the threshold. A value that
disagrees with some of them is not thereby wrong; a value that disagrees with all
four needs to say why.

## What we already have, and do not need re-derived

- The exact minimal polynomials in the table above, certified here.
- The census results: no primitive integer polynomial of degree ≤ 4 and height ≤ 100,
  and none of degree ≤ 6 and height ≤ 4, has a root in the Jacobsen 2015 interval
  `[0.59274605079208, 0.59274605079212]`. Both are exhaustive and certified, not
  searches.
- The standard theory of planar duality, matching lattices, the star-triangle
  relation and the Scullard–Ziff self-dual-cell criterion.
- That the square-site problem is not self-dual and not self-matching. We are asking
  whether that is provably fatal, not whether it is true.

## Do not spend output on

- Re-deriving the known thresholds in the table.
- A survey of percolation, or of exact solvability in statistical mechanics.
- Numerical estimation of the threshold. Ours is at 13 digits and the bottleneck is
  not precision.
- Arguing that the threshold is "probably transcendental" without a mechanism. We
  make no transcendence claim and do not want one asserted.

## Provenance of the framing above

`docs/manuscripts/p2-algebraic-exclusion/manuscript.md` §1.1 and §4.2, whose Table 6
is generated from `results/pslq-lattice-native-candidates/latest.json`; the A-lattice
row from `results/ziff-a-lattice-complexity/latest.json`; the census artifacts under
`results/pslq-degree4-*`, `results/pslq-degree6-low-height-*` and
`results/pslq-degree6-height4-*`; and the interval table in
`analysis/pslq_search_contract.json`.
