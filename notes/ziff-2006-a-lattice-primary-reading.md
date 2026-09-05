# Ziff 2006 A-lattice: primary reading

**Date:** 2026-09-05
**Ticket:** #566 (P2 blocking)
**Claim level:** C1 on the polynomial identity; theory input, not lattice evidence
**This note does not enter** `docs/STATUS.md`

The P2 manuscript's height bound is load-bearing and was held as
`CORROBORATED_NOT_PRIMARY`. The three `[LIT: primary text not verified]`
sources named in #566 were read from the arXiv versions that the journals
identify with the published papers. The polynomial we hold is the one in
the paper. Nothing in the census needs to be re-run.

The JSON artifact and the `verification_status` flip live on PR #564, which
is not this PR. This note is the reading. Once #564 is on `main`, the
status field should move with this note as the citation, and the test
`test_the_sourcing_status_is_not_quietly_upgraded` should move with it.
Do not copy the artifact onto this branch: that would duplicate #564.

## What was read

| Source | arXiv | Journal | What we needed |
|---|---|---|---|
| Ziff, *Generalized cell/dual-cell transformation for percolation, and new exact thresholds* | [cond-mat/0510245](https://arxiv.org/abs/cond-mat/0510245) | Phys. Rev. E **73**, 016134 (2006) | A-lattice bond polynomial; square-site exclusion |
| Scullard, *Exact site percolation thresholds using the site-to-bond and star-triangle transformations* | [cond-mat/0507392](https://arxiv.org/abs/cond-mat/0507392) | Phys. Rev. E **73**, 016107 (2006) | martini site polynomial |
| Suding–Ziff, *Site percolation thresholds for Archimedean lattices* | [cond-mat/9811416](https://arxiv.org/abs/cond-mat/9811416) | Phys. Rev. E **60**, 275 (1999) | square site is numerical, not exact |

The arXiv HTML of each was read in full, not a search-engine summary of it.

## 1. Ziff 2006 — A-lattice bond threshold

Ziff's A-lattice is the martini cell with the upper bond removed (or occupied
with probability 1). For equal bond occupation probability `p` he equates
the cell and dual-cell three-point connectivities and obtains, verbatim,

```text
p^5 - 4 p^4 + 3 p^3 + 2 p^2 - 1 = 0
```

"whose solution gives `p_c` (bond) = 0.625457 …"

That is exactly the polynomial in `results/ziff-a-lattice-complexity/latest.json`
on #564: ascending `(-1, 0, 2, 3, -4, 1)`, degree 5, coefficient height 4,
unique root in `(0,1)` agreeing with the quoted decimal at six places. The
primary text does **not** disagree with us.

A-lattice **site** percolation is a different number: `p_c = 1/√2`, which
Scullard already had. We do not use that row.

## 2. Square site is named as unsolved, in the same paper

Verbatim, after the new exact thresholds:

> Unfortunately, the method does not appear to work for some of the more
> notorious unsolved systems: site percolation on the square and honeycomb
> lattices, and bond percolation on the kagomé lattice.

This is the primary-text version of the Q4 second half. It is Ziff's own
statement that square site is outside the cell/dual-cell (and, in that
paper, star-triangle) reachable set. It is not a theorem that no later
mechanism can reach it. It is the 2006 author saying the method that
produced the A-lattice quintic does not produce square site.

## 3. Scullard 2006 — martini site

PRE abstract and body: martini site threshold is the real root in `(0,1)` of

```text
p^4 - 3 p^3 + 1 = 0
```

giving `p_c = 0.764826…`. The other two exact site values in that paper are
`(√5 - 1)/2` and `1/√2`. These match the manuscript table. Height of the
martini quartic is 3, so it does not move the height bound.

## 4. Suding–Ziff 1999 — square site is not an exact threshold in this paper

They compute eight Archimedean **site** thresholds numerically, errors
about `± 3×10^{-6}`. They explicitly skip square, triangular and kagomé:

> We did not consider the square, triangular, and Kagomé lattices, as
> `p_c` (site) is either known exactly (triangular and Kagomé) or has
> already been measured to a high degree of precision (square).

Square site is cited as a high-precision numerical value (`0.592746`), not
as an algebraic form. The one exact Archimedean site value they confirm is
`(3,12²)`: `[1 - 2 sin(π/18)]^{1/2}`, consistent with their Monte Carlo.
That is the degree-6 height-3 form already in the census. Nothing here
gives square site a polynomial.

## 5. What this does to P2

| Claim in the draft | After this reading |
|---|---|
| A-lattice bond = that quintic, height 4 | **confirmed from the paper** |
| "every exactly-known planar threshold has height ≤ 3" | **false**, as already corrected on #564; the primary text is why |
| Result F′ (`C(≤6, ≤4)`) is the class that covers the record | **stands** |
| square site is outside the 2006 cell/dual-cell method | **Ziff says so**; this is Q4's second half as of 2006, not a proof for all future methods |
| Scullard martini site polynomial | **confirmed** |
| Suding–Ziff as a source of an exact square-site form | **they do not claim one** |

No census artifact moves. If #564's `verification_status` is flipped to a
primary reading, cite this note and the arXiv/journal pair above.

## Not established

- that height 4 is a bound on the mechanisms (it is a data point; see the
  companion literature brief);
- that square site cannot be solved by a later, larger-cell construction;
- anything about our lattice matching-odd observable.
