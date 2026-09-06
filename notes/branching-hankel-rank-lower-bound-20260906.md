# Branching-Hankel rank lower bound: depth is the complexity resource

**Status:** THEOREM (exact, rational) — upgrades #549 from a class-count lower
bound to a genuine **response-matrix rank lower bound**, and pins the
depth-vs-rank hierarchy `r_d(k) = min(d+1, k+1)` on the parallel-gadget family.
This realises both named major outcomes of the cut-network probe:

* **A.** construct `k+1` admissible branching experiments whose response
  matrix on the `k+1` #549 predictive classes has exact rank `k+1`;
* **B.** for every fixed depth `d`, rank `<= d+1` independent of `k`, while
  the predictive class count `k+1` grows without bound — so branching depth
  is a genuine complexity resource, not a cosmetic re-parametrisation.

Script: `scripts/rank_lower_bound.py` (exact `Fraction` arithmetic);
data: `results/rank-vs-depth/latest.json`.

---

## 1. Protocol (declared composition of the #549 fork)

A **single #549 fork** on a group of `g` gadgets is: one shared update vertex
(uniform among the `8g` future vertices), then two independently continued
clones, success = both safe.  Its exact conditional probability is affine:

```
F(g, a_g) = [343 g^3 - 182 g^2 + 25 g + 4 a_g] / [8g (8g-1)^2],
slope = 4/[8g(8g-1)^2] > 0.
```

A **grouped-fork experiment** `E_{(g_1,...,g_j)}` partitions the `k` gadgets
into `j` future-disjoint groups and runs one independent fork per group,
requiring all `j` to succeed.  Given the per-group A-counts the groups are
independent, and the per-group counts are hypergeometric in the total count
`a`.

## 2. Lemma (the polynomial algebra)

For `E_{(g_1,...,g_j)}`, the success probability is a polynomial of exact
degree `j` in `a` with positive leading coefficient.

*Proof sketch.*  Each `F(g_l, a_l)` is affine in `a_l` with positive slope, so
the product `∏ F(g_l, a_l)` is a degree-`j` polynomial in the `a_l` with the
`j`-th mixed monomial coefficient `∏ slope(g_l) > 0`.  Averaging over the
hypergeometric distribution of `(a_1,...,a_j)` given `a` replaces each
monomial by a falling-factorial moment of `a`, a degree-`j` polynomial whose
leading coefficient is `∏_l (g_l / k)` times the mixed coefficient — nonzero.
`QED`

Concretely the experiment `E_j` with group sizes `(1,1,...,1,k-j+1)` gives
degree exactly `j` (verified by finite differences for `k = 2..8`).

## 3. The two theorems

**Theorem 1 (rank upgrade).**  Let `E_0` be the trivial experiment (constant 1)
and `E_1,...,E_k` the grouped-fork experiments of depth `1,...,k`.  The
`(k+1) x (k+1)` response matrix `M_{a,j} = P(E_j | x_a)` has exact rank
`k+1` for every `k >= 1`.

*Proof.*  By the lemma, column `j` lies in `span{1,a,...,a^j}` and has a
nonzero component on `a^j`.  Hence the columns span
`span{1,a,...,a^k}`, whose dimension is `k+1` (the values `0..k` give a
Vandermonde system).  `QED` — machine-verified by exact elimination at
`k = 2..8`.

**Theorem 2 (fixed depth has bounded rank).**  Let `L_d` be the linear span of
all grouped-fork experiments with at most `d` groups (depth `<= d`).  Then

```
r_d(k) := dim L_d = d+1        for 1 <= d <= k,
```

independent of `k`; only full depth `d = k` reaches the class count.  In
particular `r_0 = 1` (unbranched coordinate `a` invisible), `r_1 = 2` (the
#549 single-fork language, matching the earlier affine-obstruction result),
and the hierarchy is exactly

```
1 = r_0 < 2 = r_1 < 3 = r_2 < ... < k+1 = r_k = (class count).
```

So: **the class count grows like k+1, but every fixed depth sees only a
constant-size response space.**  Depth, not width, is the resource that
exposes the hidden multiplicity.

## 4. Relation to the repository facts

* #549's `F_{k,a}` is the `d = 1` member: one non-constant coordinate, rank 2.
* The affine obstruction recorded in the earlier probe (single-root fork
  language caps rank at 2) is exactly `r_1 = 2`.
* #435's gap `1/98` is the `k = 1` slope difference `F(1,1) - F(1,0)`.
* The `L_0 ⊆ L_1 ⊆ ... ⊆ L_k` equivalence-relation hierarchy therefore
  strictly refines at every step on this family, and stabilises at `L_k` with
  exactly `k+1` classes.

## 5. Boundaries and the next questions

* The protocol is the **declared composition of the published #549 fork**;
  the per-group events are independent because groups are future-vertex
  disjoint.  Lifting to the *full* N16 network (where two forks may touch the
  same future region or share nonlocal safe-set correlations) needs the N16
  site data and is not claimed here — the same caveat as the earlier probe's
  two-root note.
* This gives ordinary real rank.  **Nonnegative/positive rank** of the same
  response matrix is the next object (Program E): the columns here are not
  entrywise proportional, so nonnegative rank may exceed 2 even at depth 1;
  determining `r_+(L_d)` and whether it also equals `d+1` (or grows faster)
  is the immediate follow-up.
* "Depth" here counts independent fork groups.  Sequential shared updates
  *within* one group (a deeper single-lineage tree) are a different dimension
  and are left explicit rather than conflated with group depth.

Files: `scripts/rank_lower_bound.py`,
`results/rank-vs-depth/latest.json`.
