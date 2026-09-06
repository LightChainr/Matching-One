# Synthesis: what "state" means for rank-one continuation, and how it grows with the language

**Status:** the north-star answer, consolidated from the exact results of this
program.  On the one family where the experiment semantics is fully closed
(the #549 parallel-gadget family), the answer is complete and exact; the
boundaries to general cut networks and to positive realization are stated
separately.

---

## 0. The question, restated

> Given a rank-one torus continuation state, what is the minimal object that
> predicts every admissible shared-prefix / branching future, and how does its
> dimension grow as the experiment language becomes richer?

The program's exact work answers this by first *fixing the language*, then
computing the minimal predictive object **as a function of the language**, and
finally showing that no local/coarse summary can replace it.

## 1. The predictive-state hierarchy (complete, on the #549 family)

Let the hidden coordinate be `a = #{A-type gadgets}` in a bank of `k`
gadgets.  Define the depth-`d` language `L_d` = all experiments built from at
most `d` independent #549 forks (disjoint future groups).  The minimal linear
predictive state and the predictive-class count are then:

| language | minimal linear state | dimension | predictive classes |
|---|---|---|---|
| `L_0` (unbranched) | constant 1 | 1 | 1 |
| `L_1` (one fork) | `span{1, a}` | 2 | 2 |
| `L_d` (depth d) | `span{1, a, ..., a^d}` | `d+1` | `d+1` |
| `L_k` (full) | `span{1, a, ..., a^k}` = all functions of `a` | `k+1` | `k+1` |

Equivalently the **branching Hankel rank** of `L_d` equals the minimal linear
predictive-state dimension `d+1`, and the hierarchy refines strictly at every
depth and stabilises at depth `k`.  This is the exact content of the two
theorems in `branching-hankel-rank-lower-bound-20260906.md`:

```
1 = r_0 < 2 = r_1 < 3 = r_2 < ... < k+1 = r_k = (# classes).
```

**Answer.**  On this family the canonical predictive state is **the truncated
moment vector of the hidden type-count `(1, a, a^2, ..., a^d)`**; its dimension
is exactly the language depth plus one.  "State" is not a scalar attached to
the network — it is a ladder indexed by experiment depth, and each rung has an
exact finite dimension.  (A single real coordinate `a` carries all rungs, which
is why "k+1 classes" does not by itself force high dimension at shallow depth.)

## 2. What cannot replace the state (the no-go side)

* **No bounded-radius quotient** (`no-bounded-radius-quotient-20260906.md`):
  for every fixed `r`, the radius-r terminal neighbourhood plus the complete
  unbranched survival law `S(z)^k` is identical across all `k+1` classes while
  the fork probability separates them by `1/[2k(8k-1)^2]`.  The branching
  signal is carried by future vertices arbitrarily far from the terminals.
* **One witness ≠ rank** (`distinguishing-dimensions-and-positive-rank`):
  on the same `k+1` classes, one experiment already assigns `k+1` distinct
  values, yet the linear rank is 2; only `k+1` experiments give linear
  independence.  `{separating, affine, linear} = {1, 2, k+1}`.
* **Bounded one-branch memory, unbounded branching complexity**
  (`memory-branching-and-moment-algebra`): the `L_0` class is a single class
  (identical unbranched law), while `L_k` splits it into `k+1` classes.
  Shared-prefix branching measures counterfactual-future correlations that the
  marginal time series integrates out.

## 3. The algebraic mechanism (depth = tensor power, not moment order)

The single fork depends only on the **second** successor-hazard moment
(`sum x^2` = 29 vs 25), the first moment being equal.  Higher depth does not
summon higher hazard moments; it takes the **d-th tensor power of the
second-moment difference**, which projects onto `a^d`.  So depth is an
*exponent in a commutative (independent-fork) algebra*, not a derivative order
of a reliability functional (`memory-branching-and-moment-algebra`).

## 4. Structural programs — where they stand

* **Compositionality (J).**  Disjoint union is a congruence for every depth
  (success probabilities multiply).  Series composition / terminal gluing /
  vertex substitution are **open**; if they fail to be congruences, the
  missing state is the required counterexample.
* **Myhill–Nerode (G).**  The branching-Hankel-rank = minimal-linear-state
  correspondence holds on the parallel family (the only one with closed
  semantics).  A general weighted-tree-automaton statement needs the general
  tree semantics (Program A) closed first — open.
* **Positive vs signed (E).**  On the #549 response spaces nonnegative rank
  equals ordinary rank at every depth (a *negative* result: this family does
  not separate them, unlike P398).  A repository-native family with ordinary
  rank bounded but nonnegative rank unbounded remains **open**; the abstract
  slack-matrix mechanism (rank ≤ 3, nonnegative rank Θ(n)) is the known
  template to instantiate.
* **Canonical quotient existence (D, positive direction).**  The no-go rules
  out bounded-radius and reliability-polynomial summaries; the coarsest
  *exact* quotient preserving every branching future is still **open**, and is
  now the single sharpest remaining theorem target.

## 5. The north-star, in one sentence

For rank-one continuation, the canonical predictive object is **a
depth-indexed moment state `(1, a, ..., a^d)` whose dimension is the language
depth plus one**, it cannot be compressed by any bounded-radius summary, and
its only remaining open question is whether positivity (Program E) or a
general compositional semantics (Programs A/G/J) adds structure beyond this
linear ladder.

## 6. Boundaries

All results in Sections 1–3 are exact statements about the declared
composition of the published #549 fork protocol (per-group independence from
future-vertex disjointness) and about the occupied-bridge embedding
convention; the full N16 network lift needs the N16 site data and is not
claimed.  Ordinary real rank throughout; nonnegative rank equals it on the
only family computed.

Files: see `README.md` for the per-program matrix and the five notes.
