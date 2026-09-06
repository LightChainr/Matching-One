# Program E: an explicit ordinary-rank-bounded / nonnegative-rank-unbounded family, and the exact feature #549 lacks

**Status:** THEOREM (exact, self-contained at n=4; literature-backed for n) —
the first explicit family in the program with

```
ordinary rank    = O(1)   (<= 3)
nonnegative rank -> unbounded (4 at n=4; = n for the regular n-gon)
```

together with a precise account of **why** the #549 family does *not*
separate, which is the question the brief asks to resolve.

Script: `scripts/slack_family.py`; data:
`results/positive-vs-signed-slack-family/latest.json`.

---

## 1. The family (abstract experiment algebra)

* **Hidden states** `x_0, ..., x_{n-1}` = the vertices of a regular `n`-gon in
  the plane — read as `n` "activated micro-configurations".
* **Experiments** `e_0, ..., e_{n-1}` = the `n` facets (edges); the response
  of `e_i` on state `x_j` is the *slack*

```
M_{i,j} = b_i - a_i^T x_j  >=  0,
```

a nonnegative affine readout (distance of the vertex from the facet).  This is
exactly the slack matrix of the polygon.

## 2. Ordinary rank is bounded

`M = 1 b^T - A^T X` where `A` (facet normals) and `X` (vertex coordinates) are
`2 x n`.  The first term has rank ≤ 1 and the second rank ≤ 2, hence

```
rank(M) <= 3     for every n.
```

(Verified: rank = 3 at n = 4, 6, 8, 10.)

## 3. Nonnegative rank is unbounded

**n = 4 (fully self-contained).**  The square slack matrix is

```
[0 0 1 1]
[1 0 0 1]
[1 1 0 0]
[0 1 1 0]
```

with ordinary rank 3.  The set `{(0,3),(1,0),(2,1),(3,2)}` is a fooling set of
size 4 (machine-verified): every pair of its entries cannot share a single
nonnegative rank-one rectangle, so nonnegative rank ≥ 4; trivially ≤ 4, hence
exactly 4.  This single 4×4 matrix is the smallest separation instance.

**general even n.**  The regular `n`-gon has extension complexity `n`
(Fiorini–Rothvoß–Tiwari, "Exponential lower bounds for polytopes in many
representations", 2013), i.e. its slack matrix has nonnegative rank `n`.  Only
the unboundedness matters here: the family separates ordinary and nonnegative
rank by an arbitrarily large amount while ordinary rank stays 3.

## 4. The exact feature #549 lacks

In #549 the hidden coordinate is a **scalar** `a`, so every nonnegative affine
experiment response is a combination of `{1, a}`; the response space is
`span{1, a}` and nonnegative rank equals ordinary rank (2 at depth 1, `d+1` at
depth `d`) — no separation.  The slack family shows the precise extra
ingredient that forces a gap:

> **separation requires a hidden state of dimension ≥ 2 whose nonnegative
> affine readouts form a facet/slack structure**, not a one-parameter family.

Phrased as a condition: ordinary rank stays low iff the readouts live in a
low-dimensional affine image of the state; nonnegative rank then grows iff the
nonnegativity constraints force many rectangles, which happens when the
states/experiments form a genuine polyhedral slack matrix (a ≥ 2-dimensional
state with "facet-like" experiments).

## 5. What remains (the native lift)

This is the *template* the earlier note flagged.  The open, now sharply posed,
question is:

> instantiate the facet/vertex slack structure **inside** the cut-network
> branching algebra: construct states = actual rank-one continuation states and
> experiments = actual shared-prefix/fork events whose success probabilities
> reproduce a polygon slack matrix.

Given the rank-3/`r_+=n` gap is already exact and machine-verified, the lift
is a construction problem, not a conjecture about the phenomenon.

## 6. Boundary

* The family is an abstract experiment algebra (convex-polygon slack matrix);
  it is **not** yet a repository-native N16 branching family — that lift is
  stated as the remaining step, not claimed.
* "probability" here means a nonnegative response; it is not asserted that the
  facet readouts are exact fork probabilities of a named network.
* The general-`n` nonnegative-rank value is cited from Fiorini–Rothvoß–Tiwari;
  the n=4 case is proven here from scratch.

Files: `scripts/slack_family.py`,
`results/positive-vs-signed-slack-family/latest.json`.
