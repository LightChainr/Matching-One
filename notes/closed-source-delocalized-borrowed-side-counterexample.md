# A one-west contour need not contain a local peelable ear

## Result

The proposed geometric theorem is false.  A simple axis contour with
exactly one west step need not contain a bounded orthogonal ear whose
removal turns it into a directed contour.  The west step can borrow its
compensating east side across an arbitrarily long transverse run.

For integers

\[
 L\ge p+q+2,\qquad 1\le a\le L-1,
\]

consider the lifted path

\[
 \boxed{\Gamma_{L,a,p,q}
 =E^aU^pWU^qE^{L+1-a}D^{p+q}.}                    \tag{1}
\]

It runs from `(0,0)` to `(L,0)`, so it closes to an axis-essential loop
on the `L x L` torus.  It is simple, occupied-corner legal, and contains
one west step.  Any local surgery eliminating that west step must borrow
an east edge from one of its two sides; the nearer such side is at
transverse distance `min(p,q)`.  Hence no packet of uniformly bounded
diameter peels (1) when `p,q` diverge.

The smallest member with no unit orthogonal ear adjacent to the west
step is

\[
                         p=q=2,                    \tag{2}
\]

which exists already for `L=6` (for example `a=2`).  Its excess contour
length over a directed straight loop is ten.  More generally

\[
 |\Gamma|-L=2+2(p+q).                              \tag{3}
\]

Thus a topological injection

```text
one-west contour -> directed contour x fixed local packet
```

does not exist without recording an unbounded borrowed-side span.

This is a minimal obstruction, not a restoration of the old crude
`beta` bound.  Long borrowed sides pay their actual vertical source
weight: increasing `p+q` by one adds one up and one down edge and costs
another factor `m^-2`.  Therefore
the counterexample prevents a purely local ear proof, while leaving
open a weighted injection into

```text
directed contour x (span, finite packet),
```

whose span tail may still improve the asymptotic gate.

## 1. Simplicity of the lifted path

The six pieces of (1) occupy

\[
\begin{array}{c|c}
\text{piece}&\text{support}\\ \hline
E^a&[0,a]\times\{0\}\\
U^p&\{a\}\times[0,p]\\
W&[a-1,a]\times\{p\}\\
U^q&\{a-1\}\times[p,p+q]\\
E^{L+1-a}&[a-1,L]\times\{p+q\}\\
D^{p+q}&\{L\}\times[0,p+q].
\end{array}                                      \tag{4}
\]

Their interiors are disjoint.  Consecutive pieces share one endpoint;
the last endpoint `(L,0)` is the period translate of the first.  Since
`p+q<L`, the last vertical segment does not meet a period translate of
either middle vertical segment.  Therefore the projected loop is simple.

Every dual vertex used by the loop has degree two.  In particular no
alternating four-cut-edge vertex needs an ambiguous smoothing.  Pair
`Gamma` with the straight axis loop at height `p+q+1`; the two loops are
disjoint and bound an annulus.  Filling either annular side produces an
occupied NN configuration whose cut contour resolves to (1) under the
occupied-corner convention.  Thus (1) is an actual legal digital
boundary, not merely a formal lattice walk.

## 2. Why the west step has no bounded ear

Let `e_W` be the unique west edge, from `(a,p)` to `(a-1,p)`.  Removing
`e_W` leaves its two endpoints on different pieces.  A surgery that
eliminates the west displacement while preserving a simple axis loop
must pair it with an east edge crossing the same vertical strip
`[a-1,a]`.

There are exactly two relevant east sides:

* the edge on the lower horizontal piece at height zero, reached through
  the run `U^p`;
* the first edge on the upper horizontal piece at height `p+q`, reached
  through the run `U^q`.

Consequently every rectangle or orthogonal packet containing `e_W` and
a compensating east side has transverse diameter at least

\[
                         \min(p,q).                 \tag{5}
\]

If `p=1`, the lower side and `e_W` form a unit orthogonal ear; if `q=1`,
the upper side does.  Conversely, when `p,q>=2`, neither side produces
a unit ear or an L-triomino peel.  This proves both the minimality of
(2) and the absence of any uniform bound in the family `p,q->infinity`.

A local ear elsewhere on the contour is irrelevant to the proposed
injection: deleting it does not remove the unique west edge.  Any
one-step factorization that actually reduces the west count must include
one of the two compensating sides above and therefore obey (5).

## 3. Source cost and the remaining weighted route

The step counts are

\[
 E=L+1,\qquad W=1,qquad U=D=p+q.
\]

Hence (3) follows.  At fixed vertical skeleton, inserting the reversal
still carries the familiar relative `m^-2` west/east cost.  But a large
borrowed-side span requires the vertical stacks appearing in (1);
relative to the straight contour their full weight is
`m^[-2-2(p+q)]`.  Thus the obstruction has an exponentially decaying
span tail in the strong-source regime.

The correct possible replacement for the failed local-ear theorem is a
span-resolved injection.  One records `(p,q)`, removes the rectangular
borrowed side, and retains the directed remainder.  Whether the sum of
that tail weakens the `beta` gate is an analytic transfer question, not
a consequence of planar topology alone.

## Scientific boundary

This note disproves only the uniform local-ear/L-triomino claim.  It does
not show that delocalized borrowed sides survive the physical placement
average, compute their signed coefficient, or overturn the vanishing
single-hairpin coefficient in
`closed-source-single-hairpin-ward-identity.md`.
