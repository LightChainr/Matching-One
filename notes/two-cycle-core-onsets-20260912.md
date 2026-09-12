# Two-cycle cores settle the axis and diamond rank-two onsets

Date: 2026-09-12. Direct completion of the remaining geometric question in
#673 / draft #692. This is a proof, with small existing-size checks, not a
proposal to enumerate another torus. No novelty claim is made.

## Result and conventions

On the nearest-neighbour square-site graph with period lattice

* axis: `Lambda = <(L,0),(0,L)>`,
* diamond: `Lambda = <(L,L),(L,-L)>`,

let `r(B)` be the rank over Q of the occupied graph's image in torus H1.
For every integer `L >= 2` (parallel lifted edges at period two retained),

| Geometry | Least occupied mass with r=2 | All minimizers | Number |
|---|---:|---|---:|
| axis | 2L-1 | a full physical horizontal row union a full vertical column | L^2 |
| diamond | 3L-1 | a full physical horizontal/vertical line of 2L vertices, plus the L-1 internal vertices of one straight transverse length-L arc | 4L^2 |

The existing digital-Alexander theorem identifies these as the `b x n`
(exclusive cross) cell: black rank two forces complementary matching rank
zero. It is NOT the rank-one spiral / `both x both` cell.

In the `(u,v)=(x+y,y-x)` convention of #692, a physical horizontal or vertical
line is a straight diagonal line. A physical transverse arc is the straight
plug described empirically there. Thus this theorem addresses exactly the
existing conjecture, not a different lattice or event.

## 1. The reusable two-cycle-core lemma

If a graph embedded in a torus has ambient rank two, it contains a connected
subgraph of cycle rank two whose two graph-homology generators have independent
ambient images. Indeed, different connected components cannot carry independent
ambient classes: their disjoint embedded cycles have intersection number zero,
whereas independent classes in H1(T^2;Q) have nonzero intersection. In a component
carrying rank two choose a spanning tree and two off-tree edges whose fundamental
cycles have independent ambient images. Retain the tree and those two edges,
then prune leaves. The graph cycle rank is two and the ambient map is injective.

After suppressing degree-two vertices, the connected core is either a wedge
of two circles, a theta (three internally disjoint paths between two vertices),
or a dumbbell (two circles connected by a path). This classification follows
from `sum(deg(v)-2)=2` after leaves are removed. A dumbbell is impossible: its
two disjoint circles would have independent ambient images. Thus only wedge
and theta remain. All their simple cycles have nonzero ambient class; for a
theta the three classes are, up to signs, `a,b,a-b`, with a,b independent.

This lemma concerns a SUBGRAPH of the occupied induced graph. Extra occupied
edges cause no problem: a lower bound on this core is already a lower bound
on occupied vertices. At equality there can be no additional occupied vertex.

## 2. Axis proof, including the missing inequality in #692

A simple cycle of class `(a,b)` has at least `L(|a|+|b|)` edges, by its physical
integer displacement. For a simple cycle the edge count equals its number of
vertices, including a two-edge periodic circle.

A wedge has two essential circles, each of length at least L, meeting in
one vertex; hence it has at least `2L-1` vertices. Equality forces both circles
to have length L. Independence forces one horizontal and one vertical, and
geodesic equality forces each to be straight. This gives precisely a row-column
cross.

For a theta, among `a,b,a-b` at least one class has both coordinates nonzero.
Otherwise the independent a,b would have to lie on the two different axes,
and a-b would be mixed. That simple circle has at least 2L vertices. A theta
therefore cannot occur at occupied mass `2L-1` or below.

The L horizontal rows and L vertical columns give L^2 different crosses,
each with rank two. This proves the onset, equality classification and count.

The original #692 inference, “the two cycles share at least one vertex, so
union size >= 2L-1”, has its inequality in the wrong direction: a lower bound
on intersection size gives an UPPER bound on a union when sizes are fixed.
The core lemma repairs the proof without invalidating the observed onset.
It also avoids assuming that every rank-two graph contains separately chosen
simple cycles in the two prescribed coordinate classes.

## 3. Diamond lower bound

For `Lambda=<(L,L),(L,-L)>`, a period vector is

    L(a+b, a-b),

whose l1 length is `2L max(|a|,|b|)`. Every nontrivial simple cycle therefore
has at least 2L edges.

A wedge consequently needs at least `4L-1` vertices. In a theta let the three
path lengths be `l1,l2,l3`. Each pair is an essential simple circle, so

    l1+l2 >= 2L,   l1+l3 >= 2L,   l2+l3 >= 2L.

Adding gives `E=l1+l2+l3 >= 3L`. A connected theta has `V=E-1`, hence
`V >= 3L-1`. This proves the lower bound on every occupied configuration.

At equality the core must be a theta and all three pair bounds are equalities,
so `l1=l2=l3=L`. The remaining issue is whether bent geodesic paths create
additional equality cases. The period arithmetic rules them out.

## 4. Diamond equality: three straight arms, not arbitrary bent plugs

Orient the three length-L paths from the same branch vertex u to v. Lift them
from the same lattice representative of u and write their physical displacement
vectors as d1,d2,d3. Their pairwise differences are nonzero period vectors,
so

    2L <= ||di-dj||_1 <= ||di||_1+||dj||_1 <= 2L.

Every inequality is an equality. Thus `||di||_1=L`, and no two displacement
vectors use the same coordinate with the same nonzero sign: otherwise the
triangle inequality for their difference would be strict.

There are only four signed coordinate slots: x+, x-, y+, y-. Three nonzero
vectors cannot each use two slots without overlap. At least one vector uses
only one coordinate and hence equals `(+-L,0)` or `(0,+-L)`. All di are congruent
coordinatewise modulo L, because their differences lie in Lambda. Therefore
all their coordinates are multiples of L. With l1 norm L, EACH di must be
one of these four axial vectors. They are distinct, so exactly two are an
opposite pair and the third is perpendicular.

A nearest-neighbour path of length L with displacement `(L,0)`, for example,
uses L positive horizontal steps; there is no room for a detour. All three
paths are straight. The opposite pair forms a full line of 2L vertices; the
remaining arm contributes exactly L-1 internal vertices of a transverse plug.

Conversely, every such line-plus-plug is a theta with independent essential
cycles, so it has rank two. This proves the full equality classification.

## 5. Exact count 4L^2

There are two physical line directions and L disjoint full lines of each
direction. For a given full line, choose one endpoint on its 2L vertices and
one of two transverse directions. Reversing the same length-L arc counts it
twice, so there are `(2L*2)/2=2L` distinct plugs per full line.

No occupied minimizer is counted under two full lines. Parallel full lines
are disjoint. A horizontal and a vertical full line on the diamond quotient
meet in two vertices; their union has `4L-2 > 3L-1` vertices for L>=2. The full
line in a minimizer is therefore unique. The count is

    2 * L * 2L = 4L^2.

In particular the previously proposed diamond L=5 onset is now a theorem:
`k=14`, count 100. No `2^50` enumeration, or even a new low-k census, is needed.

## 6. Other #692 statements that the proof repairs or narrows

At diamond mass 2L, rank two is impossible by the new bound. A rank-one cycle
with both generator coordinates nonzero must attain the systolic bound with
both `(u,v)` displacements of magnitude 2L. Every step then has the same
sign pattern, forcing one of the 2L full straight diagonal lines. The
`both x both` spiral count 2L at this mass follows, using the existing exact
rank-one label map for the complementary graph. This replaces #692's other
use of the incorrect shared-vertex union argument.

Two local wording errors in its high-k proof should also be corrected:
there are L lines PER diagonal family (2L total), not 2L per family. For the
sharpness example with exactly L white sites on the diamond, use an NNN
chain with physical steps `(1,1)` (or `(1,-1)`) closing after L steps, not a
full NN straight line of 2L sites. The high-k conclusion itself survives.

This note does not assert a general classification of all higher-mass cells,
Galois properties of matching polynomials, or any critical scaling exponent.

## 7. Executed checks and source trail

`python scripts/rank_two_onset.py` exhaustively checks all occupied subsets
AT OR BELOW the predicted onset on axis L=2,3,4 and diamond L=2,3: 133,711
configurations total. A direct physical integer-lift graph traversal finds
zero rank-two configurations below the bound and exactly the predicted sets
at equality. Counts are 4/9/16 on the axis and 16/36 on the diamond.

Every predicted minimizer at L=2,...,8 was also constructed and its rank
checked, without enumerating the surrounding configuration space. These are
checks of the proof and constructor, not new stochastic evidence.

Source read through the connector: draft #692,
`notes/wrapping-five-cell-onset-proofs-20260908.md`, blob
`6a3ac9dcf6e1c7f8bbef853379a0e02ea8ef4cad`, especially sections 3b and 5c.
The geometry there agrees with the physical periods used here. The cell
translation uses merged #702, not the superseded #690 directional conjecture.
Full repository CI for this new file has not been run.
