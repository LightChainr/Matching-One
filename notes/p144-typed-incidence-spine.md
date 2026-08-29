# A binary typed incidence-spine state sum for square-site matching

## Result

The post-#269 local-state route closes, but not as two pairwise smoothings of
one four-valent vertex.  It closes as two **typed four-way junctions** on a
fixed doubled-lattice ribbon/incidence spine.

Place three node types on the torus:

- original sites at even-even doubled coordinates;
- edge-midpoint hubs at odd-even/even-odd coordinates;
- face-center hubs at odd-odd coordinates.

At every original site choose one of exactly two local states:

```text
B = J_edge: join the four incident edge-midpoint ports;
W = J_face: join the four incident face-center ports.
```

The `B` spine retracts to the black NN induced graph.  In each square, the
`W` face hub replaces the clique among its white corners by an embedded star.
It therefore has the same components and ambient-H1 image as white
NN+NNN connectivity.  The 16 face patterns are checked exactly.

This is the combinatorial spine of the two complementary subsurfaces in the
digital Alexander proof.  With the #269 cellwise thickening, the black spine
retracts from `U`, the white spine is a 1-skeleton of `V=closure(T^2\\U)`, and
the subsurfaces `U,V` share the boundary.  This statement is about that typed
cellwise thickening, not arbitrary metric neighborhoods of the bare graphs.
Complementing the site mask exchanges `B` and `W` locally.
It does not imply that the axis-lattice matching polynomial is complement-odd,
because the edge-hub and face-hub types are distinct.

## Rank-graded vertex/surface polynomial

The fixed-spine state sum is

\[
\Phi_L(p;x,y)=\sum_{\omega}
p^{|B|}(1-p)^{N-|B|}x^{r_b(\omega)}y^{r_w(\omega)}.
\]

By #269, `r_b+r_w=2`, and the matching polynomial is the exact specialization

\[
M_L(p)=\left.\frac12(x\partial_x-y\partial_y)\Phi_L(p;x,y)
\right|_{x=y=1}.
\]

On the honest axis `L=3` torus, exhaustive enumeration of all 512 masks shows:

- the typed spine and the original NN/NN+NNN graphs have identical component
  counts and ambient ranks for every mask;
- `q=(r_b-r_w)/2` is read without an event-channel convention;
- differentiating the rank-graded state sum reproduces every Bernstein
  coefficient of the direct matching polynomial.

The common boundary descriptor has only three cases.  Essential parallel
boundary curves give `(r_b,r_w)=(1,1)` and `q=0`.  Contractible boundary curves
have one torus side: black gives `(2,0),q=+1`, white gives `(0,2),q=-1`.

## Why literal two-smoothing language fails

Two smoothings of a four-port medial vertex give only the partitions

```text
(NE)(SW),  (NW)(ES).
```

A black site with four incident black arms requires the one-block terminal
partition `(NESW)`.  No pairing represents this matrix element.  The minimal
missing partition type is a four-way junction `J4`.  White 8-connectivity needs
the same junction class on the alternating face-port set.

Thus the useful exact statement is sharper than “two smoothings”: the mask is
binary, but its two values select `J_edge` and `J_face`.  A pure transition
polynomial restricted to pairings discards connectivity information needed by
the matching charge.

## Minimal transfer rule

For a fixed-width scan, keep two typed frontier partitions, integer lift
potentials, and the accumulated ambient-H1 generators.  At a new site branch:

```text
B: union the site with its four edge hubs;
W: union the site with its four face hubs.
```

When the periodic frontier closes, emit
`p^|B|(1-p)^(N-|B|) x^r_b y^r_w`.  This is a finite recursion at every fixed
width.  No named vertex polynomial, smaller deletion/contraction algebra, or
width-independent bond dimension is claimed.
