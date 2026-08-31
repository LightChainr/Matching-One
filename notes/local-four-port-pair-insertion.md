# A local four-port pair insertion has an exact two-component occupation mark

**Result.** Replace the vacant-site tensor at x by the specified pair
projector, and average its four C4 rotations. Its exact relative
occupation weight is

\[
 \boxed{\beta_x(Q)=\mathbf1_{x\notin A}\frac{d_2(Q)}{Q^2}
 \left[\mathbf1_{NS\mid EW}
       +\tfrac12\mathbf1_{NE\mid SW}
       +\tfrac12\mathbf1_{NW\mid ES}\right],
 \quad d_2(Q)=\frac{Q(Q-3)}2.}
 \tag{1}
\]

The three events are **exact exterior component partitions** of the
four incident edge-nodes, not merely the indicated pairwise
connections. Each requires all four neighbours of the vacant site
to be occupied and exactly two distinct occupied NN components to
meet those four neighbours, twice each. Other occupied components
away from x are unrestricted.

This supplies an actual local insertion interface. Unlike the completed
global one-seam `[2]` trace, it can have nonzero direct q/E numerators.
It does not by itself identify a continuum four-leg primary or the
cause of a normalized global-U response.

Base: `bea717e826df5a22518774b1725ae7bcbe2cb801`.
Conventions are those of
[the hypergraph law](closed-source-hypergraph-rc-twist-projection.md)
and [finite torus pair closure](closed-source-finite-torus-pair-closure.md).
All arguments below are finite algebra or explicit lattice drawings in
coordinates; no Monte Carlo or enumeration is used.

## 1. Ordered ports and the unmodified physical configuration

At a site x let the incident edge colours be N,E,S,W. The original
hypergraph tensor is `1+v delta_all4`. In its vacant term, replace 1
by `Pi2=i P_[Q-2,2] i^dagger`, acting from the ordered input pair
(N,E) to output pair (S,W), where

\[
 i|\{a,b\}\rangle=(|a,b\rangle+|b,a\rangle)/\sqrt2,
 \qquad a\ne b.
 \tag{2}
\]

The projector acts as zero on the equal-colour diagonal. Thus
`Pi2_(a,b),(c,d)=P2_{ {a,b},{c,d} }/2` when both ordered pairs have
distinct colours, and zero otherwise. Work first at integer Q>=4;
the contracted formulas below give their rational continuation.

For every term, x is **vacant in the original A**. Its hypergraph
cluster count c_H, winding subgroups, ambient rank r, q=r-1 and E=q²
are those of A. In particular, projector diagram lines must not be
treated as new occupied NN edges when computing r. The following
contractions do not include colour twists.

## 2. All fifteen exterior partitions, without a numerical table

Let p be the partition of the four edge-nodes by their exterior
hypergraph connectivity. Distinct blocks get independent colours;
their colours are allowed to coincide. If p has b blocks, the other
hypergraph clusters supply `Q^(c_H-b)`. Let C(p) be the remaining
sum of the projector tensor with p's equality constraints.

Write B for the pair-to-colour incidence matrix,
`B_{ {a,b},c}=1_{c in {a,b}}`. Its range is singlet plus standard,
so `P2 B=0`. Also `P2 1=0`. These two identities settle all cases:

| exterior partition | contraction C(p) | reason |
|---|---|---|
| one block; any 3+1 | 0 | some input or output pair is equal |
| 2+1+1, pair NE or SW | 0 | equal input or output pair |
| 2+1+1, a cross input-output pair | 0 | the sum is an incidence-vector contraction of P2 |
| four singleton blocks | 0 | invariant pair vector is killed by P2 |
| NE \| SW | 0 | both same-side pairs are equal |
| NS \| EW | d2(Q) | `Tr Pi2` |
| NW \| ES | d2(Q) | `Tr(Pi2 Swap)=Tr Pi2` |

For example, N=S=a with E and W otherwise unconstrained gives
`(1/2) sum_a B_a^T P2 B_a=0`. This also shows why a 2+1+1 partition
cannot be counted as a successful pairing just because one cross
connection is present.

The two surviving partitions have b=2. Their full colour contraction
is therefore exactly

\[
 d_2(Q)Q^{c_H-2}.
 \tag{3}
\]

The physical factor `v^K Q^(-r/2) O(r)` is unchanged, with no extra
factor v for the vacant marked site. Relative to the original weight
`v^K Q^(c_H-r/2)`, each successful unrotated pairing has factor
`d2(Q)/Q²`.

The opposite pairing NS|EW survives all four rotations. Each adjacent
pairing survives two rotations and vanishes in the other two. The
average is precisely (1); the adjacent events are disjoint, so the
notation `1_(NE|SW or NW|ES)/2` is equivalent to the displayed sum.

An incident edge-node at a vacant x is isolated if its other endpoint
is vacant. If that neighbour is occupied, the edge-node is a dangling
leaf of that neighbour's NN component. Hence a 2+2 partition is
equivalent to the four-neighbour/two-component condition stated above.
No nearby colour-only cluster has been substituted for a physical NN
component.

## 3. Rank support: adjacent pairings occur at every original rank

Here are explicit configurations on an LxL square torus, L>=7. All
coordinates are modulo L, x=(0,0) is vacant, and sites not listed are
vacant. Define

\[
 C_{NE}=\{(0,1),(1,1),(1,0)\},\qquad
 C_{SW}=\{(0,-1),(-1,-1),(-1,0)\}.
\]

Take the following occupied sets:

| original rank | occupied A | two components at x |
|---|---|---|
| 0 | `C_NE union C_SW` | two contractible three-site paths |
| 1 | preceding A, plus the full row y=2 | NE belongs to the horizontal winding component; SW stays contractible |
| 2 | preceding A, plus the full column x=2 | NE belongs to a rank-two component; SW stays contractible |

In every case the components remain distinct, the exterior partition
is NE|SW, and the bracket in (1) is 1/2. Thus the insertion has support
at q=-1,0,+1. In particular its unnormalized E numerator is strictly
positive at integer Q>=4 and positive activity on these tori: the
rank-zero and rank-two examples have positive weights, and all mark
weights are nonnegative. No corresponding sign for a normalized
covariance or global-U derivative follows from that fact.

### The opposite pairing has an additional topological meaning

For NS|EW choose a simple occupied path joining N to S in one
component and a simple occupied path joining E to W in the other.
Temporarily adding x closes two simple cycles. They meet only at x,
with alternating incident directions, and their algebraic intersection
is ±1. They therefore generate the full ambient torus homology:

\[
 NS\mid EW\quad\Longrightarrow\quad r(A\cup\{x\})=2.
 \tag{4}
\]

Moreover r(A) cannot already be two. A rank-two occupied component
would be disjoint from at least one of those two essential cycles:
if it is one of the two marked components use the other cycle, and
if it is an unmarked component use either. But the complement of a
rank-two embedded graph contains no essential curve; equivalently,
two independent classes in that graph cannot both have zero
intersection with a nonzero disjoint class. Hence r(A) is 0 or 1.

Both occur concretely. Occupy the row y=0 and column x=0 except for
their common missing site x. They are two disjoint open paths,
giving original rank zero and partition NS|EW. Adding the full row
y=2 joins the vertical path to a horizontal winding cycle while
leaving the other path separate; the original rank is then one and
the same partition persists. These examples work for L>=5.

Equation (4) classifies the exterior event; it does **not** instruct
the insertion calculation to use the rank after adding x. Its
physical q/E weight still uses the original rank.

## 4. A genuine separation between the global seam and a local pair mark

Occupy exactly two complete horizontal rows, with vacant gaps between
them. They are two essential components of winding (1,0). For the
primitive horizontal seam their global pair trace is

\[
 v^{2L}Q^{c_0-1/2}d_2(Q),
 \tag{5}
\]

which is nonzero for integer Q>=4, and its rational continuation is
not identically zero. Yet (1) vanishes **at every site**: at an
occupied site the vacant insertion is absent, while at every vacant
site both horizontal neighbours are vacant. This already proves the
global/local separation for the actual unit-site experiment.

The separation also persists for a fixed-size, single-pair local
cut. Fix R and place the rows at torus distance greater than 2R from
each other, with L large enough that a radius-R disc is injective.
A disc contains at most one row segment. For the same one-projector
construction on four distinct exterior edge ports, with all exterior
connections supplied by the unmodified occupied configuration, the
exterior partition has at most one nonsingleton block: the disc cuts
that row in one interval, whose remaining arc is connected outside.
The other row does not meet the ports, and unused edge-nodes are
isolated. None of the required 2+2 partitions occurs. The local
contraction is zero by the table in section 2.

This last statement concerns a simply connected disc cut, one copy
of Pi2, four distinct physical edge ports, and no extra spectator
wiring that creates a second cluster. It does not cover a support
wrapping the torus, two spatially separated insertions, a deliberately
added colour-loop closure, or arbitrary local fields. In particular it
does not assert that every local observable has zero coupling to the
two-row configuration. It proves that the completed seam signal alone
does not establish a nonzero coupling for this specified local pair
insertion.

## 5. Finite occupation packet and the Q=1 endpoint

Let `kappa_x` be the bracket in (1), including the vacancy indicator.
The complete first-order local packet is the finite occupation sum

\[
 \mathcal N^{loc}_O(Q,v)=
 \sum_A v^KQ^{c_H-r/2}O(r)\sum_x\beta_x(Q).
 \tag{6}
\]

Its denominator packet is obtained by O=1. A normalized source
response must subtract the denominator contribution, rather than
identify (6) alone with that response. Since

\[
 \beta_x(Q)=\frac{Q-3}{2Q}\kappa_x,\qquad
 \beta_x(1)=-\kappa_x,\qquad
 \partial_Q\beta_x(1)=\tfrac32\kappa_x,
 \tag{7}
\]

the contracted endpoint and its explicit operator derivative are
regular. A derivative of the entire Q-family additionally differentiates
its declared occupation weight and activity; it is not just the last
term of (7). This is a four-port trace contraction through two exterior
components, so it is fully compatible with the regular one-endpoint
identity `ell P_[2]=0`. A nonzero finite local packet does not settle
its continuum normalization, field identity, or global anomaly loading.
