# Every embedded rank-one prefix has a two-terminal birth representation

The earlier two-port maps are not confined to two fortunate configurations.
For a graph genuinely embedded on the torus, a fixed rank-one occupied
prefix admits a constructive, exact two-terminal representation of its
entire subsequent rank-two birth event. **There may be many parallel
two-port components, not just one.** Their safety polynomials multiply.
Small treewidth is a separate computational property and is not asserted.

This statement applies to the black nearest-neighbor graph used by P334,
including its periodic multigraph edges with their actual lifted gains.
It does not automatically apply to a drawing with unjoined edge crossings,
such as adding both square-face diagonals while retaining them as crossing
abstract graph edges.

## Precise setting and statement

Let G be a finite graph embedded without crossings in
`T=R^2/P Z^2`, where P is an invertible integer period matrix and its columns
are the period vectors. Fix an occupied subgraph B whose homology image
has rank one. Subsequent occupation of a subset S of the vacant vertices
adds every available nearest-neighbor edge between occupied endpoints.
Write d for the number of vacant vertices.

Choose any essential occupied connected component K. The homology image
of K contains a primitive class ell and is exactly `Z ell`; all occupied
components have their homology image in that same line. Define

`q(v)=det(P ell,v)`, so `q(P h)=det(P) det(ell,h)`.

The following conclusions hold simultaneously for every S:

1. Every cycle of the full graph G minus K has q-gain zero.
2. Each connected component H of G minus K has, after a vertex gauge,
   at most two distinct attachment addresses to K. If there are two,
   their difference has magnitude `abs(det(P))`.
3. Split K into two fixed terminals for each two-address H, assigning
   each attachment to its address. Rank-two birth occurs exactly when
   at least one such network connects its two terminals.
4. Identifying the lower terminals of all these networks as s and the
   upper terminals as t produces one ordinary two-terminal connectivity
   event. Distinct components share only the fixed terminals; their
   selectable vertex sets are disjoint.

Zero/one-address components cannot cause this birth and are free variables
for safety counting. Other occupied components may be contracted because
all their cycles have q-gain zero. Contraction is used only after the
topological argument on the original embedded graph.

## Proof: the annular cover exposes the two sides

**Primitive essential cycle.** A nonzero closed walk in the finite embedded
component K decomposes into simple graph cycles. At least one is essential.
An essential embedded simple circle on the torus has primitive homology;
call its class ell after choosing sign. Since B has rank one, the image of
every cycle in B is collinear with ell. In particular K's image is exactly
`Z ell`, not a proper multiple of it.

**Balanced complement.** Choose that simple circle C inside K. Any graph
cycle avoiding K is disjoint from C. Its algebraic intersection with C is
therefore zero. In period coordinates this intersection is `det(ell,h)`.
Consequently q vanishes on every such cycle. This is a statement about the
fully available complement graph; it automatically holds for every
selected-vertex subgraph, without enumerating insertion subsets.

**Ordered lifts.** Pass to the annular cover
`A=R^2/(P ell)Z` of T. Its transverse deck group is Z. Because K's homology
image is exactly `Z ell`, its inverse image consists of disjoint connected
lifts K_j, indexed by that deck group. Each K_j contains an essential
simple circle C_j in the annulus. These circles are disjoint, ordered, and
separate the infinite annulus. A transverse deck translation sends j to
j+1; choose its sign to agree with this order.

A connected lift of G minus K cannot cross any C_j, so it lies between
two consecutive circles C_j,C_(j+1). It can attach only to K_j and
K_(j+1). Indeed a different K lift and any edge to it would have to cross
one of those separating circles; an embedded edge cannot do this without
meeting K. The fact that K can have branches and additional parallel
cycles does not remove its separating C_j. It can create pockets or
several disconnected channels between neighboring lifts, but not a third
neighboring lift for one connected complement component.

**Addresses.** The zero cycle-gain property lets us integrate q-gains to
a scalar potential phi on each finite complement component. Attachments
to a common K lift have the same gauge-adjusted address. Moving from K_j
to K_(j+1) adds `det(P)` to q, since a transverse deck generator h can be
chosen with `det(ell,h)=1`. Thus the possible addresses are a singleton
or two consecutive values separated by `abs(det(P))`. Different finite
components may have different additive gauge offsets; only their internal
two-address distinction is used.

**Birth equivalence.** A nonzero transverse cycle in B plus S must visit K,
because the complement is balanced. Split it into excursions outside K
and connecting portions inside K. A nonzero total q-gain requires an
excursion between different attachment addresses. That excursion is an
occupied path between the two terminals of one complement component.
Conversely, any such path, completed by an occupied path inside K, has
q-gain `+/-det(P)` and is independent of ell. It creates rank two.
This proves all four statements.

The proof also explains the familiar primitive cross topology: ell and
the resulting transverse cycle generate all of H1(T), not a proper
finite-index rank-two subgroup. This concerns the homology image of this
embedded occupied graph, not every quantity called an index in the
repository.

## Constructive gain formula

For each occupied component choose a root and a spanning-tree lifted
displacement p(v). Set p(v)=0 at a vacant singleton. An oriented original
edge u->v with displacement e becomes a root-to-root edge with gain

`g(u,v)=q(p(u)+e-p(v))`.

All contracted occupied cycles have zero q, so a different spanning tree
changes only the vertex gauge. On a complement component set
`g(u,v)=phi(v)-phi(u)`. An edge from the root of K to v has address

`a(v)=g(K,v)-phi(v)`.

An excursion between attachments u and v has gain `a(u)-a(v)`. Equal-gain
parallel edges and zero loops may be removed. Different-gain parallel
attachments must remain distinct until the two terminals have been split.
Deleting off-terminal block-cut-tree branches is then an ordinary exact
connectivity simplification, not a further topological assumption.

For the N425 family, `P ell=(8,-19)` and `q(dx,dy)=19dx+8dy`.
Hence the observed address difference 425 is forced by the construction.
One two-port component in each original witness was incidental: the new
frozen twelve-prefix analysis has between one and fourteen such components.
Those finite observations are illustrations, not premises of this proof.

## The full physical clock factors, without a hyperedge hierarchy

Let the relevant terminal networks have disjoint variable sets V_i and
safety polynomials

`F_i(z)=sum_{S subset V_i: s not connected to t} z^|S|`.

Let r other vacant vertices be irrelevant after the exact simplifications.
Then the complete physical safety polynomial is

\[
\boxed{F_B(z)=(1+z)^r\prod_i F_i(z).}
\]

This includes every trigger order. It is not necessary to enumerate
minimal pairs, triples, quartics, or higher hyperedges before obtaining
the full waiting-time law. With `f_k=[z^k]F_B` and uniform fresh ordering,

`Pr(T>k)=f_k/binom(d,k)`.

At fixed k the component survivals must **not** be multiplied directly:
their occupancy counts share the fixed total k. Polynomial multiplication
performs the required convolution exactly. Under independent Bernoulli
occupation with probability p, the component safety probabilities do
multiply directly.

For a marked site v in component i, its pivotal generating polynomial in
the full system equals its component pivotal polynomial multiplied by all
other component safety polynomials and the free-variable factor. A first
birth attributed to v requires all competing components still to be safe.
This gives the local-to-global marked-birth bridge without duplicating
independent propagation paths as evidence.

### Which parallel channel wins the birth race?

There is an aggregate marked law requiring **no extra reliability solves**.
Let n_i be the number of random vertices assigned to component i, not the
degree of its safety polynomial, and form

`B_i(z)=n_i F_i(z)-(1+z)F_i'(z)`.

Its coefficient at k is `(n_i-k) f_i,k-(k+1) f_i,k+1`: the number of
safe-set/next-vertex pairs that cross into connectivity at that step.
These coefficients are nonnegative. The full-system contribution of
component i is

\[
H_i(z)=(1+z)^r B_i(z)\prod_{j\ne i}F_j(z),\qquad
\boxed{\Pr(T=k+1,\ I_{\rm birth}=i)=
\frac{[z^k]H_i(z)}{d\binom{d-1}{k}}.}
\]

The product rule gives `sum_i H_i=dF_B-(1+z)F_B'`, so summing the component
laws reproduces the entire waiting-time distribution. In the full NN torus,
where complete occupancy necessarily creates rank two, their integrated
probabilities sum to one. Long-tail channel shares use the same coefficients
restricted to k+1 above the declared tail step; no site-by-site force-on/off
calculation is required.

Equivalently, with Bernoulli safety
`S_i(p)=(1-p)^(n_i) F_i(p/(1-p))`, the probability that channel i produces
the first birth is `integral_0^1 [-S_i'(p)] product_(j!=i) S_j(p) dp`.
This is a competing-channel attribution, not a claim that channel waiting
times are independent under a fixed total insertion count.

## What is new here, and what is not

The established torus homology classification is the backdrop, not a new
CFT result. Arguin, *Homology of Fortuin-Kasteleyn clusters of Potts models
on the torus*, section 2, uses the trivial, primitive cyclic, and full-cross
homology sectors. Its finite embedded-graph topology applies to the NN
occupied subgraph here, without importing FK probabilities into site
percolation. The constructive annular-cover/address proof and the
conditional reliability factorization above are supplied explicitly here.

Primary reference, retrieved with the arXiv skill:
[Arguin, hep-th/0111193v2](https://arxiv.org/abs/hep-th/0111193v2),
J. Stat. Phys. 109, 301--310 (2002), DOI 10.1023/A:1019979326380.

- Exact scope: finite graph embedded on a torus, fixed rank-one prefix,
  original lifted NN incidences retained, arbitrary subsequent site subset.
- No claim: general small treewidth, a cheap evaluator for every prefix,
  applicability to unjoined crossing edges, a new exponent, or a CFT
  operator identification.
- Computational consequence: topology reduces the full event to parallel
  ordinary reliability networks; the remaining cost is their actual width
  and connectivity, not the largest minimal-trigger cardinality.
- Source lineage: the instance maps and clocks in `6358ba49`; new-prefix
  selection fixed at `b9cbe13e`. No new random data or repeated network
  computations were used to derive this statement.
