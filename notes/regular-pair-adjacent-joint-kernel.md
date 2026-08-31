# Adjacent vacant vertices: one shared edge plus an outside bypass

The existing Bell8 kernel applies to adjacent marked vertices without
alteration, provided their common undirected edge is represented by the
**same** hypergraph-node ID. That node is an isolated outside component
and is counted once in the normalization. Adjacency alone gives zero:
a nonzero kernel requires at least one additional outside component
connecting the two vertices. One such bypass gives a nonnegative
contrast product; multiple bypass components can give negative contact
terms. The homogeneous contact source is not a thermal or one-site-source
alias, as an occupied-domino example proves exactly.

Base `a237968f1d7a82d26b46e83c58179dbba7f1a908`. This uses the fixed
canonical Kreg and original occupation/rank observers. No Bell8 table,
N25 score, double-insertion norm, or Monte Carlo run was repeated.

## 1. The shared physical edge and its normalization

Take x=(0,0) and y=(1,0); other directions follow by rotation. The port
order remains

```
(xN,xE,xS,xW,yN,yE,yS,yW).
```

When both marked sites are vacant, the edge-node e={x,y} is present
once in the original hypergraph. Neither endpoint activates its
four-way join, and no other site can join that edge-node. Consequently
the external partition contains the block

```
{xE,yW},
```

with no other port in that block. It is a singleton physical component
represented twice in the port list, not two independent singletons.
In particular its colour is one shared index, not two summed colours.

Let pi' be the partition of the six other ports and b'=|pi'|. The full
eight-port partition has b=b'+1 blocks. Writing their colours as
`(a,e,c,d,f,g,h,e)`, the exact normalized closure is

```
B_pi(Q) = Q^(-(b'+1)) sum_(six external colours,e)
              Kreg(a,e,c,d) Kreg(f,g,h,e) D_pi'
        = Q^(-b') sum_(six external colours) L_Q D_pi',

L_Q(a,c,d;f,g,h)
        = Q^(-1) sum_e Kreg(a,e,c,d) Kreg(f,g,h,e).           (1)
```

Thus the usual `Q^(-|pi|)` is exactly correct: |pi| counts **unique
outside components**, not port appearances. Equation (1) is a
six-external-port contact tensor with one averaged internal edge colour;
there is no extra vacancy probability or extra factor Q to insert.
Duplicating the ID would change the contraction itself. Merely changing
the normalizing power by one, while retaining the correct contraction,
would happen not to change its first Q derivative because B_pi(1)=0;
that coincidence does not justify the wrong finite-Q normalization.

All other external factors cancel as in the Bell8 construction. The
original `Q^(-r(A)/2)`, q(A), and E(A) use the **unmodified occupation A**.
Do not infer a new topological rank from the virtual joins in L_Q.
The joint insertion is zero if either marked site is occupied.

## 2. Adjacency alone is zero; a bypass is necessary

The block {xE,yW} supplies one shared outside component. If it is the
only shared component, colour invariance makes each local conditional
factor independent of that shared colour. The normalized pair closure
therefore factors into two one-site closures. Both vanish at Q1, giving

```
no other shared outside component  =>  g_pi=0.               (2)
```

Any other component meeting both vertices contains an occupied path
between their other neighbouring sites, avoiding x and y. For an honest
simple square quotient there is no second physical edge with the same
two endpoints. Thus a nonzero adjacent kernel requires a genuine
**outside bypass**, not simply the presence of the marked contact edge.

If precisely one additional component is shared, there are two shared
colours in total. Let d_x,d_y be the equal-minus-unequal conditional
colour contrasts from the existing spatial-kernel derivation. Then

```
g_pi=d_x(1)d_y(1).                                         (3)
```

The contact-edge component occupies exactly one port at each vertex.
It cannot participate in the positive local 2+2 shared-component case.
Every possible local contrast here is therefore in

```
{0, -1/2, -3/4, -1, -3/2}.
```

This proves the sharper adjacent selection rule

```
exactly one outside bypass component => g_pi>=0,             (4)
```

with zero still possible. For the simplest north bypass, the two
neighbours xN and yN belong to one occupied component, while the
remaining four outside ports are private. Its canonical labels are
`01230451`. The shared edge is adjacent to the bypass port at each
vertex, so d_x=d_y=-1/2 and

```
g_pi=1/4.                                                  (5)
```

This pattern is realized by occupying just the two sites (0,1),(1,1)
above the marked vacant edge. It is the shortest possible bypass.

## 3. Multiple bypass components: a realizable negative contact

There is no universal nonnegative-contact theorem. Consider the actual
adjacent partition

```
labels = 01020201,
B = {xN,xS,yN,yS},   e={xE,yW},   C={xW,yE}.                 (6)
```

It has two outside bypass components B,C in addition to e. Its exact
colour coarsenings give a short independent derivation. With all three
colours distinct, each local pattern has one opposite pair, so
`(4K_x)(4K_y)=36`; the three-colour falling-factorial derivative is -1.
Merging e and C gives two opposite pairs at each vertex and product4;
the two-colour derivative is +1. Either other merger gives a 3+1
pattern and zero. Thus

```
g16=-36+4=-32,             g_pi=-2.                         (7)
```

This is realizable on an ordinary rectangular torus with Lx,Ly>=5:
leave x=(0,0),y=(1,0) vacant and occupy exactly

```
B={0,1} x {1,...,Ly-1},
C={2,...,Lx-1} x {0}.
```

B and C are individually NN-connected and have no NN edge between
them. B meets the four north/south ports; C meets the west/east ports;
e remains isolated. Both occupied components have rank0, since each
is missing the closing segment of its potential periodic direction.
Their open paths nevertheless connect the required ports around the
torus. Thus (7) is a genuine contact geometry, not an inadmissible
abstract partition or an artefact of giving the contact edge two IDs.

For comparison, splitting the C connection into two private outside
components leaves both one-site partitions unchanged and gives
`01020301`: now there is only one bypass B and d_x=d_y=-1, hence g=1.
The same two one-site connectivity patterns can therefore yield
different signs depending on their cross-site gluing. The contact
kernel is not a function of the two marginal one-site partitions.

Only these three specified rows were read from the existing lookup:
`01230451 -> g16=4`, `01020301 ->16`, and `01020201 ->-32`.
The latter two values also follow from (3) and the coarsening proof
above; no table or configuration population was regenerated.

## 4. A homogeneous contact source is not a thermal/source alias

Let

```
S_adj(A) = sum_(unordered nearest-neighbour {x,y}, x,y vacant) g_(pi_xy(A)).
```

On a rectangular torus with side lengths at least5, compare two
rank0 occupations with exactly K=2:

- Two isolated occupied sites: every occupied component has one site.
  Adjacent vertices have no common NN neighbour, so no adjacent vacant
  pair has an outside bypass. Equation (2) gives S_adj=0.
- One occupied domino: there is one occupied component of two adjacent
  sites. Exactly two adjacent vacant pairs meet this component at both
  endpoints: the opposite edges of the two incident plaquettes. Each
  has the clean bypass (5), and every other adjacent pair has zero.
  Therefore S_adj=1/2.

Both occupations have the same K, r=0, q=-1, E=1, and the same original
closed source Sstar. Indeed both occupied graphs are forests, so
`Sstar=2N+1-3K=2N-5`. They also have the same homogeneous *one-site*
regular-pair activation: each vacant site sees four different outside
components. A site cannot be adjacent to both endpoints of an occupied
domino on the triangle-free square graph. Thus the one-site activation
is `S_one=N-2` in both occupations.

Nevertheless their S_adj values differ. This rules out a pointwise
reduction of the homogeneous contact source to any function of these
named one-site/source/thermal coordinates. In particular it is not a
common Bernoulli-clock score proportional to K plus a normalization.
This is stronger than observing that an abstract contact tensor is
nonzero. It does **not** preclude an accidental zero after projection
to the particular pooled-root slope U at a particular size.

For the prescribed homogeneous local amplitude epsilon/N,

```
partial_logQ partial_epsilon^2 (partition ratio)|_(Q1,epsilon0)
   has occupation insertion (2/N^2) sum_(unordered vacant {x,y}) g_pi.
```

The adjacent part consequently uses `H_adj=2 S_adj/N^2`; the domino
has H_adj=1/N^2 and the isolated pair has zero. Equivalently an ordered
pair sum uses coefficient 1/N^2. This factor2 comes from differentiating
the finite product of local vertex factors, not from the physical edge
colour normalization. Original q/E numerators insert the same H_adj
against the fixed original observer.

## Consequence for the planned J2 split

Adjacent terms are valid entries of the existing lookup and can be
nonzero, signed, and independent of the old one-site activation. They
measure a contact edge **together with outside bypass organization**.
They are neither automatically zero nor automatically removable by a
thermal reparametrization. Conversely these local proofs do not claim
that the measured total J2 is predominantly nonadjacent or that its U
projection is nonzero. The declared adjacent/nonadjacent split is
therefore a substantive mechanism discriminator. No conclusion about
non-contact global transmission should be drawn from total J2 alone.
