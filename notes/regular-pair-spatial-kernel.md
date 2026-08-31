# Canonical regular-pair spatial kernel: two shared components are necessary

The canonical local completion `Kreg=average i(I-P1)i^dagger` has an
exact two-site Q1 kernel on all 4140 external eight-port connectivities.
It vanishes whenever at most one outside component meets both vertices.
Two shared components suffice, but their contribution is signed: a
simple realizable square-lattice pattern gives -1/8. With exactly two
shared components the kernel factors into two named local colour
contrasts, explaining both the selection rule and the negative sign.

This supplies the production lookup requested from base `2ba8863f`.
The coefficient of K0 remains identically one. No alternate completion,
N25 population, previous double-insertion norm, or Monte Carlo result
was evaluated.

## 1. Observable and exact integer construction

Order the ports as

```
(xN,xE,xS,xW,yN,yE,yS,yW).
```

Let pi be their exact outside **graph-component partition**, with b
blocks. For two distinct vacant vertices define

```
B_pi(Q) = Q^(-b) sum_colours Kreg_x Kreg_y D_pi,
g_pi = partial_Q B_pi(Q)|_(Q=1).                            (1)
```

The original occupation and any of its rank labels remain fixed.
Additional outside components not meeting the eight ports cancel from
the ratio. Shared components include all hypergraph components: if the
two vacant vertices are adjacent, their common isolated edge-port also
counts. One must not count only occupied-site clusters.

For each exact colour-equality partition sigma coarser than pi, let
k=|sigma|. It has `(Q)_k=Q(Q-1)...(Q-k+1)` distinct colour assignments.
The finite-network theorem gives B_pi(1)=0. More explicitly, the k=1
term has Kreg_x=Kreg_y=0 identically, while every k>=2 has `(1)_k=0`.
Consequently neither derivatives of the entries nor of Q^(-b) survive
in the first derivative. The only surviving derivative is

```
partial_Q (Q)_k|1 = (-1)^(k-2)(k-2)!,       k>=2.            (2)
```

For an exact local colour-equality pattern rho, put
`L(rho)=4 Kreg(rho;1)`. Direct substitution into the regular entries
gives the small local alphabet

| Local colour equality | L |
|---|---:|
| four distinct | -8 |
| one opposite pair NS or EW | -6 |
| one adjacent pair | -3 |
| two opposite pairs NS\|EW | -2 |
| two adjacent pairs | -1 |
| 3+1 or all4 | 0 |

This table concerns **colour equality after coarsening**, not the
external connectivity pi before the colour sum. The integer production
kernel is exactly

```
g16(pi) = 16 g_pi
  = sum_(sigma coarser than pi, |sigma|>=2)
      (-1)^(|sigma|-2)(|sigma|-2)! L(sigma_x)L(sigma_y).       (3)
```

The script implements (3) directly with integer arithmetic. It does not
fit a polynomial or evaluate a sequence of numerical Q values.

## 2. A shared-component selection theorem

Let s be the number of blocks of pi meeting both vertices. At integer
Q>=4, condition on the colours of these shared blocks and sum the
colours of all private blocks at each vertex. The two remaining local
factors are independent conditional on those colours.

If s=0, the full normalized closure factors as

```
B_pi(Q)=B_(pi_x)(Q) B_(pi_y)(Q).
```

If s=1, global colour invariance makes each local conditional factor
independent of the one specified shared colour, so precisely the same
normalized factorization holds. Both one-vertex factors vanish at Q1.
Differentiating proves the exact rule

```
s<=1  =>  g_pi=0.                                           (4)
```

The argument is an identity in the rational completion and does not
assume an independent physical occupation law between x and y. It is
a statement about the colour contraction for a *fixed* outside graph.

## 3. Exactly two shared components: a local contrast product

Call the two shared colours a,b. After averaging all private component
colours, write the local conditional factors as F_x(a,b;Q), F_y(a,b;Q).
Each has just two possible values, equal or unequal shared colours.
Define

```
d_x(Q)=F_x(equal;Q)-F_x(unequal;Q),
d_y(Q)=F_y(equal;Q)-F_y(unequal;Q).
```

The indicator that two independent uniform shared colours coincide
has variance `(Q-1)/Q^2`. Therefore the exact conditional decomposition is

```
B_pi(Q)=B_(pi_x)(Q)B_(pi_y)(Q)+(Q-1)d_x(Q)d_y(Q)/Q^2,
g_pi=d_x(1)d_y(1),                  s=2.                    (5)
```

This explains why a positive-looking product of local kernel entries
does not imply a positive Q-derivative. The two local contrasts can
have different signs. Their values can be compressed into the following
contact rules. Put kappa=2 for an opposite repeated pair and kappa=1
for an adjacent repeated pair.

| Four ports at one vertex, with exactly two shared components | d(1) |
|---|---:|
| shared components occupy 2+2 ports | kappa/4 |
| one shared component repeats; other shared and one private each once | -kappa/2 |
| one private component repeats; both shared occur once | -3 kappa/4 |
| four different components; shared ports adjacent | -1/2 |
| four different components; shared ports opposite | 0 |
| a component occupies 3 or 4 ports | 0 |

These follow by summing the remaining private colour(s). For example,
in the second row, fixing unequal shared colours leaves one private
colour. Matching it to the repeated shared colour gives zero; matching
it to the other shared colour gives `-kappa/4`; each other colour gives
`-3kappa/4`. Continuing their count Q-2 to Q1 yields the unequal factor
`kappa/2`, hence the contrast `-kappa/2`. The equal-shared factor is zero
at Q1. The other rows follow by the same finite colour sums.

## 4. Minimum support and a physically realizable negative example

The smallest nonzero shared count is two. A two-block external pattern

```
labels = 00110011,
blocks = {xN,xE,yN,yE}, {xS,xW,yS,yW}
```

has adjacent 2+2 contact at each vertex. Equation (5) gives
`g=(1/4)(1/4)=1/16`, so its stored integer is g16=1.

Adding one private outside component changes the sign:

```
labels = 00110012,
blocks = {xN,xE,yN,yE}, {xS,xW,yS}, {yW},
g=(1/4)(-1/2)=-1/8,                 g16=-2.                 (6)
```

This is not merely an abstract Bell partition. An explicit square-lattice
realization uses vacant x=(0,0), y=(6,0), on an honest torus large enough
that the following finite sets do not meet their periodic copies:

```
A = {(i,2): 0<=i<=7}
    union {(0,1),(1,0),(1,1),(6,1),(7,0),(7,1)},
B = {(i,-2): -1<=i<=6}
    union {(-1,0),(-1,-1),(0,-1),(6,-1)}.
```

Occupy exactly A union B. Each set is NN-connected and the two have no
NN adjacency. A carries xN,xE,yN,yE; B carries xS,xW,yS. Both endpoints
of the y-west edge are vacant, so yW is its own isolated hypergraph
component. This realizes (6), with both occupied components contractible.
A translation places the picture within, for example, a 20-by-12 torus.
Thus **even realizable ordinary square-lattice partitions need not give
nonnegative spatial kernel**. No continuum or topology-rank attribution
is needed for this sign statement.

## 5. Complete finite result and production format

One exact pass over Bell8 gives:

| Shared components s | Partitions | Nonzero | Positive | Negative |
|---:|---:|---:|---:|---:|
| 0 | 225 | 0 | 0 | 0 |
| 1 | 1369 | 0 | 0 | 0 |
| 2 | 1922 | 1250 | 986 | 264 |
| 3 | 600 | 600 | 384 | 216 |
| 4 | 24 | 24 | 24 | 0 |
| total | 4140 | 1874 | 1394 | 480 |

These are counts over the complete abstract outside-partition space,
not frequencies in a lattice ensemble. The nonnegative s=4 class does
not imply nonnegativity of the full kernel. No assertion that every
listed partition is realizable in a specified torus is made.

The production file is
[`analysis/regular_pair_spatial_kernel.tsv`](../analysis/regular_pair_spatial_kernel.tsv),
with two integer columns, `key` and `g16`, and 1874 data rows. Exactly
zero rows are omitted. **An absent valid canonical key means zero.**

Its SHA256 is
`36ae069d370b1d7a4398861c928afb41aa76885c8895c696b1bc0c97e9c314fd`.

Canonicalization is first-appearance restricted-growth labelling in
the fixed port order of Section 1: the first component receives 0,
the next previously unseen component 1, and so on. Encode it as

```
key = sum_(i=0..7) (label[i] << (3*i)),
g = g16/16.
```

Labels range from 0 to 7 and the key fits in 24 bits. Raw union-find
root IDs must be canonicalized before packing. Noncanonical keys are
outside the interface, not evidence that the physical kernel is zero.
The two sites must each use N,E,S,W in the stated common orientation.

Regenerate with the managed research Python:

```
/Users/lc/python-envs/research-py311/bin/python \
  scripts/regular_pair_spatial_kernel.py
```

The script prints the exact shared-component summary and writes only
this lookup. It retains all signed entries. A caller evaluating the
named local vacant-vertex intervention must separately impose that
both sites are vacant; that occupation indicator is not in the key.

## 6. Connection to a connected Q response

For this regular completion, each one-insertion and two-insertion
coefficient is zero at Q1. Hence the Q derivative of the disconnected
product of one-insertion coefficients is also zero. At this derivative
order the connected two-site insertion has the same g_pi in (1);
subtracting the product of their *activation derivatives* would instead
mix in a second-order-in-Q quantity. Background occupation derivatives
also multiply a zero two-insertion coefficient at Q1. This states the
kernel's order and normalization, without defining an additional source
or carrying out a new ensemble score.

The scientific content is the shared-component selection rule and its
signed two-component contrast mechanism. The lookup supplies the exact
canonical local interaction, not a generic colour catalogue. It does
not establish a continuum field, a universal positive correlator, or
an independent evidence block beyond the declared occupation ensemble.
