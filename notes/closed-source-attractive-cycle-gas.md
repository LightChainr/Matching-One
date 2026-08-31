# The closed source is an attractive cycle gas for nonnegative coupling

**Finite-lattice result.** The already fixed source `S*=C+F+Bv` is
supermodular on occupied vertex subsets. Consequently its Bernoulli tilt
`exp(t S*)` satisfies the FKG lattice condition for every finite `t>=0`.
On an unrestricted nonalias square torus, this coupling range is sharp
for the lattice condition: two three-vertex paths around one unit face
give an exact failure at every `t<0`.

This note starts at `bc17b81d`. It uses the configurationwise action in
[the closed cluster-gas note](decimation-closed-cluster-gas-action.md),
with no numerical response, new enumeration or sampling.

## 1. The ambient correction exposes two supermodular dimensions

Fix the ordinary NN graph of an honest square-cell torus with N vertices.
For an occupied subset A, let G[A] be the vertex-induced occupied graph,
K(A)=|A|, and work over the rationals throughout. Write

```text
beta1(A) = dim H1(G[A];Q),
r(A) = rank[H1(G[A];Q) -> H1(torus;Q)],
beta_null(A) = dim ker[H1(G[A];Q) -> H1(torus;Q)].
```

These are graph cycles and their ambient image; filled elementary faces
are not attached to G[A] as two-cells. Rank-nullity gives
`beta_null=beta1-r`. The proven source identity, with `q=r-1`, is therefore

```text
S*(A) = 2 beta1(A)-3K(A)-q(A)+2N
      = 2 beta1(A)-r(A)-3K(A)+2N+1
      = beta1(A)+beta_null(A)-3K(A)+2N+1.              (1)
```

The constant is **2N+1**, not 2N. For example the empty configuration has
one vacant matching component and 2N vacant-vacant bonds, hence S*=2N+1.
The apparently adverse `-r` term has become the second nonnegative
cycle-space dimension in (1).

Choose one orientation for every edge of the full graph, so all cycle
spaces are subspaces of the same edge-chain space. Denote them by V_A,
and let Phi be the fixed ambient-homology map on the full graph cycle
space. Denote `W_A=V_A intersect ker Phi`. For I=A intersect B and
U=A union B, induced edge supports imply exactly

```text
V_I = V_A intersect V_B,       V_A+V_B subset V_U,
W_I = W_A intersect W_B,       W_A+W_B subset W_U.    (2)
```

The union graph may contain additional edges joining A\B to B\A; these
only enlarge the right-hand spaces. The kernel equalities use one common
ambient map, rather than separately chosen winding coordinates.

Apply the dimension formula to each row of (2). It gives

```text
beta1(U)+beta1(I) >= beta1(A)+beta1(B),
beta_null(U)+beta_null(I) >= beta_null(A)+beta_null(B).
```

Since K is modular, their sum proves

```text
S*(U)+S*(I)-S*(A)-S*(B) >= 0.                      (3)
```

More precisely the slack is the sum of two nonnegative integers,
`dim[V_U/(V_A+V_B)] + dim[W_U/(W_A+W_B)]`. This identifies the attractive
interaction: cycles available collectively across the two subsets, and
additional ambient-null cycles, supply the log-weight gain. No graph-class
catalog or separate proposed supermodularity of r is needed.

## 2. FKG law and its conditioning scope

For independent site parameters `0<p_v<1`, define the normalized weight

```text
mu_t(A) = Z^(-1) exp[t S*(A)]
          product_(v in A) p_v product_(v notin A) (1-p_v).
```

The logarithm of the Bernoulli factor is modular. Equation (3) implies

```text
mu_t(A union B) mu_t(A intersect B) >= mu_t(A) mu_t(B),   t>=0.  (4)
```

Normalization cancels from this inequality. The finite FKG theorem
therefore gives positive association:
`Cov_mu_t(f,g)>=0` for any two increasing real functions of the occupied
subset. This is a statement about the full interacting law, not only its
linearization at t=0.

The classical implication is Proposition1 of
[Fortuin, Kasteleyn and Ginibre, *Correlation inequalities on some partially
ordered sets*, CMP22,89--103(1971)](https://math.bme.hu/~balint/oktatas/perkolacio/percolation_papers/fortuin_kasteleyn_ginibre.pdf),
DOI10.1007/BF01651330. The cycle-space supermodularity above is the
model-specific input. The paper's lattice condition is sufficient, not
necessary, for association; this distinction matters for negative t below.

Conditioning a specified set of sites to be occupied and another specified
set to be vacant restricts the Boolean lattice to an interval. Its meet
and join stay in that interval, so (3)--(4) and positive association remain
valid. Thus hard checkerboard saturation and deterministic site boundary
conditions are included; probabilities 0 or 1 can be handled on their
remaining free-site support. Arbitrary conditioning on K, rank, a crossing
event or a sampled-prefix class need not produce such a lattice interval:
no conditional FKG claim is made for those ensembles.

The cycle-space argument itself preserves any genuine edge/gain incidence
convention. To identify (1) with the repository source C+F+Bv and use the
unit-face witness below, retain the honest nonalias square-cell convention.
An injective 3x3 stencil is a sufficient scope; for a Gaussian torus N>8
suffices. Silently suppressing quotient self-loops or distinct gain edges
is not an allowed change of that convention.

## 3. A four-site certificate makes the coupling range sharp

Take one contractible unit face, with corners cyclically a,b,c,d, and put
every site outside it in the vacant state. Use the two occupied sets

```text
A={a,b,c},     B={a,c,d},
I={a,c},      U={a,b,c,d}.
```

Each of A and B induces a three-site path, I consists of two isolated
sites, and U induces one contractible four-cycle. Thus

| Set | K | beta1 | beta_null | r | S* |
|---|---:|---:|---:|---:|---:|
| A or B | 3 | 0 | 0 | 0 | 2N-8 |
| I | 2 | 0 | 0 | 0 | 2N-5 |
| U | 4 | 1 | 1 | 0 | 2N-9 |

The supermodular slack is exactly 2. For every strictly positive Bernoulli
reference law, including unequal site parameters,

```text
[mu_t(U) mu_t(I)]/[mu_t(A) mu_t(B)] = exp(2t).        (5)
```

It is less than one for every t<0. Hence, for the unrestricted torus in
the stated scope, the FKG lattice condition holds **if and only if t>=0**.
Some restrictions may remove this obstruction: on an induced forest the
action is modular and all real t give a product law after changing site
activities. Failure of (4) at negative t is not a proof of negative
association, nor by itself a proof that a particular increasing pair has
negative covariance.

## 4. What attraction does and does not imply for the global observables

At fixed t>=0 and homogeneous p, r and K are increasing. Therefore

```text
partial_p E_mu[r] = Cov_mu(r,K)/[p(1-p)] >= 0,
partial_p E_mu[q] >= 0.                              (6)
```

This supplies a monotone finite matching-mean thermal curve. It does not
by itself establish a simple root. The sector observable E=q^2 is not
increasing, and neither the orientation projector nor its root/slope
normalization preserves the increasing-function class; (4) supplies no
sign prediction for their responses or for the global U.

Also, S* is not an increasing function of occupation: adding an isolated
occupied site changes it by -3. Attraction at a fixed nonnegative coupling
must not be replaced by an assertion of monotonic occupation or rank as t
increases. The result identifies an exact attractive microscopic cycle
gas for the already closed source, without a continuum-field identity or
an additional assumption about an RG flow.
