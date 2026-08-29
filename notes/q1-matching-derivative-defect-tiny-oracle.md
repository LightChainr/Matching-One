# A finite-graph Q=1 matching derivative defect

## Exact object

For a declared pairing of a primal configuration `A` with a complement
configuration, define the integer cluster defect

```text
d(A) = k_G(A) - k_H(A^c) - ell(A)
```

and the normalized Laurent amplitude

```text
Xi_G,H(Q,p) = E_p[Q^d(A)].
```

This is finite and exact even when no positive local generic-Q model is known.
It obeys

```text
Xi(1,p) = 1,
(Q d/dQ) log Xi | Q=1 = E_p[d],
(Q d/dQ)^2 log Xi | Q=1 = Var_p(d).
```

Thus the matching cluster observable is already a well-defined Q tangent.  A
local Potts/interface realization would be an additional theorem: it must
factor this Laurent amplitude through a bounded local state space and satisfy
pull-through/fusion.  The construction does not assume that theorem.

## Tiny exact outcomes

For the planar edge-FK control `C3` and its two-vertex, three-parallel-edge
dual, Euler gives

```text
k_G(A) - k_Gstar(A^c) = |V_G| - |A| - 1 = 2 - |A|.
```

Consequently `d(A)=0` for every one of the eight configurations, not merely
after averaging.  Replacing the dual by the nondual three-edge path `P4`
while retaining the same local subtraction gives

```text
E_p[d] = p^3 - 3 p.
```

This is a minimal obstruction certificate: equal configuration counts and a
complement bijection do not produce a duality defect.

For the site-matching toy pair `C4 -> K4`, with occupied clusters on `C4` and
vacant clusters on `K4`, the raw derivative is

```text
E_p[k_C4(A)-k_K4(A^c)] = 2 p^4 - 4 p^2 + 4 p - 1,
```

and equals `1/8` at `p=1/2`.  This is a finite matching tangent, but it is not
by itself evidence for a local generic-Q seam.

## Why S/D can be only a projection

Given two scalar responses `(R_G,R_H)`, the combinations `S=R_G+R_H` and
`D=R_G-R_H` are always the two covectors associated with the formal exchange
matrix

```text
J = [[0,1],[1,0]],   P_plus/minus = (I plus/minus J)/2.
```

Calling them parity *fields* requires more: a physical two-way identification
between the two state spaces whose composite is the identity (or a declared
projector).  A one-way matching augmentation such as `C4 -> K4` supplies no
such inverse.  In that case S/D are projections of a doubled observable
vector; they do not diagonalize an intrinsic symmetry of either model.

The next exact step for Issue #233 is therefore sharply separated: search for
a bounded-locality factorization of `Xi` and compute its fusion tangent.  The
tiny oracle fixes what that construction must reproduce and supplies a cheap
negative control.

## Claim boundary

- Exact: the finite sums, the two Q-derivative identities, and all tiny-graph
  polynomials in the committed oracle.
- Interpretation: `Xi` is a candidate scalar matrix element of a derivative
  interface.
- Conjecture: the square-site matching observable admits a bounded-locality
  generic-Q defect whose tangent realizes this object.
