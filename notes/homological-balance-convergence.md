# Homological balance roots converge to the site threshold

Status: theorem reduction for honest `L x L` periodic square-cell tori.  The
only infinite-volume inputs are standard subcritical exponential decay on the
square and matching lattices and the planar matching-threshold relation.

## Statement

Let `G_L` be the square nearest-neighbour site torus and let `Ghat_L` be its
nearest-plus-next-nearest matching torus.  For a black configuration at density
`p`, put

```text
r_B = rank im[H1(black NN complex) -> H1(T^2; Q)].
```

For the white matching configuration at density `1-p`, define `r_W`
analogously.  On every honest periodic square-cell torus, digital Alexander
duality gives

```text
r_B + r_W = 2.
```

Consequently the finite matching function is

```text
M_L(p) = E_p[r_B] - 1,
```

and its unique zero `p_L` is the homological balance point `E[r_B]=1`.
If

```text
p_c(G) + p_c(Ghat) = 1,
```

then

```text
p_L -> p_c(G).
```

This conclusion uses no finite-size scaling ansatz, CFT field assignment, or
root-shift exponent.

## 1. Strict finite-volume monotonicity

The ambient-image rank is monotone under site insertion.  Russo's formula for
the increasing integer-valued observable gives

```text
M_L'(p) = sum_v E_p[r(omega^{v=1}) - r(omega^{v=0})] >= 0.
```

It is strictly positive for `0<p<1`.  Fix a horizontal primitive row and a
vertex `v` on it.  Occupy every other vertex of that row and no off-row
vertex.  Without `v` the occupied graph is a contractible path; inserting `v`
closes one primitive torus cycle, so `Delta_v r=1`.  The conditioning event
has positive probability at every interior `p`.  Hence `M_L` is strictly
increasing.  Since `M_L(0)=-1` and `M_L(1)=1`, its zero exists and is unique.

## 2. One subcritical lemma handles both sides

For a fixed finite-range planar site lattice `H`, standard subcritical
sharpness implies that for every `u<p_c(H)` there are positive constants
`c(u),C(u)` such that

```text
P_u(0 connected to distance n) <= C exp(-c n).
```

If the occupied subgraph of an `L x L` torus has nonzero ambient `H1`, a lift
of a noncontractible cycle contains an occupied connection through distance at
least `L/2`.  A union bound over the `L^2` vertices therefore gives

```text
P_u(r_H > 0) <= C' L^2 exp(-c' L) -> 0.                 (1)
```

No supercritical giant-component theorem is needed.

### Below the square-site threshold

Fix `epsilon>0` and set `p=p_c(G)-epsilon`.  Apply (1) directly to the black
NN lattice:

```text
P_p(r_B=0) -> 1,
M_L(p) -> -1.
```

More quantitatively,

```text
0 <= M_L(p)+1 = E_p[r_B] <= 2 C L^2 exp(-cL).
```

### Above the square-site threshold

Set `p=p_c(G)+epsilon`.  The white matching density is

```text
1-p = p_c(Ghat)-epsilon,
```

so (1) applies to the white matching graph.  Thus `r_W=0` with probability
tending to one.  The exact finite identity `r_B+r_W=2` now gives

```text
P_p(r_B=2) -> 1,
M_L(p) -> +1,
```

and

```text
0 <= 1-M_L(p) = E_p[r_W] <= 2 C_hat L^2 exp(-c_hat L).
```

The upper half of the proof is therefore just the lower half applied to the
matching complement.  In particular, it does not require separately proving
that a supercritical black cluster wraps in two independent directions.

## 3. Root convergence

For every fixed `epsilon>0`, the preceding bounds imply, for all sufficiently
large `L`,

```text
M_L(p_c-epsilon) < 0 < M_L(p_c+epsilon).
```

Strict monotonicity then yields

```text
p_c-epsilon < p_L < p_c+epsilon.
```

Since `epsilon` was arbitrary, `p_L -> p_c`.

The argument also proves exponential saturation of the homology rank at every
fixed distance from criticality.  It does **not** give a convergence rate for
`p_L-p_c`: that requires near-critical, rather than fixed-subcritical,
information.

## 4. Scope and scientific card

- Mechanism changed: the matching root is a consistent topological
  phenomenological-renormalization coordinate, not merely an empirically good
  polynomial zero.
- Exact layer: `M_L=E[r_B]-1`, `r_B+r_W=2`, finite-volume monotonicity.
- Imported probability layer: subcritical exponential decay on each of the
  two finite-range planar lattices and `p_c(Ghat)=1-p_c(G)`.
- Observer/sector/source/geometry: ambient `H1` rank / topology / site
  occupation / honest square-cell torus.
- Dependency group: Issues #269, #275, #276 and the digital-Alexander proof.
- Not proved: a root-shift exponent, a CFT operator identity, optimal estimator
  variance, or coverage of self-identifying short-period quotient cells.
- Next promotion observation: derive a near-critical two-sided bound for
  `P(r_B>0)` and `P(r_W>0)` at `p_c +/- lambda L^(-3/4)`; this would turn
  consistency into a scaling-window theorem.

## Literature anchors

- H. Duminil-Copin and V. Tassion, *A new proof of the sharpness of the phase
  transition for Bernoulli percolation and the Ising model*, arXiv:1502.03050.
- P. Duncan, M. Kahle and B. Schweinhart, *Homological percolation on a torus:
  plaquettes and permutohedra*, arXiv:2011.11903.

The first supplies the subcritical exponential-decay architecture.  The
second supplies the closest ambient-image homology analogue; the present
4/8-site argument is shorter because exact digital Alexander duality converts
the supercritical side into a second subcritical estimate.
