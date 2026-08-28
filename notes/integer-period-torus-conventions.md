# Integer-period torus state and edge conventions (C00)

Status: implementation contract for the general 2x2 homology engine.
This is the single correctness path shared by later Gaussian, sheared and
modular tori.  It extends the PR #21 rank-0/1/2 classifier; it does not
replace it.

## Period matrix

A quotient of `Z^2` is given by a 2x2 integer matrix `P` with `det(P) != 0`.
Columns are the two declared generators:

```text
P = [[P00, P01],
     [P10, P11]]
v0 = (P00, P10)
v1 = (P01, P11)
N = |det(P)|
```

Special cases used as regressions:

| name | stored coordinates | P | N |
|---|---|---|---|
| axis `L` | `(x, y)` | `diag(L, L)` | `L^2` |
| diamond `L` | `(u, v) = (x+y, y-x)` | `diag(2L, 2L)` | `2 L^2` |
| diamond-xy `L` | `(x, y)` | columns `(L, L)`, `(-L, L)` | `2 L^2` |
| Gaussian `(a, b)` | `(x, y)` | columns `(a, b)`, `(-b, a)` | `a^2+b^2` |

The axis and uv-diamond embeddings are exactly the PR #21 geometries.
`diamond-xy` is the same lattice in original coordinates; exhaustive
rank and generator-relative counts match diamond `L=2`.

A 2-vector `(px, py)` is accepted as the diagonal matrix `diag(px, py)`.

## Exact winding map

A closed cover displacement `d = (dx, dy)` is a torus winding iff it lies
in the column lattice of `P`.  The generator coefficients are

```text
w = adj(P) d / det(P)
```

computed with integer adjugate, integer determinant, divisibility check,
and integer division.  Floating inversion is forbidden.  Equivalent form:

```text
P w = d
w0 = ( P11 dx - P01 dy ) / det(P)
w1 = (-P10 dx + P00 dy ) / det(P)
```

For Gaussian `P = [[a, -b], [b, a]]` this reproduces the circulant formula
`(a X + b Y)/N`, `(-b X + a Y)/N`.

## Sites and edges

Sites of a general integer quotient are the unique `Z^2` points of the
half-open parallelogram `{P s | 0 <= s0, s1 < 1}`.  There are `N = |det(P)|`
of them.  Ordering is `(y, x)` so `P = diag(L, L)` matches `axis_geometry(L)`.

Each stored edge keeps:

- `i`, `j`: indices of reduced representatives;
- `(dx, dy)`: the unreduced lattice step in the universal cover, never the
  reduced chord inside the parallelogram.

Primal steps are `+(1,0)` and `+(0,1)`.  Matching steps add `+(1,1)` and
`+(1,-1)`.  Destination indices are obtained by reducing `position + step`
through `adj(P)/det(P)` (Python floor division).

The uv-diamond embedding keeps its historical step set
`(1,-1)`, `(1,1)` and matching extras `(2,0)`, `(0,2)`.  Those displacements
are already in uv cover coordinates, with diagonal periods `2L`.

## Homology channels

Per occupied component the union-find retains a rational basis of at most
two primitive winding vectors.  Configuration channels:

| channel | meaning |
|---|---|
| rank 0 | no wrapping |
| rank 1 | one-dimensional or spiral wrapping |
| rank 2 / `cross` | two independent windings in one component |
| `direction_0` | some winding has nonzero first generator coefficient |
| `direction_1` | some winding has nonzero second generator coefficient |
| `either` | rank > 0 (Boolean wrap of PR #21) |
| `both` | both generator-relative flags |

`horizontal` / `vertical` remain aliases of `direction_0` / `direction_1`
for axis tori.

The stored basis is a Q-basis of primitive directions, enough for rank and
generator-relative flags.  It is not a Hermite normal form of the Z-module.

## Unimodular invariance

If `U` has `det(U) = +/-1` and `P' = P U`, then `P` and `P'` generate the
same lattice.  Windings transform as `w' = U^{-1} w`.  Invariant:

- component sizes and roots;
- homology rank;
- `either` and `cross`.

Not invariant, by design: `direction_0`, `direction_1`, and `both`.  Those
flags are generator-relative.  A shear that mixes the declared basis can
turn a purely first-generator wrap into a spiral in the new coordinates.

## PR #21 regressions that must remain exact

Black-primal exhaustive counts:

```text
axis L=3:    rank0=259, rank1=162, rank2=91, d0=175, d1=175
diamond L=2: rank0=143, rank1=68,  rank2=45, d0=81,  d1=81
```

`either` must continue to match the Boolean wrap detector in
`matched_torus_reference.cluster_stats`.

## Small Gaussian regression

Non-axis quotient `(a,b)=(2,1)`, `N=5`, `P = [[2,-1],[1,2]]`:

```text
rank0=16, rank1=10, rank2=6, d0=11, d1=11
```

Cyclic labels `j = a x + b y (mod N)` on the same sites match the
circulant geometry reference.
