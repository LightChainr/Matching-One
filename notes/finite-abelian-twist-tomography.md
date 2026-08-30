# Finite-abelian twists close the unmarked source and open line tomography

## Result

Let `Lambda` be the global ambient integral `H1` image of a black NN or white
matching configuration. The unrestricted theorem in #269 says that every
connected carrier image is saturated in `Z^2`. This globalizes: disjoint
essential components of an embedded torus carrier must be parallel, so a
rank-one union still has one primitive image line; if any component has rank
two, its image is already all of `Z^2`. Thus the global `Lambda` is saturated
as well. Hence `Z^2/Lambda` is free of rank `2-r`, where `r=rank Lambda`.

For any finite abelian group `A` of order `n`, define a flat twist

```text
alpha in Hom(Z^2,A)
```

and let `T_alpha(p)` be the probability that `alpha` vanishes on the occupied
ambient image. Configurationwise,

```text
#{alpha : alpha|Lambda=0}
  = |Hom(Z^2/Lambda,A)|
  = n^(2-r).
```

Therefore the sum over every flat `A`-twist is

```text
S_n(p) := sum_alpha T_alpha(p)
        = n^2 P0(p) + n P1(p) + P2(p),
```

and its normalized average is exactly the intrinsic source at a discrete
negative value:

```text
S_n/n = n P0 + P1 + n^-1 P2
      = Z_top(s=-log n).
```

The aggregate depends only on `|A|`, not on the decomposition of `A`.

## Two discrete source values are complete

The unmarked source has only the three probabilities `(P0,P1,P2)` with one
normalization. Consequently the two twist traces `S_2,S_3` already reconstruct
it exactly:

```text
P0 = (S3 - 2 S2 + 1)/2,
P1 = S2 - 1 - 3 P0,
P2 = 1 - P0 - P1.
```

Thus Phase E of #337 does not need an infinite source scan. Order-2 and
order-3 flat-twist averages are an exact finite tomography of the whole
unmarked rank-source functional.

## What individual prime twists see

For `A=F_q` and nonzero `alpha`, `ker(alpha)` is a projective line in
`P^1(F_q)`. A saturated rank-one image has a primitive integral generator
`ell`, whose reduction modulo `q` is nonzero. It follows that

```text
T_alpha = P0 + L_ker(alpha),
```

where `L_line` is the probability that the carrier has rank one and its
primitive winding direction reduces to `line` modulo `q`. Every projective
line is the kernel of exactly `q-1` nonzero twists.

If `P0` is available from the ordinary rank archive, each nonzero twist
therefore returns one modular projective-line bin by subtraction. This is an
exact bridge to #334: flat twists are a Fourier/twist implementation of the
projective line mark, coarsened modulo `q`.

Several primes can be combined by CRT when an external bound on the primitive
line is available. A finite prime list alone does not determine an unbounded
integer line.

## Correction to the proposed arithmetic interpretation

The finite-field rank cannot supply saturation tomography on the
theorem-supported carriers. Saturation gives

```text
r_q = r
```

for every prime. There is no hidden `q`-primary Smith state to find. The new
information lives instead in the individual twist constraint sectors, which
resolve the already-primitive line modulo `q`.

This strengthens rather than weakens the intrinsic-source program:

```text
aggregate twists at |A|=2,3  -> complete unmarked source,
individual F_q twists         -> modular projective-line tomography.
```

It also identifies a practical continuum target. One should compare the
finite set of twist-sector amplitudes, not a hypothetical finite-field rank
defect, with affine-TL or torus defect partition functions.

## Boundary

This is an exact finite-volume cohomology transform. It does not yet identify
the twist constraints with a local CFT field or prove a particular affine-TL
module decomposition.
