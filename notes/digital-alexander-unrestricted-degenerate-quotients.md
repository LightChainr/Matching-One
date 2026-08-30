# Unrestricted digital Alexander duality on finite square quotients

Status: unrestricted finite-index theorem for Issue 269 over both rational
homology and the integral graph-image lattice.

## Theorem

Let `L=P Z^2` be any finite-index sublattice, with `P` an arbitrary
nonsingular integer `2x2` matrix. Let `B` be any `L`-periodic subset of
square-lattice sites. On `T_L=R^2/L`, let `G_B` be the black induced NN graph
and let `G_W^m` be the matching graph of the white complement, with NN and both
diagonal adjacencies. Loops, parallel edges and quotient presentations in
which unit-face corners identify are allowed.

Write `I_B` and `I_W` for their ambient images in `H_1(T_L;Q)`. Then

```text
I_W = I_B^perp,
rank(I_B)+rank(I_W)=2.
```

The maximum component ranks used by the repository are preserved by the
covering argument below. Therefore, configuration by configuration,

```text
r_black+r_white=2,
q=r_black-1=1-r_white=(r_black-r_white)/2.
```

Moreover, for every connected black or white matching component `C`, the
integral image

```text
H_C=im[H_1(C;Z)->H_1(T_L;Z)=L]
```

is saturated. It is zero at rank zero, `Z ell` for a primitive `ell in L` at
rank one, and all of `L` at rank two. Consequently every rank-one plateau has
actual saturation index one; `iota` is not an additional state coordinate.

## 1. Canonical finite regularization

Set `L'=2L`. Inclusion `2L < L` gives a four-sheeted regular cover

```text
p:T_(2L) -> T_L.
```

Every vector of `2L` has even coordinates in the ambient square-lattice
basis. Every nonzero difference of two corners of one unit square belongs to
`{-1,0,1}^2` and has at least one odd coordinate. Thus no two distinct
corners of a unit square agree modulo `2L`. Every square face in `T_(2L)` has
four distinct corners, for every nonsingular `P` and without any index bound.

This is the useful regularization: it is finite, canonical and stays inside
the exact lattice model. No limiting or generic-position argument is needed.

## 2. Apply the honest-cell theorem upstairs

Pull the coloring and both graphs back to `T_(2L)`. They are exactly the NN
and complementary matching graphs of the lifted coloring. Because all faces
upstairs are honest, the merged 16-pattern pruning and regular-neighbourhood
proof applies there. With primes denoting the lifted ambient images,

```text
I'_W=(I'_B)^perp,
r'_black+r'_white=2.
```

The local chain certificate remains useful here. In every lifted face, a
removed white diagonal and its white NN replacement have the same relative
boundary and exactly the same integer displacement. Their difference is a
closed lifted 1-chain of displacement zero. This proves that the pruning does
not alter ambient homology before or after either quotient, including regimes
that appear downstairs as loops, repeated edges or self-identified corners.

## 3. Rational graph images descend exactly

Let `C` be a connected component of either graph downstairs and let `C'` be a
component of its full inverse image. The restriction `C' -> C` is a finite
cover and is surjective. Projection sends every upstairs cycle to a
downstairs cycle, so

```text
p_* im_Q H1(C') subset im_Q H1(C).
```

For the reverse inclusion, take a loop `gamma` in `C`. A lift of `gamma` may
end at a different point in its finite deck fibre. The corresponding deck
monodromy has finite order, so some positive iterate `gamma^m` lifts to a
closed loop in `C'`. Therefore `m[gamma]` belongs to the projected upstairs
image. Division by `m` over `Q` proves equality.

Equivalently in the universal cover, if `H_C < L` is the displacement
stabilizer of a lifted component, then the stabilizer at level `2L` is
`H_C intersect 2L`, a finite-index subgroup of `H_C`. The two groups have the
same rational span and rank. This argument is componentwise, so it preserves
the repository's maximum-component ranks, not only the union-image rank.

Consequently

```text
p_* I'_B=I_B,
p_* I'_W=I_W
```

as rational subspaces.

## 4. Orthogonality descends

Use the columns of `P` as the period basis downstairs and the columns of
`2P` upstairs. In these bases

```text
p_*=2 I_2.
```

The oriented torus intersection forms consequently obey

```text
omega_L(p_*u,p_*v)=4 omega_(2L)(u,v).
```

Thus `p_*` is an invertible conformally symplectic map over `Q` and carries
orthogonal complements to orthogonal complements. Applying it to the
honest-cell identity upstairs gives

```text
I_W=p_*I'_W
   =p_*[(I'_B)^perp]
   =(p_*I'_B)^perp
   =I_B^perp.
```

This closes the unrestricted step. Self-identification is a defect of a
small quotient CW presentation, not a topological singularity of the torus.

## 4.1 Integral saturation requires all honest covers, not one

The single map `p_*=2I` is sufficient over `Q` but cannot by itself remember
an even integral index. The integral proof therefore uses the family

```text
T_(qL) -> T_L,  q>=2.
```

Every such cover is honest: vectors in `qL=qP Z^2` have both ambient
coordinates divisible by `q`, whereas a nonzero difference of unit-square
corners has a coordinate `+1` or `-1`.

First prove the carrier lemma upstairs. A connected embedded black NN graph
has a connected closed regular neighbourhood `U` which deformation retracts
to it. For white matching, each removed diagonal is replaced by an integral
same-face NN chain with the same boundary and lifted displacement. This
preserves connected components and the ambient image over `Z`; the pruned
component is the embedded one-skeleton of its complementary carrier.

For any connected compact subsurface `U subset T^2`, its integral ambient
image is saturated:

1. if `U` has genus one, it contains an intersection-one pair, whose images
   form a unimodular pair in `H_1(T^2;Z)`, so the image is all `Z^2`;
2. if `U` has genus zero, its first homology is generated by boundary curves;
   every essential embedded boundary circle is primitive, and disjoint
   essential circles on a torus are parallel, so the image is either zero or
   one primitive line.

It remains to descend this integral conclusion without using `p_*=qI` as an
inverse. Fix a component `C_tilde` of the universal periodic graph and define

```text
H={ell in L : C_tilde+ell=C_tilde}.
```

This stabilizer is exactly the downstairs integral ambient image: loops lift
to paths with endpoint displacement in `H`, and connectedness supplies a path
for every element of `H`. In the `qL` cover, the chosen lifted component has
stabilizer and integral image exactly

```text
H intersect qL.
```

Indeed, a loop upstairs has endpoint displacement in both groups; conversely,
every `h in H intersect qL` closes after projection to `T_(qL)`.

Suppose `H` were not saturated. Put `S=Sat_L(H)` and
`d=[S:H]>1`. Since `S` is primitive, it is a direct summand of `L`. Choose a
prime `q` not dividing `d` and a basis `e_1,...,e_r` of `S`, extended to `L`,
in which Smith normal form gives

```text
H = d_1 Z e_1 + ... + d_r Z e_r,
d = product_i d_i.
```

In these direct-summand coordinates, `qL intersect S=qS`, and coordinatewise

```text
d_i Z intersect q Z = lcm(d_i,q) Z = q d_i Z.
```

Therefore

```text
H intersect qL = qH,
[qS : H intersect qL] = product_i d_i = d > 1.
```

But the `qL` cover is honest, so its connected carrier image must be
saturated. This contradiction forces `d=1`. The argument covers rank one and
rank two; rank zero is saturated tautologically.

## 5. Birth, reconstruction and reflection

Along a site-addition filtration the black ambient image is nested and its
rank is nondecreasing. Define

```text
K1=min{k:r_black(k)>=1},
K2=min{k:r_black(k)=2}.
```

Rank complementarity makes the historical reverse matching birth equal to
`N-K1+1`. Therefore `K_minus=K1` and `K_plus=K2`, including a direct `0->2`
birth where the two thresholds coincide. Since ranks take only `0,1,2`,

```text
r_black(k)=1[k>=K_minus]+1[k>=K_plus].
```

Applying the same identity after swapping NN/matching and reversing the site
order gives both reflections:

```text
K_minus^G(pi)+K_plus^Ghat(reverse(pi))=N+1,
K_plus^G(pi)+K_minus^Ghat(reverse(pi))=N+1.
```

On a nonempty rank-one plateau, nested one-dimensional rational images cannot
change line. In a two-dimensional symplectic space a line is its own
orthogonal complement, so the complementary white line agrees. Its canonical
primitive representative is therefore constant as well.

## 6. Integral consequence for the persistent state

The existing `_rank_line_safe` gcd/SNF index remains useful as an
implementation and convention audit, but the theorem fixes it to

```text
iota=1
```

on every rank-one component and hence on every rank-one filtration state.
There is no independent arithmetic index process to infer. A returned index
greater than one is now a counterexample to the discrete carrier bridge or an
implementation error, not a new continuum state.

## Machine certificate and scientific boundary

The executable certificate verifies:

1. all 16 lifted local face-chain identities;
2. the canonical degree-four honest cover, computed for all 86 existing HNF
   regression representatives;
3. projection through the 35 self-identifying representatives, exercising
   loops and repeated endpoints;
4. all 31,068 cached subset states through index 10, including 17,248
   rank-one states, in the original local certificate;
5. all 16 integral white face-chain patterns and 56 rank-one/rank-two Smith
   defect rows in the integral certificate;
6. the merged index-2-through-13 frontier's 140 HNF representatives and
   101,140,028,118 paths, now classified as regression only;
7. the two-dimensional symplectic-line identity and all two-threshold paths
   through 12 sites.

The finite checks are regression and implementation diagnostics, not the
logical basis of the unrestricted theorem. The rational statement uses one
finite honest cover and conformally symplectic descent. The integral statement
uses honest carriers on all `qL` covers and the coprime-prime Smith
contradiction.
It concerns ambient homology images, not graph cyclomatic numbers, and makes
no continuum, critical-threshold, finite-size-scaling or CFT-field claim.
