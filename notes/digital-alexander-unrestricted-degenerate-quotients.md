# Unrestricted digital Alexander duality on finite square quotients

Status: unrestricted finite-index theorem for Issue 269, with the integral
saturation refinement separated from the proved rational statement.

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

When both ranks are one, black and white have the same rational winding line.
That line has a unique primitive integral representative up to sign, hence a
canonical repository direction after the fixed sign convention. This theorem
does **not** assert that the actual integral graph-image subgroup has
saturation index one.

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

## 6. Exact boundary of the integral claim

Finite-cover descent is a rational statement. It cannot by itself preserve
integral saturation because `p_*=2I`. For a rank-one component with integral
cycle windings `v_i=n_i ell`, where `ell` is the canonical primitive direction,

```text
im_Z H1(C)=gcd_i(|n_i|) Z ell.
```

Thus saturation index one is equivalent to the exact additional condition
`gcd_i(|n_i|)=1`. This condition is necessary and sufficient and is directly
executable with the existing `_rank_line_safe` classifier. It is not required
for rank complementarity, the Q-score, threshold reconstruction, reflection,
or the primitive projective-line gate.

No nonsaturated rank-one state occurs in the exact HNF regression through
quotient index 10, but that remains a finite diagnostic rather than an
unrestricted theorem. The counterexample locator is the existing subset-DP:
its first state with returned saturation index greater than one would be a
minimal witness, while leaving all rational conclusions intact.

## Machine certificate and scientific boundary

The executable certificate verifies:

1. all 16 lifted local face-chain identities;
2. the canonical degree-four honest cover, computed for all 86 existing HNF
   regression representatives;
3. projection through the 35 self-identifying representatives, exercising
   loops and repeated endpoints;
4. all 31,068 cached subset states through index 10, including 17,248
   rank-one states;
5. the two-dimensional symplectic-line identity and all two-threshold paths
   through 12 sites.

The finite checks are regression and counterexample localization, not the
logical basis of the unrestricted theorem. The theorem is the finite
honest-cover argument plus rational image and conformally symplectic descent.
It concerns ambient homology images, not graph cyclomatic numbers, and makes
no continuum, critical-threshold, finite-size-scaling or CFT-field claim.
