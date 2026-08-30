# Unrestricted degenerate-quotient digital Alexander theorem

Status: theorem for every finite-index period sublattice of `Z^2`.

## Statement

For every nonsingular integer period matrix `P` and every periodic black-site configuration,

```text
im_Q H1(white matching) = im_Q H1(black NN)^perp,
r_black + r_white = 2,
q = r_black - 1 = 1 - r_white = (r_black-r_white)/2.
```

No honest quotient-cell hypothesis is needed. Every connected integral carrier image is
saturated: rank zero gives 0, rank one gives one primitive line, and rank two gives all Z2.

## Why self-identifying faces are harmless

Replace `L` by `2L`. The four-sheeted regular cover `T_(2L) -> T_L` has honest square faces:
vectors in `2L` have even ambient coordinates, whereas every nonzero difference between two
unit-square corners has a coordinate of absolute value one. The existing honest-cell theorem
therefore applies upstairs even if the base presentation has loops, repeated edges or identified corners.
Every removed diagonal has the same lifted displacement as its white NN replacement; the difference
therefore remains ambient-null after any quotient projection.
Rational graph-image homology descends exactly because a downstairs loop's finite deck monodromy
is killed by a positive iterate. In period bases `p_*=2I`, so the intersection form scales by four
and orthogonal complements descend with it.

## Proof chain

1. **canonical finite honest cover.** For every L=P Z^2, the sublattice L'=2L defines a four-sheeted regular cover T_(2L)->T_L. Every vector in 2L has even ambient coordinates, so no nonzero difference of unit-square corners lies in 2L; all quotient faces upstairs have four distinct corners.
2. **honest-cell theorem upstairs.** Lift the coloring and both graphs to T_(2L). The existing 16-pattern pruning and complementary-subsurface proof applies without change, giving I'_W=(I'_B)^perp and r'_black+r'_white=2.
3. **lifted-chain replacement.** A removed white diagonal and its same-face NN replacement have the same relative boundary and identical lifted displacement. Their difference has ambient H1 class zero after every quotient projection, including loops and repeated edges.
4. **complementary-subsurface duality.** For complementary subsurfaces U,V of an oriented torus, exact sequence, excision and Poincare-Lefschetz duality give im H1(V;Q)=im H1(U;Q)^perp.
5. **finite-cover rational image descent.** For the full inverse-image graph, every upstairs cycle projects downstairs. Conversely, a downstairs loop has finite deck monodromy, so a positive iterate lifts closed; after tensoring with Q its class is in the projected upstairs image. The same argument componentwise preserves the repository maximum-component ranks.
6. **symplectic rational descent.** In compatible period bases p_*=2I, hence omega(p_*u,p_*v)=4omega(u,v). It carries rational orthogonal complements to orthogonal complements.
7. **honest integral carrier classification.** On every honest qL cover, black regular neighborhoods and the integrally pruned white complementary carriers have the same Z-H1 images as their graphs. A connected torus subsurface has image 0, one primitive line, or all H1(T2;Z), according to rank 0,1,2.
8. **exact component stabilizer on qL covers.** For a universal-lift component with stabilizer H<=L, a loop in its qL-cover component has endpoint displacement in H intersect qL; connectedness gives a loop for every such displacement. Thus the integral upstairs image is exactly H intersect qL.
9. **coprime-prime Smith descent.** If S=Sat_L(H) and d=[S:H]>1, choose a prime q not dividing d. In a basis of the direct summand S extended to L, Smith coordinates give H=sum d_i Z e_i, H intersect qL=qH, and index d in qS. This contradicts honest-carrier saturation upstairs, so d=1.
10. **filtration corollaries.** Monotonicity plus rank complementarity gives both births, full two-threshold reconstruction and swapped/reversed reflection. In the rank-one case a line equals its symplectic orthogonal in Q2; nested one-dimensional images keep that rational line and its canonical primitive representative constant.

## Machine certificates

- universal face patterns: 16; all pass
- projected HNF representatives: 86 (35 self-identifying)
- canonical four-sheeted honest covers: 86; failures: zero
- consistent projected face patterns: 1020
- cached subset states: 31068; rank-one states: 17248
- rank-sum, rank-mark, primitive-line and projection failures: zero
- nonsaturated rank-one states in the finite regression: 0 (theorem-predicted regression)
- symbolic threshold pairs: 364; failures: zero
- all machine gates: `True`

The finite HNF layer is a regression certificate, not the unrestricted inference. Rational
duality uses the four-sheeted cover; integral saturation uses all honest qL covers and a
coprime-prime Smith contradiction.

## Boundary

- The theorem concerns ambient H1 images of the NN and complementary matching graphs, not their graph cyclomatic numbers.
- Every connected black NN or white matching carrier has saturated integral ambient-H1 image; at rank one the actual graph-image subgroup is generated by the primitive direction.
- The index-2-through-13 saturation census is now regression and implementation-convention audit only.
- No continuum, threshold, CFT-field or finite-size-scaling consequence is claimed.
- The rational theorem uses the four-sheeted cover; the integral theorem uses the family of qL covers and a coprime-prime Smith contradiction.
