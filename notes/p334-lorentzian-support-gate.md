# Why the observed fixed-line ULC is not a standard Lorentzian consequence

The previous exact atlas left a tempting proof target: homogenize the fixed-line
rank-one support as

`F_ell(x0,x1,...,xN) = sum_S x0^(N-|S|) prod_(i in S) xi`

and try to prove that `F_ell` is Lorentzian. This route already fails at the
smallest honest quotient. For `diag(2,2)` and line `(0,1)`, the support consists
of the two masks `{0,1}` and `{2,3}`. Their exponent vectors are

`(2,1,1,0,0)` and `(2,0,0,1,1)`.

Move one unit from the first occupied site of the left vector toward either
site occupied only on the right. One exchange produces rank-zero diagonal
pairs; the other produces rank-one pairs in the orthogonal line. Neither
symmetric exchange stays in the fixed-line support. Hence the support is not
M-convex. Since M-convex support is necessary for a nonnegative homogeneous
Lorentzian polynomial, no Hessian calculation can rescue this `F_ell`.

The rank sequence remains much better behaved than its support. Nevertheless,
two natural real-rooted strengthenings are also false:

- At `diag(2,3)`, line `(1,0)`, the normalized coefficients on their support
  are `1/5, 3/5, 3/5`; the quadratic is proportional to `1+3z+3z^2`, whose
  discriminant is `-3`.
- At `diag(2,4)`, line `(1,0)`, the raw count polynomial after removing `z^2`
  is `4+24z+54z^2+48z^3+12z^4`. Its exact Sturm chain has only two real roots
  for squarefree degree four.

Even the normalized-matching property is too strong. The first full layer-flow
failure in increasing HNF order is `((11,3),(0,1))`, line `(1,-3)`, from layer
7 to layer 8. The maximum flow is `605/726`. The violating lower cut contains
11 of 66 states and has no neighbor among the 11 upper states. Mask `471` is
already a singleton witness: adding any of its four missing sites jumps the
homology rank directly from one to two.

This does **not** refute the ULC observation. It sharpens it. Through `N<=12`,
every fixed-line sequence remains strictly ULC on positive support, yet three
standard mechanisms that would imply it are exactly false. The surviving
conjecture is therefore rank-sequence-only:

`A_k^2 C(N,k-1) C(N,k+1) >= A_(k-1) A_(k+1) C(N,k)^2`.

A useful next proof must attack this two-step inequality directly—perhaps by a
weighted two-step injection or an orbitwise boundary double count—rather than
impose a global exchange or layer-expansion axiom that the geometry does not
have.
