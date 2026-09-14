# The generic-Q graph-polynomial root is an exact topological-source balance

Date: 2026-09-14. Algebraic bridge between the rank-source program, Jacobsen's torus graph polynomial and the graded Q-tangent program. Addendum to draft #773 / #636 / #746.

## 1. Three topological sectors and a closure source

Let `Z_0,Z_1,Z_2` denote the finite torus partition sums in the rank/homology sectors 0,1,2 for the declared Potts/FK or percolation model. Introduce a bounded topological source coupled to `X=r-1`:

    Z_top(s)=e^(-s) Z_0 + Z_1 + e^s Z_2.                 (1)

The source does not modify the local transfer weights. It only deforms the closure functional among the three topological sectors.

The two extreme sourced contributions are balanced when

    e^s Z_2 = e^(-s) Z_0,

or

    Z_2/Z_0 = e^(-2s).                                   (2)

## 2. Jacobsen's graph-polynomial condition is exactly (2)

The toroidal graph-polynomial criterion uses

    P_B(Q,v)=Z_2-Q Z_0,

so its root satisfies

    Z_2/Z_0=Q.                                           (3)

Comparing (2) and (3) gives the exact chemical potential

    boxed: s_Q = -1/2 log Q.                             (4)

Therefore the finite graph-polynomial estimator is precisely a **sourced topological balance root**.

For Q=1, `s_Q=0` and one recovers the ordinary rank-0/rank-2 balance used by square-site Matching One. For generic Q the source is nonzero but known exactly.

## 3. Canonical rank coordinate

Use

    b=1/2 log(Z_0/Z_2).

Equation (3) is simply

    boxed: b = -1/2 log Q = s_Q.                         (5)

Thus the natural generic-Q odd coordinate is the shifted variable

    b_tilde = b + 1/2 log Q,                             (6)

and the graph-polynomial root is the section `b_tilde=0`.

The even rank-shape coordinate

    d=log[Z_1/sqrt(Z_0 Z_2)]                             (7)

is invariant under the source translation and therefore carries topological-shape information not fixed by the graph-polynomial balancing convention.

## 4. Q tangent decomposes into known source motion plus genuine shape response

At fixed microscopic couplings, differentiating (5) with respect to `log Q` shows that the explicit topological-source part of the Q tangent is

    d s_Q / d log Q = -1/2.                              (8)

Hence a generic-Q tangent analysis should separate:

1. the exact, kinematic odd translation `b -> b + (1/2) log Q`;
2. the response of the source-invariant even coordinate d;
3. any additional response from Q-dependence of the microscopic/local weights and representation sector.

Failing to subtract item 1 mixes a known closure chemical potential with genuine continuum/source information.

This gives a concrete continuum benchmark for the graded Q->1 machinery of #746.

## 5. Long-cylinder sector crossing

At fixed width and large longitudinal size m, suppose

    Z_0 ~ c_0 lambda_0^m,
    Z_2 ~ c_2 lambda_2^m.

The sourced balance condition gives

    m log(lambda_2/lambda_0)
      + log(c_2/c_0)
      + 2s = 0.                                          (9)

For the graph-polynomial source `2s=-log Q`,

    m log(lambda_2/lambda_0)
      + log(c_2/c_0)
      - log Q =0.                                        (10)

Thus the familiar leading-eigenvalue crossing is the extensive limit of the same sourced balance equation; coefficient/prefactor mismatches produce the known `1/m` root displacement unless cancelled.

The exact width-two bridge is one explicit realization of (10).

## 6. Finite-size correction interpretation

At a continuum critical point, define the shifted odd mismatch

    B_L = b_L(p_c,Q) + 1/2 log Q.                         (11)

The finite graph-polynomial root displacement is driven by B_L divided by the thermal slope `partial_p b`. Therefore a proposed `L^-4` correction should be interpreted as the scaling of a **residual matching/topological odd coordinate after the exact Q chemical potential has been removed**, not as the scaling of the raw sector odds.

This is the generic-Q version of the Q=1 8-arm/matching-odd question.

## 7. Interfaces

- #746: report the explicit `-1/2` topological-source piece of the Q tangent separately from the even shape d response.
- #636: a topological source should be implemented as a closure-functional deformation of the same local transfer operator, not as a new microscopic local field.
- #681/#705: width-two open/closed spectral crossing and finite-length coefficient corrections are sourced-balance effects of the same equation.
- #768: any continuum selection rule for the odd correction should apply to the residual coordinate (11).
- #585: map-resolved torus solution spaces can be tested on both the odd sourced-balance coordinate and the independent even shape coordinate.

## 8. Boundaries

- The algebra above assumes the declared `Z_0,Z_1,Z_2` are exactly the sectors entering the graph-polynomial definition; no unrelated rank convention may be substituted.
- Equation (4) identifies the closure-source component of Q. It does not say all Q-dependence of a Potts/FK model is a topological source; local Boltzmann weights and representation data also vary with Q.
- No LCFT/Jordan identification follows from the sourced balance itself.
