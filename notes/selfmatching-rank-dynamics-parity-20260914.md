# Self-matching rank dynamics splits into exactly two parity relaxation functions

Date: 2026-09-14. Exact finite two-time algebra plus triangular L=3,4 Fourier controls. Addendum to draft #773 / #784 / #776.

## 1. Odd/even basis for the full rank sigma-algebra

At a self-matching symmetric point let

    X=r-1 in {-1,0,1},
    P(X=-1)=P(X=1)=a,
    P(X=0)=r0=1-2a.

Define the centered even observable

    Y=X^2-2a.

Then

    E X=E Y=0,
    E[XY]=0.

The three functions `{1,X,Y}` span every observable of the three-state rank variable. Hence any two-time rank law is completely determined by its stationary one-time probabilities plus the 2x2 matrix of correlations among X and Y.

## 2. Complement symmetry diagonalizes the two-time problem

For a self-matching model at the symmetric point, configuration complement acts by

    X -> -X,
    Y ->  Y.

Any stationary dynamics which commutes with complement therefore gives, at every time t,

    E[X_0 Y_t]=E[Y_0 X_t]=0.                              (1)

Thus the entire nontrivial two-time rank law is controlled by only two scalar functions

    C_o(t)=E[X_0 X_t],
    C_e(t)=E[Y_0 Y_t].                                    (2)

No assumption that the projected rank process is Markov is needed.

## 3. Explicit reconstruction of the 3x3 joint table

Reversibility plus complement symmetry imply the stationary two-time table has the form

             t: -1    0     +1
    0:-1       A     B      C
      0        B     D      B
     +1        C     B      A

with row sums `A+B+C=a` and `2B+D=r0`.

The two correlations are

    C_o=2(A-C),                                           (3)

and, since `Y=X^2-2a`,

    C_e=2(A+C)-4a^2.                                     (4)

Therefore

    A=[C_e+4a^2+C_o]/4,
    C=[C_e+4a^2-C_o]/4,
    B=a r0-C_e/2,
    D=r0^2+C_e.                                          (5)

Equations (5) reconstruct the full two-time rank kernel from the two parity relaxation functions.

At t=0, `C_o=2a`, `C_e=2a r0`, giving the diagonal stationary table. At infinite separation both correlations vanish and (5) becomes the independent product law.

## 4. Fourier parity gives the two relaxation functions directly

Under independent p=1/2 resampling noise with correlation rho, the Walsh expansion gives

    C_o(rho)/Var(X)
      = sum_(S odd) P_X(S) rho^|S|,                       (6)

because the self-matching theorem forces every even Fourier coefficient of X to vanish.

For Y, complement-even symmetry gives the dual statement:

    Yhat(S)=0 for every odd |S|,                          (7)

after removing its level-zero mean. Hence

    C_e(rho)/Var(Y)
      = sum_(S even, S!=empty) P_Y(S) rho^|S|.            (8)

So odd and even rank dynamics have disjoint Fourier supports at every finite self-matching size.

## 5. Exact triangular L=3,4 controls

A full integer Walsh-Hadamard transform of triangular-site rank at p=1/2 gives exact parity separation.

For X, all nonconstant even-level Walsh numerators vanish. The spectral mass fractions are:

L=3:

    level 1: 0.712444
    level 3: 0.238328
    level 5: 0.045745
    level 7: 0.003389
    level 9: 0.000094.

L=4:

    level 1: 0.615471
    level 3: 0.278853
    level 5: 0.082393
    level 7: 0.019589
    level 9: 0.003252
    level11: 0.000411
    level13: 0.000030
    level15: 0.000001.

For `Y=X^2-E X^2`, all odd-level Walsh numerators vanish exactly. Its nonzero spectral mass is:

L=3:

    level 2: 0.260241
    level 4: 0.419277
    level 6: 0.298795
    level 8: 0.021687.

L=4:

    level 2: 0.164173
    level 4: 0.317839
    level 6: 0.310667
    level 8: 0.155414
    level10: 0.041543
    level12: 0.008525
    level14: 0.001624
    level16: 0.000216.

These are exact finite Fourier controls, not scaling laws.

## 6. A continuum dynamical target

If the topological spectral samples have nondegenerate critical scaling on the `L^(3/4)` cardinality scale, then under time/resampling scale

    epsilon=t L^(-3/4)

the two normalized correlations should have universal fixed-modulus limits

    F_o(t;tau),
    F_e(t;tau).                                          (9)

Equations (5) would then give a universal complete two-time rank law without constructing a Markov rank process.

This target is stronger than a single autocorrelation but lower-dimensional than the full microscopic dynamical-percolation process.

## 7. Near-critical extension

A further conditional target is

    F_o(lambda,t;tau), F_e(lambda,t;tau),                (10)

where lambda is the near-critical thermal coordinate. Matching/complement symmetry should exchange `(lambda,X)` with `(-lambda,-X)` while leaving the even sector invariant. Establishing (10) rigorously requires an interface from the near-critical/dynamical scaling limit to torus rank events; it is not supplied by the finite algebra alone.

## 8. Interfaces

- #784: compute both X and Y spectral samples; X alone does not reconstruct the full two-time rank law.
- #776: triangular self-matching gives exact finite parity regressions for all Fourier levels.
- #778: paired birth-time data describe monotone-parameter persistence, while (5)--(10) describe dynamical resampling time; these are distinct process directions and should not be conflated.
- #585: the pair `F_o,F_e` is a source-explicit dynamical continuum fingerprint which does not require naming a CFT field first.

## 9. Boundaries

- Equations (1)--(8) are exact finite statements under self-matching/complement-symmetric stationary dynamics.
- The universal scaling functions (9)--(10) are conjectural until the topological-rank spectral sample is controlled.
- No Markov closure of the projected rank process is claimed.
