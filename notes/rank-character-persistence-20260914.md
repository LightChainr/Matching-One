# Rank character and two-birth persistence: an exact two-dimensional topology coordinate

Date: 2026-09-14. Additive continuation of draft PR #773. All algebra in §§1–5 is finite-system exact once the digital-Alexander rank variable `r in {0,1,2}` and the standard monotone site coupling are fixed. Scaling interpretations are separated in §6.

## 1. One complex character contains the entire rank law

Let

    X = r-1 in {-1,0,1},
    omega = exp(2 pi i/3),
    chi = E[omega^X].

Writing `P_j=P(r=j)`, direct discrete Fourier inversion gives

    chi = P0 omega^(-1) + P1 + P2 omega,
    Re chi = (3P1-1)/2,
    Im chi = (sqrt(3)/2)(P2-P0).                         (1)

Conversely, for `chi=x+iy`,

    P0 = (1-x)/3 - y/sqrt(3),
    P1 = (1+2x)/3,
    P2 = (1-x)/3 + y/sqrt(3).                           (2)

Hence the admissible `chi` are exactly the closed equilateral triangle with vertices `1,omega,omega-bar`. This is not a Potts/parafermion field identification; it is simply the nontrivial Z3 Fourier character of the three rank values.

The matching/complement involution sends `X -> -X`, therefore

    chi -> conjugate(chi).                               (3)

The balance condition `P2=P0` is exactly `Im chi=0`: the physical rank-law path crosses the real axis in the character triangle.

## 2. Canonical rank-simplex coordinates and source orbits

Define

    b = (1/2) log(P0/P2),
    d = log[P1/sqrt(P0 P2)].                             (4)

Equivalently, after normalization,

    (P0,P1,P2) proportional to (exp(b), exp(d), exp(-b)). (5)

Thus `b` is matching-odd and `d` is matching-even. The previous `c=P1/[2 sqrt(P0P2)]` is `c=exp(d)/2`.

A topological source `sX` changes the weights by `(e^{-s},1,e^s)`. Therefore exactly

    b -> b-s,
    d -> d.                                              (6)

So `d` is a source invariant and the topological source is a pure translation in the odd natural coordinate. In the complex character plane, constant-d source orbits are conics. Eliminating `P_j` from (2),(4),

    (1+2x)^2 = exp(2d)[(1-x)^2-3y^2].                   (7)

This supplies an exact geometric picture of the rank-source deformation: it moves along one fixed conic until the real axis is reached.

At balance `b=0`, write `P0=P2=a`, `P1=1-2a`. The natural-parameter Fisher information for `(b,d)` is the covariance of the sufficient statistics `(1,0,-1)` and `(0,1,0)`, hence

    I_bd = [[2a, 0],
            [0, 2a(1-2a)]].                             (8)

The matching-odd and matching-even rank directions are therefore exactly Fisher-orthogonal at balance. This is an exact finite statement, not an asymptotic CFT orthogonality claim.

## 3. Source zeros depend only on the even coordinate

At balance,

    Z_X(s) = a e^{-s} + (1-2a) + a e^s
           = 2a[cosh s + exp(d)/2].                     (9)

Thus the nearest complex source zero is determined entirely by `d` (or `c`). When `c<1`,

    s = +/- i theta,    theta=arccos(-c).                (10)

The critical square-torus Pinson/Newman-Ziff control used elsewhere in #773 has

    a_* = 0.30952627542983132276969363335...,
    P1_* = 0.38094744914033735446061273330...,
    c_* = 0.615371746084117...,
    d_* = 0.207618451946263...,
    chi_* = 1-3a_* = 0.071421173710506... .              (11)

These are topology-law constants at the square-torus critical continuum point. They do not identify a microscopic field.

## 4. The rank process is exactly a two-bar persistence object

Use the standard monotone coupling: give each of the N sites an iid uniform label and occupy all sites whose label is <=p. Since adding occupied sites cannot reduce the ambient homology-image rank,

    r(p): 0 -> 1 -> 2

with at most two birth times

    T1 = inf{p:r(p)>=1},
    T2 = inf{p:r(p)=2},    T1<=T2.                       (12)

Thus the whole ambient-rank filtration is a two-bar persistence barcode: two classes are born at `T1,T2` and never die in the ambient `H1(T^2)` image.

The one-time rank probabilities are exactly

    P0(p)=P(T1>p),
    P1(p)=P(T1<=p<T2),
    P2(p)=P(T2<=p).                                     (13)

The static rank law therefore gives both birth-time marginals, but not their copula.

Because `b(p)` is strictly decreasing, define canonical birth coordinates

    B1=b(T1),  B2=b(T2),   B1>=B2.                      (14)

Then

    P1(b) = P(B2 < b <= B1).                            (15)

Consequently the intrinsic curve `d=D_L(b)` is simultaneously a static rank equation of state and a one-point coverage function of the random canonical persistence interval `[B2,B1]`.

The canonical gap has the exact first-moment identity

    E[B1-B2]
      = integral_R P1(b) db
      = integral_0^1 [-b'(p)] P1(p) dp.                 (16)

Unlike the raw `p` gap, this quantity is invariant under every monotone thermal reparameterization used to define the same rank-law curve.

Higher gap moments require multi-time rank-one persistence. For `k>=2`,

    E[(B1-B2)^k]
      = k(k-1) integral_{u<v} (v-u)^(k-2)
          P(B2<u<v<=B1) du dv.                          (17)

So the next process-level object after the one-time curve is the two-time rank-one persistence kernel, not another one-point exponent.

## 5. Discrete permutation births and an exact static formula for the mean gap

Let a uniformly random site permutation be revealed one site at a time, and let

    J1=min{j:r_j>=1},
    J2=min{j:r_j=2},
    D=J2-J1 >=0.                                        (18)

Conditional on `J1,J2`, the continuous birth times are order statistics. In particular, conditional only on `D=d`,

    T2-T1 ~ Beta(d, N+1-d),                             (19)

with a point mass at zero for `d=0`. Therefore

    E[T2-T1 | D=d] = d/(N+1).                           (20)

If `C_{1,k}` is the exact number of k-site configurations with rank one, then during a random permutation the number of discrete occupation levels spent at rank one is exactly D. Hence

    E D = sum_{k=0}^N C_{1,k}/binom(N,k),                (21)
    E(T2-T1)=E D/(N+1).                                 (22)

This is useful because `E D` requires only the static rank-sector occupation polynomials requested in #775; no paired permutation archive is needed.

For the exact controls committed with this note:

    square L=3:      E D = 3/2,
    square L=4:      E D = 14122/6435 = 2.19456099456...,
    triangular L=3:  E D = 3/2,
    triangular L=4:  E D = 991/455 = 2.17802197802....   (23)

The scaled values `E D/L^(5/4)` are respectively about `0.37992, 0.38795, 0.37992, 0.38502`. This close small-size agreement is a positive control only; it is not an exponent or amplitude fit.

A simultaneous two-rank birth has `D=0`. If `Delta_v X=2` denotes a rank jump from 0 directly to 2 when site v is inserted, exchangeability and beta integration give the exact identity

    P(D=0) = N integral_0^1 P_p(Delta_0 X=2) dp.         (24)

This separates two questions that had previously been conflated:

- the absolute geometry of simultaneous two-rank births (a 6-arm/8-arm fusion diagnostic), and
- the matching-odd one-point correction controlling the finite balance-root shift.

They need not have the same exponent if the lower-arm geometry cancels in the matching-odd projection.

## 6. Scaling directions opened by this representation

### 6.1 Self-normalized continuum equation of state

If a near-critical rank law exists, eliminate the nonuniversal thermal coordinate and ask directly for

    d -> D_*(b;tau).                                    (25)

For a self-matching model `D_*` is even. Square/matching finite-size odd parts measure matching-asymmetric irrelevant corrections without paying an extra thermal exponent.

### 6.2 Birth-gap scaling

If `T2-T1=O(L^{-3/4})`, then `D=O(L^(5/4))` because N=L^2. Thus `D/L^(5/4)` is a natural process-level observable for #775 or existing paired threshold-rank archives. Its full law contains information not present in `P_j(p)` alone.

### 6.3 Character trajectory as a topology-only validation object

`chi(p)` is a two-dimensional path inside a fixed equilateral triangle. Complement matching is conjugation, topological source trajectories have the exact conics (7), and balance is real-axis crossing. A continuum candidate for the rank law should reproduce this entire trajectory after thermal reparameterization, not merely one root or one source-zero angle.

## Reproduction

`scripts/rank_character_persistence_controls.py` postprocesses exact L=3,4 rank-sector coefficient tables for square and triangular controls. It checks the character inversion and conic identity and records the exact mean permutation gap (21). No Monte Carlo is used.
