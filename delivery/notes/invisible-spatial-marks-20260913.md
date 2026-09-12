# Entire uniform thermal/rank-source data do not determine spatial-source response

Date: 2026-09-13. Exact finite countermodels with strictly positive weights.
These are marked, dependent site measures, NOT competing assertions about the known
unmodified Bernoulli percolation law. Their purpose is to delimit what a restricted
observer dictionary can identify.

## 1. Two invisible marks on the actual 4x4 torus

Let O_square be the translation/D4 orbit of the occupied 2x2 square, row masks
[3,3,0,0]. It has 16 elements. Consider two other contractible four-site shapes:

    O_L: [3,1,1,0],  orbit size 128;
    O_T: [7,2,0,0],  orbit size 64.

Every support configuration has r=0 and total occupancy K=4. The three orbits are
pairwise disjoint. Define

    f_L = 8*1_{O_square} - 1_{O_L},
    f_T = 4*1_{O_square} - 1_{O_T}.

Each is translation/D4 invariant. Both sums over all configurations vanish within
EACH (r,K) cell, since 8*16-128=0 and 4*16-64=0.

For either f choose W_+/-=1+/-delta*f, with delta_L=1/16 and delta_T=1/8.
Every weight lies between 1/2 and 3/2. They are strictly positive, p-independent
marks; they do not introduce negative probabilities.

At uniform p, every four-site configuration has the same Bernoulli weight
p^4(1-p)^12. Hence exactly, for every p in [0,1],

    E_p W_+=E_p W_-=1,
    P_{W+}(r=j,K=k)=P_{W-}(r=j,K=k)=P_Bernoulli(r=j,K=k).  (1)

Therefore all of the following also agree, at EVERY p:

- M, E_top, every rank probability, and every homogeneous thermal derivative;
- the entire joint generating function E[exp(s*(r-1)+u*K)];
- every mixed homogeneous thermal/unmarked-rank-source jet;
- the finite unperturbed matching root and its homogeneous thermal slope.

Knowing the WHOLE rank/K curve, rather than a few moments, still cannot distinguish
these marked finite families from one another using that restricted dictionary.

## 2. Two spatial sources distinguish two independent hidden directions

Apply an additive site-probability field p_i=p+epsilon*h_i, where h has zero mean
and +/-1 values. The original rank observable is retained. Put H=sum_i h_i n_i.
On the support K=4, so separate differentiation of the product measure gives

    d^2/d epsilon^2 E_{p+epsilon*h} f |0
       = p^2(1-p)^10 sum_omega f(omega) H(omega)^2.          (2)

The diagonal score term is constant on this K=4 support and cancels by sum f=0.
The first derivative cancels by translation symmetry. Let q_delta(epsilon h) be
the weighted matching root. At the unperturbed root p_*, the normalizer is 1 and
M=0, so it drops from the root equation only. Since X=-1 on the support,

    q_delta''(h)-q_0''(h)
      = delta * p_*^2(1-p_*)^10
        * sum f H^2 / M'(p_*).                             (3)

For the stripe h=(-1)^x and checkerboard h=(-1)^(x+y), the integer contrast matrix
(rows=sources, columns=L/T marks) is

              f_L     f_T
    stripe   -512    -128
    checker     0    -256

with determinant 131072 != 0. Thus the two marks are independently distinguishable
by these spatial second-order responses, despite complete equality in (1).

At p_*=0.5906721123310283..., the plus-L mark changes the stripe root curvature by
-0.00029607493663376017..., and the checkerboard curvature by zero. The plus-T mark
changes the stripe curvature by -0.00014803746831688008... and checkerboard curvature
by -0.00029607493663376017.... Opposite marks reverse these changes.

These small numbers are not a practical sample-size recommendation. The result is
exact distinguishability, not a claim that a noisy experiment separates the marks
cheaply. Finite-amplitude roots at 1/64 and 1/128 validate (3) at high precision.
For the L mark, the checkerboard response is unchanged for every amplitude: all
support configurations have two occupied sites of each checkerboard colour, so the
weighted sum f vanishes even before differentiation. This is a genuine null control.

## 3. Boundary of the inference

The Bernoulli ensemble is mathematically specified, so there is no ambiguity about
what its response actually is. The countermodels instead disprove the inference:

    identical complete uniform rank/thermal data
        => identical response to a newly addressed spatial source.

A known microscopic source can be evaluated directly, as in the full Hessian work.
An unknown candidate mechanism must supply that source map; additional derivatives
of the same spatially unmarked generating function do not supply it.

This does not supply the two named #275 continuum candidates or replace original-U.
No source definition is changed in a historical freeze. The code uses PR708's
verified rank table, enumerates only these finite symmetry orbits, and stores every
orbit member and integer contrast in invisible-spatial-marks.json.
