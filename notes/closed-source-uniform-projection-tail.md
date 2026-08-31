# The two fixed laws have opposite original-U signs throughout m >= 64

**Uniform finite-coupling result.** For the locked N25 axis/tilted pair,
`Sstar=C+F4+Bvac` and `Sdrop=Sstar+r`, put `lambda=1/m`,
`A=25^(13/8)/2` and `Delta=1152/625`. For **every real m>=64**,

\[
-0.618102\,m^{-11}
 < U_*/A < -0.454124\,m^{-11}<0,
\]
\[
0<1.376734\,m^{-42/5}
 < U_{\rm drop}/A <1.844309\,m^{-42/5}.
\tag{1}
\]

These decimal constants are outward-rounded rationals, not estimates.
The exact-rational calculation also supplies classical uniform remainders:

\[
\left|U_*/A+\frac{625}{1152}\lambda^{11}\right|
 \le378\lambda^{13},
\qquad
\left|U_{\rm drop}/A-\frac{625}{384}\lambda^{42/5}\right|
 \le1.322966\lambda^{44/5}.
\tag{2}
\]

This strengthens the existing eventual opposite-tail statement to a
single predeclared half-line. In particular neither law can cross zero
again within this interval. It does not locate either crossover or the
smallest valid endpoint. No coupling points, configurations, random
samples, fits or new source parameters were introduced.

## Inputs and preserved observer

The complete `(K,g,q)` histograms are the already locked
`results/p337-closed-source-finite-coupling/{axis,tilted}.csv`, with hashes
in [score.json](../results/p337-uniform-projection-tail/score.json).
The projected tail and moving-root cancellations are those of
`09042093`/`9b88a49b`; the fixed deletion law and its stripe tail are
`fbbaa2aa:notes/topological-projection-reverses-global-u-tail.md`.
This is another deterministic consequence of those same exact populations.

For star use `h=y/m` and weights `h^K lambda^g`. For drop use
`d=y m^(-23/25)` and weights `d^K x^e`, where

\[
x=\lambda^{1/25},\qquad e=25(g-r)+2K,
\qquad r=q+1.
\]

Each geometry is normalized separately before pooling. If its partition,
q numerator and rank-one numerator are Z, Qn and R, define

\[
F=Qn_a Z_b+Qn_b Z_a,\quad
M_i=Qn_{i,h}Z_i-Qn_i Z_{i,h},\quad
P_i=R_{i,h}Z_i-R_i Z_{i,h},
\]
\[
\mathcal N=P_bZ_a^2-P_aZ_b^2,\qquad
G=M_a Z_b^2+M_b Z_a^2.
\]

Here h denotes the relevant common thermal coordinate h or d, and every
derivative is taken **before** root substitution. At the root `F=0`,

\[
\boxed{U/A=2\mathcal N/(\Delta G).}\tag{3}
\]

Indeed `E=1-R/Z`, and the pooled q slope is `G/(2 Z_a² Z_b²)`.
The coordinate-to-p Jacobian cancels between these slopes. No fixed-root
or unnormalized pooled numerator is substituted for the original U.

## One uniform implicit-root tube, not a mesh of root fits

All coefficients are integers from the histograms. For a polynomial
`P(x)=sum_(j>=s) c_j x^j`, the calculation repeatedly uses the exact bound

\[
|P(x)|\le x^s\sum_{j\ge s}|c_j|X^{j-s},\quad 0<x\le X.\tag{4}
\]

For star set `X=1/64`; for drop set `X=17/20` in x. The rational
inequality `(17/20)^25 >= 1/64` encloses the required whole domain.
It is a majorant radius, not an added physical coupling point or a search
for a different cutoff.

The ground part of F is `2(h^50-1)`, with derivative `100 h^49`.
Bounding all excited terms gives the following uniform certificates:

| Quantity | Star | Drop |
|---|---:|---:|
| Root tube | `1<=h<=1+2 lambda²` | `1<=d<=1+2 x^52` |
| Lower bound on `F_h` in tube | `99.3550300391` | `99.4362179252` |
| Leading term of `F(1)` | `-100 lambda²` | `-100 x^52` |
| Remainder divided by that power, absolute bound | `0.576067754` | `0.582044098` |

Thus F is negative at the left endpoint and positive at the right, and
has exactly one root in each tube. The existing unique-root result for
these two fixed laws identifies this as the original matching root.
Mean-value control then proves

\[
|h_0-(1+\lambda^2)|\le1.018855\lambda^6,
\qquad 0<d_0-1\le1.011524x^{52}.
\tag{5}
\]

The first approximation uses the previously established root series;
its residual is bounded as a complete polynomial, not as an unspecified
Taylor remainder. This is why the large fixed-h/root-motion cancellation
does not invalidate the sign certification.

## Angular cancellation, root motion and denominator bounds

Substituting `hbar=1+lambda²` for star and `dbar=1` for drop into the
complete angular numerator gives leading coefficients `-200 lambda^11`
and `+600 x^210`, respectively. The exact polynomial has no lower term.
Equation (4) bounds all remaining terms; (5) times a uniform bound on
`|partial_h N|` controls the additional root displacement.

| Bound after division by the leading monomial | Star | Drop |
|---|---:|---:|
| Polynomial remainder | `26.364845198` | `58.043590979` |
| Root-displacement remainder | `0.022357794` | `18.009012703` |
| Enclosure of the complete numerator | `[-226.387203,-173.612797]` | `[523.947396,676.052604]` |
| Enclosure of G | `[397.420120,414.824433]` | `[397.744871,412.947556]` |

Every denominator lower bound is positive. Inserting these intervals into
(3), before rounding outward, yields (1). For reproducibility the script
stores rational enclosing endpoints at resolution `10^-12`; all deciding
operations are exact `Fraction` arithmetic, not floating-point intervals.

For (2), the angular remainders start two lambda powers later for star,
and ten x powers later for drop. The root-motion errors start still later.
The positive denominator majorants minus their ground value400 have
nonnegative coefficients beginning at `lambda²` or `x^52`. Factoring
these powers in (4) before division gives the exact upper constants
`377.987005108782` and `1.322965501785`, rounded outward to (2).

## What the finite-coupling comparison now establishes

Dropping the rank projection does not merely alter an inaccessible
formal leading coefficient: it reverses the original root-normalized
angular response throughout the explicit interval `t>=log64`. The
pressure-density equivalence at fixed t and increasing N cannot erase
this finite-volume observer distinction. The proof keeps N25 fixed;
it is not a uniform-in-size estimate or a continuum-field identification.

The contract was committed before the bound calculation; the selected
half-line never changed. Code:
[p337_uniform_projection_tail.py](../scripts/p337_uniform_projection_tail.py).
Final exact-bound calculation took about0.45 seconds in the local research
Python environment. No enumeration, MC, cloud operation or test suite was
run. The only intermediate execution issue was serialization of very large
rational integers; output was changed to outward rational enclosures.

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p337_uniform_projection_tail.py
```
