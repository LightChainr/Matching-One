# Signed two-cloud interface transfer closes the mesoscopic axis window

## Result

Put `tau=m^-1`, `N=L^2`, and `gamma=2L-1`.  The absolute marked
contour estimate in
[`closed-source-axis-capillary-uniformity.md`](closed-source-axis-capillary-uniformity.md)
gave the sufficient condition `L^2/m ->0`.  That condition is not the
true scale: it marked each record by `K<=N` before complement pairing.

The signed transfer proves the fixed-volume law uniformly under

\[
                         \boxed{L/m\longrightarrow0}.              \tag{1}
\]

Namely,

\[
 E_{1,h}(h_0)=
 -{L^2(L^2-6L+6)\over2}m^{-(2L+1)}\{1+o(1)\},     \tag{2}
\]

and, because the companion quotient has a longer winding barrier,

\[
 {U\over A_N}=-{L^2-6L+6\over\Delta}
                  m^{-(2L+1)}\{1+o(1)\}.           \tag{3}
\]

Thus the whole formerly unresolved window `L << m <= L^2` retains the
old coefficient.  There is no intermediate `L^2/m` mechanism.

The finite capillary scale also closes.  Set

\[
 c={L\over m},\qquad \mathcal I_d=I_d(2c),\qquad
 J_1(c)=\mathcal I_0^2-\mathcal I_1^2.             \tag{4}
\]

For bounded `c`, the complete two-cloud directed-interface limit is

\[
 \boxed{{U\over A_N}=-{L^2\over\Delta}
        m^{-(2L+1)}\{J_1(c)+o(1)\}.}               \tag{5}
\]

Since `I_0(x)>I_1(x)>0` for finite positive `x`, `J_1(c)>0`.
Capillary fluctuations renormalize the magnitude but cannot reverse
the strong-source sign.  A one-cloud calculation produces a spurious
width/area scaling function, even a plausible finite zero; the
rank-one configuration necessarily supports both black singletons in
the exterior and white singleton holes in the interior.  At the exact
root their two pressures cancel every area tilt.  The only surviving
finite-`c` factor is the positive nonintersection determinant (4).

No enumeration, fitting or sampling is used.  The fixed-volume input
is the completed first-two-shell classification at `762dbaf4`.

## 1. Exact complement-paired mark

Let `A` be a connected rank-one occupied annulus whose complement is
also connected rank one.  Write

* `K=|A|`;
* `a_i(A)` for the number of occupied vertices having an NN neighbour
  outside `A`;
* `a_o(A)` for the number of vacant vertices having an NN neighbour
  in `A`.

At the next source order an isolated occupied site can be inserted in
`M(A)=N-K-a_o(A)` positions; for the complement the count is
`M(A^c)=K-a_i(A)`.  Insert `h^K+h^(N-K)` and the two decorations into
the exact root/normalization functional of `762dbaf4`,

\[
 \mathcal C[S,T]=
 -\left({S\over1+h^N}\right)''_{h=1}
 -\left({T\over1+h^N}\right)'_{h=1}
 +\left({NhS\over(1+h^N)^2}\right)'_{h=1}.          \tag{6}
\]

Direct differentiation gives

\[
 \boxed{\mathcal C(A,A^c)=
 { (N+2-2K)a_i+(2K-N+2)a_o\over4}.}                \tag{7}
\]

The bulk mark `K=O(N)` survives only multiplied by the boundary
imbalance `a_i-a_o`.  Summing (7) over straight stripe/complement pairs
recovers exactly

\[
 \sum\mathcal C(A,A^c)
       =-{L^2(L^2-6L+6)\over2}.                    \tag{8}
\]

This identity already shows why the earlier absolute `N` mark was too
large.

## 2. The true expansion parameter is `L/m`

Suppose the two essential boundaries have total excess `2j` over
`2L`.  Resolve alternating dual vertices by turning around occupied
corners.  If the boundary layers are disjoint, only turns and one-cell
necks change `a_i-a_o`, and

\[
 |a_i-a_o|\le2j,\qquad a_i+a_o\le4L+4j.            \tag{9}
\]

If the layers overlap the first bound can be `O(L)`, but overlap is
possible for only `O(j)` widths.  After transverse translation and the
width sum, (7) is `O_j(L^4)` before selecting turn columns.  An
excess-`2j` directed record has at most `2j` turn columns, contributing
at most `C^jL^(2j)` choices.  Hence the full complement-paired shell is
bounded by

\[
                         C^j(1+j)L^{2j+4}.          \tag{10}
\]

Relative to (8), its weighted parameter is `(CL/m)^2`.  Zero-image
nonsingleton polymers have aggregate activity `O(N/m^4)` by
`closed-source-mesoscopic-black-gas.md`; under the same signed boundary
mark they are `O((L/m)^4)+o(1)`.  The growing long-contour cutoff absorbs
the remainder.  Summing (10) proves (1)--(3).

## 3. Both reference clouds are mandatory

Let a flat stripe have width `w`.  A black singleton can be placed in
the vacant bulk at NN distance at least one from the occupied annulus;
a white singleton hole can be made by deleting an occupied site at NN
distance at least one from either boundary.  Their counts are

\[
 M_w=L\max(L-w-2,0),\qquad
 H_w=L\max(w-2,0).                                  \tag{11}
\]

The first operation raises the source exponent by two and has activity
`ah`; the second raises it by four and has activity `a^2/h`, where
`a=m^-2`.  The latter is legal precisely because a site counted by
`H_w` lies in the interior of the occupied annulus: deleting it leaves
the occupied component connected with the same ambient rank and creates
one contractible white singleton.  Adjacent interior holes are excluded
from this reference cloud and have total activity `O(Na^4)=o(1)`.
Deleting a boundary site is not a white-hole gas event: it changes the
cut contour and is counted by the interface transfer below.  Its
non-directed local remainder is at most `O(La)=O(c^2/L)` for bounded
`c`.

The flat rank-one weight is therefore

\[
 W_w(h)=h^{Lw}(1+a^2/h)^{H_w}(1+ah)^{M_w}.          \tag{12}
\]

The two rank-zero reference gases are

\[
 Z_0=(1+ah)^N,qquad Z_2=(h+a^2)^N.                 \tag{13}
\]

At their exact common root

\[
 h_*=1+a,qquad B:=1+ah_*=h_*+a^2,                 \tag{14}
\]

the bulk widths `2<=w<=L-2` obey

\[
 W_w(h_*)=h_*^{2L}B^{N-4L},                        \tag{15}
\]

independent of `w`.  Equation (15) is the decisive cancellation.  A
black-only cloud leaves an artificial `q=h/(1+ah)` width score; the
white-hole factor replaces the interior `h` powers by the same base
`h+a^2=B` and removes it.

Differentiate the finite width sum, retaining the two narrow endpoints.
The bulk scores form an arithmetic progression whose centered sum is
zero, while the endpoint terms give

\[
 {\sum_w W_w'(h_*)\over\sum_wW_w(h_*)}
 -{N\over2B}(1+a)
 ={c^2\over2}+O_c(L^{-1}).                         \tag{16}
\]

No `-c^4/12` width term survives.  The same equality holds for a rough
annulus: changing the signed bridge area transfers sites between
`H(A)` and `M(A)`, but both carry the common base `B` at (14).  Thus the
area mark cancels before the interface sum.  This also explains why a
one-boundary `cosh` factor cannot be thermally marked in isolation.

## 4. Strict two-boundary transfer

In the resolved digital graph the two cut boundaries are edge-disjoint
simple loops and cannot cross or touch.  At `tau=c/L`, a horizontal
reversal uses two excess edges but loses a freely placeable column, so
it is `O(1/L)` relative to a directed record.  A surviving boundary is
a directed bridge.  For transverse displacement `d`, its kernel is

\[
 \sum_{j\ge0}{c^{2j+d}\over j!(j+d)!}=I_d(2c)=\mathcal I_d. \tag{17}
\]

The argument is `2c`, not `c`: every vertical cut edge carries `1/m`,
and choosing its column supplies `L`.  For two bridges initially
separated by `d`, strict nonintersection gives the two-path
Lindstrom--Gessel--Viennot determinant

\[
 J_d(c)=\det\begin{pmatrix}\mathcal I_0&\mathcal I_d\\
                            \mathcal I_d&\mathcal I_0
          \end{pmatrix}
       =\mathcal I_0^2-\mathcal I_d^2.              \tag{18}
\]

The determinant is valid only after occupied-corner resolution; using
the unresolved degree-four dual vertex would permit false contacts.
The `O(L)` bulk widths have transfer `I_0(2c)^2+o(1)`.  The thermal
response is an endpoint effect: width one and full-minus-one carry
occupation leverage `O(L^2)`.  Their strict survival fraction relative
to two independent bridges is

\[
 \rho(c)={J_1\over\mathcal I_0^2}
         =1-\left({\mathcal I_1\over\mathcal I_0}\right)^2.         \tag{19}
\]

The rank-one amplitude itself carries `I_0(2c)^2`; multiplying it by
the endpoint score `(c^2/2)rho(c)` from (16),(19) leaves exactly
`(c^2/2)J_1(c)`.  Restoring the fixed barrier and pooled slope proves
(5).

Finally, positivity is elementary from the integral representation

\[
 I_0(x)-I_1(x)={1\over\pi}\int_0^\pi
 e^{x\cos\theta}(1-\cos\theta)d\theta>0,           \tag{20}
\]

and `I_0+I_1>0`.  Hence `J_1=(I_0-I_1)(I_0+I_1)>0` for every finite
`c`; the putative capillary sign reversal is excluded.

## Scientific boundary

The proved finite-lattice consequences are:

* `L/m->0` is sufficient for the old derivative and original-U
  coefficient; the former `L^2/m->0` gate was nonsharp.
* bounded `L/m` is governed by the positive determinant (4), so the
  strong-source sign cannot change through capillary roughening.

This does not interchange a fixed-`m` thermodynamic limit, identify a
continuum field, or extend the directed-interface reduction to
unbounded `L/m`.
