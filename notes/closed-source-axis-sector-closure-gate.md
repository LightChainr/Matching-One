# Axis topology and sector odds add no gate beyond `beta=L/m^2`

## Result

Let the equal-area pair be the axis quotient of side `L=5k` and its
Gaussian `(4k,3k)` companion.  Put

\[
 r=m^{-1},\qquad c={L\over m},\qquad \beta={L\over m^2},
 \qquad J(c)=I_0(2c)^2-I_1(2c)^2.                 \tag{1}
\]

Assume `m -> infinity` and `beta -> 0`.  Combining the proved complete
one-carrier estimate in `closed-source-axis-one-carrier-beta-gate.md` with
`closed-source-axis-beta-cloud-closure.md`, the remaining topology and
rank-sector terms are uniformly negligible under the same gate.  For the
existing positive angular denominator `Delta`,

\[
 {U\over A_N}=-{L^2\over\Delta}\,m^{-(2L+1)}J(c)\{1+o(1)\}<0. \tag{2}
\]

In particular there is no additional `L^3/m^5` gate.  Such a scale appears
if one moves the bare stripe polynomial to the corrected root but omits the
same common polymer pressure from the stripe transfer.  In the effective
two-phase thermal coordinate the two changes are one reparametrization and
cancel from slope-normalized original `U`.

No enumeration, fit or new coupling point is used.

## 1. The primary vertical direction is already included

For one essential occupied component the winding barrier is

\[
                         \gamma_{axis}=2L-1.       \tag{3}
\]

The primitive `(1,0)` and `(0,1)` directions are co-leading; the latter is
not an omitted vertical sector.  Their translations and both orientations
are already in the prefactor of (2).  Relative to the exact directed
arbitrary-run endpoint determinant, the complete physical one-carrier
transfer, including every horizontal reversal, self-avoidance deletion and
fixed-radius collar, obeys

\[
 {Z_{1c}\over J_{run}}
 =1+O(\beta+m^{-1}+L/m^3)=1+o(1)                 \tag{4}
\]

uniformly in `c`.  The directed factor itself contains the full capillary
entropy and is asymptotic to the Bessel endpoint `J(c)>0`, with
`J(c)=e^(4c)/(8 pi c^2){1+O(c^-1)}` for large `c`.

## 2. An extra essential component is automatically suppressed

If a rank-one configuration has `e` essential occupied components and
`c0` zero-image components, the component barrier gives

\[
 g\ge (2L-1)+(e-1)(2L-2)+2c_0.                  \tag{5}
\]

Thus `e>=2` pays at least `2L-2` beyond the one-carrier sector.  Relax all
new boundaries to arbitrary lattice walks and retain every translation and
reversal.  The contour count is at most exponential in its length and the
vertical-run pressure is `exp[O(c)]`.  Even after division by the endpoint
determinant, the relative contribution is bounded by

\[
 \operatorname {poly}(L,m,c)
 \exp\{-(2L-2)\log m+C(L+c)\}=o(1).              \tag{6}
\]

Indeed `c/L=1/m`, so `C(L+c)=O(L)` while `L log m` diverges.  This closes
the explicit extra-carrier boundary left open by (4).

## 3. Transverse winding and the companion do not compete

For a relaxed no-west boundary, let `D` be net transverse displacement.
Its normalized moment generating function is

\[
 E e^{\theta D}
 =\left({1-2r\over1-2r\cosh\theta}\right)^{L+1}. \tag{7}
\]

Taking `theta=log(m/4)` gives the uniform large-deviation estimate

\[
                 P(|D|\ge L)\le e^{-L\log m+O(L)}.            \tag{8}
\]

The all-reversal multiplier is `1+O(beta)` with an exponentially summable
displacement kernel, so it preserves (8).  Hence a primary carrier which
also winds transversely is negligible.

Every non-primary primitive slope on the axis torus has Manhattan period
at least `2L`, hence two boundaries cost an extra `2L`.  The `(4k,3k)`
companion has systole `7k=(7/5)L`, so its first barrier exceeds the axis
barrier by `4k=(4/5)L`.  Staircase entropy is only `exp[O(L)]`; therefore

\[
 {Z_{other\ slope}\over Z_{axis}}
 \le e^{-2L\log m+O(L+c)}=o(1),\qquad
 {Z_{1,tilted}\over Z_{1,axis}}
 \le e^{-(4/5)L\log m+O(L+c)}=o(1).              \tag{9}
\]

Thermal differentiation adds only polynomial factors and does not change
these conclusions.

## 4. The two-cloud root is a coordinate change, not a new gate

Factor the rank-zero and rank-two laws into common local pressures and
quotient-sensitive remainders.  Use the effective thermal coordinate

\[
 \xi(h,m)={1\over N}\log {Z_2^{loc}(h,m)\over Z_0^{loc}(h,m)}. \tag{10}
\]

Near the two-phase root, `xi_h>0`.  The independent bases are
`B=1+m^-2 h` and `C=h+m^-4`; connected contractible corrections merely
change the analytic map `h -> xi`.  Equal-area quotients have the same
embedding density for every sub-systolic animal.

For any such common analytic coordinate change, original `U` is exactly
invariant:

\[
                 {Y_h\over Q_h}={Y_\xi\xi_h\over Q_\xi\xi_h}
                                ={Y_\xi\over Q_\xi}.          \tag{11}
\]

This is why moving the numerical `h`-root while leaving the bare carrier
unchanged produces a false susceptibility enhancement.  The same pressure
dresses the carrier, and (11) removes it.  The complete cloud calculation
leaves only the geometry-dependent collar remainder

\[
 {\delta Y_{collar}\over Y_{carrier}}
 =O(\beta/m)+O(\beta/m^3)+O(e^{-\kappa L\log m})=o(1).       \tag{12}
\]

## 5. Rank-zero/rank-two odds cannot bypass rank one

Configurationwise,

\[
 q=rank-1\in\{-1,0,1\},\qquad E=q^2=1-1_{rank=1}.            \tag{13}
\]

Thus the angular `E` numerator is exactly minus the rank-one probability
projector.  Redistributing weight between ranks zero and two cannot enter it
directly; it acts only through normalization, the pooled root, or its slope.

After the common pressure in (10) is removed, a rank-zero/rank-two animal
which distinguishes the quotients must see a deck translate.  Its boundary
cost is at least `2 ell_1-2`, so the difference of restricted log odds obeys

\[
 |\epsilon_{axis}-\epsilon_{tilted}|
 \le \operatorname {poly}(L,m)
 e^{-(2L-2)\log m+C(L+c)}=o(1).                 \tag{14}
\]

This includes long almost-winding rank-zero paths.  Their raw source order
may precede the signed rank-one derivative, but (13) prevents direct entry
into `Y`; their root correction multiplies the already exponentially small
carrier.  The pooled root centers the common odds, rank-one mass is
exponentially small, and consequently

\[
                      Q_h={N\over2}\{1+o(1)\}>0               \tag{15}
\]

in the effective chart, or times the positive common Jacobian in `h`.

## 6. Closed error budget

The full relative remainder is

\[
 O(\beta+m^{-1}+\beta/m)
 +\operatorname {poly}(L,m,c)
  \left[e^{-(2L-2)\log m+C(L+c)}
       +e^{-(4/5)L\log m+C(L+c)}\right]=o(1).     \tag{16}
\]

Equations (4)--(16) close the extra-carrier, homology, cloud and sector-odds
boundaries and give (2).  This is a strong-source theorem for the specified
axis/tilted sequence and existing `S*`, `q`, `E`, pooled root and angular
normalization.  It does not assert a fixed-`m` phase or continuum exponent.
The minimal conclusion is sharper: topology supplies no hidden condition
beyond

\[
                              \boxed{L/m^2\to0}.
\]
