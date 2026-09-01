# #275 restricted-trace transport law

Status: theory-derived forward contract; zero new samples.  This note fixes the
same-source, same-observer and same-normalizer object required by the current
P0.  It does not identify Q4, Jordan, H8 or a continuum primary.

## Restricted source coordinates

For geometry `g` and homology rank `r=0,1,2`, write

\[
Z_{r,g}(p,\varepsilon)=Z_{r,g}(p)+\varepsilon B_{r,g}(p)+O(\varepsilon^2),
\qquad s_{r,g}=B_{r,g}/Z_{r,g}.
\]

Use the two independent normalized sector coordinates

\[
\eta={1\over2}\log{Z_2\over Z_0},\qquad
\xi=\log{Z_1\over2\sqrt{Z_0Z_2}},
\]

so that one microscopic source has tangents

\[
\eta_B={s_2-s_0\over2},\qquad
\xi_B=s_1-{s_0+s_2\over2}.
\]

`eta_B` is the rank-0/rank-2 bias.  `xi_B` is rank-1 loading relative to the
geometric mean of the two even sectors.  They are typed source coordinates;
neither may be replaced by a static C3 character or by a different observer's
amplitude.

With `q=r-1`, `E=q^2`, and normalized sector probabilities `P_r`, direct
differentiation gives

\[
q_B=(E-q^2)\eta_B-P_1q\,\xi_B,
\qquad
E_B=P_1q\,\eta_B-P_1E\,\xi_B.
\]

These identities retain the rank-1 denominator response that is absent from a
numerator-only selection rule.

## Moving-root original-U transport

Let the pooled matching condition be `M(p_0)=0`, let `D=M_p(p_0)`, and let the
same source move the root by

\[
w=h_M/D.
\]

Define the comoving source tangents

\[
a_g=\eta_{B,g}-w\eta_{p,g},\qquad
b_g=\xi_{B,g}-w\xi_{p,g}.
\]

For the repository's original `E` readout and H4 geometry projector
`P_4`, the complete response is

\[
\boxed{
{\partial_B U_g\over A_N}
={1\over D}\,\partial_p\,
\mathcal P_4\!\left[P_1q\,a_g-P_1E\,b_g\right]_{p=p_0}.}
\]

This is the compact form of the existing four-term root/slope formula.  It
shows exactly why a critical identity for `Z2-Z0`, or a zero of one
unnormalized numerator, does not by itself determine original `U`: `xi_B`, its
thermal derivative, and root transport remain.

## What the Arguin relation removes

At `Q=1`, the continuum homology identity gives `Z2=Z0`.  If the *same* source
insertion obeys the compatible relation `B2=B0`, then

\[
\eta_B=0.
\]

This removes rank-0/rank-2 bias, but it does not remove

\[
\xi_B=s_1-(s_0+s_2)/2.
\]

Only `s0=s1=s2` is a pure common normalizer and forces the complete normalized
response to zero.  Moreover, a relation proved only on the critical surface
does not determine its transverse thermal derivative.  If `B2=B0` is known in
a neighbourhood of the root, then `eta_B=eta_B,p=0`; otherwise `eta_B,p` must
remain in the contract.

## Modulus and phase law for an ordinary weight-4 column

For a candidate whose bottom-projected response is an ordinary weight-4 Q4
column, write

\[
\widehat E_4(\tau)=(\operatorname{Im}\tau)^2E_4(\tau).
\]

Under `gamma*tau=(a*tau+b)/(c*tau+d)`,

\[
\widehat E_4(\gamma\tau)=
\left({c\tau+d\over c\bar\tau+d}\right)^2\widehat E_4(\tau).
\]

Thus the modulus is invariant and the phase changes by
`4 arg(c*tau+d)`.  The same law may be used for a Jordan bottom-projected
`log N` slope; it must not be applied to an unprojected fixed-`N` Jordan-top
intercept, whose additive part is connection/gauge dependent.

## Minimal scoreable data object

The complete missing object is the aligned three-sector source--thermal jet

\[
\mathcal J_B=\{s_r,\partial_p s_r\}_{r=0,1,2}.
\]

For each batch, geometry, orientation and `p` point, retain at least

- `count_r`, `sum_B_r`, `sum_K_r`, `sum_KB_r`;
- `sum_dB_dp_r` when the source explicitly depends on `p`;
- the common batch/dependency identifier, orientation and H4 normalization;
- one shared partition denominator and one aligned delete-one unit for all
  three sectors.

Then

\[
\partial_p s_r=\operatorname{Cov}(K,B\mid r)
+E[\partial_pB\mid r].
\]

Existing K1/K2 full curves determine the baseline probabilities, `p` jets,
root and geometry transport.  The fixed-`p` rho-child archive also provides a
six-coordinate rank-0/rank-2 covariance.  Neither archive contains an
independent same-`B` sector insertion, so neither alone closes this contract.

## Decision use

The next candidate comparison should supply two actual, constrained columns in
this data object: for example a vacuum/Ward column and a thermal Q4/Jordan
column.  After the full normalizer and allowed nuisance amplitudes are retained,
score their images in the existing covariance.

- Different images permit one existing-data model-elimination score.
- Identical images certify non-identifiability and identify the one missing
  physical relation.
- Another untyped character, source, angle or radial fit does not substitute
  for this column.

Priority is attention allocation, not a task lock.  Parallel theory remains
welcome when it supplies a typed source column or a relation that changes this
rank calculation.
