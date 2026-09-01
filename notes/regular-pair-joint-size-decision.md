# Canonical joint U: one conditional macroscopic-size decision

**Next output: a fixed-window, fixed-source original-U dilation
comparison.** The completed [N25 joint readout](../results/regular-pair-joint-u/REPORT.md)
gives

```
J2_total = -0.0055194314248394015,
J2_NN    = -0.001751074454402799,
J2_nonNN = -0.0037683569704366022.
```

These reject the stated exact global additive and NN-contact-only
closures for canonical `Kreg=K2bar+K0bar`. Non-NN on N25 still includes
microscopic separations. It supplies no macroscopic scaling window.
The useful next question keeps the canonical tensor and asks which
specified field-loading model predicts its **macroscopically separated
joint transmission to original U**. No K3 source, adjustable singlet
coefficient, or fitted counterterm is introduced here.

Under the explicit single-field assumptions below, the extensive
window response obeys

\[
 T_N:=N^2J_{2,\mathrm{macro}}(N)
       \sim C_{\mathcal W}N^{29/8-x}.
\]

Thus an area dilation by four predicts `2^(-5/4)` for x=17/4 and
`2^(-13/4)` for x=21/4. These ratios differ by a factor of four. They
are conditional model predictions, not measured exponents or proved
field identities. This note runs no computation or tests and freezes
no execution contract.

## 1. The observer and window must be fixed together

Use the same physical quotient families and source as the old U:

```
axis: (5k,0), tilted: (4k,3k), N=25 k^2, L=sqrt(N),
Delta4=1152/625, Kreg=K2bar+K0bar,
each vacant vertex: T0 -> T0 + (epsilon/N) Kreg.
```

The continuum tori have the same square shape. The physical local
lattice spacing, common N/E/S/W port orientation, and microscopic
normalization of Kreg stay fixed as k changes. No Q1 rank, q/E mark,
pooled-root definition, or slope normalization is changed.

Choose one symmetric displacement window before reading new results,
for example

\[
 \mathcal W_N(x,y)
  =\mathbf1\{c\le d_{\mathbb T_N}(x,y)/L\le d\},\qquad
                  0<c<d<\tfrac12.                         \tag{1}
\]

Here `d_T` is Euclidean shortest displacement on the physical torus,
not graph distance measured in an arbitrarily chosen quotient labeling.
The constants, boundary convention, and any angular weight are identical
at every dilation and in the two geometries; a radial window needs no
angular weight. The bound d<1/2 is a simple example avoiding the torus
cut locus, not a newly selected numerical window. A nonempty window
with positive continuum area has order N allowed displacements per
site, hence order N² ordered site pairs.

Extend the delivered signed g_xy by zero if either endpoint is occupied.
Define the window source and its original-U response by

\[
 S_{2,\mathcal W}(A)=\frac1{N^2}
       \sum_{x\ne y}\mathcal W_N(x,y)g_{xy}(A),\qquad
 J_{2,\mathrm{macro}}=\mathcal L_{U,N}[S_{2,\mathcal W}].      \tag{2}
\]

This is a fixed projection of the second derivative's **pair kernel**.
For a proper window it is not the unfiltered second derivative of the
homogeneous one-parameter vertex family. Equivalently, first allow
separate per-vertex couplings g_x, take their mixed functional
derivatives, and then sum them with the prescribed window. No full
finite-strength windowed family is needed to define this coefficient.

The complete linear functional L_U,N includes the separately centered
q/E moments, original pooled-root motion, and both thermal-slope terms,
as derived in [the joint-U interface](regular-pair-joint-u-functional.md).
Translation reduction gives the exact source-moment scale
`sum_y W_N(0,y) g16_0y/(16*N)`, under the full original probability law.
Neither a pair-count average nor a vacancy-conditioned normalization
may be substituted silently.

The corresponding per-vertex coupling is g=epsilon/N. Every mixed
second derivative therefore brings N² when converted from epsilon
units to g units, even after the linear window projection:

\[
              T_N=N^2J_{2,\mathrm{macro}}.                 \tag{3}
\]

This factor is exact at each finite N; scaling assumptions enter only
below. The macroscopic window removes collisions of the two *marked*
vertices because their minimum separation grows as cL. Deleting only
the four NN terms does not do this: fixed separations 2, 3, and so on
remain possible contact/OPE contributions at all larger N.

## 2. Derive the exponent through the intrinsic U interface

Write, for each source family,

\[
 M=\tfrac12(m_a+m_b),\quad
 Y=(e_a-e_b)/\Delta_4,\quad z=M(p),\quad
 X_N(z)=Y(p(z)),\quad U_N=A_N\partial_zX_N(0).              \tag{4}
\]

The repository convention is `A_N=N^(13/8)/2`; omitting its constant
1/2 changes amplitudes only. A simple pooled root with M_p nonzero is
required. The fixed-z coordinate removes a common thermal
reparameterization. Differentiating with respect to z contributes no
additional size power when the intrinsic scaling function is regular:
the p derivative and division by M_p carry the same thermal factor.

Let F_N(Q,z) denote the window-projected second response of this
intrinsic X curve to **per-vertex** couplings, before the Q derivative.
The following assumptions define the mechanism being compared:

1. The same specified continuum field of dimension x(Q) loads each
   insertion with an N-independent microscopic coefficient. Its scaling
   dimension and coefficients admit the stated Q continuation near 1.
2. With the fixed torus shape, angles, and window, that channel dominates
   the intrinsic response in a neighborhood of z=0:

   \[
   F_N(Q,z)=L^{4-2x(Q)}B(Q,z;\mathcal W)+R_N(Q,z).           \tag{5}
   \]

   The remainder is smaller at the required Q and z derivative order,
   uniformly enough to differentiate this expansion at Q1.
3. The activated intrinsic slope
   `partial_z partial_logQ B(1,0;W)` is nonzero. Mere existence of a
   field or an unmarked two-point function does not establish this.
4. No other same-order field, logarithmic partner, contact term outside
   the window definition, or size-dependent source normalization is
   supplying or canceling the claimed leading response.

The exponent in (5) follows from two integrations over a fixed fraction
of the torus, giving L^4, and two local fields of dimension x, giving
L^(-2x). At r/L fixed the coefficient is a **torus scaling function**;
one need not replace it by the plane formula r^(-2x). The q/E observer
and the intrinsic derivative are part of B, not optional factors added
after measuring an ordinary correlator.

Section 3 shows that a single regular channel's exact Q1 zero removes
the exponent-derivative logarithm. Applying (4) to the activated (5)
therefore gives

\[
 \begin{aligned}
 T_N&\sim\tfrac12\partial_z\partial_{\log Q}B(1,0;\mathcal W)
                    N^{13/8}L^{4-2x(1)}\\
    &=C_{\mathcal W}N^{13/8+2-x(1)}
      =C_{\mathcal W}N^{29/8-x(1)},\\
 J_{2,\mathrm{macro}}&\sim C_{\mathcal W}N^{13/8-x(1)}.
 \end{aligned}                                             \tag{6}
\]

For N to 4N, equivalently L to 2L:

| Specified common insertion dimension | T_N | T_(4N)/T_N | J2_macro(4N)/J2_macro(N) |
|---|---|---|---|
| x=17/4 | N^(-5/8) | 2^(-5/4) | 2^(-21/4) |
| x=21/4 | N^(-13/8) | 2^(-13/4) | 2^(-29/4) |

These are not the older one-insertion ratios for `N V_av`: two integrated
insertions change the power, and the fixed canonical regular source is
distinct from the old bounded-occupation tangent. The signed amplitude
C_W is not predicted here. A stable nonzero amplitude has the same sign
at successive sufficiently dominant scales; its ratio is positive even
when both responses are negative.

## 3. A Q derivative does not automatically force log L here

Every nonempty finite-network canonical Kreg insertion vanishes at Q1.
This includes the windowed two-insertion response, its normalizations,
and its full intrinsic q/E response. Under the single-leading-channel
and differentiated-remainder assumptions, (5) consequently requires
`B(1,z;W)=0` throughout the relevant z neighborhood.

For that one amplitude the exact derivative of the leading expression
is

\[
 \left.\partial_{\log Q}
  [L^{4-2x(Q)}B(Q,z)]\right|_1
 =L^{4-2x(1)}\left[B_{\log Q}(1,z)
       -2x_{\log Q}(1)B(1,z)\log L\right].                \tag{7}
\]

The second term is zero. In the factorized notation
`B(Q,z)=A(Q)Phi(Q,z)`, a regular single amplitude with A(1)=0 likewise
has only `A'(1)Phi(1,z)` at this order. A nonzero x'(1) alone does not
generate a surviving logarithm.

This reasoning does not remove every possible logarithm. For example,
two channels can satisfy `B_1(1,z)=-B_2(1,z) != 0` and
`x_1(1)=x_2(1)`, while their dimension derivatives differ. Their finite
Q1 sum vanishes, but its Q derivative contains
`-2 B_1(1,z)(x'_1(1)-x'_2(1)) L^(4-2x) log L`.
Jordan mixing or a Q-to-1 limit nonuniform in L also lies outside (5).
Such mechanisms are alternatives to the single regular channel, not
automatic repair terms to fit after a failed fixed-power comparison.

## 4. What the comparison can reject

The next bounded output is a comparison of the **same**
`T_N` and `T_(4N)` under the two declared models, not a free-exponent fit.
When the reference response is nonzero one may use their ratio. A
denominator-free alternative reports both fixed-model contrasts

\[
 D_{17}(N)=T_{4N}-2^{-5/4}T_N,\qquad
 D_{21}(N)=T_{4N}-2^{-13/4}T_N.                             \tag{8}
\]

The finite-window dominance assumption and admissible correction/error
criterion must be stated before readout; a bare asymptotic claim gives
no finite-N rejection threshold. A further predetermined dilation can
be held out when assessing approach to the ratios. Window, counterterm,
angle, correction exponent, and microscopic normalization should not
be tuned to restore a failed ratio.

- A resolved response with incompatible scaling rejects the stated
  single-field **loading model on that window**, not all occurrences of
  that dimension in the theory.
- A resolved sign reversal similarly rejects a same-leading-amplitude
  description over that declared scale range.
- A vanishing or unresolved macro response does not support either
  nonzero-loading power model. The finite N25 non-NN result still holds,
  but it cannot be promoted to macroscopic original-U transmission.
- Passing one ratio supports the specified power/loading combination;
  another field with the same dimension or a different mechanism with
  the same power remains possible. No universal normalization or unique
  local operator follows from one dilation.

The completed [spatial observer](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-observer.md)
is `C(x,y)=<g_xy>` at selected separations, without q/E source
crossmoments, the directional pair, pooled-root motion, or U-slope
response. Its positive noncontact signal and its L64/L32 ratio cannot
be inserted into (6) or (8). Even a correlator at fixed r/L would require
its own loading analysis; unmarked C can be governed by a sector that
the directional intrinsic U response removes.

C4 averaging allows spins 0, +/-4, +/-8, and higher multiples of four;
it is not a pure-spin or RG-eigenfield projection. Products of spin-four
components can themselves carry spin zero or eight, and the two-angle
Delta4 quotient can alias higher harmonics. Thus a nonzero U loading
at the claimed two-field order is a substantive assumption. If symmetry
annihilates its intrinsic slope, or if directional transmission first
requires an additional irrelevant anisotropy insertion, (6) does not
describe the leading power. Removing macroscopic pair collisions also
does not remove each local tensor's pre-existing microscopic field
mixing. These are reasons to compare specified mechanisms, not to call
an observed ratio proof of a pure four-leg field.

## Provenance and scope

The exact definitions and factor N² come from
[`regular-pair-joint-u-functional.md`](regular-pair-joint-u-functional.md)
(`7557da5271f85a69ea5426b61ce7e67b94ee8ff2`). The N25 producer and fixed
NN split were implemented in
`99b58fc18666cfa6d35b96b52bb84c78dec43a55`. The completed numerical
anchor is linked at the beginning; it was not recalculated for this
note. The old intrinsic-U normalization and dimension conventions are
documented in
[`local-pair-size-response-predictions.md`](local-pair-size-response-predictions.md),
whose first-insertion powers are kept separate here.

Equations (2)--(4) specify exact finite response identities. Equations
(5)--(8) are hypothesis-level, conditional finite-window predictions.
No source parameter was fitted, no kernel or occupation population was
regenerated, and no new numerical score, test, or server job was run.
