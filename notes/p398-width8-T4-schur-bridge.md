# P398: T4 is a fourth-order, directionally asymmetric geometric bridge

**The first source-correlation change is fourth order, not third.** The
existing hierarchy keeps both G psi and G^2 psi inside the old 7-space.
The T4 residual is first reached on the third forward step but has a
nonzero one-step left return to the source. That asymmetric four-step path
has opposite signs on the two protected rays, matching the observed
directions of their tail repairs.

This uses only the same 7->8 named extension as `30eef34a`, with no new
source, span, lag, or fit. Full left/right couplings, their declared-column
coordinates, and moments k=0..6 are in
`results/p398-width8-T4-schur-bridge/latest.json`.

## One pole, with separate left and right couplings

Let e=psi/||psi|| and let w be the physical T4 residual orthogonal to the
old 7-space, divided by its positive real norm, without rephasing. In an
orthonormal old-space basis whose first column is e, the extended backward
generator is

\[
H_8=\begin{pmatrix}A&b\\c&d\end{pmatrix}.
\]

Eliminating only w gives the matrix-valued memory bridge

\[
K_w(t)=b e^{dt}c,\qquad \Sigma_w(z)=b(z-d)^{-1}c.
\]

With R_7(z)=(z-A)^(-1), the exact scalar source-resolvent correction is

\[
\widehat u_8(z)-\widehat u_7(z)
=\frac{(e^*R_7 b)(cR_7e)}{z-d-cR_7b}.
\]

This is a fixed Schur complement, not a fitted decay. Its bare poles are
d_minus=-9.2923363936 and d_plus=-9.2991171155. These fast bridge poles
are not the full process's slow physical eigenmodes near 2; their coupling
can nevertheless shift the slow source-visible poles.

## Why the cubic term must vanish

The already established exact identities are

\[
G\psi_s=-3\psi_s+\sqrt2\,s\zeta\,P_sT_2,
\qquad GT_2=2T_3-2T_2+S_{11}-B_2.
\]

Every term on the second right side, and its Kreweras projection, is in
the old span. Consequently

\[
ce=\langle w,Ge\rangle=0,\qquad
cAe=\langle w,G^2e\rangle=0.
\]

For source moments m_k=e*H^k e, the quadratic increment is
(e*b)(ce)=0. The cubic increment is

\[
\Delta m_3=(e^*Ab)(ce)+(e^*b)(cAe)+(e^*b)d(ce)=0.
\]

All quartic paths also contain one of those zero factors **except**

\[
\boxed{\Delta m_4=(e^*b)(cA^2e)
=\langle e,Gw\rangle\langle w,G^3e\rangle.}
\]

Equivalently, cR_7e first begins at z^(-3), while e*R_7b begins at z^(-1);
the resolvent correction first begins at z^(-5), hence the fourth moment.
The bridge dwell rate d does not enter this leading coefficient.

## The actual named path

The phase below is fixed by T4, not selected from the result.

| Ray | One-step left return e*b | Three-step right entry cA^2e | Delta m4 |
|---|---|---|---:|
| Minus | -.0313928810(1+i) | -.8949377808(1-i) | **+.05618935045** |
| Plus | -.0347034680(1+i) | +3.9037079787(1-i) | **-.27094440991** |

Thus the two rays have similar direct left-return magnitudes, but opposite
signs and very different strengths in the three-step entry. The new plus
response initially shifts downward, whereas minus shifts upward:

\[
u_8-u_7=\Delta m_4\,t^4/24+O(t^5).
\]

These signs agree with the earlier finite-distance repair directions:
plus's overestimated tail is reduced, while minus's underestimated tail
is increased. The leading coefficient alone is not used to predict a
finite-distance error magnitude.

This is a concrete non-selfadjoint chain: the forward source cannot reach
w in one or two steps, but the adjoint side has an immediate return.
It is not an assertion that individual configurations follow a four-jump
microscopic path, nor a claim about morphism memory.

## Deleting current exposes a different, tiny second-order bridge

Under the same projection, S replaces H_8 by (H_8+H_8*)/2. Since ce=0 for
G, the reversible source couplings become e*b_S=(e*b)/2 and
c_S e=conj(e*b)/2. Therefore

\[
\Delta m_{2,S}=|e^*b|^2/4
=|\langle w,Je\rangle|^2.
\]

These are exactly the small amounts of current-source squared norm added
by T4 in `33c6028f`: .000492756488 on minus and .000602165346 on plus.
Thus deleting current does not simply suppress the same bridge: it removes
the forward cancellation and creates a tiny reciprocal source gate.

The complete moment increments are:

| Dynamics / ray | Delta m0..1 | Delta m2 | Delta m3 | Delta m4 | Delta m5 | Delta m6 |
|---|---:|---:|---:|---:|---:|---:|
| G minus | 0 | 0 | 0 | +.0561893505 | -1.63976114 | +30.46109673 |
| G plus | 0 | 0 | 0 | -.2709444099 | +3.94124952 | -18.48371409 |
| S minus | 0 | +.0004927565 | -.0076196669 | +.1147566401 | -1.70203523 | +24.02958170 |
| S plus | 0 | +.0006021653 | +.0089067507 | -.2078557318 | +1.45745534 | +11.43345509 |

The exact old-space inclusions, not a floating-point zero threshold,
establish the G zeros. Direct k<=6 matrix moments and the separately saved
seven quartic path terms exhibit the nonzero coefficients. No repeated
validation campaign was run.

This sharpens the previous conclusion: T4's principal plus repair is an
indirect geometric-transport bridge with a delayed forward entrance and
immediate adjoint return. It is neither a missing instantaneous current
observer nor a new slow state added by fitting. Everything belongs to the
same deterministic width-eight finite-process dependency block.
