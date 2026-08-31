# A common-label tangent preserves both instantaneous laws and exposes cross-geometry birth response

The common-label experiment has a sharply separable new mechanism:
**cross-orientation first-birth susceptibility is supported only on R0/R0,
whereas R0/R1 permits an R0 contact mark to load the other orientation's
completion and lifetime.** Both statements are exact consequences of the
prefix information, not assumed exchange symmetry. The existing S/D fork
tensor can recover these channels without new tails.

## One actual joint next-label policy

Let Z include the common ordered prefix and both labeled NN torus geometries
f,s. The remaining label set has size d. Write r_i for the prefix rank,
r_i(u) for the rank after adding u, e_i(u) for its occupied NN contact degree,
and c_i(u) for the number of touched occupied components.

For a joint degree a=(e_f,e_s), define

\[
 A_a=\{u:r_f(u)=r_f,\ r_s(u)=r_s,\ (e_f(u),e_s(u))=a\},
 \qquad \pi_a=|A_a|/d.
\]

On these jointly rank-preserving labels set

\[
 g_f=1_{\{r_f=0\}}(e_f-c_f),\quad
 g_s=1_{\{r_s=0\}}(e_s-c_s),\quad
 g_+=(g_f+g_s)/2,\quad g_-=(g_f-g_s)/2.
\]

The R0 restriction is substantive: a safe R0 insertion's new graph cycles
are contractible; this interpretation is not assigned to e-c in R1/R2.
Use the two-parameter policy

\[
 q_{t_+,t_-}(u\mid Z)=\pi_a
 \frac{\exp\{\pi_a(t_+g_+(u)+t_-g_-(u))\}}
 {\sum_{v\in A_a}\exp\{\pi_a(t_+g_+(v)+t_-g_-(v))\}},\quad u\in A_a.
\]

All other labels keep probability 1/d; empty classes are omitted. After
the chosen label, draw the original uniform remaining permutation, shared
by the two orientations. This is one common label policy, not the average
of two orientation-specific interventions in9ce53a5a/7c60b8a7.

Every class retains mass pi_a at every finite parameter. On that class both
new ranks and both **NN graph Euler increments** `Delta chi_i=1-e_i` are
fixed. All other label probabilities are unchanged. Thus the *entire joint*
distribution of `(r_f(u),r_s(u),Delta chi_f,Delta chi_s)` is exactly invariant,
conditional on Z. No claim about an Euler convention that also attaches
filled plaquettes is needed or implied.

For any tail output Y let `m_Y(u)=E[Y|Z,u]`. Differentiation at zero gives

\[
 R_{Y,\alpha}(Z)=\partial_{t_\alpha}E_t[Y\mid Z]\big|_0
 =\sum_a\pi_a^2\operatorname{Cov}_{A_a}(g_\alpha,m_Y),\quad \alpha=+,-.
\]

Indeed, the conditional score is `pi_a(g_alpha-bar g_alpha,a)`; integrating
it against class mass pi_a produces the squared weight. Every function of
Z and the immediate joint tuple has zero response, at all finite t as well
as to first order. A nonzero future response demonstrates failure of
conditional-mean closure on that tuple within the same old prefix.

For iid uniform U,V, including U=V, the existing quartet estimator is

\[
 \widehat R_{Y,\alpha}=
 1_{\{U,V\text{ in the same }A_a\}}
 \frac{(g_\alpha(U)-g_\alpha(V))
 [Y_{U0}+Y_{U1}-Y_{V0}-Y_{V1}]}4.
\]

Its conditional expectation is the boxed-class covariance sum above.
There is no additional orientation-mixture half: g+/g- and S already
contain their declared halves. Keep the full original20k-prefix denominator
and original20 batches, including zero-contribution cells. Products of
class estimates or division by sampled small class frequencies are not
required.

## Exchange covariance, not automatic zeros at a fixed geometry pair

Let J exchange the two geometries and all their marks in Z. It maps joint
degree class `(e_f,e_s)` to `(e_s,e_f)`, preserves pi_a, and sends
`(g_+,g_-)` to `(g_+,-g_-)`. The exact finite-parameter identity is

\[
 q^{JG}_{t_+,t_-}(u\mid JZ)=q^G_{t_+,-t_-}(u\mid Z).
\]

For any same-semantic orientation output Y define
`S=(Y_f+Y_s)/2` and the **raw** half-difference `D0=(Y_f-Y_s)/2`.
They have exchange parities + and -. With rows `(S,D0)` and columns `(+,-)`,
the susceptibility matrix obeys

\[
 R(JG,JZ)=P\,R(G,Z)\,P,\qquad P=\operatorname{diag}(1,-1).
\]

Hence `R_S,+` and `R_D0,-` are exchange-even, while `R_S,-` and `R_D0,+`
are exchange-odd. The odd entries vanish only after a genuinely
exchange-invariant ensemble average (or an explicitly symmetrized one).
A fixed ordered pair of different Gaussian geometries is not automatically
such an ensemble just because its permutations are uniform and paired.
An actual label/geometry bijection preserving the source law would have to
establish that stronger symmetry. Relabeling and averaging a record with
its exchanged copy can manufacture the odd zeros; it supplies no extra
physical evidence or independent sample.

The production difference is instead
`D=(Y_f-Y_s)/delta`, `delta=cos(4theta_f)-cos(4theta_s) != 0`.
Thus `D0=(delta/2)D`. Under exchange **delta also changes sign**, so this
normalized D is exchange-even. With both output rows `(S,D)`, its matrix
transforms as `R_norm(JG,JZ)=R_norm(G,Z)P`: plus-column responses are even
and minus-column responses odd. Holding the numerical delta fixed while
swapping names is a different convention. Neither normalized plus response
may be set to zero by raw-D oddness.

## Recover the physically named cross channels from the saved tensor

For an orientation output Y define, with the same joint-class weights,

\[
 C_{ij}[Y]=\sum_a\pi_a^2\operatorname{Cov}_{A_a}(g_i,m_{Y_j}),
 \qquad i,j\in\{f,s\}.
\]

These are source-mark/response susceptibilities, not a symmetric covariance
matrix of interchangeable quantities. There is no general reciprocity
requiring `C_fs=C_sf`. In the production S/D normalization,

\[
 \boxed{C_{fs}=R_{S,+}+R_{S,-}
              -\tfrac\delta2(R_{D,+}+R_{D,-}),}
\]
\[
 \boxed{C_{sf}=R_{S,+}-R_{S,-}
              +\tfrac\delta2(R_{D,+}-R_{D,-}).}
\]

In particular

```
C_fs + C_sf = 2 R_S,+ - delta R_D,-,
C_fs - C_sf = 2 R_S,- - delta R_D,+.
```

The following source supports are strict at each prefix, and their nulls
hold in each saved quartet, not just asymptotically:

| Prefix ranks | Cross first-birth response | Cross completion/lifetime response |
|---|---|---|
| 00 | C_fs[F1] and C_sf[F1] both allowed | Both directions allowed |
| 01 | Both zero | f-source to s-completion allowed; reverse zero |
| 10 | Both zero | s-source to f-completion allowed; reverse zero |
| 02 or20 | Both zero | Both zero: receiving R2 has both births in the prefix |
| Both ranks >=1 | Entire policy tangent zero | Entire policy tangent zero |

The reason is structural: g_i vanishes outside R0, while K1 of an R1/R2
orientation is already known; both K1,K2 are known in R2. Thus the two
displayed cross combinations for F1 are **pure00 channels**. They cannot be
generated by mixing the one-R0 prevalence cells. In mixed cells, g- equals
g+ for0r and -g+ forr0, so the input columns are collinear. For F1 the raw
`(S,D0)` matrix there is proportional to `[[1,1],[1,1]]` or
`[[1,-1],[-1,1]]`, respectively. The vanishing cross combinations therefore
need no orientation-exchange assumption.

For the canonical/integrated responses this same reconstruction commutes
with all linear birth identities:

```
C_ij[K_j] = -(N+1) C_ij[F_j,integral],
C_ij[Wclock] = -(N+1)(C_ij[F2,integral]-C_ij[F1,integral]),
R_A - R_E = 2 R_F1,   R_A + R_E = 2 R_F2,
E = 1-F1+F2,         Wclock=K2-K1.
```

Most directly, **cell01 gives `C_fs[Wclock]=C_fs[K2]`**, since K1s is fixed;
cell10 gives `C_sf[Wclock]=C_sf[K2]`. A nonzero value means a common safe
label weighted by the R0 side's contact structure changes the other
side's remaining completion time, although neither immediate rank/Euler
outcome law changed. It is cross-geometry response under a real shared-label
selection policy, not evidence of temporal path memory or physical contact
between the two separately embedded systems. Cell00 separately tests
whether first-birth information transmits through that same source mark.

## A local source check with genuine sign content

Taking the *chosen-label mark itself* as output gives

\[
 T_{\beta\alpha}=\partial_{t_\alpha}E_t[g_\beta]\big|_0
  =\sum_a\pi_a^2\operatorname{Cov}_{A_a}(g_\alpha,g_\beta).
\]

This two-coordinate contact-response matrix is symmetric PSD. A zero
quadratic form means that source combination is constant in every
positive-mass joint class and its policy is exactly uniform for every
finite tilt along that combination: **all** future responses in that
direction must vanish. In00 the classes already fix both degrees, so plus
is a tilt in minus total touched-component count and minus is a tilt in
their imbalance. Distinguishing these directions identifies a finite
geometric source; it does not identify a continuum field. Nonzero source
rank can coexist with zero response of a particular future observer, and
positivity of T imposes no sign on a birth susceptibility.

## Scientific scope and handoff

This strengthens the earlier orientation-mixture result by conserving the
entire paired instantaneous law under one actual next-label rule. A nonzero
normalized-D response is directly a paired-loading sensitivity; a nonzero
S response alone need not imply one. Neither can be inferred by averaging
the old single-orientation tangent curves: the common safe/joint-degree
classes and cross terms are different.

The proposed decisive readout is the pair of pure00 cross-first-birth
channels alongside the directional01/10 completion/lifetime channels,
canonical and integrated, using exactly the same source/batch covariance.
This selects finite mechanisms without assuming zeros from geometry names,
claiming a unique loop-edit cause, or identifying a CFT spin/thermal/Jordan
operator from the label-exchange Z2 alone. It does not change the
unperturbed global mean or supply independent evidence.

Sources: the contact interpretation at e67d9b90; nested shared-label fork
semantics at84018e19 and source e32a8593/contact959a7fa2; the
orientation-specific policy at9ce53a5a and its full shape at7c60b8a7.
The root's common-policy collection/readout definition is frozen at
ffb70969. This note adds exact identities only: no sampling, DP, replay,
empirical re-scoring, PR or issue comment. Production and the unique shared
batch covariance remain with the root and its coordinator.
