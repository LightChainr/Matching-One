# Mixed curvature in the fixed Euler-preserving source coordinates

The four saved finite policies give an exact interaction contrast in the
named physical source coordinates. It equals a **box average of the mixed
Hessian**, not the Hessian at zero. Both are supported only by original00
prefixes. A nonzero value excludes additive response in these coordinates;
it does not demonstrate noncommuting operations, time ordering, a second
independent mode, or a continuum field identity.

The allocation6bace935 and frozen
`notes/p334-common-label-euler-tangent-definition.md` supply class
conservation and orientation relabeling. They supply **no independent
source-reflection involution**. None is assumed below. In particular,
negating a parameter is an alternative probability law, not a proven
symmetry of the labeled geometry or its response.

## The exact two-source family and conserved instantaneous law

Fix the complete original prefix Z, its d vacant labels and joint-safe
degree classes A_a. Write n_a=|A_a| and pi_a=n_a/d. The physical marks are
the R0-only loop coordinates
`L_i=1[old_rank_i=0](e_i-c_i)`, i=f,s. The natural parameters are fixed by

\[
 q_{\boldsymbol t}(u\mid Z)=\pi_a
 {\exp\{\pi_a(t_f L_f(u)+t_s L_s(u))\}
  \over\sum_{v\in A_a}\exp\{\pi_a(t_f L_f(v)+t_s L_s(v))\}},
 \quad u\in A_a.
\]

Every other label retains probability1/d. Empty classes are omitted.
The suffix after the selected label stays uniform and shared by the two
geometries. Each A_a retains mass pi_a for all finite parameters. Its two
immediate ranks and graph Euler increments `(1-e_f,1-e_s)` are fixed.
Therefore the entire joint instantaneous rank/Euler law is unchanged,
conditional on every original prefix, throughout the parameter plane.

Let m_F(Z,u) be the unchanged conditional suffix mean of a fixed observer.
All derivatives below act on q alone, not on m_F or on the population of Z.
The finite family is analytic for real parameters and its partial
derivatives commute: `H_fs F=H_sf F`.

## Density Hessian, class centering and its integer numerator

At zero, define

\[
 s_i(u)=\pi_a(L_i(u)-\mu_{ai}),\qquad
 T_{ij}(u)=s_i(u)s_j(u)-\pi_a^2
                    \operatorname{Cov}_a(L_i,L_j).
\]

Then `partial_i q|0=q0 s_i` and
`partial_i partial_j q|0=q0 T_ij`, with q0=1/d. T is the second
**density** score: the Hessian of log q alone is only the negative
covariance term and is not a response estimator.

For exact class sums
`S_i=sum_a L_i`, `Q_ij=sum_a L_i L_j`, n=n_a,

\[
 s_i(u)={nL_i(u)-S_i\over d},\qquad
 \boxed{T_{ij}(u)=
 {(nL_i(u)-S_i)(nL_j(u)-S_j)
                  -(nQ_{ij}-S_iS_j)\over d^2}.}
\]

All numerator entries are integers. Their centering is exact class by
class: the sum of the product term is `n(nQ_ij-S_iS_j)`, exactly canceled
by summing the subtracted constant over n labels. Thus both sum_a s_i and
sum_a T_ij vanish. Outside active classes both scores are zero. Every
function of the conserved immediate tuple consequently has zero first
and second response, in addition to the finite-law invariance above.

For a future observer,

\[
 H_{ij}F=E_Z E_{U\,{\rm uniform}}[T_{ij}(U)m_F(Z,U)]
        =E_Z\sum_a\pi_a^3\,
             \kappa_a(L_i,L_j,m_F),
\]

where the last quantity is the ordinary third joint cumulant within the
uniform class. The same cumulant formula holds at a finite parameter with
the tilted within-class law. Mixed curvature measures a mark-mark-response
cumulant, not merely the mark Gram entry. In particular,
`Cov_a(L_f,L_s)=0` alone does not force H_fs F=0.

The exact-centered two-label estimator used in this round is

```
(T_ij(U)-T_ij(V))*(Fbar_U-Fbar_V)/2,
```

with iid uniform U,V and two suffixes averaged per label. No same-class
pair mask is applied. Its expectation is H_ij F because E_U T_ij=0.
The response Hessian need not be PSD: its ff/ss entries can have either
sign, and `H_fs^2<=H_ff H_ss` is not a valid generic restriction. The PSD
Hessian of a log partition function is a different object.

## Why the mixed term has exact00 support

If either original rank is nonzero, its L_i is identically zero on every
label. The family is then independent of that parameter for every finite
t. Therefore all derivatives involving that source vanish, and in
particular T_fs and H_fs vanish outside original00. This is a structural
zero, not an ensemble cancellation or an orientation symmetry.

The whole-source mixed response is thus
`P(G=00) E[H_fs F|G=00]`. The new64-quartet00 extension must retain the
original full-prefix/batch denominator, not replace that contribution by
an unweighted conditional00 mean. Old8 and new64 are separately retained
estimators with shared original prefixes, not independent populations.

## The saved plus/minus points form the physical rectangle exactly

Because `g_plus=(L_f+L_s)/2` and `g_minus=(L_f-L_s)/2`, the saved scalar
policies at parameter lambda have physical coordinates

```
plus(lambda):  (t_f,t_s)=(lambda/2, lambda/2),
minus(lambda): (t_f,t_s)=(lambda/2,-lambda/2).
```

Their already saved lambda=+/-1 points are the four corners of
`[-1/2,1/2]^2`. For Phi(t_f,t_s)=E_t F, their contrast is

\[
 \mathcal R=\Phi(1/2,1/2)+\Phi(-1/2,-1/2)
             -\Phi(1/2,-1/2)-\Phi(-1/2,1/2)
 =\int_{-1/2}^{1/2}\!\int_{-1/2}^{1/2}
          \partial_f\partial_s\Phi(x,y)\,dy\,dx.
\]

The rectangle has area1, so the saved contrast is numerically its **average
mixed curvature**. For a general symmetric half-width h, the contrast
divided by4h^2 is the average. Repeated one-dimensional fundamental
theorems of calculus prove the identity without a small-parameter
approximation.

Locally, with both plus/minus parameters available as coordinates,

```
H_++ = (H_ff+2H_fs+H_ss)/4,
H_-- = (H_ff-2H_fs+H_ss)/4,
H_+- = (H_ff-H_ss)/4,
H_fs = H_++-H_--.
```

These are basis identities, not independent evidence. In general R is not
H_fs(0): its small-box expansion is
`R(h)/(4h^2)=H_fs(0)+h^2*(H_fffs+H_fsss)/6+O(h^4)`.
For h=1/2 the leading correction has coefficient1/24. A zero local
Hessian can coexist with a nonzero finite contrast, and a nonzero local
Hessian can cancel on box averaging.

## Additivity, source-coordinate changes and a one-mode counterexample

On a connected parameter rectangle the additive form
`Phi(t_f,t_s)=c+A(t_f)+B(t_s)` is equivalent to H_fs identically zero there.
A nonzero R rejects that form somewhere in its box; a nonzero H_fs(0)
rejects local additivity at the origin. One zero R, or one zero H_fs(0),
does not establish additivity.

This is a statement in the specified physical source axes. For a nonlinear
coordinate change t=t(z), the response Hessian transforms as

\[
 \partial_{z_a}\partial_{z_b}\Phi
 =\sum_{ij}(\partial_{z_a}t_i)H_{ij}(\partial_{z_b}t_j)
   +\sum_i(\partial_{t_i}\Phi)\partial_{z_a}\partial_{z_b}t_i.
\]

The second term can create apparent curvature from first response alone;
even a linear mixing of source axes mixes diagonal and mixed curvature.
Separate smooth reparameterizations `t_f=a(z_f), t_s=b(z_s)` do preserve
the zero/nonzero mixed criterion locally through the factor a'b', with
nonzero derivatives. Arbitrary mixing does not. No such reparameterization
is performed in this analysis.

Nor does mixed response require two effective modes. The abstract finite
class with three labels, `L_f=L_s=F=(0,0,1)` and pi=1 has

```
Phi(t_f,t_s)=exp(t_f+t_s)/(2+exp(t_f+t_s)),
H_fs(0)=2/27 != 0.
```

It depends on only one collective coordinate and has rank-one source
variation. This is a finite-probability counterexample, not a new lattice
configuration or production sample. Nonlinear collective response can be
nonadditive in the named axes without implying independent fields.

## Which symmetries would actually force the mixed response to vanish?

**Orientation exchange.** Swapping the two labeled geometries exchanges
t_f,t_s and physical receiver labels. The mixed derivative is unchanged
as a source tensor component. S=(F_f+F_s)/2 is exchange-even; the raw
half-difference D0 is odd. Thus H_fs S is even and H_fs D0 is odd under
that relabeling. Only a separately established exchange-invariant ensemble
would force the latter average to vanish. A fixed ordered geometry pair
does not provide that assumption. For the actual normalized difference
`D=(F_f-F_s)/delta_cos4`, the denominator also changes sign, so D and its
mixed response are exchange-even: exchange does not force them to zero.

**Independent source reflections, if they existed.** Suppose an actual
measure-preserving involution preserved each class, flipped its centered
L_f and left centered L_s unchanged, with a covariant suffix map. Then
the source covariance is zero by oddness, T_fs is odd under that
involution, and any even observer would have H_fs F=0. Its finite response
would also be even in t_f, forcing every centered rectangle contrast to
vanish. With two independent reflections, a mixed response is permitted
only for an observer odd under both. A simultaneous flip of both sources
instead leaves T_fs even and gives different restrictions.

No such independent reflections are supplied by the current frozen family.
Centering L_i, using a positive exponential tilt, or evaluating both signs
of a parameter does not construct them. They cannot be used to assign
zeros to A/E or to any physical receiver in the current data.

Scientific card: this note fixes the exact probability meaning of the
second-response tensor and existing four-corner contrast, their00 support,
and their coordinate/symmetry limits. It follows allocation6bace935 and
Issue334's bounded mixed-curvature work. No finite weights or trajectories
were recomputed, no literature or test campaign was run, and no PR/comment
was made. A single covariance join should retain the original20 batches
and old8/new64 dependence when comparing the finite box with the local
Hessian; differing values are nonlinear response information, not an
estimator-independence claim.
