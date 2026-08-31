# The intrinsic norm-4 readout as a thermal-clock source quotient

This note defines a physical source-response question on the archived P40
N65/N85 configurations. It does not claim a completed norm-4 transport test:
the original lineages are 65→130→260 and 85→170→340. The primary source is
the same explicit density in both directions,

\[
S=(C_B+C_W)/N.
\]

The new object below is a response of the **original intrinsic readout**,
not a replacement of that readout by a covariance with a similar name.
The source response is at \(\lambda=0\); it is not a continuum field identity.

## 1. The fixed readout and its source derivative

At a given N, write g=f,s for the stored first and second orientations and
define

\[
\mathcal P X=(X_f-X_s)/\Delta\cos4\theta,\qquad
\overline X=(X_f+X_s)/2.
\]

Let \(A_g(p)=\langle q\rangle_g\),
\(e_g(p)=\langle E\rangle_g\), and \(E=q^2\). The pooled matching root
\(p_0\) satisfies \(\overline A(p_0)=0\). With every quantity evaluated there,

\[
d=\overline{A'},\qquad B=\tfrac12\mathcal P e',\qquad
U=N^{13/8}B/d.
\]

This is the frozen U in
[`norm4_two_generator_transfer_20260829.yaml`](../predictions/norm4_two_generator_transfer_20260829.yaml),
introduced at `2236d36c80c8a466d9317c929bc33e92a7ca9d33`.
The legacy parity coordinate is \(S_{\rm legacy}=1-e/2\), and its stored
projector uses second−first. Thus the legacy `P4_S_prime` equals
\(\mathcal P e'/2\), with the positive sign used here. Do not confuse
\(S_{\rm legacy}\) with the absolute cluster source S.

For the positive measure

\[
P_{p,\lambda,g}(\eta)\propto P_{p,g}(\eta)e^{\lambda S_g(\eta)},
\]

write \(J_{f,g}=\operatorname{Cov}_g(f,S_g)\). All source functions and any
declared counterterm coefficients are held fixed when taking p or lambda
derivatives. In the Bernoulli family,

\[
J'_{f,g}=\frac{\kappa_g(f,S_g,K)}{p(1-p)}.
\]

The moving-root derivative and the full readout derivative are

\[
t=\dot p_0=-\overline{J_q}/d,
\qquad
\dot B=\tfrac12\mathcal P(J'_E+t e''),
\qquad
\dot d=\overline{J'_q}+t\overline{A''},
\]
\[
\boxed{\ \mathcal L_N[S]:=\dot U
=N^{13/8}\left(\frac{\dot B}{d}-\frac{B\dot d}{d^2}\right).\ }
\]

The baseline coefficients do not depend on S, while J, J′ and t are linear
in S. Consequently \(\mathcal L_N\) is a linear functional of the pair of
source functions. A simple pooled root and nonzero d are assumed.

## 2. Exact common-clock invariance

For a pure source \(S=a+bK\), with a,b allowed to depend on N but b common
to the two directions,

\[
P_p(\eta)e^{\lambda(a+bK)}\propto P_{\widetilde p}(\eta),
\qquad
\operatorname{logit}\widetilde p=\operatorname{logit}p+b\lambda.
\]

Re-finding the pooled matching root holds \(\widetilde p=p_0\). Both B and d
acquire the same Jacobian \(\partial_p\widetilde p\), which cancels in B/d.
Therefore

\[
\mathcal L_N[1]=\mathcal L_N[K]=0,\qquad
\boxed{\ \mathcal L_N[S+a+bK]=\mathcal L_N[S].\ }
\]

The pure-clock identity even holds at finite lambda wherever the transformed
root is defined; only its first derivative is needed here. By linearity the
boxed statement holds for any admissible microscopic source. It is a readout
on source functions modulo the **common Bernoulli thermal clock**.

This does not quotient out every occupation-dependent function: K², an Euler
polynomial, and separate coefficients \(b_fK,b_sK\) are not generally null.
Per-geometry fitted counterterms must not be described as one common physical
source. Additive constants are invisible through normalization.

## 3. What changes relative to the fixed-matching covariance

The single-p quantity

\[
C_{\rm even\mid q}=J_E-
\frac{\operatorname{Cov}(E,q)}{\operatorname{Var}(q)}J_q
\]

keeps the matching mean fixed by adding a **relative q-fugacity**. It asks
whether a source has a topology-even marginal tangent beyond that q tangent.
It does not remove the common Bernoulli thermal clock. In contrast,
\(\mathcal L_N[S]\) compares the angular even thermal slope with the pooled
matching slope after the root has moved, and removes the common clock exactly.

Thus a resolved nonzero raw \(\dot U\) would directly reject the explanation
that this source merely shifts the Bernoulli thermal coordinate. This is a
stronger discriminator for that specific mechanism, not a logical ordering
of the two statistics: either statistic can vanish while the other responds.
The thermal derivatives J′ carry information absent from C at one p.

A nonzero \(\dot U\) alone does not exclude a q-only source, identify an energy
operator, establish a Jordan module, or establish norm-4 transport. Its kernel
also contains nonthermal sources: one scalar readout can miss other sectors
or contain cancellations. A null does not identify the whole source with K.

## 4. One next physical quantity: a within-topology thermal-jet response

The K-stratified archive permits a more demanding, precisely defined
secondary question without a new sample or S² mark:

> After removing a common clock and every source term determined by the
> current three-state topology, does the absolute cluster source still change
> the intrinsic angular thermal slope?

All conditional expectations below are taken at the baseline pooled root,
separately in each geometry, and then **frozen as functions of q**. Define

\[
K_{\perp,g}=K-\mathbb E_g[K\mid q],\qquad
S_{\perp,g}=S-\mathbb E_g[S\mid q],
\]
\[
b_*=
\frac{\sum_g\mathbb E_g[S_{\perp,g}K_{\perp,g}]}
     {\sum_g\mathbb E_g[K_{\perp,g}^2]},\qquad
R_g=S_{\perp,g}-b_*K_{\perp,g}.
\]

The coefficient b* is common to both geometries; equal geometry weights are
part of this definition. If the denominator is exactly zero, K is already
determined by q on the support and the redundant clock term can be omitted.
Then

\[
\mathbb E_g[R_g\mid q]=0,\qquad
\sum_g\operatorname{Cov}_g(R_g,K)=0.
\]

Moreover R is unchanged under
\(S_g\mapsto S_g+bK+f_g(q)\), where b is common and f_g is any fixed
three-state function. This construction therefore defines a second, explicit
source quotient. It is not a change to the primary raw source.

The proposed output is simply

\[
\boxed{\ W_N=\mathcal L_N[R].\ }
\]

Although \(J_q(R)=J_E(R)=0\) at the baseline by construction, W need not vanish:
the frozen compensators do not force their p derivatives to vanish. Put
\(p_{g,r}=P_g(q=r)\) and

\[
\gamma_{g,r}=\operatorname{Cov}_g(S,K\mid q=r)
              -b_*\operatorname{Var}_g(K\mid q=r).
\]

Then the exact derivative formula is

\[
J'_f(R)=\frac{1}{p_0(1-p_0)}
 \sum_{r=-1}^1 p_{g,r}\bigl(f(r)-\langle f\rangle_g\bigr)\gamma_{g,r},
\qquad f=q,E.
\]

Here \(t(R)=0\), so

\[
W_N=\frac{N^{13/8}}d
\left[\tfrac12\mathcal P J'_E(R)-\frac Bd\overline{J'_q(R)}\right].
\]

A resolved W would locate an observed response in the coupling between
cluster structure and occupancy **within topological sectors**, beyond a
common clock plus any single-point q-only source. It would not prove a
new continuum field or escape every function of the joint variables (q,K).
A null would constrain this particular thermal-jet readout, not prove source
closure. Report raw \(\mathcal L_N[S]\) first and W as the declared compensated
diagnostic; do not choose among compensators after seeing their scores.

## 5. Exact input sufficiency and finite-estimator boundary

For each batch, orientation and K, the proposed records
`count,sum_q,sum_e,sum_s,sum_qs,sum_es` suffice. Here s is the integer cluster
sum and must be divided by N. Sector counts and source sums are recovered as

\[
n_\pm=(\Sigma e\pm\Sigma q)/2,\quad n_0=n-\Sigma e,
\]
\[
s_\pm=(\Sigma es\pm\Sigma qs)/2,\quad s_0=\Sigma s-\Sigma es.
\]

Multiplying these records by K or K² yields every conditional moment above.
S² is not needed for \(\mathcal L_N[S]\), b*, W, or their aligned batch
delete-one uncertainty. It would be needed to report the residual source
variance or a source-variance-normalized efficiency, which are not proposed
outputs here. Unobserved topology sectors cannot supply conditional means.

With reference probability \(p_*\), use

\[
w_K(p)=(p/p_*)^K[(1-p)/(1-p_*)]^{N-K},\qquad
\langle f\rangle_p=\frac{\sum_iw_{K_i}(p)f_i}{\sum_iw_{K_i}(p)}.
\]

Every derivative includes this denominator. The common-clock invariance is
also exact for this finite empirical family because multiplying its weights
by \(e^{\lambda bK}\) has the same log-odds reparametrization. For W, estimate
and freeze the conditional functions and b* at each full/leave-one-out root;
do not differentiate freshly refitted conditional functions as p varies.

These are self-normalized importance estimates of the same underlying
Bernoulli observables. They are not identical to the old threshold-integrated
estimator and are not finite-sample unbiased. General finite-lambda reweighting
cannot be recovered from these first source moments: exponentiating a
conditional mean s would change the measure. Only the lambda-zero tangent
is claimed. Root location/slope, effective sample size and weight concentration
should accompany the output; no new production or repeated test suite is
implied by this note.

## 6. Complete-lineage prediction: source rigidity versus a moving generator

### The source coordinate must agree across N

The current source density S=s/N, where s=CB+CW, uses the tilt
\(e^{\lambda s/N}\). A common lambda across sizes therefore means the
size-dependent cluster fugacity \(Q_N=e^{\lambda/N}\). For a common
microscopic absolute-cluster fugacity \(Q=e^t\), the measure instead contains
\(e^{t s}\), and

\[
u_N=\partial_\lambda U_N\big|_{\lambda=0},\qquad
\boxed{\ v_N=\partial_t U_N\big|_{t=0}=N u_N.\ }
\]

The same conversion applies to W by linearity. Both coordinates are useful,
but a density-source residual cannot be presented as a common-Q residual.
This Q weights the total black-NN plus white-matching cluster count; it is
not an unproved identification with a different Potts/FK source interface.

### Two falsifiable fixed-coefficient source hypotheses

On the original two lineages, with each pooled matching root relocated under
the same bulk coupling t, define

\[
R_{q2,N}(t)=U_N(t)-3U_{2N}(t)+2U_{4N}(t),
\qquad
R_{J,N}(t)=U_N(t)-2U_{2N}(t)+U_{4N}(t).
\]

The explicit **source-rigidity** hypotheses are, separately,

\[
H_{q2}^{\rm src}:\quad v_N-3v_{2N}+2v_{4N}=0,
\]
\[
H_J^{\rm src}:\quad v_N-2v_{2N}+v_{4N}=0,
\qquad N=65,85.
\]

These hold if the corresponding transfer coefficients remain fixed and the
closure residual has no first-order response to the same microscopic source.
A sufficient, stronger family assumption would be
\(U_N(t)=A(t)+C(t)/N\) or \(U_N(t)=A(t)+B(t)\log N\) throughout a neighborhood
of t=0. Neither assumption follows from an unperturbed fit. In particular,
compatibility of a source-tangent law cannot rehabilitate a rejected t=0 law:
\(R'(0)=0\) is not \(R(0)=0\). The tangent hypotheses can also be posed with a
nonzero but source-rigid baseline closure defect, stated explicitly.

Using the reported density derivatives, the correct common-Q residuals are

\[
R'_{q2,N}(0)=N\,[u_N-6u_{2N}+8u_{4N}],\qquad
R'_{J,N}(0)=N\,[u_N-4u_{2N}+4u_{4N}].
\]

This is a direct two-lineage physics output from source marks on the six
existing-size permutation archives. It requires the joint uncertainty of
the actual reused streams; neither a fresh independent sample nor a claim
that the old law automatically transports the source is implied.

### A sharper alternative with one common moving generator

Write the normalized recurrence as

\[
U_{4N}(t)-[1+\kappa(t)]U_{2N}(t)+\kappa(t)U_N(t)=\epsilon_N(t),
\]

where \(\kappa_0=1/2\) gives q2 and \(\kappa_0=1\) gives Jordan. A concrete
next hypothesis allows one **common** \(\dot\kappa\) across both lineages while
holding the baseline defect source-rigid, \(\epsilon'_N(0)=0\). Differentiation
then predicts

\[
r_N:=v_{4N}-(1+\kappa_0)v_{2N}+\kappa_0v_N
=\dot\kappa\,(U_{2N}-U_N).
\]

Use this normalized r in the determinant: for q2,
\(r_N=\tfrac12R'_{q2,N}(0)\); for Jordan,
\(r_N=R'_{J,N}(0)\). The reported legacy q2 residual is twice the normalized
recurrence residual, and its covariance scales accordingly. If the baseline
law is not satisfied, the source calculation is a diagnostic of its proposed
rigid extension (with the stated baseline defect), not certification of the
original model.

Thus the two response residuals must be parallel to the two unperturbed
increments. One directly scoreable, gain-free contrast is

\[
\boxed{\ D_{\rm src}
=r_{65}(U_{170}-U_{85})-r_{85}(U_{130}-U_{65}).\ }
\]

All quantities can come from the same K-stratified estimates and their joint
aligned-delete-one vectors; do not insert independent-source error bars or
divide by an unresolved increment. The determinant is a derived uncertain
quantity, not an exact probability certificate.

This yields three concrete readings: unresolved r leaves source rigidity
compatible; resolved r consistent with one common slope supports that
particular moving-generator description; resolved non-collinearity rejects
one common generator drift with a source-rigid defect. The latter does not
by itself count fields: source-dependent defects, additional finite-size
terms, or a richer transfer can also produce it. No new exponent, new
counterterm election, or general framework is needed to state this test.
