# Root-comoving rank-1 occupancy under an absolute cluster source

The proposed readout asks a direct physical question: **after the common
Bernoulli parameter has moved enough to restore matching balance, does cluster
fugacity change the probability of having exactly one homological direction?**
It uses source covariances, not an additional thermal derivative, and is a
distinct finite-volume observable from the angular thermal-slope response Udot.
Lower derivative order is a reason to measure it, not a promise of lower variance.

## 1. Geometry, source and normalization

For each of the two fixed geometries g, let h be black homology rank and define

\[
q=h-1,\quad E=q^2,\quad P_{r,g}=P_g(h=r),\quad
A_g=\langle q\rangle_g,\quad e_g=\langle E\rangle_g.
\]

Use a bar for the equal average of the **separately normalized** geometry
expectations. The common microscopic paired-cluster fugacity is
\(Q=e^t\), with measure proportional to
\(P_{p,g}(\eta)e^{t(C_B+C_W)}\). If the stored source is
\(S=(C_B+C_W)/N\), put

\[
J_{q,g}=\operatorname{Cov}_g(q,S),\qquad
J_{E,g}=\operatorname{Cov}_g(E,S),\qquad
D=\overline{\partial_p A_g}.
\]

The pooled matching root satisfies \(\overline A(p_0(t),t)=0\). At t=0,

\[
\dot p_0=-N\overline{J_q}/D.
\]

The required mean covariance is \(\overline{\operatorname{Cov}_g(f,S)}\),
not the covariance formed after pooling the geometry populations. Indeed,

\[
\operatorname{Cov}_{\rm mixture}(f,S)
=\overline{\operatorname{Cov}_g(f,S)}
+\tfrac14(\langle f\rangle_f-\langle f\rangle_s)
          (\langle S\rangle_f-\langle S\rangle_s).
\]

The extra term would change the mixture weights under the source. The
experiment instead keeps the two geometries equally weighted.

## 2. Root motion and rank composition are different responses

Since \(\overline P_1=1-\overline e\), the comoving rank-1 response is

\[
\boxed{\ \dot P_1^{\rm root}
=-N\left[\overline{J_E}
-\frac{\overline{E_p}\,\overline{J_q}}D\right].\ }
\]

Here \(\overline{E_p}=\overline{\partial_p\langle E\rangle_g}\). The first
term changes the rank composition at fixed p; the second accounts for moving
the common Bernoulli parameter to the new matching root. This is not the
fixed-q-fugacity projection using \(\operatorname{Cov}(E,q)/\operatorname{Var}(q)\).

At the pooled root, \(\overline P_0=\overline P_2\), throughout the source
family. Consequently

\[
\dot{\overline P}_0^{\rm root}
=\dot{\overline P}_2^{\rm root}
=-\tfrac12\dot P_1^{\rm root}.
\]

A positive rank1dot transfers probability from both extreme sectors into
rank 1; a negative value transfers it in the other direction. These
equalities are for direction averages. A pooled matching root does not
require \(P_{0,g}=P_{2,g}\) separately in each geometry.

For a common thermal source \(S=a+bK\), Bernoulli differentiation gives
\(J_f=b\,p(1-p)\partial_p\langle f\rangle\). Therefore

\[
\dot p_0=-Nb\,p_0(1-p_0),\qquad
\boxed{\ \dot P_1^{\rm root}=0.\ }
\]

A resolved nonzero rank1dot would thus exclude a source response that is
only a common affine-K thermal-clock shift. A zero does not prove that model:
this one readout can miss other deformations. K², geometry-dependent clocks,
or q-only source functions are not all removed by this null. The response
remains within the three-state topology algebra; E is not thereby identified
with a continuum energy operator.

## 3. First/second cumulative activation: an exact root budget

Use distinct notation for the cumulative indicators

\[
H_1=\mathbf1\{h\ge1\}=1-(E-q)/2,\qquad
H_2=\mathbf1\{h=2\}=(E+q)/2.
\]

Their fixed-p source responses are \((J_q-J_E)/2\) and
\((J_q+J_E)/2\). Hence the root-displacement allocation is

\[
r_1=-\frac{N(\overline{J_q}-\overline{J_E})}{2D},\qquad
r_2=-\frac{N(\overline{J_q}+\overline{J_E})}{2D},\qquad
r_1+r_2=\dot p_0.
\]

These are additive terms in a **static susceptibility budget**. They are
not the comoving changes of H1/H2. The latter obey

\[
\dot{\overline H}_1^{\rm root}=\tfrac12\dot P_1^{\rm root},\qquad
\dot{\overline H}_2^{\rm root}=-\tfrac12\dot P_1^{\rm root},\qquad
\dot{\overline H}_1^{\rm root}+\dot{\overline H}_2^{\rm root}=0.
\]

If a scorer calls these cumulative events `I1` and `I2`, this is the meaning
of `movingI1+movingI2=0`. It would be incorrect for I1/I2 interpreted as the
two mutually exclusive rank-1/rank-2 indicators.

Neither r1/r2 nor a rank-sector allocation measures the source at the instant
of a particular K1/K2 birth. The same-K archive gives \(\Delta(SH_i)\), but
separating \(S_K\Delta H_i\) from \(H_i(K+1)\Delta S\) requires cross-step
marks such as \(S_Kq_{K+1}\) and \(S_KE_{K+1}\), which are not currently stored.

## 4. Mechanism conjecture: sector population can move without angular deformation

A strong rank1dot together with weak Udot would motivate the conjecture that
absolute cluster fugacity predominantly changes the **population of the
one-direction sector**, while the conditional spatial organization responsible
for the global H4 thermal slope changes little or cancels in that observer.
This is a hypothesis, not an inference from two significance labels; the
quantities have different units and their uncertainty must accompany comparison.

The spatial discriminator is the **unoriented primitive winding line within
rank 1**. If its generator has physical displacement
\(v=m\omega_1+n\omega_2\), use
\(O_4=(v/|v|)^4\), which is unchanged by reversing the generator. Define it
in the same physical lattice frame across geometries. The comoving source
response of \(\mathbb E[O_4\mid h=1]\) separates a population-only change from
conditional winding-line rearrangement. A change of this conditional harmonic
with weak global Udot would expose observer cancellation rather than establish
the absence of angular structure. Current q/E/cluster-count moments do not
contain this winding-line mark, so no such spatial result is asserted here.

## Appendix: existing K profiles already resolve state-support contributions

For this appendix use \(r=q\in\{-1,0,1\}\). For the previously defined W
quotient, let \(n_{gkr},s_{gkr}\) be the sector
count and sum of source density at K=k, and let \(m^K_{gr},m^S_{gr},b_*\) be
its fixed-root conditional means and common clock coefficient. With T
permutations and normalized Binomial weight w_k,

\[
\xi_{gkr}=\frac{w_k(k-m^K_{gr})}{T p_0(1-p_0)}
\left[s_{gkr}-n_{gkr}\{m^S_{gr}+b_*(k-m^K_{gr})\}\right],
\]
\[
W=\sum_{g,k,r}\frac N D
\left[\frac{N^{13/8}\sigma_g}{2\Delta\cos4}(r^2-e_g)
-\frac U2(r-A_g)\right]\xi_{gkr},\qquad \sigma_f=1,\ \sigma_s=-1.
\]

Grouping this identity by the three actual q states and three fixed K bands
gives a nine-cell support description without redefining W. The bands can
use pooled complement K1/K2 medians, frozen independently of the source means.
Complete paired batch trajectories supply the cell covariance V; the signed
noise allocation is \(\nu_a=\sum_bV_{ab}=\operatorname{Cov}(W_a,W)\), whose sum
is Var(W). Negative allocations describe covariance cancellation, not negative
variances. These are support and mean-influence descriptions, not individual
birth-event attribution or independent evidence rows. The quotient derivation
is in [the thermal-clock source note](p40-thermal-clock-source-quotient.md).
