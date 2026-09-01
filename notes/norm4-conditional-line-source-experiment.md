# A conditional winding-line experiment after the resolved rank-population response

The physical question is now **whether absolute cluster fugacity rearranges
winding inside rank1 beyond the change predicted by a common topological
source and thermal clock**. The existing root-comoving rank1 response already
resolves a population change. This experiment adds a spatial readout; it does
not repeat that marginal measurement or assume E_top is an energy operator.

## Readout and the moving-root prediction

At each fixed geometry g, use black homology rank h, q=h−1, E=q² and R=1{h=1}.
For its unoriented primitive winding vector in physical lattice coordinates,
v=mω1+nω2, define O=(v_x+i v_y)^4/|v|^4. Reversing v leaves O unchanged.
Use the same physical frame in both geometries, not the bare integer pair (m,n).
The source is the **bulk** count s=C_B+C_W, conjugate to t=log Q; a stored
density-source derivative must first be multiplied by N.

All following conditional quantities are formed **inside each geometry first**:

\[
\mu_g=\frac{\langle RO\rangle_g}{\langle R\rangle_g},\qquad
\mu_{p,g}=\frac{\partial_p\langle RO\rangle_g
 -\mu_g\partial_p\langle R\rangle_g}{\langle R\rangle_g},\qquad
C_{s,g}=\frac{\langle ROs\rangle_g-\mu_g\langle Rs\rangle_g}
 {\langle R\rangle_g}.
\]

Bars denote fixed equal averages of separately normalized geometries. At the
pooled matching root, let D=mean_g ∂p〈q〉_g. Then

\[
\dot p_s=-\overline{\operatorname{Cov}_g(q,s)}/D,\qquad
\nu_{s,g}=C_{s,g}+\dot p_s\mu_{p,g}.
\]

Conditioning after pooling the two geometries would instead introduce
source-dependent geometry weights. It is not this estimand. Full-K batch
profiles of R, RO, Rs and ROs suffice, with Binomial differentiation for the
p jets; same-batch uncertainty does not require an additional s² moment.

Consider the finite first-order **common E-plus-clock** response model,
s∼a_g+b_N K+d_N E. Constants may differ by geometry; b_N,d_N are common to
the pair. For the dimensionless E source,

\[
\dot p_E=-\overline{\operatorname{Cov}_g(q,E)}/D,\qquad
\dot r_E=-\overline{\operatorname{Var}_g(E)}
             -\dot p_E\,\overline{\partial_p\langle E\rangle_g},\qquad
d_N=\dot r_s/\dot r_E,
\]

where r=mean_g〈R〉. E vanishes inside rank1, so its fixed-p conditional
response is zero. A common K source gives C_K=p(1−p)μ_p, exactly canceled
by its matching-root displacement. The model therefore predicts

\[
\boxed{\nu_{s,g}^{(E+clock)}=d_N\dot p_E\mu_{p,g}},\qquad
e_g=\nu_{s,g}-d_N\dot p_E\mu_{p,g}.
\]

Root transport is essential: **nonzero ν_s alone does not reject E-plus-clock**.
The primary residual is the four-real vector
(Re e_first, Im e_first, Re e_second, Im e_second). Recompute the conditional
ratios, root, d_N and all prediction terms together in each aligned delete-one
view, including the declared root-anchor uncertainty. Its joint quadratic
score uses the supported covariance rank (normally four); do not subtract a
further degree of freedom for d_N, which was determined by the additional
rank-population equation. If rdot_E is poorly resolved, label this ratio
prediction unstable; its division-free equivalent is
rdot_E ν_s−rdot_s pdot_E μ_p, not evidence that the model closes.

## A stronger, division-free common-topology-plus-clock comparison

For s∼a_g+b_N K+c_N q+d_N E, both topological terms vanish conditionally on
rank1. Define four-real vectors C from C_s and T from p(1−p)μ_p, in the same
geometry/quadrature order. The necessary prediction is C=b_N T. Thus all six

\[
w_{ij}=C_iT_j-C_jT_i\quad(i<j)
\]

must vanish, independently of c_N,d_N and without the rank-population ratio.
These six numbers are **not six independent constraints**. At a regular null
C=bT with T≠0, write (L_T x)_{ij}=x_iT_j−x_jT_i. The first-order perturbation
is δw=L_T(δC−bδT), whose rank is at most three. From the saved eight-real
(C,T) covariance Σ, construct the null-tangent covariance

\[
\Omega=L_T[\Sigma_{CC}-\hat b\Sigma_{CT}
 -\hat b\Sigma_{TC}+\hat b^2\Sigma_{TT}]L_T^\mathsf T.
\]

A consistent choice is b_hat=(T·C)/(T·T). Its first-order estimation error is
parallel to T and is annihilated by L_T. Under this regular null, wᵀΩ⁺w has
a nominal chi-square comparison with the supported rank r≤3. Report r and
the component vector. A direct nonlinear delete-one wedge covariance can
acquire higher-order extra eigenvalues; counting them as five/six degrees
of freedom would misrepresent the null. When T is small or noisy, this
first-order approximation is weakly identified: report the measured C,T and
wedges without treating unresolved products as source closure. Exact known
symmetries can further reduce the supported rank.

## What the outcomes would change

- A resolved primary e rejects this named E-plus-clock prediction in the
  conditional line readout. A q contribution or another topological map may
  still explain it; this alone does not prove an irreducible spatial source.
- A resolved wedge rejects the **common** clock-plus-q/E response map here.
  A within-geometry Re/Im wedge detects a response not parallel to that
  geometry's thermal line direction. Only cross-geometry wedges changing can
  instead expose unequal effective clock coefficients, an informative but
  narrower failure of sharing.
- If rank population moves while both comparisons remain unresolved, keep
  the measured population response and spatial uncertainty separately. This
  does not prove population-only dynamics or identify a continuum field.

This tests the conjecture that strong population transfer can coexist with
a nearly rigid conditional winding distribution; a resolved conditional
residual would expose spatial rearrangement hidden by weak global Udot.
The four cyclic sizes share a counter domain, while N260 and N340 form their
own groups; paired deletions and cross-N covariance remain mandatory for
combined interpretation. Additional marks on the same old configurations,
real/imaginary components and E/wedge views are correlated observations,
not independent replications. Coefficients may depend on N: finite agreement
does not establish a universal source map, a scaling law or the original
norm-4 operator identity. This note defines the experiment; it claims no
conditional-line result before its archived-data calculation.
