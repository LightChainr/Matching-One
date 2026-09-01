# Exact finite landing transfer for the root-conditioned mixed Hessian

## Result

This note fixes the finite matrix whose minors test the ordinary-block lemma
in Issue #537.  The operations must occur in this order:

1. compute the global root, means, `R`, and Schur coefficients in each full
   geometry;
2. disintegrate the complete Schur bilinear into exact two-state fibres at
   the thermal site `z`;
3. group fibres by source and thermal landing labels and C4 orbits;
4. combine separately normalized axis and tilted matrices with P4;
5. only then take `2 x 2` minors.

A nonzero minor of the final matrix disproves the finite statement that the
ordinary-four-arm block is a pure thermal reparametrization.  A minor of
`(D_z a,D_z H)` before these steps is only a raw support check.

## 1. Global quantities and exact source columns

Let `g` be axis or tilted, with pool weight `1/2`, and put

\[
 c_{axis}={1\over\Delta},\qquad c_{tilted}=-{1\over\Delta},
 \qquad y_g=2c_gE.                                            \tag{1}
\]

Then `1/2 sum_g E_g y_g` is the declared P4 mean `Y`.  Decompose the one
fixed physical source into components

\[
                              a=\sum_\lambda a^\lambda,       \tag{2}
\]

with all normalization weights absorbed into `a^lambda`.  For the canonical
spatial source, `lambda` may be an ordered pair `(x,y)` carrying its `N^-2`
weight.

Translation reduction does not change that normalization silently.  If one
fixes `x=0` and sums over displacements `y`, the physical component is
`a^y=N^-1 g_0y`; the other factor `N` has already been supplied by the
translation orbit.  There is no second `N^-1` attached to the thermal site
`z`.

At the common finite root `p0`, define

\[
\begin{aligned}
 M&={1\over2}\sum_gE_gq=0,\\
 M_t&={1\over2}\sum_g\operatorname {Cov}_g(q,S),\\
 Y_t&={1\over2}\sum_g\operatorname {Cov}_g(y_g,S),
       &R&={Y_t\over M_t},\\
 jM^\lambda&={1\over2}\sum_g\operatorname {Cov}_g(q,a^\lambda),
       &\beta_\lambda&={jM^\lambda\over M_t},\\
 H_g&=y_g-Rq,
       &\beta&=\sum_\lambda\beta_\lambda.
\end{aligned}                                                \tag{3}
\]

Here `t=logit(p)` and `S=K-Np0`.  The geometry-specific means are

\[
 \mu_{a^\lambda,g}=E_ga^\lambda,qquad
 \mu_{H,g}=E_gH_g.                                           \tag{4}
\]

All quantities in (3)--(4) are global background data; none is re-estimated
inside a landing cell.  Linearity gives the exact allocation

\[
 (a-Ea)S-\beta B
 =\sum_\lambda\{(a^\lambda-Ea^\lambda)S-\beta_\lambda B\},
 \quad B=S^2-Np_0(1-p_0).                                   \tag{5}
\]

Thus the target is

\[
 T_t={1\over2}\sum_g\sum_\lambda E_g\left[
 (H_g-\mu_{H,g})
 \{(a^\lambda-\mu_{a^\lambda,g})S-\beta_\lambda B\}
 \right].                                                    \tag{6}
\]

Equation (5) is important: assigning the full `-beta B` independently to
every pair column would overcount the root/slope projection.

## 2. Exact one-site allocation and midpoint fibre

Put `u_z=X_z-p0`, so `S=sum_z u_z`.  The exact symmetric allocation

\[
 B=\sum_zb_z,qquad b_z=u_zS-p_0(1-p_0)                     \tag{7}
\]

follows from `sum_z u_z S=S^2`.  It keeps the second thermal score inside
the same landing matrix without selecting an artificial partner site.

Fix all sites except `z` and call that background `eta`.  With

\[
 K_-=K(\eta),\qquad S_-=K_--(N-1)p_0,                        \tag{8}
\]

use `i=0,1` for `X_z` and define

\[
\begin{aligned}
 w_0&=1-p_0,&w_1&=p_0,\\
 u_i&=i-p_0,&S_i&=S_-+u_i,\\
 b_i&=u_iS_i-p_0(1-p_0),\\
 \widetilde H_i&=2c_gE_i-Rq_i-\mu_{H,g},\\
 A_i^\lambda&=a_i^\lambda-\mu_{a^\lambda,g}.
\end{aligned}                                                \tag{9}
\]

The exact contribution assigned to `(lambda,z)` in this fibre is

\[
 \Phi_{g,\lambda}(\eta,z)=
 \sum_{i=0}^1w_i\widetilde H_i
 \{A_i^\lambda u_i-\beta_\lambda b_i\}.                    \tag{10}
\]

Summing (10) over `lambda`, the off-`z` law, `z`, and `g/2` recovers (6).
Its source part has the exact midpoint form

\[
 \sum_iw_i\widetilde H_iA_i^\lambda u_i
 =p_0(1-p_0)\left\{
 \widetilde H^{mid}D_za^\lambda
 +(a^{\lambda,mid}-\mu_{a^\lambda,g})D_z\widetilde H
 \right\}.                                                   \tag{11}
\]

Thus both kernel reconnection and readout pivotality are present.  The Schur
piece

\[
             -\beta_\lambda\sum_iw_i\widetilde H_i b_i       \tag{12}
\]

must remain in the same fibre.  Equal `K_-` makes `b_0,b_1` common but does
not cancel (12), because `H_0,H_1` can differ.

## 3. Rows, columns, and geometry matrices

Choose fixed landing cuts around the marked pair and thermal site.  A label
scheme is admissible only if it is determined from the two-state fibre, C4
covariant, and identical in both geometries.

* A **source column** `alpha` contains `lambda` (or a symmetry class of
  source indices), the pair convention, the eight source-port landing
  partition, and its `z=0 -> 1` transition.
* A **thermal/readout row** `tau` contains the alternating four-arm landing,
  `(rank_0,rank_1)`, and the occupied/matching landing identifications needed
  for `q_i,E_i`.

The source label set includes an explicit absent/endpoint-occupied state.
Such a fibre can have `a_0^lambda=a_1^lambda=0` and still contributes through
`-beta_lambda b_i` in (10); dropping it changes the Schur Hessian.  Likewise,
an `ordinary-four-arm` row name is licensed only when the label records two
distinct occupied branches and two distinct vacant matching separators at
the landing cut, together with the off-port extra-contact flag.  Local
occupied/vacant alternation by itself is only a near-block label.

Let `Omega_g(tau,alpha;z)` contain the ordinary-four-arm off-`z`
backgrounds with these labels, and let

\[
 \nu_{g,-z}(\eta)=p_0^{K_-}(1-p_0)^{N-1-K_-}.                \tag{13}
\]

The complete geometry matrix is

\[
 L_g[\tau,\alpha]
 =\sum_{(\lambda,z):\lambda\in\alpha}
  \sum_{\eta\in\Omega_g(\tau,\alpha;z)}
  \nu_{g,-z}(\eta)\Phi_{g,\lambda}(\eta,z).                 \tag{14}
\]

Rows are thermal/readout landing states and columns are source landing
states.  Every entry is already a scalar contribution to the full Schur
Hessian, not a vector of jumps.  Pair normalization belongs in
`a^lambda`; pair multiplicity belongs in the sums in (14).

## 4. C4 and P4 are distinct projections

Canonicalize labels under simultaneous quarter-turns, or form

\[
 L_g^{C4}[\tau,\alpha]={1\over4}\sum_{j=0}^3
 L_g[\mathcal R^j\tau,\mathcal R^j\alpha].                   \tag{15}
\]

For spin four the quarter-turn character is `+1`.  This local orbit average
does not perform the axis-minus-tilted projection.

If an implementation fixes one nearest-neighbour direction for `z`, it has
computed only one representative of a four-element physical orbit.  It must
restore the four directions before (15) (equivalently use the correctly
normalized orbit sum).  This multiplicity is independent of the `N^-1`
fixed-origin source weight above.

The final matrix is

\[
 L^{P4,Schur}={1\over2}
 \left(L_{axis}^{C4}+L_{tilted}^{C4}\right),                 \tag{16}
\]

because the P4 signs and `Delta^-1` are already in `y_g=2c_gE`, while
`-Rq` and the Schur columns use their pooled coefficients.  Equation (16)
requires separately normalized geometry laws and the same global
`p0,R,beta_lambda`; an axis C4 orbit cannot substitute for the tilted term.

Summing all entries of (16) recovers the ordinary-four-arm portion of `T_t`.
Labels requiring extra arms belong to a separate remainder matrix.

## 5. The finite rank test

The precise finite pure-thermal lemma is the factorization

\[
 L^{P4,Schur}[\tau,\alpha]=u_\tau v_\alpha
 \quad\hbox{on the ordinary-four-arm label set}.              \tag{17}
\]

It implies rank at most one.  Hence any two row labels `tau_1,tau_2` and
two column labels `alpha_1,alpha_2` with

\[
 \boxed{\det
 \begin{pmatrix}
 L[\tau_1,\alpha_1]&L[\tau_1,\alpha_2]\\
 L[\tau_2,\alpha_1]&L[\tau_2,\alpha_2]
 \end{pmatrix}\ne0}                                         \tag{18}
\]

for `L=L^{P4,Schur}` kill that lemma.  The determinant is taken only after
(10)--(16), with signed entries preserved.

The following do not suffice: a `(D_za,D_zH)` minor; a single-geometry
minor; a C4 orbit without the tilted companion; fibre-specific `R` or
`beta`; or channelwise absolute values.  Conversely, vanishing of all
finite minors is only the algebraic gate, not the asymptotic cancellation
theorem.

## Minimal computation contract

A legitimate four-entry certificate retains:

1. both geometry IDs and their P4 coefficients;
2. exact/certified `p0,mu_a^lambda,g,mu_H,g,R,beta_lambda` from the complete
   finite laws;
3. for every fibre,
   `(lambda,K_-,q_0,q_1,E_0,E_1,a^lambda_0,a^lambda_1)` and its weight;
4. shared C4-canonical `tau,alpha` labels;
5. separate axis and tilted contributions to the four final entries;
6. the exact determinant before decimal rendering.

This is finite algebra; it requires neither an arm-probability estimate nor
a new Monte Carlo stream.
