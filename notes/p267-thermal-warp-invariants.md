# When is thermal redistribution a coordinate change?

## Result

There are two distinct hypotheses. A **scalar observable** changes as
`G = D composed with phi`; a **density profile** changes as
`G = phi' * (D composed with phi)`. Only the second gives the infinitesimal
ansatz `R = (v D)'` and preserves the integrated clock. Rejecting a few
polynomial velocities does not reject general finite density transport.

The sharp next question is not whether one single-sign A profile can be
warped—it always has a monotone cumulative-coordinate construction—but
whether **E is transported by the same A-defined coordinate**. This gives
an exact, root-free two-channel null:

\[
\boxed{\Omega(U)-r_A r_E\Omega(D)=0},\qquad
\Omega(D)=\int_0^1[F_A(p)D_E(p)-F_E(p)D_A(p)]\,dp,
\]

where `F_j(p)=integral_0^p D_j` and the two amplitudes may differ. A stronger
sequence uses the normalized cumulative A coordinate without dividing by
the possibly weak E mass. These are finite-warp invariants, not a new
velocity fit or a field identification.

This is an exact-theory companion to [PR #484](https://github.com/LightChainr/Matching-One/pull/484),
source commit `7b30648be558df0652a7ff22143cc87ed399d042`,
`notes/etop-clock-redistribution-n100.md`, and the clock quotient in
[PR #485](https://github.com/LightChainr/Matching-One/pull/485), science commit
`6543c0e`. It does not rerun their Monte Carlo or rescore their production.

## 1. Fix the transformation law before fitting

Let `phi:[0,1]->[0,1]` be orientation-preserving, endpoint-fixed, with
strictly positive derivative. For `phi(p)=p+epsilon v(p)+O(epsilon^2)`:

\[
\begin{array}{c|c|c}
 &G(p)&\partial_\epsilon G|_0\\
\text{scalar observable}&D(\phi(p))&vD'\\
\text{density/profile}&\phi'(p)D(\phi(p))&(vD)'
\end{array}
\]

Thus the Jacobian is a physical coupling assumption. An integrated
clock identity makes density transport a useful hypothesis, but does not
automatically turn an ordinary wrapping probability into a density.
Both preserve ordered sign/zero structure; only density transport
preserves signed mass and signed lobe areas.

For the three N100 shapes use

\[
D_j=Y_j(4i)-Y_j(2i),\quad U_j=Y_j(\tfrac12+i)-Y_j(2i).
\]

The original common-amplitude clock ansatz is
`U_j = r_C phi' D_j composed with phi`. Its normalized target is `U_j/r_C`,
not `D_j + (U_j-r_C D_j)`. A nonzero negative `r_C` reverses the raw overall
sign and is removed before comparing signed profiles. The more permissive
common-coordinate hypothesis is

\[
U_A=r_A\phi'D_A\circ\phi,\qquad U_E=r_E\phi'D_E\circ\phi.
\tag{1}
\]

Independent amplitudes test a common coordinate without already insisting
on common A/E coupling. A failure of (1) is therefore stronger than a
failure of the equal-amplitude version.

## 2. The order-free regular tangent gate

Write `R=alpha(vD)'`, `alpha != 0`, and `H(p)=integral_0^p R`. Fixed
endpoints imply `v(0)=v(1)=0`, hence

\[
v(p)=\frac{H(p)}{\alpha D(p)}.
\tag{2}
\]

Consequences for analytic profiles, including finite Bernstein polynomials:

- `integral R=0` is necessary, not sufficient.
- At an interior zero `zeta` of D of order m, H must vanish to order at
  least m. For a simple zero, `H(zeta)=0` is the condition and
  `v(zeta)=R(zeta)/(alpha D'(zeta))`.
- If `H(zeta)!=0` at a simple zero, the required velocity has pole residue
  `H(zeta)/(alpha D'(zeta))`. Increasing polynomial degree cannot remove it.
- At an endpoint zero of D of order m, H needs order at least m+1 for the
  endpoint velocity to vanish, not merely to remain bounded.
- For rational polynomial inputs, cancel the exact polynomial gcd of H
  and `alpha D`; the reduced denominator must have no roots on `[0,1]`
  and the reduced velocity must vanish at both endpoints. This is an exact
  criterion, independent of any chosen velocity-polynomial degree.

If one only considers a local interval with free boundary flux, (2)
becomes `(H+c)/(alpha D)`. In a zero-free interval, **every smooth residual
has such a representation**. Across several source zeros, all H values
at those zeros must agree with the same `-c`, together with the higher
vanishing-order conditions. One unrestricted local velocity for one curve
therefore has no identifying power.

For two channels with the same velocity and common alpha, the zero-flux
condition is

\[
H_A D_E-H_E D_A\equiv0,
\tag{3}
\]

plus the regularity conditions above. For different fixed amplitudes use
`(H_A/r_A)D_E-(H_E/r_E)D_A=0`. This is stronger than fitting each channel
with its own velocity.

The familiar moment relation follows by integration by parts. With
`z=s(p-p_ref)`, `D_j=integral z^j D`, and `v=a+b(p-p_ref)`,

\[
R_j=-\alpha a s jD_{j-1}-\alpha b jD_j.
\]

Absorbing alpha into the fit coefficients gives the existing convention.
These are tangent constraints, not finite-warp moment invariants.

### Why the tangent pole is not a finite-warp no-go

For a genuine finite density warp `G=phi' D composed with phi`,

\[
\int_0^\zeta(G-D)\,dp=F_D(\phi(\zeta))-F_D(\zeta).
\]

At a source zero zeta this is generically nonzero at second order in the
displacement. So `G-D` can fail the regular **linear** generator gate while
G is exactly a finite pullback. The certificate explicitly constructs this
case: `D=p(1-p)(p-1/3)`, `phi=p+p(1-p)/4`, and the old-zero flux is
`151/419904 != 0`. No finite-warp rejection is inferred from that pole.

## 3. Finite density transport: signed lobes, not ordinary moments

After amplitude normalization, a regular finite density warp preserves:

1. total signed mass;
2. the ordered sign pattern and finite-order zero multiplicities;
3. the signed area of every corresponding lobe;
4. consequently positive/negative masses and total variation `integral |D|`.

The total-variation gate is necessary but weaker than the full ordered
lobe vector. For finitely many isolated finite-order zeros, matching sign
order, zero multiplicities, and corresponding lobe masses permits a
piecewise cumulative-inverse construction

\[
\phi=(F_D|_{\text{corresponding lobe}})^{-1}\circ F_G.
\]

Matching zero orders gives finite positive limiting derivatives across
the matched zeros. Endpoint orders must also match if an endpoint-regular
diffeomorphism, rather than merely a monotone homeomorphism, is required.
For example, at corresponding simple zeros the derivative obeys
`phi' = sqrt(G'/D')` with source derivative evaluated at the matched zero.

Ordinary p moments need not be preserved:

\[
\int f(p)G(p)\,dp=\int f(\phi^{-1}(q))D(q)\,dq.
\]

The certificate gives equal-mass linear profiles with a nonzero TV gap
`17/72`; these cannot be finite density pullbacks of each other. Conversely,
when both amplitude-normalized profiles have one strict sign in the
interior, matching mass always yields a monotone cumulative warp. For
such a pair TV adds no information. Rejection of a low-degree tangent then
only rejects that restricted tangent family.

## 4. The A-defined clock makes the shared-warp hypothesis identifiable

Assume D_A and U_A have one strict sign in the interior and nonzero mass.
Let `M_A=integral D_A`, `M_UA=integral U_A`, and

\[
q_D(p)=F_A(p)/M_A,\quad q_U(p)=F_{U_A}(p)/M_{UA}.
\]

Both increase from 0 to 1 even if their unnormalized A profiles are
negative. Their only possible monotone coordinate transport is

\[
\boxed{\phi=q_D^{-1}\circ q_U.}
\tag{4}
\]

There is no remaining velocity fit. In this coordinate, E defines a
finite signed measure, not necessarily a probability measure:

\[
d\nu_D(q)=D_E(p)\,dp,
\quad \frac{d\nu_D}{dq}=
M_A\frac{D_E(p)}{D_A(p)},\quad p=q_D^{-1}(q).
\]

Equation (1) holds exactly iff `nu_U = r_E nu_D` under the A clock (with
the regularity qualification for (4)). In particular, define

\[
J_m(D)=\int_0^1 q_D(p)^mD_E(p)\,dp.
\]

Then

\[
\boxed{J_m(U)-r_EJ_m(D)=0,\quad m=0,1,2,\ldots .}
\tag{5}
\]

All moments are sufficient for equality of finite signed measures on
`[0,1]`, because polynomials are uniformly dense in continuous functions.
A finite selection such as `m=1,...,6` is necessary, not sufficient. It is
nonetheless a precise common-coordinate test with **no division by E's
weak mass, no estimated zero locations, and no polynomial velocity order**.
If r_E is itself estimated as the E mass ratio, retain that uncertainty;
do not normalize by U_E's weak mass or treat the fitted amplitude as fixed.

The smallest nontrivial invariant is the oriented cumulative-path area:

\[
\Omega(D)=2\int F_A D_E-M_A M_E.
\]

Under (1), `F_UA=r_A F_A composed with phi` and similarly for E, which
proves `Omega(U)=r_A r_E Omega(D)` by change of variable. If the two mass
nulls vanish,

\[
J_1(U)-r_EJ_1(D)=
\frac{\Omega(U)-r_A r_E\Omega(D)}{2r_A M_A}.
\tag{6}
\]

Thus J1 is not independent of the area test. Individually valid A and E
warps need not pass it: the exact certificate preserves each channel's
mass using different polynomial phis but produces a nonzero Omega gap.
This is a genuine shared-coordinate obstruction, not evidence for a
particular additional field or state count.

The first variation connects the finite invariant to (3):

\[
\delta\Omega=2\int(H_A D_E-H_E D_A)
+M_A H_E(1)-M_E H_A(1).
\]

For mass-preserving variations only the integrated flux determinant
remains. A zero integrated determinant is weaker than the pointwise
tangent identity (3).

## 5. What the existing full-p archives can now test

The archive can propagate its same-batch covariance or leave-one-batch-out
jackknife through area amplitudes, Omega, and a fixed vector of J moments.
No new random samples are required. Polynomial Bernstein data permit
exact integration; rational arithmetic is useful for small certificates,
while numerically stable Bernstein/Beta integration is appropriate for
the existing large-degree production histogram.

For simple stable lobe boundaries, their random displacement contributes
no first-order term to a lobe integral because the integrand is zero at
the boundary. In contrast the old-source-zero tangent flux has influence

\[
\delta H(\zeta)=\int_0^\zeta\delta R
-\frac{R(\zeta)}{D'(\zeta)}\delta D(\zeta),
\]

so treating a data-dependent source zero as fixed can miss covariance.
Weak/unresolved zeros do not justify a plug-in exact no-go. Likewise an
exact nonzero functional of empirical mean polynomials is not an
uncertainty-calibrated rejection. This note makes no new production
significance claim.

## Reproduce the exact certificates

```sh
python3 scripts/p267_thermal_warp_invariants.py \
  --json results/p267-thermal-warp-invariants/certificate.json
python3 -m unittest discover -s tests -p test_p267_thermal_warp_invariants.py -v
```

The script uses only the standard library and exact rational arithmetic.
An optional input JSON supplies arrays `D_A,D_E,U_A,U_E`, amplitudes
`r_A,r_E` (or common `alpha`), and `basis: "power"` (low-order first) or
`"bernstein"`. It computes mass nulls, Omega null, and `J_0,...,J_6` nulls.
It does not estimate covariance or choose a production hypothesis.

## Scientific card

- **Mechanism space changed:** scalar relabeling, regular infinitesimal
  density transport, single-channel finite transport, and shared A/E
  finite transport are now separate, executable hypotheses. A's cumulative
  coordinate removes arbitrary-warp freedom when A has one sign.
- **Not proved:** no N100 rejection is asserted here; no critical field,
  LCFT identity, new state number, or cross-scale law is inferred.
- **Observer / sector / geometry:** thermal p profiles of ordinary
  P4-projected A_top and E_top contrasts; N100 shapes `2i,4i,1/2+i`.
- **Source / dependency:** exact rational certificates are new algebra;
  any application to PR484 uses its existing common-stream N100 block,
  not independent evidence from PR485's clock reuse.
- **Next discriminant:** finite shared-clock `J_m(U)-r_EJ_m(D)` with the
  source and amplitude uncertainty retained. Rejection means the E profile
  needs a deformation beyond the A-defined thermal coordinate; survival
  of finitely many moments is not a proof of common transport.
