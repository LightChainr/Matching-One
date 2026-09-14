# Rare topological charge: an exact algebra for balance without concentration

2026-09-14.  Finite identities intended as an explanatory layer for the geometric-balance manuscript.  No asymptotic theorem beyond those already on the branch is needed for Sections 1--3.

## 1. Charge and neutral gas

For an occupied configuration put

\[
D=r-1\in\{-1,0,+1\}.
\]

By the wrapping-component identity this is also

\[
D=W_4-W_8.
\]

The exact projective-homology classification may be viewed as

- `D=-1`: black rank zero / white rank two, no neutral paired winding gas;
- `D=0`: black and white rank one, with a common slope and `K>=1` paired essential components;
- `D=+1`: black rank two / white rank zero, no neutral paired winding gas.

Thus `D` is an `O(1)` topological charge, whereas the rank-one count `K` is the neutral component gas.

The matching observable is simply

\[
M(p)=E_pD=P_2(p)-P_0(p).
\]

The finite matching root `p_*` is therefore the zero-mean-charge condition

\[
\boxed{E_{p_*}D=0.}
\]

It does not require the neutral gas to be small, concentrated, or absent.

## 2. Charge susceptibility

The exact charge variance is

\[
\chi(p):=\operatorname{Var}_p(D)
=P_0(p)+P_2(p)-M(p)^2.
\]

At the matching root, write

\[
P_0(p_*)=P_2(p_*)=\varepsilon_*.
\]

Then

\[
\boxed{
\chi(p_*)=2\varepsilon_*=1-P_1(p_*).}
\]

Hence on a long rank-one plateau the topological susceptibility can be exponentially small even though the finite zero-charge point remains unique.

## 3. Exact factorization of the median slope

Let

\[
H(p)=\frac{P_2(p)}{P_0(p)+P_2(p)}
\]

be the endpoint-sector conditional odds, and let

\[
F(p)=\frac12(1+M(p))
\]

be the two-birth mixture CDF.

Put

\[
L(p)=\log\frac{P_2(p)}{P_0(p)}.
\]

Then `H'=H(1-H)L'`.  At the matching root `H=1/2`, so `H'=L'/4`.

Also, because `P_0=P_2=epsilon_*`,

\[
M'(p_*)
=P_2'(p_*)-P_0'(p_*)
=\varepsilon_*L'(p_*).
\]

Using `chi=2 epsilon_*` gives the exact identities

\[
\boxed{M'(p_*)=2\chi(p_*)H'(p_*),}
\]

and therefore

\[
\boxed{F'(p_*)=\chi(p_*)H'(p_*).}
\]

This is a literal factorization of the finite median density into

1. the probability mass left in the two charged endpoint sectors, and
2. the speed at which their **conditional odds** rotate through equality.

The root theorem in the main manuscript works by controlling item 2 through a ratio comparison; full-law concentration is governed by whether item 1 stays substantial away from `p_c`.  The two questions are algebraically different already at finite size.

## 4. Bernoulli score form

For `N` sites let `X` be the occupied-site count and

\[
S_p=\frac{X-Np}{p(1-p)}
\]

be the product-Bernoulli score.  For each sector with positive probability,

\[
\frac d{dp}\log P_j(p)=E[S_p\mid r=j].
\]

Thus

\[
\boxed{
L'(p)
=E[S_p\mid r=2]-E[S_p\mid r=0].}
\]

At the root,

\[
F'(p_*)
=\frac{\chi(p_*)}{4}
\left(E[S_{p_*}\mid r=2]-E[S_{p_*}\mid r=0]\right).
\]

This representation may be useful for finite exact transfers or rare-event samplers because it says precisely which two conditional ensembles determine the condition number of the matching root.

## 5. Asymptotic interpretation in elongated geometries

In a separated two-window geometry, the interior of the interval between the births is overwhelmingly rank one.  At the matching root one therefore expects

\[
\chi(p_*)=1-P_1(p_*)\to0.
\]

The root can nevertheless remain sharply selected because the endpoint **odds ratio** `P_2/P_0` changes exponentially with geometry.  The factorization above makes this compatible with an increasingly flat mixture CDF:

- the mixture may converge to a broad two-atom law or an interior plateau rather than concentrate;
- the finite zero-charge point still converges to `p_c` under the much weaker shortest-period condition of the main root theorem.

This is the algebraic core of the phrase **balance without concentration**.

In the extreme-elongation limit the mixture approaches `1/2` throughout every fixed interior `p`, so `F'` vanishes there.  There is no contradiction with finite uniqueness: the limiting plateau loses the information carried by the exponentially rare charged sectors.

## 6. Numerical implication

Estimating the root by treating `F` as an ordinary well-conditioned median can become inefficient when `chi` is tiny.  The identity

\[
F'=\chi H'
\]

says that a tiny absolute error in `F` must be compared to a tiny derivative, whereas a targeted endpoint-sector odds calculation works directly with the object that selects the root.

This does not authorize an unweighted ratio estimator on rare sectors; it only identifies the correct conditioning.  Any computational method must retain the covariance/rare-event error model.

## 7. Source-field generating function

The exact topological source may be written

\[
Z_p(h,z;\chi_{slope})
=P_2e^h+P_0e^{-h}+P_1H_p(z;\chi_{slope}),
\]

as in `projective-homology-gas-20260914.md`.  Then `D` is the charge coupled to `h` and the paired component count/slope live in the neutral sector.  At `h=0,z=1`,

\[
\partial_hZ=M,
\qquad
\partial_h^2\log Z=\chi.
\]

No continuum/CFT interpretation is required for this source decomposition.
