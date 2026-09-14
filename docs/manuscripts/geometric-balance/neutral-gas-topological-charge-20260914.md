# Neutral component gas versus bounded topological charge

2026-09-14.  Exact finite topology plus asymptotic consequences for transfer/count descriptions.  This note explains why extensive winding-count thermodynamics and the Matching-One balance observable live at different orders.

## 1. Exact decomposition

At one common parameter on one honest torus define

\[
D=W_4-W_8=r_4-1\in\{-1,0,+1\}.
\]

By the winding-component classification there is a nonnegative integer `K` such that

\[
\boxed{
W_4=K+1_{\{D=+1\}},
\qquad
W_8=K+1_{\{D=-1\}}.}                                          \tag{1.1}
\]

Here `K>=1` in the neutral rank-one sector `D=0`, whereas `K=0` in the two charged endpoint sectors.

The symmetric/antisymmetric coordinates are

\[
S=W_4+W_8=2K+|D|,
\qquad
D=W_4-W_8.                                                     \tag{1.2}
\]

Thus `K` is the paired neutral component gas and `D` is a **bounded topological charge**.

## 2. Exact source decomposition

For real sources `s,t`,

\[
sW_4+tW_8=(s+t)K+s1_{\{D=+1\}}+t1_{\{D=-1\}}.                \tag{2.1}
\]

Hence configuration by configuration

\[
\left|sW_4+tW_8-(s+t)K\right|\le |s|+|t|.                    \tag{2.2}
\]

On a height-`m` torus, whenever the neutral-gas pressure exists,

\[
\psi_K(u)=\lim_{m\to\infty}\frac1m\log E e^{uK},
\]

(2.2) gives immediately

\[
\boxed{
\psi_{4,8}(s,t)
:=\lim_{m\to\infty}\frac1m\log E e^{sW_4+tW_8}
=\psi_K(s+t).}                                                 \tag{2.3}
\]

This recovers the pressure identity previously obtained from the bounded count difference, but gives it a typed interpretation: **the thermodynamic winding-count pressure has only one extensive source direction.**

The antisymmetric source couples only to a bounded charge and disappears after division by height.

## 3. Rank-one covariance at the extensive scale

Because `|D|<=1`,

\[
\operatorname{Var}(W_4-W_8)\le1.                              \tag{3.1}
\]

Therefore if the individual count variances are asymptotically linear in height,

\[
\frac1m\operatorname{Var}(W_4-W_8)\to0.                       \tag{3.2}
\]

Using complementary equality of the two marginal count pressures, their variance rates agree; call the common rate `v_w(p)` at fixed width.  Then

\[
\boxed{
\frac1m
\begin{pmatrix}
\operatorname{Var}W_4 & \operatorname{Cov}(W_4,W_8)\\
\operatorname{Cov}(W_4,W_8) & \operatorname{Var}W_8
\end{pmatrix}
\longrightarrow
v_w(p)
\begin{pmatrix}1&1\\1&1\end{pmatrix}.}                        \tag{3.3}
\]

In particular, when `v_w(p)>0`,

\[
\boxed{\operatorname{Corr}(W_4,W_8)\to1.}                     \tag{3.4}
\]

This is the same-parameter long-cylinder limit.  It deliberately contrasts with the independent Poisson limits in the two **separated parameter windows**.

Any simulation or transfer implementation of common-parameter component counts has a strong structural control: the antisymmetric variance is at most one at every finite size.

## 4. LDP consequence

Suppose `K_m/m` satisfies an LDP with speed `m` and rate `I_K`.  Since

\[
\left|W_4-K\right|\le1,
\qquad
\left|W_8-K\right|\le1,                                      \tag{4.1}
\]

all three normalized variables are exponentially equivalent:

\[
\frac{W_4}m,\quad\frac{W_8}m,\quad\frac Km.
\]

Therefore they have the **same** speed-`m` rate function.  More generally the joint normalized vector is supported asymptotically on the diagonal:

\[
(W_4/m,W_8/m)\approx(k,k).
\]

The topological charge has no extensive rate degree of freedom.  A bivariate extensive count model that assigns a nonzero rate to `W_4/m-W_8/m` violates exact topology.

## 5. Why the matching observable is invisible to bulk count pressure

The matching observable is

\[
M(p)=E D=P_2(p)-P_0(p).                                       \tag{5.1}
\]

It is an `O(1)` response of the bounded topological sector.  By contrast the count pressure in (2.3) is an `O(m)` bulk quantity.  Therefore no derivative of the **per-row bulk pressure** in the antisymmetric source can recover `M`:

\[
\lim_{m\to\infty}\frac1m E D=0                               \tag{5.2}
\]

regardless of the finite sign of `M`.

This gives a structural explanation for a phenomenon already seen in finite transfer calculations: slow/leading modes common to topological sectors can cancel from the matching difference.  Such cancellation is not, by itself, evidence for an exotic field or Jordan block.  The balance observable is designed to remove the neutral bulk contribution and retain the subextensive topological charge.

A transfer-matrix paper should therefore distinguish:

1. **bulk Perron pressure / neutral component gas**;
2. **sector amplitudes and subleading topological-charge response**;
3. **the finite zero-charge condition `M=0`.**

These are not interchangeable spectral quantities.

## 6. Integrated birth coordinates recover the limiting two-atom centres

The exact birth identities give

\[
E G=\int_0^1P_1(p)\,dp,
\qquad
E C=-\frac12\int_0^1M(p)\,dp,                                \tag{6.1}
\]

where

\[
G=T_2-T_1,\qquad C=\frac{T_1+T_2-1}{2}.
\]

Suppose along an exponential-aspect sequence the two individually sharp births converge to deterministic centres `(a,b)`.  Then

\[
E G\to b-a,
\qquad
E C\to\frac{a+b-1}{2}.                                       \tag{6.2}
\]

Hence the centres can be reconstructed from **integrated rank probabilities** without fitting quantiles:

\[
\boxed{
a=\frac{1-J-A}{2},
\qquad
b=\frac{1-J+A}{2},}                                           \tag{6.3}
\]

where

\[
A=\lim\int_0^1P_1(p)\,dp,
\qquad
J=\lim\int_0^1M(p)\,dp.                                      \tag{6.4}
\]

The strict matching mass gap proved in `matching-enhancement-mass-gap-20260914.md` gives `a+b>1`, and therefore the new sign prediction

\[
\boxed{
\int_0^1M_n(p)\,dp\longrightarrow1-a-b<0.}                    \tag{6.5}
\]

So the limiting complement-odd centre shift has an equivalent **area under the matching curve** interpretation.

This is useful for archives that store rank counts over occupation number: the two centre combinations can be estimated by exact beta/binomial integrations rather than repeated root/quantile solves.

## 7. Root balance is a defect-amplitude problem inside the plateau

In a geometry with a macroscopic rank-one plateau, the neutral sector has overwhelming probability over most of the plateau and both charged sectors are rare.  The matching root is not located by the extensive neutral gas: it is the point where the two rare defect weights are equal,

\[
P(D=+1)=P(D=-1).                                               \tag{7.1}
\]

The exact factorization from `rare-charge-balance-20260914.md`,

\[
F'(p_*)=\chi(p_*)H'(p_*),                                     \tag{7.2}
\]

makes the hierarchy quantitative.  The charged susceptibility `chi=1-P_1` may be tiny, while the conditional defect odds `H` vary extremely rapidly.  This is precisely how a well-defined finite balance root survives while the unconditioned limiting CDF becomes flat.

## 8. Claim boundary

Sections 1--5 are deterministic/asymptotic algebra once the fixed-width pressure/variance limits exist.  Section 6 additionally uses the already-proved sharpness of each birth around its finite centre and the strict centre theorem on this continuation branch.  No continuum-field interpretation is made.
