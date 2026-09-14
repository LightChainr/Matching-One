# One-dimensional collapse of same-parameter two-colour count fluctuations

2026-09-14.  Exact consequence of the bounded topological charge

\[
D=W_4-W_8\in\{-1,0,1\}.
\]

This sharpens `neutral-gas-topological-charge-20260914.md`: not only the pressure and variance rates, but **every growing-scale fluctuation limit and every extensive LDP collapse to the diagonal**.

## 1. Uniform exponential-moment comparison

Write

\[
W_4=K+1_{\{D=1\}},\qquad
W_8=K+1_{\{D=-1\}}.                                           \tag{1.1}
\]

For real sources `s,t`, the source correction relative to the neutral gas is one of `0,s,t`.  Hence with

\[
c(s,t)=\max(|s|,|t|),
\]

we have pathwise

\[
\left|sW_4+tW_8-(s+t)K\right|\le c(s,t).                      \tag{1.2}
\]

Therefore at every finite size

\[
\boxed{
e^{-c(s,t)}E e^{(s+t)K}
\le E e^{sW_4+tW_8}
\le e^{c(s,t)}E e^{(s+t)K}.}                                  \tag{1.3}
\]

The error is multiplicative by a size-independent constant.  No asymptotic estimate is needed.

If the system has longitudinal size `m`, division by `m` gives

\[
\left|
\frac1m\log Ee^{sW_4+tW_8}
-
\frac1m\log Ee^{(s+t)K}
\right|
\le \frac{c(s,t)}m.                                           \tag{1.4}
\]

Thus any limiting scaled cumulant generating function satisfies

\[
\boxed{\psi_{4,8}(s,t)=\psi_K(s+t).}                           \tag{1.5}
\]

The antisymmetric source direction is exactly subextensive.

## 2. Joint LDP is infinite off the diagonal

For every realization,

\[
|W_4-K|\le1,\qquad |W_8-K|\le1.                               \tag{2.1}
\]

Hence

\[
\left\|
(W_4/m,W_8/m)-(K/m,K/m)
\right\|_\infty\le1/m.                                       \tag{2.2}
\]

Suppose `K_m/m` satisfies an LDP with speed `m` and good rate `I_K`.  Deterministic exponential equivalence gives the full joint rate function

\[
\boxed{
I_{4,8}(x,y)=
\begin{cases}
I_K(x),&x=y,\\
+\infty,&x\ne y.
\end{cases}}                                                   \tag{2.3}
\]

So a same-parameter bivariate count theory with a finite extensive cost for `x-y != 0` is structurally impossible.

This is stronger than saying the correlation tends to one: there is no second extensive large-deviation coordinate at all.

## 3. Every diverging-scale fluctuation limit is diagonal

Let `a_m->infinity` be any deterministic scale.  Then

\[
\frac{W_4-EW_4}{a_m}
-
\frac{K-EK}{a_m}
\to0,                                                         \tag{3.1}
\]

and the same for `W_8`, because each difference is bounded by a constant divided by `a_m`.

Therefore if

\[
\frac{K_m-EK_m}{a_m}\Rightarrow Z                             \tag{3.2}
\]

for **any** nondegenerate limit `Z`, not necessarily Gaussian, then

\[
\boxed{
\left(
\frac{W_4-EW_4}{a_m},
\frac{W_8-EW_8}{a_m}
\right)
\Rightarrow (Z,Z).}                                           \tag{3.3}
\]

Consequences include:

- a CLT for one count automatically gives the same CLT jointly for both;
- a stable/non-Gaussian scaling limit also transfers diagonally;
- moderate deviations with a growing normalization inherit the same one-dimensional collapse.

No separate two-colour fluctuation theorem is required once the neutral count limit is known.

## 4. Covariance matrix and principal directions

If `Var(K_m)->infinity`, then

\[
\operatorname{Var}(W_4)=\operatorname{Var}(K)+o(\operatorname{Var}K),
\]

\[
\operatorname{Var}(W_8)=\operatorname{Var}(K)+o(\operatorname{Var}K),
\]

and

\[
\operatorname{Cov}(W_4,W_8)=\operatorname{Var}(K)+o(\operatorname{Var}K).\tag{4.1}
\]

Thus the normalized covariance matrix tends to the rank-one projector on `(1,1)` and

\[
\boxed{\operatorname{Corr}(W_4,W_8)\to1.}                     \tag{4.2}
\]

The antisymmetric principal component is exactly `D`, whose variance is at most one at every finite size.

This yields a severe implementation check: if a long-cylinder common-parameter simulation reports an antisymmetric variance growing with length, the count dictionary is wrong.

## 5. Separation of bulk and root physics

The neutral gas `K` carries all extensive pressure, variance and LDP information.  The matching observable

\[
M=E(W_4-W_8)=ED                                                   \tag{5.1}
\]

lives entirely in the bounded topological defect sector.

Therefore the matching root is a **subextensive sector-amplitude problem sitting on top of an extensive one-dimensional gas**.  This explains why common transfer/Perron modes can cancel in `M` even when they dominate each individual count.

The correct spectral hierarchy is:

1. one bulk `(1,1)` count mode;
2. bounded endpoint/rank charge sectors;
3. the zero-charge condition selecting the matching root.

A second extensive count mode should not be introduced to explain root motion.

## 6. Generating-function fingerprint

At finite size define

\[
\Phi_m(s,t)=Ee^{sW_4+tW_8}.                                   \tag{6.1}
\]

Equation (1.3) means that along two source pairs with the same sum `s+t`, the ratio of generating functions is bounded uniformly in `m`:

\[
\boxed{
\left|\log\frac{\Phi_m(s,t)}{\Phi_m(s',t')}\right|
\le c(s,t)+c(s',t')
\quad\text{if }s+t=s'+t'.}                                   \tag{6.2}
\]

Thus any transfer calculation in which this log-ratio grows linearly with cylinder height contradicts exact topology before any asymptotic interpretation is attempted.

## 7. Claim boundary

All statements are deterministic/probabilistic consequences of the exact finite decomposition (1.1), conditional only on existence of the one-dimensional limit invoked in each asymptotic corollary.  No Poisson assumption is used.