# Exact charge-neutral coordinates for the common black/white crossover

2026-09-14.  Exact finite reparameterization of the same-parameter joint law.  This is a direct answer to the “minimal topology-compatible family” part of #767 before any near-critical scaling hypothesis is introduced.

## 1. Exact support

At one common parameter, with black NN and complementary white matching on the same labels,

\[
(W_4,W_8)\in\{(0,1),(1,0),(K,K):K\ge1\}.                     \tag{1.1}
\]

Equivalently define the bounded topological charge

\[
D=W_4-W_8=r_4-1\in\{-1,0,1\}.                                \tag{1.2}
\]

The charged sectors `D=-1,+1` have `K=0`; the neutral sector `D=0` has `K>=1`.

## 2. Charged susceptibility and charge fugacity

Write

\[
\boxed{\chi(p)=P_0(p)+P_2(p)=1-P_1(p)}                        \tag{2.1}
\]

for the total charged-sector probability, and

\[
\boxed{\theta(p)=\log\frac{P_2(p)}{P_0(p)}}                   \tag{2.2}
\]

for the log charge odds.

For every interior `0<p<1`, both endpoint probabilities are positive, so `theta` is finite.  Since `P_2` is strictly increasing and `P_0` strictly decreasing,

\[
\boxed{\theta(p)\text{ is strictly increasing}.}              \tag{2.3}
\]

Solving (2.1)--(2.2),

\[
\boxed{
P_2
=\chi\frac{e^{\theta/2}}{2\cosh(\theta/2)},
\qquad
P_0
=\chi\frac{e^{-\theta/2}}{2\cosh(\theta/2)}.}               \tag{2.4}
\]

Thus the matching observable is exactly

\[
\boxed{
M=P_2-P_0
=\chi\tanh(\theta/2).}                                       \tag{2.5}
\]

This identity holds at every finite size and every interior parameter.

## 3. The matching root is the zero-fugacity point

The unique matching root satisfies

\[
P_2=P_0.                                                       \tag{3.1}
\]

Therefore

\[
\boxed{p_*=\theta^{-1}(0).}                                   \tag{3.2}
\]

The location of the root is controlled by the **sign change of the charge fugacity**, while `chi` controls how much total charged mass is available to reveal that sign.

This is the exact algebra behind “balance without concentration”:

- in a macroscopic rank-one plateau, `chi` may be extremely small;
- `theta` can still run monotonically from large negative to large positive values;
- the finite root remains the unique point `theta=0` even though the unconditioned birth CDF is nearly flat at `1/2`.

## 4. Exact slope factorization at the root

Differentiate (2.5):

\[
M'
=\chi'\tanh(\theta/2)
+\frac{\chi\theta'}2\operatorname{sech}^2(\theta/2).          \tag{4.1}
\]

At the root `theta=0`,

\[
\boxed{M'(p_*)=\frac12\chi(p_*)\theta'(p_*).}                 \tag{4.2}
\]

Since

\[
F=(1+M)/2,                                                     \tag{4.3}
\]

we obtain

\[
\boxed{F'(p_*)=\frac14\chi(p_*)\theta'(p_*).}                 \tag{4.4}
\]

This is the log-odds version of the earlier rare-charge factorization.  If

\[
H=P_2/(P_0+P_2)=\frac1{1+e^{-\theta}},                        \tag{4.5}
\]

then `H'(p_*)=theta'(p_*)/4`, recovering `F'=chi H'`.

The nonlinear expansion is

\[
M=\chi\left(\frac\theta2-\frac{\theta^3}{24}+O(\theta^5)\right).\tag{4.6}
\]

So a common-window theory can separately model the charged mass `chi` and the charge field `theta`; conflating them loses the exact root mechanism.

## 5. Minimal joint PGF after forgetting slope

Let

\[
H_p(z)=E[z^K\mid D=0]                                        \tag{5.1}
\]

be the neutral rank-one component-count PGF.  It is supported on positive integers and satisfies `H_p(1)=1`.

The exact same-parameter joint PGF is

\[
\boxed{
G_p(s,t)
=\chi
 \frac{e^{\theta/2}s+e^{-\theta/2}t}{2\cosh(\theta/2)}
 +(1-\chi)H_p(st).}                                           \tag{5.2}
\]

Thus the entire two-colour law is characterized by only three objects:

1. one scalar charged susceptibility `chi(p)`;
2. one scalar strictly monotone charge fugacity `theta(p)`;
3. one **one-dimensional** positive-integer neutral-count law `H_p`.

No bivariate count copula is needed or allowed by the exact support.

At the root,

\[
\boxed{
G_{p_*}(s,t)
=\frac{\chi_*}{2}(s+t)
 +(1-\chi_*)H_{p_*}(st).}                                    \tag{5.3}
\]

So a root-window limit is determined by the single charged mass `chi_*` plus the neutral count law.

## 6. Slope-marked minimal law

Let `L` be the unique projective slope in the neutral sector and let

\[
H_{p,\ell}(z)
=E[z^K\mid D=0,L=\ell],                                      \tag{6.1}
\]

with slope weights

\[
w_{p,\ell}=P(L=\ell\mid D=0).                               \tag{6.2}
\]

Then

\[
\boxed{
G_p(s,t;\{y_\ell\})
=\chi
 \frac{e^{\theta/2}s+e^{-\theta/2}t}{2\cosh(\theta/2)}
 +(1-\chi)
 \sum_\ell w_{p,\ell}y_\ell H_{p,\ell}(st).}                 \tag{6.3}
\]

This is the fully marked topology-compatible family.  Different slopes cannot be independently active in the neutral sector; that global hard-core constraint is already encoded by the single sum over `ell`.

## 7. Exact moment identities useful for data reduction

Let

\[
\mu=EK,
\qquad
\sigma_K^2=\operatorname{Var}K,                              \tag{7.1}
\]

where `K=0` in charged sectors as in the exact decomposition.  Because

\[
KD=0                                                          \tag{7.2}
\]

pathwise, several joint count moments simplify.

First,

\[
EW_4=\mu+P_2,
\qquad
EW_8=\mu+P_0.                                                 \tag{7.3}
\]

Also

\[
W_4W_8=K^2                                                     \tag{7.4}
\]

pathwise, so

\[
\boxed{
\operatorname{Cov}(W_4,W_8)
=\sigma_K^2-\mu\chi-P_0P_2.}                                 \tag{7.5}
\]

Moreover

\[
\boxed{
E(W_4-W_8)^2=\chi,}                                           \tag{7.6}
\]

and hence

\[
\boxed{
\operatorname{Var}(W_4-W_8)
=\chi-M^2.}                                                    \tag{7.7}
\]

At the matching root `M=0`, the antisymmetric variance is exactly `chi_*`.

These are cheap exact controls for any same-parameter joint-count implementation.

## 8. Boundary conditions from the separated Poisson windows

In the lower separated birth window, the existing component theorem gives

\[
K_4\Rightarrow\operatorname{Poi}(\lambda).                   \tag{8.1}
\]

Topology converts this into

\[
D=-1\quad\text{when }K=0,
\qquad
D=0\quad\text{when }K\ge1,                                  \tag{8.2}
\]

with `D=+1` negligible at that scale.  Therefore

\[
\chi\to e^{-\lambda},
\qquad
\theta\to-\infty,                                             \tag{8.3}
\]

while the neutral conditional count is zero-truncated Poisson.

The upper separated window is the reflected limit

\[
\chi\to e^{-\lambda},
\qquad
\theta\to+\infty.                                             \tag{8.4}
\]

Thus a genuine common-window crossover must interpolate the charge fugacity from `-infinity` to `+infinity` while matching these one-sided neutral-count boundary laws.

This is much more restrictive than choosing an arbitrary correlated bivariate Poisson family.

## 9. Three geometric regimes in these coordinates

The coordinates separate the regimes requested in #767.

### Fixed positive exponential aspect `d>0`

The two birth windows stay separated.  In the plateau between them,

\[
\chi\to0,                                                     \tag{9.1}
\]

while `theta` still passes through zero at the finite matching root.  The common-parameter law is overwhelmingly neutral.

### Fixed macroscopic aspect ratio

There is no exponentially long rank-one plateau forced by geometry.  `chi` may remain `O(1)` through the critical crossover, and the neutral count law need not be a rare-component Poisson law.

### Intermediate `rho->infinity`, `log rho=o(w)`

This is the natural merging regime for a nontrivial crossover: `chi`, `theta`, and `H_p` may all have nondegenerate scaling limits.  Any proposed scaling function should be stated for these three objects (or the slope-marked extension), not for two unconstrained colour intensities.

## 10. A concrete crossover target

A minimal conditional programme is therefore:

\[
\chi_w(p)\to\Chi(x,\rho),
\qquad
\theta_w(p)\to\Theta(x,\rho),
\qquad
H_{w,p}(z)\to\mathcal H_{x,\rho}(z),                           \tag{10.1}
\]

under a declared near-critical scaling variable `x` and aspect parameter `rho`.

The exact constraints are

\[
0\le\Chi\le1,
\qquad
\Theta\text{ monotone in the thermal parameter},              \tag{10.2}
\]

`H` is a PGF on positive integers, and the separated-window limits are (8.3)--(8.4).

The matching root in the scaling theory is simply

\[
\boxed{\Theta=0.}                                             \tag{10.3}
\]

This is the smallest state space in which a common-window crossover can live without violating finite topology.

## 11. Claim boundary

Sections 1--9 are exact finite identities plus already-proved separated-window boundary conditions.  Section 10 is a conditional scaling programme.  No square-site near-critical universality function is asserted, and no fixed-p prefactor is analytically continued to criticality here.