# A dual-odd universal thermal scaling function for the charge free energy

2026-09-14.  Exact cumulant identities plus a finite-width scaling conjecture.

The first-derivative charge slope was already identified with safe-sector pivotal density.  Going to higher logit derivatives reveals a stronger structure: the leading NN/matching thermal scaling function appears to be ODD under the complementary-sector exchange, so even derivatives cancel from the charge difference while odd derivatives add.

## 1. Use the Bernoulli logit coordinate

Put

\[
h=\log\frac p{1-p},\qquad p=\frac{e^h}{1+e^h}.              \tag{1.1}
\]

Complementation sends

\[
p\leftrightarrow1-p
\quad\Longleftrightarrow\quad
h\leftrightarrow-h.                                          \tag{1.2}
\]

For a safe transfer kernel, write the unnormalized row matrix as

\[
A_G(h)=\sum_{k=0}^w e^{hk}A_{G,k},                            \tag{1.3}
\]

and let `Lambda_G(h)` be its Perron root.  The Bernoulli-normalized safe Perron root is

\[
\lambda_G(h)=\frac{\Lambda_G(h)}{(1+e^h)^w},                  \tag{1.4}
\]

so

\[
I_G(h)=-\log\lambda_G(h)
=w\log(1+e^h)-\log\Lambda_G(h).                              \tag{1.5}
\]

The charge free-energy difference is

\[
\boxed{
\Theta_w(h)=I_{4,w}(h)-I_{8,w}(-h).}                         \tag{1.6}
\]

Its zero is the semi-infinite charge root.

## 2. Exact Q-process cumulant identities

The Perron transform of `A_G(h)` defines a stationary finite-state Q-process.  Let `K_j` be the number of occupied sites added in row `j` under that process.

Standard Perron pressure differentiation gives

\[
\partial_h\log\Lambda_G
=\bar K_G,                                                     \tag{2.1}
\]

and

\[
\partial_h^2\log\Lambda_G
=\sigma_G^2,                                                   \tag{2.2}
\]

where

\[
\sigma_G^2
=\operatorname{Var}(K_0)
+2\sum_{j\ge1}\operatorname{Cov}(K_0,K_j)                    \tag{2.3}
\]

is the Green--Kubo asymptotic variance per row.  More generally higher derivatives are asymptotic cumulant rates of the additive row occupancy.

Therefore

\[
\boxed{
I_{G,h}=wp-\bar K_G,}                                         \tag{2.4}
\]

\[
\boxed{
I_{G,hh}=wp(1-p)-\sigma_G^2,}                                 \tag{2.5}
\]

and

\[
I_{G,hhh}
=wp(1-p)(1-2p)-\kappa_{3,G}^{\rm asym}.                       \tag{2.6}
\]

Equation (2.4) is the already-used safe occupation deficit / pivotal identity.  Equations (2.5)--(2.6) show that higher charge derivatives measure **deficits of thermal cumulants** under topological survival conditioning.

## 3. Exact parity of charge derivatives

Differentiate (1.6):

\[
\Theta_h
=I_{4,h}(h)+I_{8,h}(-h),                              \tag{3.1}
\]

\[
\Theta_{hh}
=I_{4,hh}(h)-I_{8,hh}(-h),                                   \tag{3.2}
\]

\[
\Theta_{hhh}
=I_{4,hhh}(h)+I_{8,hhh}(-h).                                 \tag{3.3}
\]

Thus odd thermal cumulant responses ADD between the complementary sectors, while even responses SUBTRACT.

This exact algebra is the natural place to look for a continuum dual-odd scaling function.

## 4. Finite transfer evidence

At the charge root `h_w`, the transparent transfer gives:

| `w` | `w^(1/4) Theta_h` | `w^(1/2) Theta_hh` | `w^(-5/4) Theta_hhh` |
|---:|---:|---:|---:|
| 4 | 0.84215 | -0.08115 | 0.06078 |
| 5 | 0.83264 | -0.08124 | 0.06223 |
| 6 | 0.82724 | -0.08291 | 0.06306 |
| 7 | 0.82384 | -0.08499 | 0.06355 |
| 8 | 0.82156 | -0.08714 | 0.06385 |

The half-step finite-difference controls in `fixed-width-charge-spectrum-derivatives-w4-w8-20260914.json` agree at the displayed precision.

The individual second derivatives behave very differently from their difference:

\[
I_{4,hh}/\sqrt w
=0.1766,0.1765,0.1765,0.1766,0.1767,                        \tag{4.1}
\]

while

\[
I_{8,hh}/\sqrt w
=0.1969,0.1927,0.1903,0.1887,0.1876.                        \tag{4.2}
\]

Each sector separately has the expected thermal second-derivative scale `sqrt(w)`, but their difference is only about `w^-1/2`:

\[
\Theta_{hh}\sqrt w\approx-0.08\text{ to }-0.09.             \tag{4.3}
\]

So the leading `w^(+1/2)` even thermal response is cancelling very strongly between primal and complementary matching sectors.

## 5. Dual-odd scaling-function conjecture

The simplest continuum organization is

\[
\boxed{
\Theta_w(h)
=\frac1w\,\mathcal F(X)
+\text{subleading sector-odd lattice corrections},}           \tag{5.1}
\]

with

\[
X=c_h(h-h_c)w^{3/4},                                         \tag{5.2}
\]

and

\[
\boxed{\mathcal F(-X)=-\mathcal F(X).}                       \tag{5.3}
\]

Here `c_h` is the thermal metric factor fixed by the matching/complement convention.  At the level of this branch, the especially strong amplitude matching in the safe pivotal deficits suggests that the Bernoulli logit is already close to the natural common normalization.

If (5.3) holds, then

\[
\Theta_h\asymp w^{-1/4},                                     \tag{5.4}
\]

\[
\Theta_{hh}=0\times w^{1/2}+\text{subleading},               \tag{5.5}
\]

\[
\Theta_{hhh}\asymp w^{5/4}.                                  \tag{5.6}
\]

These are exactly the three patterns visible in section 4.

The semi-infinite root displacement `p_w-p_c=O(w^-4)` is far smaller than the thermal coordinate scale `w^-3/4`, so evaluating derivatives at `p_w` rather than exactly `p_c` does not affect the leading scaling function.

## 6. First nonlinear coefficients

Use the finite-width scaled derivatives as rough controls for

\[
\mathcal F(X)=a_1X+\frac{a_3}{6}X^3+O(X^5).                  \tag{6.1}
\]

In the uncalibrated logit metric, the data suggest

\[
a_1\approx0.8\text{--}0.82,
\qquad
a_3\approx0.064.                                             \tag{6.2}
\]

The ratio

\[
\frac{a_3}{6a_1}\approx1.3\times10^{-2}                      \tag{6.3}
\]

is small but clearly nonzero in the current widths.

These numbers should not be called universal until the thermal metric factor is fixed.  The **parity** prediction (5.3), and dimensionless ratios after a declared normalization, are the robust targets.

## 7. Fixed-aspect charge crossover

For a torus with aspect ratio

\[
\rho=m/w,                                                     \tag{7.1}
\]

the periodic charge log-odds at leading scaling order should be

\[
\boxed{
\theta_{w,m}(h)
\approx \rho\,\mathcal F(X)}                                 \tag{7.2}
\]

before exponentially small longitudinal trace corrections are added.

This gives #767 a much sharper common-window object than two unconstrained colour intensities:

- `rho` controls how strongly the dual-odd scaling function is amplified;
- `X` is the ordinary thermal coordinate;
- the matching root is `X=0` at leading universal order;
- sector-odd lattice anisotropy shifts that zero only at the much smaller `w^-4` scale on the square lattice.

Thus the charge-neutral crossover coordinates and the fixed-width transfer now fit into one two-variable scaling picture.

## 8. A hierarchy of derivative tests

The conjecture predicts alternating behavior:

\[
\Theta^{(2j+1)}_h\asymp
w^{-1+(2j+1)3/4},                                             \tag{8.1}
\]

while

\[
\Theta^{(2j)}_h
\]

has its leading universal thermal contribution cancelled and is controlled by irrelevant/metric corrections.

So the next discriminating transfer test is not another root fit but the fourth and fifth logit derivatives, with careful numerical differentiation or exact Perron perturbation.  A fifth derivative at scale `w^(11/4)` and a strongly suppressed fourth derivative would be powerful evidence for (5.3).

## 9. Claim boundary

The Perron/cumulant identities and parity signs (3.1)--(3.3) are exact.  The finite derivative table is a reproducible small-width control.  The odd universal scaling function, the cancellation of all leading even derivatives, and the interpretation of `rho F(X)` as the full fixed-aspect charge scaling law are conjectural scaling statements requiring either a near-critical primal/matching scaling-limit proof or substantially wider transfer evidence.
