# A common primal/dual safe-sector pivotal amplitude

2026-09-14.  Strong scaling conjecture extracted from the exact Perron/Russo identities.

The fixed-width charge slope is already exactly typed as a sum of two conditioned pivotal intensities.  The new observation is that the NN-safe and complementary-matching-safe contributions appear to approach the SAME leading Bernoulli-score amplitude.

## 1. Exact finite-width quantities

At the charge-coexistence root `p_w`, put

\[
q_w=1-p_w.                                                     \tag{1.1}
\]

Let

\[
\bar K^0_{4,w},\qquad \bar K^0_{8,w}                         \tag{1.2}
\]

be the mean occupied sites in one row under the two Perron/Doob safe phases.  Define their occupation deficits relative to an unconditioned Bernoulli row:

\[
\boxed{
d_{4,w}=wp_w-\bar K^0_{4,w},
\qquad
d_{8,w}=wq_w-\bar K^0_{8,w}.}                             \tag{1.3}
\]

The exact score/pivotal identities give

\[
d_{4,w}=p_w r^0_{4,w},
\qquad
d_{8,w}=q_w r^0_{8,w},                                 \tag{1.4}
\]

where `r^0` is the safe-conditioned topological pivotal density per row.

Also

\[
\boxed{
p_wq_w\Theta'_w(p_w)=d_{4,w}+d_{8,w}.}                 \tag{1.5}
\]

No scaling assumption enters (1.3)--(1.5).

## 2. Transfer data

For `w=4,...,8`:

| `w` | `d4 w^(1/4)` | `d8 w^(1/4)` | `d4/d8` |
|---:|---:|---:|---:|
| 4 | 0.4152842 | 0.4268705 | 0.972858 |
| 5 | 0.4130767 | 0.4195635 | 0.984539 |
| 6 | 0.4116079 | 0.4156285 | 0.990327 |
| 7 | 0.4105738 | 0.4132612 | 0.993497 |
| 8 | 0.4098310 | 0.4117281 | 0.995392 |

The approach to a common amplitude is substantially cleaner than either raw derivative alone.

The corresponding conditional pivotal intensities satisfy

| `w` | `r4 w^(1/4)` | `r8 w^(1/4)` | `r4/r8` |
|---:|---:|---:|---:|
| 4 | 0.7021850 | 1.0447587 | 0.672103 |
| 5 | 0.6974869 | 1.0289366 | 0.677872 |
| 6 | 0.6946883 | 1.0199656 | 0.681090 |
| 7 | 0.6928117 | 1.0144357 | 0.682953 |
| 8 | 0.6914963 | 1.0108043 | 0.684105 |

while

\[
\frac{1-p_c}{p_c}
=0.6870631169\ldots.                                         \tag{2.1}
\]

## 3. Common-amplitude conjecture

The natural scaling statement is

\[
\boxed{
d_{4,w}w^{1/4}\to C_{\rm piv},
\qquad
d_{8,w}w^{1/4}\to C_{\rm piv}.}                        \tag{3.1}
\]

Equivalently,

\[
\boxed{
\frac{r^0_{4,w}}{r^0_{8,w}}
\longrightarrow\frac{q_c}{p_c}.}                             \tag{3.2}
\]

The finite data suggest `C_piv` around `0.40--0.41`; no precision claim is made from these widths.

Then (1.5) forces one thermal-slope amplitude:

\[
\boxed{
\Theta'_w(p_w)w^{1/4}
\longrightarrow
\frac{2C_{\rm piv}}{p_cq_c}.}                                \tag{3.3}
\]

For `C_piv≈0.41`, the right side is about `3.40`, matching the direct transfer sequence.

## 4. Continuum interpretation

The two safe sectors converge to the same magnetic/topological primary.  Their derivatives with respect to complementary Bernoulli thermal coordinates probe opposite sides of the same thermal perturbation.

The occupation deficit `wp-Kbar` is the Bernoulli score response of that sector.  Unlike the raw pivotal count, it already includes the microscopic parameter metric factor.  Therefore a common leading amplitude in (3.1) is exactly what one expects if the complement map `q=1-p` fixes the relative normalization of the two thermal coordinates in the scaling limit.

This is stronger than merely saying both sides have exponent `1/4`: it asserts a primal/dual amplitude relation.

## 5. Four-arm reading

The exact identity

\[
d_{G,w}=p\,r^0_{G,w}                                        \tag{5.1}
\]

and the critical four-arm mechanism suggest

\[
r^0_{G,w}\asymp C_G w^{-1/4}.                               \tag{5.2}
\]

Equation (3.1) is then the amplitude relation

\[
p_c C_4=q_c C_8=C_{\rm piv}.                                \tag{5.3}
\]

So the two raw pivotal amplitudes are NOT predicted equal; their Bernoulli-weighted amplitudes are.

This is an important distinction for any direct pivotal sampling experiment.

## 6. A stronger diagnostic than numerical differentiation

Future wider transfer calculations should report all three objects:

\[
d_4w^{1/4},\qquad d_8w^{1/4},\qquad
p q\Theta'_w w^{1/4}.                                       \tag{6.1}
\]

The exact identity demands the third equal the sum of the first two at every width.  The scaling conjecture further demands the first two coalesce.

A direct Q-process pivotal sampler can test the same statement without taking eigenvalue finite differences.

## 7. Relation to the sector-odd correction

This common thermal amplitude is a **sector-even normalization statement about the derivative**.  It does not imply the critical sector energies are equal at finite width.  Their tiny difference is the separate sector-odd irrelevant correction responsible for the `w^-4` pseudo-critical shift.

Thus the picture is:

- thermal response: common primal/dual amplitude, order `w^-1/4`;
- critical mismatch: much smaller rotational sector-odd amplitude, order `w^-17/4` on the square lattice.

## 8. Claim boundary

Equations (1.3)--(1.5) are exact.  The convergence to a common amplitude and the ratio (3.2) are scaling conjectures strongly supported by the transparent widths.  A rigorous proof would require a near-critical primal/matching scaling limit with the relative thermal metric fixed by complement duality, or an equivalent arm-amplitude theorem.
