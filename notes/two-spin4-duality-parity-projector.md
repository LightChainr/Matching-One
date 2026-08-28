# Two spin-4 sectors and the matching/duality-parity projector

Status: synthesis of high-statistics server evidence with a CFT operator-classification precedent. The operator identification remains conjectural; the two scaling sectors are empirical hypotheses with falsifiable controls.

## 1. High-statistics evidence now supports two angular sectors

For the same-N Gaussian tori, write

\[
S=\frac{R_G+R_{\hat G}}2,\qquad
D=\frac{R_G-R_{\hat G}}2=\frac{M}{2}.
\]

For two orientations define the spin-4 projector

\[
P_4[X]=\frac{X(\theta_1)-X(\theta_2)}
{\cos4\theta_1-\cos4\theta_2}.
\]

The independent 100M-replica P31 run gives, for the `either` wrapping convention, the matching-even amplitudes scaled for an `N^-1=L^-2` law:

| N | `N P4[S]` | SE |
|---:|---:|---:|
| 65 | 0.0102593 | 0.0014433 |
| 85 | 0.0130056 | 0.0018008 |
| 130 | 0.0124497 | 0.0030798 |
| 145 | 0.00718869 | 0.00239126 |
| 170 | 0.00816274 | 0.00375246 |

A simple inverse-variance constant diagnostic gives

\[
\boxed{A_I\approx0.01060\pm0.00094}
\]

with chi-square about `4.66/4`. This is a strong nonzero matching-even spin-4 signal compatible with `L^-2` over the present sizes.

Independently, the matching-function sector satisfies

\[
\Delta M\approx A_M\,\Delta\cos4\theta\,N^{-13/8},
\]

with P31 pooled

\[
\boxed{A_M=0.7885\pm0.0352}
\]

and chi-square `1.53/4`. P32 then selected the fixed `13/8` H4 law out of sample over zero effect and found no predictive need for H8, a logarithm, or a free exponent at current precision.

Thus the data no longer suggest that only one angular correction is present. They are consistent with **two different spin-4 sectors with different matching parity and scaling dimensions**.

## 2. CFT precedent: square lattices admit two such families

Caselle, Hasenbusch, Pelissetto & Vicari (J. Phys. A 35, 4861; arXiv:cond-mat/0106372) classify irrelevant operators in the 2D Ising CFT. Their square-lattice list contains

- identity-family spin-4 `Q_4^I + Qbar_4^I`, with scaling dimension `x=4`;
- energy-family spin-4 `Q_4^epsilon + Qbar_4^epsilon`, whose dimension is `x_epsilon+4`.

They also emphasize that square-lattice scalar quantities admit continuum spins `4j`, and that identity- and energy-family operators have opposite duality parity in the Ising model.

For percolation,

\[
x_t=5/4,
\]

so the direct energy-family analogue has

\[
\boxed{x_t+4=21/4}.\]

The exact `c=0,h=5/8` Virasoro check in this repository confirms that the level-2 null relation does not kill the level-4 quasiprimary: a non-null, non-total-derivative direction exists with bulk spin `+/-4` and `x=21/4`.

## 3. Working duality-parity identification

The strongest current structural hypothesis is:

### Sector I: identity-family / matching-even

\[
S_4(L)\sim A_I L^{-2}.
\]

Candidate continuum field: the square-lattice identity-family spin-4 anisotropy, schematically `T^2+Tbar^2` (with the usual `c=0` caveats/logarithmic partners).

### Sector epsilon: thermal-family / matching-odd

\[
D_4(L)\sim A_\epsilon L^{-13/4}.
\]

Candidate continuum field: the level-4 spin-4 quasiprimary in the thermal `h=hbar=5/8` family.

The site-matching/complement transformation is conjectured to act in the continuum as a duality-like involution on these irrelevant couplings:

- identity-family spin-4 is even and survives the sum `S`;
- thermal-family spin-4 is odd and survives the difference `D`.

This is an identification hypothesis, not a theorem. Matching parity must ultimately be derived from the lattice-to-CFT map or validated by independent self-dual controls.

## 4. Why the faster-decaying D sector is currently numerically larger

Because `D=M/2`, its present scaled amplitude is approximately

\[
A_D=A_M/2\approx0.3943.
\]

Compared with `A_I≈0.01060`, the raw ratio is predicted to behave as

\[
\frac{D_4}{S_4}
\approx
\frac{A_D}{A_I}N^{-5/8}
\approx37.2\,N^{-5/8}.
\]

Therefore `D` can be larger at the current small/moderate systems even though it decays faster. The two contributions cross near

\[
N\sim3.3\times10^2,\qquad L\sim18.
\]

This resolves the apparent contradiction between the early observation that the matching-odd signal was larger and the asymptotic expectation that the `x=4` identity-family correction decays more slowly.

## 5. Root shift follows only from the matching-odd sector

P35 verifies configuration-independent amplitude closure for the finite curves:

\[
-\Delta p^*\frac{\overline{M'}}{\Delta M}\approx1
\]

to `3e-4` or better at all five threshold-rank sizes, while

\[
B=N^{-3/8}\overline{M'}
\]

is nearly constant (`1.7514 -> 1.7462`).

Hence the angular root bias implied by the thermal-family sector is

\[
\Delta p^*
\sim
-\frac{A_M}{B}\,\Delta\cos4\theta\,N^{-2}
\]

or equivalently

\[
\boxed{
-\frac{N^2\Delta p^*}{\Delta\cos4\theta}
\to\frac{A_M}{B}.
}
\]

The angular factor is essential. A test of bare `N^2 Delta p` across different orientation pairs is not the correct H4 radial invariant.

Using the current high-stat `A_M≈0.7885` and `B≈1.748` gives the rough no-new-fit prediction

\[
A_p\equiv-\frac{N^2\Delta p^*}{\Delta\cos4\theta}
\approx0.451.
\]

The current 10M threshold-rank root estimates are still noisy enough that this is not a sharp test; a higher-stat threshold-rank production should use `A_p`, not bare `N^2 Delta p`, as the primary root-amplitude closure statistic.

## 6. Strong independent controls

### Square-bond self-dual control

On a genuinely self-dual square-lattice critical model, a duality-odd thermal-family spin-4 coupling should be forced to vanish at the self-dual point (or obey a correspondingly stronger cancellation), while the duality-even identity-family anisotropy may remain.

Therefore run the same orientation/thermal-parity decomposition for square-bond percolation at exact `p_c=1/2`.

Working prediction:

- an `x=4`, spin-4 even anisotropy can be present;
- the central duality-odd `x=21/4` analogue should be absent/suppressed.

This is a much stronger control of the parity identification than merely checking that square bond has orientation-dependent finite-size corrections.

### Self-matching triangular-site control

Triangular site percolation is self-matching at `p=1/2`, but its `C6` microscopic symmetry does not admit a generic spin-4 lattice scalar. It is useful as a negative symmetry/parity control, not as the closest positive C4 control.

## 7. Next numerical tests

1. Use P33 threshold-rank curves to report the thermal-even/odd parts of `P4[S]` and `P4[D]` explicitly.
2. Increase threshold-rank statistics until the angular-normalized root amplitude `A_p` can distinguish constant/slow correction models.
3. Run the square-bond self-dual C4 parity control.
4. Execute issue #38 with a **new post-registration seed** for the exact Gaussian-doubling ratio.
5. On held-out sizes, compare the two-family fixed model

\[
S_4=A_I N^{-1},\qquad M_4=A_M N^{-13/8}
\]

against free exponents and log alternatives without using root data for model selection.

## Falsification

The projector interpretation should be rejected or revised if any of the following persist with higher-stat held-out data:

- `S_4` does not approach an `N^-1` law;
- `M_4` fails the frozen H4 `N^-13/8` prediction;
- the self-dual square-bond control carries a comparable central duality-odd thermal-family signal;
- the Gaussian doubling sign/ratio test fails prospectively;
- the thermal-even component does not contain the central root-moving signal.
