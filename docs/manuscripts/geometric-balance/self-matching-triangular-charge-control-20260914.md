# Triangular-site self-matching as an exact zero control for sector-odd charge corrections

2026-09-14.  Exact topological control; no scaling theory is needed for the principal statement.

## 1. Self-matching triangulations

For independent site percolation on a planar triangulation, every face already has all pairs of its vertices adjacent.  The site matching graph therefore adds no new face diagonal.  In particular the standard triangular lattice is self-matching:

\[
G^\star_{match}=G_{tri}.                                      \tag{1.1}
\]

On an honest periodic triangular triangulation, digital Alexander duality gives configurationwise

\[
\boxed{
r_G(\omega)+r_G(\omega^c)=2.}                               \tag{1.2}
\]

At `p=1/2`, `omega` and `omega^c` have the same law.  Hence

\[
\boxed{P_0(1/2)=P_2(1/2)}                                    \tag{1.3}
\]

for every finite honest torus.

Thus the matching/root observable has the exact finite value

\[
\boxed{p^*_{\Lambda}=1/2}                                    \tag{1.4}
\]

with no finite-size drift.

## 2. Fixed-width charge free energy

Let `I^0_{tri,w}(p)` be the no-horizontal-homology strip free energy.  Because the primal and matching graphs are literally the same graph,

\[
\Theta_w^{tri}(p)
=I^0_{tri,w}(p)-I^0_{tri,w}(1-p).                             \tag{2.1}
\]

Therefore

\[
\boxed{\Theta_w^{tri}(1/2)=0}                                \tag{2.2}
\]

for every width `w`, and strict monotonicity gives the unique semi-infinite charge root

\[
\boxed{p_w^{ch}=1/2\quad\text{for all }w.}                   \tag{2.3}
\]

This is the probability/digital-Alexander version of the size-independent factor that appears in exactly solvable/self-dual cases of the graph-polynomial/eigenvalue method.

## 3. Consequence for sector-odd irrelevant channels

At the exact self-matching point, the full charge free-energy difference vanishes:

\[
I^0_{primal,w}(1/2)-I^0_{matching,w}(1/2)=0                   \tag{3.1}
\]

for every finite width.  Therefore every contribution that is ODD under primal/matching sector exchange must cancel in the sum of finite-size corrections.

In particular, if the square-lattice `w^-17/4` mismatch is produced by a sector-odd spin-four thermal-family anisotropy, its analogue has exactly zero amplitude in the self-matching triangular-site problem.  But the statement is stronger: **all** sector-odd amplitudes vanish, regardless of spin.

This makes triangular site a control for sector parity, not by itself a way to distinguish spin four from a scalar odd operator.

## 4. Three-lattice discrimination strategy

The clean comparison is therefore:

| lattice/model | rotation | self-matching? | expected charge-shift pattern |
|---|---|---|---|
| square site | `C4` | no | spin-4 odd correction allowed; observed `Delta=4` |
| kagome bond / hexagonal Bravais symmetry | `C6` (Jacobsen discusses three-fold basis symmetry) | no | spin-4 forbidden; observed leading `Delta=6` |
| triangular site | `C6` | yes | every sector-odd amplitude zero; root exactly `1/2` |

This separates two mechanisms that would otherwise be conflated:

1. **rotational selection**, which removes particular spins;
2. **exact primal/matching self-duality**, which removes the entire odd sector.

A scalar `x=21/4` explanation of the square `Delta=4` has no natural account of the square/kagome contrast, whereas the spin-four hypothesis does.

## 5. Symmetry-breaking experiment

The sharpest future test would start from a six-fold non-self-matching model (kagome bond is the existing literature example) and add a controlled anisotropy that lowers rotational symmetry while retuning to criticality.

The spin-four hypothesis predicts that once `C6/C3` no longer forbids spin four, a `Delta=4` charge-root correction should reappear with amplitude proportional to the appropriate spin-four component of the anisotropy for weak perturbations.

By contrast, changing a self-matching triangular-site model while preserving exact self-matching would still keep the entire sector-odd charge difference zero.  To activate the test there one must break self-matching as well as rotational symmetry.

## 6. Claim boundary

The finite root `1/2` and fixed-width charge-free-energy identity are exact consequences of self-matching plus complement symmetry.  The square/kagome/triangular comparison as evidence for a spin-four thermal-family correction is a research interpretation.  No new triangular transfer computation is required to establish the exact zero control.
