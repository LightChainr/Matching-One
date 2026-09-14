# The first charge-sector relaxation gap appears to be 2 pi / w

2026-09-14.  Small-width spectral control for the corrected finite-aspect charge crossover.

## 1. Observable

At the fixed-width charge-coexistence point `p_w^ch`, let

\[
\lambda_{1,G,w}>|\lambda_{2,G,w}|\ge\cdots                  \tag{1.1}
\]

be the leading spectral moduli of the transparent safe-homology transfer kernel for graph `G`.  Define the dimensionless relaxation gap

\[
\boxed{
g_{G,w}
=w\log\frac{\lambda_{1,G,w}}{|\lambda_{2,G,w}|}.}             \tag{1.2}
\]

If critical transfer scaling applies, `g_{G,w}` should converge to `2 pi Delta x`, where `Delta x` is the first scaling-dimension gap inside the corresponding magnetic/topological sector.

## 2. Direct spectrum

Sparse Arnoldi evaluation of the same safe kernels used for the charge roots gives

| `w` | `g_4,w` | `g_8,w` |
|---:|---:|---:|
| 4 | 7.1077281 | 5.5172789 |
| 5 | 6.7138550 | 5.9713148 |
| 6 | 6.5406250 | 6.1426213 |
| 7 | 6.4512010 | 6.2165936 |
| 8 | 6.4004937 | 6.2516436 |

The two sequences approach from opposite sides and strongly bracket

\[
\boxed{2\pi=6.283185307179586\ldots}.                        \tag{2.1}
\]

The simplest interpretation is

\[
\boxed{\Delta x=1.}                                          \tag{2.2}
\]

That is exactly the spacing of a level-one conformal descendant above a primary magnetic state on the cylinder.

This is a much less speculative identification than the separate Kac-(4,2) hypothesis for the **sector-odd irrelevant correction**: the present gap concerns ordinary relaxation *within* each magnetic sector, not the tiny difference between the two sectors.

## 3. Consequence for finite aspect

For a periodic transfer trace, subleading contamination at aspect ratio

\[
\rho=m/w                                                       \tag{3.1}
\]

is therefore expected to scale as

\[
\boxed{e^{-2\pi\rho}}                                        \tag{3.2}
\]

up to amplitudes and further gaps.

Combining this with the charge thermal slope

\[
\Theta'_w\asymp w^{-1/4}                                     \tag{3.3}
\]

and the derivative factor `m`, the corrected finite-length root displacement has scale

\[
\boxed{
|p^*_{w,m}-p_w^{ch}|
\sim \frac{e^{-2\pi\rho}}{\rho\,w^{3/4}}}                    \tag{3.4}
\]

unless an additional symmetry cancels the first descendant contribution from the **difference** of the two periodic traces.

Equation (3.4) is a falsifiable target, not a theorem of the present safe kernel alone.

## 4. When is the intrinsic w^-4 shift visible?

The semi-infinite pseudo-critical displacement is

\[
p_c-p_w^{ch}\asymp w^{-4}.                                   \tag{4.1}
\]

Demanding the finite-length term (3.4) be asymptotically smaller gives

\[
\frac{e^{-2\pi\rho}}\rho\ll w^{-13/4}.                       \tag{4.2}
\]

If

\[
\rho=C\log w,                                                 \tag{4.3}
\]

then a sufficient leading-power condition is

\[
\boxed{C>\frac{13}{8\pi}=0.517253\ldots.}                    \tag{4.4}
\]

So only logarithmic aspect growth may be needed to expose the semi-infinite `w^-4` correction on periodic tori.

This is qualitatively different from a generic open-boundary amplitude mismatch, which would give a `1/m` correction without the exponential factor.

## 5. Possible extra cancellation

The primal and complementary matching sectors approach the same continuum magnetic module.  It is therefore possible that not only their leading primary eigenvalues but also the first descendant trace amplitudes agree in the **sector difference** at criticality.  If so, the actual periodic-root correction would begin at `e^{-4 pi rho}` or another higher gap.

This should be tested directly by constructing the exact periodic topological transfer blocks or by comparing finite-`m` roots at fixed `w`.  The current safe-kernel spectrum shows that a `2 pi/w` relaxation mode exists; it does not by itself prove that this mode survives subtraction of the two torus traces.

## 6. New hierarchy of spectral questions

The charge programme now separates three distinct spectra.

1. **Primary magnetic gap:**
   \[
   w I^0_G(p_c)\to2\pi(5/48).
   \]
2. **Within-sector relaxation:**
   \[
   w\log(\lambda_1/|\lambda_2|)\to2\pi.
   \]
3. **Primal/dual sector-odd mismatch:**
   conjecturally
   \[
   I^0_4(p_c)-I^0_8(1-p_c)\asymp w^{-17/4},
   \]
   with the tentative `(4,2)` LCFT interpretation recorded separately.

Confusing these three is exactly what makes the fast `w^-4` pseudo-critical convergence look mysterious.

## 7. Claim boundary

The finite spectral numbers in section 2 are direct outputs of the transparent safe transfer.  The limit `2 pi`, the descendant interpretation, and the use of that gap in the exact periodic topological trace difference are scaling conjectures to be checked at larger widths / with the precise periodic sector transfer.  No claim is made that the safe-kernel second eigenvalue is already a certified eigenvalue of Jacobsen's reduced TL block, although the matching leading eigenvalue and root sequence strongly support the shared sector semantics.
