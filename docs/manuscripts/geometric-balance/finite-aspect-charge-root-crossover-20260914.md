# From the thermal w^-3/4 window to the semi-infinite w^-4 charge root

2026-09-14.  Corrected two-parameter scaling note for #767.

A first draft treated the two charged sectors as generic boundary-to-boundary Perron amplitudes and obtained an `O(1/m)` root displacement.  That is a valid generic transfer situation but is **not the right default for the periodic torus matching observable**.  At `q=1`, the graph-polynomial/topological sectors are periodic transfer traces.  A simple Perron eigenvalue appears in a trace with coefficient one, so at the semi-infinite eigenvalue crossing the two leading trace terms cancel exactly.  Finite-`m` convergence is then controlled by subleading spectral gaps and is exponential in `m` at fixed `w`, as already observed in Jacobsen's finite-`m` calculations.

This correction materially changes the aspect scale needed to see the intrinsic `w^-4` pseudo-critical shift.

## 1. Semi-infinite charge crossing

For circumference `w`,

\[
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p),                     \tag{1.1}
\]

and

\[
p_w^{ch}:\quad\Theta_w(p_w^{ch})=0.                         \tag{1.2}
\]

The transparent site transfer reproduces the Jacobsen `w x infinity` root sequence.  Critical scaling gives

\[
\Theta'_w(p_w^{ch})\asymp w^{-1/4},                           \tag{1.3}
\]

while

\[
p_c-p_w^{ch}\asymp w^{-4}.                                  \tag{1.4}
\]

## 2. Why a generic amplitude formula would give 1/m

For two arbitrary open-boundary transfer observables one can have

\[
Z_i(m)=A_i\lambda_i^m[1+o(1)].                               \tag{2.1}
\]

At `lambda_1=lambda_2`, unequal `A_i` produce a root displacement

\[
\delta p\sim\frac{\log(A_2/A_1)}{m\Theta'_w}.                \tag{2.2}
\]

This remains a useful warning for nonperiodic boundary conditions or modified observables.  It was the source of the discarded `m~w^(17/4)` estimate in the first draft.

## 3. Periodic torus: trace coefficients cancel

For the actual periodic topological matching observable, the graph-polynomial/eigenvalue formulation decomposes the periodic transfer into open and closed topological blocks.  At `q=1`, schematically

\[
Z_{2D}=\operatorname{Tr}(T_{open}^m),\qquad
Z_{0D}=\operatorname{Tr}(T_{closed}^m),                       \tag{3.1}
\]

up to the common local normalization converting Boltzmann weights to Bernoulli probabilities.

Let the dominant eigenvalues be `lambda_o,lambda_c`, and write the next spectral moduli as `mu_o,mu_c`.  For simple Perron roots,

\[
\operatorname{Tr}(T_i^m)
=\lambda_i^m\left[1+O\left((|\mu_i|/\lambda_i)^m\right)\right].\tag{3.2}
\]

The leading coefficient is **one**, not an arbitrary boundary overlap.

At `p=p_w^{ch}`,

\[
\lambda_o=\lambda_c=\lambda_w,                               \tag{3.3}
\]

so the leading terms in `Z_{2D}-Z_{0D}` cancel exactly.  Define

\[
\Delta_w
=\min_i\log\frac{\lambda_w}{|\mu_i|}>0.              \tag{3.4}
\]

Then the finite-`m` mismatch at the semi-infinite root is

\[
Z_{2D}-Z_{0D}
=\lambda_w^m O(e^{-m\Delta_w}).                              \tag{3.5}
\]

Differentiating the leading eigenvalue difference contributes the factor

\[
m\Theta'_w(p_w^{ch}).                                        \tag{3.6}
\]

Therefore the natural periodic-root correction is

\[
\boxed{
|p^*_{w,m}-p_w^{ch}|
=O\left(\frac{e^{-m\Delta_w}}
{m\Theta'_w(p_w^{ch})}\right),}                              \tag{3.7}
\]

provided the subleading eigenvalues remain separated and no equal-modulus oscillatory degeneracy changes the prefactor.  This is exponential in `m` for every fixed `w`.

The exact site-rank transfer should be audited directly if (3.7) is promoted as a theorem; the key correction here is that periodic trace structure removes the generic `O(1/m)` amplitude term.

## 4. Critical scaling of the subleading gap

At criticality a transfer spectral gap within a fixed CFT sector has the form

\[
\Delta_w\sim\frac{g}{w},                                     \tag{4.1}
\]

where `g=2 pi Delta x` times the lattice velocity/geometry factor for the relevant next state.  Put

\[
\rho=m/w.                                                     \tag{4.2}
\]

Combining (1.3), (3.7), and (4.1),

\[
\boxed{
|p^*_{w,m}-p_w^{ch}|
\lesssim
\frac{e^{-g\rho}}{\rho\,w^{3/4}}}                           \tag{4.3}
\]

at the level of critical finite-size scaling.

This is the corrected aspect crossover.

## 5. Fixed aspect recovers the ordinary thermal window

If `rho` is fixed,

\[
|p^*_{w,m}-p_w^{ch}|
=O(w^{-3/4})                                                   \tag{5.1}
\]

with an aspect-dependent coefficient `e^{-g rho}/rho`.

Since the intrinsic semi-infinite displacement is only `w^-4`, the finite-aspect root is governed at leading order by the standard near-critical thermal scale.  Thus a square or fixed-aspect torus should not be expected to reveal the Jacobsen `w^-4` shift directly.

## 6. Growing aspect: only logarithmic growth may be needed

To resolve the semi-infinite width shift before finite-length corrections dominate, require

\[
\frac{e^{-g\rho}}{\rho w^{3/4}}\ll w^{-4},                   \tag{6.1}
\]

or

\[
\boxed{
\frac{e^{-g\rho}}\rho\ll w^{-13/4}.}                        \tag{6.2}
\]

Thus a logarithmically growing aspect ratio

\[
\rho=C\log w                                                   \tag{6.3}
\]

is already sufficient when

\[
gC>13/4                                                       \tag{6.4}
\]

(up to logarithmic factors and the actual sector gap constant).

This replaces the discarded polynomial requirement `rho >> w^(13/4)`.  The semi-infinite regime is reached much sooner because periodic traces cancel their leading amplitudes.

## 7. Exponential aspect

When

\[
\log m/w\to d>0,                                              \tag{7.1}
\]

`rho` itself is exponentially large, so finite-length transfer corrections to the charge root are beyond all algebraic orders in `w`.  The fixed-width charge crossing is then effectively exact long before the lower/upper winding-component birth windows are considered.

This cleanly coexists with the main geometric phenomenon: the matching root can track a very sharp charge-sector eigenvalue crossing while the two individual rank births converge to macroscopically separated parameters.

## 8. Interface to #767

The corrected picture gives three compatible resolutions of the same charge field.

1. **Fixed aspect `rho=O(1)`:** root motion occurs on the standard `w^-3/4` near-critical scale.
2. **Growing aspect:** finite-length corrections are additionally suppressed by the transfer factor `e^{-g rho}`.
3. **Semi-infinite / exponential aspect:** the intrinsic charge root is `p_w^{ch}=p_c+O(w^-4)` and finite-length corrections are negligible.

The crossover is controlled by a **spectral gap in aspect ratio**, not by an arbitrary Perron amplitude ratio.

## 9. A new concrete target

The next useful computation is not a larger threshold table.  It is the subleading spectrum of the transparent safe/topological transfer at `p_w^{ch}`:

\[
g_w=w\Delta_w.                                                \tag{9.1}
\]

If `g_w` approaches a nonzero limit, then (4.3) becomes quantitatively testable.  Comparing that limit with CFT candidate gaps will identify which descendant/excitation controls approach to the semi-infinite charge criterion.

## 10. Claim boundary

- The generic amplitude formula (2.2) is correct for arbitrary boundary overlaps but is **not** used as the periodic-torus conclusion.
- Exponential fixed-`w` convergence follows from periodic trace blocks with simple separated Perron eigenvalues; Jacobsen's eigenvalue formulation supplies the corresponding graph-polynomial structure and finite-`m` controls.
- The simultaneous scaling `Delta_w~g/w` and `Theta'_w~w^-1/4` is CFT/critical-scaling input, not yet a rigorous square-site theorem.
- The earlier `m~w^(17/4)` crossover claim is withdrawn for the periodic rank observable.
