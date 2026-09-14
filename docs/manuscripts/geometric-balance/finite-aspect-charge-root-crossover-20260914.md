# From the thermal w^-3/4 window to the semi-infinite w^-4 charge root

2026-09-14.  Conditional two-parameter scaling consequence of the fixed-width charge free energy.  This note links the fixed-aspect near-critical regime in #767 to the `m->infinity` Jacobsen/charge-coexistence root.

The exponent algebra is robust; the sharp finite-`m` statement requires control of the two Perron sector amplitudes, which is kept explicit below.

## 1. Long-cylinder charge fugacity with amplitudes

At fixed circumference `w`, analytic finite-state transfer gives, away from accidental subleading degeneracies,

\[
P_0^{(4)}(p;w,m)
=A_{4,w}(p)\,e^{-m I^0_{4,w}(p)}[1+O(e^{-m\Delta_{4,w}})],    \tag{1.1}
\]

and by digital Alexander

\[
P_2^{(4)}(p;w,m)
=P_0^{(8)}(1-p;w,m)
=A_{8,w}(1-p)\,e^{-m I^0_{8,w}(1-p)}[1+O(e^{-m\Delta_{8,w}})].\tag{1.2}
\]

Therefore the charge log-odds is

\[
\theta_{w,m}(p)
=m\Theta_w(p)
+\log\frac{A_{8,w}(1-p)}{A_{4,w}(p)}+o(1),                   \tag{1.3}
\]

with

\[
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p).                     \tag{1.4}
\]

The semi-infinite charge root `p_w^ch` satisfies `Theta_w(p_w^ch)=0`.

## 2. Finite-length displacement from the semi-infinite root

Assume `Theta'_w(p_w^ch)>0` and write

\[
R_w=\log\frac{A_{8,w}(1-p_w^{ch})}{A_{4,w}(p_w^{ch})}.        \tag{2.1}
\]

Linearizing (1.3), the finite-torus matching root obeys

\[
\boxed{
p^*_{w,m}-p_w^{ch}
=-\frac{R_w}{m\Theta'_w(p_w^{ch})}
+o\!\left(\frac{|R_w|+1}{m\Theta'_w}\right).}               \tag{2.2}
\]

If `R_w=O(1)`, the width of the charge-sign crossover is

\[
\boxed{\delta p_{charge}\asymp\frac1{m\Theta'_w}.}           \tag{2.3}
\]

If the amplitude ratio is only polynomial in `w`, then `R_w=O(log w)` and the same formula acquires only a logarithm.

## 3. Critical scaling of the slope

The transfer control and the pivotal identity give

\[
\Theta'_w(p_w^{ch})\asymp w^{-1/4},                           \tag{3.1}
\]

which is the ordinary thermal scaling `y_t-1` with `y_t=3/4`.  Hence

\[
\boxed{
\delta p_{charge}\asymp\frac{w^{1/4}}m.}                    \tag{3.2}
\]

For aspect ratio

\[
\rho=m/w,                                                     \tag{3.3}
\]

this is

\[
\boxed{
\delta p_{charge}\asymp\frac{w^{-3/4}}\rho.}                \tag{3.4}
\]

This is the natural interpolation between the fixed-aspect thermal window and the semi-infinite eigenvalue crossing.

## 4. Three aspect regimes

### Fixed macroscopic aspect ratio

If `m~rho w` with fixed positive `rho`,

\[
\delta p_{charge}\asymp w^{-3/4}.                             \tag{4.1}
\]

So the balance root lives in the standard near-critical thermal window.  The extremely accurate semi-infinite pseudo-critical point is not the relevant finite-size centre at this order.

### Growing but sub-semi-infinite aspect

If `rho->infinity`,

\[
\delta p_{charge}\asymp w^{-3/4}/\rho.                       \tag{4.2}
\]

The charge crossing sharpens continuously as the cylinder grows longer.  This is a natural intermediate regime for #767.

### Exponential aspect

If `log m / w -> d>0`, then

\[
\frac{w^{1/4}}m
\]

is exponentially small.  Finite-length amplitude corrections to the charge root are negligible compared with any algebraic width correction.  This explains why the balance root can lock onto the fixed-width charge free-energy crossing even while the two rank births split to distinct fixed parameters.

## 5. When does the Jacobsen w^-4 shift become visible on a finite torus?

The semi-infinite sequence satisfies

\[
p_c-p_w^{ch}\asymp w^{-4}.                                   \tag{5.1}
\]

To resolve this intrinsic width shift before finite-length displacement dominates, require

\[
\frac{w^{1/4}}m\ll w^{-4}.                                   \tag{5.2}
\]

Thus

\[
\boxed{m\gg w^{17/4}.}                                       \tag{5.3}
\]

Equivalently

\[
\boxed{\rho=m/w\gg w^{13/4}.}                                \tag{5.4}
\]

If `R_w=O(log w)` rather than `O(1)`, insert the corresponding logarithmic factor.

This is a strong warning for numerical comparisons: a width sequence on square or mildly elongated tori should not be expected to display the semi-infinite `w^-4` estimator shift cleanly.  One must first be beyond the much larger finite-length charge window.

## 6. Polynomial aspect phase diagram

Let

\[
m=w^a.                                                        \tag{6.1}
\]

Ignoring possible logarithmic amplitude factors, (3.2) gives

\[
p^*_{w,m}-p_w^{ch}=O(w^{1/4-a}).                             \tag{6.2}
\]

Therefore:

- `a=1`: ordinary `w^-3/4` thermal scale;
- `1<a<17/4`: finite-length displacement is smaller than the thermal window but still larger than the intrinsic `w^-4` shift;
- `a=17/4`: the two corrections compete;
- `a>17/4`: the semi-infinite `w^-4` shift dominates.

This gives a concrete meaning to “take `m->infinity` first”: the aspect need not literally be infinite, but it must grow past a high polynomial crossover if the `w^-4` pseudo-critical structure is the target.

## 7. Relation to #767 variables

The common-window programme in `charge-neutral-crossover-coordinates-20260914.md` uses the charge field `theta` rather than two independent colour intensities.  Equation (1.3) supplies its large-aspect matching form:

\[
\theta\approx m\Theta_w+R_w.                                 \tag{7.1}
\]

Thus a crossover theory should match

1. the fixed-aspect thermal scaling of `theta` at `m=O(w)`;
2. the linear-in-`m` charge free energy for `m/w->infinity`;
3. the separated birth-component Poisson windows in exponential aspect.

These are not three unrelated ansatzes; they are different resolutions of the same bounded topological charge.

## 8. Claim boundary

The fixed-`w` Perron expansion and the definition of `Theta_w` are standard finite-state transfer facts.  The power `Theta'_w~w^-1/4` is currently an empirically strong/CFT-and-pivotal-supported scaling statement for square-site percolation.  The crossover threshold `m~w^(17/4)` additionally assumes the sector amplitude log-ratio grows at most subleadingly (for example `O(log w)`).  No uniform theorem in the simultaneous `w,m->infinity` limit is claimed here.
