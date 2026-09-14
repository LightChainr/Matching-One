# Exponential aspect should lock the finite balance root to the semi-infinite charge sequence

2026-09-14.  Conditional finite-size prediction combining the directional birth theory with the fixed-width charge-sector transfer.

The conclusion is deliberately counterintuitive:

> at fixed positive exponential aspect parameter `d`, the two individual homology births move to `d`-dependent macroscopic centres, but the finite matching/balance root should have the SAME leading `w^-4` width correction as the semi-infinite Jacobsen sequence, independent of `d`.

## 1. Two different centres coexist

For axial square tori

\[
C_w\times C_m,
\qquad
\frac{\log m}{w}\to d>0,                                    \tag{1.1}
\]

the individual NN rank births satisfy

\[
T_1\to a(d),
\qquad
T_2\to b(d),                                                  \tag{1.2}
\]

with

\[
\kappa_4(a(d))=d,
\qquad
\kappa_8(1-b(d))=d.                                          \tag{1.3}
\]

The two centres can stay a fixed positive distance apart and move strongly with `d`.

The matching root, however, is the zero of

\[
P_2(p)-P_0(p),                                                \tag{1.4}
\]

or equivalently of the charge fugacity.  At fixed width and `m->infinity`, that root tends to

\[
p_w^{ch}:\quad
\lambda^0_{4,w}(p_w^{ch})
=\lambda^0_{8,w}(1-p_w^{ch}).                                \tag{1.5}
\]

This is the Jacobsen semi-infinite-cylinder eigenvalue sequence.

## 2. Longitudinal finite-size correction

For the periodic torus, the two charged topological sectors are transfer traces.  At the semi-infinite crossing their leading eigenvalue contributions cancel.  If

\[
\Delta_w
=\min\log\frac{\lambda_1}{|\lambda_2|}               \tag{2.1}
\]

is the relevant within-sector relaxation gap, then the finite-length root correction has the form

\[
|p^*_{w,m}-p_w^{ch}|
\lesssim
\frac{e^{-m\Delta_w}}
{m\Theta'_w(p_w^{ch})}                                       \tag{2.2}
\]

up to subleading trace amplitudes and possible extra cancellations.

The transparent safe spectrum gives strong evidence for

\[
\Delta_w\sim\frac{2\pi}{w},                                  \tag{2.3}
\]

while the thermal slope obeys

\[
\Theta'_w\asymp w^{-1/4}.                                    \tag{2.4}
\]

## 3. Exponential aspect annihilates longitudinal corrections

Under (1.1),

\[
\frac mw
=\frac{e^{dw+o(w)}}w.                                \tag{3.1}
\]

Equations (2.2)--(2.4) then give schematically

\[
|p^*_{w,m}-p_w^{ch}|
\le
\exp[-c e^{dw+o(w)}]                                         \tag{3.2}
\]

up to polynomial factors in `w`.

In particular the finite-length error is smaller than every algebraic power of `w`:

\[
\boxed{
p^*_{w,m}=p_w^{ch}+o(w^{-K})
\quad\text{for every fixed }K.}                              \tag{3.3}
\]

This is conditional on a uniform `Delta_w>=c/w` spectral-gap bound (or an equivalent periodic-transfer mixing statement).  The current small-width spectrum strongly supports that hypothesis but does not prove it.

## 4. Consequence: the root shift forgets d at leading order

The semi-infinite square-site sequence has

\[
p_c-p_w^{ch}\asymp A_\square w^{-4},                         \tag{4.1}
\]

with the published Jacobsen widths and modern `p_c` reference giving a coefficient near `0.29` after the leading `1/w^2` extrapolation of `(p_c-p_w)w^4`.

Combining with (3.3),

\[
\boxed{
p_c-p^*_{w,m}
\sim A_\square w^{-4}}                                       \tag{4.2}
\]

for **every fixed `d>0`**, with the same leading amplitude as the semi-infinite cylinder.

Thus the aspect parameter is absent from the leading balance-root correction even though it completely controls the separated birth centres in (1.2)--(1.3).

## 5. A striking two-scale statement

At fixed `d>0`, the finite system should simultaneously exhibit

\[
T_1=a(d)+O_P(w^{-1}),
\qquad
T_2=b(d)+O_P(w^{-1}),                                        \tag{5.1}
\]

at regular mass points (Gumbel windows), while

\[
p^*_{w,m}=p_c-A_\square w^{-4}+o(w^{-4}).                   \tag{5.2}
\]

So the typical birth fluctuations are `1/w`, the birth-centre separation is `O(1)`, and the balance-root displacement is only `w^-4`.

This is a particularly sharp version of **balance without concentration**:

```text
birth locations        : geometry-sensitive, d-dependent, O(1) apart
birth fluctuations     : O(1/w)
balance root correction: d-independent to leading order, O(w^-4)
```

These are three different spectral/probability mechanisms.

## 6. Why d=0 is singular for this statement

The claim is for fixed positive `d`.  At fixed macroscopic aspect ratio (`m=O(w)`), longitudinal trace corrections are only `e^{-O(1)}` and the root lies in the ordinary `w^-3/4` near-critical thermal window before the much smaller semi-infinite `w^-4` correction can be isolated.

Therefore one must not take the `d>0` prediction continuously to `d=0` without resolving the simultaneous aspect crossover.

The order/scale distinction is exactly the one encoded in `finite-aspect-charge-root-crossover-20260914.md`.

## 7. Falsification test using existing-style archives

No new large critical simulation is needed in principle.  For several fixed `d>0` values and available widths, compute the exact/MC balance roots and form

\[
R_w(d)=w^4[p_c-p^*_{w,m(d)}].                                \tag{7.1}
\]

The prediction is

\[
\boxed{R_w(d_1)-R_w(d_2)\to0}                                \tag{7.2}
\]

for any two positive fixed `d_1,d_2`, and all should approach the same semi-infinite coefficient.

By contrast the lower/upper birth medians should visibly separate by different `d`-dependent amounts.  Measuring both in the same data block would be an especially clean control.

## 8. Claim boundary

The fixed-`d` birth centres are author-level results already on the branch.  The semi-infinite root sequence and its `w^-4` correction are established numerical/literature facts of the eigenvalue method.  Equation (4.2) additionally assumes a periodic topological-sector relaxation gap of order at least `1/w`; the transparent safe spectrum points strongly to the specific constant `2 pi`, but a uniform proof is still missing.
