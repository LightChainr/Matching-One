# Pell-paired axis/diamond estimator

Status: high-priority, low-tuning experiment.

## Core observation

The axis square torus of integer side `a` has

\[
N_A=a^2,
\qquad L_A=a.
\]

The `pi/4`-rotated diamond torus with integer parameter `d` and periods

\[
(d,d),(-d,d)
\]

has

\[
N_D=2d^2,
\qquad L_D=\sqrt2\,d.
\]

Choose `(a,d)` from Pell's equation

\[
a^2-2d^2=\pm1.
\]

Then simultaneously

\[
|N_A-N_D|=1
\]

and

\[
L_A/L_D\to1
\]

extremely rapidly.

The first useful pairs are

\[
(3,2),\ (7,5),\ (17,12),\ (41,29),\ (99,70),\ldots
\]

Thus we get two square tori with almost identical physical scale and area, but with the microscopic square lattice rotated by `pi/4` between them.

This is unusually clean for testing a spin-4 correction: size/area matching is fixed by number theory rather than fitted after seeing the data.

## Hypothesis

Let `p_A(a)` and `p_D(d)` be the unique roots of the same finite matching function `M(p)=0` in the axis and diamond quotients.

Assume the leading orientation-sensitive root bias is a spin-4 contribution,

\[
p_A(a)=p_c-Ca^{-w}+\cdots,
\]

\[
p_D(d)=p_c+C(\sqrt2 d)^{-w}+\cdots,
\]

with `w` expected near 4 for the ordinary matching root.

Under the `pi/4` rotation, the sign change follows the spin-4 factor `exp(i4 pi/4)=-1`; equality of magnitudes is the nontrivial part to test.

## Why Pell pairing improves cancellation

Write `L=\sqrt2 d`. Pell gives

\[
a^2=L^2\pm1,
\]

so

\[
a=L\left(1\pm L^{-2}\right)^{1/2}
=L\left[1\pm\tfrac12L^{-2}+O(L^{-4})\right].
\]

Therefore

\[
a^{-w}-L^{-w}=O(L^{-w-2}).
\]

The simple orientation mean

\[
\mu(a,d)=\frac{p_A(a)+p_D(d)}{2}
\]

cancels an exactly equal-and-opposite leading amplitude down to the scale-mismatch residual

\[
O(L^{-w-2}).
\]

If `w=4`, the *Pell mismatch itself* only leaves an `L^-6` contribution from the spin-4 sector.

This does not remove rotation-even or logarithmic corrections. Those become the next object to identify.

## Threshold-free diagnostic

The orientation gap

\[
\Delta(a,d)=p_D(d)-p_A(a)
\]

should satisfy

\[
\Delta\sim2C L^{-w}
\]

if the sign-flipped field dominates.

This is important: the exponent and amplitude can be tested **without inserting any assumed value of `p_c`**.

Across successive Pell pairs, estimate

\[
w_{\rm eff}
=-\frac{\log|\Delta_{k+1}/\Delta_k|}
{\log(L_{k+1}/L_k)}.
\]

Only after the orientation law is established should `mu_k` be used as an accelerated threshold estimator.

## Fixed-power exact cancellation

If `w=4` is preregistered from independent evidence, one can cancel the sign-flipped term exactly even though `L_A != L_D` by using

\[
\widehat p_{4}
=w_Ap_A+w_Dp_D,
\]

with

\[
w_A=\frac{L_A^4}{L_A^4+L_D^4},
\qquad
w_D=\frac{L_D^4}{L_A^4+L_D^4}.
\]

No amplitude is fitted. Only the transformation law and exponent are assumed.

The simple mean is scientifically cleaner as a first test; the fixed-power weights are a second, preregistered estimator.

## Existing exact tiny-size evidence

Our brute-force reference data already include the Pell pair `(a,d)=(3,2)`:

| geometry | parameter | sites | physical side | matching root |
|---|---:|---:|---:|---:|
| axis | `a=3` | 9 | 3 | `0.586511455112676...` |
| diamond | `d=2` | 8 | `2 sqrt(2)` | `0.604563277853507...` |

The two roots straddle the accepted threshold neighborhood. This pair is much too small for an asymptotic cancellation claim, but it confirms that the implementation supplies the desired opposite-sign finite-size bias.

A second, non-Pell comparison `(axis a=4, diamond d=3)` gives roots

```text
axis:    0.590672112331028...
diamond: 0.594252321168569...
```

again on opposite sides.

Using a reference threshold only as a diagnostic, the apparent `L^4`-scaled biases for this second pair are approximately

```text
axis:    -0.531
diamond: +0.488
```

which are already surprisingly close in magnitude at very small sizes. This numerical observation is exploratory, not an asymptotic fit.

## Production sequence

Use these Pell pairs in increasing order:

| axis `a` | diamond `d` | `N_A` | `N_D` | use |
|---:|---:|---:|---:|---|
| 3 | 2 | 9 | 8 | exact regression |
| 7 | 5 | 49 | 50 | CPU exact/very-high-statistics if feasible |
| 17 | 12 | 289 | 288 | CPU Monte Carlo |
| 41 | 29 | 1681 | 1682 | CPU Monte Carlo |
| 99 | 70 | 9801 | 9800 | CPU/GPU |
| 239 | 169 | 57121 | 57122 | GPU later |
| 577 | 408 | 332929 | 332928 | GPU later |

The site-count difference stays exactly one.

## Acceptance criteria

Call the spin-4/Pell mechanism supported only if:

1. `Delta_k` has a stable sign and power law over several Pell pairs;
2. `w_eff` approaches the independently expected exponent rather than being fitted arbitrarily;
3. axis and diamond scaled amplitudes approach equal magnitude with opposite sign after physical-length normalization;
4. the simple mean `mu_k` converges predictively faster on held-out larger pairs;
5. the same behavior appears in an exactly solved square-bond control;
6. a higher-symmetry control does not show the same square-specific orientation pattern.

## Potential payoff

If this works, it gives a finite-size acceleration mechanism with almost no continuous fitting:

- topology chooses the observable;
- rotation chooses the sign;
- Pell arithmetic matches the physical sizes;
- only the asymptotic operator content remains to be identified.

That is much stronger evidence than finding another inverse-power polynomial that happens to fit one threshold sequence.
