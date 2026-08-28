# Orientation Fourier selection rules for square-lattice corrections

Status: analytic scaffold for interpreting Gaussian-integer torus data.

## Symmetry before CFT

Consider a continuum square torus of fixed modulus `tau=i`, while the microscopic square lattice is rotated by angle `theta` relative to the torus periods.

Let `A(theta)` be the amplitude of any scalar finite-size correction in a reflection-symmetric observable such as the matching function evaluated at a fixed thermal coordinate.

The microscopic square lattice has `C4` rotational symmetry, so

\[
A(\theta+\pi/2)=A(\theta).
\]

Reflection implies

\[
A(-\theta)=A(\theta).
\]

Therefore, without invoking any specific field theory,

\[
\boxed{
A(\theta)=A_0+\sum_{m\ge1} A_{4m}\cos(4m\theta).
}
\]

Sine harmonics are forbidden by reflection. Harmonics not divisible by four are forbidden by microscopic `C4` symmetry.

Thus `cos(4 theta)` is not merely a heuristic fit. It is the **first symmetry-allowed anisotropic harmonic**.

## Connection to conformal spin

A continuum irrelevant field with conformal spin `s=h-\bar h` transforms under a microscopic rotation by

\[
\Phi_s\mapsto e^{is\theta}\Phi_s.
\]

A real reflection-even combination gives a contribution proportional to

\[
\cos(s\theta).
\]

Compatibility with the square lattice requires `s` to be a multiple of four for a scalar lattice observable. Therefore a leading `cos(4 theta)` harmonic is naturally interpreted as a spin-4 correction, while `cos(8 theta)` can arise from a spin-8 field, products of spin-4 fields, or nonlinear mixing.

This identification is stronger if the **finite-size exponent** of the same harmonic agrees with the scaling dimension of the proposed spin-4 field.

## Two independent axes of information

For each orientation-dependent term fit both:

1. angular harmonic `cos(4m theta)`;
2. size exponent `L^{-omega_m}`.

A convincing spin-4 identification needs both angular and radial evidence.

For example, if

\[
M_L(p_c;\theta)
=L^{-13/4}\left[A_0+A_4\cos4\theta\right]+
L^{-q}B(\theta)+\cdots,
\]

then the root bias is

\[
p_L^*(\theta)-p_c
\simeq -\frac{A_0+A_4\cos4\theta}{B_t}\,L^{-4}
\]

when `M'_L(p_c) ~ B_t L^{3/4}`.

The earlier axis/diamond sign flip suggests `|A_4| > |A_0|` on small sizes, but the Gaussian-orientation experiment is required to determine whether the asymptotic `A_0` vanishes, remains nonzero, or is merely smaller.

## Critical distinction: pure spin-4 versus mixed scalar + spin-4

There are three qualitatively different possibilities.

### H4-pure

\[
A_0=0,
\qquad A(\theta)=A_4\cos4\theta+\cdots.
\]

Consequences:

- axis and 45-degree amplitudes become exactly equal and opposite asymptotically;
- the simple axis/diamond average removes the leading correction;
- `theta=pi/8` is an asymptotic zero-amplitude geometry.

### H4-mixed

\[
A_0\neq0,
\qquad |A_4|\gtrsim |A_0|.
\]

Consequences:

- axis/diamond may still have opposite signs;
- their magnitudes need not match;
- the zero occurs at

\[
\cos4\theta_*=-A_0/A_4,
\]

not necessarily at `pi/8`;
- fixed-N angular tomography can estimate `theta_*` without knowing `p_c` precisely.

### H4-false

The leading angular behavior is not stable `cos4`, or its size exponent is inconsistent with the claimed leading correction. Then the axis/diamond sign flip is preasymptotic or produced by a different mechanism.

## Threshold-free extraction of A4

For two orientations at identical `N`, `tau`, and physical scale,

\[
D_{12}(p)=M_{\theta_1}(p)-M_{\theta_2}(p).
\]

The orientation-even thermal displacement and scalar correction cancel to leading order. Under the Fourier expansion,

\[
D_{12}
=L^{-13/4}A_4\left(\cos4\theta_1-\cos4\theta_2\right)
+L^{-13/4}A_8\left(\cos8\theta_1-\cos8\theta_2\right)+\cdots.
\]

This is why same-N Gaussian-integer pairs are the cleanest spin-analysis experiment.

## Fixed-N discrete Fourier regression

For an `N` with `r>=3` inequivalent sum-of-two-squares representations, construct design columns

\[
1,\quad \cos4\theta,\quad \cos8\theta,\ldots
\]

and fit only as many harmonics as the number of independent angles permits.

Recommended hierarchy:

1. constant only;
2. constant + `cos4`;
3. constant + `cos4` + `cos8`.

Use leave-one-orientation-out prediction when `r>=4` rather than relying on in-sample residuals.

The first useful `N` with four primitive orientations in the current design list is `N=1105`.

## Estimating the magic angle rather than assuming pi/8

If `constant + cos4` is accepted, estimate

\[
r_0=-A_0/A_4.
\]

A leading-amplitude zero exists only when `|r_0|<=1`, with

\[
\theta_*={1\over4}\arccos r_0
\]

in the fundamental interval `0<=theta<=pi/4`.

Freeze `theta_*` using smaller same-N data, then choose rational Gaussian-integer approximants to `tan theta_*` and test larger systems. This is a more rigorous successor to the preselected `pi/8` sequence.

## Why orientation and modular shape must be separated

Rotating the microscopic lattice at fixed `tau=i` probes lattice anisotropy / conformal-spin sectors.

Changing `tau` while keeping microscopic orientation fixed probes continuum shape dependence of torus one-point functions and scaling functions.

A `cos4 theta` zero and a modular-form zero such as a possible `E4(tau)` structure are conceptually distinct. They may interact, but fitting both simultaneously before each has been isolated would be underdetermined.
