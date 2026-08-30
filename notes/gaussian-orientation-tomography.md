# Gaussian-integer orientation tomography on an exactly square torus

Status: high-priority exploratory design. This strengthens the earlier axis/diamond Pell-pair experiment by eliminating the residual size mismatch entirely.

## Core construction

For the square lattice, choose torus period vectors

\[
v_1=(a,b),\qquad v_2=(-b,a),\qquad a,b\in\mathbb Z.
\]

They satisfy

\[
v_1\cdot v_2=0,\qquad |v_1|=|v_2|=\sqrt{a^2+b^2}.
\]

Therefore every such quotient is an **exact square torus** with continuum modulus

\[
\tau=i,
\]

physical side length

\[
L_{\rm phys}=\sqrt{N},
\]

and number of lattice sites

\[
N=|\det(v_1,v_2)|=a^2+b^2.
\]

The only intended change is the orientation of the microscopic square lattice relative to the torus axes,

\[
\theta=\arctan(b/a).
\]

The axis and 45-degree diamond geometries are the endpoints `(a,b)=(L,0)` and `(d,d)`.

## Spin-4 prediction

If the leading orientation-sensitive correction is dominated by a spin-4 field, its amplitude at fixed torus modulus should have the angular form

\[
A(\theta)=A_0+A_4\cos(4\theta)+A_8\cos(8\theta)+\cdots.
\]

Reflection symmetry removes sine terms in the simplest version. The exact algebraic value of the first harmonic is

\[
\cos(4\theta)
=\frac{a^4-6a^2b^2+b^4}{(a^2+b^2)^2}.
\]

For the square-site matching function, a sharpened working model is

\[
M_{N,\theta}(p_c)
=N^{-13/8}\left[A_0+A_4\cos(4\theta)+\cdots\right],
\]

if the ordinary matching-root bias is asymptotically `L^-4` and `M'_L(p_c)~L^(3/4)`.

This model is a hypothesis, not an assumption to be forced in fits.

## Same-N comparisons: the decisive feature

Many integers have multiple representations as a sum of two squares. Then distinct microscopic orientations have **identical**:

- number of sites `N`,
- physical side length `sqrt(N)`,
- torus modulus `tau=i`,
- aspect ratio and area.

Examples with primitive representations include:

| N | representation 1 | theta 1 | cos(4 theta 1) | representation 2 | theta 2 | cos(4 theta 2) |
|---:|---|---:|---:|---|---:|---:|
| 65 | `(8,1)` | 7.125 deg | +0.8788166 | `(7,4)` | 29.745 deg | -0.4844970 |
| 85 | `(9,2)` | 12.529 deg | +0.6412457 | `(7,6)` | 40.601 deg | -0.9532180 |
| 145 | `(12,1)` | 4.764 deg | +0.9452081 | `(9,8)` | 41.634 deg | -0.9725089 |
| 205 | `(14,3)` | 12.095 deg | +0.6641999 | `(13,6)` | 24.775 deg | -0.1581678 |
| 425 | `(19,8)` | 22.834 deg | -0.02329135 | `(16,13)` | 39.094 deg | -0.9161910 |

Larger `N` can have four or more inequivalent primitive orientations. For example `N=1105` has `(33,4)`, `(32,9)`, `(31,12)`, `(24,23)`.

These comparisons are substantially cleaner than comparing axis `a` with diamond `d` at merely similar physical lengths.

## Threshold-independent orientation difference

At a common probability `p_ref` near `p_c`, define

\[
D_N(\theta_1,\theta_2;p_{ref})
=M_{N,\theta_1}(p_{ref})-M_{N,\theta_2}(p_{ref}).
\]

Because the two systems have identical continuum shape and physical scale, the common leading thermal displacement caused by `p_ref-p_c` should cancel to leading order. Under the spin-4 model,

\[
N^{13/8}D_N
\to A_4\,[\cos(4\theta_1)-\cos(4\theta_2)].
\]

Therefore the ratio

\[
\mathcal A_4(N)=
\frac{N^{13/8}D_N}
{\cos(4\theta_1)-\cos(4\theta_2)}
\]

should approach an orientation-independent constant without requiring the disputed last digits of `p_c`.

This is the preferred first test.

## Fixed-N angular regression

When one `N` has at least three inequivalent representations, fit at fixed `N`

\[
M_{N,\theta}(p)=\alpha_N+\beta_N\cos(4\theta)
\]

using no size extrapolation at all.

Then test whether adding `cos(8 theta)` is demanded out of sample. Across increasing `N`, inspect the scaling of `alpha_N` and `beta_N` separately.

This separates orientation-even and orientation-odd sectors more directly than an axis/diamond average.

## A magic orientation near 22.5 degrees

For a pure leading spin-4 correction,

\[
\cos(4\theta)=0
\]

at

\[
\theta=\pi/8=22.5^\circ,
\qquad \tan\theta=\sqrt2-1.
\]

Rational approximants to `sqrt(2)-1` give exact square tori whose spin-4 amplitude should be strongly suppressed. A convenient sequence is

```text
(a,b) = (5,2), (12,5), (29,12), (70,29), (169,70), ...
```

with

```text
N = 29, 169, 985, 5741, 33461, ...
```

and rapidly decreasing `|cos(4 theta)|`.

The first two are CPU-scale. If their matching-root bias is already anomalously small compared with nearby orientations at comparable `N`, this is a strong signal before any GPU campaign.

## Stronger prediction than a sign flip

The original axis/diamond result tests only

\[
A(0)\quad\text{versus}\quad A(\pi/4).
\]

A genuine spin-4 mechanism predicts the full angular law. The hypothesis should be considered substantially strengthened only if:

1. same-N orientation differences have the predicted signs;
2. normalized differences collapse when divided by `Delta cos(4 theta)`;
3. a `pi/8` approximant suppresses the leading amplitude;
4. fixed-N multi-angle data prefer `cos(4 theta)` over arbitrary orientation labels;
5. higher harmonics decrease systematically with size.

A mere axis/diamond sign flip is not sufficient.

## Exact-control models

Run the identical quotient geometries first for square **bond** percolation at the exact threshold `p_c=1/2`. This is a same-shape, same-lattice-symmetry positive control with no threshold uncertainty.

Then apply the frozen geometry conventions and angular analysis to square site percolation and its matching lattice.

## Computational notes

The quotient should be implemented for a general integer period matrix, not with special-case coordinate formulas for axis and diamond systems. The displacement-potential union-find must recover winding vectors in the chosen period basis.

CPU discovery stage:

- `N=29` magic-orientation pilot;
- same-N pairs `N=65,85,145,205,425`;
- exact/small-system regression tests where feasible;
- square-bond control before square-site interpretation.

GPU stage is justified only after a stable `cos(4 theta)` signal is seen. A GPU campaign should process many independent replicas for many orientations simultaneously.

## Relation to the modular-shape conjecture

This experiment keeps `tau=i` fixed and rotates the **microscopic lattice**. It therefore isolates lattice-spin anisotropy from modular shape dependence.

A later experiment can keep microscopic orientation fixed and vary `tau`; that tests the separate torus one-point/shape factor. Do not conflate the two effects.
