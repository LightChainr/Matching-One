# Magic-torus conjecture: target the hexagonal modular point first

Status: deliberately speculative. This is a geometry experiment motivated by conformal spin and modular covariance, not a claimed formula for the percolation correction amplitude.

## Distinguish two controls

There are two different geometric operations in the research program.

### Microscopic orientation at fixed torus shape

The axis torus with periods

\[
(L,0),(0,L)
\]

and the diamond torus with periods

\[
(L,L),(-L,L)
\]

both have square modular shape `tau=i`; the second simply rotates the microscopic square lattice by `pi/4` relative to the torus cycles (and rescales the physical period by `sqrt(2)`).

This is the clean test for a spin-4 sign flip.

### Torus modular shape at fixed microscopic lattice

Choose two integer period vectors

\[
v_1=(m,0),\qquad v_2=(a,b),
\]

so that

\[
\tau=(a+ib)/m.
\]

Changing `(a,b,m)` varies the aspect ratio and shear/twist of the torus without rotating the microscopic square lattice itself.

Percolation wrapping/excess-cluster quantities are known to depend universally on aspect ratio and periodic twist, so this is a physically meaningful control parameter rather than a coordinate trick.

## CFT motivation

Square-lattice rotational symmetry permits continuum operators whose spin is a multiple of four. In other two-dimensional lattice CFTs, and in explicit finite-size studies, the leading lattice-anisotropy correction is often associated with

\[
T^2+\bar T^2,
\]

which has conformal spin `+4` and `-4` components and RG exponent `y=-2`.

A rotation by `pi/4` multiplies a pure spin-4 component by

\[
e^{i4\pi/4}=-1,
\]

which is precisely the sign structure sought in the axis/diamond experiment.

For percolation specifically, Feng, Deng and Blöte found a strong orientation dependence and sign change in the non-logarithmic `X_t2=4` correction amplitude on the square lattice. This does not by itself identify the field uniquely, but it makes the spin-4 mechanism plausible.

## Aggressive modular conjecture

Suppose the leading anisotropic finite-size correction of a dimensionless torus observable can be schematically factorized as

\[
\delta Q_4(L,\tau)
= g_4\,L^{-2}\,F_4(\tau,\bar\tau)+\cdots,
\]

where `g_4` is the microscopic square-lattice coupling to the spin-4 irrelevant field and `F_4` is a torus shape factor with weight-4 modular covariance (possibly nonholomorphic and mixed with logarithmic descendants in the `c=0` theory).

The simplest holomorphic weight-4 modular form for the full modular group is the Eisenstein series `E_4(tau)`. It has a simple zero at the hexagonal/equianharmonic elliptic point

\[
\tau_\hex = e^{i\pi/3}
=\tfrac12+i\tfrac{\sqrt3}{2},
\]

up to the conventional modular-equivalent choice `e^{2pi i/3}`.

This motivates the **targeted conjecture**:

> A torus whose modular parameter is near `tau_hex` may strongly suppress the leading spin-4 finite-size amplitude of square-lattice percolation observables.

This is not the statement `F_4=E_4`. In percolation the continuum theory is logarithmic (`c=0`), the observable is not the partition function, and scalar/logarithmic fields can mix with descendants. The value of the conjecture is that it gives a precise first geometry to test before a blind two-dimensional shape scan.

## Integer square-lattice approximants

A square lattice with periodic integer vectors can approximate `tau_hex` arbitrarily well.

Use

\[
v_1=(m,0),\qquad
v_2=(m/2,b)
\]

for even `m`, with

\[
b\approx \frac{\sqrt3}{2}m.
\]

Examples:

| `m` | `a=m/2` | `b` | `tau` imaginary part | sites `N=mb` |
|---:|---:|---:|---:|---:|
| 30 | 15 | 26 | 0.8666666667 | 780 |
| 52 | 26 | 45 | 0.8653846154 | 2340 |
| 82 | 41 | 71 | 0.8658536585 | 5822 |
| 112 | 56 | 97 | 0.8660714286 | 10864 |

The exact Diophantine approximants should be generated systematically rather than selected by hand; these examples only show that useful sizes are modest.

## Experiment M1: exact controls

Before square-site production runs:

1. implement a general integer-period quotient `(v1,v2)` with displacement-potential homology tracking;
2. verify the finite matching identity by exhaustive enumeration on tiny sheared quotients;
3. test square **bond** percolation at exact `p_c=1/2`;
4. compare the signed leading correction amplitude at:
   - square shape `tau=i`,
   - several paths approaching `tau_hex`,
   - nearby shapes on both sides of the candidate zero.

A real zero should manifest as a sign change, not merely as one accidentally small fitted coefficient.

## Experiment M2: square-site matching root

For each modular shape, measure the same matching observable and fit a joint model

\[
p_L^*(\tau)=p_c+A_4(\tau)L^{-w_4}+A_{next}(\tau)L^{-w_{next}}+\cdots
\]

with shared `p_c` and shape-dependent amplitudes.

Do not allow independent intercepts for each shape in the decisive fit. The point is to test whether a *single* threshold is approached with an amplitude that crosses zero in `tau`.

## Experiment M3: combine modular zero with orientation projection

The most aggressive version uses both knobs:

- microscopic axis versus `pi/4` orientation at fixed modular shape;
- modular shape near `tau=i` versus `tau_hex`.

If the leading correction is genuinely a spin-4 lattice field, these operations constrain its transformation in different ways. A data cube over

\[
(\text{matching parity})\times
(\text{orientation parity})\times
(\tau\text{ shape})
\]

could separate anisotropic, scalar, and logarithmic sectors much more cleanly than one long width sequence.

## Falsification

Demote this conjecture if any of the following occurs:

1. no stable orientation-odd correction is visible in the exact square-bond control;
2. the leading amplitude does not change systematically with torus shape;
3. no minimum/sign change appears near the hexagonal modular point or any nearby shape after physical-length normalization;
4. apparent cancellation moves substantially with the fitted size window;
5. the effect cannot be reproduced in at least two observables.

A negative result is still useful: it would separate microscopic orientation effects from torus modular-shape effects and constrain the operator interpretation.

## Why this is computationally attractive

The experiment does not require a transfer-matrix frontier. Moderate tori with thousands to millions of sites are straightforward for Monte Carlo, and many modular shapes can be batched.

It is therefore a good CPU discovery problem and an excellent later GPU many-replica problem. The expensive part is statistical discrimination of a small correction amplitude, not memory growth of a connectivity-state basis.

## References / conceptual anchors

- Feng, Deng, Blöte (2008): square-lattice percolation correction amplitude depends strongly on orientation; `X_t2=4` correction and logarithmic companion.
- Ziff, Lorenz, Kleban (1999): critical percolation excess-cluster/cross-configuration quantities depend on torus aspect ratio and periodic twist.
- Standard modular-form fact: the normalized weight-4 Eisenstein series `E_4` vanishes at the hexagonal/equianharmonic elliptic point.

Again: the final bullet is a **heuristic analogy**, not an identification of the percolation correction with `E_4`.
