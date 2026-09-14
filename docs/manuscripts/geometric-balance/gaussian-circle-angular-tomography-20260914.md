# Gaussian-circle angular tomography for square-lattice finite-size corrections

Date: 2026-09-14

Status: exact angular / arithmetic algebra plus a finite-size-scaling programme.  The linear-algebra statements below do not identify any continuum field and do not assume the matching-root correction is exhausted by finitely many harmonics.

## 1. Fixed-norm periods remove size, modulus and Smith confounds

Let

```text
u=(a,b),
gcd(a,b)=1,
N=a^2+b^2.
```

Use the square torus period basis

```text
u=(a,b),
nu_perp=(-b,a).
```

Then

```text
det(nu,nu_perp)=N,
```

and the physical modulus is exactly square.  Because `gcd(a,b)=1`, the Smith invariants of the 2x2 period matrix are

```text
(1,N).
```

Therefore two primitive representations of the same integer `N` as a sum of two squares have simultaneously

```text
same site count,
same physical circumference sqrt(N),
same continuum square modulus,
same cyclic Smith class,
```

while differing only in the embedding angle of the microscopic square lattice relative to the period cycle.

This is the natural setting for angular finite-size tomography.

## 2. D4/reflection symmetry reduces the angular basis to Chebyshev polynomials

For a reflection-even square-lattice scalar quantity, orientation dependence is periodic under `theta -> theta+pi/2` and even under `theta ->-theta`.  Hence its angular Fourier expansion has the form

```text
y_N(theta)
 = A_0(N)
 + sum_{m>=1} A_{4m}(N) cos(4m theta).
```

Put

```text
h=cos(4 theta).
```

Then exactly

```text
cos(4m theta)=T_m(h),
```

where `T_m` is the Chebyshev polynomial of the first kind.  Thus at fixed `N`, angular tomography is ordinary polynomial interpolation in the scalar coordinate `h`:

```text
y_N(h)
 = A_0+A_4 T_1(h)+A_8 T_2(h)+A_12 T_3(h)+... .
```

No continuum-field assumption enters this change of basis.

## 3. Exact invertibility theorem

Suppose the same Gaussian circle has `k` inequivalent primitive orientations with distinct

```text
h_i=cos(4 theta_i), i=1,...,k.
```

Consider the `k x k` matrix

```text
V_{i,m}=T_m(h_i),  m=0,...,k-1.
```

Because `T_m` has degree `m`, with leading coefficient

```text
1, 1, 2, 4, 8, ...,
```

the determinant is a nonzero constant times the ordinary Vandermonde determinant:

```text
det V
 = 2^((k-1)(k-2)/2)
   product_{i<j}(h_j-h_i).
```

Hence:

> **Gaussian-circle tomography theorem.**  Distinct `h_i` imply that the first `k` square/reflection-even angular irreps
>
> ```text
> H0,H4,...,H_{4(k-1)}
> ```
>
> are exactly identifiable from `k` same-circle measurements, with no cross-size, modulus or Smith transfer law.

This is purely algebraic.

## 4. Closure residuals are higher-irrep detectors

With `k` orientations, fit only the first `r<k` angular columns.  The left nullspace of the resulting `k x r` design provides exact contrasts annihilating those irreps.

For example, four N1105 orientations and the `H0/H4/H8` model give the integer null contrast

```text
-121 y_(4,33)
+342 y_(9,32)
-247 y_(12,31)
+ 26 y_(23,24).
```

It annihilates `H0,H4,H8` exactly.  Its unit-norm gain on `H12` is about `-0.705`, so it is a genuine higher-harmonic closure test rather than a nearly blind residual direction.

Conversely, using all four columns `H0,H4,H8,H12` gives an invertible design with condition number about `4.27` for N1105.

## 5. Exact irrep-improved estimators

Let `p_i` be orientation-dependent pseudo-critical roots at one fixed Gaussian circle.  If one wants the angular-scalar component after removing the first `k-1` nontrivial harmonics, solve weights `w_i` satisfying

```text
sum_i w_i = 1,
sum_i w_i H4(theta_i)=0,
...
sum_i w_i H_{4(k-1)}(theta_i)=0.
```

Then

```text
p_scalar^(<=4(k-1) removed)
 = sum_i w_i p_i
```

is exactly insensitive to those angular irreps at that finite `N`.

For N1105, annihilating `H4,H8,H12` gives approximately

```text
w = (
 +0.4149325453,
 -0.4015794286,
 +0.7330946817,
 +0.2535522016
).
```

The remaining effective `H16` gain is about `0.9041`; thus the projector is not accidentally blind to all higher harmonics.

This is an angular analogue of Symanzik improvement: known lattice irreps are removed algebraically before any residual size extrapolation.

## 6. Two-angle H4 projector and its residual H8 geometry

For two orientations with

```text
h_i=H4(theta_i),
```

the H4-annihilating scalar projector leaves an H8 coefficient

```text
C8
 = [h1 H8(theta2)-h2 H8(theta1)]/(h1-h2)
 = -(1+2 h1 h2).
```

Therefore the residual sensitivity to H8 can itself be designed arithmetically.

Examples:

```text
N=65:  (1,8) vs (4,7),   C8=-0.1484319457
N=85:  (2,9) vs (6,7),   C8=+0.2224938303
N=377: (4,19) vs (11,16), C8=+0.0036146395
```

The N377 pair is thus an `H4/H8` double-notch design: after exact H4 cancellation it retains only about `0.36%` of the H8 coupling.  Its projected gains on later harmonics are not simultaneously tiny (`C12≈-0.1378`, `C16≈-0.9811`), so it can distinguish an H8 residual from scalar/H12/higher structure.

A Python state-cap pilot shows the N377 width-one safe automata exceed 200k states in both graph families; this is therefore a targeted CPU candidate, not a default small calculation.

## 7. Why the N85 threshold coincidence should not be overinterpreted

The H4-projected root estimates currently include

```text
N=65  : pc_hat - pc_ref ≈ +5.66e-8,
N=85  : pc_hat - pc_ref ≈ +4.33e-10,
ell=10: pc_hat - pc_ref ≈ -7.57e-8.
```

The N85 value is striking, but these three pairs have different residual H8 weights `C8`.  Treating their physical circumferences `8.06--10` only as a narrow exploratory band, the three residuals are nearly linear in `C8`:

```text
r ≈ S0 + C8 S8,
S0 ≈ +3.4e-8,
S8 ≈ -1.59e-7,
```

with remaining discrepancies of order `1e-9` in this crude local linearization.

Thus a plausible explanation of the exceptionally accurate N85 estimator is **accidental cancellation between a small scalar residual and an H8 residual**, not a proven asymptotic error of `4e-10`.

This is precisely why same-circle `H0/H4/H8` tomography (N1105) is more informative than celebrating one unusually accurate two-angle estimate.

## 8. Number-theoretic design principle

The useful resource is not a large site count by itself, but the number and geometry of primitive sum-of-two-squares representations of `N`.

- two representations permit exact H4 cancellation;
- three permit same-circle `H0/H4/H8` identification;
- four permit one held-out closure test for that three-column model or exact inclusion of H12.

N1105 is the first small circle in the current search with four useful primitive orientations and therefore occupies a qualitatively different role from ordinary larger-width production.

The design problem can be phrased before any simulation:

```text
choose N and primitive representations
 to optimize conditioning / target-irrep gain / state-space cost.
```

This converts angular acquisition from an ad hoc angle list into an exact arithmetic design problem.

## 9. Connection to Matching One mechanisms

For the matching charge root, the current finite evidence suggests:

```text
A4 dominates the raw finite-size displacement;
thermal slope is leading-angle-scalar;
root/slope-normalized charge shape is nearly angle-invariant.
```

Gaussian-circle tomography now asks the next correct question:

```text
after removing the dominant H4 regularization mismatch,
is the residual angular-scalar, H8, H12, or something outside the D4 harmonic truncation?
```

That question is more informative than another power-law fit of the raw axial root.

For other observables, the same algebra applies provided the observable is declared reflection-even and the compared tori genuinely share N, modulus and Smith class.

## 10. Claim boundary

- The Gaussian-circle/Smith statements, Chebyshev basis, determinant formula and exact null/projector constructions are algebraic facts.
- Finite roots quoted from the current transfer programme are deterministic controls but not rigorous threshold enclosures.
- The decomposition of the H4-projected residual into scalar/H8 pieces is exploratory until same-circle multi-angle data such as N1105 are computed.
- No angular coefficient by itself identifies a CFT/LCFT field.
- Arithmetic design does not authorize a production run; state-space cost and source-specific signal/noise remain separate decisions.
