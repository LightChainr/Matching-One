# Modular-ring grading as a three-child correction spectrometer

The existing degree-2 hexagonal oracle gives the generator images

```text
E4 -> A(1,zeta,zeta^2),
E6 -> B(1,1,1).
```

Because the ordinary holomorphic ring is `C[E4,E6]`, this already determines every monomial correction on the same three geometries.

## 1. A C3 grading of the modular ring

For a monomial `E4^a E6^b`, the degree-2 child vector has character

```text
chi_child(E4^a E6^b)=zeta^(a mod 3).              (1)
```

Thus the Hecke child orbit defines a ring grading

```text
C[E4,E6] -> Z/3Z,
deg E4=1,
deg E6=0.
```

At the parent hexagonal point, every monomial with `a>0` vanishes to order `a`, while a pure `E6^b` survives.  The child character and parent zero order therefore give two independent discrete labels.

## 2. Modular derivatives do not add low-weight shapes

With `D_k=D-(k/12)E2`, the Ramanujan--Serre identities give

```text
D_4 E4=-E6/3,
D_6 E6=-E4^2/2,
6 D_6 D_4 E4=E4^2,
D_8(E4^2)=-2E4E6/3.                              (2)
```

Equation (2) is the useful syzygy: a modular-covariant “derivative correction” through weight 8 is not an extra fit direction.  `D4E4` is exactly the E6 shape, and `D6E6` or the second derivative of E4 is exactly the E4-squared shape.  Treating them as independent amplitudes would be an algebraic overparameterization.

An ordinary derivative `D` without the Serre connection is quasimodular.  If it produces child-vector leakage beyond the character prescribed by (1), that leakage is a positive quasimodular/Jordan diagnostic rather than another member of the ordinary ring.

## 3. Three children separate the first three shapes

For a complex response `y=(y0,y1,y2)` on

```text
(2omega, omega/2, (omega+1)/2),
```

the zero-parameter predictions are

```text
E4:       y/y0=(1,zeta,zeta^2), DFT support r=1,
E6:       y/y0=(1,1,1),         DFT support r=0,
E4^2:     y/y0=(1,zeta^2,zeta), DFT support r=2.  (3)
```

By (2), the same rows classify `D4E4` with E6 and `D6E6` with E4 squared.  More than one significant DFT component is a direct mixture/tangent/nonlocal-sector signal.

The imaginary channel is essential.  The real projections of the E4 and E4-squared rows are both

```text
1:-1/2:-1/2,
```

so real-only data aliases the two conjugate characters.  Their sine quadratures have opposite signs.

## 4. Missing correction powers

The ring contains no weight-2 modular form and no odd-weight modular forms.  Its low-weight spectrum is

```text
weight 4:  E4,
weight 6:  E6,
weight 8:  E4^2,
weight 10: E4E6,
weight 12: span(E4^3,E6^2).
```

For a dimensionless torus response with `x=weight`, these correspond to

```text
N^-1, N^-2, N^-3, N^-4, N^-5, ... .
```

Therefore an ordinary holomorphic vacuum correction cannot produce a half-integer power of `N`.  The shapes at `N^-1` through `N^-4` are one-dimensional; the first ordinary modular ambiguity occurs only at `N^-5`/weight 12.  Square and hexagonal values resolve it exactly:

```text
E4^3: visible at square, zero at hex,
E6^2: zero at square, visible at hex.
```

This does not forbid fractional corrections from thermal, charged, logarithmic, or nonlocal modules.  It says they cannot be relabeled as ordinary vacuum-ring descendants.

## 5. Frozen score

Keep one typed complex/chiral observable over the three degree-2 children and compute its three-point DFT.  Score the declared pure supports in (3), with full covariance.  No amplitude fit is required.

The parent hexagonal value is a secondary discriminator:

- nonzero parent plus constant child character selects an E6-containing direction;
- zero parent plus constant child character is compatible with `E4^3` rather than E6;
- the zero order distinguishes powers of E4 inside a shared child character.

## Boundary

The identities are exact for ordinary modular forms and Serre derivatives.  A lattice bridge requires a common complex frame and a declared chiral response.  Quasimodular, Jordan, thermal, charged, or topological corrections are allowed to violate the ordinary-ring grading; such a violation is the intended diagnostic.
