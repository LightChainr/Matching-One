# Angular-irrep improvement of the matching charge root

Date: 2026-09-14

Status: finite deterministic / exploratory transfer analysis plus a new finite-size-improvement conjecture.  The transfer semantics are the committed safe lifted-gain semantics; local reimplementations were checked against committed roots to about `1e-11` or better.  No new infinite-volume theorem is claimed.

## 1. Motivation

The oblique safe-transfer programme shows that the dominant finite-width matching-root displacement on the square lattice transforms as the real spin-four harmonic

```text
H4(theta)=cos(4 theta).
```

The usual presentation is

```text
p_ch(ell,theta)-pc
 = -A4 H4(theta) ell^-4 + ... .
```

A stronger use of this fact is possible.  Different orientations at the **same physical circumference** can be combined as an exact D4-irrep projector which removes the entire observed spin-four component without knowing `pc` or `A4` in advance.

This turns angular dependence from a field-identification diagnostic into a finite-size-improvement mechanism.

## 2. Two-orientation projector

Suppose two short-period directions have the same physical circumference `ell`, with spin-four harmonics `h1!=h2` and charge roots `p1,p2`.

At leading order

```text
p_i = pc - A4 h_i ell^-4 + higher.
```

Then

```text
A4_hat = ell^4 (p2-p1)/(h1-h2),

pc_hat^(4-perp)
 = (h1 p2 - h2 p1)/(h1-h2).
```

Neither formula uses an external threshold value.  The second is the component of the two-root vector orthogonal to the `H4` irrep.

More generally, write a D4 harmonic expansion

```text
p_ch(ell,theta)
 = pc + S0(ell)
      + H4(theta) S4(ell)
      + H8(theta) S8(ell)
      + ... .
```

The two-angle projector removes **all** terms proportional to `H4(theta)`, not only the first `ell^-4` coefficient.  If the axial `ell^-6` correction belongs to the same spin-four tower, it is removed as well.  The remainder directly probes scalar and higher-irrep contamination.

This is the finite-size analogue of a Symanzik/improved-observable projection.

## 3. Exact equal-circumference pair at ell=5

The primitive direction

```text
u=(3,4)
```

has Euclidean length `5`.  Thus

```text
axis:       (1,0), n=5, ell=5, H4=1,
oblique:    (3,4), n=1, ell=5, H4=-527/625.
```

The `n=1` oblique quotient is honest: the primitive period has length five, larger than the local NN/matching interaction range, and the transformed edges have strictly positive longitudinal displacement.  The lifted-gain transfer has only

```text
G4 states = 45,
G8 states = 147.
```

Roots:

```text
p_axis5 = 0.5922358232050258,
p_34,1  = 0.5931582013380546.
```

The spin-four projector gives

```text
pc_hat_ell5 = 0.5927362453692130,
A4_hat_ell5 = 0.3127638526162537.
```

For comparison only after the fact, the diagnostic high-precision threshold reference `0.59274605079` differs by

```text
-9.8054e-6.
```

The individual raw roots are displaced by roughly `5e-4`; the irrep projection has already removed most of that finite-width error.

## 4. Exact equal-circumference pair at ell=sqrt(65)

The integer circle `a^2+b^2=65` contains two inequivalent primitive directions

```text
u1=(1,8),
u2=(4,7),
```

with

```text
H4(u1)= 3713/4225 = 0.8788165680...,
H4(u2)=-2047/4225 =-0.4844970414....
```

Both use `n=1`, hence

```text
ell=sqrt(65)=8.062257748... .
```

Safe-transfer state counts remain small:

```text
(1,8): G4/G8 = 1394 / 2153,
(4,7): G4/G8 = 1134 / 6216.
```

The roots are

```text
p_18 = 0.5926839084948207,
p_47 = 0.5927803979513810.
```

The exact harmonic projector is

```text
pc_hat_sqrt65
 = (3713 p_47 + 2047 p_18)/5760
 = 0.5927461073406902,

A4_hat_sqrt65
 = 0.2990272752626905.
```

Again the reference threshold is not used in the estimator.  Its after-the-fact difference is

```text
+5.66e-8.
```

This pair is especially useful because neither direction is an axis: the improvement is therefore not an axis-specific cancellation.

## 5. Exact equal-circumference pair at ell=10

Use

```text
axis:       (1,0), n=10, ell=10, H4=1,
oblique:    (3,4), n=2,  ell=10, H4=-527/625.
```

The oblique transfer has

```text
G4/G8 states = 7259 / 75541,
row memory    = 4 / 7,
p_34,2        = 0.5927709164472407.
```

A dedicated C++ realization of the committed axial lifted-gain algorithm was checked at `w=9` against the repository root to below `1e-12`, and gives

```text
w=10 safe states G4/G8 = 8953 / 8953,
p_axis10               = 0.5927163956300879.
```

Therefore

```text
pc_hat_ell10
 = (625 p_34,2 + 527 p_axis10)/1152
 = 0.5927459750664773,

A4_hat_ell10
 = 0.2957943638928858.
```

After the fit, comparison to the diagnostic reference gives

```text
-7.57e-8.
```

The ell=5 and ell=10 projected errors differ in magnitude by about a factor `129`.  With only two related scales this must **not** be promoted to an exponent fit.  It does show that the algebraic angular cancellation removes far more than the raw `ell^-4` error.

## 6. Strong constraint on a leading scalar x=21/4 explanation

A scalar correction of the same root exponent four would enter `S0(ell)`, not `H4 S4`, and would survive the projector.

The projected residuals at `ell≈8--10` are of order `1e-7`, while the observed spin-four raw correction at the same sizes is of order `3e-5`--`7e-5` in `p`.  Therefore an orientation-independent exponent-four component large enough to explain the original root displacement is incompatible with these angular projections.

This does not prove that every scalar `x=21/4` amplitude is exactly zero.  It says that such a scalar cannot be the dominant mechanism behind the observed matching-root `ell^-4` shift.

## 7. Thermal denominator is also orientation-scalar at equal ell

Define the physical thermal-slope indicator

```text
D(theta,ell)
 = |u| d_p Theta_row(p_root) ell^(1/4)
 = Q_p(p_root)/ell^(3/4),

Q = n |u|^2 Theta_row.
```

For the exact `ell=10` pair,

```text
D_axis10 = 3.39179932,
D_34,2   = 3.39187874,
relative difference = 2.34e-5.
```

For the equal `ell=sqrt(65)` pair,

```text
D_18 = 3.40273574,
D_47 = 3.40296413,
```

again differing only at the `1e-4` level.

Thus the root angular law is cleanly factorized at these sizes into

```text
spin-four critical numerator / scalar leading thermal denominator.
```

## 8. Root/slope-normalized shape: the TANGENT_SPIN4 diagnostic

For each orientation define

```text
y = Q_p(p_root)(p-p_root),
Psi(y)=Q(p).
```

If the leading spin-four correction acts only by translating the thermal scaling coordinate, all orientation dependence should disappear from `Psi` to first order.

At `ell=10`, axis versus `(3,4)` gives

```text
             axis10             (3,4),n2          difference

y=-0.50   -0.4995465515       -0.4995187757       +2.78e-5
y=-0.25   -0.2495311293       -0.2495230387       +8.09e-6
y=+0.25   +0.2512046459       +0.2512132325       +8.59e-6
y=+0.50   +0.5063352579       +0.5063680241       +3.28e-5
```

Local normalized coefficients

```text
Psi(y)=y+a2 y^2+a3 y^3+...
```

are

```text
axis10:     a2≈0.0133250, a3≈0.0235507,
(3,4),n2:   a2≈0.0134625, a3≈0.0235648.
```

The equal `sqrt(65)` pair behaves similarly.  Its slope-normalized curves differ by only `~1e-5--5e-5` on `|y|<=0.5`, despite roots lying on opposite sides of the threshold.

This is strong finite-width evidence for

```text
Q_{ell,theta}(X)
 = F(X + delta X_4(theta,ell))
   + smaller normal correction,
```

rather than an orientation-dependent deformation of the full scaling-function shape at leading order.

## 9. Angular-improvement hierarchy conjecture

The natural stronger conjecture is

```text
p_ch(ell,theta)-pc
 = H4(theta) [a4 ell^-4+a6 ell^-6+...]
   + H8(theta) [b8 ell^-8+...]
   + S0(ell)
   + ... .
```

If the known axial `4,6,...` correction ladder initially belongs to one spin-four family, the two-orientation projector removes the whole bracket, not just `a4 ell^-4`.  The first surviving term could then be an `H8` channel or a genuinely scalar dual-odd correction.

For the `(3,4)` angle,

```text
H8=cos(8 theta)=164833/390625,
```

and the H4-cancelled axis/(3,4) combination carries effective

```text
H8_eff = 429/625 = 0.6864.
```

The current ell=5 and 10 residuals are not sufficient to decide whether the next asymptotic term is scalar `ell^-6`, spin-eight `ell^-8`, a logarithmic collision term, or a mixture.  Their main use is to demonstrate that the leading spin-four tower can be removed by symmetry before fitting the next mechanism.

## 10. Practical consequence

The angular projector provides two outputs from the same small transfer calculations:

1. a mechanism test for the D4 representation content of the root correction;
2. an improved threshold estimator which can converge much faster than either orientation root separately.

This suggests a new order of operations for future high-precision work:

```text
first project out known lattice irreps,
then fit / identify the residual correction.
```

Do not spend larger widths to fit a correction that an exact symmetry projector can remove at smaller widths.

## 11. Claim boundary

- Integer directions, physical circumferences, D4 harmonics, safe-state counts and finite roots listed above are deterministic finite-transfer quantities.
- `pc_hat^(4-perp)` is a reference-free finite-size estimator conditional only on using the observed D4 harmonic as the projected component; comparison with the external threshold is diagnostic after the fact.
- The statement that the entire `ell^-6` axial correction belongs to the same spin-four tower is a conjecture.
- The apparent superconvergence of the projected estimator is not yet assigned an exponent.
- The functional `TANGENT_SPIN4` interpretation remains a scaling conjecture, albeit now supported by equal-circumference full-curve controls rather than root exponents alone.
