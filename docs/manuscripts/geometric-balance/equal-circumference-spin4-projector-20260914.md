# Equal-circumference angular projector for the sector-odd correction

Date: 2026-09-14

Status: analysis/design note motivated by `docs/research-compass-beyond-exactness-20260914.md`.  It uses only already committed safe-transfer outputs for the numerical compression below.  The proposed `(4,3)` computation is a discriminating next experiment, not an automatic width/angle scan.

## 1. Question after the oblique spin-four result

The current oblique safe-transfer evidence is already substantially stronger than an exponent fit:

- the critical charge-sector mismatch changes sign with `cos(4 theta)`;
- after physical metric normalization the free-energy amplitude collapses across axis, diagonal and several oblique directions;
- the charge-root shift shows the same angular sign/magnitude law.

The useful remaining question is therefore no longer simply

> is the leading correction compatible with spin four?

but rather

> after the leading spin-four piece is projected out, is the first residual still the same angular sector with an ordinary radial `ell^-2` dressing, or is there a genuinely angular-orthogonal lower competitor (scalar / spin eight / another map sector)?

This distinction is directly aligned with the research compass: identify what breaks the common sector, rather than extending a solved finite model for its own sake.

## 2. Angular Fourier bookkeeping

Let

```text
F(ell,theta) = ell^(17/4) E_charge^phys(pc;ell,theta),
```

where `E_charge^phys=|u| Theta_row` for primitive circumference direction `u` and `ell=n|u|`.

Square `C4` symmetry plus reflection allows the angular expansion

```text
F(ell,theta)
  = A0(ell)
  + A4(ell) cos(4 theta)
  + A8(ell) cos(8 theta)
  + A12(ell) cos(12 theta) + ... .
```

The current leading hypothesis is `A4(ell)->B != 0`, with the other harmonics subleading after the `ell^(17/4)` rescaling.

For axis and diagonal orientations,

```text
theta_axis = 0,
theta_diag = pi/4,
```

so

```text
F_odd  = [F(axis)-F(diag)]/2 = A4 + A12 + ...,
F_even = [F(axis)+F(diag)]/2 = A0 + A8 + A16 + ... .
```

Thus a scalar `x=21/4` competitor belongs to the angular-even projector, while the proposed spin-four term belongs to the angular-odd projector.

The same construction can be applied to the scaled root response

```text
R(ell,theta) = -(p_root-pc) ell^4.
```

## 3. Existing no-new-compute Pell projector

There is already a near-equal-circumference axis/diagonal pair:

```text
axis:     w=7, ell_a=7,
diagonal: n=5, ell_d=5 sqrt(2),
ell_d^2-ell_a^2 = 1.
```

The relative circumference mismatch is only about one percent and is forced by the Pell relation `7^2-2*5^2=-1`.

### Critical free-energy mismatch

From the committed axis and diagonal controls,

```text
F_axis(w=7)   = +1.0358025515351343,
F_diag(n=5)   = -1.0206765596904341.
```

Therefore

```text
F_odd  = 1.0282395556127844,
F_even = 0.0075629959223501,
|F_even/F_odd| = 0.0073553.
```

So the angular-even leakage is only about `0.74%` of the angular-odd component at this already-available pair.

Because the physical circumferences are not exactly equal, this number must **not** be called a scalar-amplitude measurement: ordinary radial finite-size drift can itself generate an `O(1%)` mismatch between the two members.  It is nevertheless a strong indication that a same-order scalar term is not competing at comparable amplitude.

### Root response

Using the same pair,

```text
R_axis(w=7) = +0.3035281536402238,
R_diag(n=5) = -0.2991554482270620,
```

hence

```text
R_odd  = 0.3013418009336429,
R_even = 0.0021863527065809,
|R_even/R_odd| = 0.0072554.
```

The free-energy and root projectors independently put the even leakage at essentially the same `~0.7%` level.

This is more informative than another fit of the exponent four: it says that the leading angular-orthogonal contamination is already small before any new production.

## 4. A stronger compression already visible: one spin-four tower with an `ell^-2` dressing

For three directions with two usable circumferences each, fit only the two-term form

```text
A_est(ell) = A_inf + A_2 / ell^2,
```

where `A_est=-(p_root-pc)ell^4/cos(4theta)`.

The resulting two-point coefficients are

| direction | `A_inf` | `A_2` |
|---|---:|---:|
| axis `(1,0)` from `w=8,9` | `0.2899430` | `0.656237` |
| diagonal `(1,1)` from `n=4,5` | `0.2863246` | `0.641543` |
| `(2,1)` from `n=3,4` | `0.2921425` | `0.690466` |

The agreement is unexpectedly tight for such small widths: the inferred asymptotic amplitude varies only at the percent level, and the `ell^-2` coefficient varies by only a few percent.

The direct critical free-energy amplitude shows the same pattern.  Writing

```text
B_est(ell)=B_inf+B_2/ell^2,
```

for the corresponding two-width pairs gives approximately

| direction | `B_inf` | `B_2` |
|---|---:|---:|
| axis `(1,0)` | `0.9774` | `2.830` |
| diagonal `(1,1)` | `0.9653` | `2.769` |
| `(2,1)` | `0.9847` | `2.965` |

The axis thermal-slope controls (`w=4..8`) are also well summarized by

```text
C(ell)=C_inf+C_2/ell^2,
C_inf ~= 3.377,
C_2   ~= 1.74,
```

where `C` is the physical `p` derivative coefficient in `partial_p E_charge^phys ~ C ell^-1/4`.

The independent leading amplitudes then satisfy

```text
B_inf/C_inf ~= 0.289,
```

matching the root-amplitude estimates above.

This motivates a sharper, falsifiable hypothesis:

> **Factorized spin-four tower.**  Over the first two visible correction orders, the critical sector-odd response is dominated by one `cos(4 theta)` angular sector with ordinary radial `ell^-2` dressing, while the leading thermal response is scalar.  The observed root ladder `ell^-4, ell^-6, ...` can therefore arise without assigning a new angular mechanism to every even power.

Schematically,

```text
E_charge^phys(pc;ell,theta)
 = cos(4theta) ell^(-17/4)
     [B0 + B2 ell^-2 + O(ell^-4)]
   + R_perp(ell,theta),

partial_p E_charge^phys(pc;ell,theta)
 = ell^(-1/4)
     [C0 + C2 ell^-2 + O(ell^-4)]
   + S_perp(ell,theta),
```

and therefore

```text
p_root-pc
 = -cos(4theta) ell^-4
     [A0 + A2 ell^-2 + O(ell^-4)]
   + Q_perp(ell,theta).
```

This is a data compression / conjecture, not an operator theorem.  The `ell^-2` dressing could arise from nonlinear scaling fields, descendant mixing, geometry, or another correction with the same angular character.  It should not yet be named as a specific Virasoro descendant.

## 5. One decisive next computation: exact equal circumference

The cleanest next experiment is **not** another angle ladder.  Use a primitive Pythagorean direction

```text
u=(4,3), |u|=5,
n=2,
ell=10,
cos(4theta)= -527/625 = -0.8432,
```

and compare it with the ordinary axis transfer at

```text
u=(1,0), w=10, ell=10.
```

These two cylinders have **exactly the same physical circumference** in the same microscopic square-site model.  Radial finite-size corrections therefore cancel from the leading angular ratio instead of being estimated or interpolated.

Compute only three objects at the fixed diagnostic `pc_ref` and at the charge root:

1. physical critical mismatch

```text
E = |u| [I4(pc)-I8(1-pc)];
```

2. physical thermal derivative

```text
D = |u| partial_p[I4(p)-I8(1-p)] at pc,
```

preferably by left/right Perron Feynman--Hellmann differentiation rather than a coarse finite difference;

3. charge root `p_root`.

The factorized hypothesis predicts the same-length ratios

```text
E_(4,3) / E_axis       ~= -527/625,
D_(4,3) / D_axis       ~= 1,
(p_root_(4,3)-pc)/(p_root_axis-pc) ~= -527/625,
```

up to genuinely angular-orthogonal residuals and numerical error.

The most useful reported quantities are therefore the **residuals**, not another fitted amplitude:

```text
r_E = E_(4,3) - (-527/625) E_axis,
r_D = D_(4,3) - D_axis,
r_p = (p_root_(4,3)-pc)
      - (-527/625)(p_root_axis-pc).
```

Normalize them by the corresponding leading axis signal and propagate Perron/eigensolver error into the differences.

### Diagnostic numerical scale, not an acceptance target

The existing two-width compression gives roughly

```text
A(ell=10) ~ 0.296,
B(ell=10) ~ 1.00,
C(ell=10) ~ 3.39.
```

So one expects root shifts of order

```text
axis ell=10:      p_root-pc ~ -2.96e-5,
(4,3), n=2:       p_root-pc ~ +2.50e-5,
```

if the factorized description is approximately right.  These estimates are only sanity checks; the test is the same-length angular ratio, not agreement with the extrapolated digits.

## 6. Compute boundary and stop rule

The generic oblique safe-transfer code already handles arbitrary primitive directions.  For `(4,3)` the matching row memory is larger than the diagonal case, but `n=2` is deliberately chosen so the first task is bounded.

Execution order:

1. build `(4,3), n=2` state spaces and report exact state counts / memory before any extension;
2. reproduce a previously committed oblique direction with the same executable environment;
3. compute `E,D,p_root` and certified/residual numerical errors;
4. compare only to axis `w=10` at the same physical circumference;
5. **stop** after this pair unless the residual is large enough to require identification.

Decision:

- if `r_E,r_D,r_p` are small at the percent level or below, do **not** open an angle ladder; the next high-value work is the module/matrix-element origin of the already isolated spin-four sector;
- if a residual is clearly larger than solver error and the known radial correction scale, identify its angular character before adding widths;
- only if a real angular-orthogonal residual is found should a near-node direction such as `(5,2)` receive a second width.

This experiment directly tests two explanations and therefore meets the research-compass criterion for additional compute.

## 7. Claim boundary

The existing numerical projector and two-width coefficient compression are algebraic re-analyses of committed deterministic safe-transfer outputs.  They are not asymptotic fits with enough widths to establish correction exponents.

The factorized spin-four tower is a new working conjecture.  The equal-circumference `(4,3)` design is proposed because it can falsify that conjecture while controlling the principal radial/metric confounder in one computation.