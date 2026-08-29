# Dimensionless amplitude ratios

Status: C0/C1 derivation freeze for issue #118. Descriptive P48 reconstructions are development-only and are not a numerical target for issue #57.

## Metric cancellation

Write the frozen P48 leading laws as

```text
P4[S]    = A_I(0) F_I(0) N^{-1}            + ...
P4[D]    = u_T     F_T(0) N^{-13/8}        + ...
P4[D']   = u_I b   F_I'(0) N^{-5/8}        + ...
P4[S']   = u_T b   F_T'(0) N^{-5/4}        + ...
Mbar'    = b       M0'(0) N^{3/8}          + ...
```

Here `b` is the nonuniversal thermal metric `dt/dp` and `u_I`, `u_T` are lattice couplings to the two spin-4 sectors. The combinations

```text
R_I = P4[D'] / (P4[S] * Mbar') = A_Dp / (A_S * B)
R_T = P4[S'] / (P4[D] * Mbar') = A_Sp / (A_D * B)
```

have net exponent `N^0` and cancel `b`. At fixed torus modulus the putative limits are response ratios of the even/odd scaling functions,

```text
R_I -> F_I'(0) / (F_I(0) M0'(0))
R_T -> F_T'(0) / (F_T(0) M0'(0))
```

The raw ratio `A_D / A_S` has net exponent `N^{-5/8}` and does **not** cancel two independent couplings. It is not a universality candidate.

`R_T` inherits the prospective failure of pure `P4[S'] ~ N^{-5/4}`. Its finite-size drift is itself a discriminator between the frozen q=2 and Jordan-log corrections. `R_I` uses only channels whose pure laws survived N=185/265, so it is the cleaner first target.

## Cross-model test is not implied

These ratios become universal only after the lattice observable, wrapping channel, and torus modulus are matched across microscopic realizations (square-site/matching, C4 self-matching, square-bond, later isoradial controls). This note does not claim that agreement. It only freezes the metric-cancelling monomials.

A stronger factorization proposed on #118,

```text
A_H4(model, tau) = lambda4(model) * C_observable * g2(tau)
```

requires the #106 microscopic tensor and is out of scope here.

## Descriptive reconstruction

Committed P48 scaled amplitudes plus the last P35 `B` row (N=170, `B≈1.74619`, common-metric approximation) give

```text
N=185:  R_I ~ 1.4218    R_T ~ 5.2645
N=265:  R_I ~ 1.2009    R_T ~ 5.1199
```

These numbers are not jackknife ratios inside a shared delete-one replicate, and `B` is not measured at N=185/265. They exist only to check that the scaled formula reproduces the hand estimates already recorded on the issue. They must not be copied into an issue #57 prediction file.

Relative drift of `R_T` between these two sizes is smaller than that of the raw scaled `P4[S']` amplitude, which is the comparison recorded on the issue. This does not make `R_T` universal.

## Required next measurement

Recompute `Mbar'`, `P4[D]`, `P4[S']`, `P4[S]`, `P4[D']` inside each delete-one replicate of an existing full-curve block so the ratio covariance is exact. Freeze any parent/child `R_I` or `R_T` target for norm-5 *before* reading the N=325/425 children.

## Oracle

```bash
python3 scripts/dimensionless_amplitude_ratios.py
```
