# Dimensionless amplitude ratios

Status: C0/C1 derivation for issue #118. This definition was integrated after the Issue #57 production reveal; it is not a numerical target for that block.

## Metric cancellation

Write the P48 leading laws as

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
R_T -> F_T'(0) / (F_T(0) M0'(0)).
```

The raw ratio `A_D/A_S` has net exponent `N^{-5/8}` and does not cancel two independent couplings. It is not a universality candidate.

`R_T` inherits the prospective failure of pure `P4[S'] ~ N^-5/4`; its finite-size drift can be used as another coordinate for the q=2/Jordan question. `R_I` uses only channels whose pure laws survived N=185/265, so it is the cleaner first diagnostic.

## Cross-model comparison is a separate claim

These ratios become universality candidates only after the lattice observable, wrapping channel, and torus modulus are matched across microscopic realizations. The derivation itself only identifies the metric-cancelling monomials.

## Descriptive reconstruction

Committed P48 scaled amplitudes plus the last P35 `B` row (N=170, `B≈1.74619`, common-metric approximation) give

```text
N=185:  R_I ~ 1.4218    R_T ~ 5.2645
N=265:  R_I ~ 1.2009    R_T ~ 5.1199
```

These numbers are not jackknife ratios inside a shared delete-one replicate, and `B` is not measured at N=185/265. They are development checks only.

## Fast next use

Recompute `Mbar'`, `P4[D]`, `P4[S']`, `P4[S]`, and `P4[D']` inside each delete-one replicate of existing full-curve blocks so the ratio covariance is exact. N=325/425 may be inspected retrospectively now that Issue #57 is revealed. If a held-out ratio score is useful, freeze it against an unrevealed block such as the N=145->290 full curve rather than retroactively labeling norm-5 as prospective.

## Oracle

```bash
python3 scripts/dimensionless_amplitude_ratios.py
```
