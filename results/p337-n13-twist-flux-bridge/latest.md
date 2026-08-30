# Exact N13 bridge: F3 twist H to primitive A4 flux

All state, source/sink, odd-sector and quarter-turn transport gates pass.

## Fixed character bridge

For `P=[[3,-2],[2,3]]`, both axis lines have

`z_axis=-119/169 + (120/169) i`.

Both diagonal lines have `-z_axis`. Therefore, coefficient by coefficient,

```text
H_F3 = T_01+T_10-T_12-T_11
     = L_axis_x+L_axis_y-L_diag_plus-L_diag_minus,
A4   = z_axis H_F3.
```

The unit-norm convention in `a7cb19a` is `H_unit=H_F3/2`, so in that stored coordinate `A4=2 z_axis H_unit`. This factor is normalization, not a fitted amplitude.

## Source/sink derivative

At every lower size `k`, with degree-12 Bernstein normalization,

```text
dH_F3/dp = J_H,birth1 - J_H,exit2,
J4_birth1 = z_axis J_H,birth1,
J4_exit2  = z_axis J_H,exit2,
dA4/dp    = z_axis dH_F3/dp.
```

These identities pass independently at every coefficient; no evaluation point or path sampling is used.

## Odd sectors

Both `A=L_axis_x-L_axis_y` and `D=L_diag_plus-L_diag_minus` vanish coefficientwise in the rank-one state curve, birth source and exit sink. They are exact `S`-quarter-turn-odd symmetry zeros of the `3+2i` quotient, not projective sectors forbidden on general geometries.

## At p_ref

- `H_F3=0.299772544239`
- `A4=-0.211082442393+0.21285624443 i`
- `dH_F3/dp=0.125800093888`
- `dA4/dp=-0.088581131199+0.0893255104528 i`

This closes the exact state-to-current dictionary between #337 and #334. It does not identify a continuum field.
