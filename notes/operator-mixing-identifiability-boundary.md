# Operator-mixing identifiability boundary

Status: exact algebraic slice of Issue 125.

For one scaling field of matching parity `eta`, the involution gives

```text
S^(n) allowed iff (-1)^n = eta,
D^(n) allowed iff (-1)^n = -eta.
```

Together with `L^(2-x+n y_t)` and `y_t=3/4`, the two proposed fields imply the
four familiar powers exactly:

```text
identity-even x=4:   P4[S]  ~ N^-1,    P4[D'] ~ N^-5/8,
thermal-odd x=21/4: P4[D]  ~ N^-13/8, P4[S'] ~ N^-5/4.
```

The important identifiability point is negative. A field contributes an analytic scaling function
`F(z)=f0+f1 z+...`; parity selects which combination sees `f0` and which sees `f1`, but does not
relate them. The primary four-channel map therefore has four independent columns

```text
[I:f0, I:f1, T:f0, T:f1]
```

and exact rank four. Calling this a two-amplitude model silently assumes two additional dynamical
relations, namely fixed `f_I1/f_I0` and `f_T1/f_T0`. Those relations must come from a continuum
calculation or be frozen from training data before a joint held-out prediction.

This audit does not dispute the zero pattern or powers. It prevents the symmetry selection rule from
being overinterpreted as an amplitude closure relation.
