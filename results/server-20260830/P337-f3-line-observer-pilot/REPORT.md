# P337 F3 projective-line observer pilot

The exact `P1(F3)` observer does not lift the `W_line/JS` coupling rank at
20k per size.

The two frozen rows were the existing finite-abelian line characters

```text
F3_H = +1 axes, -1 diagonals
F3_D = +1 diag+, -1 diag-, 0 axes.
```

They were evaluated at fixed pre-insertion `k0` and coupled in the same batch
to columns `W_line` and `JS`.  Each orientation supplied a complex 2x2
determinant; the four real determinant components were tested jointly within
size.

| size | joint chi2(4) | frozen gate | first det | second det |
|---|---:|---:|---:|---:|
| N325 | 2.269 | 13.277 | 0.01485-0.00667i | -0.02417-0.01140i |
| N425 | 2.312 | 13.277 | -0.01861+0.02029i | 0.00235+0.00011i |

The null has a useful structural diagnosis.  `F3_H` couples strongly to both
sources, whereas every `F3_D` coupling is close to zero.  The second row is a
reflection/projective-odd channel, while `W_line` and `JS` are the even
spin-four sources used here.  Thus this attempt did not fail because W_line
collapsed back into JS; it chose a symmetry-null second observer row.

The smallest non-adaptive continuation is `P1(F5)`.  Its axes, diagonals and
oblique lines form three D4-even orbits, hence two independent zero-sum even
characters.  Those two rows can test coupling rank without borrowing a
reflection-odd channel.  This continuation requires a separately frozen
counter block; the present F3 gate remains recorded as null.

N325 ran on TgFr7R in 6.44 seconds and N425 on XPk2PZ in 9.25 seconds.  All
exact audits were zero.  No field, exponent or continuum defect is assigned.
