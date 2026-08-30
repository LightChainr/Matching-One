# P337 P1(F5) D4-even line-observer follow-up

The first line-addressed coupling-rank lift is resolved independently at N325
and N425.

The preceding frozen `P1(F3)` attempt had only one D4-even contrast; its
second row was reflection/projective odd and the determinant stayed null.
`P1(F5)` is the smallest prime projective line space with three D4-even
orbits:

```text
axes      = {0, infinity}
diagonal  = {+1, -1}
oblique   = {+2, -2}.
```

This gives two independent integer zero-sum rows,

```text
F5_X = axes - diagonal
F5_Y = axes + diagonal - 2 oblique.
```

At fixed pre-insertion occupancy, each row was coupled in the same batch to
`W_line` and `JS`.  The two complex determinants, one per orientation, were
tested jointly with their full aligned delete-one covariance.

## Frozen gate

| size | joint chi2(4) | gate | first det | second det |
|---|---:|---:|---:|---:|
| N325 | **149.93** | 13.277 | 0.30654-0.13469i | -0.26458-0.12725i |
| N425 | **246.93** | 13.277 | -0.20503+0.22254i | 0.33078+0.01526i |

The per-orientation two-component statistics are also individually resolved:
N325 gives 70.00 and 118.75; N425 gives 140.92 and 49.89.  No size pooling or
cross-size amplitude fit was used.  The 20k stop rule therefore fires and no
sample expansion is warranted.

The source-rank conclusion is narrower and stronger than a generic new-signal
claim: the physical continuous line phase `W_line` and the even birth source
`JS` have different response vectors across exact projective winding orbits.
The earlier `O_far/O_sep4` null is consequently observer blindness, not source
equality.  The F3 null is also explained: F3 lacks a second D4-even orbit
contrast, whereas F5 supplies one.

This is still a finite global-topology statement.  It does not name a CFT
field, establish an exponent, or make the F5 character an independent data
source from the ambient-H1 filtration.  Its next useful role is as a frozen
readout on another geometry family, not as another request for larger N325 or
N425 samples.

N325 ran on TgFr7R in 6.63 seconds and N425 on XPk2PZ in 8.89 seconds.  All
exact audits were zero.
