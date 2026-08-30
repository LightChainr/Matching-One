# Projective essential-birth Phase-A certificate

Overall exact gates: **PASS**.

| geometry | N | paths | direct 0->2 | ell adds direction beyond (K1,K2) |
|---|---:|---:|---:|---|
| axis-L2 | 4 | 24 | 8 | True |
| gaussian-2-1 | 5 | 120 | 0 | True |
| c4-self-matching-3-1 | 10 | 3628800 | 518400 | True |

The primitive line is canonical up to sign and covariant under the tested unimodular basis changes.  On the axis quotient, D4 rotations leave `chi4` fixed and reflections conjugate it.  Complement/Alexander duality maps `(K1,K2,ell,site1,site2)` to `(N+1-K2,N+1-K1,ell,site2,site1)` exactly.

Integral saturation fixes `iota=1`: after the first line-bearing birth the subgroup is `Z ell`, and after the second birth it is `Z^2`.

The exact controls sharpen the word *independent*.  The line is not determined by `(K1,K2)`, yet conditional on a line-bearing birth it factorizes exactly from `(K1,K2)` in all three controls.  Moreover the two supported lines are related by a quarter-turn and have the same `chi4`, so these tiny quotients show a projective direction degree of freedom but no additional spin-4 value.  A larger quotient is needed to test directional spin-4 bias.

There is also an exact boundary to the proposed mark: axis-L2 and the C4 control contain direct `0->2` births.  Such paths have no canonical projective line.  Consequently the full `K1` histogram is recovered by summing over projective lines **plus** the typed `DIRECT_RANK2` atom; summing over `ell` alone recovers only line-bearing births.

## Crosswalk to Issue 156

Configuration by configuration, the plateau character

```text
A4(p)=E[1{tau1<=p<tau2} chi4(ell1)]
```

is exactly the fixed-p rank-one primitive-sector character.  It is not a new observable.  The genuinely new information is the source/sink split

```text
dA4/dp = j4_birth1(p) - j4_exit2(p),
```

where the exit at the second birth retains the line born at the first.  The certificate verifies this identity coefficient by coefficient in the degree-`N-1` Bernstein basis.  Saving only `A4` repeats Issue 156; saving both births localizes timing versus direction.
