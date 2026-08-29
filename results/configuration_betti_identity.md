# Configuration Euler–Poincaré / Betti identity

Source: `scripts/configuration_betti_identity.py`.
Claim level: C5 finite identity. Equivalent to P34; not a new control variate.

On every enumerated square torus configuration

```text
chi = V - E + F0 = beta0_black - beta0_white - q
```

The wrapping-difference variable `q` is the common event difference in
`{-1,0,+1}` already used by P34. It is not the homology-rank difference
`r_black - r_white`. Cyclomatic numbers bound wrapping rank from above:

```text
kappa_black = E_primal - V_black + beta0_black >= r_black
kappa_white = E_matching_white - V_white + beta0_white >= r_white
```

The issue-#111 variance-reduction branch therefore closes: every proposed
Betti control is an algebraic rewrite of `(V, E, F0, q, C_black, C_white)`.

## Exhaustive tiny quotients

| geometry | N | configs | identity fail | cyclo fail | `q = r_b-r_w` | empty `(q,r_b,r_w)` |
|---|---:|---:|---:|---:|---:|---|
| axis | 4 | 16 | 0 | 0 | 4/16 | `(-1,0,2)` |
| axis | 9 | 512 | 0 | 0 | 162/512 | `(-1,0,2)` |
| gaussian-2-1 | 5 | 32 | 0 | 0 | 10/32 | `(-1,0,2)` |
| diamond | 8 | 256 | 0 | 0 | 68/256 | `(-1,0,2)` |

Empty-mask counterexample (every listed quotient): `q=-1`, `r_black=0`,
`r_white=2`, `beta0_black=0`, `beta0_white=1`. The Euler identity still
holds because `chi = 0 = 0 - 1 - (-1)`.

## Boundary

This does not introduce a production Newman–Ziff Betti statistic and does
not claim a variance reduction relative to P34 motif controls. Averaging
the configuration identity recovers the expected Mertens–Ziff relation
already used in P34.
