# Exact projection-leakage commutator control

This control implements the finite matrix warning in Issue 400. For a declared
orthogonal projection `P`, `Q=I-P`, and arbitrary maps `U,V`, exact
multiplication gives

```text
[PUP,PVP] = P[U,V]P + PVQUP - PUQVP.
```

The frozen rational witness uses commuting orthogonal involutions

```text
U = diag(1,-1,1),
V = diag(1,1,-1),
P = (1/3) [[2,1,-1],[1,2,1],[-1,1,2]].
```

Although `[U,V]=0`, the compressed commutator is

```text
(4/9) [[0,1,1],[-1,0,1],[-1,-1,0]],
```

with squared Frobenius norm `32/27`. Entrywise, it equals the leakage term
`PVQUP-PUQVP`; the intrinsic projected commutator is exactly zero. A diagonal
projection commuting with both maps provides an independent zero control.

Therefore a nonzero compressed rectangle can measure omitted-sector
propagation without witnessing microscopic noncommutation. Intermediate
Q-sector rows are required to subtract that ambiguity in production.

## Boundary

This exact finite-algebra certificate does not inspect P250 data, estimate
Q-sector excursions, claim microscopic memory or magnetic translations, alter
frozen rank exclusions, or authorize a new simulation. Issue 400 remains open.
