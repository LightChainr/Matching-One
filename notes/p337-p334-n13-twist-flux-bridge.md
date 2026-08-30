# Exact bridge from F3 flat twists to projective A4 current

Status: exact N13 state and directed-boundary census on the Gaussian quotient
`3+2i`.  This connects the #337 constraint-sector coordinates to the #334
primitive state/source/sink coordinates without path enumeration or sampling.

## Normalization dictionary

Use one nonzero F3 twist representative for each projective kernel:

```text
T_01 -> axis x,      T_10 -> axis y,
T_12 -> diagonal +, T_11 -> diagonal -.
```

The common rank-zero term cancels from

```text
H_F3 = T_01+T_10-T_12-T_11
     = L_axis_x+L_axis_y-L_diag_plus-L_diag_minus.
```

For the period matrix

```text
P=[[3,-2],[2,3]],
```

the physical primitive fourth character is constant on each orbit:

```text
z_axis = (-119+120 i)/169,
z_diagonal = -z_axis.
```

It follows at every microcanonical subset size, hence as an exact polynomial
identity in `p`, that

```text
A4(p)=z_axis H_F3(p).
```

The #337/a7cb19a unit-vector convention used

```text
H_unit=(L_axis_x+L_axis_y-L_diag_plus-L_diag_minus)/2,
```

so its exact dictionary is `A4=2 z_axis H_unit`.  The factor two is entirely
the chosen projective-character normalization; `z_axis` is a fixed geometric
character, not a fitted coupling.

## The derivative is the same current in two coordinates

Let `h_k` denote the degree-13 Bernstein coefficient of `H_F3`.  At every
lower size `k=0,...,12`, divide directed edge counts by `C(12,k)` and define

```text
J_H,birth(k) = axes birth edges - diagonal birth edges,
J_H,exit(k)  = axes exit edges  - diagonal exit edges.
```

The exact state-boundary identity is

```text
13 [h_(k+1)-h_k] = J_H,birth(k)-J_H,exit(k).
```

Line by line the physical character then gives

```text
J4_birth(k) = z_axis J_H,birth(k),
J4_exit(k)  = z_axis J_H,exit(k),
dA4/dp      = z_axis [J_H,birth-J_H,exit].
```

All four equalities pass separately for all 13 state coefficients and all 13
degree-12 derivative coefficients.  Thus the F3 twist derivative and the #334
birth1/exit2 current are not two mechanisms on this quotient; they are the
same exact current written in constraint and physical-character bases.

At `p_ref` the exact values are

```text
H_F3       = 0.299772544239,
A4         = -0.211082442393 + 0.212856244430 i,
dH_F3/dp   = 0.125800093888,
dA4/dp     = -0.088581131199 + 0.089325510453 i.
```

These numbers only illustrate the polynomial identity; the proof is
coefficientwise and does not depend on `p_ref`.

## The odd coordinates are exact nulls here

The complete state census shows

```text
A = L_axis_x-L_axis_y = 0,
D = L_diag_plus-L_diag_minus = 0
```

at every subset size.  More strongly, their birth-source and exit-sink edge
differences vanish separately at every lower size.

This is not an N13 accident.  Every Gaussian-ideal period matrix has

```text
P(a,b)=[[a,-b],[b,a]],
R=[[0,-1],[1,0]],
R P=P R.
```

The quarter-turn `R` is a quotient-graph automorphism and exchanges the two
lines inside each orbit.  Since both `A` and `D` have `R` eigenvalue `-1`, any
unweighted scalar site ensemble on this Gaussian-ideal geometry has zero
`A/D` state, birth and exit expectation coefficientwise.

Therefore the roughly two-sigma N65 diagonal-odd marginal in `a7cb19a` is now
assigned to an exact symmetry-null direction, not retained as a candidate
unweighted physical response.  A nonzero `A/D` experiment must explicitly
insert an `R`-odd/charged source, or move to a quotient geometry that does not
have this automorphism.  This is precisely where an explicit defect source
becomes scientifically nonredundant.

## Boundary

The bridge is an exact finite-volume identity and a general Gaussian-symmetry
selection rule.  It does not identify the continuum field carried by the
surviving H current, nor predict its scaling exponent.

Reproduce with:

```bash
python3 scripts/p337_n13_twist_flux_bridge.py \
  --json results/p337-n13-twist-flux-bridge/latest.json \
  --markdown results/p337-n13-twist-flux-bridge/latest.md

python3 -m unittest discover -s tests \
  -p 'test_p337_n13_twist_flux_bridge.py'
```
