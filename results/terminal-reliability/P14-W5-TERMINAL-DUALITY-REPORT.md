# P14 W5 two-orbit terminal-duality certificate

## Decision

The structurally gated W5 candidate does **not** close as a four-terminal
disk-self-dual cell.

There is a genuine but weaker fact: the spherical planar dual is another W5,
exchanges rim and spoke orbits and gives the complement map

```text
(r,s) -> (1-s,1-r),    fixed line r+s=1.
```

That fact does not induce a planar-complement permutation of the W5 terminal
partition vector.  The correct disk-relative dual has four extra boundary
terminals and is not W5.

## Exact two-parameter distribution

All `2^8=256` edge configurations were enumerated in edge order
`R_0,...,R_3,S_0,...,S_3`.  For each of the 15 four-terminal set partitions,
the certificate stores the exact table

```text
C_pi[a,b] r^a (1-r)^(4-a) s^b (1-s)^(4-b),  0<=a,b<=4.
```

Every `(a,b)` table slice sums to `binom(4,a) binom(4,b)`.  Exactly one
partition is identically absent: `01|23`, the crossing pairing for cyclic
terminal order `(0,2,1,3)`.  The full tables and all 256 configuration rows
are in `p14-w5-terminal-duality.json`.

## Why spherical self-duality is insufficient

Write the four triangular faces as `F_i` and the outer face as `O`.  The
spherical dual has

```text
R_i* = (O,F_i),
S_i* = (F_(i-1),F_i).
```

The reflection `O -> h, F_i -> q_-i` makes complement-duality an exact
involution on the eight edge bits.  It maps open-edge counts as

```text
(a,b) -> (4-b,4-a)
```

and preserves every configuration weight on `r+s=1`.

However, its output terminal partition is not a function of the input
terminal partition.  The frozen minimal witness is:

```text
bits          primal partition   spherical-dual partition
00000011      0|13|2             0123
00100011      0|13|2             023|1
```

The second configuration adds `R_2` to open spokes `S_2,S_3`.  This closes an
internal triangle without changing the primal terminal partition, but it
changes the spherical-dual terminal partition.  Thus an internal cycle datum
lost by the 15-state terminal projection is already decisive at three open
edges.

## The actual boundary dual

For disk terminal duality the outer face must be split into four boundary
dual terminals `B_i`, one in each boundary arc.  The relative dual is

```text
R_i* = (F_i,B_i),
S_i* = (F_(i-1),F_i).
```

It is a four-cycle of internal face vertices with one terminal leaf at each
vertex.  Configuration by configuration, its terminal partition is a unique
planar/Kreweras complement of the primal partition; in particular

```text
0123 <-> 0|1|2|3.
```

But it cannot be boundary-terminal-preserving isomorphic to W5:

```text
                         W5 primal       disk-relative dual
vertex count             5               8
degree multiset           3,3,3,3,4       1,1,1,1,3,3,3,3
terminal degrees          3,3,3,3         1,1,1,1
```

This is the smallest structural obstruction; no probability choice repairs
it.

As a scalar check rather than an assumed criterion, substitution `s=1-r`
gives

```text
P(0123)-P(0|1|2|3)
= 1 - 4r + 14r^2 - 36r^3 + 76r^4 - 112r^5
  + 108r^6 - 60r^7 + 14r^8,
```

which is not an identity.  At `r=s=1/2` it equals exactly `67/128`.

## Correct continuation for Issue #14

The next exact finite-cell object is the joint primal/disk-relative-dual
boundary connectivity, equivalently the alternating-boundary medial or
Temperley--Lieb connectivity state.  It must compare W5 with the four-cycle
plus terminal leaves, not identify the spherical face centres with boundary
dual terminals.  A rigorous-bound route would then still need an explicit
periodic tiling and stochastic/comparison map.

This certificate makes no threshold, critical-polynomial, periodic tiling or
rigorous-bound claim.

## Reproduction

```bash
python3 scripts/p14_w5_terminal_duality.py \
  --output results/terminal-reliability/p14-w5-terminal-duality.json
uv run --with pytest python -m pytest -q \
  tests/test_p14_w5_terminal_duality.py \
  tests/test_p14_four_terminal_balance_roots.py
```
