# Spin-4 direction-character alias gate

This is exact direction-character algebra, not percolation data.

| direction orbit | `exp(-4 i theta)` | scalar/spin-4 character rank |
|---|---:|---:|
| axial | `['1', '0']` | `1` |
| diagonal | `['-1', '0']` | `1` |

A constant defect of `7` therefore has naive axial spin-4 average `[7, 0]`, exactly equal to
its scalar average. On the diagonal orbit the same readout is only the negative scalar average.

Using both orbit averages gives the exact response matrix

`[['1', '1'], ['1', '-1']]`

with determinant `-2` and rank `2`. Hence

- scalar: `(axis_average + diagonal_average)/2`;
- spin 4: `(axis_average - diagonal_average)/2`.

## Decision

A direction-only defect on one C4 orbit cannot distinguish scalar from spin 4. At least two C4 orbits with unequal exp(-4 i theta0) phases, or additional typed/internal edge information, are required.

## Boundary

This gate does not implement Zhou edge observables, prove discrete harmonicity, measure percolation, establish L^-2 scaling, or identify a matching/KdV operator.
