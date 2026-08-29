# Discrete-holomorphic spin-4 alias gate

Status: exact observable-design slice of Issue 109.

For a direction-only scalar defect sampled on one C4 orbit,

```text
theta_j = theta_0 + j*pi/2,
exp(-4 i theta_j) = exp(-4 i theta_0).
```

The nominal spin-4 character is therefore constant on that orbit. On the four axial directions it is
`+1`, so the naive fourth angular moment is exactly the scalar sum. On the diagonal orbit it is `-1`,
which only changes the sign of the same scalar sum. Translation or reflection within either orbit cannot
restore identifiability.

The smallest direction-only repair uses two C4 orbits with different fourth-angle phases. Axis and
diagonal averages obey

```text
[A_axis]     [1  1] [a_scalar]
[A_diag ] =  [1 -1] [a_spin4 ].
```

The matrix has determinant `-2` and rank two, giving

```text
a_scalar = (A_axis + A_diag)/2,
a_spin4  = (A_axis - A_diag)/2.
```

This is a fail-fast design constraint, not an implementation of Zhou's two edge observables. A typed or
internally complex edge observable may carry additional transformation data and must be analyzed from its
actual definition. No percolation result, decay law, or matching/KdV identification follows here.
