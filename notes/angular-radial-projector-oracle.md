# Angular/radial projector commutation oracle

This bounded Issue #8 gate treats orientation and scale as separate two-dimensional
vector spaces.  The angular scalar row `(1/2, 1/2)` kills the declared H4
orientation `(1, -1)`.  Independently, the radial row `(-1/8, 1)` kills the
declared H4 size response `(8, 1)`.

The checked-in Fraction oracle verifies that applying these rows in either order
gives the same scalar on a pure H4 sector, a pure scalar sector, and their sum.
This is the finite-dimensional tensor identity `(a^T tensor r^T)X`; it does not
depend on floating-point tolerances.  Both individual filters kill the pure
separable H4 sector, so double filtering is exactly redundant there.  A mixed
sector retains a nonzero scalar residue, preventing that narrow statement from
being misread as equivalence on arbitrary data.

Run:

```text
python3 scripts/angular_radial_projector_oracle.py
python3 -m unittest tests/test_angular_radial_projector_oracle.py -v
```

## Boundary

The ratios in the contract are synthetic exact coordinates, not measured
finite-size amplitudes.  This gate does not show that a production observable is
separable or pure H4, does not estimate 13/8, and does not replace covariance-
aware projector construction on real data.  Issue #8 therefore remains open.
