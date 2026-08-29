# A local complement-odd pivotal readout for the C4 tangent

## Result

The global wrapping rows tested in PR #208 were statistically rank one at
`N=130` and `N=170`.  A different readout is therefore needed before adding
more samples.  The smallest local candidate is

```text
O_local = [H4_pivotal(black) - H4_pivotal(complement)] / 2,
```

where `H4_pivotal` is the fixed-even-root cross-pivotal indicator marked by
axis-minus-diagonal four-arm landings on a square annulus.  The checkerboard
triangulation is its own matching graph, so occupation complement exchanges
the two terms and makes `O_local` exactly matching odd.  A 45-degree change of
the landing registry flips its sign.

The complete `N=10`, `(a,b)=(3,1)`, radius-one enumeration gives the response
matrix with rows `(global cross half-difference, O_local)` and columns
`(t, lambda)`:

```text
R_10 = [[ 15/8,   5/4  ],
        [-3/64,  11/64 ]]

det(R_10) = 195/512.
```

Thus the local row is exactly linearly independent of the global wrapping
row.  Among all 1024 configurations the twice-observable has counts
`{-1: 88, 0: 848, +1: 88}`.  Complement oddness, the exact C4 rotation, and
the 45-degree registry sign rule have zero violations.  Reflection is not
claimed inside the chiral `(3,1)` quotient; it maps to the conjugate period
representation and belongs in a future paired-geometry check.

## Scientific use

This closes only the microscopic observability gate.  It does not show that a
second continuum RG direction remains resolved at large size.  The next run
should evaluate the full `2 x 2` score-response matrix at `N=130` and `N=170`
inside aligned batches and recompute its generalized eigenvalues and
condition number in each delete-one replicate.  Continue to a third size only
if the second singular direction is statistically resolved.

The primary value is structural: this observable directly crosses the
rank-one bottleneck of the previously tested topological projections without
changing the exact `(t,lambda)` UV tangent or assuming the `x=21/4` answer.

## Reproduction

```bash
python3 scripts/c4_local_odd_pivotal.py \
  --json results/local-20260829/P155-c4-local-odd-pivotal/exact-n10-r1.json
python3 -m unittest tests.test_c4_local_odd_pivotal -v
```
