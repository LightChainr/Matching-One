# N112 rho-child Alexander-even complex C3 production

The fresh common-field production used 2,000,000 configurations in 100 batches
on Zy and completed in 223.15 seconds.  All three homology invariant-failure
counts are zero and the production stderr is empty.

## Primary result

The continuum-subtracted Alexander-even coordinate has a strong nontrivial
child character:

```text
delta E_top,r1 = -0.001744631 - 0.000047198 i
chi-square = 142.199 / 2
p = 1.324e-31
```

The scalar row is smaller and does not pass the frozen 0.01 gate:

```text
delta E_top,r0 = -0.000506428
chi-square = 5.716 / 1
p = 0.01681
```

This is a resolved Alexander-even rank-redistribution C3 response, not an
`A_top` signal and not another H4/H8/H12 vote.

## Independent observer direction

The same stream retains primitive H4 only to test whether the two named
observers occupy one child-character ray.  The amplitude-free determinant is

```text
D = H4_r0 Etop_r1 - H4_r1 Etop_r0
  = -3.29715e-6 - 4.79724e-7 i
chi-square = 21.202 / 2
p = 2.489e-5
```

Thus a common ray is rejected at the frozen threshold.  Primitive rank-one
polarization and Alexander-even rank redistribution are distinct topology
observer directions on the rho-child triple.

## Boundary

This result concerns one finite N112 square-bond common-field block.  It does
not name a continuum field, fit an exponent, identify a unique state count, or
transport the result to square-site matching.  The strong claim is the exact
observer/sector distinction plus its covariance-scored finite-lattice readout.

The first remote launch stopped during import because the fresh image lacked
`mpmath`; no replica was generated.  Its log is preserved.  The identical
frozen stream completed after installing `mpmath 1.4.1`.
