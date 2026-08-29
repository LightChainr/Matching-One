# PR267 Target 2: named non-A_top complex C3 sector

The selected observer is the full primitive-homology automorphic harmonic

```text
delta H4(tau) = E[1_rank1 (P ell/|P ell|)^4] - H4_Pinson/Arguin(tau).
```

It is explicitly not `A_top=P2-P0`: it varies inside the rank-one sector and
retains the ambient winding line.  The physical lifted vector `P ell` fixes
one laboratory complex frame without using raw period coordinates.

The exact alias gate `83e98fc` does not invalidate this row.  That gate has
rank one only for a direction statistic restricted to a single C4 orbit.
Here the sum runs over every primitive winding line.  In each frozen child,
the two included lines `(1,0)` and `(0,1)` already have unequal exact fourth-
angle phases, so the scalar/spin-4 character matrix has rank two before the
remaining primitive orbits are added.  The runner now checks this fail-fast.

The three index-two children of the N56 positive Pell parent are frozen in
the exact Hecke order `2 tau`, `tau/2`, `(tau+1)/2`.  The runner gives all
three the same counter-derived 224-bit Bernoulli bond vector, so the six real
coordinates have one measured covariance block.  Pinson/Arguin baselines are
evaluated at the actual finite Pell-child moduli before any character score.

The three pure competitors are not ordinary H4/H8/H12 angular harmonics.
They are the low-weight modular-ring supports `E4` (C3 r=1), `E6` (r=0), and
`E4^2` (r=2), transported to the actual child moduli.  The null and every
two-character mixture are frozen opponents.  A mixture outcome says that
this finite topology-typed row is not a single ordinary modular-ring
character; it does not refute a continuum field family.

## Scientific card

1. **Question:** which modular C3 support carries a named non-`A_top` rank-one polarization response?
2. **Observer:** continuum-subtracted full primitive winding `H4`, not a global rank imbalance.
3. **Geometry:** all three degree-two children of one N56 rho Pell parent, one physical complex frame.
4. **Dependency:** one common 224-bit field and a full 6x6 covariance; no coordinate is counted as a separate vote.
5. **Decision:** compare null, three pure modular supports and frozen two-support mixtures; never reopen H4/H8/H12 ordinary-channel voting.

## Independent production outcome

The preregistered three-machine run used 2,000,000 samples per child with
independent seeds and disjoint counter intervals.  The primary clean E4/r1
row is decisively rejected (`chi2=235.2105/4`, `p=9.97e-50`); E6/r0 and
E4-squared/r2 are also decisively rejected.  The secondary E4+E6 subspace is
the unique near-fit but narrowly fails the frozen `.05` criterion
(`chi2=6.0645/2`, `p=.04821`).

The more informative mechanism statement is character support: the ideal C3
r0 and r1 rows are independently resolved, while r2 is not.  Thus the named
topology-typed observable is genuinely multi-character at N112.  This closes
the proposed clean-row interpretation without reopening ordinary angular
H4/H8/H12 selection.
