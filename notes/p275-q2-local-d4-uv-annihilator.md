# P275 paired-q2 local D4 selector

The nine-geometry global-line experiment rejected its scaling map because
`chi4(ell)/B` is a scale-zero projective H1 polarization.  This follow-up
removes `ell` completely.  It reads a local landing character around the same
site that creates each rank birth.

The implementation preserves the legacy R=1 field from `6899b11` and adds the
Euclidean-disk landing convention used by PR #247 at R=2 and R=4.  Both disks
are injective on all four frozen N65/N130 Gaussian quotients.  Each preinsertion
microcanonical row stores `S,D,qS,qD,I02` at both radii, so delete-one scoring can
retain the complete common-field covariance and independently verify the exact
rank-gate contact identities from `8bef10b`, including the local mark carried by
a direct `0->2` birth.

## Frozen mechanism selector

The primary is the mean odd source

```text
U = [local_D4(R=4)-local_D4(R=2)]/log(2).
```

The shell removes a radius-independent UV/contact term.  It is never replaced
by `Cov(A_top,U)`: that connected channel is exact contact algebra.  The even
`local_S4` shell is the explicit thermal-gate nuisance.

The two N65 raw orientations train exactly two coefficients: one common
thermal-shell coefficient and one H4 amplitude multiplying the registered
`cos(4 theta)`.  Nothing is tuned on N130.  The two N130 orientations are then
predicted with the fixed q2 transfer `-2^(-13/8)`, giving a two-coordinate
heldout score.  Recomputing roots, shells, nuisance coefficients and residuals
inside delete-one blocks propagates the training uncertainty.

## Prereveal gate

ARM64 GCC 10.3.1 produced binary
`12f3fc0daf7709f518f812038d5a991fc532841be01b2e64fdd1c08fe95a5c77`
from runner commit `6ecb339f93b878389301a1dc978ae4a38c522b5c`.  The exhaustive
R2 oracle, all four R2/R4 injectivity checks, and the two contact identities
passed.  On ZyTrST, 10,000 replicas took 1.344 seconds at N65 and 2.589 seconds
at N130 with eight threads each.  The now-authorized production target remains
one million replicas per size in disjoint seed/counter domains.

## Scientific card

1. MECHANISM SPACE: local landing D4 after an R4-R2 UV annihilator, with even thermal-shell nuisance.
2. NOT PROVED: survival nominates a local H4 covector but does not identify `Q4 epsilon`.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: no primary global observer | Alexander gate odd | mean local D4 | exact q2 Gaussian parent/child.
4. DEPENDENCY GROUP: orientations/radii/channels share batches within N; N65 and N130 use disjoint seed/counter domains.
5. UPWEIGHT OBSERVATION: the preregistered two-coordinate N130 residual after N65-only nuisance training.
