# Frozen radius-one contact residual extension

Status: `FROZEN_BEFORE_OUTPUT`

This extension is applied only to the already frozen kernel-changing diagonal
edge set.  It partitions the signed mass by the producer's four fixed contact
masks.  Mask `0` means that neither source cut contacts either local occupied
thermal arm; masks `1/2` contact one arm and mask `3` contacts both.

There is one decision:

- exact nonzero mask-0 Schur sum:
  `RADIUS_ONE_CONTACT_ONLY_CLOSURE_REJECTED`;
- otherwise: `CONTACT_ZERO_RESIDUAL_UNRESOLVED`.

No contact radius, corner word, source orbit, rank stage, root coefficient, or
sample is selected after output.  The full 2x6 mask-0 matrix and all four mask
sums are retained so the surviving functional can be read without another
scan.
