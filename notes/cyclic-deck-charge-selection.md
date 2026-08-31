# Cyclic deck-charge selection certificate

This certificate closes a bounded algebraic slice of parent issue #244.  For a
cyclic deck group `C_q`, a product of character rows with charges `r_i`
transforms by charge `sum r_i mod q`.  Averaging over a complete deck orbit
therefore vanishes exactly unless the total charge is zero.

The checked contract specializes this rule to the active `C2` and `C5` cover
fibers.  It proves:

- every nontrivial linear detail row is absent from an invariant scalar;
- the `C2` quadratic detail row is algebraically allowed;
- at norm five only conjugate pairs `(1,4)` and `(2,3)` are neutral;
- a charge-4 marked observable paired with a charge-1 score is neutral;
- five copies of any nontrivial `C5` charge are the first same-charge neutral
  tensor power.

All cancellations use exact modular exponents and cyclotomic orbit sums, not
floating roots of unity.

## Boundary

Charge neutrality is a necessary selection rule, not evidence that an allowed
channel has nonzero lattice overlap.  The certificate reads no measured
response and makes no production or operator-identification claim.  Parent
issue #244 remains open.
