# P250 common-counter radius selector freeze

Before evaluating any new radius, freeze `R={1,2,3,4}` for the existing N325
norm-5 axis-minus-diagonal landing-H4 insertion.  Every candidate uses the same
4,000 replicas, batch boundaries, translation schedule, probability and seed.
`R=1` is reused from `c5bd025`; the remaining radii are evaluated once.

The selector does not optimize a phase p-value.  A separation is usable only
when all four charged two-point denominators have `|z|>=2` and the aligned
eight-real local-variance-normalized cubic covariance is numerically
nondegenerate.  A candidate advances only if at least two of `d={1,2,3}` are
usable.  If more than one advances, choose the first radius in the frozen
order, preserving the smallest local modification.

This is a geometry screen within one local insertion family.  It cannot reject
a continuum charged OPE or a different leg-defect/charged insertion.
