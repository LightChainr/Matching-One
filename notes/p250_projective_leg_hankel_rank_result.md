# P250 result: five states survive within each hand, not across hands

The model-free Hankel test localizes the earlier shared-rank failure.  Each
plus/minus hand separately is compatible with a five-dimensional commuting
moment state.  A single raw state shared across both hands is not.

## Frozen rank table

The table reports finite-400-batch Hotelling probabilities for the complete
Schur/vanishing-minor chart.  The test assumes no exponential roots and
includes defective and Jordan commuting transfers.

| group | rank<=3 | rank<=4 | rank<=5 | lower bound at 0.01 |
|---|---:|---:|---:|---:|
| plus charge 1 | `0.117` | `0.693` | `0.400` | 3 |
| plus charge 2 | `1.22e-7` | `1.71e-4` | `0.389` | 5 |
| minus charge 1 | `0.215` | `0.114` | `0.0529` | 3 |
| minus charge 2 | `7.24e-7` | `7.00e-6` | `0.0126` | 5 |
| plus two-charge block | `1.67e-36` | `2.80e-11` | `0.0543` | 5 |
| minus two-charge block | `3.83e-17` | `0.00235` | `0.0655` | 5 |
| all-channel shared block | `4.04e-137` | `1.80e-13` | `2.84e-5` | 6 |

The result explains the prior free rank-five miss without invoking magnetic
translation.  Its shared roots forced plus and minus into one state.  Once the
hands are separated, rank five survives for both.  The information gradient
therefore points to **sector-specific five-state quotients and a nontrivial
Alexander/complement map between them**, not to a single rank-six object and
not to noncommuting Tx/Ty.

The charge asymmetry inside each hand is also structured: charge 1 needs only
rank at least three at this resolution, while charge 2 needs at least five.
The two-charge block needs five, so charge 2 supplies the extra observable
directions rather than merely repeating charge 1.

## Flat-extension boundary

Order-one flatness at rank at most three fails for both hand blocks.  That is a
state-count result.  Order-two flatness cannot be tested: a degree-three
monomial matrix requires moments through total degree six, while the current
diamond ends at degree four.

The endpoint decomposition gate is exactly zero, but this is constructional:
the stream records one `G(a,b)` and reuses it for every decomposition
`u+v=(a,b)`.  It never records ordered `xy` and `yx` histories separately.
Consequently current data neither show nor refute path/context memory.

## Scientific card

- **Mechanism space changed:** one common raw plus/minus transfer state of
  dimension at most five is removed; hand-specific dimension-five states
  remain compatible.
- **Not proved:** exact rank five, a unique cross-hand dimension six, a Jordan
  block, noncommutation, or path dependence.
- **Observer/sector/source/geometry:** C4-gauged neutral projective-leg pair,
  charges 1/2, plus/minus Gaussian children, N505 radius-four diamond.
- **Dependency group:** the same fresh seed `25050510120261130`, 80k samples,
  400 batches; no new simulation.
- **Next discriminator:** extract each hand block's covariance-aware
  degree-two annihilating null line and test whether Alexander/C4 maps the plus
  polynomial to the minus polynomial.  This tests the sector bridge without
  forcing identical raw eigenpairs.
