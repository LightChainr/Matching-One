# P250 result: an Alexander bridge survives, but is not identified

The two hand-specific quadratic annihilator lines can be transported into one
another by every frozen Alexander reflection-plus-conjugation convention.  No
new simulation was used.

| primary map | finite-batch Hotelling p |
|---|---:|
| reflection + conjugation + `R^0` | `0.1163` |
| reflection + conjugation + `R^1` | `0.3863` |
| reflection + conjugation + `R^2` | `0.3542` |
| reflection + conjugation + `R^3` | `0.0612` |

The internal C4 calibration is also clean: all six within-hand conjugating
line maps have `p=0.690..0.991`.  The C4-covariant gauge is therefore behaving
as required at the level of the truncated annihilator.

This is compatibility, not an identification.  The orientation-preserving
identity/conjugation comparator has `p=0.986`, and most of the other declared
maps also survive.  The current null-line covariance is too broad to decide
which intertwiner is physical.  Selecting the largest p-value would simply
rename uncertainty.

## What shared rank six actually meant

Each hand block has six degree-two columns and one compatible annihilator line:

```text
ker(H_plus)  = <q_plus>,
ker(H_minus) = <q_minus>.
```

The raw shared block is their vertical stack.  Its kernel is the intersection
of those two lines.  Unless `q_plus` and `q_minus` are literally the same in
the identity column chart, the intersection is zero and the six-column stack
has rank six.  Thus the prior shared lower bound was a valid rejection of the
**identity intertwiner**, but it was not a physical six-state lower bound.

The minimal truncated two-sector realization is instead

```text
Q_plus  = C^6 / <q_plus>   (dimension 5),
Q_minus = C^6 / <q_minus>  (dimension 5),
Q = Q_plus direct-sum Q_minus  (dimension 10),
```

with Alexander reflection/conjugation retained as a candidate anti-linear
intertwiner between the two five-dimensional sectors.  Deck-generator changes
do not create extra candidates: charge permutations and fifth-root phases are
invertible row operations and leave the common right line unchanged.

The exact finite geometries remain distinct.  The C4 fiber multipliers `3`
and `2` are inverse modulo five, but `(10+i)(2+i)` and `(10+i)(2-i)` are not D4
reflections of one another.  Any surviving bridge is therefore a truncated
Hecke-sector morphism, not an exact graph isomorphism.

## Scientific card

- **Mechanism space changed:** the raw shared-rank-six result is reclassified
  as failure of the identity sector map; a `5+5` two-sector realization with
  an Alexander-compatible bridge survives.
- **Not proved:** a unique Alexander intertwiner, exact dimension five beyond
  the truncation, or a closed multiplication algebra.
- **Observer/sector/source/geometry:** degree-two projective-leg annihilator;
  charge-1/2 blocks, plus/minus N505 Gaussian children, C4-covariant gauge.
- **Dependency group:** the same 80k/400-batch fresh stream; all mapping scores
  are correlated reanalyses, not independent votes.
- **Next discriminator:** collect the C4-closed Manhattan radius-five shell.
  Its six new first-quadrant degree-five endpoints extend every degree-three
  shift of both annihilators.  Add degree six only when testing the full next
  flat-extension plateau.
